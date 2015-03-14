// Simple byte-aligned binary arithmetic coder (Ilya Muravyov's variant) - public domain - Fabian 'ryg' Giesen 2015
//
// Written for clarity not speed!

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>

// Probabilities are expressed in fixed point, with kProbBits bits of
// resolution. No need to go overboard with this.
static int const kProbBits = 12;
static uint32_t const kProbMax = 1u << kProbBits;

// Type used for buffers.
typedef std::vector<uint8_t> ByteVec;

// Binary arithmetic encoder (Ilya Muravyov's variant)
// Encodes/decodes a string of binary (0/1) events with
// probabilities that are not 1/2.
//
// This code is written for clarity, not performance.
class BinArithEncoder
{
    uint32_t lo, hi;
    ByteVec &bytes;

    // noncopyable
    BinArithEncoder(BinArithEncoder const &);
    BinArithEncoder &operator =(BinArithEncoder const &);

public:
    // Initialize
    explicit BinArithEncoder(ByteVec &target) : lo(0), hi(~0u), bytes(target) { }

    // Finish encoding - flushes remaining codeword
    ~BinArithEncoder()
    {
        for (int i = 0; i < 4; ++i)
        {
            bytes.push_back(lo >> 24);
            lo <<= 8;
        }
    }

    // Encode a binary symbol "bit" with the probability of a 1 being "prob".
    // Note that prob=0 (or prob=1<<kProbBits) really mean that a 1 (or 0,
    // respectively) cannot occur!
    void encode(int bit, uint32_t prob)
    {
        // Midpoint of active probability interval subdivided via prob
        uint32_t x = lo + ((uint64_t(hi - lo) * prob) >> kProbBits);

        if (bit)
            hi = x;
        else
            lo = x + 1;

        // Renormalize: when top byte of lo/hi is same, shift it out.
        while ((lo ^ hi) < (1u << 24))
        {
            bytes.push_back(lo >> 24);
            lo <<= 8;
            hi = (hi << 8) | 0xff;
        }
    }
};

// Corresponding decoder.
class BinArithDecoder
{
    uint32_t code, lo, hi;
    ByteVec const &bytes;
    size_t read_pos;

    // noncopyable
    BinArithDecoder(BinArithDecoder const &);
    BinArithDecoder &operator =(BinArithDecoder const &);

public:
    // Start decoding
    explicit BinArithDecoder(ByteVec const &source)
        : lo(0), hi(~0u), bytes(source), read_pos(0)
    {
        code = 0;
        for (int i = 0; i < 4; ++i)
            code = (code << 8) | bytes[read_pos++];
    }

    // Decode a binary symbol with the probability of a 1 being "prob".
    int decode(uint32_t prob)
    {
        int bit;

        // Midpoint of active probability interval subdivided via prob
        uint32_t x = lo + ((uint64_t(hi - lo) * prob) >> kProbBits);

        if (code <= x)
        {
            hi = x;
            bit = 1;
        }
        else
        {
            lo = x + 1;
            bit = 0;
        }

        // Renormalize
        while ((lo ^ hi) < (1u << 24))
        {
            code = (code << 8) | bytes[read_pos++];
            lo <<= 8;
            hi = (hi << 8) | 0xff;
        }

        return bit;
    }
};

// ---- A few basic models

// NOTE: Again, this is written for clarity and ease of tinkering.
// In practice, you will write more direct code for these once you've
// figured out your coding structure.

// Adaptive binary model. These are pretty good!
// Lower Inertia = faster.
//
// You typically build more sophisticated models out of these
// by having lots of them and choosing the active model based on
// context.
template<int Inertia>
struct BinShiftModel
{
    uint16_t prob;

    BinShiftModel() : prob(kProbMax / 2) {}

    void encode(BinArithEncoder &enc, int bit)
    {
        enc.encode(bit, prob);
        adapt(bit);
    }

    int decode(BinArithDecoder &dec)
    {
        int bit = dec.decode(prob);
        adapt(bit);
        return bit;
    }

    void adapt(int bit)
    {
        // Note prob never his 0 or kProbMax with this update rule!
        if (bit)
            prob += (kProbMax - prob) >> Inertia;
        else
            prob -= prob >> Inertia;
    }
};

// BitTree model. A tree-shaped cascade of BinShiftModels.
// This is the de-facto standard way to build a multi-symbol coder
// (values with NumBits bits) out of binary models.
//
// LZMA (as in 7zip/xz) uses this type of model (backed by a BinShiftModel
// as above) for its literals.
template<typename BitModel, int NumBits>
struct BitTreeModel
{
    static size_t const kNumSyms = 1 << NumBits;
    static size_t const kMSB = kNumSyms / 2;

    BitModel model[kNumSyms - 1];

    void encode(BinArithEncoder &enc, size_t value)
    {
        assert(value < kNumSyms);

        // The first bit sent is the MSB of the value and coded without context
        // Second bit is the bit below the MSB, using the value of the MSB as context
        // and so forth.
        //
        // 1 + 2 + 4 + ... = 2^NumBits - 1 contexts.
        // Numbering the MSB context 1 and then shifting in the coded bits from the
        // bottom is a convenient way to index them. (So ctx is 1-based)
        size_t ctx = 1;
        while (ctx < kNumSyms)
        {
            int bit = (value & kMSB) != 0;
            value += value; // shift value by 1 for next iter
            model[ctx - 1].encode(enc, bit);
            ctx += ctx + bit; // shift in "bit" into context
        }
    }

    size_t decode(BinArithDecoder &dec)
    {
        // Corresponding decoder is nice and easy:
        size_t ctx = 1;
        while (ctx < kNumSyms)
            ctx += ctx + model[ctx - 1].decode(dec);

        return ctx - kNumSyms;
    }
};

// ---- Random utility code

static double log_2(double x)
{
    return log(x) / log(2.0);
}

// ---- Some examples

static void example_static()
{
    // A static binary source with known probability of 1 being 1/5.
    ByteVec source;
    uint32_t const kProbOne = kProbMax / 5;

    srand(1234);
    for (size_t i = 0; i < 10000; ++i)
        source.push_back(rand() < (RAND_MAX/5));

    // Encode it
    ByteVec coded;
    {
        BinArithEncoder coder(coded);
        for (size_t i = 0; i < source.size(); ++i)
            coder.encode(source[i], kProbOne);
    }

    // Print actual and expected size (based on order-0 entropy)
    {
        double p = kProbOne / (double)kProbMax;
        double entropy_bits_per_sym = -p * log_2(p) - (1.0 - p) * log_2(1.0 - p);
        printf("static size: %d bytes - entropy: %.2f bytes\n", coded.size(), source.size() * entropy_bits_per_sym / 8.0);
    }

    // Decode it
    ByteVec decoded;
    {
        BinArithDecoder coder(coded);
        for (size_t i = 0; i < source.size(); ++i)
            decoded.push_back((uint8_t) coder.decode(kProbOne));
    }

    if (decoded != source)
        printf("error decoding!\n");
    else
        printf("decodes ok!\n");
}

static void example_dynamic()
{
    // A binary source that keeps changing its probability of 1 regularly
    // in a way opaque to the coder.
    // Use this as example for an adaptive model.
    static int const kInertia = 4;
    ByteVec source;

    srand(2345);
    for (size_t chunk = 0; chunk < 50; ++chunk)
    {
        int threshold = rand();
        for (size_t i = 0; i < 200; ++i)
            source.push_back(rand() < threshold);
    }

    // Encode it
    ByteVec coded;
    {
        BinArithEncoder coder(coded);
        BinShiftModel<kInertia> model;
        for (size_t i = 0; i < source.size(); ++i)
            model.encode(coder, source[i]);
    }

    printf("dynamic size: %d bytes\n", coded.size());

    // Decode it
    ByteVec decoded;
    {
        BinArithDecoder coder(coded);
        BinShiftModel<kInertia> model;
        for (size_t i = 0; i < source.size(); ++i)
            decoded.push_back((uint8_t) model.decode(coder));
    }

    if (decoded != source)
        printf("error decoding!\n");
    else
        printf("decodes ok!\n");
}

static void example_multisymbol()
{
    // Example for a multi-symbol alphabet - bytes in this case.
    // Let's get meta and use this source code as our source!
    typedef BitTreeModel<BinShiftModel<5>, 8> ByteModel;
    ByteVec source;

    {
        FILE *f = fopen("main.cpp", "rb");
        if (!f)
            return;

        fseek(f, 0, SEEK_END);
        source.resize(ftell(f));
        fseek(f, 0, SEEK_SET);
        fread(&source[0], 1, source.size(), f);
        fclose(f);
    }

    // Encode it
    ByteVec coded;
    {
        BinArithEncoder coder(coded);
        ByteModel model;
        for (size_t i = 0; i < source.size(); ++i)
            model.encode(coder, source[i]);
    }

    printf("multisymbol size: %d bytes\n", coded.size());

    // Decode it
    ByteVec decoded;
    {
        BinArithDecoder coder(coded);
        ByteModel model;
        for (size_t i = 0; i < source.size(); ++i)
            decoded.push_back((uint8_t) model.decode(coder));
    }

    if (decoded != source)
        printf("error decoding!\n");
    else
        printf("decodes ok!\n");
}

int main()
{
    example_static();
    example_dynamic();
    example_multisymbol();
    return 0;
}

// vim:et:sts=4:sw=4
