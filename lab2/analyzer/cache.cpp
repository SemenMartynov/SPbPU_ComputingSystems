#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <cmath>
#include <utility>
#include <stdexcept>
#include <iomanip>

class cache32
{
public:
    cache32(size_t cache_lines, size_t words_number, size_t associativity);
    bool check(size_t addr, signed char size);

    static const long cache_max_size = 512 * 1024 * 8; // in bits
    friend std::ostream& operator<<(std::ostream &os, const cache32& other);
private:
    const size_t cache_lines;
    const size_t words_number;
    const size_t associativity;

    // Cache is a vector of 'associativity' banks.
    // Each bank contains 'cache_lines' / 'associativity' lines.
    // Each cach line stores record with counter for LFU and address tag.
    typedef std::pair<signed char, size_t> cache_record;
    typedef std::vector<cache_record> cache_line;
    std::vector<cache_line> cache;

    bool isPowerOf2(size_t a) const;
    long wordsToBits(size_t a) const;
    long wordsToBytes(size_t a) const;
};

// It's for debug purpose
std::ostream& operator<<(std::ostream &os, const cache32& obj)
{
    std::vector<cache32::cache_line>::const_iterator bank;
    for (bank = obj.cache.begin(); bank != obj.cache.end(); ++bank)
    {
        os << "{ ";
        cache32::cache_line::const_iterator record;
        for (record = bank->begin(); record != bank->end(); ++record)
        {
            os << "(" << static_cast<size_t>(record->first) << "," << std::hex
                    << record->second << ") ";
        }
        os << "}" << std::endl;
    }
    return os;
}

cache32::cache32(size_t cache_lines, size_t words_number, size_t associativity) :
        cache_lines(cache_lines), words_number(words_number), associativity(
                associativity)
{
    // check correctness
    if (!cache_lines || !words_number || !associativity)
    {
        throw std::runtime_error(
                "Can't create a cache model -- one or more params are equal to zero");
    }
    if (!isPowerOf2(words_number))
    {
        throw std::runtime_error(
                "Number of words in the cache line in not a power of 2.");
    }
    if (!isPowerOf2(associativity))
    {
        throw std::runtime_error("Associativity index in not a power of 2.");
    }
    if (cache_lines % associativity != 0)
    {
        throw std::runtime_error(
                "Cache line number does not correlate with associativity index.");
    }

    // check cache size
    size_t tag_size = static_cast<size_t>(32 - log2(wordsToBits(words_number)));
    size_t flag_size = 8; // I hope 8 bits is enough for LFU
    long cache_size = (tag_size + wordsToBits(words_number) + flag_size)
            * cache_lines;
    if (cache_size > cache_max_size)
    {
        throw std::range_error("Cache size is greater than available space.");
    }

    // alocate cache
    cache.reserve(associativity);
    for (size_t bank = 0; bank != associativity; ++bank)
    {
        cache.push_back(
                cache_line(cache_lines / associativity, std::make_pair(0, 0)));
    }
}

bool cache32::isPowerOf2(size_t a) const
{
    return a > 0 && (a & (a - 1)) == 0;
}

long cache32::wordsToBits(size_t a) const
{
    // We assume, that 1 machine word is equal to 32 bits
    return a * 32;
}

long cache32::wordsToBytes(size_t a) const
{
    // We assume, that 1 machine word is equal to 4 bytes
    return a * 4;
}

bool cache32::check(size_t addr, signed char size)
{
    // line number in memory
    size_t N_memory = static_cast<size_t>(addr / wordsToBytes(words_number));

    // line number in cache (but we don't know the bank number)
    size_t N_bank_line = N_memory % (cache_lines / associativity);

    // Effective address (tag)
    size_t cache_line_tag = addr
            - addr % static_cast<size_t>(addr / wordsToBytes(words_number));

    // Check for cache hit
    for (std::vector<cache32::cache_line>::iterator bank = cache.begin();
            bank != cache.end(); ++bank)
    {
        if ((*bank)[N_bank_line].second == cache_line_tag)
        {
            // Cache hit!!
            // Increase LFU counter!
            ++(*bank)[N_bank_line].first;
            // Now we need to check size and either return result or split request!
            size_t cach_line_end = cache_line_tag
                    + static_cast<size_t>(wordsToBytes(words_number));
            size_t data_reqt_end = addr + static_cast<size_t>(size);
            return (data_reqt_end > cach_line_end) ?
                    check(cach_line_end,
                            static_cast<unsigned char>(data_reqt_end
                                    - cach_line_end)) :
                    true;
        }
    }

    // Cache miss!!
    // Update cache, search for minimum LFU counter.
    size_t lfu_num = 0;
    size_t lfu_val = static_cast<size_t>(cache[lfu_num][N_bank_line].first);
    for (size_t bank = 0; bank != associativity; ++bank)
    {
        if (static_cast<size_t>(cache[bank][N_bank_line].first) > lfu_val)
        {
            lfu_val = static_cast<size_t>(cache[bank][N_bank_line].first);
            lfu_num = bank;
        }
    }
    cache[lfu_num][N_bank_line] = std::make_pair(static_cast<signed char>(1),
            cache_line_tag);
    return false;
}

int main(int argc, char* argv[])
{
    /*
     * READ LOG
     */

    //log_record: arrd, len
    typedef std::pair<size_t, signed char> log_record;
    std::vector<log_record> log;
    log.reserve(48000);

    std::ifstream pinatrace("livermorec-mxopt.out");
    //std::ifstream pinatrace("livermorec-noopt.out");
    if (pinatrace.is_open())
    {
        while (pinatrace.good())
        {
            size_t itrash;
            char ctrash;
            size_t addr;
            size_t size;

            pinatrace >> std::hex >> itrash;
            pinatrace >> ctrash;
            pinatrace >> ctrash;
            pinatrace >> std::hex >> addr;
            pinatrace >> size;
            pinatrace >> std::hex >> itrash;

            log.push_back(std::make_pair(addr, static_cast<signed char>(size)));
        }
        pinatrace.close();
    }

    /*
     * Check loop
     */

    std::cout << "/---------------------------------------------------\\"
            << std::endl;
    std::cout << "| C.Lines |  Words  |  Assoc. | Miss ctr |  Hit ctr |"
            << std::endl;
    std::cout << "|---------|---------|---------|----------|----------|"
            << std::endl;

    for (size_t associativity = 2; associativity != 16; associativity *= 2)
    {
        for (size_t cache_lines = associativity * 2;
                cache_lines < associativity * 1024; cache_lines *= 2)
        {
            for (size_t words_number = 1; words_number != 16; words_number *= 2)
            {
                long long miss_counter = 0;
                long long hit_counter = 0;

                try
                {
                    // Create cache
                    cache32 cache(cache_lines, words_number, associativity);

                    // Begin processing
                    std::vector<log_record>::const_iterator record;
                    for (record = log.begin(); record != log.end(); ++record)
                    {
                        cache.check(record->first, record->second) ?
                                ++hit_counter : ++miss_counter;
                    }
                } catch (std::exception& ex)
                {
                    std::cout << "Error: " << ex.what() << std::endl;
                    return -1;
                }

                // Show results
                std::cout << "|" << std::setw(8) << cache_lines;
                std::cout << " |" << std::setw(8) << words_number;
                std::cout << " |" << std::setw(8) << associativity;
                std::cout << " |" << std::setw(9) << hit_counter;
                std::cout << " |" << std::setw(9) << miss_counter;
                std::cout << " |" << std::endl;

            }
        }
    }

    std::cout << "\\---------------------------------------------------/"
            << std::endl;
    return 0;
}
