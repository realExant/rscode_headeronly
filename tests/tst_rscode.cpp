#include <catch2/catch_test_macros.hpp>
#include "../rscode_headeronly.h"

template<typename T = uint8_t, uint32_t t = 2, bool gf_static = true>
bool tst_rscoder(const std::vector<T>& data,
                 const std::vector<T>& noise,
                 const std::vector<T>& tst_remainder,
                 const std::vector<T>& tst_syndromes,
                 const std::vector<T>& tst_sigma,
                 const std::vector<T>& tst_lambda,
                 const std::vector<T>& tst_locator,
                 const std::vector<T>& tst_error) {
    std::vector<T> remainder;
    std::vector<T> syndromes;
    std::vector<T> sigma;
    std::vector<T> lambda;
    std::vector<T> locator;
    std::vector<T> error;

    const auto coder{ ZLib::ErrCorr::RSCoder<T, 2, true>() };

    auto encoded{ coder.encode(data.begin(), data.end(), &remainder) };

    for (size_t i = 0; i < noise.size(); ++i)
        encoded[i] = noise[i];

    const auto& decoded{ coder.decode(encoded, &syndromes, &sigma, &lambda, &locator, &error) };

    return (data == decoded) && (remainder == tst_remainder) && (syndromes == tst_syndromes) &&
           (sigma == tst_sigma) && (lambda == tst_lambda) && (locator == tst_locator) && (error == tst_error);
}

TEST_CASE("T=uint8 t=2") {
    using type = uint8_t;

    const std::vector<type> data{ 240, 229, 234, 224, 213 };
    const std::vector<type> noise{ 203, 236 };
    const std::vector<type> remainder{ 228, 81, 61, 5 };
    const std::vector<type> syndromes{ 1, 72, 225, 174, 48 };
    const std::vector<type> sigma{ 1, 3, 2 };
    const std::vector<type> lambda{ 1, 75, 59 };
    const std::vector<type> locator{ 0, 1 };
    const std::vector<type> error{ 47, 189 };

    SECTION("gf_static=true") { REQUIRE(tst_rscoder<type, 2, true>(data, noise, remainder, syndromes, sigma, lambda, locator, error)); }
    SECTION("gf_static=false") { REQUIRE(tst_rscoder<type, 2, false>(data, noise, remainder, syndromes, sigma, lambda, locator, error)); }
}
