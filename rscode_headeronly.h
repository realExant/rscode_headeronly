#pragma once

#include <algorithm>
#include <cassert>
#include <tuple>
#include <vector>

namespace ZLib::ErrCorr {
    /*!
    * \brief Класс реализует помехозащищенное кодирование кодом Рида-Соломона (РС).
    * Поле Галуа строится над примитивным членом 2.
    * \param[in] T минимальый размер кодового слова (uint8_t .. uint64_t).
    * \param[in] t количество кодовых слов для исправления в искаженном закодированном сообщении.
    * Размер закодированного сообщения n = k + 2*t, где k - длина исходного сообщения.
    * \param[in] gf_static указывает как вычисялются элементы поля Галуа:
    *	\li \c gf_static == true - элемент поля формируется динамически при каждом запросе, что увеличивает время работы алгоритма,
    *		но практически не накладывает расходы на память.
    *	\li \c gf_static == false - элемент поля запрашивается из статической таблицы сформированной на шаге инициализации класса,
    *		Это приводит к увеличению динамической памяти, но уменьшает время работы алгоритма.
    *		Размер таблиц прямого и обратношо полей Галуа равен 2 ^ (sizeof(T) * 8).
    * \code{.cpp}
    * 	//Данные для кодирования
    *	std::vector<uint8_t> data = { 240, 229, 234, 224, 213 };
    *
    *	//Кодирование данных кодом РС.
    *	//Определение кодера:
    *	// - 8 бин на кодовое слово. Поле Галуа на 256 слов (2 ^ 8);
    *	// - Возможность исправления до 2 кодовых слов;
    *	// - Прямая и обратная таблицы кодирования в поле Галуа предварительно вычислены.
    *	auto coder = errcorr::RSCoder<uint8_t, 2, true>();
    *
    *	//Кодирование данных
    *	auto encoded_data = coder.encode(data);
    *
    *	//Намереное искажение данных
    *	encoded_data[0] = 203; encoded_data[1] = 236;
    *
    *	//Декодирование кодом РС.
    *	auto decoded_data = coder.decode(encoded_data);
    * \endcode
    */
    template<typename T = uint8_t, uint32_t t = 2, bool gf_static = true>
    struct RSCoder {
    private:
        static constinit inline const T zero_{ 0 };
        static constinit inline const T one_{ 1 };
        /// Кодовое расстояние
        static constinit inline const size_t d_{ 2 * t };
        /// Размер поля Галуа
        static constinit inline const size_t gf_size_{ std::numeric_limits<T>::max() };

        /*!
        * \brief Операция сложения над полем Галуа.
        */
        [[nodiscard]] inline const T add(const T op1, const T op2) const noexcept
        { return op1 ^ op2; }

        /*!
        * \brief Операция вычитания над полем Галуа. Операция идентична сложению, введена для удобства.
        */
        [[nodiscard]] inline const T sub(const T op1, const T op2) const noexcept
        { return add(op1, op2); }

        /*!
        * \brief Операция умножения над полем Галуа.
        */
        [[nodiscard]] inline const T mul(const T op1, const T op2) const noexcept {
            if (op1 == 0 || op2 == 0)
                return 0;

            const T v1{ gf_rev(op1) };
            const T v2{ gf_rev(op2) };

            return (v1 >= gf_size_ - v2)
                ? gf(v1 - (v2 ^ gf_size_))
                : gf(v1 + v2);
        }

        /*!
        * \brief Операция деления над полем Галуа.
        */
        [[nodiscard]] inline const T div(const T op1, const T op2) const noexcept {
            if (op1 == zero_ && op2 != zero_)
                return zero_;

            assert(op2 != zero_);

            const T v1{ gf_rev(op1) };
            const T v2{ gf_rev(op2) };

            return (v1 < v2)
                ? gf((v2 - v1) ^ gf_size_)
                : gf(v1 - v2);
        }

        /*!
        * \brief Операция деления полиномов над полем Галуа.
        * return tuple(quotient, remainder)
        */
        [[nodiscard]] inline std::tuple<std::vector<T>, std::vector<T>>
        div(const std::vector<T>& dividend, const std::vector<T>& divisor) const noexcept {
            const size_t divisor_len{ divisor.size() };
            const size_t len{ dividend.size() - divisor_len + 1 };

            assert(len > 0);

            std::vector<T> quotient(len);
            std::vector<T> remainder(divisor_len - 1);

            std::vector<T> tmp_dividend = dividend;
            std::vector<T> tmp_divisor(divisor_len);

            size_t i{ len - 1 };
            std::for_each(tmp_dividend.rbegin(), tmp_dividend.rbegin() + len, [&](T& dividend) {
                quotient[i] = div(dividend, divisor.back());

                for (size_t z = 0; z < divisor_len; ++z)
                    tmp_divisor[z] = mul(divisor[z], quotient[i]);

                for (size_t z = 0; z < divisor_len; ++z)
                    tmp_dividend[i + z] = sub(tmp_dividend[i + z], tmp_divisor[z]);

                i--;
                }
            );

            std::copy(tmp_dividend.begin(), tmp_dividend.begin() + remainder.size(), remainder.begin());

            erase_lead_zeros(remainder);

            return std::make_tuple(quotient, remainder);
        }

        /*!
        * \brief Операция умножения полиномов над полем Галуа.
        */
        [[nodiscard]] inline std::vector<T> mul(const std::vector<T>& multiplicanda,
                                                const std::vector<T>& multiplier) const noexcept {
            const size_t multiplicanda_len{ multiplicanda.size() };
            const size_t multiplier_len{ multiplier.size() };

            if (0 == multiplicanda_len || 0 == multiplier_len)
                return { zero_ };

            const size_t len{ multiplicanda_len + multiplier_len - 1 };
            std::vector<T> product(len);

            for (size_t i = 0; i < multiplicanda_len; ++i)
                for (size_t z = 0; z < multiplier_len; ++z)
                    product[i + z] ^= mul(multiplicanda[i], multiplier[z]);

            return product;
        }

        /*!
        * \brief Операция сложения полиномов над полем Галуа.
        */
        [[nodiscard]] inline const std::vector<T> add(const std::vector<T>& v1,
                                                      const std::vector<T>& v2) const noexcept {
            const size_t min_len{ std::min(v1.size(), v2.size()) };
            const size_t max_len{ std::max(v1.size(), v2.size()) };

            std::vector<T> res(max_len);

            for (size_t i = 0; i < min_len; ++i)
                res[i] = v1[i] ^ v2[i];

            if (v1.size() > v2.size())
                std::copy(v1.begin() + min_len, v1.end(), res.begin() + min_len);

            if (v1.size() < v2.size())
                std::copy(v2.begin() + min_len, v2.end(), res.begin() + min_len);

            return res;
        }

        /*!
        * \brief Операция вычитания полиномов над полем Галуа.
        */
        [[nodiscard]] inline const std::vector<T> sub(const std::vector<T>& v1,
                                                      const std::vector<T>& v2) const noexcept {
            return sum(v1, v2);
        }

        /*!
        * \brief Возвращает неприводимый полином по умолчанию в замисимости от размера кодового слова.
        * \note В качестве оптимизации старший бит полинома не используется, так полином 0x11D, соответствет 0x1D.
        *  ______________________________________________
        * | Размер кодового слова | Неприводимый полином |
        * |----------------------------------------------|
        * |           8           | x^8+x^4+x^3+x^2+1    |
        * |          16           | x^16+x^12+x^3+x+1    |
        * |          32           | x^32+x^22+x^2+x+1    |
        * |          64           | x^64+x^4+x^3+x+1     |
        *  ----------------------------------------------
        */
        [[nodiscard]] static inline const T generate_irreducible_poly() noexcept {
            switch (sizeof(T)) {
            case 1: return static_cast<T>(0x1D);
            case 2: return static_cast<T>(0x100B);
            case 4: return static_cast<T>(0x400007);
            case 8: return static_cast<T>(0x1B);
            }
            return T{};
        }

        /*!
        * \brief Генерирует прямое поле Галуа GF[2 ^ x], где x - размер кодового слова (sizeof(T) * 8).
        * Последний элемент поля всегда равен 1 (GF[std::numeric_limits<T>::max() + 1] == 1),
        * поэтому не учитываем его в массиве, что позволяет использовать размер массива размерности std::numeric_limits<T>::max().
        * \return прямое поле Галуа GF[2 ^ x].
        */
        [[nodiscard]] inline const std::vector<T> generate_gf() const noexcept {
            if constexpr (gf_static) {
                std::vector<T> gf(gf_size_);
                gf[0] = 1;
                gf[1] = 2;

                const T gf_half_size{ static_cast<T>((gf_size_ / 2) + 1) };

                for (T i = 2; i < gf_size_; ++i) {
                    const T p{ gf[i - 1] };
                    gf[i] = p << 1; // если переполнение, то при присваивании последний бит срезается, ...

                    if (p >= gf_half_size)
                        gf[i] ^= irreducible_poly_; //... что позволяет убрать последний бит у неприводимого полинома
                }
                return gf;
            }
            else
                return {};
        }

        /*!
        * Генерирует обратное поле Галуа GF^-1[2 ^ x], где x - размер кодового слова (sizeof(T) * 8).
        * Первый элемент поля всегда соотвествет NAN (GF^-1[0] == NAN),
        * поэтому не учитываем его в массиве, что позволяет использовать размер массива размерности std::numeric_limits<T>::max().
        */
        [[nodiscard]] inline const std::vector<T> generate_gf_rev() const noexcept {
            if constexpr (gf_static) {
                std::vector<T> gf_rev(gf_size_);

                for (size_t i = 0; i < gf_size_; ++i)
                    gf_rev[gf(static_cast<T>(i)) - one_] = static_cast<T>(i);

                return gf_rev;
            }
            else
                return {};
        }

        /*!
        * \brief Динамически вычислет элемент прямого поля Галуа.
        * \param[in] idx индекс в массиве.
        * \param[in] gf_prev элемент прямого поля Галуа вычисленный на предыдущем шаге.
        * Используется для повышения скорости посика текущего элемента.
        */
        [[nodiscard]] inline const T gf_dyn(const T idx, T gf_prev = 0) const noexcept {
            switch (idx) {
            case 0:
            case gf_size_: return 1;
            case 1:		   return 2;
            }

            const size_t gf_half_size{ (gf_size_ / 2) + 1 };

            T gf_curr{ zero_ };

            if (0 == gf_prev) {
                gf_prev = 2;

                for (T i = 2; i < gf_size_; ++i) {
                    gf_curr = gf_prev << 1; // если переполнение, то при присваивании последний бит срезается, ...

                    if (gf_prev >= gf_half_size)
                        gf_curr ^= irreducible_poly_; //... что позволяет убрать последний бит у неприводимого полинома

                    if (idx == i)
                        return gf_curr;

                    gf_prev = gf_curr;
                }
            }
            else {
                gf_curr = gf_prev << 1; // если переполнение, то при присваивании последний бит срезается, ...
                return (gf_prev < gf_half_size)
                    ? gf_curr
                    : gf_curr ^ irreducible_poly_; //... что позволяет убрать последний бит у неприводимого полинома
            }
            return T{};
        }

        /*!
        * \brief Возвращает элемент прямого поля Галуа.
        * \param[in] idx индекс в массиве.
        */
        [[nodiscard]] inline const T gf(const T idx) const noexcept {
            if constexpr (gf_static)
                return gf_size_ == idx ? one_ : gf_[idx];
            else
                return gf_dyn(idx);
        }

        /*!
        * \brief Динамически вычислет элемент обратного поля Галуа.
        * \param[in] idx индекс в массиве.
        */
        [[nodiscard]] inline const T gf_rev_dyn(const T idx) const noexcept {
            T gf{ zero_ };
            for (T i = 0; i < gf_size_; ++i) {
                gf = gf_dyn(i, gf);
                if (idx == gf)
                    return i;
            }
            return T{};
        }

        /*!
        * \brief Возвращает элемент обратного поля Галуа.
        * \param[in] idx индекс в массиве.
        */
        [[nodiscard]] inline const T gf_rev(const T idx) const noexcept {
            if constexpr (gf_static) {
                assert(idx != 0);
                return idx == 0
                    ? 0
                    : gf_rev_[idx - one_];
            }
            else
                return gf_rev_dyn(idx);
        }

        /*!
        * \brief Возвращает порождающий полином степени d=2t.
        */
        [[nodiscard]] inline const std::vector<T> generate_generic() const noexcept {
            std::vector<T> generic{ gf(1), 1 };
            for (T i = 1; i < d_; ++i)
                generic = mul(generic, std::vector<T>{ gf(i + 1), one_ });
            return generic;
        }

        /// Неприводимы полином.
        static inline const T irreducible_poly_{ generate_irreducible_poly() };

        /// Прямое поле Галуа (используется только при агрументе шаблона класса is_static == true).
        const std::vector<T> gf_;
        /// Обратное поле Галуа (используется только при агрументе шаблона класса is_static == true).
        const std::vector<T> gf_rev_;

        /*!
        * \brief Стирает нулевые коэффициенты в старших индексах.
        * \param[in, out] poly полином.
        */
        inline void erase_lead_zeros(std::vector<T>& poly) const noexcept {
            const size_t poly_len{ poly.size() };
            for (size_t i = 0; i < poly_len; ++i) {
                if (0 != poly[poly_len - 1 - i]) {
                    if (0 < i && i < poly_len)
                        poly.erase(poly.end() - i, poly.end());
                    break;
                }
            }
        }

        /*!
        * \brief Приводит полином к числу.
        */
        [[nodiscard]] inline const T eval_by_pow(const std::vector<T>& v, const T pow) const noexcept {
            T sum{ v[0] };

            const size_t v_size{ v.size() };
            for (size_t i = 1; i < v_size; ++i) {
                size_t idx = pow * i;
                sum ^= mul(v[i], gf(static_cast<T>(idx > gf_size_
                    ? idx % gf_size_
                    : idx)));
            }

            return sum;
        }

        /*!
        * \brief Приводит полином к числу.
        */
        [[nodiscard]] inline const T eval(const std::vector<T>& v, const T x) const noexcept
        { return eval_by_pow(v, gf_rev(x)); }

        /*!
        * \brief Оператор интертирования.
        */
        [[nodiscard]] inline const T invert(const T val) const noexcept
        { return static_cast<T>(gf_size_ - val); }

        /*!
        * \brief Вычисялет полином синдромов.
        * \param[in] encoded_data данные в коде РС.
        * \result полином синдромов.
        */
        [[nodiscard]] inline const std::vector<T> syndromes(const std::vector<T>& encoded_data) const noexcept {
            std::vector<T> syndromes(d_ + 1);
            syndromes[0] = 1;

            for (size_t i = 0; i < d_; ++i)
                syndromes[i + 1] = eval_by_pow(encoded_data, static_cast<T>(i + 1));

            erase_lead_zeros(syndromes);

            return syndromes;
        }

        /*!
        * \brief Вычисляет полином локаторов ошибок и полином велечин ошибок. Используется
        * алгоритм Евклида нахождения НОД двух челых чисел.
        * \param[in] syndromes полином синдромов.
        * \result tuple(sigma, lambda), где
        * sigma полином локаторов ошибок,
        * lambda полином велечин ошибок.
        */
        [[nodiscard]] inline const std::tuple<const std::vector<T>, const std::vector<T>>
        euclid_gcd(const std::vector<T>& syndromes) const noexcept {
            const size_t len{ d_ + 1 };
            std::vector<T> r0(len + 1);
            r0[len] = 1;
            auto r1{ syndromes };
            std::vector<T> b0{ zero_ };
            std::vector<T> b1{ one_ };

            while (true) {
                const auto& [q, r2] { div(r0, r1) };
                std::vector<T> b2{ add(b0, mul(q, b1)) };

                if (r2.size() - 1 <= t) {
                    const auto& [sigma, _] { div(b2, { b2[0] }) };
                    const auto& [lambda, __] { div(r2, { r2[0] }) };

                    return { sigma, lambda };
                }
                else {
                    r0 = std::move(r1);
                    r1 = r2;
                    b0 = std::move(b1);
                    b1 = std::move(b2);
                }
            }
        }

        /*!
        * \brief Осуществляет поиск по методу Ченя: поиск позиций ошибок кодовых слов путем
        * перебора элементов поля Галуа.
        * \param[in] sigma полином локаторов ошибок
        * \result список позиций ошибок кодовых слов в порядке возрастания.
        */
        [[nodiscard]] inline const std::vector<T> chen_search(const std::vector<T>& sigma) const noexcept {
            std::vector<T> locator;

            for (size_t i = gf_size_; i > 0; --i) {
                if (0 == eval_by_pow(sigma, static_cast<T>(i)))
                    locator.push_back(invert(static_cast<T>(i)));
            }

            return locator;
        }

        /*!
        * \brief Вычисляет формальную производную полинома.
        * \param[in] poly полином
        * \result формальная производная.
        */
        [[nodiscard]] inline const std::vector<T> formal_derivative(const std::vector<T>& poly) const noexcept {
            const size_t size{ poly.size() - (poly.size() % 2 == 0 ? 1 : 2) };
            std::vector<T> derivative(size);
            for (size_t i = 0; i < size; ++i) {
                derivative[i] = ((i + 1) % 2 == 0)
                    ? 0
                    : poly[i + 1];
            }
            return derivative;
        }

        /*!
        * \brief Вычисляет полином ошибок по алгоритму Форни.
        * \param[in] tuple(sigma, lambda), где sigma - полином локаторов ошибок; lambda - полином велечин ошибок.
        * \param[in] locator список позиций ошибок кодовых слов.
        * \result полином ошибок.
        */
        [[nodiscard]] inline const std::vector<T> algorithm_forni(const std::tuple<const std::vector<T>,
                                                                  const std::vector<T>>&euclid_gcd,
                                                                  const std::vector<T>& locator) const noexcept {
            const size_t locator_len{ locator.size() };
            std::vector<T> error(static_cast<T>(locator[locator_len - 1] + 1), zero_);

            const auto& [sigma, lambda] { euclid_gcd };
            const auto& derivative_sigma{ formal_derivative(sigma) };

            for (size_t i = 0; i < locator_len; ++i) {
                const auto& w{ mul(gf(locator[i]), eval_by_pow(lambda, invert(locator[i]))) };
                const auto& l{ eval_by_pow(derivative_sigma, invert(locator[i])) };

                error[locator[i]] = div(w, l);
            }

            return error;
        }

        /*!
        * \brief Вычисляет полином ошибок по алгоритму Форни.
        * \param[in] sigma полином локаторов ошибок.
        * \result locator список позиций ошибок кодовых слов.
        * \result полином ошибок.
        */
        [[nodiscard]] inline const std::vector<T> algorithm_forni(const std::vector<T>& syndromes,
                                                                  const std::vector<T> sigma,
                                                                  const std::vector<T>& locator) const noexcept {
            std::vector<T> error(locator[locator.size() - 1] + 1, 0);

            const auto& derivative_sigma{ formal_derivative(sigma) };

            std::vector<T> s(syndromes.size() - 1);
            std::copy(syndromes.begin() + 1, syndromes.end(), s.begin());
            std::vector<T> omega{ mul(s, sigma) };

            if (omega.size() > d_)
                std::fill(omega.begin() + d_, omega.end(), 0);

            erase_lead_zeros(omega);

            const size_t locator_size{ locator.size() };
            for (size_t i = 0; i < locator_size; ++i) {
                const auto& w{ eval_by_pow(omega, invert(locator[i])) };
                const auto& l{ eval_by_pow(derivative_sigma, invert(locator[i])) };

                error[locator[i]] = div(w, l);
            }

            return error;
        }

    public:
        inline RSCoder()
            : gf_{ generate_gf() }
            , gf_rev_{ generate_gf_rev() } {
        }

        /*!
        * \brief Выполяет кодирование кодом РС.
        * \param[in] data исходные данные.
        * \return кодированные данные.
        */
        template<class _It>
        [[nodiscard]] inline std::vector<T> encode(const _It begin, const _It end,
                                                   std::vector<T>* const out_remainder = nullptr) const {
            std::vector<T> encoded(std::distance(begin, end) + d_);
            std::copy_backward(begin, end, encoded.end());
            const std::vector<T> generic_{ generate_generic() };
            const auto& [_, remainder] = div(encoded, generic_);
            std::copy(remainder.begin(), remainder.end(), encoded.begin());
            if (out_remainder) *out_remainder = remainder;
            return encoded;
        }

        /*!
        * \brief Выполяет кодирование кодом РС.
        * \param[in] data указатель на массив исходных данные.
        * \param[in] size размер массива исходных данных.
        * \return кодированные данные.
        */
        [[nodiscard]] inline std::vector<T> encode(const T* const data, const size_t size,
                                                   std::vector<T>* const out_remainder = nullptr) const noexcept {
            return encode(data, data + size, out_remainder);
        }

        /*!
        * \brief Выполяет декодирование. В каческе алгоритма вычисления полинома локаторов ошибок используется
        * алгоритм Евклида нахождения НОД двух челых чисел.
        * \param[in] encoded_data данные кодированне кодом РС.
        * \param[out] out_syndromes полином синдрома ошибок. При полиноме равном 1 ошибки в кодированном
        * сообщении не обнаружены.
        * \param[out] out_sigma полином локаторов ошибок.
        * \param[out] out_lambda полином величин ошибок.
        * \param[out] out_locator локатор позиций ошибок.
        * \param[out] out_error полином ошибок.
        * \return декодированные данные.
        */
        [[nodiscard]] inline std::vector<T> decode(const std::vector<T>& encoded,
                                                   std::vector<T>* const out_syndromes = nullptr,
                                                   std::vector<T>* const out_sigma = nullptr,
                                                   std::vector<T>* const out_lambda = nullptr,
                                                   std::vector<T>* const out_locator = nullptr,
                                                   std::vector<T>* const out_error = nullptr) const noexcept {
            const auto& syndromes{ RSCoder::syndromes(encoded) };

            if (1 == syndromes.size()) {
                if (out_syndromes) *out_syndromes = syndromes;
                return encoded;
            }

            const auto& [sigma, lambda] { euclid_gcd(syndromes) };
            const auto& locator{ chen_search(sigma) };
            const auto& error{ algorithm_forni({ sigma, lambda }, locator) };
            std::vector<T> decoded{ add(error, encoded) };

            decoded.erase(decoded.begin(), decoded.begin() + d_);
            decoded.shrink_to_fit();

            if (out_syndromes) *out_syndromes = syndromes;
            if (out_sigma) *out_sigma = sigma;
            if (out_lambda) *out_lambda = lambda;
            if (out_locator) *out_locator = locator;
            if (out_error) *out_error = error;

            return decoded;
        }
    };
}
