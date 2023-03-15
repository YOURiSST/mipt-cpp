#pragma once

#include <vector>
#include <sstream>
#include <cmath>

class BigInteger {
private:
    std::vector<int64_t> digits_;
    bool is_positive_;

    static size_t GetSizeFromString(size_t string_size) {
        if (string_size == 0 || string_size == 1) {
            return 1;
        }
        return (string_size - 1) / GetBaseLen() + 1;
    }

    static size_t GetLen(int64_t in_int) {
        if (in_int == 0) {
            return 1;
        }
        size_t len = 0;
        while (in_int > 0) {
            len++;
            in_int /= 10;
        }
        return len;
    }

    static const int64_t BASE_ = 10'000'000;
    static const int64_t BASE_len_ = 7;

    void KillLeadingZeros() {
        while (GetSize() > 1 && digits_.back() == 0) {
            digits_.pop_back();
        }
        if (IsNull()) {
            is_positive_ = true;
        }
    }

    bool ComparationRun(const BigInteger& other) const;


public:
    BigInteger();

    BigInteger(int64_t x);

    explicit BigInteger(const std::string& s);

    BigInteger(const BigInteger& other) : digits_(other.digits_), is_positive_(other.is_positive_) {}

    void BuildByString(const std::string& s);

    std::string toString() const;

    bool IsPositive() const;

    static int64_t GetBase();

    static size_t GetBaseLen();

    size_t GetSize() const;

    bool IsNull() const;

    explicit operator bool() const;

    BigInteger& IdentifyModuleOperation(const BigInteger& other, bool sum);

    void MultByBase();

    void ModuleSum(const BigInteger& other);

    void ModuleDif(const BigInteger& other, bool our_greater);

    BigInteger operator-() const;

    friend bool operator==(const BigInteger& our, const BigInteger& other);

    friend bool operator<(const BigInteger& our, const BigInteger& other);

    BigInteger& operator=(const BigInteger& other) = default;

    BigInteger& operator+=(const BigInteger& other);

    BigInteger& operator-=(const BigInteger& other);

    BigInteger& operator++();

    BigInteger& operator--();

    BigInteger operator++(int);

    BigInteger operator--(int);

    BigInteger& operator*=(const BigInteger& other);

    BigInteger& operator/=(const BigInteger& other);

    BigInteger& operator%=(const BigInteger& other);

    int64_t Front() const;
};

bool BigInteger::IsPositive() const {
    return is_positive_;
}

int64_t BigInteger::GetBase() {
    return BigInteger::BASE_;
}

size_t BigInteger::GetBaseLen() {
    return BigInteger::BASE_len_;
}

size_t BigInteger::GetSize() const {
    return digits_.size();
}

bool BigInteger::IsNull() const {
    return digits_[0] == 0 && GetSize() == 1ull;
}

BigInteger::operator bool() const {
    return !IsNull();
}

int64_t BigInteger::Front() const {
    return digits_[0];
}

std::string BigInteger::toString() const {
    std::stringstream ss;
    if (!is_positive_) {
        ss << "-";
    }
    ss << digits_.back();
    for (int i = static_cast<int>(GetSize()) - 2; i > -1; --i) {
        size_t len = GetLen(digits_[i]);
        for (size_t j = 0; j < GetBaseLen() - len; ++j) {
            ss << "0";
        }
        ss << digits_[i];
    }

    return ss.str();
}

void BigInteger::BuildByString(const std::string& s) {
    digits_.clear();
    if (s == "0" || s == "-0" || s.empty()) {
        digits_.emplace_back(0);
        is_positive_ = true;
        return;
    }
    is_positive_ = true;
    size_t s_len = s.size();
    if (s.front() == '-') {
        is_positive_ = false;
        s_len--;
    }
    digits_.resize(GetSizeFromString(s_len));
    for (size_t i = 0; i < s_len; ++i) {
        digits_[i / GetBaseLen()] += (s[s.size() - i - 1] - '0') *
                                     static_cast<int64_t>(powl(10, i % GetBaseLen()));
    }
    KillLeadingZeros();
}

// constructors

BigInteger::BigInteger() {
    digits_.emplace_back(0);
    is_positive_ = true;
}

BigInteger::BigInteger(int64_t x) {
    is_positive_ = true;
    if (x == 0) {
        digits_.emplace_back(0);
    }
    if (x < 0) {
        is_positive_ = false;
        x *= -1;
    }
    while (x) {
        digits_.emplace_back(x % GetBase());
        x /= GetBase();
    }
}

BigInteger::BigInteger(const std::string& s) {
    BuildByString(s);
}

bool operator==(const BigInteger& our, const BigInteger& other) {
    if (our.is_positive_ != other.is_positive_) {
        return false;
    }
    if (other.GetSize() != our.GetSize()) {
        return false;
    }
    for (size_t i = 0; i < our.GetSize(); ++i) {
        if (our.digits_[i] != other.digits_[i]) {
            return false;
        }
    }
    return true;
}

// operators

bool BigInteger::ComparationRun(const BigInteger& other) const {
    if (GetSize() < other.GetSize()) {
        return true;
    }
    if (GetSize() > other.GetSize()) {
        return false;
    }
    for (int i = static_cast<int>(GetSize()) - 1; i > -1; --i) {
        if (digits_[i] == other.digits_[i]) {
            continue;
        }
        if (digits_[i] > other.digits_[i]) {
            return false;
        }
        return true;
    }
    if (is_positive_) {
        return false;
    }
    return true;
}

bool operator!=(const BigInteger& our, const BigInteger& other) {
    return !(our == other);
}

bool operator<(const BigInteger& our, const BigInteger& other) {
    if (our.is_positive_ && !other.is_positive_) {
        return false;
    }
    if (!our.is_positive_ && other.is_positive_) {
        return true;
    }
    if (our.is_positive_) {
        return our.ComparationRun(other);
    } else {
        return !our.ComparationRun(other);
    }
}

bool operator>(const BigInteger& our, const BigInteger& other) {
    return (other < our);
}

bool operator<=(const BigInteger& our, const BigInteger& other) {
    return !(our > other);
}

bool operator>=(const BigInteger& our, const BigInteger& other) {
    return !(our < other);
}

BigInteger operator "" _bi(const char* c_string, size_t) {
    BigInteger res(c_string);
    return res;
}

BigInteger operator "" _bi(unsigned long long int_input) {
    BigInteger res(int_input);
    return res;
}


void BigInteger::ModuleSum(const BigInteger& other) {
    size_t result_len = std::max(GetSize(), other.GetSize()) + 1;
    digits_.resize(result_len);
    for (size_t i = 0; i < other.GetSize(); ++i) {
        digits_[i] += other.digits_[i];
        if (digits_[i] >= GetBase()) {
            digits_[i] -= GetBase();
            ++digits_[i + 1];
        }
    }
    KillLeadingZeros();
}

void BigInteger::ModuleDif(const BigInteger& other, bool our_greater) {
    size_t result_len = std::max(GetSize(), other.GetSize()) + 1;
    digits_.resize(result_len);
    for (size_t i = 0; i < other.GetSize(); ++i) {
        if (our_greater) {
            digits_[i] -= other.digits_[i];
        } else {
            digits_[i] = other.digits_[i] - digits_[i];
        }
    }
    for (int i = 0; i < static_cast<int>(result_len) - 1; ++i) {
        if (digits_[i] < 0) {
            digits_[i] += GetBase();
            --digits_[i + 1];
        }
    }
    KillLeadingZeros();
}

BigInteger& BigInteger::IdentifyModuleOperation(const BigInteger& other, bool sum) {
    if ((is_positive_ == other.is_positive_ && sum) ||
        (is_positive_ != other.is_positive_ && !sum)) {
        ModuleSum(other);
        return *this;
    }
    bool our_greater = !ComparationRun(other);
    ModuleDif(other, our_greater);
    if (!our_greater && !IsNull()) {
        is_positive_ = !is_positive_;
    }
    return *this;
}

BigInteger& BigInteger::operator+=(const BigInteger& other) {
    return IdentifyModuleOperation(other, true);
}

BigInteger& BigInteger::operator-=(const BigInteger& other) {
    return IdentifyModuleOperation(other, false);
}

BigInteger BigInteger::operator-() const {
    BigInteger temp(*this);
    if (IsNull()) {
        return temp;

    }
    temp.is_positive_ = !temp.is_positive_;
    return temp;
}

BigInteger operator-(const BigInteger& left, const BigInteger& right) {
    BigInteger res(left);
    return res -= right;
}

BigInteger operator+(const BigInteger& left, const BigInteger& right) {
    BigInteger res(left);
    return res += right;
}

BigInteger& BigInteger::operator++() {
    *this += 1;
    return *this;
}


BigInteger& BigInteger::operator--() {
    return *this -= 1;
}

BigInteger BigInteger::operator++(int) {
    *this += 1;
    return (*this - 1);
}


BigInteger BigInteger::operator--(int) {
    *this -= 1;
    return (*this + 1);
}

BigInteger& BigInteger::operator*=(const BigInteger& other) {
    if (IsNull()) {
        return *this;
    }
    if (other.IsNull()) {
        digits_.clear();
        digits_.emplace_back(0);
        is_positive_ = true;
        return *this;
    }
    if (other == 1) {
        return *this;
    }
    if (other == -1) {
        is_positive_ = !is_positive_;
        KillLeadingZeros();
        return *this;
    }
    std::vector<int64_t> result(GetSize() + other.GetSize() + 4);
    for (size_t i = 0; i < GetSize(); ++i) {
        for (size_t j = 0; j < other.GetSize(); ++j) {
            result[i + j] += digits_[i] * other.digits_[j];
        }
    }
    for (size_t i = 0; i < result.size() - 1; ++i) {
        if (result[i] >= GetBase()) {
            result[i + 1] += result[i] / GetBase();
            result[i] %= GetBase();
        }
    }
    digits_ = result;
    if (is_positive_ == other.is_positive_) {
        is_positive_ = true;
    } else {
        is_positive_ = false;
    }
    KillLeadingZeros();
    return *this;
}

BigInteger operator*(const BigInteger& left, const BigInteger& right) {
    BigInteger res(left);
    return res *= right;
}

void BigInteger::MultByBase() {
    digits_.insert(digits_.begin(), 0);
}

BigInteger& BigInteger::operator/=(const BigInteger& other) {
    if (IsNull()) {
        return *this;
    }

    BigInteger result;
    std::vector<int64_t> result_digits_(GetSize() + 1);
    BigInteger current(0);
    BigInteger other_mod = other;
    other_mod.is_positive_ = true;
    for (int i = static_cast<int>(GetSize()) - 1; i > -1; --i) {
        current.MultByBase();
        current.digits_[0] = digits_[i];
        current.KillLeadingZeros();
        int64_t res_x = 0, left = 0, right = GetBase();
        while (left < right + 1) {
            int64_t mid = (left + right) >> 1;
            BigInteger now = other_mod * mid;
            if (now <= current) {
                res_x = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        digits_[i] = res_x;
        current -= other_mod * res_x;
    }
    is_positive_ = (is_positive_ == other.is_positive_);
    KillLeadingZeros();
    return *this;
}

BigInteger operator/(const BigInteger& left, const BigInteger& right) {
    BigInteger res(left);
    return res /= right;
}

BigInteger& BigInteger::operator%=(const BigInteger& other) {
    return *this -= (*this / other * other);
}

BigInteger operator%(const BigInteger& left, const BigInteger& right) {
    BigInteger res(left);
    return res %= right;
}

std::ostream& operator<<(std::ostream& os, const BigInteger& our) {
    os << our.toString();
    return os;
}

std::istream& operator>>(std::istream& is, BigInteger& our) {
    std::string s;
    is >> s;
    our.BuildByString(s);
    return is;
}


/*/////////////////////////////////RATIONAL\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/


class Rational {
private:
    static const size_t default_precision = 18;
    BigInteger numerator_;
    BigInteger denumerator_;

    static BigInteger GCD(BigInteger& left, BigInteger& right);

    void SignControl();

    void Normilize();

    bool IsNull() const {
        return numerator_.IsNull();
    }

public:
    Rational() : numerator_(0), denumerator_(1) {}

    Rational(const BigInteger& numerator, const BigInteger& denumerator) :
        numerator_(numerator), denumerator_(denumerator) {
        Normilize();
    }

    Rational(const BigInteger& numerator) : numerator_(numerator), denumerator_(1) {}

    Rational(int numerator) : numerator_(numerator), denumerator_(1) {}

    Rational operator-() const;

    friend bool operator==(const Rational& our, const Rational& other);

    friend bool operator<(const Rational& our, const Rational& other);

    Rational& operator+=(const Rational& other);

    Rational& operator-=(const Rational& other);

    Rational& operator*=(const Rational& other);

    Rational& operator/=(const Rational& other);

    std::string toString() const;

    std::string asDecimal(size_t precision) const;

    explicit operator double() const;
};

BigInteger Rational::GCD(BigInteger& left, BigInteger& right) {
    if (left == 0 || right == 0) {
        return 1;
    }
    if (left == 1 || right == 1) {
        return 1;
    }


    while (right != 0) {
        left %= right;
        std::swap(left, right);
    }
    return right + left;

}

void Rational::SignControl() {
    if ((numerator_.IsPositive() != denumerator_.IsPositive() && numerator_.IsPositive()) ||
        (numerator_.IsPositive() == denumerator_.IsPositive() && !numerator_.IsPositive())) {
        numerator_ *= -1;
        denumerator_ *= -1;
    }
}

void Rational::Normilize() {
    if (IsNull()) {
        denumerator_ = 1;
        return;
    }
    if (numerator_ == denumerator_) {
        numerator_ = denumerator_ = 1;
        return;
    }
    if (numerator_ == 1 || numerator_ == -1) {
        SignControl();
        return;
    }
    if (numerator_ >= denumerator_ && (numerator_ % denumerator_ == 0)) {
        numerator_ /= denumerator_;
        denumerator_ = 1;
        SignControl();
        return;
    }
    if (numerator_ < denumerator_ && (denumerator_ % numerator_ == 0)) {
        denumerator_ /= numerator_;
        numerator_ = 1;
        SignControl();
        return;
    }
    BigInteger numerator_abs = numerator_, denumerator_abs = denumerator_;
    if (!numerator_abs.IsPositive()) {
        numerator_abs *= -1;
    }
    if (!denumerator_abs.IsPositive()) {
        denumerator_abs *= -1;
    }
    BigInteger cur_gcd = GCD(numerator_abs, denumerator_abs);

    numerator_ /= cur_gcd;
    denumerator_ /= cur_gcd;
    SignControl();
}

std::string Rational::toString() const {
    if (denumerator_ == 1) {
        return numerator_.toString();
    }
    return numerator_.toString() + '/' + denumerator_.toString();
}

bool operator<(const Rational& our, const Rational& other){
    return (our.numerator_ * other.denumerator_ < our.denumerator_ * other.numerator_);
}

bool operator==(const Rational& our, const Rational& other) {
    return (our.numerator_ * other.denumerator_ == our.denumerator_ * other.numerator_);
}

bool operator!=(const Rational& left, const Rational& right) {
    return !(left == right);
}

bool operator>(const Rational& left, const Rational& right) {
    return right < left;
}

bool operator<=(const Rational& left, const Rational& right) {
    return !(left > right);
}

bool operator>=(const Rational& left, const Rational& right) {
    return !(left < right);
}

Rational& Rational::operator*=(const Rational& other) {
    if (IsNull()) {
        return *this;
    }
    numerator_ *= other.numerator_;
    denumerator_ *= other.denumerator_;
    Normilize();
    return *this;
}

Rational& Rational::operator+=(const Rational& other) {
    if (denumerator_ == other.denumerator_) {
        numerator_ += other.numerator_;
        Normilize();
        return *this;
    }
    numerator_ *= other.denumerator_;
    numerator_ += (other.numerator_ * denumerator_);
    denumerator_ *= other.denumerator_;
    Normilize();
    return *this;
}

Rational& Rational::operator-=(const Rational& other) {
    if (denumerator_ == other.denumerator_) {
        numerator_ -= other.numerator_;

        Normilize();
        return *this;
    }
    numerator_ *= other.denumerator_;
    numerator_ -= (other.numerator_ * denumerator_);
    denumerator_ *= other.denumerator_;
    Normilize();
    return *this;
}

Rational& Rational::operator/=(const Rational& other) {
    if (*this == 0) {
        return *this;
    }
    numerator_ *= other.denumerator_;
    denumerator_ *= other.numerator_;
    Normilize();
    return *this;
}

Rational operator+(const Rational& left, const Rational& right) {
    Rational res(left);
    return res += right;
}

Rational operator-(const Rational& left, const Rational& right) {
    Rational res(left);
    return res -= right;
}

Rational operator*(const Rational& left, const Rational& right) {
    Rational res(left);
    return res *= right;
}

Rational operator/(const Rational& left, const Rational& right) {
    Rational res(left);
    return res /= right;
}

Rational Rational::operator-() const {
    Rational res(*this);
    return res *= -1;
}

std::ostream& operator<<(std::ostream& os, const Rational& our) {
    os << our.toString();
    return os;
}

std::string Rational::asDecimal(size_t precision = 0) const {
    if (precision == 0) {
        return (numerator_ / denumerator_).toString();
    }
    BigInteger numerator_copy = numerator_;
    if (!numerator_copy.IsPositive()) {
        numerator_copy *= -1;
    }
    BigInteger integer_part = (numerator_copy / denumerator_);
    numerator_copy %= denumerator_;
    std::vector<int64_t> temp_digits;
    for (size_t i = 0; i < precision + 1; ++i) {
        if (numerator_copy) {
            numerator_copy *= 10;
            for (; numerator_copy < denumerator_;) {
                numerator_copy *= 10;
                temp_digits.emplace_back(0);
            }
            temp_digits.emplace_back((numerator_copy / denumerator_).Front());
            numerator_copy %= denumerator_;
        } else {
            temp_digits.emplace_back(0);
        }
    }

    if (temp_digits[precision] > 4) {
        ++temp_digits[precision - 1];
        {
            int i = precision - 1;
            while (temp_digits[i] >= 10 && i > 0) {
                temp_digits[i] -= 10;
                ++temp_digits[i - 1];
                --i;
            }
        }

        if (temp_digits[0] >= 10) {
            temp_digits[0] -= 10;
            ++integer_part;
        }
    }
    std::string ans;
    if (!numerator_.IsPositive()) {
        ans += "-";
    }
    ans += integer_part.toString();
    ans += '.';
    for (size_t i = 0; i < precision; ++i) {
        ans += std::to_string(temp_digits[i]);
    }
    return ans;
}

Rational::operator double() const {
    std::stringstream ss;
    ss << asDecimal(Rational::default_precision);
    double to_ret;
    ss >> to_ret;
    return to_ret;
}
