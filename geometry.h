#pragma once

#include <array>
#include <iostream>
#include <cmath>
#include <vector>

class Line;

struct Point {
    constexpr static const double EPS = 1e-7;
    constexpr static const double PI = M_PI;

    double x = 0, y = 0;

    Point() = default;
    Point(double x, double y) : x(x), y(y) {}
    Point& operator+=(const Point& other);
    Point& operator-=(const Point& other);
    Point& operator*=(double coeff);
    Point& operator/=(double coeff);

    void radianRotation(double ang) {
        double tmp_x_ = x * cos(ang) - y * sin(ang);
        double tmp_y_ = x * sin(ang) + y * cos(ang);
        x = tmp_x_;
        y = tmp_y_;
    }
    void degreeRotation(double ang) {
        ang = (ang * Point::PI / 180.0);
        radianRotation(ang);
    }
    void reflectByPoint(const Point& p);
    void reflectByAxis(const Line& a);
    bool isOnVec(const Point& a, const Point& b) const;
    double abs() const;
};

bool isEqual(double a, double b) {
    return std::abs(a - b) < Point::EPS;
}

bool operator==(const Point& left, const Point& right) {
    return isEqual(left.x, right.x) && isEqual(left.y, right.y);
}

bool operator!=(const Point& left, const Point& right) {
    return !(left == right);
}

Point& Point::operator+=(const Point& other) {
    x += other.x;
    y += other.y;
    return *this;
}

Point& Point::operator-=(const Point& other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

Point operator+(const Point& left, const Point& right) {
    Point res = left;
    return res += right;
}

Point operator-(const Point& left, const Point& right) {
    Point res = left;
    return res -= right;
}

Point& Point::operator*=(double coeff) {
    x *= coeff;
    y *= coeff;
    return *this;
}

Point& Point::operator/=(double coeff) {
    x /= coeff;
    y /= coeff;
    return *this;
}

Point operator*(const Point& left, double right) {
    Point res = left;
    return res *= right;
}

Point operator/(const Point& left, double right) {
    Point res = left;
    return res /= right;
}

void Point::reflectByPoint(const Point& p) {
    *this = (p * 2.0 - *this);
}

double Point::abs() const {
    return sqrt(x * x + y * y);
}

double dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y;
}

double cross(const Point& a, const Point& b) {
    return a.x * b.y - b.x * a.y;
}

double vecLen(const Point& ft, const Point& sc) {
    return (ft - sc).abs();
}

double getCos(const Point& a, const Point& b) {
    return dot(a, b) / (a.abs() * b.abs());
}

double getSin(const Point& a, const Point& b) {
    return cross(a, b) / (a.abs() * b.abs());
}

bool Point::isOnVec(const Point& a, const Point& b) const {
    return x + Point::EPS >= std::min(a.x, b.x) && x - Point::EPS <= std::max(a.x, b.x) &&
           y + Point::EPS >= std::min(a.y, b.y) && y - Point::EPS <= std::max(a.y, b.y);
}

double getzCrossProduct(const Point& ft, const Point& sc, const Point& trd) {
    double dx1 = sc.x - ft.x;
    double dy1 = sc.y - ft.y;
    double dx2 = trd.x - sc.x;
    double dy2 = trd.y - sc.y;
    return dx1 * dy2 - dy1 * dx2;
}

class Line {
private:
public:
    double a_, b_, c_;
    Point getNormal() const {
        return Point{a_, b_};
    }
    double getC() const {
        return c_;
    }

    double getB() const {
        return b_;
    }

    double getA() const {
        return a_;
    }

    Line(double a, double b, double c) : a_(a), b_(b), c_(c) {}
    Line(const Point& ft, const Point& sc) {
        a_ = ft.y - sc.y;
        b_ = sc.x - ft.x;
        c_ = ft.x * sc.y - sc.x * ft.y;
    }

    Line(double k, double b) : a_(-k), b_(1), c_(-b) {}
    Line(const Point& p, double k) {
        a_ = -k;
        b_ = 1;
        c_ = p.x * (p.y + k) - (p.x + 1) * p.y;
    }

    double insertPoint(const Point& p) {
        return p.x * a_ + p.y * b_ + c_;
    }

    Point crossPoint(const Line& other) {
        Point a = {a_, other.a_};
        Point b = {b_, other.b_};
        Point c = {-c_, -other.c_};
        return {cross(c, b) / cross(a, b), cross(a, c) / cross(a, b)};
    }
    bool containsPoint(const Point& p) const;
};


bool operator==(const Line& our, const Line& other) {
    std::vector<double> ftLine = {our.a_, our.b_, our.c_};
    std::vector<double> scLine = {other.a_, other.b_, other.c_};
    int notNull = 4;
    for (int i = 0; i < 3; ++i) {
        if (!isEqual(ftLine[i], 0.0)) {
            notNull = i;
            break;
        }
    }
    if (notNull != 4) {
        if (isEqual(ftLine[notNull], 0.0)) {
            return false;
        }
        for (int i = 0; i < 3; ++i) {
            ftLine[i] /= ftLine[notNull];
            scLine[i] /= scLine[notNull];
        }
    }
    for (int i = 0; i < 3; ++i) {
        if (!isEqual(ftLine[i], scLine[i])) {
            return false;
        }
    }
    return true;
}
bool operator!=(const Line& our, const Line& other) {
    return !(our == other);
}


bool Line::containsPoint(const Point& p) const {
    return isEqual((a_ * p.x + b_ * p.y + c_), 0.0);
}

void Point::reflectByAxis(const Line& a) {
    Point normalVec = a.getNormal();
    double dist = std::abs(normalVec.x * x + normalVec.y * y + a.getC());
    normalVec /= (normalVec.abs() * normalVec.abs());

    normalVec *= dist;
    *this += normalVec;
    if (isEqual(a.getNormal().x * x + a.getNormal().y * y + a.getC(), 0.0)) {
        *this += normalVec;
    } else {
        *this -= normalVec * 3;
    }
}

class Shape {
public:
    virtual bool operator==(const Shape& other) const = 0;
    virtual bool operator!=(const Shape& other) const = 0;
    virtual bool isCongruentTo(const Shape& other) const = 0;
    virtual bool isSimilarTo(const Shape& other) const = 0;
    virtual bool containsPoint(const Point& point) const = 0;

    virtual Shape& rotate(const Point& center, double angle) = 0;
    virtual Shape& reflect(const Point& center) = 0;
    virtual Shape& reflect(const Line& axis) = 0;
    virtual Shape& scale(const Point& center, double coefficient) = 0;

    virtual double perimeter() const = 0;
    virtual double area() const = 0;

    virtual ~Shape() = default;
};

class Polygon : public Shape {
private:
    std::array<Point, 6> congruentHelper(size_t ft_i, size_t sc_i_pos, size_t sc_i_neg,
                                       const std::vector<Point>& second_vertices_) const;
protected:
    std::vector<Point> vertices_;
public:
    explicit Polygon(const std::vector<Point>& vertices) : vertices_(vertices) {}

    template<class... Args>
    explicit Polygon(Args&& ... vertices) : vertices_{std::forward<Args>(vertices)...} {}

    const std::vector<Point>& getVertices() const {
        return vertices_;
    }
    size_t verticesCount() const {
        return vertices_.size();
    }
    bool isConvex() const;
    double perimeter() const override;
    double area() const override;
    bool isEqualVertices(const std::vector<Point>& second_vertices_) const;

    bool operator==(const Shape& other) const override;
    bool operator==(const Polygon& other) const;
    bool operator!=(const Shape& other) const override;
    bool operator!=(const Polygon& other) const;

    bool isSimilarTo(const Shape& other) const override;
    bool isCongruentTo(const Shape& other) const override;
    bool containsPoint(const Point& point) const override;

    Polygon& rotate(const Point& center, double angle) override;
    Polygon& reflect(const Point& center) override;
    Polygon& reflect(const Line& axis) override;
    Polygon& scale(const Point& center, double coefficient) override;
};

bool Polygon::isConvex() const {
    bool first_sign = (getzCrossProduct(vertices_[0], vertices_[1], vertices_[2]) + Point::EPS) > 0;
    for (size_t i = 1; i < vertices_.size(); ++i) {
        if ((getzCrossProduct(vertices_[i], vertices_[(i + 1) % vertices_.size()],
                              vertices_[(i + 2) % vertices_.size()]) + Point::EPS > 0) != first_sign) {
            return false;
        }
    }
    return true;
}

double Polygon::perimeter() const {
    double ans = 0.0;
    for (size_t i = 0; i < vertices_.size(); ++i) {
        ans += (vertices_[(i + 1) % vertices_.size()] - vertices_[i]).abs();
    }
    return ans;
}

double Polygon::area() const {
    double ans = 0.0;
    for (size_t i = 0; i < vertices_.size(); ++i) {
        double deltax = vertices_[(i + 1) % vertices_.size()].x - vertices_[i].x;
        double deltay = vertices_[(i + 1) % vertices_.size()].y + vertices_[i].y;
        ans += deltax * deltay;
    }
    ans /= 2.0;
    return std::abs(ans);
}

bool Polygon::isEqualVertices(const std::vector<Point>& second_vertices_) const {
    if (second_vertices_.size() != vertices_.size()) {
        return false;
    }
    size_t similar_to_first_ind = vertices_.size();
    for (size_t i = 0; i < second_vertices_.size(); ++i) {
        if (vertices_[0] == second_vertices_[i]) {
            similar_to_first_ind = i;
            break;
        }
    }


    if (similar_to_first_ind == vertices_.size()) {
        return false;
    }
    bool is_found = false;
    for (auto& flag: {-1, 1}) {
        is_found = true;
        for (size_t i = 0; i < vertices_.size(); ++i) {
            size_t sec_ind;
            if (flag == 1) {
                sec_ind = (similar_to_first_ind + i) % second_vertices_.size();
            } else {
                sec_ind = (similar_to_first_ind + second_vertices_.size() - i) % second_vertices_.size();
            }

            if (vertices_[i] != second_vertices_[sec_ind]) {
                is_found = false;
                break;
            }
        }
        if (is_found) {
            return true;
        }
    }
    return is_found;
}

bool Polygon::operator==(const Shape& other) const {
    if (!dynamic_cast<const Polygon*>(&other)) {
        return false;
    }
    std::vector<Point> second_vertices_ = dynamic_cast<const Polygon*>(&other)->vertices_;
    return isEqualVertices(second_vertices_);
}

bool Polygon::operator==(const Polygon& other) const {
    return isEqualVertices(other.vertices_);
}

bool Polygon::operator!=(const Shape& other) const {
    return !(*this == other);
}

bool Polygon::operator!=(const Polygon& other) const {
    return !(*this == other);
}

std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << p.x << " " << p.y;
    return os;
}

std::array<Point, 6> Polygon::congruentHelper(size_t ft_i, size_t sc_i_pos, size_t sc_i_neg,
                                            const std::vector<Point>& second_vertices_) const {
    Point our_vec1 = vertices_[ft_i] - vertices_[(ft_i + 1) % vertices_.size()];
    Point our_vec2 = vertices_[(ft_i + 2) % vertices_.size()] - vertices_[(ft_i + 1) % vertices_.size()];
    Point other_vec1 = second_vertices_[sc_i_pos % vertices_.size()] -
                       second_vertices_[(sc_i_pos + 1) % vertices_.size()];
    Point other_vec2 = second_vertices_[(sc_i_pos + 2) % vertices_.size()] -
                       second_vertices_[(sc_i_pos + 1) % vertices_.size()];
    Point other_vec1_neg = second_vertices_[sc_i_neg % vertices_.size()] -
                           second_vertices_[(sc_i_neg - 1) % vertices_.size()];
    Point other_vec2_neg = second_vertices_[(sc_i_neg - 2) % vertices_.size()] -
                           second_vertices_[(sc_i_neg - 1) % vertices_.size()];
    return {our_vec1, our_vec2, other_vec1, other_vec2, other_vec1_neg, other_vec2_neg};
}

bool Polygon::isSimilarTo(const Shape& other) const {
    if (!dynamic_cast<const Polygon*>(&other)) {
        return false;
    }

    std::vector<Point> second_vertices_ = dynamic_cast<const Polygon*>(&other)->vertices_;
    if (vertices_.size() != second_vertices_.size()) {
        return false;
    }

    double ftCos = getCos(vertices_[0] - vertices_[1], vertices_[2] - vertices_[1]);
    for (size_t i = 0; i < second_vertices_.size(); ++i) {
        Point ft = second_vertices_[i];
        Point sc = second_vertices_[(i + 1) % second_vertices_.size()];
        Point trd = second_vertices_[(i + 2) % second_vertices_.size()];
        Point vec1 = ft - sc;
        Point vec2 = trd - sc;
        bool is_found_min = true, is_found_plus = true;
        if (isEqual(ftCos, getCos(vec1, vec2))) {
            for (size_t ft_i = 0, cnt = 0, sc_i_pos = i, sc_i_neg = i + vertices_.size() + 2; cnt < vertices_.size();
                 ++ft_i, ++sc_i_pos, --sc_i_neg, ++cnt) {
                auto [
                        our_vec1, our_vec2, other_vec1,
                        other_vec2, other_vec1_neg, other_vec2_neg
                     ] = congruentHelper(ft_i, sc_i_pos, sc_i_neg, second_vertices_);
                if (!(isEqual(getCos(our_vec1, our_vec2), getCos(other_vec1, other_vec2)) &&
                      isEqual(our_vec1.abs() / our_vec2.abs(), other_vec1.abs() / other_vec2.abs()))) {
                    is_found_plus = false;
                }
                if (!(isEqual(getCos(our_vec1, our_vec2), getCos(other_vec1_neg, other_vec2_neg)) &&
                      isEqual(our_vec1.abs() / our_vec2.abs(), other_vec1_neg.abs() / other_vec2_neg.abs()))) {
                    is_found_min = false;
                }

                if (!is_found_min && !is_found_plus) { break; }
            }
            if (is_found_min || is_found_plus) {
                return true;
            }
        }
    }
    return false;
}

bool Polygon::isCongruentTo(const Shape& other) const {
    if (!dynamic_cast<const Polygon*>(&other)) {
        return false;
    }

    std::vector<Point> second_vertices_ = dynamic_cast<const Polygon*>(&other)->vertices_;
    if (vertices_.size() != second_vertices_.size()) {
        return false;
    }

    double ftCos = getCos(vertices_[0] - vertices_[1], vertices_[2] - vertices_[1]);
    for (size_t i = 0; i < second_vertices_.size(); ++i) {
        Point ft = second_vertices_[i];
        Point sc = second_vertices_[(i + 1) % second_vertices_.size()];
        Point trd = second_vertices_[(i + 2) % second_vertices_.size()];
        Point vec1 = ft - sc;
        Point vec2 = trd - sc;
        bool is_found_min = true, is_found_plus = true;
        if (isEqual(ftCos, getCos(vec1, vec2))) {
            for (size_t ft_i = 0, cnt = 0, sc_i_pos = i, sc_i_neg = i + vertices_.size() + 2; cnt < vertices_.size();
                 ++ft_i, ++sc_i_pos, --sc_i_neg, ++cnt) {
                auto [
                    our_vec1, our_vec2, other_vec1,
                    other_vec2, other_vec1_neg, other_vec2_neg
                ] = congruentHelper(ft_i, sc_i_pos, sc_i_neg, second_vertices_);
                if (!(isEqual(getCos(our_vec1, our_vec2), getCos(other_vec1, other_vec2)) &&
                      isEqual(our_vec1.abs(), other_vec1.abs()) &&
                      isEqual(our_vec2.abs(), other_vec2.abs()))) {
                    is_found_plus = false;
                }

                if (!(isEqual(getCos(our_vec1, our_vec2), getCos(other_vec1_neg, other_vec2_neg)) &&
                      isEqual(our_vec1.abs(), other_vec1_neg.abs()) &&
                      isEqual(our_vec2.abs(), other_vec2_neg.abs()))) {
                    is_found_min = false;
                }

                if (!is_found_min && !is_found_plus) { break; }
            }
            if (is_found_min || is_found_plus) {
                return true;
            }
        }
    }
    return false;
}

bool Polygon::containsPoint(const Point& point) const {
    double angle = 0.0;
    auto getRad = [&](double nowCos, double nowSin) -> double {
        return acos(nowCos) * ((nowSin + Point::EPS >= 0) ? 1 : -1);
    };
    for (size_t i = 0; i < vertices_.size(); ++i) {
        Point ft = vertices_[i];
        Point trd = vertices_[(i + 1) % vertices_.size()];
        if (point.isOnVec(ft, trd) && Line(ft, trd).containsPoint(point)) {
            return true;
        }
        Point vec1 = ft - point;
        Point vec2 = trd - point;
        double nowCos = getCos(vec1, vec2);
        double nowSin = getSin(vec1, vec2);
        double radAngle = getRad(nowCos, nowSin);
        angle += radAngle;
    }
    return !isEqual(angle, 0.0);
}

Polygon& Polygon::rotate(const Point& center, double angle) {
    for (auto& vert: vertices_) {
        vert -= center;
        vert.degreeRotation(angle);
        vert += center;
    }
    return *this;
}

Polygon& Polygon::reflect(const Line& axis) {
    for (auto& vert: vertices_) {
        vert.reflectByAxis(axis);
    }
    return *this;
}

Polygon& Polygon::reflect(const Point& center) {
    for (auto& vert: vertices_) {
        vert.reflectByPoint(center);
    }
    return *this;
}

Polygon& Polygon::scale(const Point& center, double coefficient) {
    for (auto& point: vertices_) {
        point -= center;
        point *= coefficient;
        point += center;
    }
    return *this;
}

class Rectangle : public Polygon {
public:
    Rectangle(const Point& p1, const Point& p2, const Point& p3, const Point& p4) : Polygon(p1, p2, p3, p4) {}
    Rectangle(const Point& ft, const Point& sc, double k);

    Point center() const {
        return (vertices_[1] + vertices_[3]) / 2.0;
    }
    std::pair<Line, Line> diagonals() const {
        return {Line(vertices_[0], vertices_[2]), Line(vertices_[1], vertices_[3])};
    }
};

Rectangle::Rectangle(const Point& ft, const Point& sc, double k) {
    if (k < 1) {
        k = 1 / k;
    }
    Point diffVec = (sc - ft) * cos(atan(k));
    diffVec.radianRotation(atan(k));
    *this = Rectangle(ft, ft + diffVec, sc, sc - diffVec);
}

class Ellipse : public Shape {
protected:
    Point f1_, f2_;
    double fDistSum_ = 0;
    Ellipse() = default;
public:
    Ellipse(const Point& f1, const Point& f2, double fDistSum) : f1_(f1), f2_(f2), fDistSum_(fDistSum) {}

    std::pair<Point, Point> focuses() const {
        return {f1_, f2_};
    }

    double eccentricity() const {
        return vecLen(f1_, f2_) / fDistSum_;
    }

    Point center() const {
        return (f1_ + f2_) / 2.0;
    }

    std::pair<Line, Line> directrices() const;
    double perimeter() const override;
    double area() const override;
    bool isEqualEllipse(const Ellipse& other) const {
        return fDistSum_ == other.fDistSum_ &&
               ((f1_ == other.f1_ && f2_ == other.f2_) || (f1_ == other.f2_ && f2_ == other.f1_));
    }

    bool operator==(const Shape& other) const override;
    bool operator!=(const Shape& other) const override;
    bool operator==(const Ellipse& other) const;
    bool operator!=(const Ellipse& other) const;

    bool isSimilarTo(const Shape& other) const override;
    bool isCongruentTo(const Shape& other) const override;
    bool containsPoint(const Point& point) const override;

    Ellipse& rotate(const Point& center, double angle) override;
    Ellipse& reflect(const Point& center) override;
    Ellipse& reflect(const Line& axis) override;
    Ellipse& scale(const Point& center, double coefficient) override;
};

std::pair<Line, Line> Ellipse::directrices() const {
    Point fVec = f2_ - f1_;
    Point nVec = fVec;
    nVec.radianRotation(Point::PI / 2.0);

    Point d1 = center() + fVec * fDistSum_ / (fVec.abs() * 2.0 * eccentricity());
    Point d2 = center() - fVec * fDistSum_ / (fVec.abs() * 2.0 * eccentricity());

    return {Line{d1, d1 + nVec}, Line{d2, d2 + nVec}};
}

bool Ellipse::operator==(const Shape& other) const {
    if (!dynamic_cast<const Ellipse*>(&other)) {
        return false;
    }
    const Ellipse& other_casted = dynamic_cast<const Ellipse&>(other);
    return isEqualEllipse(other_casted);
}

bool Ellipse::operator==(const Ellipse& other) const {
    return isEqualEllipse(other);
}

bool Ellipse::operator!=(const Shape& other) const {
    return !(*this == other);
}

bool Ellipse::operator!=(const Ellipse& other) const {
    return !(*this == other);
}

bool Ellipse::isCongruentTo(const Shape& other) const {
    if (!dynamic_cast<const Ellipse*>(&other)) {
        return false;
    }
    const Ellipse& other_casted = dynamic_cast<const Ellipse&>(other);
    return isEqual((f2_ - f1_).abs(), (other_casted.f2_ - other_casted.f1_).abs()) &&
           isEqual(fDistSum_, other_casted.fDistSum_);
}

bool Ellipse::isSimilarTo(const Shape& other) const {
    if (!dynamic_cast<const Ellipse*>(&other)) {
        return false;
    }
    const Ellipse& other_casted = dynamic_cast<const Ellipse&>(other);
    return isEqual((f2_ - f1_).abs() / fDistSum_, (other_casted.f2_ - other_casted.f1_).abs() / other_casted.fDistSum_);
}

bool Ellipse::containsPoint(const Point& point) const {
    return vecLen(f1_, point) + vecLen(f2_, point) - Point::EPS <= fDistSum_;
}

Ellipse& Ellipse::rotate(const Point& center, double angle) {
    f1_ -= center;
    f2_ -= center;
    f1_.degreeRotation(angle);
    f2_.degreeRotation(angle);
    f1_ += center;
    f2_ += center;
    return *this;
}

Ellipse& Ellipse::reflect(const Line& axis) {
    f1_.reflectByAxis(axis);
    f2_.reflectByAxis(axis);
    return *this;
}

Ellipse& Ellipse::reflect(const Point& center) {
    f1_.reflectByPoint(center);
    f2_.reflectByPoint(center);
    return *this;
}

Ellipse& Ellipse::scale(const Point& center, double coefficient) {
    f1_ -= center;
    f2_ -= center;
    f1_ *= coefficient;
    f2_ *= coefficient;
    f1_ += center;
    f2_ += center;
    fDistSum_ *= coefficient;
    return *this;
}

double Ellipse::area() const {
    Point fVec = f2_ - f1_;
    double res = sqrt((fDistSum_ / 2.0) * (fDistSum_ / 2.0) - (fVec.abs() / 2.0) * (fVec.abs() / 2.0));
    return Point::PI * fDistSum_ * res / 2.0;
}

double Ellipse::perimeter() const {
    Point fVec = f2_ - f1_;
    double halfDistSum = fDistSum_ / 2.0;
    double fVecLenHalf = fVec.abs() / 2.0;
    double res = sqrt(halfDistSum * halfDistSum - fVecLenHalf * fVecLenHalf);
    double tmpCf = (halfDistSum - res) * (halfDistSum - res) / ((halfDistSum + res) * (halfDistSum + res));
    tmpCf *= 3.0;
    return (halfDistSum + res) * (1.0 + tmpCf / (10.0 + sqrt(4.0 - tmpCf))) * Point::PI;
}

class Circle : public Ellipse {
public:
    Circle(const Point& c, double rad) : Ellipse(c, c, 2.0 * rad) {}
    Circle(const Point& ft, const Point& sc, const Point& trd) {
        Point f12 = Line((ft + sc) / 2, (ft + sc) / 2 + Line(ft, sc).getNormal()).crossPoint(
            Line((ft + trd) / 2, (ft + trd) / 2 + Line(ft, trd).getNormal())
        );
        f1_ = f2_ = f12;
        fDistSum_ = 2.0 * (f12 - ft).abs();

    }
    double radius() const {
        return fDistSum_ / 2.0;
    }
};

class Square : public Rectangle {
public:
    Square(const Point& ft, const Point& sc) : Rectangle(ft, sc, 1.0) {}
    Circle inscribedCircle() const {
        return {center(), (vertices_[1] - vertices_[0]).abs() / 2.0};
    }
    Circle circumisedCircle() const {
        return {center(), (vertices_[1] - vertices_[0]).abs() / sqrt(2.0)};
    }
};

class Triangle : public Polygon {
private:
public:
    Triangle(const Point& ft, const Point& sc, const Point& trd) : Polygon(ft, sc, trd) {}
    Circle circumscribedCircle() const {
        return {vertices_[0], vertices_[1], vertices_[2]};
    }
    Point centroid() const {
        return (vertices_[0] + vertices_[1] + vertices_[2]) / 3.0;
    }
    Point orthocenter() const {
        return vertices_[0] + vertices_[1] + vertices_[2] - (circumscribedCircle().center() * 2.0);
    }
    Circle ninePointsCircle() const {
        return {(circumscribedCircle().center() + orthocenter()) / 2.0, circumscribedCircle().radius() / 2.0};
    }
    Line EulerLine() const {
        return {ninePointsCircle().center(), centroid()};
    }
    Circle inscribedCircle() const;
};

Circle Triangle::inscribedCircle() const {
    Point p1 =
        vertices_[0] + (vertices_[1] - vertices_[0]) * (perimeter() / 2.0 - (vertices_[1] - vertices_[2]).abs()) /
                       (vertices_[1] - vertices_[0]).abs();
    Point p2 =
        vertices_[0] + (vertices_[2] - vertices_[0]) * (perimeter() / 2.0 - (vertices_[1] - vertices_[2]).abs()) /
                       (vertices_[2] - vertices_[0]).abs();
    Point p3 =
        vertices_[1] + (vertices_[2] - vertices_[1]) * (perimeter() / 2.0 - (vertices_[2] - vertices_[0]).abs()) /
                       (vertices_[2] - vertices_[1]).abs();
    return {p1, p2, p3};
}
