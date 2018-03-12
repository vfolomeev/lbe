#pragma once
#include "math_api.h"


//Vector2d

Vector2d operator+(Vector2d a,Vector2d b){
    double x=a.x+b.x;
    double y=a.y+b.y;
    return Vector2d(x,y);
}
Vector2d operator*(double a,Vector2d b){
    return Vector2d(a*b.x,a*b.y);
}
Vector2d operator*(Vector2d a,double b){
    return Vector2d(b*a.x,b*a.y);
}
double operator*(Vector2d a,Vector2d b){
    return a.x*b.x+a.y*b.y;
}
Vector2d::Vector2d(){
    x=0;y=0;
}
Vector2d::Vector2d(double a){
    x=a;y=a;
}
Vector2d::Vector2d(double a,double b){
    x=a;y=b;
}
double Vector2d::mag(){
    return sqrt(x*x+y*y);
}
Vector2d Vector2d::operator %(Vector2d a){
    return Vector2d();
}
std::ostream& operator <<(std::ostream& cout_,const Vector2d& r){
    cout_<<" "<<r.x<<" "<<r.y;
    return cout_;
}
QTextStream &operator <<(QTextStream &s, const Vector2d &r){
    return s << " "<<r.x<<" "<<r.y;
}



