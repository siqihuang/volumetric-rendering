/*class point{
private:
	float x,y,z;
public:
	point(){}
	point(float a,float b,float c){
		x=a;
		y=b;
		z=c;
	}
	void setValue(float a,float b,float c){
		x=a;
		y=b;
		z=c;
	}
	float getX(){
		return x;
	}
	float getY(){
		return y;
	}
	float getZ(){
		return z;
	}
	point operator +(vector v){
		point temp;
		temp.x=x+v.getX();
		temp.y=y+v.getY();
		temp.z=z+v.getZ();
		return temp;
	}
	vector operator -(point p){
		vector v(x-p.getX(),y-p.getY(),z-p.getZ());
		return v;
	}
	void print(){
		cout<<x<<" "<<y<<" "<<z<<endl;
	}
};*/

#include"vector.h"
/*
class point{
	private:
	float x,y,z;
public:
	point();
	point(float a,float b,float c);
	void setValue(float a,float b,float c);
	float getX();
	float getY();
	float getZ();
	point operator +(vector v);
	vector operator -(point p);
	void print();
};
*/