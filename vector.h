/*class vector{
private:
	float x,y,z;
public:
	vector(float a,float b,float c){
		x=a;
		y=b;
		z=c;
	}
	vector(){}
	float getX(){
		return x;
	}
	float getY(){
		return y;
	}
	float getZ(){
		return z;
	}
	float getLength(){
		return sqrt(x*x+y*y+z*z);
	}
	void setValue(float a,float b,float c){
		x=a;
		y=b;
		z=c;
	}
	vector operator +(vector v){
		float a,b,c;
		a=x+v.getX();
		b=y+v.getY();
		c=z+v.getZ();
		vector temp(a,b,c);
		return temp;
	}
	vector operator *(vector v){//cross product
		float a,b,c;
		a=y*v.getZ()-z*v.getY();
		b=z*v.getX()-x*v.getZ();
		c=x*v.getY()-y*v.getX();
		vector temp(a,b,c);
		return temp;
	}
	vector operator %(float f){//dot product
		vector temp(x*f,y*f,z*f);
		return temp;
	}
	vector unitVecotr(){
		float res,a,b,c;
		res=sqrt(x*x+y*y+z*z);
		a=x/res;
		b=y/res;
		c=z/res;
		vector temp(a,b,c);
		return temp;
	}
	void print(){
		cout<<"x: "<<x<<endl;
		cout<<"y: "<<y<<endl;
		cout<<"z: "<<z<<endl;
	}
};
#ifndef _vector_
#define _vector_
class vector{
	private:
	float x,y,z;
public:
	vector(float a,float b,float c);
	vector();
	float getX();
	float getY();
	float getZ();
	float getLength();
	void setValue(float a,float b,float c);
	vector operator +(vector v);
	vector operator *(vector v);
	vector operator %(float f);
	vector unitVecotr();
	void print();
};
#endif
*/