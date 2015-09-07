#include<iostream>

/*class camera{
private:
	input inFile;
	string fileName;
	BMP SampleOutput;
	vector A,B,H,V;
	point M;
public:
	camera(string fileName){
		inFile.setFile(fileName);
		SampleOutput.SetSize(inFile.width,inFile.height);
		A=inFile.VDIR*inFile.UVEC;
		B=A*inFile.VDIR;
		V=B.unitVecotr();
		V=V%tan(inFile.FOVY);
		V=V%inFile.VDIR.getLength();
		H=A.unitVecotr();
		H=H%tan(inFile.FOVY);
		H=H%(inFile.width/inFile.height);
		H=H%inFile.VDIR.getLength();
		M.setValue(0,0,0);
		M=M+inFile.VDIR;

		point **group=new point*[inFile.height];
		for(int i=0;i<inFile.height;i++)
			group[i]=new point[inFile.width];

		for(int i=0;i<inFile.height;i++){
			for(int j=0;j<inFile.width;j++){
				float sx,sy;
				sx=1.0*j/(inFile.width-1);
				sy=1.0*i/(inFile.height-1);
				group[i][j]=M+H%(2*sx-1)+V%(2*sy-1);

				vector v=group[i][j]-inFile.EYEP;
				v=v.unitVecotr();
				SampleOutput(j,i)->Red=(int)abs(v.getX()*256);
				SampleOutput(j,i)->Green=(int)abs(v.getY()*256);
				SampleOutput(j,i)->Blue=(int)abs(v.getZ()*256);
			}
		}

		SampleOutput.WriteToFile("SAMPLE.bmp");
	}

};

camera::camera(string fileName){
	inFile.setFile(fileName);
	SampleOutput.SetSize(inFile.width,inFile.height);
	A=inFile.VDIR*inFile.UVEC;
	B=A*inFile.VDIR;
	V=B.unitVecotr();
	V=V%tan(inFile.FOVY);
	V=V%inFile.VDIR.getLength();
	H=A.unitVecotr();
	H=H%tan(inFile.FOVY);
	H=H%(inFile.width/inFile.height);
	H=H%inFile.VDIR.getLength();
	M.setValue(0,0,0);
	M=M+inFile.VDIR;

	point **group=new point*[inFile.height];
	for(int i=0;i<inFile.height;i++)
		group[i]=new point[inFile.width];

	for(int i=0;i<inFile.height;i++){
		for(int j=0;j<inFile.width;j++){
			float sx,sy;
			sx=1.0*j/(inFile.width-1);
			sy=1.0*i/(inFile.height-1);
			group[i][j]=M+H%(2*sx-1)+V%(2*sy-1);
	
			vector v=group[i][j]-inFile.EYEP;
			v=v.unitVecotr();
			SampleOutput(j,i)->Red=(int)abs(v.getX()*256);
			SampleOutput(j,i)->Green=(int)abs(v.getY()*256);
			SampleOutput(j,i)->Blue=(int)abs(v.getZ()*256);			
		}
	}

	SampleOutput.WriteToFile("SAMPLE.bmp");
}

color::color(){}

color::color(float red,float green,float blue){
	r=red;
	g=green;
	b=blue;
}
void color::setValue(float red,float green,float blue){
	r=red;
	g=green;
	b=blue;
}

point::point(){}

point::point(float a,float b,float c){
	x=a;
	y=b;
	z=c;
}

void point::setValue(float a,float b,float c){
	x=a;
	y=b;
	z=c;
}

float point::getX(){
	return x;
}

float point::getY(){
	return y;	
}

float point::getZ(){
	return z;
}

point point::operator +(vector v){
	point temp;
	temp.x=x+v.getX();
	temp.y=y+v.getY();
	temp.z=z+v.getZ();
	return temp;
}
	
vector point::operator -(point p){
	vector v(x-p.getX(),y-p.getY(),z-p.getZ());
	return v;
}

void point::print(){
	cout<<x<<" "<<y<<" "<<z<<endl;
}

input::input(){}

input::input(string fileName){
	inFile.open(fileName);
	readFile();
}

void input::setFile(string fileName){
	inFile.open(fileName);
	readFile();
}

void input::print(){
	cout<<FOVY<<endl;
}
	
void input::readFile(){
	while(true){
		string str;
		inFile>>str;
		if(str=="STEP"){
			inFile>>STEP;
		}
		else if(str=="XYZC"){//need adjustment
			float a,b,c;
			inFile>>a;
			inFile>>b;
			inFile>>c;
			XYZC.setValue(a,b,c);
		}
		else if(str=="BRGB"){
			float red,green,blue;
			inFile>>red;
			inFile>>green;
			inFile>>blue;
			BRGB.setValue(red,green,blue);
		}
		else if(str=="MRGB"){
			float red,green,blue;
			inFile>>red;
			inFile>>green;
			inFile>>blue;
			MRGB.setValue(red,green,blue);
		}
		else if(str=="FILE"){
			inFile>>FileOutput;
		}
		else if(str=="RESO"){
			inFile>>width;
			inFile>>height;
		}
		else if(str=="EYEP"){
			float a,b,c;
			inFile>>a;
			inFile>>b;
			inFile>>c;
			EYEP.setValue(a,b,c);
		}
		else if(str=="VDIR"){
			float a,b,c;
			inFile>>a;
			inFile>>b;
			inFile>>c;
			VDIR.setValue(a,b,c);
		}
		else if(str=="UVEC"){
			float a,b,c;
			inFile>>a;
			inFile>>b;
			inFile>>c;
			UVEC.setValue(a,b,c);
		}
		else if(str=="FOVY"){
			inFile>>FOVY;
			FOVY=FOVY/180*3.1416;
			cout<<FOVY<<endl;
		}
		else if(str=="LPOS"){
			float a,b,c;
			inFile>>a;
			inFile>>b;
			inFile>>c;
			LPOS.setValue(a,b,c);
		}
		else if(str=="LCOL"){
			float a,b,c;
			inFile>>a;
			inFile>>b;
			inFile>>c;
			LCOL.setValue(a,b,c);
		}
		else break;
	}
}

vector::vector(float a,float b,float c){
	x=a;
	y=b;
	z=c;
}
vector::vector(){}

float vector::getX(){
	return x;
}

float vector::getY(){
	return y;
}

float vector::getZ(){
	return z;
}

float vector::getLength(){
	return sqrt(x*x+y*y+z*z);
}

void vector::setValue(float a,float b,float c){
	x=a;
	y=b;
	z=c;
}

vector vector::operator +(vector v){
	float a,b,c;
	a=x+v.getX();
	b=y+v.getY();
	c=z+v.getZ();
	vector temp(a,b,c);
	return temp;
}

vector vector::operator *(vector v){//cross product
	float a,b,c;
	a=y*v.getZ()-z*v.getY();
	b=z*v.getX()-x*v.getZ();
	c=x*v.getY()-y*v.getX();
	vector temp(a,b,c);
	return temp;
}

vector vector::operator %(float f){//dot product
	vector temp(x*f,y*f,z*f);
	return temp;
}

vector vector::unitVecotr(){
	float res,a,b,c;
	res=sqrt(x*x+y*y+z*z);
	a=x/res;
	b=y/res;
	c=z/res;
	vector temp(a,b,c);
	return temp;
}

void vector::print(){
	cout<<"x: "<<x<<endl;
	cout<<"y: "<<y<<endl;
	cout<<"z: "<<z<<endl;
}

*/