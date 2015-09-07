#include"vector.h"
#include"point.h"
#include"EasyBMP.h"
#include"perlin.h"
#include<fstream>
#include<string>
using namespace std;


class vector{
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

class point{
	private:
	float x,y,z;
public:
	point(){};
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
};

class Quaternion{
private:
	float group[4][4];
public:
	Quaternion(){
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				group[i][j]=0;
		group[0][0]=1;
		group[1][1]=1;
		group[2][2]=1;
		group[3][3]=1;
	}
	Quaternion(float x,float y,float z){//displacement
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				group[i][j]=0;
		group[0][0]=1;
		group[1][1]=1;
		group[2][2]=1;
		group[3][3]=1;
		group[0][3]=x;
		group[1][3]=y;
		group[2][3]=z;
	}
	Quaternion(char axe,float degree){
		float angle=degree*3.1415/180;
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				group[i][j]=0;
		group[3][3]=1;
		if(axe=='x'){
			group[0][0]=1;
			group[1][1]=cos(angle);
			group[1][2]=-sin(angle);
			group[2][1]=sin(angle);
			group[2][2]=cos(angle);
		}
		else if(axe=='y'){
			group[1][1]=1;
			group[0][0]=cos(angle);
			group[0][2]=sin(angle);
			group[2][0]=-sin(angle);
			group[2][2]=cos(angle);
		}
		else{//axe=='c'
			group[2][2]=1;
			group[0][0]=cos(angle);
			group[0][1]=-sin(angle);
			group[1][0]=sin(angle);
			group[1][1]=cos(angle);
		}
	}
	float getValueAt(int i,int j){
		return group[i][j];
	}
	void setValueAt(int i,int j,float value){
		group[i][j]=value;
	}
	Quaternion operator *(Quaternion q){
		Quaternion temp;
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				float value=0;
				for(int k=0;k<4;k++)
					value+=group[i][k]*q.getValueAt(k,j);
				temp.setValueAt(i,j,value);
			}
		}
		return temp;
	}
	point operator %(point v){
		float a,b,c;
		a=group[0][0]*v.getX()+group[0][1]*v.getY()+group[0][2]*v.getZ()+group[0][3];
		b=group[1][0]*v.getX()+group[1][1]*v.getY()+group[1][2]*v.getZ()+group[1][3];
		c=group[2][0]*v.getX()+group[2][1]*v.getY()+group[2][2]*v.getZ()+group[2][3];
		point temp(a,b,c);
		return temp;
	}
};

class color{
private:
	float r,g,b;
public:
	color(){}
	color(float red,float green,float blue){
		r=red;
		g=green;
		b=blue;
	}
	void setValue(float red,float green,float blue){
		r=red;
		g=green;
		b=blue;
	}
	float getR(){
		return r;
	}
	float getG(){
		return g;
	}
	float getB(){
		return b;
	}
	color operator *(color c){
		float red,green,blue;
		red=r*c.getR();
		green=g*c.getG();
		blue=b*c.getB();
		return color(red,green,blue);
	}
	color operator +(color c){
		float red,green,blue;
		red=r+c.getR();
		green=g+c.getG();
		blue=b+c.getB();
		return color(red,green,blue);
	}
	color operator %(float f){
		float red,green,blue;
		red=r*f;
		green=g*f;
		blue=b*f;
		return color(red,green,blue);
	}
	bool check(){
		if(r>1||g>1||b>1) return true;
		else return false;
	}
};

class input{
private:
	float STEP,FOVY;
	vector XYZC,VDIR,UVEC,LPOS[4];
	point EYEP;
	color BRGB,MRGB,LCOL[4];
	char file[100];
	char *FileOutput;
	int width,height,*type,length,lightNum;
	ifstream inFile;
	float *buffer,*position,*radis;
public:
	friend class camera;
	friend class voxel;
	input(){
		FileOutput=file;
		type=new int[5];
		position=new float[15];
		radis=new float[5];
		lightNum=0;
	}
	input(string fileName){
		FileOutput=file;
		inFile.open(fileName);
		readFile();
	}
	void setFile(string fileName){
		inFile.open(fileName);
		readFile();
	}

	void print(){
		cout<<FOVY<<endl;
	}
	
	void readFile(){
		int count=0;
		string str;
		char group[30];
		while(!inFile.eof()){
			inFile>>str;
			if(str=="STEP"){
				inFile>>STEP;
			}
			else if(str=="XYZC"){//need adjustment
				int a,b,c;
				inFile>>a;
				inFile>>b;
				inFile>>c;
				buffer=new float[a*b*c];
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
			}
			else if(str=="LPOS"){
				float a,b,c;
				inFile>>a;
				inFile>>b;
				inFile>>c;
				LPOS[lightNum].setValue(a,b,c);
			}
			else if(str=="LCOL"){
				float a,b,c;
				inFile>>a;
				inFile>>b;
				inFile>>c;
				LCOL[lightNum++].setValue(a,b,c);
			}
			else if(str!=""){
				length=atoi(str.c_str());
				for(int i=0;i<length;i++){
					inFile>>str;
					if(str=="sphere") type[i]=0;
					else if(str=="cloud") type[i]=1;
					else type[i]=2;

					inFile>>position[3*i];inFile>>position[3*i+1];inFile>>position[3*i+2];
					inFile>>radis[i];
				}
			}
		}
	}
};

class voxel{
private:
	float **lightBuffer,*densityBuffer,step;
	int sca,count,m,pri;
	float *position,*radis;
public:
	friend class camera;
	voxel(){
		lightBuffer=nullptr;
		densityBuffer=nullptr;
	}
	void setParameter(int scale,float s,float *buffer,int mode,int primitive,float *pos,float *rad){
		count=0;
		sca=scale;
		step=s;
		m=mode;
		pri=primitive;
		position=pos;
		radis=rad;
		if(lightBuffer==nullptr){
			lightBuffer=new float*[3];
			for(int i=0;i<3;i++)
				lightBuffer[i]=new float[sca*sca*sca];
			densityBuffer=buffer;
		}
	}

	void initLightBuffer(int n){
		for(int i=0;i<sca*sca*sca;i++){
			lightBuffer[n][i]=-1;
		}//i
	}

	void setLightBuffer(vector *p,int n){
		if(m>2){
			for(int i=0;i<sca;i++){
				for(int j=0;j<sca;j++){
					for(int k=0;k<sca;k++){
						setLightBufferValue(p[n],n,i,j,k);
					}//k
				}//j
			}//i
		cout<<"light buffer computation complete"<<endl;
		}
	}

	void initBuffer(float fre,float amp,int *type,int length){
		int num=sca*sca*sca;
		for(int i=0;i<num;i++) densityBuffer[i]=0;
		for(int i=0;i<length;i++){
			if(type[i]==0)
				initSphere(i,sca,position,radis);
			else if(type[i]==1)
				initCloud(i,sca,position,radis,fre,amp);
			else
				initPyroclastic(i,sca,position,radis,fre,amp);
		}
	}

	float quadProportion(float x,float y,float z,float r){
		float rate,max;
		if(fabs(x)>fabs(y)) max=fabs(x);
		else max=fabs(y);
		if(fabs(max)<fabs(z)) max=fabs(z);
		rate=max/r;
		return rate;
	}

	float cylinderProportion(float x,float y,float z,float r){
		float rate,length,depth=0.3;
		rate=fabs(y)/depth;
		length=sqrt(x*x+z*z);
		if(rate==0) return length/r;
		if(length/rate>1) return length/r;
		else return fabs(y)/depth/r;
	}

	float ovalProportion(float x,float y,float z,float r){
		float rate,length,a=1,b=0.4,c=0.8;
		length=sqrt(x*x/(a*a)+y*y/(b*b)+z*z/(c*c));
		rate=length/r;
		return rate;
	}

	void initPyroclastic(int n,int scale,float *position,float *radis,float fre, float amp){
		int num=scale*scale*scale;
		int x,y,z,rx,ry,rz;
		float d,fbm,result;
		rx=scale*position[3*n];
		ry=scale*position[3*n+1];
		rz=scale*position[3*n+2];
		Perlin p(1,fre,amp,0.5);
		for(int i=0;i<num;i++){
			x=i/(scale*scale);
			y=(i%(int)(scale*scale))/scale;
			z=i-x*scale*scale-y*scale;
			z=-z;
			d=sqrt((x-rx)*(x-rx)+(y-ry)*(y-ry)+(z-rz)*(z-rz));
		
			fbm=p.Get(1.0*(x-rx)/scale,1.0*(y-ry)/scale,1.0*(z-rz)/scale);
			if(pri==1)
				result=(radis[n]*scale-fabs(d/radis[n]/scale)+fabs(fbm))*scale;
			else if(pri==2)
				result=(radis[n]*scale-fabs(quadProportion(1.0*(x-rx)/scale,1.0*(y-ry)/scale,1.0*(z-rz)/scale,radis[n]))
				+fabs(fbm))*scale;
			else if(pri==3)
				result=(radis[n]*scale-fabs(cylinderProportion(1.0*(x-rx)/scale,1.0*(y-ry)/scale,1.0*(z-rz)/scale,radis[n]))
				+fabs(fbm))*scale;
			else
				result=(radis[n]*scale-fabs(ovalProportion(1.0*(x-rx)/scale,1.0*(y-ry)/scale,1.0*(z-rz)/scale,radis[n]))
				+fabs(fbm))*scale;
			if(result>0)
				densityBuffer[i]+=result;
			if(densityBuffer[i]>scale) densityBuffer[i]=scale;

			//if(i%10000==0) cout<<i/10000<<endl;
		}
	}

	void initCloud(int n,int scale,float *position,float *radis,float fre,float amp){
		int num=scale*scale*scale;
		int x,y,z,rx,ry,rz;
		float d,fbm,result;
		float oct=log(scale)/log(2);
		rx=scale*position[3*n];
		ry=scale*position[3*n+1];
		rz=scale*position[3*n+2];
		Perlin p(oct,fre,amp,3.5);//2-12,1;7,0-2
		for(int i=0;i<num;i++){
			x=i/(scale*scale);
			y=(i%(int)(scale*scale))/scale;
			z=i-x*scale*scale-y*scale;
			z=-z;
			d=sqrt((x-rx)*(x-rx)+(y-ry)*(y-ry)+(z-rz)*(z-rz));
			
			fbm=p.Get(1.0*(x-rx)/scale,1.0*(y-ry)/scale,1.0*(z-rz)/scale);
			if(pri==1)
				result=(fbm+1-d/(scale*radis[n]))*scale;
			else if(pri==2)
				result=(fbm+1-quadProportion(1.0*(x-rx)/scale,1.0*(y-ry)/scale,1.0*(z-rz)/scale,radis[n]))*scale;
			else if(pri==3)
				result=(fbm+1-cylinderProportion(1.0*(x-rx)/scale,1.0*(y-ry)/scale,1.0*(z-rz)/scale,radis[n]))*scale;
			else 
				result=(fbm+1-ovalProportion(1.0*(x-rx)/scale,1.0*(y-ry)/scale,1.0*(z-rz)/scale,radis[n]))*scale;
			if(result>0)
			densityBuffer[i]+=result;
			if(densityBuffer[i]>scale) densityBuffer[i]=scale;

			//if(i%10000==0) cout<<i/10000<<endl;
		}
	}

	void initSphere(int n,int scale,float *position,float *radis){
		int num=scale*scale*scale;
		int x,y,z,rx,ry,rz,dx,dy,dz;
		float d;
		rx=scale*position[3*n];
		ry=scale*position[3*n+1];
		rz=scale*position[3*n+2];
		for(int i=0;i<num;i++){
			x=i/(scale*scale);
			y=(i%(int)(scale*scale))/scale;
			z=i-x*scale*scale-y*scale;
			z=-z;
			d=sqrt((x-rx)*(x-rx)+(y-ry)*(y-ry)+(z-rz)*(z-rz));

			if(d<scale*radis[n])
				densityBuffer[i]+=scale*(1-d/(scale*radis[n]));
			if(densityBuffer[i]>scale) densityBuffer[i]=scale;
			//if(i>999000) cout<<i<<endl;
		}
	}
	
	float TrilinearInterpolation(float rx,float ry,float rz,float *densityBuffer){
		int mx,my,mz;//integer form of rxyz
		float px,py,pz;//centter of mxyz
		float dx,dy,dz;//distance between rxyz and pxyz
		int width,height,depth;
		width=height=depth=sca;
		mx=floor(rx*width);
		my=floor(ry*height);
		mz=floor(rz*depth);
		px=(0.5+mx)/width;
		py=(0.5+my)/height;
		pz=(0.5+mz)/depth;
		dx=fabs(px-rx);
		dy=fabs(py-ry);
		dz=fabs(pz-rz);
		mx+=width/2;
		my+=height/2;
		mz+=depth/2;

		if((dx*1000*width>=1)&&(dy*1000*height>=1)&&(dz*1000*depth>=1)){//common form
			if(my==0||my==height-1||mx==0||mx==width-1||mz==0||mz==depth-1) 
				return densityBuffer[mx*height*depth+my*depth+mz];
			float wy1,wy2,wx1,wx2,wz1,wz2,wdensity,sum;
			if(rx>px&&ry>py&&rz>pz){
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				sum=wx1*wy1*wz1+wx1*wy1*wz2+wx1*wy2*wz1+wx1*wy2*wz2+wx2*wy1*wz1+wx2*wy1*wz2+wx2*wy2*wz1+wx2*wy2*wz2;
				wdensity=wy1*wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy1*wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz+1]+
					wy2*wx1*wz1/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz]+wy2*wx1*wz2/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz+1]+
					wy1*wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz]+wy1*wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz+1]+
					wy2*wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+(my+1)*depth+mz]+wy2*wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+(my+1)*depth+mz+1];		
			}
			else if(rx>px&&ry>py&&rz<=pz){
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				sum=wx1*wy1*wz1+wx1*wy1*wz2+wx1*wy2*wz1+wx1*wy2*wz2+wx2*wy1*wz1+wx2*wy1*wz2+wx2*wy2*wz1+wx2*wy2*wz2;
				wdensity=wy1*wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy1*wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz-1]+
					wy2*wx1*wz1/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz]+wy2*wx1*wz2/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz-1]+
					wy1*wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz]+wy1*wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz-1]+
					wy2*wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+(my+1)*depth+mz]+wy2*wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+(my+1)*depth+mz-1];		
			}
			else if(rx>px&&ry<=py&&rz>pz){
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				sum=wx1*wy1*wz1+wx1*wy1*wz2+wx1*wy2*wz1+wx1*wy2*wz2+wx2*wy1*wz1+wx2*wy1*wz2+wx2*wy2*wz1+wx2*wy2*wz2;
				wdensity=wy1*wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy1*wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz+1]+
					wy2*wx1*wz1/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz]+wy2*wx1*wz2/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz+1]+
					wy1*wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz]+wy1*wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz+1]+
					wy2*wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+(my-1)*depth+mz]+wy2*wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+(my-1)*depth+mz+1];
			}
			else if(rx>px&&ry<=py&&rz<=pz){
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				sum=wx1*wy1*wz1+wx1*wy1*wz2+wx1*wy2*wz1+wx1*wy2*wz2+wx2*wy1*wz1+wx2*wy1*wz2+wx2*wy2*wz1+wx2*wy2*wz2;
				wdensity=wy1*wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy1*wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz-1]+
					wy2*wx1*wz1/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz]+wy2*wx1*wz2/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz-1]+
					wy1*wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz]+wy1*wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz-1]+
					wy2*wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+(my-1)*depth+mz]+wy2*wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+(my-1)*depth+mz-1];
			}
			else if(rx<=px&&ry>py&&rz>pz){
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				sum=wx1*wy1*wz1+wx1*wy1*wz2+wx1*wy2*wz1+wx1*wy2*wz2+wx2*wy1*wz1+wx2*wy1*wz2+wx2*wy2*wz1+wx2*wy2*wz2;
				wdensity=wy1*wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy1*wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz+1]+
					wy2*wx1*wz1/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz]+wy2*wx1*wz2/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz+1]+
					wy1*wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz]+wy1*wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz+1]+
					wy2*wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+(my+1)*depth+mz]+wy2*wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+(my+1)*depth+mz+1];
			}
			else if(rx<=px&&ry>py&&rz<=pz){
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				sum=wx1*wy1*wz1+wx1*wy1*wz2+wx1*wy2*wz1+wx1*wy2*wz2+wx2*wy1*wz1+wx2*wy1*wz2+wx2*wy2*wz1+wx2*wy2*wz2;
				wdensity=wy1*wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy1*wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz-1]+
					wy2*wx1*wz1/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz]+wy2*wx1*wz2/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz-1]+
					wy1*wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz]+wy1*wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz-1]+
					wy2*wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+(my+1)*depth+mz]+wy2*wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+(my+1)*depth+mz-1];
			}
			else if(rx<=px&&ry<=py&&rz>pz){
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				sum=wx1*wy1*wz1+wx1*wy1*wz2+wx1*wy2*wz1+wx1*wy2*wz2+wx2*wy1*wz1+wx2*wy1*wz2+wx2*wy2*wz1+wx2*wy2*wz2;
				wdensity=wy1*wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy1*wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz+1]+
					wy2*wx1*wz1/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz]+wy2*wx1*wz2/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz+1]+
					wy1*wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz]+wy1*wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz+1]+
					wy2*wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+(my-1)*depth+mz]+wy2*wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+(my-1)*depth+mz+1];
			}
			else{
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				sum=wx1*wy1*wz1+wx1*wy1*wz2+wx1*wy2*wz1+wx1*wy2*wz2+wx2*wy1*wz1+wx2*wy1*wz2+wx2*wy2*wz1+wx2*wy2*wz2;
				wdensity=wy1*wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy1*wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz-1]+
					wy2*wx1*wz1/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz]+wy2*wx1*wz2/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz-1]+
					wy1*wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz]+wy1*wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz-1]+
					wy2*wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+(my-1)*depth+mz]+wy2*wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+(my-1)*depth+mz-1];
			}
			return wdensity;
		}

		if((dx*1000*width<1)&&(dy*1000*height<1)&&(dz*1000*depth<1)){//nearly coincidence
			return densityBuffer[mx*height*depth+my*depth+mz];
		}

		if((dx*1000*width>=1)&&(dy*1000*height<1)&&(dz*1000*depth<1)){//only coincide x axe
			if(mx==0||mx==width-1) return densityBuffer[mx*height*depth+my*depth+mz];
			float wx1,wx2,wdensity;
			if(rx>px){
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				wdensity=wx1*densityBuffer[mx*height*depth+my*depth+mz]+wx2*densityBuffer[(mx+1)*height*depth+my*depth+mz];
				return wdensity;
			}
			else{
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				wdensity=wx1*densityBuffer[mx*height*depth+my*depth+mz]+wx2*densityBuffer[(mx-1)*height*depth+my*depth+mz];
				return wdensity;
			}
		}

		if((dx*1000*width<1)&&(dy*1000*height>=1)&&(dz*1000*depth<1)){//only coincide y axe
			if(my==0||my==height-1) return densityBuffer[mx*height*depth+my*depth+mz];
			float wy1,wy2,wdensity;
			if(ry>py){
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wdensity=wy1*densityBuffer[mx*height*depth+my*depth+mz]+wy2*densityBuffer[mx*height*depth+(my+1)*depth+mz];
				return wdensity;
			}
			else{
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wdensity=wy1*densityBuffer[mx*height*depth+my*depth+mz]+wy2*densityBuffer[mx*height*depth+(my-1)*depth+mz];
				return wdensity;
			}
		}

		if((dx*1000*width<1)&&(dy*1000*height<1)&&(dz*1000*depth>=1)){//only coincide z axe
			if(mz==0||mz==depth-1) return densityBuffer[mx*height*depth+my*depth+mz];
			float wz1,wz2,wdensity;
			if(rz>pz){
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				wdensity=wz1*densityBuffer[mx*height*depth+my*depth+mz]+wz2*densityBuffer[mx*height*depth+my*depth+mz+1];
				return wdensity;
			}
			else{
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				wdensity=wz1*densityBuffer[mx*height*depth+my*depth+mz]+wz2*densityBuffer[mx*height*depth+my*depth+mz-1];
				return wdensity;
			}
		}

		if((dx*1000*width<1)&&(dy*1000*height>=1)&&(dz*1000*depth>=1)){//only coincide yz plain
			if(my==0||my==height-1||mz==0||mz==depth-1) return densityBuffer[mx*height*depth+my*depth+mz];
			float wy1,wy2,wz1,wz2,wdensity,sum;
			if(ry>py&&rz>pz){
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				sum=wy1*wz1+wy2*wz1+wy1*wz2+wy2*wz2;
				wdensity=wy1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy2*wz1/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz]+
					wy1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz+1]+wy2*wz2/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz+1];
				return wdensity;
			}
			else if(ry>py&&rz<=pz){
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				sum=wy1*wz1+wy2*wz1+wy1*wz2+wy2*wz2;
				wdensity=wy1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy2*wz1/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz]+
					wy1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz-1]+wy2*wz2/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz-1];
				return wdensity;
			}
			else if(ry<=py&&rz>pz){
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				sum=wy1*wz1+wy2*wz1+wy1*wz2+wy2*wz2;
				wdensity=wy1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy2*wz1/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz]+
					wy1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz+1]+wy2*wz2/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz+1];
				return wdensity;
			}
			else{
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				sum=wy1*wz1+wy2*wz1+wy1*wz2+wy2*wz2;
				wdensity=wy1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy2*wz1/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz]+
					wy1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz-1]+wy2*wz2/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz-1];
				return wdensity;
			}
		}

		if((dx*1000*width>=1)&&(dy*1000*height<1)&&(dz*1000*depth>=1)){//only coincide xz plain
			if(mx==0||mx==width-1||mz==0||mz==depth-1) return densityBuffer[mx*height*depth+my*depth+mz];
			float wx1,wx2,wz1,wz2,wdensity,sum;
			if(rx>px&&rz>pz){
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				sum=wx1*wz1+wx2*wz1+wx1*wz2+wx2*wz2;
				wdensity=wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz]+
					wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz+1]+wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz+1];
				return wdensity;
			}
			else if(rx>px&&rz<=pz){
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				sum=wx1*wz1+wx2*wz1+wx1*wz2+wx2*wz2;
				wdensity=wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wx2*wz1/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz]+
					wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz-1]+wx2*wz2/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz-1];
				return wdensity;
			}
			else if(rx<=px&&rz>pz){
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				wz1=(1-depth*(rz-pz));
				wz2=1-wz1;
				sum=wx1*wz1+wx2*wz1+wx1*wz2+wx2*wz2;
				wdensity=wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz]+
					wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz+1]+wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz+1];
				return wdensity;
			}
			else{
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				wz1=(1-depth*(pz-rz));
				wz2=1-wz1;
				sum=wx1*wz1+wx2*wz1+wx1*wz2+wx2*wz2;
				wdensity=wx1*wz1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wx2*wz1/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz]+
					wx1*wz2/sum*densityBuffer[mx*height*depth+my*depth+mz-1]+wx2*wz2/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz-1];
				return wdensity;
			}
		}

		if((dx*1000*width>=1)&&(dy*1000*height>=1)&&(dz*1000*depth<1)){//only coincide xy plain
			if(my==0||my==height-1||mx==0||mx==width-1) return densityBuffer[mx*height*depth+my*depth+mz];
			float wy1,wy2,wx1,wx2,wdensity,sum;
			if(ry>py&&rx>px){
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				sum=wy1*wx1+wy2*wx1+wy1*wx2+wy2*wx2;
				wdensity=wy1*wx1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy2*wx1/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz]+
					wy1*wx2/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz]+wy2*wx2/sum*densityBuffer[(mx+1)*height*depth+(my+1)*depth+mz];
				return wdensity;
			}
			else if(ry>py&&rx<=px){
				wy1=(1-height*(ry-py));
				wy2=1-wy1;
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				sum=wy1*wx1+wy2*wx1+wy1*wx2+wy2*wx2;
				wdensity=wy1*wx1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy2*wx1/sum*densityBuffer[mx*height*depth+(my+1)*depth+mz]+
					wy1*wx2/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz]+wy2*wx2/sum*densityBuffer[(mx-1)*height*depth+(my+1)*depth+mz];
				return wdensity;
			}
			else if(ry<=py&&rx>px){
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wx1=(1-width*(rx-px));
				wx2=1-wx1;
				sum=wy1*wx1+wy2*wx1+wy1*wx2+wy2*wx2;
				wdensity=wy1*wx1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy2*wx1/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz]+
					wy1*wx2/sum*densityBuffer[(mx+1)*height*depth+my*depth+mz]+wy2*wx2/sum*densityBuffer[(mx+1)*height*depth+(my-1)*depth+mz];
				return wdensity;
			}
			else{
				wy1=(1-height*(py-ry));
				wy2=1-wy1;
				wx1=(1-width*(px-rx));
				wx2=1-wx1;
				sum=wy1*wx1+wy2*wx1+wy1*wx2+wy2*wx2;
				wdensity=wy1*wx1/sum*densityBuffer[mx*height*depth+my*depth+mz]+wy2*wx1/sum*densityBuffer[mx*height*depth+(my-1)*depth+mz]+
					wy1*wx2/sum*densityBuffer[(mx-1)*height*depth+my*depth+mz]+wy2*wx2/sum*densityBuffer[(mx-1)*height*depth+(my-1)*depth+mz];
				return wdensity;
			}
		}
	}

	void setLightBufferValue(vector p,int n,int xi,int yj,int zk){
		float x,y,z,rx,ry,rz,dx,dy,dz,result=1,px,py,pz;
		int height,width,depth;
		height=width=depth=sca;
		//count++;
		px=p.getX();py=p.getY();pz=p.getZ();
		x=1.0/width*(xi+0.5)-0.5;
		y=1.0/height*(yj+0.5)-0.5;
		z=1.0/depth*(zk+0.5)-0.5;//xyz center coordinate of the voxel to be calculate
		dx=1.0/width/2;
		dy=1.0/height/2;
		dz=1.0/depth/2;//dxyz represent half length of a voxel
		vector v=(point(px,py,pz)-point(x,y,z)).unitVecotr();//light vector
		rx=x+step*v.getX();
		ry=y+step*v.getY();
		rz=z+step*v.getZ();//path along the light
		while(rx>-0.5&&rx<0.5&&ry>-0.5&&ry<0.5&&rz>-0.5&&rz<0.5){
			int mx,my,mz;//mxyz represent the current int form position of rxyz
			mx=floor(rx*width)+width/2;
			my=floor(ry*height)+height/2;
			mz=floor(rz*depth)+depth/2;
			if(m<=2){//simple Interpolation
				result*=pow(2.718,-densityBuffer[mx*height*depth+my*depth+mz]*step);
			}
			else 
				result*=pow(2.718,-TrilinearInterpolation(rx,ry,rz,densityBuffer)*step);
			rx+=step*v.getX();
			ry+=step*v.getY();
			rz+=step*v.getZ();
		}//while
		lightBuffer[n][height*depth*xi+depth*yj+zk]=result;
	}
};

class camera{
private:
	input inFile;
	BMP SampleOutput,mapping;
	vector A,B,H,V,**background;
	point M,PA,PB,PM,**group,nEYEP;
	voxel vox;
	Quaternion Q,QX,QY,QZ,D;
public:
	camera(string fileName,int mode,int primitive){
		initCamera(fileName);
		initVoxel(inFile.XYZC.getX(),inFile.STEP,inFile.buffer,mode,primitive);
	}

	void initCamera(string fileName){
		inFile.setFile(fileName);
		cout<<"file import compelete"<<endl;
		cout<<"initiating camera"<<endl;

		SampleOutput.SetSize(inFile.width,inFile.height);
		A=inFile.VDIR*inFile.UVEC;
		B=A*inFile.VDIR;
		V=B.unitVecotr();
		V=V%tan(inFile.FOVY);
		V=V%inFile.VDIR.getLength();
		H=A.unitVecotr();
		H=H%tan(inFile.FOVY);
		H=H%(1.0*inFile.width/inFile.height);
		H=H%inFile.VDIR.getLength();
		M=inFile.EYEP+inFile.VDIR;

		/*Quaternion qx('x',-rx),qy('y',-ry),qz('z',-rz),d(-x,-y,-z);
		Quaternion qInverse=qx*qy*qz*d;
		point nEYEP=qInverse%inFile.EYEP;
		point nA=qInverse%(inFile.EYEP+A);
		point nB=qInverse%(inFile.EYEP+B);
		point nM=qInverse%(inFile.EYEP+inFile.VDIR);*/
		
		group=new point*[inFile.height];
		for(int i=0;i<inFile.height;i++)
			group[i]=new point[inFile.width];

		background=new vector*[inFile.height];
		for(int i=0;i<inFile.height;i++)
			background[i]=new vector[inFile.width];

		for(int i=0;i<inFile.height;i++){
			for(int j=0;j<inFile.width;j++){
				float sx,sy;
				sx=1.0*j/(inFile.width-1);
				sy=1.0*i/(inFile.height-1);
				group[inFile.height-i-1][j]=M+H%(2*sx-1)+V%(2*sy-1);

				vector v=group[inFile.height-i-1][j]-inFile.EYEP;
				v=v.unitVecotr();
				SampleOutput(j,i)->Red=(int)abs(v.getX()*255);
				SampleOutput(j,i)->Green=(int)abs(v.getY()*255);
				SampleOutput(j,i)->Blue=(int)abs(v.getZ()*255);

				background[inFile.height-1-i][j]=v;
			}
		}

		SampleOutput.WriteToFile("SAMPLE.bmp");
		nEYEP=inFile.EYEP;

		for(int i=0;i<inFile.height;i++){
			for(int j=0;j<inFile.width;j++){
				float sx,sy;
				sx=1.0*j/(inFile.width-1);
				sy=1.0*i/(inFile.height-1);
				group[inFile.height-i-1][j]=M+H%(2*sx-1)+V%(2*sy-1);
				vector v=group[inFile.height-i-1][j]-nEYEP;
				v=v.unitVecotr();

				background[inFile.height-1-i][j]=v;
			}
		}
		//inFile.EYEP=nEYEP;

		cout<<"camera initated"<<endl;
	}

	void initVoxel(int w,float s,float *buffer,int mode,int primitive){
		vox.setParameter(w,s,inFile.buffer,mode,primitive,inFile.position,inFile.radis);
	}

	void move(float x,float y,float z,float rx,float ry,float rz){
		Quaternion qx('x',-rx),qy('y',-ry),qz('z',-rz),d(-x,-y,-z);
		QX=QX*qx;QY=QY*qy;QZ=QZ*qz;D=D*d;
		Q=QX*QY*QZ*D;
		nEYEP=Q%inFile.EYEP;
		point nA=Q%(inFile.EYEP+A);
		point nB=Q%(inFile.EYEP+B);
		point nM=Q%(inFile.EYEP+inFile.VDIR);

		vector VA,VB;
		VA=nA-nEYEP;
		VB=nB-nEYEP;
		V=VB.unitVecotr();
		V=V%tan(inFile.FOVY);
		V=V%inFile.VDIR.getLength();
		H=VA.unitVecotr();
		H=H%tan(inFile.FOVY);
		H=H%(1.0*inFile.width/inFile.height);
		H=H%inFile.VDIR.getLength();
		M=nM;

		for(int i=0;i<inFile.height;i++){
			for(int j=0;j<inFile.width;j++){
				float sx,sy;
				sx=1.0*j/(inFile.width-1);
				sy=1.0*i/(inFile.height-1);
				group[inFile.height-i-1][j]=nM+H%(2*sx-1)+V%(2*sy-1);

				vector v=group[inFile.height-i-1][j]-nEYEP;
				v=v.unitVecotr();
				SampleOutput(j,i)->Red=(int)abs(v.getX()*255);
				SampleOutput(j,i)->Green=(int)abs(v.getY()*255);
				SampleOutput(j,i)->Blue=(int)abs(v.getZ()*255);

				background[inFile.height-1-i][j]=v;
			}
		}
	}

	void outputMapping(char *fileName,int pixesX,int pixesY,float fre,float amp,int n){
		mapping.ReadFromFile(fileName);
		output(fre,amp,n,fileName,pixesX,pixesY);
	}

	void output(float fre,float amp,int n,char *fileName,int pixesX,int pixesY){
			vox.initBuffer(fre,amp,inFile.type,inFile.length);
				
			for(int i=0;i<inFile.lightNum;i++)
				vox.initLightBuffer(i);

			for(int i=0;i<inFile.lightNum;i++)
				vox.setLightBuffer(inFile.LPOS,i);
			int p,q,r;
			p=n/100;q=(n-100*p)/10;r=n-100*p-10*q;
			outputImage(p+3,q,r,n,fileName,pixesX,pixesY);
	}

	void outputImage(int p,int q,int r,int n,char *fileName,int pixesX,int pixesY){
		int u=0,v=0;
		float uscale,vscale;
		for(int i=0;i<inFile.height;i++){
			//if(i%10==0) cout<<"generating the "<<i<<"'s row out of "<<inFile.height<<endl;
			for(int j=0;j<inFile.width;j++){
				if(checkXIntersection(nEYEP,background[i][j],-0.5,0.5,-0.5,0.5)&&
				checkYIntersection(nEYEP,background[i][j],-0.5,0.5,-0.5,0.5)&&
				checkZIntersection(nEYEP,background[i][j],-0.5,0.5,-0.5,0.5)){
					color c=accumulatePixes(nEYEP,background[i][j],-0.5,0.5,-0.5,0.5,-0.5,0.5,v,u,fileName,pixesX,pixesY);
					SampleOutput(j,i)->Red=(int)(c.getR()*255);
					SampleOutput(j,i)->Green=(int)(c.getG()*255);
					SampleOutput(j,i)->Blue=(int)(c.getB()*255);
					v++;
				}
				else{
					SampleOutput(j,i)->Red=(int)abs(inFile.BRGB.getR()*255);
					SampleOutput(j,i)->Green=(int)abs(inFile.BRGB.getG()*255);
					SampleOutput(j,i)->Blue=(int)abs(inFile.BRGB.getB()*255);
				}
			}
			u++;
			v=0;
		}
		if(n==-1)
			SampleOutput.WriteToFile(inFile.FileOutput);
		else{
			char c[20]="output/";
			char buf1[10],buf2[10],buf3[10];
			itoa(p,buf1,10);itoa(q,buf2,10);itoa(r,buf3,10);
			c[7]=buf1[0];c[8]=buf2[0];c[9]=buf3[0];c[10]='.';c[11]='b';c[12]='m';c[13]='p';c[14]='\0';
			SampleOutput.WriteToFile(c);
		}
	}

	color accumulatePixes(point p,vector v,float xp1,float xp2,float yp1,float yp2,float zp1,float zp2,int col,int row,char *fileName,int pixesX,int pixesY){
		float rx,ry,rz,tran=1,length;
		vector nv=v.unitVecotr();
		color c(0,0,0);
		int width,height,depth;
		width=vox.sca;
		height=vox.sca;
		depth=vox.sca;
		//rx=p.getX()+(dis*vox.step)*nv.getX();
		//ry=p.getY()+(dis*vox.step)*nv.getY();
		//rz=p.getZ()+(dis*vox.step)*nv.getZ();//path along the ray
		rx=p.getX()+(vox.step)*nv.getX();
		ry=p.getY()+(vox.step)*nv.getY();
		rz=p.getZ()+(vox.step)*nv.getZ();//path along the ray
		while(tran>0){
			int mx,my,mz;//mxyz represent the int form of coordinate
			mx=rx*width;
			mx+=width/2;
			my=ry*height;
			my+=height/2;
			mz=rz*depth;
			mz+=+depth/2;
			//if(tran<0.000001) break;
			if(rx>=xp2&&nv.getX()>=0) break;
			if(rx<xp1&&nv.getX()<=0) break;
			if(ry>=yp2&&nv.getY()>=0) break;
			if(ry<yp1&&nv.getY()<=0) break;
			if(rz>=zp2&&nv.getZ()>=0) break;
			if(rz<zp1&&nv.getZ()<=0) break;
			//if(mx>=100||mx<0||my>=100||my<0||mz>=100||mz<0)
				//break;
			if(rx>=xp2||rx<xp1||ry>=yp2||ry<yp1||rz>=zp2||rz<zp1){
				rx+=vox.step*nv.getX();
				ry+=vox.step*nv.getY();
				rz+=vox.step*nv.getZ();
				continue;
			}
			if((vox.m==1||vox.m==3)&&vox.densityBuffer[mx*height*depth+my*depth+mz]==0){
				rx+=vox.step*nv.getX();
				ry+=vox.step*nv.getY();
				rz+=vox.step*nv.getZ();
				continue;
			}

			for(int i=0;i<inFile.lightNum;i++)
				if(vox.lightBuffer[i][mx*height*depth+my*depth+mz]==-1){
					vox.setLightBufferValue(inFile.LPOS[i],i,mx,my,mz);
				//vox.setLightBufferValue(-inFile.LPOS.getX(),-inFile.LPOS.getY(),-inFile.LPOS.getZ(),mx,my,mz);
			}

			if(vox.m==1){
				float dt=pow(2.718,-vox.densityBuffer[mx*height*depth+my*depth+mz]*vox.step);
				tran*=dt;
				color clight=inFile.LCOL[0]%vox.lightBuffer[0][mx*height*depth+my*depth+mz];
				for(int i=1;i<inFile.lightNum;i++)
					clight=clight+inFile.LCOL[i]%vox.lightBuffer[i][mx*height*depth+my*depth+mz];
				float r,g,b;
				r=clight.getR();g=clight.getG();b=clight.getB();
				if(r>1) r=1;
				if(g>1) g=1;
				if(b>1) b=1;
				clight.setValue(r,g,b);
				//c=c+inFile.MRGB*inFile.LCOL%vox.lightBuffer[mx*height*depth+my*depth+mz]%tran%vox.densityBuffer[mx*height*depth+my*depth+mz]%vox.step;
				c=c+inFile.MRGB%tran%(1-dt)*clight;
			}
			else if(vox.m==2){
				float dt=pow(2.718,-vox.TrilinearInterpolation(rx,ry,rz,vox.densityBuffer)*vox.step);
				tran*=dt;
				color clight=inFile.LCOL[0]%vox.lightBuffer[0][mx*height*depth+my*depth+mz];
				for(int i=1;i<inFile.lightNum;i++)
					clight=clight+inFile.LCOL[i]%vox.lightBuffer[i][mx*height*depth+my*depth+mz];
				float r,g,b;
				r=clight.getR();g=clight.getG();b=clight.getB();
				if(r>1) r=1;
				if(g>1) g=1;
				if(b>1) b=1;
				clight.setValue(r,g,b);
				//c=c+inFile.MRGB*inFile.LCOL%tran%vox.TrilinearInterpolation(rx,ry,-rz,vox.lightBuffer)%vox.TrilinearInterpolation(rx,ry,-rz,vox.densityBuffer)%vox.step;
				//c=c+inFile.MRGB*inFile.LCOL%tran%vox.TrilinearInterpolation(rx,ry,-rz,vox.lightBuffer)%(1-dt);
				c=c+inFile.MRGB*clight%tran%(1-dt);
			}
			else if(vox.m==3){
				float dt=pow(2.718,-vox.densityBuffer[mx*height*depth+my*depth+mz]*vox.step);
				tran*=dt;
				color clight=inFile.LCOL[0]%vox.TrilinearInterpolation(rx,ry,rz,vox.lightBuffer[0]);
				for(int i=1;i<inFile.lightNum;i++)
					clight=clight+inFile.LCOL[i]%vox.TrilinearInterpolation(rx,ry,rz,vox.lightBuffer[i]);
				float r,g,b;
				r=clight.getR();g=clight.getG();b=clight.getB();
				if(r>1) r=1;
				if(g>1) g=1;
				if(b>1) b=1;
				clight.setValue(r,g,b);
				c=c+inFile.MRGB*clight%tran%(1-dt);
			}
			else{
				float dt=pow(2.718,-vox.TrilinearInterpolation(rx,ry,rz,vox.densityBuffer)*vox.step);
				tran*=dt;
				color clight=inFile.LCOL[0]%vox.TrilinearInterpolation(rx,ry,rz,vox.lightBuffer[0]);
				for(int i=1;i<inFile.lightNum;i++)
					clight=clight+inFile.LCOL[i]%vox.TrilinearInterpolation(rx,ry,rz,vox.lightBuffer[i]);
				float r,g,b;
				r=clight.getR();g=clight.getG();b=clight.getB();
				if(r>1) r=1;
				if(g>1) g=1;
				if(b>1) b=1;
				clight.setValue(r,g,b);
				c=c+inFile.MRGB*clight%tran%(1-dt);
			}
			rx+=vox.step*nv.getX();
			ry+=vox.step*nv.getY();
			rz+=vox.step*nv.getZ();
		}
		if(fileName!=nullptr){
			float r,g,b;
			r=c.getR()*mapping(pixesX+col,pixesY+row)->Red/255;
			g=c.getG()*mapping(pixesX+col,pixesY+row)->Green/255;
			b=c.getB()*mapping(pixesX+col,pixesY+row)->Blue/255;
			c=color(r,g,b);
		}
		c=c+inFile.BRGB%tran;
		return c;
	}

	bool checkZIntersection(point p,vector v,float xp1,float xp2,float yp1,float yp2){
		//check Z plant first
		if(v.getX()==0&&v.getY()==0){//if parallel with Z axes
			if(p.getX()>xp1&&p.getX()<xp2&&p.getY()>yp1&&p.getY()<yp2) return true;
			else return false;
		}
		float a,b,x,y,x1,x2,y1,y2,xn,xf,yn,yf,n,f,temp;
		a=v.getX();
		b=v.getY();
		x=p.getX();
		y=p.getY();
		if(x>=xp1&&x<=xp2&&y>=yp1&&y<=yp2) return true;
		if(a==0){
			if(x<xp1||x>xp2) return false;
			else{
				if(b>0&&y<=yp2) return true;
				else if(b>0&&y>yp2) return false;
				else if(b<=0&&y<yp1) return false; 
				else return true;
			}
		}
		if(b==0){
			if(y<yp1||y>yp2) return false;
			else{
				if(a>0&&x<=xp2) return true;
				else if(a>0&&x>xp2) return false;
				else if(a<=0&&x<xp1) return false; 
				else return true;
			}
		}
		x1=x+a/b*(yp1-y);
		x2=x+a/b*(yp2-y);
		y1=xp1;
		y2=xp2;//y1,y2 represnet their x coordinate
		xn=(x1-x)*a;
		xf=(x2-x)*a;
		yn=(y1-x)*a;
		yf=(y2-x)*a;
		if(xn>xf){
			temp=xn;
			xn=xf;
			xf=temp;
		}
		if(yn>yf){
			temp=yn;
			yn=yf;
			yf=temp;
		}

		if(xn>yn) n=xn;
		else n=yn;
		if(yf<xf) f=yf;
		else f=xf;
		
		/*
		cout<<"!"<<endl;
		cout<<xn<<endl;
		cout<<xf<<endl;
		cout<<yn<<endl;
		cout<<yf<<endl;
		*/
		if(f>n&&f>0) return true;
		else return false;
	}

	bool checkYIntersection(point p,vector v,float zp1,float zp2,float xp1,float xp2){
		if(v.getX()==0&&v.getZ()==0){
			if(p.getX()>xp1&&p.getX()<xp2&&p.getZ()>zp1&&p.getZ()<zp2) return true;
			else return false;
		}
		float a,b,x,z,x1,x2,z1,z2,xn,xf,zn,zf,n,f,temp;
		a=v.getZ();
		b=v.getX();
		x=p.getX();
		z=p.getZ();
		if(x>=xp1&&x<=xp2&&z>=zp1&&z<=zp2) return true;
		if(b==0){
			if(x<xp1||x>xp2) return false;
			else{
				if(a>0&&z<=zp2) return true;
				else if(a>0&&z>zp2) return false;
				else if(a<=0&&z<=zp1) return false; 
				else return true;
			}
		}
		if(a==0){
			if(z<zp1||z>zp2) return false;
			else{
				if(b>0&&x<=xp2) return true;
				else if(b>0&&x>xp2) return false;
				else if(b<=0&&x<=zp1) return false; 
				else return true;
			}
		}
		z1=z+a/b*(xp1-x);
		z2=z+a/b*(xp2-x);
		x1=zp1;
		x2=zp2;
		xn=(x1-z)*a;
		xf=(x2-z)*a;
		zn=(z1-z)*a;
		zf=(z2-z)*a;
		if(xn>xf){
			temp=xn;
			xn=xf;
			xf=temp;
		}
		if(zn>zf){
			temp=zn;
			zn=zf;
			zf=temp;
		}

		if(xn>zn) n=xn;
		else n=zn;
		if(zf<xf) f=zf;
		else f=xf;
		
		if(f>n&&f>0) return true;
		else return false;
	}

	bool checkXIntersection(point p,vector v,float yp1,float yp2,float zp1,float zp2){
		if(v.getZ()==0&&v.getY()==0){//if parallel with Z axes
			if(p.getZ()>zp1&&p.getZ()<zp2&&p.getY()>yp1&&p.getY()<yp2) return true;
			else return false;
		}
		float a,b,z,y,z1,z2,y1,y2,zn,zf,yn,yf,n,f,temp;
		a=v.getY();
		b=v.getZ();
		z=p.getZ();
		y=p.getY();
		if(z>=zp1&&z<=zp2&&y>=yp1&&y<=yp2) return true;
		if(b==0){
			if(z<zp1||z>zp2) return false;
			else{
				if(a>0&&y<=yp2) return true;
				else if(a>0&&y>yp2) return false;
				else if(a<=0&&y<=yp1) return false; 
				else return true;
			}
		}
		if(a==0){
			if(y<yp1||y>yp2) return false;
			else{
				if(b>0&&z<=zp2) return true;
				else if(b>0&&z>zp2) return false;
				else if(b<=0&&z<=zp1) return false; 
				else return true;
			}
		}
		y1=y+a/b*(zp1-z);
		y2=y+a/b*(zp2-z);
		z1=yp1;
		z2=yp2;
		zn=(z1-y)*a;
		zf=(z2-y)*a;
		yn=(y1-y)*a;
		yf=(y2-y)*a;
		if(zn>zf){
			temp=zn;
			zn=zf;
			zf=temp;
		}
		if(yn>yf){
			temp=yn;
			yn=yf;
			yf=temp;
		}

		if(zn>yn) n=zn;
		else n=yn;
		if(yf<zf) f=yf;
		else f=zf;
		
		if(f>n&&f>0) return true;
		else return false;
	}
};
