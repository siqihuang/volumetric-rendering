/*class input{
private:
	float STEP,FOVY;
	vector XYZC,VDIR,UVEC,LPOS,LCOL;
	point EYEP;
	color BRGB,MRGB;
	string FileOutput;
	int width,height;
public:
	ifstream inFile;
	friend class camera;
	input(){}
	input(string fileName){
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
};
#ifndef _Input_
#define _Input_

#include"vector.h"
#include"point.h"
#include"color.h"

#include<iostream>
#include<string>
#include<fstream>

using namespace std;
class input{
private:
	float STEP,FOVY;
	vector XYZC,VDIR,UVEC,LPOS,LCOL;
	point EYEP;
	color BRGB,MRGB;
	string FileOutput;
	int width,height;
public:
	ifstream inFile;
	friend class camera;
	input();
	input(string fileName);
	void setFile(string fileName);
	void print();
	void readFile();
};
#endif
*/