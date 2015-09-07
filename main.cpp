/**
 * An example program that creates a 24-bit bitmap file that is 800 x 600 and makes a blue/red checkerboard pattern.
 * Uses EasyBMP
 *
 * Cory Boatright
 * University of Pennsylvania, Fall 2011
 **/

#include "camera.h"
#include "EasyBMP.h"
#include <iostream>
#include <string>
#include "perlin.h"


using namespace std;

int main(){
	string fileName;
	char textFileName[20],*textFile=textFileName;
	float x,y,z,rx,ry,rz,fre1,fre2,amp1,amp2,perlinFre,perlinAmp;
	int mode,primitive,texture,multiOutput,rotate,batch,copy,changePerlin,pixesX,pixesY;
	cout<<"Please input file name(please note only txt file is received! and you do not need to end with '.txt'): "<<endl;
	cin>>fileName;
	cout<<"please indicate if using trilinear trapolation or not,1 for no, 2 for only density buffer,";
	cout<<"3 for only light buffr, 4 for both density and light buffer"<<endl;
	cin>>mode;
	cout<<"please indicate primitives. 1 for sphere, 2 for quad, 3 for cylinder, 4 for oval"<<endl;
	cin>>primitive;
	cout<<"please indicate perlin frequency"<<endl;
	cin>>perlinFre;
	cout<<"please indicate perlin amplitude"<<endl;
	cin>>perlinAmp;
	
	cout<<"importing files"<<endl;
	camera c(fileName+".txt",mode,primitive);
	c.output(perlinFre,perlinAmp,-1,nullptr,0,0);
	while(1){
		cout<<"Do you wish to change perlin parameter? 1 for no, 2 for yes"<<endl;
		cin>>changePerlin;
		while(changePerlin==2){
			cout<<"please indicate perlin frequency"<<endl;
			cin>>perlinFre;
			cout<<"please indicate perlin amplitude"<<endl;
			cin>>perlinAmp;
			cout<<"Do you wish to change perlin parameter? 1 for no, 2 for yes"<<endl;
			cin>>changePerlin;
			c.output(perlinFre,perlinAmp,-1,nullptr,0,0);
		}
		cout<<"Do you wish to move the camera around? 1 for no, 2 for yes"<<endl;
		cin>>rotate;
		while(rotate==2){
			cout<<"input the x displacement of the camera with respect to the WORLD"<<endl;
			cin>>x;
			cout<<"input the y displacement of the camera with respect to the WORLD"<<endl;
			cin>>y;
			cout<<"input the z displacement of the camera with respect to the WORLD"<<endl;
			cin>>z;
			cout<<"input the x rotation of the camera with respect to ITSELF"<<endl;
			cin>>rx;
			cout<<"input the y rotation of the camera with respect to ITSELF"<<endl;
			cin>>ry;
			cout<<"input the z rotation of the camera with respect to ITSELF"<<endl;
			cin>>rz;
			c.move(x,y,z,rx,ry,rz);
			c.output(perlinFre,perlinAmp,-1,nullptr,0,0);
			cout<<"Do you wish to move the camera around? 1 for no, 2 for yes"<<endl;
			cin>>rotate;
		}

		cout<<"please indicate if texture mapping the output image or not, 1 for no, 2 for yes"<<endl;
		cin>>texture;
		while(texture==2){
			cout<<"You have to notice that if the picture size is small than the size of the cloud, error may occur, so you can try many times to get a satisified result"<<endl;
			cout<<"please input the texture source file name(please note only bmp file is accepted! you do not have to end with '.bmp'"<<endl;
			cin>>textFile;
			cout<<"please input the x displacement of the picture"<<endl;
			cin>>pixesX;
			cout<<"please input the y displacement of the picture"<<endl;
			cin>>pixesY;
			int i;
			for(i=0;textFileName[i]!='\0';i++);
			textFileName[i++]='.';textFileName[i++]='b';textFileName[i++]='m';
			textFileName[i++]='p';textFileName[i]='\0';
			c.outputMapping(textFile,pixesX,pixesY,perlinFre,perlinAmp,-1);
			cout<<"please indicate if texture mapping the output image or not, 1 for no, 2 for yes"<<endl;
			cin>>texture;
		}

		cout<<"do you wish to make a batch of photos? 1 for no, 2 for yes"<<endl;
		cin>>batch;
		while(batch==2){
			cout<<"Be sure to create an 'output' folder in your current file path, the program will not create one, otherwise the program will meet error"<<endl;
			cout<<"please input the final x displacement"<<endl;
			cin>>x;
			cout<<"please input the final y displacement"<<endl;
			cin>>y;
			cout<<"please input the final z displacement"<<endl;
			cin>>z;
			cout<<"please input the final x rotation"<<endl;
			cin>>rx;
			cout<<"please input the final y rotation"<<endl;
			cin>>ry;
			cout<<"please input the final z rotation"<<endl;
			cin>>rz;
			cout<<"please input the beginning perlin frequency"<<endl;
			cin>>fre1;
			cout<<"please input the ending perlin frequency"<<endl;
			cin>>fre2;
			cout<<"please input the beginning perlin amplitude"<<endl;
			cin>>amp1;
			cout<<"please input the ending perlin amplitude"<<endl;
			cin>>amp2;
			cout<<"please input the number of copies"<<endl;
			cin>>copy;

			for(int i=0;i<copy;i++){
				float unit=1.0/copy;
				float unitPerlinFre=(fre2-fre1)/copy;
				float unitPerlinAmp=(amp2-amp1)/copy;
				c.move(unit*x,unit*y,unit*z,unit*rx,unit*ry,unit*rz);
				c.output(fre1+unitPerlinFre*i,amp1+unitPerlinAmp*i,i,nullptr,0,0);
				cout<<"copy "<<i<<"complete "<<endl;
			}
			cout<<"do you wish to make a batch of photos? 1 for no, 2 for yes"<<endl;
			cin>>batch;
		}
	}
	return 0;
}