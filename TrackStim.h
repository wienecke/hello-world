#pragma once
#include <ctime>
#include <cmath>
#include <windows.h>		// Header File For Windows
#include <math.h>			// Math Library Header File
#include <stdio.h>			// Header File For Standard Input/Output
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library

#include <NIDAQmx.h>

//#include "Animation.h"      // header file for animation class to make random dot patterns.
//#include "RandomC.h"
//#include "EnforcedCorr1D.h"

//#include "dirent.h"
//#include <sys/types.h>

//#include "H5Cpp.h"
//using namespace H5;

#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else

class TrackStim
	
	//Class for creating stimuli 
{

public:
	
	int colorChannels = 3;
	int bitDepth = 8;

	float *buffer[20] = { NULL };
	float *movie[20] = { NULL };
	float *noiseBuffer[20] = { NULL };
	long *noiseTemplate[20] = { NULL };
	float *noiseValuesBuffer[20] = { NULL };
	char *filStr[20] = { NULL };
	char *filStrNoise[20] = { NULL };

	double *meshx = NULL;
	double *meshy = NULL;
	double *meshu = NULL;
	double *meshv = NULL;
	double *meshi = NULL;

	int meshxres = 0;
	int meshyres = 0;
	
	int texType[20] = { 0 };
	int texXres[20] = { 0 };
	int texYres[20] = { 0 };
	int texTres[20] = { 0 };
	int texSpacing[20] = { 0 };
	int texDir[20] = { 0 };
	int texPrecision[20] = { 0 };

	int noiseCellNumber[20] = { 0 };
	int noiseFrameNumber[20] = { 0 };
	
	char filenames[20][255] = { "" };
	char noiseNames[20][255] = { "" };


	int fileIndex = 0;

	float Sfullx = 0;
	float Sfully = 0;
	float Sfulltheta = 0;
	float Sfulltcurr = 0;
	float *noiseValues[20] = { NULL };
	float *noiseValuesFrame[20] = { NULL };

	//float noiseValues[33] = { 0 };
	bool firstError = TRUE;
	char DAQErrtmp[500];
	uInt32 dataTwo = 0;

	long int *framenumberArray = NULL;
	float *tcurrArray = NULL;
	int *boutIndArray = NULL;
	int *epochchooseArray = NULL;
	float *xArray = NULL;
	float *yArray = NULL;
	float *thetaArray = NULL;
	uInt32 *dataFrameArray = NULL;


	struct OLPARAMS {
		float mean = 0;
		float amp = 0;
		float per = 0;
		float phase = 0;
	};

	struct STIMPARAMS { // everything attached to the stimulus; array of this to get stim sequence by index: "stim epoch"
		int stimtype = 0;
		float lum = 0;
		float contrast = 0;
		float duration = 0;
		float transgain = 0;
		float rotgain = 0;
		// below: changes fly's frame with OL
		OLPARAMS trans;  // OL on fly translation & rotation -- OL on stimulus position and angle
		OLPARAMS rot;
		// below: in fly's frame
		OLPARAMS stimtrans;  // OL on stimulus translation, where applicable: halves, bars, inf corridor
		OLPARAMS stimrot;    // OL on stim rotation, where applicable: bars
		float spacing = 0;
		float spacing2 = 0;
		float density = 0;
		float tau = 0;
		float tau2 = 0;
		float arenasize = 0; 
		float arenaheight = 0; 
		int randomize = 0;
		float x_center = 0;
		float y_center = 0;
		float base_radius = 0;
		float base_extent = 0; // radius of 3d cylinders to be drawn, and their heights.
		int USEFRUSTUM = 0;
	} Stimulus;

	int boutInd = 0; // bout index. starts at 0
	int frameIndex = 0;
	int numepochs; // number of epochs in stimulus
	int epochchoose = 0; // control of which after reading them in...
	float epochtimezero = 0;
	bool firstie = TRUE;
	int displaylist = 1;
	float maxRunTime = 0; // duration of recording in seconds
	long int framenumber = 0;

	struct VIEWPOSITIONS
	{
		int x[2],y[2],h[2],w[2];  // (x,y,h,w) coords for all 4 view windows...
	} ViewPorts;

	char stimulusDataPath[255] = "";
	float * stimulusData;
	char * outputDir;
	int windowWidth, windowHeight, bits;

private:
	
	//CAnimation Rdots; // random dots class
	//RandomC r1; // random number generator, with correlation time, etc.

	// basic
	double PI;
	float ScreenUpdateRate; // throw up at Hz
	char paramfilename_two[255] = ""; // stores param file name to be written out

	LARGE_INTEGER CounterBeginTime;
	LARGE_INTEGER CounterFreq;
	FILE *monitor_stimulus_file;   // file to write for stimulus info

	void readstring(FILE *f,char *str);
	float sign(float in);

	float getDLPColor(float DLPintensity, char channel);

	float backgroundtriple[3],foregroundtriple[3],blacktriple[3],whitetriple[3];
	void setColorTriples();

	float perspELEVATION = 0;
	float perspAZIMUTH = 0;
	float perspHEAD = 0;

	float perspAlpha = 0;
	float perspBeta = 0;
	float perspGamma = 0;
	float perspDelta = 0;

	float cylinderELEVATION = 0;
	float cylinderAZIMUTH = 0;

public:

	//////////basic
	TrackStim();  // **
	virtual ~TrackStim(); // **

	//////////////set various things in the stimulus
	void setZeroTime();        // ** also done with class first defined...
	void writeStim(TaskHandle taskHandle);    // ** writes out stimulus data, do each color?
	float queryCurrTime(); // ** gets current time
	void incrementFrameNumber();
	
	STIMPARAMS initializeStimulus();
	TrackStim::STIMPARAMS * readParams(char *szFile);

	void setBackground();  // sets background color
	void setForeground();  // sets foreground color

	/////////////making real stimuli
	//void readStimulusData(); //for noise, based on HDF5 file 
	void readViewPositions();  //**

	bool drawScene();  // ** executes the drawing set by setStimType

	void drawSineCylinder(float tc);
	void drawSineCylinder_Standing();
	void drawSineCylinder_Standing_ContrastModulated(float tc);
	void drawSquareWaveCylinder(float tc);
	void drawCylinderStripe(float tc);
	void drawCylinderStripes_Drifting(float tc);
	void drawCylinderStripe_Standing(float tc);
	void drawCylinderWedge(float angle, float w);
	void drawLocalCircle();

	void drawSineBeaters(float tc);
	void drawSquareWaveBeaters(float tc);
	void drawSineCircles(float tc);
	void drawSquareWaveCircles(float tc);
	void drawSavedStimMesh(float tc);
	void drawNoise(float tc);
	void drawSinePlane(float tc);
	void drawSineSphere(float tc);
	void drawUniformMesh();
	void drawUniformFullscreen();

	void readMesh();
	void readSavedStim();
	//void listFilesInDir();
	void recordStimInfo(TaskHandle taskHandle);


} ;
