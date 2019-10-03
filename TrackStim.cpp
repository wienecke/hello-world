#include "stdafx.h"
#include "TrackStim.h"

//const GLubyte* sssExtensions = glGetString(GL_EXTENSIONS);

/////////////// functions to set and write stimulus stuff from fly measurements //////////////////////////////

TrackStim::TrackStim()
{
	PI = 3.14159265359;
	ScreenUpdateRate = 100; // throw up at Hz

	windowWidth = 1140; // 1140whole window, PRO4500
	windowHeight = 912; // 912whole window, PRO4500
	bits = 32; // whole window, PRO4500

	monitor_stimulus_file = NULL;   // file to write for stimulus info
	paramfilename[0] = NULL;

	outputDir = "C:\\JL_stimuli\\";
	CreateDirectory(outputDir, NULL);

	TrackStim::setZeroTime();

	// initialize FMeas
	Sfull.x = Sfull.y = Sfull.theta = Sfull.xold = Sfull.yold = Sfull.thetaold = 0.0f;
	Sfull.vx = Sfull.vy = Sfull.vtheta = 0.0f;
	Sfull.tcurr = Sfull.told = Sfull.dt = 0.0f;
	Sfull.dx = Sfull.dy = Sfull.dtheta = 0.0f;

	Saccumulate = Sfull;

	boutInd = 0;

	// initialize stimulus stuff in case not read in
	numepochs = 1;
	epochchoose = 0;
	maxRunTime = 0;

	TrackStim::Stimulus = TrackStim::initializeStimulus();

	TrackStim::ViewPorts = TrackStim::initializeViewPositions();

	framenumber = 0;


	//perspELEVATION = asin(2.3/8) * 180 / PI; // ~16.7083 degrees. corresponds to rotation of stimulus screen
	//perspAZIMUTH = asin(1.6/8) * 180 / PI; // ~11.5370 degrees corresponds to rotation of stimulus screen
	//perspHEAD = 45; // corresponds to head tilt, though in matlab scripts, head tilt is 60 degrees

	//perspAZIMUTH = 13.2; // corresponds to frustum position so that the (x, 0, 0) point of the unrotated cylinder is rotated into the frustum. this minimizes distortion of the cylinder
	//perspELEVATION = 40.4; // corresponds to frustum position so that the (x, 0, 0) point of the unrotated cylinder is rotated into the frustum. this minimizes distortion of the cylinder
	//perspHEAD = 0; // no head tilt. head tilt in matlab code was to tilt ommatidia positions, not screen. here it doesn't make sense to have a head tilt since we consider everything from the perspective of columns that see the screen, and viewing axis, which is orthogonal to the plane of the screen

	// coordinates are: origin at fly's position, (1,0,0) in the optical axis of the microscope below the fly, (0,1,0) to the left of the fly's head from its own perspective, (0,0,1) above the fly's head from its own perspective
	////screen 1
	//perspAlpha = 17; // elevation drop from the microscope stage. in degrees
	//perspBeta = 11.7; // roll correction from away from (0,0,1). roll rotation is about (1,0,0). in degrees
	// screen 2
	perspAlpha = 0; // elevation drop from the microscope stage. in degrees
	perspBeta = 90; // roll correction from away from (0,0,1). roll rotation is about (1,0,0). in degrees

	perspGamma = atan((cos(perspAlpha * PI / 180) * sin(perspBeta * PI / 180)) / sin(perspAlpha * PI / 180)) * 180 / PI; // in degrees, 90 when alpha = 0 beta = 90 
	perspDelta = asin(cos(perspAlpha * PI / 180) * cos(perspBeta * PI / 180)) * 180 / PI; // in degrees, 0 (tiny number really) when alpha = 0 beta = 90

	perspAZIMUTH = perspGamma; // 90 degrees;
	perspELEVATION = perspDelta; // 0 degrees;
	perspHEAD = 0; // no head tilt. head tilt in matlab code was to tilt ommatidia positions, not screen. here it doesn't make sense to have a head tilt since we consider everything from the perspective of columns that see the screen, and viewing axis, which is orthogonal to the plane of the screen

	//// screen 1
	//cylinderAZIMUTH = perspAZIMUTH - 8; // offset empirically chosen to center the stimulus on the screen
	//cylinderELEVATION = perspELEVATION - 42; // offset empirically chosen to center the stimulus on the screen
	// screen 2
	cylinderAZIMUTH = perspAZIMUTH - 52; // offset empirically chosen to center the stimulus on the screen
	cylinderELEVATION = perspELEVATION + 15; // offset empirically chosen to center the stimulus on the screen

};

TrackStim::~TrackStim()
{
	// close files in use
	if (!(monitor_stimulus_file == NULL))
		fclose(monitor_stimulus_file);

	// make a time-stamped copy of output files
	char dateStr[9];
	char timeStr[9];
	_strdate(dateStr);
	_strtime(timeStr);
	char myDateTime[30];
	char cmd1[150];
	char cmd2[150];
	sprintf(myDateTime, "20%c%c_%c%c_%c%c_%c%c_%c%c_%c%c", dateStr[6], dateStr[7], dateStr[0], dateStr[1], dateStr[3], dateStr[4],
		timeStr[0], timeStr[1], timeStr[3], timeStr[4], timeStr[6], timeStr[7]);
	sprintf(cmd1, "copy .\\_stimulus_output.txt .\\_stimulus_output_%s.txt", myDateTime);
	sprintf(cmd2, "copy .\\_main_booleans.txt .\\_main_booleans_%s.txt", myDateTime);
	system(cmd1);
	system(cmd2);
	sprintf(cmd1, "copy .\\_stimulus_output.txt %s_stimulus_output_%s.txt", outputDir, myDateTime);
	sprintf(cmd2, "copy .\\_main_booleans.txt %s_main_booleans_%s.txt", outputDir, myDateTime);
	system(cmd1);
	system(cmd2);
};

TrackStim::STIMPARAMS TrackStim::initializeStimulus()
{
	TrackStim::STIMPARAMS Stim;   // declare it here...
	Stim.stimtype = 0;
	Stim.lum = .5;
	Stim.contrast = 1;
	Stim.duration = 3;   // in seconds
	Stim.transgain = 1;
	Stim.rotgain = 1;
	Stim.trans.mean = 0;
	Stim.trans.amp = 0;
	Stim.trans.per = 1;
	Stim.trans.phase = 0;
	Stim.rot.mean = 0;
	Stim.rot.amp = 0;
	Stim.rot.per = 1;
	Stim.rot.phase = 0;
	Stim.stimtrans.mean = 360;
	Stim.stimtrans.amp = 0;
	Stim.stimtrans.per = 1;
	Stim.stimtrans.phase = 0;
	Stim.stimrot.mean = 0;
	Stim.stimrot.amp = 0;
	Stim.stimrot.per = 1;
	Stim.stimrot.phase = 30;
	Stim.spacing = 30; //.010; //.1 or 10
	Stim.spacing2 = 5;
	Stim.density = .5;
	Stim.tau = .1;  // in seconds
	Stim.arenasize = 100; // in mm
	Stim.arenaheight = 30; // in mm
	Stim.tau2 = 0.3; // in seconds...
	Stim.randomize = 1;
	Stim.x_center = 0;
	Stim.y_center = 0;
	Stim.x_center2 = 0;
	Stim.y_center2 = 0;
	Stim.base_radius = 500;
	Stim.base_extent = 5000;
	Stim.USEFRUSTUM = 1;
	Stim.optogenetics = 0;

	return Stim;   // pass it back.
}

TrackStim::VIEWPOSITIONS TrackStim::initializeViewPositions()
{
	TrackStim::VIEWPOSITIONS temp;
	temp.x[0] = 0; temp.x[1] = 250;
	temp.y[0] = 0; temp.y[1] = 300;
	temp.w[0] = 200; temp.w[1] = 200;
	temp.h[0] = 200; temp.h[1] = 200;

	return temp;
};

int TrackStim::compareStrings(string stringyOne, string stringyTwo) {
	if (stringyTwo.compare(stringyOne)) {
		return 0;
	}
	else {
		return 1;
	}
	
}

void TrackStim::incrementFrameNumber()
{
	framenumber++;
};

void TrackStim::setZeroTime()
{
	QueryPerformanceCounter(&CounterBeginTime);  // time == 0
	QueryPerformanceFrequency(&CounterFreq);     // frequency of counter, needed to calculate time
};

float TrackStim::queryCurrTime()
{
	LARGE_INTEGER TimeNow;
	QueryPerformanceCounter(&TimeNow);
	return (float)(TimeNow.QuadPart - CounterBeginTime.QuadPart) / (float)CounterFreq.QuadPart;
};

void TrackStim::writeStim(TaskHandle taskHandle)
{
	uInt32 data;
	int         error = 0;
	char        errBuff[2048] = { '\0' };

	if (monitor_stimulus_file == NULL)
	{
		monitor_stimulus_file = fopen("_stimulus_output.txt", "w");
		fprintf(monitor_stimulus_file, "%s \n", paramfilename);
	}
	else
	{
		//fprintf(monitor_stimulus_file,"%8d %8.3f %2d %5.4f %5.4f %4.3f \n",framenumber,Sfull.tcurr,epochchoose,Sfull.x,Sfull.y,Sfull.theta);
		DAQmxErrChk(DAQmxReadCounterScalarU32(taskHandle, 1.0, &data, NULL));
		//char tmp[500];
		//sprintf(tmp,"data: %s\n",data);
		//MessageBox(NULL,tmp,"End of program",MB_YESNO|MB_ICONQUESTION)==IDNO;
		fprintf(monitor_stimulus_file, "%8d %8.3f %8d %8d %8.4f %8.4f %8.4f %8d\n", framenumber, Sfull.tcurr, boutInd, epochchoose, Sfull.x, Sfull.y, Sfull.theta, data);
		//fprintf(monitor_stimulus_file,"%8d %8.3f %2d %5.4f %5.4f %4.3f %d\n",framenumber,Sfull.tcurr,epochchoose,Sfull.x,Sfull.y,Sfull.theta,data);
		//fprintf(monitor_stimulus_file,"%8d %8.3f %2d %5.4f %5.4f %4.3f \n",framenumber,Sfull.tcurr,epochchoose,Sfull.x,Sfull.y,Sfull.theta);
		//fprintf(monitor_stimulus_file,"%8d %8.3f %2d %5.4f %5.4f %3.2f \n",framenumber,Sfull.tcurr,epochchoose,Scurr.x,Scurr.dx,Scurr.dt);
	};
Error:
	{
		puts("");
		if (DAQmxFailed(error))
			DAQmxGetExtendedErrorInfo(errBuff, 2048);
		//if( taskHandle!=0 ) {
			/*********************************************/
			// DAQmx Stop Code
			/*********************************************/
		//	DAQmxStopTask(taskHandle);
		//	DAQmxClearTask(taskHandle);
		//}
		if (DAQmxFailed(error))
		{
			char tmp[500];
			printf("DAQmx Error: %s\n", errBuff);
			sprintf(tmp, "DAQmx Error: %s\n", errBuff);
			fprintf(monitor_stimulus_file, "%8d %8.3f %2d %5.4f %5.4f %4.3f \n", framenumber, Sfull.tcurr, epochchoose, Sfull.x, Sfull.y, Sfull.theta);
			//MessageBox(NULL,tmp,"End of program",MB_YESNO|MB_ICONQUESTION)==IDNO;
			getchar();
			//return;
		}
	}
};

float TrackStim::mod(float x, float y)
{
	if (x >= 0)
		return x - y * (int)(x / y);
	else
		return x - y * (int)((x / y) - 1);
}

float TrackStim::round(float x)
{
	return floor(x + 0.5);
}

float TrackStim::sign(float in) {
	if (in > 0)
		return float(1);
	else if (in < 0)
		return float(-1);
	else if (in == 0)
		return 1;

}

//////////////////making real stimuli//////////////////////////

float TrackStim::getDLPColor(float DLPintensity, char channel)
{
	// returns the digital value to get that intensity in that channel

	//const static float gamma_gr = 2.0783; // new DLP bulb and chassis, measurements 131211
	//const static float scale_gr = 1.0155;
	//const static float gamma_bl = 2.0691;
	//const static float scale_bl = 1.0147;

	const static float gamma_re = 1; // red values are dummy values, but dummy values are theoretically correct for LightCrafter
	const static float scale_re = 1;
	const static float gamma_gr = 1; // green values are dummy values
	const static float scale_gr = 1;
	const static float gamma_bl = 1; // blue values are dummy values
	const static float scale_bl = 1;

	//const static float gamma_gr = 0.9958; // LightCrafter 4500, 6 bit depth, measurements 141110, second (frankenstein) unit
	//const static float scale_gr = 1.0088;
	//const static float gamma_bl = 1; // blue values are dummy values
	//const static float scale_bl = 1;

	//const static float gamma_gr = 0.9820; // LightCrafter 4500, 7 bit depth, measurements 140508
	//const static float scale_gr = 1.0427;
	//const static float gamma_bl = 1; // blue values are dummy values
	//const static float scale_bl = 1;

	float temp = 0;

	if (channel == 'R')
	{
		temp = pow(DLPintensity / scale_re, 1 / gamma_re);
	};

	if (channel == 'G')
	{
		temp = pow(DLPintensity / scale_gr, 1 / gamma_gr);
	};

	if (channel == 'B')
	{
		temp = pow(DLPintensity / scale_bl, 1 / gamma_bl);
	};

	//keep the output in the closed interval [0, 1]
	if (temp > 1)
	{
		temp = 1;
	}
	if (temp < 0)
	{
		temp = 0;
	}

	if (bitDepth == 8)
	{
		temp *= 255.0f / 255.0f; //maintain 8 bit depth 
	}
	else {

		static int nChannels;
		if (stimulusDataPath != (string) "NULL")
		{
			nChannels = int(Stimulus.density); // if we are presenting noise, check Stimulus.density to obtain nChannels

			if (nChannels == 0)
			{
				nChannels = 1;
			}
		}
		else
		{
			nChannels = 1; // default to one channel presentation
		}
		switch (nChannels) // Stimulus.density
		{
		case 1:
			temp *= 63.0f / 255.0f; // convert from 8 bit depth to 6 bit depth.
			break;
		case 2:
			temp *= 15.0f / 255.0f; // convert from 8 bit depth to 4 bit depth.
		case 3:
			temp *= 3.0f / 255.0f; // convert from 8 bit depth to 2 bit depth.
			break;
		}

		//temp = DLPintensity; // debug line to use if we want to turn off gamma correction
	}

	return temp;
};

void TrackStim::setColorTriples()
{
	float bgcol = Stimulus.lum * (1 - Stimulus.contrast);
	float fgcol = Stimulus.lum * (1 + Stimulus.contrast);

	backgroundtriple[0] = getDLPColor(bgcol, 'R');
	backgroundtriple[1] = getDLPColor(bgcol, 'G');
	backgroundtriple[2] = getDLPColor(bgcol, 'B');

	foregroundtriple[0] = getDLPColor(fgcol, 'R');
	foregroundtriple[1] = getDLPColor(fgcol, 'G');
	foregroundtriple[2] = getDLPColor(fgcol, 'B');

	blacktriple[0] = blacktriple[1] = blacktriple[2] = 0.0f;
	whitetriple[0] = whitetriple[1] = whitetriple[2] = 1.0f;
};

void TrackStim::setBackground() // sets background color
{
	// float color = pow((Stimulus.lum * (1 - Stimulus.contrast/2)) , (1/GammaCorrect));
	// glClearColor(color, color, color, 0.0f);
	//
	//float bgcol = Stimulus.lum * (1 - Stimulus.contrast/2);
	//glClearColor(getDLPColor(bgcol,'R'),getDLPColor(bgcol,'G'),getDLPColor(bgcol,'B'),0.0f);
	glClearColor(backgroundtriple[0], backgroundtriple[1], backgroundtriple[2], 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

};

void TrackStim::setBackgroundColor()
{
	// float color = pow((Stimulus.lum * (1 - Stimulus.contrast/2)) , (1/GammaCorrect));
	// glColor3f(color, color, color);				//

	//float bgcol = Stimulus.lum * (1 - Stimulus.contrast/2);
	//glColor3f(getDLPColor(bgcol,'R'),getDLPColor(bgcol,'G'),getDLPColor(bgcol,'B'));
	glColor3f(backgroundtriple[0], backgroundtriple[1], backgroundtriple[2]);
};

void TrackStim::setForeground() // sets foreground color
{
	// float color = pow((Stimulus.lum * (1 + Stimulus.contrast/2)) , (1/GammaCorrect));
	// glColor3f(color, color, color);

	//float col = Stimulus.lum * (1 + Stimulus.contrast/2);
	//glColor3f(getDLPColor(col,'R'),getDLPColor(col,'G'),getDLPColor(col,'B'));
	glColor3f(foregroundtriple[0], foregroundtriple[1], foregroundtriple[2]);
};

void TrackStim::readstring(FILE *f, char *string)
{
	do
	{
		fgets(string, 255, f);
	} while ((string[0] == '/') || (string[0] == '\n'));
	return;
}

void TrackStim::readNoise()   // arena size from elsewhere -- spacing in stimulus
{
	/* to get this to work, use code:
	paramfilename = GetFileName() ;
	for (int ii=0; ii<260; ii++) paramfilekeep[ii]=paramfilename[ii]; // seems very klugey, and below
	 -- and pass this paramfilekeep to this function
	*/

	float x;
	int numpoints;
	FILE *filein;
	char oneline[255];

	// add in dialog box here, use global variable to keep world name
	filein = fopen("E:/Code/trackball_final/Trackball/Debug/Data/noise.txt", "rt");				// File To Load World Data From

	readstring(filein, oneline);
	sscanf(oneline, "NUMPOINTS %d\n", &numpoints);

	for (int loop = 0; loop < min(numpoints, 14400); loop++)
	{
		readstring(filein, oneline);
		sscanf(oneline, "%f", &x);
		noiseValues[loop] = x;
	}
	fclose(filein);

};

float TrackStim::retrieveNextNoiseValue()
{
	static int count = -1;
	count++;
	return noiseValues[count / 2 % 14400]; // noise is sampled at 120Hz, 2 minutes
};
/*
void TrackStim::readStimulusData()
{
	stimulusData = new float[24 * 24 * 100000 * 3];

	// define memory dataspace
	hsize_t dims_out[2];
	dims_out[0] = 1;
	dims_out[1] = 24 * 24 * 100000 * 3;

	DataSpace memspace(2, dims_out);

	// open HDF5 file and open dataset "stimulus"
	H5File file(stimulusDataPath, H5F_ACC_RDONLY);
	DataSet dataset = file.openDataSet("stimulus"); // hard coded stimulus dataset

	// get dataset dataspace
	DataSpace dataspace = dataset.getSpace();

	hsize_t dims_in[3];
	dataspace.getSimpleExtentDims(dims_in, NULL);

	//hsize_t offset_in[3];
	//offset_in[0] = 0;
	//offset_in[1] = 0;
	//offset_in[2] = 0;

	//dataspace.selectHyperslab( H5S_SELECT_SET, dims_in, offset_in );

	dims_out[0] = 1;
	dims_out[1] = dims_in[0] * dims_in[1] * dims_in[2];

	hsize_t offset_out[2];
	offset_out[0] = 0;
	offset_out[1] = 0;

	memspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset_out);

	dataset.read(stimulusData, PredType::NATIVE_FLOAT, memspace, dataspace);
};
*/
void TrackStim::readViewPositions()
{
	FILE *filein;
	char oneline[255];

	if (GetFileAttributes("viewpositions.txt") != INVALID_FILE_ATTRIBUTES) // attempt to find world data in current directory
	{
		filein = fopen("viewpositions.txt", "rt");				// File To Load World Data From

		readstring(filein, oneline);
		sscanf(oneline, "%i %i", &ViewPorts.x[0], &ViewPorts.x[1]);
		readstring(filein, oneline);
		sscanf(oneline, "%i %i", &ViewPorts.y[0], &ViewPorts.y[1]);
		readstring(filein, oneline);
		sscanf(oneline, "%i %i", &ViewPorts.w[0], &ViewPorts.w[1]);
		readstring(filein, oneline);
		sscanf(oneline, "%i %i", &ViewPorts.h[0], &ViewPorts.h[1]);

		fclose(filein);
	}
	else if (GetFileAttributes("..\\viewpositions.txt") != INVALID_FILE_ATTRIBUTES) // attempt to find world data in parent directory
	{
		filein = fopen("..\\viewpositions.txt", "rt");				// File To Load World Data From

		readstring(filein, oneline);
		sscanf(oneline, "%i %i", &ViewPorts.x[0], &ViewPorts.x[1]);
		readstring(filein, oneline);
		sscanf(oneline, "%i %i", &ViewPorts.y[0], &ViewPorts.y[1]);
		readstring(filein, oneline);
		sscanf(oneline, "%i %i", &ViewPorts.w[0], &ViewPorts.w[1]);
		readstring(filein, oneline);
		sscanf(oneline, "%i %i", &ViewPorts.h[0], &ViewPorts.h[1]);

		fclose(filein);
	}
	else // otherwise load default world data
	{
		ViewPorts = initializeViewPositions();
	}
};

void TrackStim::writeViewPositions()
{
	FILE *fileout;

	fileout = fopen("viewpositions.txt", "w");
	fprintf(fileout, "%i %i \n", ViewPorts.x[0], ViewPorts.x[1]);
	fprintf(fileout, "%i %i \n", ViewPorts.y[0], ViewPorts.y[1]);
	fprintf(fileout, "%i %i \n", ViewPorts.w[0], ViewPorts.w[1]);
	fprintf(fileout, "%i %i \n", ViewPorts.h[0], ViewPorts.h[1]);
	fclose(fileout);
};

void TrackStim::readSavedStim() // Read a saved stim from binary file
{

	FILE * pFile;
	long lSize;
	size_t result;

	const char* prefixy = "..\\";
	const char* suffixy = filenames[fileIndex];

	//char *filStr = (char*)malloc(strlen(prefixy) + 1 + 4);
	char *filStr = NULL;
	filStr = new char[260];
	strcpy(filStr, prefixy);
	strcat(filStr, suffixy);
	
	pFile = fopen(filStr, "rb");
	if (pFile == NULL) { fputs("File error", stderr); exit(1); }
	
	char line[256];

	strcpy(line, filStr);
	char *subStringOne = strtok(line, "\_");
	subStringOne = strtok(NULL, "\_");
	char *subStringTwo = strtok(NULL, "\_");
	char *subStringThree = strtok(NULL, "\_");
	char *subStringFour = strtok(NULL, "\_");
	char *subStringFive = strtok(NULL, "\_");
	char *subStringSix = strtok(NULL, "\_");

	
	texType[fileIndex] = atoi(subStringOne);
	texXres[fileIndex] = atoi(subStringTwo);
	texYres[fileIndex] = atoi(subStringThree);
	texTres[fileIndex] = atoi(subStringFour);
	texSpacing[fileIndex] = atoi(subStringFive);
	texDir[fileIndex] = atoi(subStringSix);


	delete[] filStr;
	filStr = nullptr;
	// obtain the file size
	fseek(pFile, 0, SEEK_END);
	lSize = ftell(pFile) / sizeof(float);
	rewind(pFile);

	// allocate memory to contain the whole file 
	buffer = new float[lSize];
	//buffer = (float*)malloc(sizeof(float)*lSize);
	if (buffer == NULL) { fputs("Memory error", stderr); exit(2); }

	//copy the file into the buffer
	result = fread(buffer, sizeof(float), lSize, pFile);
	if (result != lSize) { fputs("Reading error", stderr); exit(3); }

	// terminate 
	fclose(pFile);

	// copy buffer to every third element of movie (populate RG values with zeros)
	movie[fileIndex] = new float[lSize * colorChannels];
	//movie[fileIndex] = (float*)malloc(sizeof(float)*lSize * colorChannels);
	for (int i = 0; i < texXres[fileIndex] * texYres[fileIndex] * texTres[fileIndex]; i++) {
		movie[fileIndex][3 * i] = 0;
		movie[fileIndex][3 * i + 1] = buffer[i];
		movie[fileIndex][3 * i + 2] = 0;
	}

	delete[] buffer;
	buffer = nullptr;

};

void TrackStim::readMesh() //Read a saved mesh for geometric/luminance correction

{
	int itx, ity, itt;

	int ix, iy;
	double x, y, u, v, br;

	FILE *filein;
	filein = fopen("..\\mesh_2_.txt", "rt");

	//Get mesh dimensions
	if (fscanf(filein, "%d %d", &ix, &iy) != 2) {
		fprintf(stderr, "Failed to read the mesh dimensions\n");
		fclose(filein);
	}

	int meshlength;
	meshxres = ix;
	meshyres = iy;
	meshlength = meshxres * meshyres;

	// Create new mesh
	meshx = new double[meshlength];
	meshy = new double[meshlength];
	meshu = new double[meshxres];
	meshv = new double[meshyres];
	meshi = new double[meshlength];
	//meshx = (double*)malloc(sizeof(double)*meshlength);
	//meshy = (double*)malloc(sizeof(double)*meshlength);
	//meshu = (double*)malloc(sizeof(double)*meshxres);
	//meshv = (double*)malloc(sizeof(double)*meshyres);
	//meshi = (double*)malloc(sizeof(double)*meshlength);

	for (int i = 0; i < meshlength; i++) {
		if (fscanf(filein, "%lf %lf %lf", &x, &y, &br) != 3) {
			fprintf(stderr, "Unexpected end of mesh file encountered\n");
			fclose(filein);
		}
		meshx[i] = y;// *1.6;
		meshy[i] = x;
		//meshu[i] = u;
		//meshv[i] = v;
		meshi[i] = br;
		if (br > -1)
			meshi[i] = 1;
	}

	int foo = 1;

	for (int iii = 0; iii < meshxres; iii++) {
		meshu[iii] = (double)iii / (meshxres - 1);
	}
	for (int iii = 0; iii < meshyres; iii++) {
		meshv[iii] = (double)iii / (meshyres - 1);
	}

};

/*
void TrackStim::listFilesInDir(void)
{
	FILE *fp;
	DIR           *d;
	struct dirent *dir;
	//char *filenames[100];
	int a;
	int i = 0;
	d = opendir("..\\savedStimuli");
	if (d)
	{
		while ((dir = readdir(d)) != NULL)
		{
			filenames[i] = (char*)malloc(strlen(dir->d_name) + 1);
			strcpy(filenames[i], dir->d_name);
			i++;
		}
		closedir(d);
	}
}
*/



bool TrackStim::drawScene() // executes the drawing set by setStimType
{

	glEnable(GL_SCISSOR_TEST);
	glViewport(0, 0, windowWidth, windowHeight);
	glScissor(0, 0, windowWidth, windowHeight);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // set whole screen to be black; scissor testing makes backgrounds good in viewports
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear The Screen And The Depth Buffer
	glLoadIdentity();									// Reset The View
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // still on whole window

	Sfull.tcurr = queryCurrTime();

	Stimulus.base_radius = 200.0f; // hard-coded radius falls in the frustum depth range, dynamic for each ViewPort

	glViewport(ViewPorts.x[0], ViewPorts.y[0], ViewPorts.w[0], ViewPorts.h[0]);//glViewport(0, 0, 1140, 912);

	glScissor(ViewPorts.x[0], ViewPorts.y[0], ViewPorts.w[0], ViewPorts.h[0]);//(0,0,1140,912);

	glClear(GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();// Reset The View

	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_FALSE);

	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);// Select The Projection Matrix
	glLoadIdentity();

	if (!Stimulus.USEFRUSTUM)
	{
		gluPerspective(90.0f, 2 * (float)ViewPorts.w[0] / (float)ViewPorts.h[0], 0.1f, Stimulus.arenasize); // for LightCrafter 4500 pixels

	}
	else
	{
		switch (Stimulus.stimtype)
		{
		case 165:
			glOrtho(-50.0f, 50.0f, -50.0f, 50.0f, 0.04f, 0.06f);
			//glOrtho(-5.0f, 5.0f, -5.0f, 5.0f, 0.04f, 0.06f);
			//glOrtho(-100.0f, 100.0f, -100.0f, 100.0f, 0.04f, 0.06f);
			break;
		case 70:
		case 50:
			glFrustum(110.0f, 20.0f, -51.0f, 93.0f, 50.0f, 1000.0f); // modified from "screen 2 on InFocus DLP" settings, 170624. screen is now left-right inverted.
			//glFrustum(-500.0f, 500.0f, -500.0f, 500.0f, 50.0f, 1000.0f); // modified from "screen 2 on InFocus DLP" settings, 170624. screen is now left-right inverted.
			//glFrustum(25.0f, 20.0f, 15.0f, 20.0f, 50.0f, 1000.0f);
			//gluPerspective(160.0f, 1.25, 49.9f, 50.1); // for LightCrafter 4500 pixels
			break;
		case 11:
		case 46:
		case 64:
		case 65:
		case 66:
		case 67:
		case 71:
		case 68:
			//glOrtho(-1.0f, 1.0f, -1.6f, 1.6f, 0.04f, 0.06f);
			glOrtho(-1.0f, 1.0f, -1.56f, 1.56f, 0.04f, 0.06f);
			//glOrtho(-1.0f, 1.0f, -1.6f, 1.6f, 0.04f, 0.06f); // this is the one
			//glOrtho(-1.25f, 1.25f, -2.0f, 2.0f, 0.04f, 0.06f);
			//glOrtho(-.50f, .50f, -.32f, .32f, 0.04f, 0.06f);
			//glOrtho(-1.29f, 1.29f, -2.07f, 2.07f, 0.04f, 0.06f); 
			// NO y*1.6,x,br=1/-1 (ellipse, vertical overfills window fills screen, hor full underfills window and screen)
			// NO y*1.6,x,br=1 same as above y*1.6,x,br=1/-1 
			// YES BEST y,x,br=1/-1 (circle, vertical jagged underfills window and screen, hor full underfills window and screen)
			// NO x,y,br=1/-1 (circle, vertical full fills window and screen, hor jagged underfills window and screen)
			// NO x*1.6,y,br=1/-1 (ellipse, vertical overfills window fills screen, hor jagged underfills window and screen)
			// NO x,y,br=1 (circle, vertical full fills window and screen, hor full underfills window and screen)
			//glOrtho(-2.07f, 2.07f, -1.29f, 1.29f, 0.04f, 0.06f); 
			// NO y*1.6,x,br=1/-1 (ellipse, vertical jagged underfills window and screen, hor full fills window overfills screen crescent)
			// NO y*1.6,x,br=1 (ellipse, vertical full fills window and screen, hor full fills window overfills screen)
			// NO y,x,br=1/-1 (ellipse, vertical jagged underfills window and screen, hor full fills window overfills screen)
			// NO x,y,br=1/-1 (ellipse, vertical full underfills window and screen, hor jagged underfills window and screen)
			// NO x*1.6,y,br=1/-1 (ellipse, vertical full fills window and screen, hor jagged underfills window and screen)
			//glOrtho(-1.29f, 1.29f, -1.29f, 1.29f, 0.04f, 0.06f); 
			// NO/MAYBE y*1.6,x,br=1/-1 (circle, vertical overfills window fills screen, hor full fills window overfills screen crescent)
			// NO/MAYBE y*1.6,x,br=1 same as above y*1.6,x,br=1/-1
			// NO y,x,br=1/-1 (ellipse, vertical jagged underfills window and screen, hor full fills window overfills screen)
			// NO x,y,br=1/-1 (ellipse, vertical full fills window and screen, hor jagged underfills window and screen)
			// NO x*1.6,y,br=1/-1 (circle, vertical overfills window fills screen, hor jagged underfills window and screen)
			//glOrtho(-2.07f, 2.07f, -2.07f, 2.07f, 0.04f, 0.06f); 
			// YES SECOND BEST y*1.6,x,br=1/-1 (circle, vertical jagged underfills window and screen, hor full underfills window and screen too MUCH, BUT MAYBE SAME AS ABOVE YES BEST FOR PROJECTOR EVEN FURTHER BACK?)
			// NO y*1.6,x,br=1 (circle, vertical full fills window and screen, hor full underfills window and screen)
			// NO/MAYBE y,x,br=1/-1 (ellipse, vertical jagged underfills window and screen, hor full underfills window and screen too MUCH, BUT MAYBE SAME AS ABOVE FOR PROJECTOR EVEN FURTHER BACK?)
			// NO x,y,br=1/-1 (ellipse, vertical full underfills window and screen, hor jagged underfills window and screen)
			// NO x*1.6,y,br=1/-1 (circle, vertical full fills window and screen, hor jagged underfills window and screen)
			//glOrtho(-100.0f, 100.0f, -100.0f, 100.0f, 0.0f, 500.0f); // modified from "screen 2 on InFocus DLP" settings, 170624. screen is now left-right inverted.
			// PRO4500
			//glFrustum(84.3f, 2.3f, -65.6f, 65.6f, viewInd * 1000.0f + 50.0f, viewInd * 1000.0f + 1000.0f); // modified from "screen 2 on LightCrafter 4500" settings, 170226. screen is still left-right inverted. Screen size is now the full screen for the PRO4500, which is 131.2 mm x 82 mm
			//glFrustum(41.0f, -41.0f, -65.6f, 65.6f, viewInd * 1000.0f + 50.0f, viewInd * 1000.0f + 1000.0f); // modified from "screen 2 on LightCrafter 4500" settings, 170226. screen is still left-right inverted. Screen size is now the full screen for the PRO4500, which is 131.2 mm x 82 mm					break;

			break;
		};
	}

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();  // start at origin
	/*
	float lookAtVertex[3];
	// this lookAt vector is calculated by rotating the "initial" lookAt vector (1,0,0) CW in elevation and left (CCW) in azimuth
	lookAtVertex[0] = cos(perspELEVATION * PI / 180) * cos(perspAZIMUTH * PI / 180); //0
	lookAtVertex[1] = cos(perspELEVATION * PI / 180) * sin(perspAZIMUTH * PI / 180); //1
	lookAtVertex[2] = sin(perspELEVATION * PI / 180); //0
	//// this lookAt vector is calculated by rotating the "initial" lookAt vector (0,0,1) CCW in elevation and left (CW) in roll
	//lookAtVertex[0] = sin( perspAlpha * PI / 180 );
	//lookAtVertex[1] = cos( perspAlpha * PI / 180 ) * sin( perspBeta * PI / 180 );
	//lookAtVertex[2] = cos( perspAlpha * PI / 180 ) * cos( perspBeta * PI / 180 );

	float lookAtUp[3];
	// this up vector is calculated by rotating the "initial" up vector (0,0,1) CW in elevation and left (CCW) in azimuth. corresponds to the fly lookAtUp
	lookAtUp[0] = -sin(perspELEVATION * PI / 180) * cos(perspAZIMUTH * PI / 180); //0
	lookAtUp[1] = -sin(perspELEVATION * PI / 180) * sin(perspAZIMUTH * PI / 180); //0
	lookAtUp[2] = cos(perspELEVATION * PI / 180); //1
	////this up vector is calculated by rotating the "initial" up vector (-1,0,0) CCW in elevation and left (CW) in roll. corresponds to the physical screen lookAtUp
	//lookAtUp[0] = -cos( perspAlpha * PI / 180 );
	//lookAtUp[1] = sin( perspAlpha * PI / 180 ) * sin( perspBeta * PI / 180 );
	//lookAtUp[2] = sin ( perspAlpha * PI / 180 ) * cos( perspBeta * PI / 180 );

	gluLookAt(0.0f, 0.0f, 0.0f, lookAtVertex[0], lookAtVertex[1], lookAtVertex[2], lookAtUp[0], lookAtUp[1], lookAtUp[2]); //0,1,0  0,0,1
	*/
	gluLookAt(0.0f, 0.0f, 0.0f, 0, 1, 0, 0, 0, 1); //0,1,0  0,0,1

	setColorTriples();
	setBackground();
	setForeground();

	// all set, now draw the stimulus
	switch (Stimulus.stimtype)
	{
	case 11:
	case 46:
		drawSavedStimMesh(Sfull.tcurr);
		break;
	case 50:
		drawCylinderStripe(Sfull.tcurr - epochtimezero);
		break;
	case 51:
		drawSquareWaveCylinder(Sfull.tcurr - epochtimezero);
		break;
	case 56:
		drawSineCylinder_Standing_ContrastModulated(Sfull.tcurr);
		break;
	case 60:
		drawSineCylinder_Standing();
		break;
	case 61:
		drawCylinderStripe_Standing(Sfull.tcurr);
		break;
	case 62:
		drawCylinderStripes_Drifting(Sfull.tcurr);
		break;
	case 64:
		drawSineBeaters(Sfull.tcurr);
		break;
	case 65:
		drawSquareWaveBeaters(Sfull.tcurr);
		break;
	case 66:
		drawSineCircles(Sfull.tcurr);
		break;
	case 67:
		drawSquareWaveCircles(Sfull.tcurr);
		break;
	case 68:
		drawLocalCircle();
		break;
	case 70:
		drawSineCylinder(Sfull.tcurr);
		break;
	case 71:
		drawUniformFullscreen();
		break;
	case 165:
		drawSinePlane(Sfull.tcurr);
		break;
	case 166:
		drawSineSphere(Sfull.tcurr);
		break;
	case 169:
		drawUniformMesh();
		break;
	};

	glDisable(GL_SCISSOR_TEST); // for some reason, want to do this before swap buffer operation

	return TRUE;
};

void TrackStim::drawUniformFullscreen()
{
	const int domainy = 1024;
	static float imagedata[domainy * 3];

	for (int ii = 0; ii < domainy; ii++)
	{
		imagedata[ii * 3] = 1;
		imagedata[ii * 3 + 1] = 1;
		imagedata[ii * 3 + 2] = 1;
	}

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, domainy, 0, GL_RGB, GL_FLOAT, imagedata);

	glPushMatrix(); // just to be sure...

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);

	glBegin(GL_QUADS);

	glTexCoord1f(0.5);
	glVertex3f(50, 0.05f, -50);

	glTexCoord1f(0.5);
	glVertex3f(-50, 0.05f, -50);

	glTexCoord1f(0.5);
	glVertex3f(-50, 0.05f, 50);

	glTexCoord1f(0.5);
	glVertex3f(50, 0.05f, 50);


	glEnd();

	glDisable(GL_TEXTURE_1D);
	glPopMatrix();

};

void TrackStim::drawSinePlane(float tc)
{
	const int domainy = 1024;
	static float imagedata[domainy * 3];

	float added = PI / 2;	// add the stimulus OL
							// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
							// we want to choose this value to be the temporal frequency here.
							//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																		// do tc-epochtimezero so that the start is at the same position every time
	//added = fmod(double(added - PI / 2), double(Stimulus.spacing * PI / 180)); // keep added to within a single period. this should prevent the seam from appearing. 130628 leongjcs
	Sfull.x = added - PI / 2; // save the phase

	float per = Stimulus.spacing;
	int scale;
	scale = 1;//200/per;

	for (int ii = 0; ii < domainy; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin((2 * PI * ii / domainy) - added)), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin((2 * PI * ii / domainy) - added)), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin((2 * PI * ii / domainy) - added)), 'B');
	}

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, domainy, 0, GL_RGB, GL_FLOAT, imagedata);

	float wedgeWidthDeg = 5;
	float wedgewidth = wedgeWidthDeg * PI / 180; // width of each wedge
	float phase = Stimulus.stimrot.phase * PI / 180;
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	float heighty;
	float angle;
	float angleTwo;

	//Sfull.theta = phase;
	glPushMatrix(); // just to be sure...

	//transformToScopeCoords();

	// now (then), put it into the right pitch coordinates...
	//glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

	// now (first), roll it to get right perp angle...
	float rollang = Stimulus.stimrot.mean;
	Sfull.theta = rollang;
	glRotatef(rollang, 0.0f, 1.0f, 0.0f);
	//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);

	float perspe = 10;

	glBegin(GL_QUADS);
	for (int ii = 0; ii < 100; ii++) // through each wedge...
	{
		angle = (ii - 50);// + (added * 180 / PI);
		angleTwo = ((ii - 50) + 1);// +(added * 180 / PI);

		glTexCoord1f(ii / per);
		glVertex3f(angle, 0.05f, -50);

		glTexCoord1f((ii + 1) / per);
		glVertex3f(angleTwo, 0.05f, -50);

		glTexCoord1f((ii + 1) / per);
		glVertex3f(angleTwo, 0.05f, 50);

		glTexCoord1f(ii / per);
		glVertex3f(angle, 0.05f, 50);
	}

	glEnd();

	glDisable(GL_TEXTURE_1D);
	glPopMatrix();

};

void TrackStim::drawUniformMesh()
{

	glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texYres[epochchoose], texXres[epochchoose], 0, GL_RGB, GL_FLOAT, movie[epochchoose]);

	glPushMatrix(); // just to be sure...

	//static bool firstie = TRUE;
	//static int displaylist;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);

	if (firstie)//(foist)
	{
		glDeleteLists(displaylist, 1);
		displaylist = glGenLists(1);
		glNewList(displaylist, GL_COMPILE);
		glBegin(GL_QUADS);

		int nn = meshyres;
		for (int ii = 0; ii < meshxres - 1; ii++) // through each wedge...
		{
			for (int jj = 0; jj < meshyres - 1; jj++)
			{
				int ki = nn*ii + jj;
				int kj = nn*(ii + 1) + jj;
				int kk = nn*(ii + 1) + (jj + 1);
				int kl = nn*ii + (jj + 1);

				if (meshi[ki] && meshi[kj] && meshi[kk] && meshi[kl])
				{

					glTexCoord2f(meshv[jj], meshu[ii]);
					glVertex3f(meshx[ki], 0.05f, meshy[ki]);

					glTexCoord2f(meshv[jj], meshu[ii + 1]);
					glVertex3f(meshx[kj], 0.05f, meshy[kj]);

					glTexCoord2f(meshv[jj + 1], meshu[ii + 1]);
					glVertex3f(meshx[kk], 0.05f, meshy[kk]);

					glTexCoord2f(meshv[jj + 1], meshu[ii]);
					glVertex3f(meshx[kl], 0.05f, meshy[kl]);
				}
			}
		}
		glEnd();
		glEndList();
		firstie = FALSE;
	}
	else
	{
		glCallList(displaylist);
	}

	glDisable(GL_TEXTURE_2D);
	glPopMatrix();

};

void TrackStim::drawSquareWaveBeaters(float tc)
{
	const int domainy = 512;
	static float imagedata[domainy * 3];

	float added = PI / 2;	// add the stimulus OL
							// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
							// we want to choose this value to be the temporal frequency here.
							//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																		// do tc-epochtimezero so that the start is at the same position every time
	Sfull.x = added - PI / 2; // save the phase

	float per = Stimulus.spacing / 100;
	int scale;
	scale = 8;//200/per;

	for (int ii = 0; ii < domainy; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sign(sin(scale*((2 * PI * ii / domainy) - added)))), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sign(sin(scale * ((2 * PI * ii / domainy) - added)))), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sign(sin(scale * ((2 * PI * ii / domainy) - added)))), 'B');
	}

	glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, domainy, 1, 0, GL_RGB, GL_FLOAT, imagedata);

	glPushMatrix(); // just to be sure...

	//static bool first = TRUE;
	//static int displaylist;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);

	if (firstie)//(foist)
	{
		glDeleteLists(displaylist, 1);
		displaylist = glGenLists(1);
		glNewList(displaylist, GL_COMPILE);
		glBegin(GL_QUADS);

		int nn = meshyres;
		for (int ii = 0; ii < meshxres - 1; ii++) // through each wedge...
		{
			for (int jj = 0; jj < meshyres - 1; jj++)
			{
				int ki = nn*ii + jj;
				int kj = nn*(ii + 1) + jj;
				int kk = nn*(ii + 1) + (jj + 1);
				int kl = nn*ii + (jj + 1);

				if (meshi[ki] && meshi[kj] && meshi[kk] && meshi[kl])
				{

					glTexCoord2f(meshv[jj], meshu[ii]);
					glVertex3f(meshx[ki], 0.05f, meshy[ki]);

					glTexCoord2f(meshv[jj], meshu[ii + 1]);
					glVertex3f(meshx[kj], 0.05f, meshy[kj]);

					glTexCoord2f(meshv[jj + 1], meshu[ii + 1]);
					glVertex3f(meshx[kk], 0.05f, meshy[kk]);

					glTexCoord2f(meshv[jj + 1], meshu[ii]);
					glVertex3f(meshx[kl], 0.05f, meshy[kl]);
				}
			}
		}
		glEnd();
		glEndList();
		firstie = FALSE;
	}
	else
	{
		glCallList(displaylist);
	}

	glDisable(GL_TEXTURE_2D);
	glPopMatrix();

};

void TrackStim::drawSineBeaters(float tc)
{

	const int domainy = 512;
	static float imagedata[domainy * 3];

	float added = PI / 2;	// add the stimulus OL
							// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
							// we want to choose this value to be the temporal frequency here.
							//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																		// do tc-epochtimezero so that the start is at the same position every time
	Sfull.x = added - PI / 2; // save the phase

	float per = Stimulus.spacing / 100;
	int scale;
	scale = 8;//200/per;

	for (int ii = 0; ii < domainy; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(scale*((2 * PI * ii / domainy) - added))), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(scale * ((2 * PI * ii / domainy) - added))), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(scale * ((2 * PI * ii / domainy) - added))), 'B');
	}

	glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, domainy, 1, 0, GL_RGB, GL_FLOAT, imagedata);

	glPushMatrix(); // just to be sure...

	//static bool first = TRUE;
	//static int displaylist;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);

	if (firstie)//(foist)
	{
		glDeleteLists(displaylist, 1);
		displaylist = glGenLists(1);
		glNewList(displaylist, GL_COMPILE);
		glBegin(GL_QUADS);

		int nn = meshyres;
		for (int ii = 0; ii < meshxres - 1; ii++) // through each wedge...
		{
			for (int jj = 0; jj < meshyres - 1; jj++)
			{
				int ki = nn*ii + jj;
				int kj = nn*(ii + 1) + jj;
				int kk = nn*(ii + 1) + (jj + 1);
				int kl = nn*ii + (jj + 1);

				if (meshi[ki] && meshi[kj] && meshi[kk] && meshi[kl])
				{

					glTexCoord2f(meshv[jj], meshu[ii]);
					glVertex3f(meshx[ki], 0.05f, meshy[ki]);

					glTexCoord2f(meshv[jj], meshu[ii + 1]);
					glVertex3f(meshx[kj], 0.05f, meshy[kj]);

					glTexCoord2f(meshv[jj + 1], meshu[ii + 1]);
					glVertex3f(meshx[kk], 0.05f, meshy[kk]);

					glTexCoord2f(meshv[jj + 1], meshu[ii]);
					glVertex3f(meshx[kl], 0.05f, meshy[kl]);
				}
			}
		}
		glEnd();
		glEndList();
		firstie = FALSE;
	}
	else
	{
		glCallList(displaylist);
	}

	glDisable(GL_TEXTURE_2D);
	glPopMatrix();

};

void TrackStim::drawSquareWaveCircles(float tc)
{
	const int domainy = 128;
	static float imagedata[domainy * 3];

	float added = PI / 2;	// add the stimulus OL
							// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
							// we want to choose this value to be the temporal frequency here.
							//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																		// do tc-epochtimezero so that the start is at the same position every time
	Sfull.x = added - PI / 2; // save the phase

	float per = Stimulus.spacing / 100;
	int scale;
	scale = Stimulus.spacing / PI;//200/per;

	for (int ii = 0; ii < domainy; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sign(sin(scale*((2 * PI * ii / domainy) - added)))), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sign(sin(scale * ((2 * PI * ii / domainy) - added)))), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sign(sin(scale * ((2 * PI * ii / domainy) - added)))), 'B');
	}

	glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, domainy, 0, GL_RGB, GL_FLOAT, imagedata);

	glPushMatrix(); // just to be sure...

					//static bool first = TRUE;
					//static int displaylist;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);

	if (firstie)//(foist)
	{
		glDeleteLists(displaylist, 1);
		displaylist = glGenLists(1);
		glNewList(displaylist, GL_COMPILE);
		glBegin(GL_QUADS);

		int nn = meshyres;
		for (int ii = 0; ii < meshxres - 1; ii++) // through each wedge...
		{
			for (int jj = 0; jj < meshyres - 1; jj++)
			{
				int ki = nn*ii + jj;
				int kj = nn*(ii + 1) + jj;
				int kk = nn*(ii + 1) + (jj + 1);
				int kl = nn*ii + (jj + 1);

				if (meshi[ki] && meshi[kj] && meshi[kk] && meshi[kl])
				{

					glTexCoord2f(meshv[jj], meshu[ii]);
					glVertex3f(meshx[ki], 0.05f, meshy[ki]);

					glTexCoord2f(meshv[jj], meshu[ii + 1]);
					glVertex3f(meshx[kj], 0.05f, meshy[kj]);

					glTexCoord2f(meshv[jj + 1], meshu[ii + 1]);
					glVertex3f(meshx[kk], 0.05f, meshy[kk]);

					glTexCoord2f(meshv[jj + 1], meshu[ii]);
					glVertex3f(meshx[kl], 0.05f, meshy[kl]);
				}
			}
		}
		glEnd();
		glEndList();
		firstie = FALSE;
	}
	else
	{
		glCallList(displaylist);
	}

	glDisable(GL_TEXTURE_2D);
	glPopMatrix();

};

void TrackStim::drawSineCircles(float tc)
{

	const int domainy = 128;
	static float imagedata[domainy * 3];

	float added = PI / 2;	// add the stimulus OL
							// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
							// we want to choose this value to be the temporal frequency here.
							//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																		// do tc-epochtimezero so that the start is at the same position every time
	Sfull.x = added - PI / 2; // save the phase

	float per = Stimulus.spacing / 100;
	int scale;
	scale = Stimulus.spacing/PI;//200/per;

	for (int ii = 0; ii < domainy; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(scale*((2 * PI * ii / domainy) - added))), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(scale * ((2 * PI * ii / domainy) - added))), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(scale * ((2 * PI * ii / domainy) - added))), 'B');
	}

	glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, domainy, 0, GL_RGB, GL_FLOAT, imagedata);

	glPushMatrix(); // just to be sure...

					//static bool first = TRUE;
					//static int displaylist;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);

	if (firstie)//(foist)
	{
		glDeleteLists(displaylist, 1);
		displaylist = glGenLists(1);
		glNewList(displaylist, GL_COMPILE);
		glBegin(GL_QUADS);

		int nn = meshyres;
		for (int ii = 0; ii < meshxres - 1; ii++) // through each wedge...
		{
			for (int jj = 0; jj < meshyres - 1; jj++)
			{
				int ki = nn*ii + jj;
				int kj = nn*(ii + 1) + jj;
				int kk = nn*(ii + 1) + (jj + 1);
				int kl = nn*ii + (jj + 1);

				if (meshi[ki] && meshi[kj] && meshi[kk] && meshi[kl])
				{

					glTexCoord2f(meshv[jj], meshu[ii]);
					glVertex3f(meshx[ki], 0.05f, meshy[ki]);

					glTexCoord2f(meshv[jj], meshu[ii + 1]);
					glVertex3f(meshx[kj], 0.05f, meshy[kj]);

					glTexCoord2f(meshv[jj + 1], meshu[ii + 1]);
					glVertex3f(meshx[kk], 0.05f, meshy[kk]);

					glTexCoord2f(meshv[jj + 1], meshu[ii]);
					glVertex3f(meshx[kl], 0.05f, meshy[kl]);
				}
			}
		}
		glEnd();
		glEndList();
		firstie = FALSE;
	}
	else
	{
		glCallList(displaylist);
	}

	glDisable(GL_TEXTURE_2D);
	glPopMatrix();

};

void TrackStim::drawSavedStimMesh(float tc)
{

	glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

	int rollang = Stimulus.stimrot.mean;
	const int n = texXres[epochchoose] * texYres[epochchoose] * colorChannels;

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texYres[epochchoose], texXres[epochchoose], 0, GL_RGB, GL_FLOAT, movie[epochchoose] + n*frameIndex);

	frameIndex = (frameIndex + 1) % texTres[epochchoose];

	glPushMatrix(); // just to be sure...

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);

	if (firstie)
	{
		glDeleteLists(displaylist, 1);
		displaylist = glGenLists(1);
		glNewList(displaylist, GL_COMPILE);
		glBegin(GL_QUADS);

		int nn = meshyres;
		for (int ii = 0; ii < meshxres - 1; ii++)
		{
			for (int jj = 0; jj < meshyres - 1; jj++)
			{
				int ki = nn*ii + jj;
				int kj = nn*(ii + 1) + jj;
				int kk = nn*(ii + 1) + (jj + 1);
				int kl = nn*ii + (jj + 1);

				if (meshi[ki] && meshi[kj] && meshi[kk] && meshi[kl])
				{

					glTexCoord2f(meshv[jj], meshu[ii]);
					glVertex3f(meshx[ki], 0.05f, meshy[ki]);

					glTexCoord2f(meshv[jj], meshu[ii + 1]);
					glVertex3f(meshx[kj], 0.05f, meshy[kj]);

					glTexCoord2f(meshv[jj + 1], meshu[ii + 1]);
					glVertex3f(meshx[kk], 0.05f, meshy[kk]);

					glTexCoord2f(meshv[jj + 1], meshu[ii]);
					glVertex3f(meshx[kl], 0.05f, meshy[kl]);
				}
			}
		}
		glEnd();
		glEndList();
		firstie = FALSE;
	}
	else
	{
		glCallList(displaylist);
	}

	glDisable(GL_TEXTURE_2D);
	glPopMatrix();
};

void TrackStim::drawSineSphere(float tc)
{
	const int domainy = 128;
	static float imagedata[domainy * 3];

	float added = PI / 2;	// add the stimulus OL
							// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
							// we want to choose this value to be the temporal frequency here.
							//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																		// do tc-epochtimezero so that the start is at the same position every time
	added = fmod(double(added - PI / 2), double(Stimulus.spacing * PI / 180)); // keep added to within a single period. this should prevent the seam from appearing. 130628 leongjcs
	Sfull.x = added - PI / 2; // save the phase

	for (int ii = 0; ii < domainy; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(2 * PI * ii / domainy)), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(2 * PI * ii / domainy)), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(2 * PI * ii / domainy)), 'B');
	}

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, domainy, 0, GL_RGB, GL_FLOAT, imagedata);

	float wedgeWidthDeg = 5;
	float wedgewidth = wedgeWidthDeg * PI / 180; // width of each wedge
	float phase = Stimulus.stimrot.phase * PI / 180;
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	float heighty;
	float angle;
	float per = Stimulus.spacing;

	//Sfull.theta = phase;
	glPushMatrix(); // just to be sure...

					//transformToScopeCoords();

					// now (then), put it into the right pitch coordinates...
					//glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

					// now (first), roll it to get right perp angle...
	float rollang = Stimulus.stimrot.mean;
	Sfull.theta = rollang;
	glRotatef(rollang, 0.0f, 1.0f, 0.0f);
	//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);


	glBegin(GL_QUADS);
	for (int iii = 0; iii < 40 / wedgeWidthDeg; iii++) // through each wedge...
	{
		heighty = iii * wedgewidth;
		for (int ii = 0; ii < 360 / wedgeWidthDeg; ii++) // through each wedge...
		{
			angle = ii * wedgewidth + added;

			glTexCoord1f(ii*wedgeWidthDeg / per);
			glVertex3f(radius*sin(heighty)*cos(angle), radius*sin(heighty)*sin(angle), radius*cos(heighty));

			glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
			glVertex3f(radius*sin(heighty + wedgewidth)*cos(angle + wedgewidth), radius*sin(heighty + wedgewidth)*sin(angle + wedgewidth), radius*cos(heighty));

			glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
			glVertex3f(radius*sin(heighty + wedgewidth)*cos(angle + wedgewidth), radius*sin(heighty + wedgewidth)*sin(angle + wedgewidth), radius*cos(heighty + heighty));

			glTexCoord1f(ii*wedgeWidthDeg / per);
			glVertex3f(radius*sin(heighty)*cos(angle), radius*sin(heighty)*sin(angle), radius*cos(heighty + heighty));
		}
	}
	glEnd();

	//Distort(0, 0, 1, 1, DISTORT_ASIN_R, 10);

	glDisable(GL_TEXTURE_1D);
	glPopMatrix();

	//const GLubyte* sssExtensions = glGetString(GL_EXTENSIONS);

};

void TrackStim::drawSineCylinder(float tc)
{
	const int domainy = 1024;
	static float imagedata[domainy * 3];

	float added = PI / 2;	// add the stimulus OL
						// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
						// we want to choose this value to be the temporal frequency here.
						//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																		// do tc-epochtimezero so that the start is at the same position every time
	added = fmod(double(added - PI / 2), double(Stimulus.spacing * PI / 180)); // keep added to within a single period. this should prevent the seam from appearing. 130628 leongjcs
	Sfull.x = added - PI / 2; // save the phase

	for (int ii = 0; ii < domainy; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(2 * PI * ii / domainy)), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(2 * PI * ii / domainy)), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(2 * PI * ii / domainy)), 'B');
	}

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, domainy, 0, GL_RGB, GL_FLOAT, imagedata);

	float wedgeWidthDeg = 1;
	float wedgewidth = wedgeWidthDeg * PI / 180; // width of each wedge
	float phase = Stimulus.stimrot.phase * PI / 180;
	float per = Stimulus.spacing;
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	float angle;
	//Sfull.theta = phase;

	glPushMatrix(); // just to be sure...

	//transformToScopeCoords();

	//now (then), put it into the right pitch coordinates...
	//glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

	// now (first), roll it to get right perp angle...
	float rollang = Stimulus.stimrot.mean;
	Sfull.theta = rollang;
	glRotatef(rollang, 0.0f, 1.0f, 0.0f);
	//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);

	glBegin(GL_QUADS);

	for (int ii = 0; ii < 360 / wedgeWidthDeg; ii++) // through each wedge...
	{
		angle = ii * wedgewidth + added;

		//glColor3f(0.0f, 0.2f, 0.0f);
		glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), -height / 2);

		//glColor3f(0.0f, 0.2f, 0.0f);
		glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle + wedgewidth), radius*sin(angle + wedgewidth), -height / 2);

		//glColor3f(0.0f, 0.2f, 0.0f);
		glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle + wedgewidth), radius*sin(angle + wedgewidth), height / 2);

		//glColor3f(0.0f, 0.2f, 0.0f);
		glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), height / 2);
	}

	glEnd();

	//Distort(0, 0, 1, 1, DISTORT_ASIN_R, 10);

	glDisable(GL_TEXTURE_1D);
	glPopMatrix();

	//const GLubyte* sssExtensions = glGetString(GL_EXTENSIONS);

};

void TrackStim::drawSineCylinder_Standing()
{
	static bool FIRST = true;
	static float imagedata[128 * 3];

	for (int ii = 0; ii < 128; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*sin(ii * 2 * PI / 128.0f)), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*sin(ii * 2 * PI / 128.0f)), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*sin(ii * 2 * PI / 128.0f)), 'B');
	}
	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);

	float wedgeWidthDeg = 5;
	float wedgewidth = wedgeWidthDeg * PI / 180; // width of each wedge
	float phase = Stimulus.stimrot.phase * PI / 180 + PI / 2; // extra PI / 2 offset to make sure the seam does not show
	float per = Stimulus.spacing;
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	float angle;
	Sfull.theta = phase;

	float fulltime = Sfull.tcurr - epochtimezero;

	if (fulltime < Stimulus.tau)
	{
		Sfull.x = Stimulus.lum * (1 + Stimulus.contrast);

		glPushMatrix(); // just to be sure...

		// now (then), put it into the right pitch coordinates...
		glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

		// now (first), roll it to get right perp angle...
		glRotatef(Stimulus.stimrot.mean, 1.0f, 0.0f, 0.0f);
		//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, 1);
		glBegin(GL_QUADS);
		for (int ii = 0; ii < 360 / wedgeWidthDeg; ii++) // through each wedge...
		{
			angle = ii * wedgewidth + phase;

			glTexCoord1f(ii*wedgeWidthDeg / per);
			glVertex3f(radius*cos(angle), radius*sin(angle), -height / 2);

			glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
			glVertex3f(radius*cos(angle + wedgewidth), radius*sin(angle + wedgewidth), -height / 2);

			glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
			glVertex3f(radius*cos(angle + wedgewidth), radius*sin(angle + wedgewidth), height / 2);

			glTexCoord1f(ii*wedgeWidthDeg / per);
			glVertex3f(radius*cos(angle), radius*sin(angle), height / 2);
		}
		glEnd();
		glDisable(GL_TEXTURE_1D);

		glPopMatrix();
	}
	else
	{
		Sfull.x = Stimulus.lum * (1 - Stimulus.contrast);
	}
};

void TrackStim::drawSineCylinder_Standing_ContrastModulated(float tc)
{
	static float imagedata[128 * 3];

	float added = 0;	// add the stimulus OL
						// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
						// we want to choose this value to be the temporal frequency here.
						//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero) * (Stimulus.stimtrans.mean / Stimulus.spacing);	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																						// do tc-epochtimezero so that the start is at the same position every time
	Sfull.x = added;

	//recalculate the texture with the appropriate contrast modulation on every call
	for (int ii = 0; ii < 128; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(ii * 2 * PI / 128.0f) * sin(added * 2 * PI)), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(ii * 2 * PI / 128.0f) * sin(added * 2 * PI)), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * sin(ii * 2 * PI / 128.0f)  * sin(added * 2 * PI)), 'B');
	}
	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);

	float wedgeWidthDeg = 5;
	float wedgewidth = wedgeWidthDeg * PI / 180; // width of each wedge
	float phase = (Stimulus.stimrot.phase / 360) * Stimulus.spacing * PI / 180 + PI / 2;
	float per = Stimulus.spacing;
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	float angle;
	Sfull.theta = phase;

	glPushMatrix(); // just to be sure...

	// now (then), put it into the right pitch coordinates...
	glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

	// now (first), roll it to get right perp angle...
	//Sfull.theta = Stimulus.stimrot.mean;
	glRotatef(Stimulus.stimrot.mean, 1.0f, 0.0f, 0.0f);
	//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);

	glBegin(GL_QUADS);
	for (int ii = 0; ii < 360 / wedgeWidthDeg; ii++) // through each wedge...
	{
		angle = ii * wedgewidth + phase;

		glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), -height / 2);

		glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle + wedgewidth), radius*sin(angle + wedgewidth), -height / 2);

		glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle + wedgewidth), radius*sin(angle + wedgewidth), height / 2);

		glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), height / 2);
	}
	glEnd();
	glDisable(GL_TEXTURE_1D);

	glPopMatrix();

};

void TrackStim::drawCylinderStripes_Drifting(float tc)
{
	static float imagedata[128 * 3];

	float added = PI / 2;	// add the stimulus OL
							// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
							// we want to choose this value to be the temporal frequency here.
							//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																		// do tc-epochtimezero so that the start is at the same position every time
	added = fmod(double(added - PI / 2), double(Stimulus.spacing * PI / 180)); // keep added to within a single period. this should prevent the seam from appearing. 130628 leongjcs
	Sfull.x = added - PI / 2; // save the phase

	//float added = PI / 4;	// add the stimulus OL
						// offset by PI/2 rotation to get away from the seam. at 0 offset it appears on the edge of the screen. 130628 leongjcs
						// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
						// we want to choose this value to be the temporal frequency here.
						//i believe the units of tc to be seconds, leongjcs 130611
	//added += (tc - epochtimezero)*Stimulus.stimtrans.mean * PI / 180;	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																	   // do tc-epochtimezero so that the start is at the same position every time
	//added = (double(added - PI / 2), double(Stimulus.spacing * PI / 180)) + PI / 2;	// keep added to within a single period. this should prevent the seam from appearing. 130628 leongjcs
	//Sfull.x = added - PI / 2; // save the phase

	for (int ii = 0; ii < 128; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast* -1 * sign(sin(ii * 2 * PI / 128.0f))), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast* -1 * sign(sin(ii * 2 * PI / 128.0f))), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast* -1 * sign(sin(ii * 2 * PI / 128.0f))), 'B');
	}
	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);

	float wedgeWidthDeg = 5;
	float wedgewidth = wedgeWidthDeg * PI / 180; // width of each wedge
	float phase = Stimulus.stimrot.phase * PI / 180;
	float per = Stimulus.spacing;
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	float angle;
	//Sfull.theta = phase;

	glPushMatrix(); // just to be sure...

	// now (then), put it into the right pitch coordinates...
	glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

	// now (first), roll it to get right perp angle...
	Sfull.theta = Stimulus.stimrot.mean;
	glRotatef(Stimulus.stimrot.mean, 1.0f, 0.0f, 0.0f);
	//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);

	glBegin(GL_QUADS);
	for (int ii = 0; ii < 360 / wedgeWidthDeg; ii++) // through each wedge...
	{
		angle = ii * wedgewidth + added;

		glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle - wedgewidth), radius*sin(angle - wedgewidth), -height / 2);

		glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), -height / 2);

		glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), height / 2);

		glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle - wedgewidth), radius*sin(angle - wedgewidth), height / 2);
	}
	glEnd();
	glDisable(GL_TEXTURE_1D);

	glPopMatrix();

};

void TrackStim::drawCylinderStripe_Standing(float tc)
{
	static float imagedata[128 * 3];

	float added = 0;	// add the stimulus OL
						// Stimulus.stimtrans.mean/Stimulus.spacing is the contrast frequency of drifting gratings, as mean is the temporal freq and spacing is the spat freq
						// we want to choose this value to be the temporal frequency here.
						//i believe the units of tc to be seconds, leongjcs 130611
	added += (tc - epochtimezero) * (Stimulus.stimtrans.mean / (Stimulus.spacing * 2));	// put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
																						// do tc-epochtimezero so that the start is at the same position every time
	Sfull.x = added;

	//recalculate the texture with the appropriate contrast modulation on every call
	for (int ii = 0; ii < 128; ii++)
	{
		imagedata[ii * 3] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * -1 * sign(sin(ii * 2 * PI / 128.0f) * sin(added * 2 * PI))), 'R');
		imagedata[ii * 3 + 1] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * -1 * sign(sin(ii * 2 * PI / 128.0f) * sin(added * 2 * PI))), 'G');
		imagedata[ii * 3 + 2] = getDLPColor(Stimulus.lum * (1 + Stimulus.contrast * -1 * sign(sin(ii * 2 * PI / 128.0f) * sin(added * 2 * PI))), 'B');
	}
	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);

	float wedgeWidthDeg = Stimulus.spacing;
	float wedgewidth = wedgeWidthDeg * PI / 180; // width of each wedge
	float phase = Stimulus.stimrot.phase * PI / 180;
	float per = Stimulus.spacing * 2;
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	float angle;
	Sfull.theta = phase;

	glPushMatrix(); // just to be sure...

	// now (then), put it into the right pitch coordinates...
	glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

	// now (first), roll it to get right perp angle...
	//Sfull.theta = Stimulus.stimrot.mean;
	glRotatef(Stimulus.stimrot.mean, 1.0f, 0.0f, 0.0f);
	//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);

	glBegin(GL_QUADS);
	for (int ii = 0; ii < 1; ii++) // through each wedge...
	{
		angle = ii * wedgewidth + phase;

		glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle - wedgewidth), radius*sin(angle - wedgewidth), -height / 2);

		glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), -height / 2);

		glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), height / 2);

		glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle - wedgewidth), radius*sin(angle - wedgewidth), height / 2);
	}
	glEnd();
	glDisable(GL_TEXTURE_1D);

	glBegin(GL_QUADS);
	for (int ii = 1; ii < 360 / wedgeWidthDeg; ii++) // through each wedge...
	{
		angle = ii * wedgewidth + phase;

		float currlum = 0.5;
		glColor3f(getDLPColor(currlum, 'R'), getDLPColor(currlum, 'G'), getDLPColor(currlum, 'B'));

		//glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle - wedgewidth), radius*sin(angle - wedgewidth), -height / 2);

		//glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), -height / 2);

		//glTexCoord1f((ii + 1)*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle), radius*sin(angle), height / 2);

		//glTexCoord1f(ii*wedgeWidthDeg / per);
		glVertex3f(radius*cos(angle - wedgewidth), radius*sin(angle - wedgewidth), height / 2);
	}
	glEnd();

	glPopMatrix();

};

void TrackStim::drawSquareWaveCylinder(float tc)
{

	float wedgewidth = Stimulus.spacing / 2; // width of each wedge
	float per = Stimulus.spacing;
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	static float contrast = 0;
	float angle;

	float added = 0; // add the stimulus OL
	added += tc * Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!

	static bool FIRST = TRUE;
	if (FIRST)
		contrast = Stimulus.contrast;
	FIRST = FALSE;

	glPushMatrix(); // just to be sure...

	// now (then), put it into the right pitch coordinates...
	glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

	// now (first), roll it to get right perp angle...
	float rollang = Stimulus.stimrot.mean;
	glRotatef(rollang, 1.0f, 0.0f, 0.0f);
	//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

	float time = Sfull.tcurr - epochtimezero;
	//if (time < float(epochchoose+1)/10)
	if (time > 2)
	{
		Stimulus.contrast = contrast;
		setForeground();
		for (int ii = 0; ii<int(360.001 / per); ii++)
		{
			drawCylinderWedge(added + ii * per, wedgewidth);
		};
	}
	else
		Stimulus.contrast = 0;

	glPopMatrix();

};

void TrackStim::drawCylinderStripe(float tc)
{
	float wedgewidth = Stimulus.spacing; // width of each wedge

	float added = Stimulus.stimtrans.amp; // add the stimulus OL

	if (tc <= Stimulus.tau)
	{
		added += tc * Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!

	}
	else
	{
		added += Stimulus.tau*Stimulus.stimtrans.mean;
	}

	Sfull.x = added;

	glPushMatrix(); // just to be sure...


	// now (then), put it into the right pitch coordinates...
	glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

	// now (first), roll it to get right perp angle...
	float rollang = Stimulus.stimrot.mean;
	Sfull.theta = rollang;
	glRotatef(rollang, 1.0f, 0.0f, 0.0f);
	//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

	setForeground();

	if (Stimulus.density == 0) // set Stimulus.density to 0 to get a stripe of fixed width that drifts
	{
		// drawCylinderWedge( added, wedgewidth ); // center of the stripe at the same position as the edge
		drawCylinderWedge(added - (wedgewidth / 2), wedgewidth); // leading edge of the stripe at the same position as the (leading edge of the) edge
	}
	else if (Stimulus.density == 1) // set Stimulus.density to 1 to get an edge that drifts (expands)
	{
		drawCylinderWedge((added - Stimulus.stimtrans.amp) / 2 + Stimulus.stimtrans.amp, (added - Stimulus.stimtrans.amp));
	}

	glPopMatrix();
};

void TrackStim::drawCylinderWedge(float angle, float w)
{
	const float wedgeWidthDeg = 5;
	float w_t = wedgeWidthDeg * PI / 180;
	float w_s;

	// takes angles in degrees!
	float radius = Stimulus.base_radius;
	float height = Stimulus.base_extent;
	angle *= PI / 180;
	w *= PI / 180;

	for (w_s = angle - w / 2; w_s <= angle + w / 2 - w_t; w_s = w_s + w_t)
	{
		glBegin(GL_QUADS);
		glVertex3f(radius*cos(w_s), radius*sin(w_s), -height / 2);
		glVertex3f(radius*cos(w_s + w_t), radius*sin(w_s + w_t), -height / 2);
		glVertex3f(radius*cos(w_s + w_t), radius*sin(w_s + w_t), height / 2);
		glVertex3f(radius*cos(w_s), radius*sin(w_s), height / 2);
		glEnd();
	};

	w_t = angle + w / 2 - w_s;
	glBegin(GL_QUADS);
	glVertex3f(radius*cos(w_s), radius*sin(w_s), -height / 2);
	glVertex3f(radius*cos(w_s + w_t), radius*sin(w_s + w_t), -height / 2);
	glVertex3f(radius*cos(w_s + w_t), radius*sin(w_s + w_t), height / 2);
	glVertex3f(radius*cos(w_s), radius*sin(w_s), height / 2);
	glEnd();


};

void TrackStim::drawLocalCircle()
{
	float viewDistance = 0.05;//Stimulus.base_radius; // how far away the circle is in the openGL world
	float viewAngle = Stimulus.spacing * PI / 180; // the arc angle covered by the circle
	float radius = viewDistance * tan(viewAngle / 2); // the radius of the circle in the openGL world

	float offsetAzimuth = Stimulus.x_center; //specified in degrees
	float offsetElevation = Stimulus.y_center; // specified in degrees

	Sfull.y = 0;
	Sfull.theta = 0;

	float fulltime = Sfull.tcurr - epochtimezero;

	if (fulltime < Stimulus.tau)
	{
		Sfull.x = Stimulus.lum * (1 + Stimulus.contrast);

		glPushMatrix(); // just to be sure...

		//transformToScopeCoords();

		// now (then), put it into the right pitch coordinates...
		//glRotatef(-perspHEAD, 0.0f, 1.0f, 0.0f); // rotates for head now draw however...

		// rotate the circle about the cardinal axes so that it goes to different positions on the screen
		//glRotatef(offsetAzimuth, 0.0f, 0.0f, 1.0f); // azimuth. positive values of offsetAzimuth result in a CCW rotation viewing down the Z axis of the fly, i.e. to the left, from the fly's perspective
		//glRotatef(-offsetElevation, 0.0f, 1.0f, 0.0f); // elevation. positive values of offsetElevation result in a CW rotation viewing down the Y axis of the fly, i.e. up, from the fly's perspective

		// now (first), roll it to get right perp angle...
		float rollang = Stimulus.stimrot.mean;
		//glRotatef(rollang, 1.0f, 0.0f, 0.0f);
		//glRotatef(0.0f,-cos(perspHEAD*180.0f/PI),0.0f,sin(perspHEAD*180.0f/PI));

		/*
		for (int ii = 0; ii < 20; ii++)
		{
			float da = PI / 10; // wedge size for drawing the circle

			glBegin(GL_TRIANGLES);
			glVertex3f(viewDistance, 0, 0);
			glVertex3f(viewDistance, radius * cos((ii + 1)*da), radius * sin((ii + 1)*da));
			glVertex3f(viewDistance, radius * cos(da*ii), radius * sin(da*ii));
			glEnd();
		}

		glPopMatrix();
	}
	else
	{
		Sfull.x = Stimulus.lum * (1 - Stimulus.contrast);
	}
	*/

		for (int ii = 0; ii < 20; ii++)
		{
			float da = PI / 10; // wedge size for drawing the circle

			glBegin(GL_TRIANGLES);
			//glVertex3f(viewDistance, 0, 0);
			//glVertex3f(viewDistance, radius * cos((ii + 1)*da), radius * sin((ii + 1)*da));
			//glVertex3f(viewDistance, radius * cos(da*ii), radius * sin(da*ii));
			glVertex3f(0, viewDistance, 0);
			glVertex3f(radius * cos((ii + 1)*da), viewDistance, radius * sin((ii + 1)*da));
			glVertex3f(radius * cos(da*ii), viewDistance, radius * sin(da*ii));
			glEnd();
		}

		glPopMatrix();
	}
	else
	{
		Sfull.x = Stimulus.lum * (1 - Stimulus.contrast);
	}
};

TrackStim::STIMPARAMS * TrackStim::readParams(char * szFile) // read in stimulus file, output from excel file

{
	FILE *filein;
	int numlines;
	char linename[255];
	char oneline[255];
	float pcurr[100];
	STIMPARAMS Stims[100];
	for (int jj = 0; jj < 100; jj++)
	{
		Stims[jj] = TrackStim::initializeStimulus();
	}

	string currlinename = " "; // initialize string objects once
	string complinename = " ";
	string currline = " ";

	//filein = fopen("data/feedbackparams.txt", "rt");				// File To Load Data From
	filein = fopen(szFile, "rt");				// File To Load Data From
	for (int ii = 0; ii < 260; ii++) TrackStim::paramfilename[ii] = szFile[ii];
	char *fileNameTemp = szFile;
	char *aux = fileNameTemp;
	while (*fileNameTemp++);
	while (*fileNameTemp-- != '\\' && fileNameTemp != aux);
	char *fileNameTempTemp = (aux == fileNameTemp) ? fileNameTemp : fileNameTemp + 6;

	//determine if we are using old style param file or new instruction file
	char stringyOneOne[260];
	string stringyOne = " ";
	TrackStim::readstring(filein, oneline);
	sscanf(oneline, "%s", &stringyOneOne);
	stringyOne = stringyOneOne;
	string stringyTwo = "PARAMS";
	
	if (TrackStim::compareStrings(stringyOne, stringyTwo)) {
		// if old style, read in initial values about file
		sscanf(oneline, "PARAMS %d", &numlines);
		TrackStim::readstring(filein, oneline);
		sscanf(oneline, "EPOCHS %d", &numepochs);
		TrackStim::readstring(filein, oneline);
		sscanf(oneline, "MAXRUNTIME %f", &maxRunTime);
		TrackStim::readstring(filein, oneline);
		sscanf(oneline, "STIMULUSDATAPATH %s", &stimulusDataPath);

		//if (stimulusDataPath != (string) "NULL")
		//{
		//TrackStim::readStimulusData();
		//}

		for (int loop = 0; loop < numlines; loop++)   // loop on number of parameters
		{
			// go through remaining lines...
			fscanf(filein, "%s", &linename);
			currlinename = linename;

			for (int jj = 0; jj < numepochs; jj++)    // loop on number of parameter sets
			{
				fscanf(filein, "%f", &pcurr[jj]);   // this appears to work
			};

			for (int jj = 0; jj < numepochs; jj++)
			{
				complinename = "Stimulus.stimtype";
				if ((currlinename == complinename))  // does this work, or do i have to declare a new one here?
					Stims[jj].stimtype = pcurr[jj];
				complinename = "Stimulus.lum";
				if ((currlinename == complinename))
					Stims[jj].lum = pcurr[jj];
				complinename = "Stimulus.contrast";
				if ((currlinename == complinename))
					Stims[jj].contrast = pcurr[jj];
				complinename = "Stimulus.duration";
				if ((currlinename == complinename))
					Stims[jj].duration = pcurr[jj];
				complinename = "Stimulus.contrastVal1";
				if ((currlinename == complinename))
					Stims[jj].contrastVal1 = pcurr[jj];
				complinename = "Stimulus.contrastVal2";
				if ((currlinename == complinename))
					Stims[jj].contrastVal2 = pcurr[jj];
				complinename = "Stimulus.transgain";
				if ((currlinename == complinename))
					Stims[jj].transgain = pcurr[jj];
				complinename = "Stimulus.rotgain";
				if ((currlinename == complinename))
					Stims[jj].rotgain = pcurr[jj];
				complinename = "Stimulus.trans.mean";
				if ((currlinename == complinename))
					Stims[jj].trans.mean = pcurr[jj];
				complinename = "Stimulus.trans.amp";
				if ((currlinename == complinename))
					Stims[jj].trans.amp = pcurr[jj];
				complinename = "Stimulus.trans.per";
				if ((currlinename == complinename))
					Stims[jj].trans.per = pcurr[jj];
				complinename = "Stimulus.trans.phase";
				if ((currlinename == complinename))
					Stims[jj].trans.phase = pcurr[jj];
				complinename = "Stimulus.rot.mean";
				if ((currlinename == complinename))
					Stims[jj].rot.mean = pcurr[jj];
				complinename = "Stimulus.rot.amp";
				if ((currlinename == complinename))
					Stims[jj].rot.amp = pcurr[jj];
				complinename = "Stimulus.rot.per";
				if ((currlinename == complinename))
					Stims[jj].rot.per = pcurr[jj];
				complinename = "Stimulus.rot.phase";
				if ((currlinename == complinename))
					Stims[jj].rot.phase = pcurr[jj];
				complinename = "Stimulus.stimrot.mean";
				if ((currlinename == complinename))
					Stims[jj].stimrot.mean = pcurr[jj];
				complinename = "Stimulus.stimrot.amp";
				if ((currlinename == complinename))
					Stims[jj].stimrot.amp = pcurr[jj];
				complinename = "Stimulus.stimrot.per";
				if ((currlinename == complinename))
					Stims[jj].stimrot.per = pcurr[jj];
				complinename = "Stimulus.stimrot.phase";
				if ((currlinename == complinename))
					Stims[jj].stimrot.phase = pcurr[jj];
				complinename = "Stimulus.stimtrans.mean";
				if ((currlinename == complinename))
					Stims[jj].stimtrans.mean = pcurr[jj];
				complinename = "Stimulus.stimtrans.amp";
				if ((currlinename == complinename))
					Stims[jj].stimtrans.amp = pcurr[jj];
				complinename = "Stimulus.stimtrans.per";
				if ((currlinename == complinename))
					Stims[jj].stimtrans.per = pcurr[jj];
				complinename = "Stimulus.stimtrans.phase";
				if ((currlinename == complinename))
					Stims[jj].stimtrans.phase = pcurr[jj];
				complinename = "Stimulus.stimtrans2.mean";
				if ((currlinename == complinename))
					Stims[jj].stimtrans2.mean = pcurr[jj];
				complinename = "Stimulus.stimtrans2.amp";
				if ((currlinename == complinename))
					Stims[jj].stimtrans2.amp = pcurr[jj];
				complinename = "Stimulus.stimtrans2.per";
				if ((currlinename == complinename))
					Stims[jj].stimtrans2.per = pcurr[jj];
				complinename = "Stimulus.stimtrans2.phase";
				if ((currlinename == complinename))
					Stims[jj].stimtrans2.phase = pcurr[jj];
				complinename = "Stimulus.spacing";
				if ((currlinename == complinename))
					Stims[jj].spacing = pcurr[jj];
				complinename = "Stimulus.spacing2";
				if ((currlinename == complinename))
					Stims[jj].spacing2 = pcurr[jj];
				complinename = "Stimulus.density";
				if ((currlinename == complinename))
					Stims[jj].density = pcurr[jj];
				complinename = "Stimulus.tau";
				if ((currlinename == complinename))
					Stims[jj].tau = pcurr[jj];
				complinename = "Stimulus.tau2";
				if ((currlinename == complinename))
					Stims[jj].tau2 = pcurr[jj];
				complinename = "Stimulus.arenasize";
				if ((currlinename == complinename))
					Stims[jj].arenasize = pcurr[jj];
				complinename = "Stimulus.arenaheight";
				if ((currlinename == complinename))
					Stims[jj].arenaheight = pcurr[jj];
				complinename = "Stimulus.randomize";
				if ((currlinename == complinename))
					Stims[jj].randomize = pcurr[jj];
				complinename = "Stimulus.x_center";
				if ((currlinename == complinename))
					Stims[jj].x_center = pcurr[jj];
				complinename = "Stimulus.y_center";
				if ((currlinename == complinename))
					Stims[jj].y_center = pcurr[jj];
				complinename = "Stimulus.x_center2";
				if ((currlinename == complinename))
					Stims[jj].x_center2 = pcurr[jj];
				complinename = "Stimulus.y_center2";
				if ((currlinename == complinename))
					Stims[jj].y_center2 = pcurr[jj];
				complinename = "Stimulus.base_radius";
				if ((currlinename == complinename))
					Stims[jj].base_radius = pcurr[jj];
				complinename = "Stimulus.base_extent";
				if ((currlinename == complinename))
					Stims[jj].base_extent = pcurr[jj];
				complinename = "Stimulus.USEFRUSTUM";
				if ((currlinename == complinename))
					Stims[jj].USEFRUSTUM = pcurr[jj];
				complinename = "Stimulus.optogenetics";
				if ((currlinename == complinename))
					Stims[jj].optogenetics = pcurr[jj];
			};

		}

	}
	else
	{
		//determine number of lines (number of files plus one) in instruction file 
		int countie = 0;
		string liney;
		ifstream file(paramfilename);
			while (getline(file, liney))
				countie++;

		sscanf(oneline, "MAXRUNTIME %f", &maxRunTime);
		
		int randomizeInt;
		TrackStim::readstring(filein, oneline);
		sscanf(oneline, "Stimulus.randomize %d", &randomizeInt);
		
		for (int jj = 0; jj < 100; jj++)
		{
			Stims[jj].randomize = randomizeInt;
		}
		
		float durationInt[100] = { NULL };
		int jj = 0;
		fscanf(filein, "Stimulus.duration %f", &durationInt[jj]);
		jj = 1;
		while (fscanf(filein, "%f", &durationInt[jj])) {
			jj++;
		}
		
		int durationEntries = jj-1;
		
		numepochs = countie-3;

		//read filesnames and load those files 
		for (int loop = 0; loop < numepochs; loop++)
		{
			TrackStim::readstring(filein, oneline);
			sscanf(oneline, "%s", &filenames[loop]);
			TrackStim::readSavedStim();
			fileIndex++;
		}
		
		int kk = 0;
		for (int jj = 0; jj < 100; jj++) {
			if (jj < durationEntries) { kk = jj; }
			else { kk = durationEntries; }
			Stims[jj].stimtype = texType[jj];
			Stims[jj].duration = durationInt[kk];
		}

		//copy old style param file, with default values, into new param file for writing 
		//std::ifstream inFile("..\\prototypeParamFile.txt");
		//std::ofstream outFile("..\\Oute.txt");
		//outFile << inFile.rdbuf();
		
		const char* prefixy = "C:\\JL_stimuli\\";
		const char* suffixy = fileNameTempTemp;
		char *filStrWrite = NULL;
		filStrWrite = new char[260];
		strcpy(filStrWrite, prefixy);
		strcat(filStrWrite, suffixy);

		FILE *fileout;
		fileout = fopen(filStrWrite, "w");
		delete[] filStrWrite;
		filStrWrite = nullptr;

		fprintf(fileout, "PARAMS\t%d\t\t\t\n", 31);
		fprintf(fileout, "EPOCHS\t%d\t\t\t\n", numepochs);
		fprintf(fileout, "MAXRUNTIME\t%i\t\t\t\n", maxRunTime);
		fprintf(fileout, "STIMULUSDATAPATH\t%s\t\t\t\n",  "NULL");
		////////////////////////
		fprintf(fileout, "Stimulus.stimtype\t");
		for (int jj = 0; jj < numepochs-1; jj++) { fprintf(fileout, "%d\t", texType[jj]); }
		fprintf(fileout, "%d\n", texType[numepochs-1]);
		///////////////////////
		fprintf(fileout, "Stimulus.lum\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.lum); }
		fprintf(fileout, "%.2f\n", Stimulus.lum);
		////////////////////////
		fprintf(fileout, "Stimulus.contrast\t"); 
		for (int jj = 0; jj < numepochs - 1; jj++){ fprintf(fileout, "%.2f\t", Stimulus.contrast); }
		fprintf(fileout, "%.2f\n", Stimulus.contrast);
		////////////////////////
		fprintf(fileout, "Stimulus.duration\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stims[jj].duration); }
		fprintf(fileout, "%.2f\n", Stims[numepochs - 1].duration);
		////////////////////////
		fprintf(fileout, "Stimulus.transgain\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.transgain); }
		fprintf(fileout, "%.2f\n", Stimulus.transgain);
		////////////////////////
		fprintf(fileout, "Stimulus.rotgain\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.rotgain); }
		fprintf(fileout, "%.2f\n", Stimulus.rotgain);
		////////////////////////
		fprintf(fileout, "Stimulus.trans.mean\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.trans.mean); }
		fprintf(fileout, "%.2f\n", Stimulus.trans.mean);
		////////////////////////
		fprintf(fileout, "Stimulus.trans.amp\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.trans.amp); }
		fprintf(fileout, "%.2f\n", Stimulus.trans.amp);
		////////////////////////
		fprintf(fileout, "Stimulus.trans.per\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.trans.per); }
		fprintf(fileout, "%.2f\n", Stimulus.trans.per);
		////////////////////////
		fprintf(fileout, "Stimulus.trans.phase\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.trans.phase); }
		fprintf(fileout, "%.2f\n", Stimulus.trans.phase);
		////////////////////////
		fprintf(fileout, "Stimulus.rot.mean\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.rot.mean); }
		fprintf(fileout, "%.2f\n", Stimulus.rot.mean);
		////////////////////////
		fprintf(fileout, "Stimulus.rot.amp\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.rot.amp); }
		fprintf(fileout, "%.2f\n", Stimulus.rot.amp);
		////////////////////////
		fprintf(fileout, "Stimulus.rot.per\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.rot.per); }
		fprintf(fileout, "%.2f\n", Stimulus.rot.per);
		////////////////////////
		fprintf(fileout, "Stimulus.rot.phase\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.rot.phase); }
		fprintf(fileout, "%.2f\n", Stimulus.rot.phase);
		////////////////////////
		fprintf(fileout, "Stimulus.stimtrans.mean\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.stimtrans.mean); }
		fprintf(fileout, "%.2f\n", Stimulus.stimtrans.mean);
		////////////////////////
		fprintf(fileout, "Stimulus.stimtrans.amp\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.stimtrans.amp); }
		fprintf(fileout, "%.2f\n", Stimulus.stimtrans.amp);
		////////////////////////
		fprintf(fileout, "Stimulus.stimtrans.per\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.stimtrans.per); }
		fprintf(fileout, "%.2f\n", Stimulus.stimtrans.per);
		////////////////////////
		fprintf(fileout, "Stimulus.stimtrans.phase\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.stimtrans.phase); }
		fprintf(fileout, "%.2f\n", Stimulus.stimtrans.phase);
		////////////////////////
		fprintf(fileout, "Stimulus.stimrot.mean\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%d\t", texDir[jj]); }
		fprintf(fileout, "%d\n", texDir[numepochs-1]);
		////////////////////////
		fprintf(fileout, "Stimulus.stimrot.amp\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.stimrot.amp); }
		fprintf(fileout, "%.2f\n", Stimulus.stimrot.amp);
		////////////////////////
		fprintf(fileout, "Stimulus.stimrot.per\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.stimrot.per); }
		fprintf(fileout, "%.2f\n", Stimulus.stimrot.per);
		////////////////////////
		fprintf(fileout, "Stimulus.stimrot.phase\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.stimrot.phase); }
		fprintf(fileout, "%.2f\n", Stimulus.stimrot.phase);
		////////////////////////
		fprintf(fileout, "Stimulus.spacing\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", texSpacing[jj]); }
		fprintf(fileout, "%.2f\n", texSpacing[numepochs-1]);
		////////////////////////
		fprintf(fileout, "Stimulus.density\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.density); }
		fprintf(fileout, "%.2f\n", Stimulus.density);
		////////////////////////
		fprintf(fileout, "Stimulus.tau\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.tau); }
		fprintf(fileout, "%.2f\n", Stimulus.tau);
		////////////////////////
		fprintf(fileout, "Stimulus.arenasize\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.arenasize); }
		fprintf(fileout, "%.2f\n", Stimulus.arenasize);
		////////////////////////
		fprintf(fileout, "Stimulus.arenaheight\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.arenaheight); }
		fprintf(fileout, "%.2f\n", Stimulus.arenaheight);
		////////////////////////
		fprintf(fileout, "Stimulus.tau2\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.tau2); }
		fprintf(fileout, "%.2f\n", Stimulus.tau2);
		////////////////////////
		fprintf(fileout, "Stimulus.randomize\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%d\t", randomizeInt); }
		fprintf(fileout, "%d\n", Stimulus.randomize);
		////////////////////////
		fprintf(fileout, "Stimulus.spacing2\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.spacing2); }
		fprintf(fileout, "%.2f\n", Stimulus.spacing2);
		////////////////////////
		fprintf(fileout, "Stimulus.USEFRUSTUM\t");
		for (int jj = 0; jj < numepochs - 1; jj++) { fprintf(fileout, "%.2f\t", Stimulus.USEFRUSTUM); }
		fprintf(fileout, "%.2f\n", Stimulus.USEFRUSTUM);

		fclose(fileout);
		
		TrackStim::readMesh();

	}
	
	fclose(filein);

	return Stims;

};
