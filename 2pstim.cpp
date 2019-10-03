// Trackball.cpp : Defines the entry point for the application.
//
#include "2pstim.h"

char * GetFileName()
{

	OPENFILENAME ofn;
	char szFile[255]; // Buffer for selected file name.
	szFile[0]=NULL;
	// blah, blah, more funtions...
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = hWnd;
	ofn.lpstrFile = szFile;
	ofn.nMaxFile = sizeof(szFile);
	ofn.lpstrFilter = "Text:\0*.txt\0";
	ofn.nFilterIndex = 0;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileName(&ofn))
		return ofn.lpstrFile;
	else
		return NULL;

};

void readstr(FILE *f,char *str)
{
	do
	{
		fgets(str, 255, f);
	} while ((str[0] == '/') || (str[0] == '\n'));
	return;
}

int ShowHideTaskBar ( BOOL bHide )
{
	// find taskbar window
	CWnd* pWnd = CWnd::FindWindow(_T("Shell_TrayWnd"), _T(""));
	if(!pWnd )
		return 0;
	if( bHide )
	{
		// hide taskbar
		pWnd->ShowWindow(SW_HIDE);
	}
	else
	{
		// show taskbar
		pWnd->ShowWindow(SW_SHOW);
	}
	return 1;
}

GLvoid ReSizeGLScene(GLsizei width, GLsizei height)		// Resize And Initialize The GL Window
{
	if (height==0)										// Prevent A Divide By Zero By
	{
		height=1;										// Making Height Equal One
	}

	glViewport(0,0,width,height);						// Reset The Current Viewport

	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();									// Reset The Projection Matrix

	// Calculate The Aspect Ratio Of The Window
	gluPerspective(90.0f,(GLfloat)width/(GLfloat)height,0.1f,100.0f);

	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix
}

int InitGL(GLvoid)										// All Setup For OpenGL Goes Here
{

	//glEnable(GL_TEXTURE_2D);							// Enable Texture Mapping
	//glBlendFunc(GL_SRC_ALPHA,GL_ONE);					// Set The Blending Function For Translucency
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);				// This Will Clear The Background Color To Black
	glClear(GL_COLOR_BUFFER_BIT);
	glClearDepth(1.0);									// Enables Clearing Of The Depth Buffer
	glDepthFunc(GL_LESS);								// The Type Of Depth Test To Do
	glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
	//glShadeModel(GL_SMOOTH);							// Enables Smooth Color Shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations

	// anti-aliasing stuff? DAC 080827
	glEnable (GL_BLEND); // needed!
	//glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);			// Set The Blending Function For Translucency
	//glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); // not great -- get the seams showing up
	glBlendFunc(GL_ONE,GL_ZERO);                         // aliased, but no shadows -- one solution to shadows is to make the clear color the background color again
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POINT_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	//glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);  // can be GL_NICEST, GL_FASTEST, GL_DONT_CARE
	//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	return TRUE;										// Initialization Went OK
}

GLvoid KillGLWindow(GLvoid)								// Properly Kill The Window
{
	if (fullscreen)										// Are We In Fullscreen Mode?
	{
		ChangeDisplaySettings(NULL,0);					// If So Switch Back To The Desktop
		ShowCursor(TRUE);								// Show Mouse Pointer
	}

	ShowHideTaskBar(FALSE);  // show taskbar when done...

	if (hRC)											// Do We Have A Rendering Context?
	{
		if (!wglMakeCurrent(NULL,NULL))					// Are We Able To Release The DC And RC Contexts?
		{
			MessageBox(NULL,"Release Of DC And RC Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		}

		if (!wglDeleteContext(hRC))						// Are We Able To Delete The RC?
		{
			MessageBox(NULL,"Release Rendering Context Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		}
		hRC=NULL;										// Set RC To NULL
	}

	if (hDC && !ReleaseDC(hWnd,hDC))					// Are We Able To Release The DC
	{
		MessageBox(NULL,"Release Device Context Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hDC=NULL;										// Set DC To NULL
	}

	if (hWnd && !DestroyWindow(hWnd))					// Are We Able To Destroy The Window?
	{
		MessageBox(NULL,"Could Not Release hWnd.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hWnd=NULL;										// Set hWnd To NULL
	}

	if (!UnregisterClass("OpenGL",hInstance))			// Are We Able To Unregister Class
	{
		MessageBox(NULL,"Could Not Unregister Class.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hInstance=NULL;									// Set hInstance To NULL
	}
}

/*	This Code Creates Our OpenGL Window.  Parameters Are:					*
 *	title			- Title To Appear At The Top Of The Window				*
 *	width			- Width Of The GL Window Or Fullscreen Mode				*
 *	height			- Height Of The GL Window Or Fullscreen Mode			*
 *	bits			- Number Of Bits To Use For Color (8/16/24/32)			*
 *	fullscreenflag	- Use Fullscreen Mode (TRUE) Or Windowed Mode (FALSE)	*/

BOOL CreateGLWindow(char* title, int width, int height, int bits, bool fullscreenflag)
{
	GLuint		PixelFormat;			// Holds The Results After Searching For A Match
	WNDCLASS	wc;						// Windows Class Structure
	DWORD		dwExStyle;				// Window Extended Style
	DWORD		dwStyle;				// Window Style
	RECT		WindowRect;				// Grabs Rectangle Upper Left / Lower Right Values
	WindowRect.left=(long)0;			// Set Left Value To 0
	WindowRect.right=(long)width;		// Set Right Value To Requested Width
	WindowRect.top=(long)0;				// Set Top Value To 0
	WindowRect.bottom=(long)height;		// Set Bottom Value To Requested Height

	fullscreen=fullscreenflag;			// Set The Global Fullscreen Flag

	hInstance			= GetModuleHandle(NULL);				// Grab An Instance For Our Window
	wc.style			= CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Redraw On Size, And Own DC For Window.
	wc.lpfnWndProc		= (WNDPROC) WndProc;					// WndProc Handles Messages
	wc.cbClsExtra		= 0;									// No Extra Window Data
	wc.cbWndExtra		= 0;									// No Extra Window Data
	wc.hInstance		= hInstance;							// Set The Instance
	wc.hIcon			= LoadIcon(NULL, IDI_WINLOGO);			// Load The Default Icon
	wc.hCursor			= LoadCursor(NULL, IDC_ARROW);			// Load The Arrow Pointer
	wc.hbrBackground	= NULL;									// No Background Required For GL
	wc.lpszMenuName		= NULL;									// We Don't Want A Menu
	wc.lpszClassName	= "OpenGL";								// Set The Class Name

	if (!RegisterClass(&wc))									// Attempt To Register The Window Class
	{
		MessageBox(NULL,"Failed To Register The Window Class.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;											// Return FALSE
	}

	if (fullscreen)												// Attempt Fullscreen Mode?
	{
		DEVMODE dmScreenSettings;								// Device Mode
		memset(&dmScreenSettings,0,sizeof(dmScreenSettings));	// Makes Sure Memory's Cleared
		dmScreenSettings.dmSize=sizeof(dmScreenSettings);		// Size Of The Devmode Structure
		dmScreenSettings.dmPelsWidth	= width;				// Selected Screen Width
		dmScreenSettings.dmPelsHeight	= height;				// Selected Screen Height
		dmScreenSettings.dmBitsPerPel	= bits;					// Selected Bits Per Pixel
		dmScreenSettings.dmDisplayFrequency = 100;              // DAC 080819 -- see if this works?
		dmScreenSettings.dmFields=DM_BITSPERPEL|DM_PELSWIDTH|DM_PELSHEIGHT;

		// Try To Set Selected Mode And Get Results.  NOTE: CDS_FULLSCREEN Gets Rid Of Start Bar.
		if (ChangeDisplaySettings(&dmScreenSettings,CDS_FULLSCREEN)!=DISP_CHANGE_SUCCESSFUL)
		{
			// If The Mode Fails, Offer Two Options.  Quit Or Use Windowed Mode.
			if (MessageBox(NULL,"The Requested Fullscreen Mode Is Not Supported By\nYour Video Card. Use Windowed Mode Instead?","NeHe GL",MB_YESNO|MB_ICONEXCLAMATION)==IDYES)
			{
				fullscreen=FALSE;		// Windowed Mode Selected.  Fullscreen = FALSE
			}
			else
			{
				// Pop Up A Message Box Letting User Know The Program Is Closing.
				MessageBox(NULL,"Program Will Now Close.","ERROR",MB_OK|MB_ICONSTOP);
				return FALSE;									// Return FALSE
			}
		}
	}

	if (useDLP)
	{
		if (!ShowHideTaskBar(TRUE)) // to hide the task bar!
			MessageBox(NULL,"Could not hide the taskbar","Taskbar error",MB_OK | MB_ICONINFORMATION);
	};

	if (fullscreen)												// Are We Still In Fullscreen Mode?
	{
		dwExStyle=WS_EX_APPWINDOW;								// Window Extended Style
		dwStyle=WS_POPUP;										// Windows Style
		ShowCursor(FALSE);										// Hide Mouse Pointer
	}
	else
	{
		dwExStyle=WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;			// Window Extended Style
		dwStyle=WS_OVERLAPPEDWINDOW;							// Windows Style
		ShowCursor(FALSE);                                      // DAC 080819
	}

	if (useDLP)
	{
		dwExStyle=WS_EX_APPWINDOW;								// Window Extended Style
		dwStyle=WS_POPUP;										// Windows Style
		ShowCursor(FALSE);										// Hide Mouse Pointer
	};

	AdjustWindowRectEx(&WindowRect, dwStyle, FALSE, dwExStyle);		// Adjust Window To True Requested Size

	// Create The Window
	if (!useDLP)
	{
		if (!(hWnd=CreateWindowEx(	dwExStyle,							// Extended Style For The Window
									"OpenGL",							// Class Name
									title,								// Window Title
									dwStyle |							// Defined Window Style
									WS_CLIPSIBLINGS |					// Required Window Style
									WS_CLIPCHILDREN,					// Required Window Style
									0, 0,								// Window Position -- can modify this! DAC
									WindowRect.right-WindowRect.left,	// Calculate Window Width
									WindowRect.bottom-WindowRect.top,	// Calculate Window Height
									NULL,								// No Parent Window
									NULL,								// No Menu
									hInstance,							// Instance
									NULL)))								// Dont Pass Anything To WM_CREATE
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL,"Window Creation Error.","ERROR",MB_OK|MB_ICONEXCLAMATION);
			return FALSE;								// Return FALSE
		};
	}
	else
	{
		if (!(hWnd=CreateWindowEx(	dwExStyle,							// Extended Style For The Window
									"OpenGL",							// Class Name
									title,								// Window Title
									//dwStyle |							// Defined Window Style
									//WS_CLIPSIBLINGS |					// Required Window Style
									//WS_CLIPCHILDREN |
									WS_BORDER | WS_POPUP,				// Required Window Style
									1920, 0, // 1600, 0, // 2970, 0, // 1680, 0,				// Window Position -- can modify this! DAC
									WindowRect.right-WindowRect.left,	// Calculate Window Width
									WindowRect.bottom-WindowRect.top,	// Calculate Window Height
									NULL,								// No Parent Window
									NULL,								// No Menu
									hInstance,							// Instance
									NULL)))								// Dont Pass Anything To WM_CREATE
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL,"Window Creation Error.","ERROR",MB_OK|MB_ICONEXCLAMATION);
			return FALSE;								// Return FALSE
		};
	}


	static	PIXELFORMATDESCRIPTOR pfd=				// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		bits,										// Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		24,											// 24 bit Z-Buffer (Depth Buffer)
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};

	if (!(hDC=GetDC(hWnd)))							// Did We Get A Device Context?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Create A GL Device Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(PixelFormat=ChoosePixelFormat(hDC,&pfd)))	// Did Windows Find A Matching Pixel Format?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Find A Suitable PixelFormat.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if(!SetPixelFormat(hDC,PixelFormat,&pfd))		// Are We Able To Set The Pixel Format?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Set The PixelFormat.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(hRC=wglCreateContext(hDC)))				// Are We Able To Get A Rendering Context?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Create A GL Rendering Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if(!wglMakeCurrent(hDC,hRC))					// Try To Activate The Rendering Context
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Activate The GL Rendering Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	ShowWindow(hWnd,SW_SHOW);						// Show The Window
	SetForegroundWindow(hWnd);						// Slightly Higher Priority
	SetFocus(hWnd);									// Sets Keyboard Focus To The Window
	ReSizeGLScene(width, height);					// Set Up Our Perspective GL Screen

	if (!InitGL())									// Initialize Our Newly Created GL Window
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Initialization Failed.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	return TRUE;									// Success
}

LRESULT CALLBACK WndProc(	HWND	hWnd,			// Handle For This Window
							UINT	uMsg,			// Message For This Window
							WPARAM	wParam,			// Additional Message Information
							LPARAM	lParam)			// Additional Message Information
{
	switch (uMsg)									// Check For Windows Messages
	{
		case WM_ACTIVATE:							// Watch For Window Activate Message
		{
			if (!HIWORD(wParam))					// Check Minimization State
			{
				active=TRUE;						// Program Is Active
			}
			else
			{
				active=FALSE;						// Program Is No Longer Active
			}

			return 0;								// Return To The Message Loop
		}

		case WM_SYSCOMMAND:							// Intercept System Commands
		{
			switch (wParam)							// Check System Calls
			{
				case SC_SCREENSAVE:					// Screensaver Trying To Start?
				case SC_MONITORPOWER:				// Monitor Trying To Enter Powersave?
				return 0;							// Prevent From Happening
			}
			break;									// Exit
		}

		case WM_CLOSE:								// Did We Receive A Close Message?
		{
			PostQuitMessage(0);						// Send A Quit Message
			return 0;								// Jump Back
		}

		case WM_KEYDOWN:							// Is A Key Being Held Down?
		{
			keys[wParam] = TRUE;					// If So, Mark It As TRUE
			return 0;								// Jump Back
		}

		case WM_KEYUP:								// Has A Key Been Released?
		{
			keys[wParam] = FALSE;					// If So, Mark It As FALSE
			return 0;								// Jump Back
		}

		case WM_SIZE:								// Resize The OpenGL Window
		{
			ReSizeGLScene(LOWORD(lParam),HIWORD(lParam));  // LoWord=Width, HiWord=Height
			//if (!fileout == NULL)
			//	fprintf (fileout,"WM_SIZE called\n");
			return 0;								// Jump Back
		}

	}

	// Pass All Unhandled Messages To DefWindowProc
	return DefWindowProc(hWnd,uMsg,wParam,lParam);
}

int WINAPI WinMain(	HINSTANCE	hInstance,			// Instance
					HINSTANCE	hPrevInstance,		// Previous Instance
					LPSTR		lpCmdLine,			// Command Line Parameters
					int			nCmdShow)			// Window Show State
{
	MSG		msg;									// Windows Message Structure
	char *paramfilename = "" ;
	char paramfilekeep[255] = "" ;    // both will be passed to trackstim and used to write things out...

	TrackStim::STIMPARAMS AllStims[100];
	TrackStim::STIMPARAMS * AllStimsTemp;
	float currtime = 0;
	
	// set up stimulus here...
	TrackStim Onscreen;   // instantiate class

	uInt32 data; // data from NIDAQ counter
	int         error=0;
	TaskHandle  taskHandle=0;
	char        errBuff[2048]={'\0'};

	TaskHandle  taskHandle2=0;

	static bool done = false;							// Bool Variable To Exit Loop
	static bool FIRST_EPOCH = TRUE;
	static int * index;
	static int cur_index = 0;
	static int lastDataFrame = -1;
	static float lastDataFrameStartTime = 0;
	

	fullscreen=FALSE;

	if (MessageBox(NULL,"Put it onto the DLP", "Start FullScreen?",MB_YESNO|MB_ICONQUESTION)==IDNO)
	{
		useDLP = FALSE;							// Windowed Mode
	}
	// read in parameter data...
	paramfilename = GetFileName() ;

	for (int ii=0; ii<255; ii++) paramfilekeep[ii]=paramfilename[ii];
	AllStimsTemp = Onscreen.readParams(paramfilekeep);
	for (int jjj=0; jjj<Onscreen.numepochs; jjj++) AllStims[jjj]=AllStimsTemp[jjj];
	Onscreen.epochchoose = 0;
	Onscreen.Stimulus = AllStims[Onscreen.epochchoose];
	Onscreen.boutInd = 0;
	Onscreen.frameIndex = 0;

	// read in viewport data & screen positioning data
	Onscreen.readViewPositions(); // if no file exists, just initializes it...

	// Create Our OpenGL Window
    if (!CreateGLWindow("TrackBall",Onscreen.windowWidth,Onscreen.windowHeight,Onscreen.bits,fullscreen))
	{
		return 0;											// Quit If Window Was Not Created
	}

	// Onscreen.setZeroTime(); // set time to be zero
	// Onscreen.epochtimezero=Onscreen.queryCurrTime();    // zero it...

	// DAQ SETUP FOR IMAGING SYNCHRONIZATION

	/*********************************************/
	// DAQmx Configure Code
	/*********************************************/
	DAQmxErrChk (DAQmxCreateTask("2",&taskHandle));
	DAQmxErrChk (DAQmxCreateCICountEdgesChan(taskHandle,"/Dev1/ctr1","",DAQmx_Val_Rising,0,DAQmx_Val_CountUp));
	//DAQmxErrChk (DAQmxCreateCICountEdgesChan(taskHandle,"Dev1/ctr1","",DAQmx_Val_CI,0,DAQmx_Val_CountUp));
	//DAQmxErrChk (DAQmxCreateCICountEdgesChan(taskHandle,"Dev1/pfi0","",DAQmx_Val_CI,0,DAQmx_Val_Low));

	/*********************************************/
	// DAQmx Start Code
	/*********************************************/
	DAQmxErrChk (DAQmxStartTask(taskHandle));

	bool first = TRUE;

	while(!done)									// Loop That Runs While done=FALSE
	{
		if (PeekMessage(&msg,NULL,0,0,PM_REMOVE))	// Is There A Message Waiting?
		{
			if (msg.message==WM_QUIT)				// Have We Received A Quit Message?
			{
				done = true;							// If So done=TRUE
				continue;
			}
			else									// If Not, Deal With Window Messages
			{
				TranslateMessage(&msg);				// Translate The Message
				DispatchMessage(&msg);				// Dispatch The Message
			}
		}
		else										// If There Are No Messages
		{
			if(first==TRUE)
			{
				first = FALSE;

				/*********************************************/
				// DAQmx Configure Code
				/*********************************************/
				DAQmxErrChk (DAQmxCreateTask("1",&taskHandle2));
				//DAQmxErrChk (DAQmxCreateCOPulseChanTime(taskHandle2,"Dev1/ctr0","",DAQmx_Val_Seconds,DAQmx_Val_Low,1.00,0.50,1.00));
				DAQmxErrChk (DAQmxCreateCOPulseChanTime(taskHandle2,"Dev1/ctr0","",DAQmx_Val_Seconds,DAQmx_Val_Low,0,0.05,0.05));

				/*********************************************/
				// DAQmx Start Code
				/*********************************************/
				DAQmxErrChk (DAQmxStartTask(taskHandle2));

				/*********************************************/
				// DAQmx Wait Code
				/*********************************************/
				//DAQmxErrChk (DAQmxWaitUntilTaskDone(taskHandle2,5.0));

				Onscreen.setZeroTime(); // set time to be zero
				//Onscreen.epochtimezero=Onscreen.queryCurrTime();    // zero it...

				// Ask The User Which Screen Mode They Prefer
				/*if (MessageBox(NULL,"Would You Like To Run In Fullscreen Mode?", "Start FullScreen?",MB_YESNO|MB_ICONQUESTION)==IDNO)
				{
					fullscreen=FALSE;							// Windowed Mode
				};*/

				Error:
				{
					puts("");
					if( DAQmxFailed(error) )
						DAQmxGetExtendedErrorInfo(errBuff,2048);
					//if( taskHandle!=0 ) {
						/*********************************************/
						// DAQmx Stop Code
						/*********************************************/
					//	DAQmxStopTask(taskHandle);
					//	DAQmxClearTask(taskHandle);
					//}
					if( DAQmxFailed(error) )
					{
						char tmp[500];
						printf("DAQmx Error: %s\n",errBuff);
						sprintf(tmp,"DAQmx Error 2: %s\n",errBuff);
						printf("End of program, press Enter key to quit\n");
						MessageBox(NULL,tmp,"End of program",MB_YESNO|MB_ICONQUESTION)==IDNO;
						getchar();
						return 0;
					}
				}
			}

			currtime = Onscreen.queryCurrTime();

			DAQmxErrChk (DAQmxReadCounterScalarU32(taskHandle,1.0,&data,NULL));
			if (lastDataFrame != data)
			{
				lastDataFrame = data;
				lastDataFrameStartTime = currtime;
			}

			// Draw The Scene.  Watch For ESC Key And Quit Messages From DrawGLScene()
			//if ((active && !Onscreen.drawScene()) || keys[VK_ESCAPE])	// not most recent mouse commands in draw?
			if ((keys[VK_ESCAPE]) || ((Onscreen.maxRunTime != 0) && ((currtime - lastDataFrameStartTime > 1) || (currtime >= Onscreen.maxRunTime))))
			{
				done = true;							// ESC or DrawGLScene Signalled A Quit
				continue;
			}
			else									// Not Time To Quit, Update Screen
			{
				// do all the good stuff here, updating things, looking for things, etc., including mouse update here...
				if (cur_index == 0) // if it's the very first epoch presentation or we've reached the end of a run
				{
					switch (Onscreen.Stimulus.randomize) //initialize the next run
					{
					case 0:
						index = new int[Onscreen.numepochs];
						for (int ii = 0; ii < Onscreen.numepochs; ii++)
							index[ii]=ii;
						break;
					case 1:
						srand ( time(NULL) );
						index = new int[Onscreen.numepochs-1];
						for (int ii = 0; ii< Onscreen.numepochs-1; ii++)
							index[ii]=ii+1;
						for (int ii = 0; ii < Onscreen.numepochs-1; ii++)
						{
							int randNum, temporary;
							randNum = ( rand() % (Onscreen.numepochs-1) );
							temporary = index[ii];
							index[ii] = index[randNum];
							index[randNum] = temporary;
						}
						break;
					case 2:
						srand ( time(NULL) );
						index = new int[Onscreen.numepochs];
						for (int ii = 0; ii < Onscreen.numepochs; ii++)
							index[ii]=ii;
						for (int ii = 0; ii < Onscreen.numepochs; ii++)
						{
							int randNum, temporary;
							randNum = ( rand() % Onscreen.numepochs );
							temporary = index[ii];
							index[ii] = index[randNum];
							index[randNum] = temporary;
						}
						break;
					}
				}

				if ((FIRST_EPOCH) || ((Onscreen.Stimulus.duration > 0) && (currtime - Onscreen.epochtimezero >= Onscreen.Stimulus.duration))) // if it's the first epoch, or the epoch presentation has finite duration and has ended
				{
					FIRST_EPOCH = false;
					switch (Onscreen.Stimulus.randomize) //update epochchoose
					{
					case 0:
						Onscreen.epochchoose = index[cur_index];
                        cur_index = (cur_index+1) % Onscreen.numepochs;
//						Onscreen.epochchoose = (Onscreen.epochchoose + 1) % (Onscreen.numepochs);
						break;
					case 1:
						if (cur_index % 2 == 1)
							Onscreen.epochchoose = index[(cur_index-1)/2];
						else
							Onscreen.epochchoose = 0; // set back to background...
						cur_index = (cur_index+1) % (2*(Onscreen.numepochs-1));
						break;
					case 2:
                        Onscreen.epochchoose = index[cur_index];
                        cur_index = (cur_index+1) % Onscreen.numepochs;
						break;
					}
					Onscreen.Stimulus = AllStims[Onscreen.epochchoose];
					Onscreen.epochtimezero=currtime;    // zero the epoch time
					Onscreen.boutInd++; // increment the bout index. TrackStim initializes this value to 0, but in the data it starts at 1
				}

				Onscreen.drawScene();	// has own Fcurr and Sfull update lines...
				//Onscreen.forceWait();	// kluge since vsyncing isn't working... d'oh! 081107 -- 091008 commented out to get vsyncing working...
				SwapBuffers(hDC);		// Swap Buffers (Double Buffering)
				Onscreen.recordStimInfo(taskHandle); //record stimulus info into arrays, later write to file in writeStim
				Onscreen.incrementFrameNumber(); // important to keep track of 
			}
		}
	}

	switch (Onscreen.texPrecision[0]) //only noise will have precision 8 or 16
	{
	case 8:	
		for (int iiii = 0; iiii < Onscreen.numepochs; iiii++) { 
			delete Onscreen.noiseTemplate[iiii];
			delete Onscreen.movie[iiii];
		}
		break;
	case 32:
		for (int iiii = 0; iiii < Onscreen.numepochs; iiii++) { 
			delete Onscreen.movie[iiii]; 
		}
		break;
	}
	delete[] index;
	delete[] Onscreen.meshx;
	delete[] Onscreen.meshy;
	delete[] Onscreen.meshu;
	delete[] Onscreen.meshv;
	delete[] Onscreen.meshi;

	Onscreen.writeStim(taskHandle); //record stim info to file 

	char cmd1 [1024];
	sprintf(cmd1,"copy %s %s",paramfilekeep,Onscreen.outputDir);
	system(cmd1);

    // Stop and clear NIDAQ tasks
    DAQmxErrChk( DAQmxStopTask( taskHandle ) );
    DAQmxErrChk( DAQmxClearTask( taskHandle ) );
    DAQmxErrChk( DAQmxStopTask( taskHandle2 ) );
    DAQmxErrChk( DAQmxClearTask( taskHandle2 ) );

	// Shutdown
	KillGLWindow();										// Kill The Window
	return (msg.wParam);								// Exit The Program
}
