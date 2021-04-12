//	USMCDLL.h	- header file with declorations of USMCDLL.dll functions and used structures
//	-SD
//  Note::C++
//  Use functions with reference arguments as in sample
//  Note.C
//	Use functions with address of operator in invoke expression: func(var1, &struct);
//	(For Rev. 24) Will work correctly with firmware version 24.03 or better

#pragma once

#include <windows.h>

// Please use the USMCDLL.lib to easely import functions from dll with "dllimport" (In Microsoft Visual Studio)
#ifdef USMCDLL_EXPORTS
#define USMCDLL_API __declspec(dllexport)
#else
#define USMCDLL_API __declspec(dllimport)
#endif

// Structure representing connected devices
typedef struct USMC_Devices_st{
	DWORD NOD;			// Number of the devices ready to work

	char **Serial;		// Array of 16 byte ASCII strings
	char **Version;		// Array of 4 byte ASCII strings
} USMC_Devices;			// Structure representing connected devices

// Structure representing some of divice parameters
typedef struct USMC_Parameters_st{
    float AccelT;		// Acceleration time (in ms)
    float DecelT;		// Deceleration time (in ms)
    float PTimeout;		// Time (in ms) after which current will be reduced to 60% of normal
    float BTimeout1;	// Time (in ms) after which speed of step motor rotation will be equal to the one specified at
						// BTO1P field (see below). (This parameter is used when controlling step motor using buttons)
    float BTimeout2;	//
    float BTimeout3;	//
    float BTimeout4;	//
    float BTimeoutR;	// Time (in ms) after which reset command will be performed
    float BTimeoutD;	// This field is reserved for future use
    float MinP;			// Speed (steps/sec) while performing reset operation. (This parameter is used when controlling
						// step motor using buttons)
    float BTO1P;		// Speed (steps/sec) after BTIMEOUT 1 time have passed. (This parameter is used when controlling
						// step motor using buttons)
    float BTO2P;		//
    float BTO3P;		//
    float BTO4P;		//
    WORD MaxLoft;		// Value in full steps that will be used performing backlash operation
    DWORD StartPos;		// Current Position Saved to FLASH (see Test MicroSMC.cpp)
	WORD RTDelta;		// Revolution distance – number of full steps per one full revolution
    WORD RTMinError;	// Number of full steps missed to raise the error flag
	float MaxTemp;		// Maximum allowed temperature (Celsius)
	BYTE SynOUTP;		// Duration of the output synchronization pulse
	float LoftPeriod;	// Speed of the last phase of the backlash operation.
	float EncMult;		// Should be <Encoder Steps per Evolution> / <SM Steps per Evolution> and should be integer multiplied by 0.25

	BYTE Reserved[16];	// <Unused> File padding

} USMC_Parameters;

// Structure representing start function parameters
typedef struct USMC_StartParameters_st{
	BYTE SDivisor;		// Step is divided by this factor (1,2,4,8)
	BOOL DefDir;		// Direction for backlash operation (relative)
    BOOL LoftEn;		// Enable automatic backlash operation (works if slow start/stop mode is off)
	BOOL SlStart;		// If TRUE slow start/stop mode enabled.
	BOOL WSyncIN;		// If TRUE controller will wait for input synchronization signal to start
	BOOL SyncOUTR;		// If TRUE output synchronization counter will be reset
	BOOL ForceLoft;		// If TRUE and destination position is equal to the current position backlash operation will be performed.
	BYTE Reserved[4];	// <Unused> File padding
} USMC_StartParameters;

// Structure representing some of divice parameters
typedef struct USMC_Mode_st{
    BOOL PMode;			// Turn off buttons (TRUE - buttons disabled)
    BOOL PReg;			// Current reduction regime (TRUE - regime is on)
    BOOL ResetD;		// Turn power off and make a whole step (TRUE - apply)
    BOOL EMReset;		// Quick power off
    BOOL Tr1T;			// Trailer 1 TRUE state (TRUE : +3/+5Â; FALSE : 0Â)
    BOOL Tr2T;			// Trailer 2 TRUE state (TRUE : +3/+5Â; FALSE : 0Â)
    BOOL RotTrT;		// Rotary Transducer TRUE state (TRUE : +3/+5Â; FALSE : 0Â)
    BOOL TrSwap;		// If TRUE, Trailers are treated to be swapped
    BOOL Tr1En;			// If TRUE Trailer 1 Operation Enabled
    BOOL Tr2En;			// If TRUE Trailer 2 Operation Enabled
    BOOL RotTeEn;		// If TRUE Rotary Transducer Operation Enabled
    BOOL RotTrOp;		// Rotary Transducer Operation Select (stop on error for TRUE)
    BOOL Butt1T;		// Button 1 TRUE state (TRUE : +3/+5Â; FALSE : 0Â)
    BOOL Butt2T;		// Button 2 TRUE state (TRUE : +3/+5Â; FALSE : 0Â)
    BOOL ResetRT;		// Reset Rotary Transducer Check Positions (need one full revolution before it can detect error)
    BOOL SyncOUTEn;		// If TRUE output syncronization enabled
    BOOL SyncOUTR;		// If TRUE output synchronization counter will be reset
    BOOL SyncINOp;		// Synchronization input mode:
						// True - Step motor will move one time to the destination position
						// False - Step motor will move multiple times by steps equal to the value destination position
	DWORD SyncCount;	// Number of steps after which synchronization output sygnal occures
	BOOL SyncInvert;	// Set to TRUE to invert output synchronization signal

	BOOL EncoderEn;		// Enable Encoder on pins {SYNCIN,ROTTR} - disables Synchronization input and Rotary Transducer
	BOOL EncoderInv;	// Invert Encoder Counter Direction
	BOOL ResBEnc;		// Reset <Encoder Position> and <SM Position in Encoder units> to 0
	BOOL ResEnc;		// Reset <SM Position in Encoder units> to <Encoder Position>

	BYTE Reserved[8];	// <Unused> File padding

} USMC_Mode;

// Structure representing divice state
typedef struct USMC_State_st{
	int CurPos;			// Current position (in microsteps)
	float Temp;			// Current temperature of the driver
	BYTE SDivisor;		// Step is divided by this factor
	BOOL Loft;			// Indicates backlash status
	BOOL FullPower;		// If TRUE then full power.
	BOOL CW_CCW;		// Current direction. Relatively!
	BOOL Power;			// If TRUE then Step Motor is ON.
	BOOL FullSpeed;		// If TRUE then full speed. Valid in "Slow Start" mode only.
	BOOL AReset;		// TRUE After Device reset, FALSE after "Set Position".
	BOOL RUN;			// Indicates if step motor is rotating
	BOOL SyncIN;		// Logical state directly from input synchronization PIN
	BOOL SyncOUT;		// Logical state directly from output synchronization PIN
	BOOL RotTr;			// Indicates current rotary transducer press state
	BOOL RotTrErr;		// Indicates rotary transducer error flag
	BOOL EmReset;		// Indicates state of emergency disable button (local control)
	BOOL Trailer1;		// Indicates trailer 1 logical state.
	BOOL Trailer2;		// Indicates trailer 2 logical state.
	float Voltage;		// Input power source voltage (6-39V) -=24 version 0nly=-

	BYTE Reserved[8];	// <Unused> File padding
} USMC_State;

// New For Firmware Version 2.4.1.0 (0x2410)
typedef struct USMC_EncoderState_st{
	int EncoderPos;		// Current position measured by encoder
	int ECurPos;		// Current position (in Encoder Steps) - Synchronized with request call

	BYTE Reserved[8];	// <Unused> File padding
} USMC_EncoderState;

// ~New For Firmware Version 2.4.1.0 (0x2410)

typedef struct USMC_Info_st{
	char serial[17];
	DWORD dwVersion;    
    char DevName[32];
	int CurPos, DestPos;
	float Speed;
	BOOL ErrState;

	BYTE Reserved[16];	// <Unused> File padding
} USMC_Info;


#ifdef __cplusplus	// C++
extern "C" {
// 
//	 The USMC_Init function initializes driver and returns devices information
//
USMCDLL_API DWORD USMC_Init( USMC_Devices &Str);		// OUT - Array of structures describing all divices (may be NULL) MUST NOT be deleted

//
//	The USMC_GetState function returns structure representing current state of device
//
USMCDLL_API DWORD USMC_GetState( DWORD Device,				// IN - Device number
								 USMC_State &Str			// OUT - Structure representing divice state
								 );
//
//	The USMC_SaveParametersToFlash function saves current parameters of controller in static memory
//	so thay can be loaded at start up time
//
USMCDLL_API DWORD USMC_SaveParametersToFlash( DWORD Device	// IN - Device number
											 );
//
//	The USMC_SetCurrentPosition function sets current position of controller
//
USMCDLL_API DWORD USMC_SetCurrentPosition( DWORD Device,	// IN - Device number
										   int Position		// IN - New position
										   );
//
//	The USMC_GetMode function returns USMC_Mode structure last sent to device
//
USMCDLL_API DWORD USMC_GetMode( DWORD Device,				// IN - Device number
							    USMC_Mode &Str				// OUT - Structure representing some of divice parameters
								);
//
//	The USMC_SetMode function sets some of device parameters
//
USMCDLL_API DWORD USMC_SetMode( DWORD Device,				// IN - Device number
							    USMC_Mode &Str				// IN/OUT Structure representing some of divice parameters
								);
//
//	The USMC_GetParameters function returns USMC_Parameters structure last sent to device
//
USMCDLL_API DWORD USMC_GetParameters( DWORD Device,			// IN - Device number
									  USMC_Parameters &Str	// OUT - Structure representing some of divice parameters
									  );
//
//	The USMC_SetParameters function sets some of divice parameters
//
USMCDLL_API DWORD USMC_SetParameters( DWORD Device,			// IN - Device number
									  USMC_Parameters &Str	// IN/OUT Structure representing some of divice parameters
									  );
//
//	The USMC_GetStartParameters function returns USMC_StartParameters structure last sent to device
//
USMCDLL_API DWORD USMC_GetStartParameters( DWORD Device,	// IN - Device number
										   USMC_StartParameters &Str	// OUT - Structure representing start function parameters
										   );
//
//	The USMC_Start function sets start parameters and starts motion
//
USMCDLL_API DWORD USMC_Start( DWORD Device,					// IN - Device number
							  int DestPos,					// IN - Destination position
							  float &Speed,					// IN/OUT - Speed of rotation
							  USMC_StartParameters &Str		// IN/OUT - Structure representing start function parameters
							  );
//
//	The USMC_Stop function stops device
//
USMCDLL_API DWORD USMC_Stop( DWORD Device					// IN - Device number
							 );
//
//	The USMC_GetLastErr function return string representing last error
//
USMCDLL_API void USMC_GetLastErr( char *str,				// OUT - String buffer
								  size_t len				// IN - Lenght of that string buffer in bytes
								  );

//
//	The USMC_GetDllVersion function returnes version values of USMCDLL.dll
//
USMCDLL_API void USMC_GetDllVersion( DWORD &dwHighVersion,	// OUT - High Version Value
									 DWORD &dwLowVersion);	// OUT - Low Version Value

//
//	The USMC_Close function closes virtual driver window "microsmc.exe"
//
USMCDLL_API DWORD USMC_Close( void );

//
//	The USMC_RestoreCurPos function checks AReset bit and if it is TRUE
//  restores previous CurPos value
//
USMCDLL_API DWORD USMC_RestoreCurPos(DWORD Device					// IN - Device number
									 );

// New For Firmware Version 2.4.1.0 (0x2410)

//
//	The USMC_GetEncoderState function returns structure representing current position of encoder
//
USMCDLL_API DWORD USMC_GetEncoderState( DWORD Device,			// IN - Device number
										USMC_EncoderState &Str	// IN/OUT Structure containing encoder state
										);

// ~New For Firmware Version 2.4.1.0 (0x2410)

}	// extern "C"
#else	// C
DWORD USMC_Init( USMC_Devices *Str );
DWORD USMC_GetState( DWORD Device, USMC_State *Str );
DWORD USMC_SaveParametersToFlash( DWORD Device );
DWORD USMC_SetCurrentPosition( DWORD Device, int Position );
DWORD USMC_GetMode( DWORD Device, USMC_Mode *Str );
DWORD USMC_SetMode( DWORD Device, USMC_Mode *Str );
DWORD USMC_GetParameters( DWORD Device, USMC_Parameters *Str );
DWORD USMC_SetParameters( DWORD Device, USMC_Parameters *Str );
DWORD USMC_GetStartParameters( DWORD Device, USMC_StartParameters *Str );
DWORD USMC_Start( DWORD Device, int DestPos, float *Speed, USMC_StartParameters *Str);
DWORD USMC_Stop( DWORD Device );
void  USMC_GetLastErr( char *str, size_t len );
DWORD USMC_Close( void );
DWORD USMC_GetEncoderState( DWORD Device, USMC_EncoderState *Str);
#endif	// C/C++

