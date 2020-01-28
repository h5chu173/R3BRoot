// ------------------------------------------------------------
// -----                  R3BGlobalAnalysis               -----
// -----            Created 13-04-2016 by M.Heil          -----
// -----            Modified by H.Schulte June 2019       -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BGLOBALANALYSIS
#define R3BGLOBALANALYSIS
#define N_PLANE_MAX_TOFD 4
#define N_PADDLE_MAX_TOFD 50
#define N_PADDLE_MAX_PTOF 100
#define N_PSPX 4
#define N_FIBER_PLOT 1050 // range to plot

#include "FairTask.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <array>

#include "TClonesArray.h"
#include "TMath.h"
#include "TH3F.h"
#include <cstdlib>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;
 
/**
 * This taks reads all detector data items and plots histograms 
 * for online checks. 
 */
class R3BGlobalAnalysis : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BGlobalAnalysis();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BGlobalAnalysis(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BGlobalAnalysis();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    // begin HS_I
    /**
     * A method to extract the cross talk of the fiber detectors and save it with dedicated histograms.
     * Is called from within the Exec function after processing one event's hits.
     */
    virtual void PutHitsIntoHists(std::vector<double>& x_Anode_y_Anode_Tot_MA, /*const int fNEvents,*/ int maxAnodesInRow, int maxAnodeNr, TH2F* hist2DAnodeHitMap, const char* detName, const int ifibcount/*, bool stopDueToTimeDiff*/);

    /**
     *
     *
     */
    virtual void CheckTheEventsCrossTalk(TH2F* hist2DAnodeHitMap, int& CrossTalkLabel, const int ifibcount);

    /**
     *
     *
     */
    virtual void GetHistsCrossTalkLabel1(TH2F* hist2DAnodeHitMap, const int ifibcount, const char* detName);

    /**
     *
     *
     */
    virtual void GetHistsCrossTalkLabel2(TH2F* hist2DAnodeHitMap, const int ifibcount, const char* detName, TH1F* hist1DMAPMTTimeDistEvent, TH1F* hist1DMAPMTTimeDiffsEvent); 
    
    /**
     *
     *
     */
    virtual void GetHistsCrossTalkLabel4(TH2F* hist2DAnodeHitMap, const int ifibcount, const char* detName, TH1F* hist1DMAPMTTimeDistEvent, TH1F* hist1DMAPMTTimeDiffsEvent); 
    // end HS_I

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger)
    {
        fTrigger = trigger;
    }
    
    /**
     * Methods for setting number of planes and paddles
     */
    inline void SetNofModules(Int_t planes, Int_t ppp)
    {
        fNofPlanes   = planes;
        fPaddlesPerPlane  = ppp;
    }

//   virtual void SetParContainers();

    
 
  private:
    std::vector<TClonesArray *> fMappedItems;
    std::vector<TClonesArray *> fCalItems;
    std::vector<TClonesArray *> fHitItems;

    enum DetectorInstances
    {
      DET_AMS,
      DET_CALIFA,
      DET_FI_FIRST,
      DET_FI1A = DET_FI_FIRST,
      DET_FI1B,
      DET_FI2A,
      DET_FI2B,
      DET_FI3A,
      DET_FI3B,
      DET_FI4,
      DET_FI5,
      DET_FI6,
      DET_FI7,
      DET_FI8,
      DET_FI9,
      DET_FI10,
      DET_FI11,
      DET_FI12,
      DET_FI13,
      DET_FI_LAST = DET_FI13,
      DET_L3T,
      DET_LOS,
      DET_MUSIC,
      DET_NEULAND,
      DET_PSPX,
      DET_PTOF,
      DET_ROLU,
      DET_SCI8,
      DET_STRAW,
      DET_TOFD,
      DET_MAX
    };
 
#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

    const char *fDetectorNames[DET_MAX + 1] =
    {
      "Ams",
      "Califa",
      "Fi1a",
      "Fi1b",
      "Fi2a",
      "Fi2b",
      "Fi3a",
      "Fi3b",
      "Fi4",
      "Fi5",
      "Fi6",
      "Fi7",
      "Fi8",
      "Fi9",
      "Fi10",
      "Fi11",
      "Fi12",
      "Fi13",
      "L3t",
      "Los",
      "Music",
      "Neuland",
      "Pspx",
      "Ptof",
      "Rolu",
      "Sci8",
      "Straw",
      "Tofd",
      NULL
   };

    // If FiberI is present or not:
    Int_t  ifibdet;           
    // Number of fibers per detector
    Double_t n_fiber[NOF_FIB_DET]={256.,256.,256.,256.,512.,512.,2048.,2048.,1024.,512.,512.,512.,1024.,1024.};    
    
   
	// check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header;                     /**< Event header. */
    Int_t fTrigger;                             /**< Trigger value. */
    Double_t fClockFreq;     /**< Clock cycle in [ns]. */
    UInt_t fNofPlanes;  
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */    

    Int_t fNEvents = 0;         /**< Event counter. */
    Int_t fFibEvents = 0;         /**< Event counter. */
    Double_t max_values[NOF_FIB_DET][2048];
    Int_t FibMax[NOF_FIB_DET];   
     
    TH2F *fhChargeLosTofD;
    TH2F *fh_los_pos;
    
    TH1F *fh_channels_Fib[NOF_FIB_DET];
    TH1F *fh_fibers_Fib[NOF_FIB_DET];
    TH1F *fh_fiber_Fib[NOF_FIB_DET];
    TH1F *fh_mult_Fib[NOF_FIB_DET];
    TH2F *fh_Fib_ToF[NOF_FIB_DET];
    TH1F *fh_xpos_Fib[NOF_FIB_DET];
    TH1F *fh_ypos_Fib[NOF_FIB_DET];
    
    TH2F *fh_time_Fib[NOF_FIB_DET];
    TH2F *fh_multihit_m_Fib[NOF_FIB_DET];   
    TH2F *fh_multihit_s_Fib[NOF_FIB_DET];
    TH2F *fh_ToT_m_Fib[NOF_FIB_DET];
    TH2F *fh_ToT_s_Fib[NOF_FIB_DET];
    TH2F *fh_Fib_vs_Events[NOF_FIB_DET];
    TH2F *fh_Fibs_vs_Events[NOF_FIB_DET];
    TH2F *fh_Fibs_vs_Tofd[NOF_FIB_DET];

    TH2F *fh_Fib_vs_Fib[NOF_FIB_DET][NOF_FIB_DET];
    TH2F *fh_Fib_dx[NOF_FIB_DET][NOF_FIB_DET];
    
    TH2F *fh_Cave_position;
    
    TH2F* fh_tofd_pos;
    TH1F* fh_tofd_charge;    
    
    TH1F *fh_ptof_channels; 
    TH1F *fh_ptof_channels_cut; 
    TH1F *fh_ptof_test1;   
    TH1F *fh_ptof_test2;   
    TH1F* fh_ptof_TotPm1[N_PADDLE_MAX_PTOF]; 
    TH1F* fh_ptof_TotPm2[N_PADDLE_MAX_PTOF]; 

    
    TH1F *fh_pspx_strips_psp[N_PSPX];
    TH1F *fh_pspx_energy_psp[N_PSPX];
    TH1F *fh_pspx_multiplicity_psp[N_PSPX];
    
    TH2F *fh_pspx_pos1_strips;
    TH2F *fh_pspx_pos2_strips;
    TH2F *fh_pspx_pos1_energy;
    TH2F *fh_pspx_pos2_energy;
    
    TH2F *fh_pspx_cor_x_strips;
    TH2F *fh_pspx_cor_y_strips;
    TH2F *fh_pspx_cor_x_energy;
    TH2F *fh_pspx_cor_y_energy;

    // begin HS_I
    TH2F* hist2DSmallEntrVSBigEntr[NOF_FIB_DET];
    TH2F* hist2DSmallEntrVSBigEntr_CTL4[NOF_FIB_DET];
    TH1F* hist1DAnodeNeighbours[NOF_FIB_DET];
    TH1F* hist1DAnodeNeighboursVsEventNr[NOF_FIB_DET];
    TH1F* hist1DAnodeNeighboursToTSum[NOF_FIB_DET];  
    TH2F* hist2DAnodeNeighboursToTSumVsNeighbourNr[NOF_FIB_DET];
    TH1F* hist1DMAPMTTimeDistAll[NOF_FIB_DET];
    TH1F* hist1DTimeDiffsAll[NOF_FIB_DET];
    TH2F* hist2DInvestigateMaxEntriesI[NOF_FIB_DET];
    //TH3F* hist3DInvestigateMaxEntriesII[NOF_FIB_DET];
    TH2F* hist2DMaxAnode_CTL2[NOF_FIB_DET];
    TH1F* hist1DTimeDiffs_CTL2[NOF_FIB_DET];
    TH2F* hist2DMaxAnode_CTL4[NOF_FIB_DET];
    TH1F* hist1DTimeDiffs_CTL4[NOF_FIB_DET];
    TH2F* hist2DMaxAnode_CTL1[NOF_FIB_DET];
    // end HS_I

    
 
  
  public:
    ClassDef(R3BGlobalAnalysis, 1)
};

#endif
