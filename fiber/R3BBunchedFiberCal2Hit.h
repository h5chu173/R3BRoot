/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BBUNCHEDFIBERCAL2HIT
#define R3BBUNCHEDFIBERCAL2HIT

#include "FairTask.h"
#include <R3BTCalEngine.h>

#include <list>

class TH1F;
class TH2F;

class R3BBunchedFiberCalData;
class R3BBunchedFiberHitPar;
class R3BBunchedFiberHitModulePar;

/**
 * Transforms bunched fiber Cal level data to Hit level.
 *
 * A note about the UCESB unpacker:
 * All channels in a full plane of fibers should be mapped in sequence. For
 * example, Fib4 is built out of two "sub-detectors" (left and right side),
 * where each subdet has 256 CTDC/MAPMT 2 TAMEX2/SPMT channels, which should be
 * mapped as:
 *  FIB4_MAPMT1 .. FIB4_MAPMT512
 *  FIB4_SPMT1 .. FIB4_SPMT4
 * This class will then figure out how to map the actual fibers (except for
 * non-ambiguous sorting errors).
 */
class R3BBunchedFiberCal2Hit : public FairTask
{
  public:
    enum Direction
    {
        HORIZONTAL,
        VERTICAL
    };
    struct ToT
    {
        ToT(R3BBunchedFiberCalData const*, R3BBunchedFiberCalData const*, Double_t);
        R3BBunchedFiberCalData const* lead;
        R3BBunchedFiberCalData const* trail;
        Double_t tot;
    };
    struct Channel
    {
        std::list<R3BBunchedFiberCalData const*> lead_list;
        std::list<ToT> tot_list;
    };

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BBunchedFiberCal2Hit(const char*,
                           Int_t,
                           enum R3BTCalEngine::CTDCVariant,
                           Direction,
                           UInt_t,
                           UInt_t,
                           UInt_t,
                           Bool_t,
                           Bool_t,
                           Bool_t);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BBunchedFiberCal2Hit();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
    virtual InitStatus ReInit();

    virtual void SetParContainers();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t*);

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

    R3BBunchedFiberHitModulePar* GetModuleParAt(Int_t fiber);

    /**
     * Getting all fibers correct is difficult, this function lets us fix that
     * per detector.
     */
    virtual UInt_t FixMistake(UInt_t) = 0;

    // begin HS_II_2
    /**
     * A new method to use the TofWall data for different ways to filter the fiber detectors' based on the one from A.Kelic-Heil
     * Is called in the Exec routine.
     */
    //virtual void UseTofWall(size_t& cal_num, std::vector<double>& FibHitToFs, double arrTimeStampsPMTs[2]);
    virtual void UseTofWall_new(std::vector<double>& FibDataVec);
    
    // end HS_II_2
    //
  private:
    TString fName;
    Int_t fnEvents;
    Int_t maxevent;
    Int_t fnEventsfill = 0;  // HS_II_2

    Double_t tsync_temp[2048]={0};
    Double_t gain_temp[2048]={10};
    Bool_t tofdin;

    double fClockFreq;
    Direction fDirection;
    UInt_t fSubNum;
    UInt_t fChPerSub[2];
    Bool_t fIsCalibrator;
    Bool_t fIsGain;  // HS_II_2
    Bool_t fIsTsynch; // HS_II_2
    TClonesArray* fCalItems;
    TClonesArray* fHitItems;
    TClonesArray* fTofdHitItems;  // HS_II_2
    R3BBunchedFiberHitPar* fCalPar; /**< Parameter container. */
    R3BBunchedFiberHitPar* fHitPar; /**< Hit parameter container. */
    Int_t fNofHitPars;              /**< Number of modules in parameter file. */
    Int_t fNofHitItems;
    Bool_t involveToFWall;  // HS_II_2
    Bool_t testTof;  // HS_II_2
    Bool_t oldHists; // HS_II_2
    //R3BTofdHitPar* fTofdHitPar;
    //int fNumOfTofdHitPars;
    // [0=MAPMT,1=SPMT][Channel].
    std::vector<Channel> fChannelArray[2];

    // histograms for gain matching
    TH2F* fh_ToT_MA_Fib;
    TH2F* fh_ToT_Single_Fib;
    TH2F* fh_ToT_s_Fib[4];
    TH2F* fh_ToT_ToT;
    TH2F* fh_dt_Fib;

    // histograms concerning ToF
    //HS_II_2
    TH2F* fh_ToF_FibToTofwall;
    TH2F* fh_TestTof;
    // end HS_II_2

  public:
    ClassDef(R3BBunchedFiberCal2Hit, 3)
};

#endif
