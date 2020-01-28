// ------------------------------------------------------------
// -----                  R3BGlobalAnalysis               -----
// -----          Created April 13th 2016 by M.Heil       -----
// -----          Modiefied by H.Schulte June 2019        -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow  
 * to test the detectors online
 *  
 */


#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"
#include "R3BLosHitData.h"

#include "R3BGlobalAnalysis.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BPaddleCalData.h"

#include "R3BPspxMappedData.h"
#include "R3BPspxCalData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"

#include "TClonesArray.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <array>
#include "TMath.h"
#include <cstdlib>
#include <ctime>
#include <TRandom3.h>
#include <TRandomGen.h>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

// begin HS_I

void SetUpTheCoordinates(std::vector<double>& x_Anode_y_Anode_Tot_MA, int maxAnodesInRow, int maxAnodeNr, const char* detName);                                                                 
//void PrepareDataForPlotting(std::vector<double>& x_Anode_y_Anode_Tot_MA, std::vector<int>& x_FibCoords, std::vector<int>& y_FibCoords, int iFib, double ToT_MA, int maxAnodeNr);
//void PutHitsIntoHists(std::vector<double>& x_Anode_y_Anode_Tot_MA, const int fNEvents, int maxAnodesInRow, int maxAnodeNr, TH2F* hist2DAnodeHitMap, const char* detName, int& CrossTalkLabel);
//void CheckTheEventsCrossTalk(TH2F* hist2DAnodeHitMap, int& CrossTalkLabel, const int ifibcount);
//void GetHistsCrossTalkLabel1([> put in pointer for hists <]);
//void GetHistsCrossTalkLabel2([> put in pointer for hists <]);
void GetAverageToTMAPerNeighbour(int fNEvents, int neighbours, double ToTsumAnodeNeighboursOfOneMax, const char* detName); 

// end HS_I

R3BGlobalAnalysis::R3BGlobalAnalysis()
    : R3BGlobalAnalysis("GlobalAnalysis", 1)
{}

R3BGlobalAnalysis::R3BGlobalAnalysis(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fNofPlanes(N_PLANE_MAX_TOFD)  
    , fPaddlesPerPlane(N_PADDLE_MAX_TOFD) 
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvents(0)
    , fFibEvents(0)
{
}

R3BGlobalAnalysis::~R3BGlobalAnalysis()
{
   for(int i = 0; i < NOF_FIB_DET; i++) {   
     delete fh_channels_Fib[i];
     delete fh_fibers_Fib[i];
     delete fh_fiber_Fib[i];
     delete fh_mult_Fib[i];
     delete fh_time_Fib[i];
     delete fh_xpos_Fib[i];
     delete fh_ypos_Fib[i];
     delete fh_multihit_m_Fib[i];   
     delete fh_multihit_s_Fib[i];
     delete fh_ToT_m_Fib[i];
     delete fh_ToT_s_Fib[i];
   } 
}

InitStatus R3BGlobalAnalysis::Init()
{

  // begin HS_I

  /*
    ofstream CheckInit;
    CheckInit.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckInit << " Die Init wird aufgerufen! \n";
    CheckInit.close();
    */
  // end HS_I

    // Initialize random number:
    std::srand(std::time(0)); //use current time as seed for random generator


    LOG(INFO) << "R3BGlobalAnalysis::Init " << FairLogger::endl;

    // try to get a handle on the EventHeader. EventHeader may not be 
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    FairRunOnline *run = FairRunOnline::Instance();

        // Get objects for detectors on all levels
        assert(DET_MAX + 1 == sizeof(fDetectorNames)/sizeof(fDetectorNames[0]));
        printf("Have %d fiber detectors.\n", NOF_FIB_DET);
        for (int det = 0; det < DET_MAX; det++)
        {
/*
                fMappedItems.push_back((TClonesArray *)mgr->GetObject(Form("%sMapped", fDetectorNames[det])));
                if (NULL == fMappedItems.at(det)) {
                  printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
                }
                fCalItems.push_back((TClonesArray *)mgr->GetObject(Form("%sCal", fDetectorNames[det])));
                fHitItems.push_back((TClonesArray *)mgr->GetObject(Form("%sHit", fDetectorNames[det])));
*/
                fMappedItems.push_back((TClonesArray *)mgr->GetObject(Form("%sMapped", fDetectorNames[det])));
                if (NULL == fMappedItems.at(det)) {
                  printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
                }
                fCalItems.push_back((TClonesArray *)mgr->GetObject(Form("%sCal", fDetectorNames[det])));
                if (NULL == fCalItems.at(det)) {
                  printf("Could not find Cal data for '%s'.\n", fDetectorNames[det]);
                }
                fHitItems.push_back((TClonesArray *)mgr->GetObject(Form("%sHit", fDetectorNames[det])));
                if (NULL == fHitItems.at(det)) {
                  printf("Could not find hit data for '%s'.\n", fDetectorNames[det]);
                }
 
        }


    //------------------------------------------------------------------------ 
    // create histograms of all detectors  
    //------------------------------------------------------------------------ 

    char strNameC[255];
    sprintf(strNameC, "Cave_C_position");
    fh_Cave_position = new TH2F(strNameC, "", 2100, -100., 2000., 1000, -500., 500.);              

    if(fHitItems.at(DET_LOS) || fCalItems.at(DET_LOS)){
        char strName[255];
        sprintf(strName, "LOS_X_vs_Y_MCFD");
        //     fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);
        fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);              
            
        sprintf(strName, "LosQ_vs_TofdQ");
        fhChargeLosTofD = new TH2F(strName, "", 1000, 0., 80., 1000, 0., 80.);              
        fhChargeLosTofD->GetYaxis()->SetTitle("Charge LOS");
        fhChargeLosTofD->GetXaxis()->SetTitle("Charge ToFD");                                             
    }

    //-----------------------------------------------------------------------
    // Fiber Detectors 1-NOF_FIB_DET

    char canvName[255];

    for(Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++){

        std::stringstream fibhitobj;
        std::stringstream fibmapobj;
        std::string temp;     
        std::stringstream FileNameTemp;   

        std::stringstream histName1,histName2,histName22,histName3,histName4,histName5,histName6,histName33,
                  histName7,histName8,histName9,histName10,histName11,histName12,histName13,histName14;
        std::stringstream histTitle1,histTitle2,histTitle22,histTitle3,histTitle4,histTitle5,histTitle6,histTitle33,
                  histTitle7,histTitle8,histTitle9,histTitle10,histTitle11,histTitle12,histTitle13,histTitle14;              

        if(fCalItems.at(DET_FI_FIRST + ifibcount) || fHitItems.at(DET_FI_FIRST + ifibcount)) {

            const char* chistName;
            const char* chistTitle;
            const char* detName;
            const char* detName2;
            const char* FileName;
            std::string tempTitle;
            std::string tempName;
            std::stringstream tempFibName;
            std::string tempFibNames;
            std::stringstream tempCanvName;
            
            detName = fDetectorNames[DET_FI_FIRST + ifibcount];

            cout << "I am creating canvas " << detName <<endl;

            // begin HS_I
            // additional histograms to investigate the fiber detectors' cross talk
            hist2DSmallEntrVSBigEntr[ifibcount] = new TH2F(Form("SecondEntryOfMAPMmapVSMax_%s", detName), Form("TOT at MAPM from possible fibrous cross talk at %s; second peak vs max", detName), 40, 0, 200, 40, 0, 200); 
            hist2DSmallEntrVSBigEntr[ifibcount]->GetXaxis()->SetTitle("ToT_MA of max"); 
            hist2DSmallEntrVSBigEntr[ifibcount]->GetYaxis()->SetTitle("ToT_MA of second peak"); 

            hist2DSmallEntrVSBigEntr_CTL4[ifibcount] = new TH2F(Form("SecondEntryOfMAPMmapVSMax_CTL4_%s", detName), Form("TOT at MAPM from possible fibrous cross talk at %s; second peak vs max", detName), 40, 0, 200, 40, 0, 200); 
            hist2DSmallEntrVSBigEntr_CTL4[ifibcount]->GetXaxis()->SetTitle("ToT_MA of max"); 
            hist2DSmallEntrVSBigEntr_CTL4[ifibcount]->GetYaxis()->SetTitle("ToT_MA of second peak"); 

            hist1DAnodeNeighbours[ifibcount] = new TH1F(Form("AnodeNeighboursofOneMax_%s", detName), Form("Number of neighbouring entries on the MAPM anode around one max; det: %s",detName),10, 0, 10); 
            hist1DAnodeNeighbours[ifibcount]->GetXaxis()->SetTitle("Number of hit neighbouring Anodes of max");

            //hist1DAnodeNeighboursVsEventNr[ifibcount] = new TH1F(); 
            //hist1DAnodeNeighboursVsEventNr[ifibcount]->SetNameTitle(Form("AnodeNeighboursofOneMaxVsEventNr_%s", detName), Form("Number of neighbouring entries on the MAPM anode around one max; det: %s", detName));
            //hist1DAnodeNeighboursVsEventNr[ifibcount]->GetXaxis()->SetTitle("Event number"); 
            //hist1DAnodeNeighboursVsEventNr[ifibcount]->GetYaxis()->SetTitle("Number of hit neighbouring Anodes of max"); 

            hist1DAnodeNeighboursToTSum[ifibcount] = new TH1F(Form("ToTsumAnodeNeighboursOfOneMax_%s",detName), Form("Sum of MAPM ToT of neighbouring entries around one max; det: %s", detName), 300, 0, 300  ); 
            hist1DAnodeNeighboursToTSum[ifibcount]->GetXaxis()->SetTitle("Tot_MA sum of hit neighbouring Anodes of max");

            hist2DAnodeNeighboursToTSumVsNeighbourNr[ifibcount] = new TH2F(Form("ToTsumAnodeNeighboursOfOneMaxVsNeighbourNr_%s", detName), Form("Sum of MAPM ToT of neighbouring entries around one max normalized by the number of neighbours vs corresponding Nr of neighbours; det: %s", detName), 10, 0, 10, 2000, 0, 300);
            hist2DAnodeNeighboursToTSumVsNeighbourNr[ifibcount]->GetXaxis()->SetTitle("Number of hit neighbouring Anodes of max"); 
            hist2DAnodeNeighboursToTSumVsNeighbourNr[ifibcount]->GetYaxis()->SetTitle("Sum of ToT_MA of hit neighbouring Anodes of max"); 
            hist2DAnodeNeighboursToTSumVsNeighbourNr[ifibcount]->GetZaxis()->SetTitle("Number of entries");


            hist1DMAPMTTimeDistAll[ifibcount] = new TH1F(Form("MAPMTTimeDist_%s", detName), Form("Time distribution of hits on MAPMTs of det %s", detName), 330000, 0, 16500);
            hist1DMAPMTTimeDistAll[ifibcount]->GetXaxis()->SetTitle("Time stamp (tMAPMT) [ns]");
            hist1DMAPMTTimeDistAll[ifibcount]->GetYaxis()->SetTitle("Entries");

            hist1DTimeDiffsAll[ifibcount] = new TH1F(Form("MAPMTTimeDiffs_%s", detName), Form("Differences of time stamps from Hits to time stamp from first hit of the same event in det %s", detName), 30000, 0, 3000);
            hist1DTimeDiffsAll[ifibcount]->GetXaxis()->SetTitle("Time difference (tMAPMT) [ns]");
            hist1DTimeDiffsAll[ifibcount]->GetYaxis()->SetTitle("Entries");

            hist2DMaxAnode_CTL2[ifibcount] = new TH2F(Form("MaxAnode_CTL2_%s", detName), Form("The anode's coordiantes of the maximum hit from each CTL2 event, det %s", detName), 16, 0, 16, 16, 0, 16);
            hist2DMaxAnode_CTL2[ifibcount]->GetXaxis()->SetTitle("anode's x coord");
            hist2DMaxAnode_CTL2[ifibcount]->GetYaxis()->SetTitle("anode's y coord");
            hist2DMaxAnode_CTL2[ifibcount]->GetZaxis()->SetTitle("entries");

            hist2DMaxAnode_CTL4[ifibcount] = new TH2F(Form("MaxAnode_CTL4_%s", detName), Form("The anode's coordiantes of the maximum hit from each CTL4 event, det %s", detName), 16, 0, 16, 16, 0, 16);
            hist2DMaxAnode_CTL4[ifibcount]->GetXaxis()->SetTitle("anode's x coord");
            hist2DMaxAnode_CTL4[ifibcount]->GetYaxis()->SetTitle("anode's y coord");
            hist2DMaxAnode_CTL4[ifibcount]->GetZaxis()->SetTitle("entries");

            hist2DMaxAnode_CTL1[ifibcount] = new TH2F(Form("MaxAnode_CTL1_%s", detName), Form("The anode's coordiantes of the maximum hit from each CTL1 event, det %s", detName), 16, 0, 16, 16, 0, 16);
            hist2DMaxAnode_CTL1[ifibcount]->GetXaxis()->SetTitle("anode's x coord");
            hist2DMaxAnode_CTL1[ifibcount]->GetYaxis()->SetTitle("anode's y coord");
            hist2DMaxAnode_CTL1[ifibcount]->GetZaxis()->SetTitle("entries");

            hist1DTimeDiffs_CTL2[ifibcount] = new TH1F(Form("MAPMTTimeDiffs_CTL2_%s", detName), Form("Differences of time stamps from Hits to time stamp from first hit of the same event in det %s", detName), 1000, 0, 100);
            hist1DTimeDiffs_CTL2[ifibcount]->GetXaxis()->SetTitle("Time difference (tMAPMT) [ns]");
            hist1DTimeDiffs_CTL2[ifibcount]->GetYaxis()->SetTitle("Entries");

            //hist1DMAPMTTimeDist_CTL4[ifibcount] = new TH1F(Form("MAPMTTimeDist_CTL4_%s", detName), Form("Time distribution of hits on MAPMTs of det %s", detName), 330000, 0, 16500);
            //hist1DMAPMTTimeDist_CTL4[ifibcount]->GetXaxis()->SetTitle("Time stamp (tMAPMT) [ns]");
            //hist1DMAPMTTimeDist_CTL4[ifibcount]->GetYaxis()->SetTitle("Entries");

            hist1DTimeDiffs_CTL4[ifibcount] = new TH1F(Form("MAPMTTimeDiffs_CTL4_%s", detName), Form("Differences of time stamps from Hits to time stamp from first hit of the same event in det %s", detName), 1000, 0, 100);
            hist1DTimeDiffs_CTL4[ifibcount]->GetXaxis()->SetTitle("Time difference (tMAPMT) [ns]");
            hist1DTimeDiffs_CTL4[ifibcount]->GetYaxis()->SetTitle("Entries");

            hist2DInvestigateMaxEntriesI[ifibcount] = new TH2F(Form("InvestigateMaxEntriesI_%s", detName), " ", 5, 0, 5, 500, 0, 500);
            hist2DInvestigateMaxEntriesI[ifibcount]->GetXaxis()->SetTitle("Label");
            hist2DInvestigateMaxEntriesI[ifibcount]->GetYaxis()->SetTitle("ToTMA of max of Event");

            //hist3DInvestigateMaxEntriesII[ifibcount] = new TH3F(Form("InvestigateMaxEntriesII_%s", detName), " ", 5, 0, 5, 500, 0, 500, 20, 0, 20);
            //hist3DInvestigateMaxEntriesII[ifibcount]->GetXaxis()->SetTitle("Label");
            //hist3DInvestigateMaxEntriesII[ifibcount]->GetYaxis()->SetTitle("ToTMA of max of Event");
            //hist3DInvestigateMaxEntriesII[ifibcount]->GetZaxis()->SetTitle("Entries in 'HitMap'");

            // end HS_I

            // Channels:   
			fh_channels_Fib[ifibcount] = new TH1F(Form("%s_channels",detName), Form("%s channels",detName), 520, 0., 520.);
            fh_channels_Fib[ifibcount]->GetXaxis()->SetTitle("Channel number");
            fh_channels_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Fibers:
			fh_fibers_Fib[ifibcount] = new TH1F(Form("%s_fibers",detName), Form("%s fibers",detName), N_FIBER_PLOT, 0., N_FIBER_PLOT); 
            fh_fibers_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_fibers_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Fiber:
			fh_fiber_Fib[ifibcount] = new TH1F(Form("%s_fiber",detName), Form("%s fiber",detName), N_FIBER_PLOT, 0., N_FIBER_PLOT); 
            fh_fiber_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_fiber_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Multiplicity (number of hit fibers):
			fh_mult_Fib[ifibcount] = new TH1F(Form("%s_mult",detName), Form("%s # of fibers",detName), 100, 0., 100.);	   
            fh_mult_Fib[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Multihit MAPMT:
			fh_multihit_m_Fib[ifibcount] = new TH2F(Form("%s_multihit_m",detName), Form("%s multihits MAPMT",detName), 520, 0., 520., 20, 0., 20.);
            fh_multihit_m_Fib[ifibcount]->GetXaxis()->SetTitle("MAPMT channel");
            fh_multihit_m_Fib[ifibcount]->GetYaxis()->SetTitle("Multihit");

            // Multihit SAPMT:
			fh_multihit_s_Fib[ifibcount] = new TH2F(Form("%s_multihit_s",detName), Form("%s multihits single PMT",detName), 16, 0., 16., 20, 0., 20.);   	   
            fh_multihit_s_Fib[ifibcount]->GetXaxis()->SetTitle("SAPMT channel");
            fh_multihit_s_Fib[ifibcount]->GetYaxis()->SetTitle("Multihit");

            // ToT MAPMT:  
			fh_ToT_m_Fib[ifibcount] = new TH2F(Form("%s_tot_m",detName), Form("%s ToT of MAPMT",detName), N_FIBER_PLOT, 0., N_FIBER_PLOT, 400, 0., 20.);   	   
            fh_ToT_m_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_ToT_m_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // ToT SAPMT:  
			fh_ToT_s_Fib[ifibcount] = new TH2F(Form("%s_tot_s",detName), Form("%s ToT of single PMT",detName), N_FIBER_PLOT, 0., N_FIBER_PLOT, 400, 0., 40.);   	   
            fh_ToT_s_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_ToT_s_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // Time of fiber: 
			fh_time_Fib[ifibcount]= new TH2F(Form("%s_TimevsFiber",detName), Form("%s Time vs Fiber",detName), N_FIBER_PLOT, 0., N_FIBER_PLOT, 20000, -1024., 1024.);
            fh_time_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_time_Fib[ifibcount]->GetYaxis()->SetTitle("tMAPMT-tSPMT");

            // ToF LOS -> Fiber:
			fh_Fib_ToF[ifibcount] = new TH2F(Form("%s_tof",detName), Form("%s ToF LOS to Fiber",detName), N_FIBER_PLOT, 0., N_FIBER_PLOT,10000, -5000., 5000.);	   
            fh_Fib_ToF[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF[ifibcount]->GetXaxis()->SetTitle("Fiber number");

            // Not-calibrated position:
			fh_xpos_Fib[ifibcount] = new TH1F(Form("%s_pos",detName), Form("%s Not-calibrated position",detName), 6000, -1500., 1500.);	   
            fh_xpos_Fib[ifibcount]->GetXaxis()->SetTitle("x-Position");
            fh_xpos_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // hit fiber number vs. event number:
			fh_Fib_vs_Events[ifibcount] = new TH2F(Form("%s_fib_vs_event",detName), Form("%s Fiber # vs. Event #",detName), 1000,0,5e6,1100, 0., 1100.);	   
			fh_Fib_vs_Events[ifibcount]->GetYaxis()->SetTitle("Fiber number");
			fh_Fib_vs_Events[ifibcount]->GetXaxis()->SetTitle("Event number");

            // hit MA number vs. event number:
			fh_Fibs_vs_Events[ifibcount] = new TH2F(Form("%s_fibs_vs_event",detName), Form("%s Fiber # vs. Event #",detName), 1000,0,5e6,1100, 0., 1100.);	   
			fh_Fibs_vs_Events[ifibcount]->GetYaxis()->SetTitle("Fiber number");
			fh_Fibs_vs_Events[ifibcount]->GetXaxis()->SetTitle("Event number");

            // hit fiber number vs. TofD position:
            fh_Fibs_vs_Tofd[ifibcount] = new TH2F(Form("%s_fib_vs_TofdX",detName), Form("%s Fiber # vs. Tofd x-pos",detName), 200,-100,100,1100, 0., 1100.);    
            fh_Fibs_vs_Tofd[ifibcount]->GetYaxis()->SetTitle("Fiber number");
            fh_Fibs_vs_Tofd[ifibcount]->GetXaxis()->SetTitle("Tofd x-pos number");

            // hit fiber number vs. TofD position:
            for(Int_t j = ifibcount+1; j < NOF_FIB_DET; j++){
				detName2 = fDetectorNames[DET_FI_FIRST + j];
	     	    fh_Fib_vs_Fib[ifibcount][j] = new TH2F(Form("fib%s_vs_fib%s",detName,detName2), 
	     	        Form("Fiber %s vs. Fiber %",detName,detName2), 110,0,1100,110, 0., 1100.);    
                fh_Fib_vs_Fib[ifibcount][j]->GetYaxis()->SetTitle(Form("Fiber%s",detName2));
                fh_Fib_vs_Fib[ifibcount][j]->GetXaxis()->SetTitle(Form("Fiber%s",detName));;
		    }
	
            // dx between fibers vs x
            for(Int_t j = ifibcount+1; j < NOF_FIB_DET; j++){
				detName2 = fDetectorNames[DET_FI_FIRST + j];
	     	    fh_Fib_dx[ifibcount][j] = new TH2F(Form("fib%s_fib%s_dx",detName,detName2), 
	     	        Form("dx of Fiber %s and Fiber %",detName2,detName), 1100,0,1100,1000, -500., 500.);    
                fh_Fib_dx[ifibcount][j]->GetYaxis()->SetTitle("dx");
                fh_Fib_dx[ifibcount][j]->GetXaxis()->SetTitle(Form("Fiber%s",detName));;
		    }

            // Not-calibrated position:
            histName14 << detName << "_ypos";
            tempName=histName14.str();
            chistName=tempName.c_str();
            histTitle14 << detName << " Not-calibrated y-position ";
            tempTitle=histTitle14.str();
            chistTitle=tempTitle.c_str();
            fh_ypos_Fib[ifibcount] = new TH1F(chistName, chistTitle, 6000, -1500., 1500.);     
            fh_ypos_Fib[ifibcount]->GetXaxis()->SetTitle("y-Position");
            fh_ypos_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");
            tempName.clear();
            tempTitle.clear();


        }   // end if(Mapped)
        
    }  // end for(ifibcount)


    //---------------------------------------------------------------------------------------------------
    //TofD detector

    if(fHitItems.at(DET_TOFD) || fCalItems.at(DET_TOFD)){

        char strName1[255];
        char strName2[255];
        sprintf(strName1, "TofD_X_vs_Y");
        sprintf(strName2, "TofD X vs. Y");
        fh_tofd_pos = new TH2F(strName1, strName2, 2000, -10., 10., 2000, -10., 10.);              

        sprintf(strName1, "TofD_Q");
        sprintf(strName2, "TofD Q");
        fh_tofd_charge = new TH1F(strName1, strName2, 200, 0., 20.);              


/*
        for (Int_t j = 0; j < 4; j++)
        {
            char strName1[255];
            sprintf(strName1, "tofd_channels_plane_%d", j+1);
            char strName2[255];
            sprintf(strName2, "Tofd channels plane %d", j+1);        
            fh_tofd_channels[j] = new TH1F(strName1, strName2, 90, -45., 45.);       
            fh_tofd_channels[j]->GetXaxis()->SetTitle("Channel");
            fh_tofd_channels[j]->GetYaxis()->SetTitle("Counts");

            char strName3[255];
            sprintf(strName3, "tofd_ToT_plane_%d", j+1);
            char strName4[255];
            sprintf(strName4, "Tofd ToT plane %d", j+1);        
            fh_tofd_TotPm[j] = new TH2F(strName3, strName4, 90,-45,45,3000,0.,300.);
            fh_tofd_TotPm[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_TotPm[j]->GetYaxis()->SetTitle("ToT / ns");

            char strName7[255];
            sprintf(strName7, "tofd_ToF_plane_%d", j+1);
            char strName8[255];
            sprintf(strName8, "TOFD-LOS ToF plane %d", j+1); 
            fh_tofd_ToF[j] = new TH2F(strName7, strName8, 45,0,45, 1e4, 0, 1e3);
            fh_tofd_ToF[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_ToF[j]->GetYaxis()->SetTitle("ToF / ns");

            char strName9[255];
            sprintf(strName9, "tofd_multihit_plane_%d", j+1);
            char strName10[255];
            sprintf(strName10, "Tofd multihit plane %d", j+1);        
            fh_tofd_multihit[j] = new TH2F(strName9, strName10, 45, 0., 45., 10, 0, 10);         
            fh_tofd_multihit[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_multihit[j]->GetYaxis()->SetTitle("Multihit");


          
            for(Int_t p = 0;p<N_PADDLE_MAX_PTOF; p++){
                char strName11[255];
                sprintf(strName11, "tofd_ToT_vs_pos_paddle_%d_plane_%d", p+1,j+1);
                char strName12[255];
                sprintf(strName12, "Tofd ToT vs. Pos bar %d plane %d", p+1, j+1);        
                fh_tofd_TotvsPos[j][p] = new TH2F(strName11, strName12, 100, -50., 50., 300, 0, 300);        
                fh_tofd_TotvsPos[j][p]->GetXaxis()->SetTitle("position");
                fh_tofd_TotvsPos[j][p]->GetYaxis()->SetTitle("ToT / ns");
            }


            
            if(j<3)
            {
                int jk = j + 1;
                char strName5[255];
                sprintf(strName5, "tofd_ToF_plane_%d_%d", jk, jk+1);
                char strName6[255];
                sprintf(strName6, "Tofd ToF plane %d and %d ", jk, jk+1); 
                fh_tofd_dt[j] = new TH2F(strName5, strName6, 50,0,50,8000, -40., 40);
                fh_tofd_dt[j]->GetXaxis()->SetTitle("Bar number");
                fh_tofd_dt[j]->GetYaxis()->SetTitle("dt / ns");
            }
        }



        char strName11[255];
        sprintf(strName11, "tofd_ToT_vs_pos_paddle_plane");
        char strName12[255];
        sprintf(strName12, "Tofd ToT vs. Pos bar plane");        
        fh_tofd_TotvsPos = new TH2F(strName11, strName12, 1000, -50., 50., 3000, 0, 300);        
        fh_tofd_TotvsPos->GetXaxis()->SetTitle("position");
        fh_tofd_TotvsPos->GetYaxis()->SetTitle("ToT / ns");

*/


    }

    // -------------------------------------------------------------------------


    return kSUCCESS;
}


void R3BGlobalAnalysis::Exec(Option_t* option)
{
  //begin HS_I 
  /*
  if(fNEvents > 200) {
    //printf("The event number is to low! (%d) \n", fNEvents);
    fNEvents++;
    return;
  }
  */
  //end HS_I

    //printf("die Exec wird aufgerufen! \n");
    ofstream CheckExec;
    CheckExec.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckExec << "Die Exec wird aufgerufen! \n";
    CheckExec.close();


    if(fNEvents/1000000.==(int)fNEvents/1000000) cout<<"Events: "<<fNEvents<<flush<<'\r';


    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
      LOG(fatal) << "FairRootManager not found";
    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
    return;

    //----------------------------------------------------------------------
    // LOS detector
    //----------------------------------------------------------------------
    Double_t timeLos[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Double_t LosQ=0.;

    if(fHitItems.at(DET_LOS)){
        auto det = fHitItems.at(DET_LOS);
        Int_t nHits = det->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nHits; ihit++)     
        {
            R3BLosHitData *hitData = (R3BLosHitData*)det->At(ihit);
            timeLos[ihit] = hitData->fTime_ns;

            fh_los_pos->Fill(hitData->fX_cm,hitData->fY_cm); 

            fh_Cave_position->Fill(0.,hitData->fX_cm);

            if (hitData->fZ>LosQ) LosQ=hitData->fZ;
        }
    }


    //----------------------------------------------------------------------
    // TOFD
    //----------------------------------------------------------------------
    Double_t timeTofd[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Double_t xTofd[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Double_t yTofd[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Double_t qTofd[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Double_t TofdQ=0.;
    Double_t TofdX=0.;

    if(fHitItems.at(DET_TOFD))
    {
	    auto det = fHitItems.at(DET_TOFD);
	    Int_t nHits = det->GetEntriesFast();
//	    cout<<"ToFD hits: "<<nHits<<endl;
        if (nHits==1) {
	    for (Int_t ihit = 0; ihit < nHits; ihit++)     
	    {
	        R3BTofdHitData* hit = (R3BTofdHitData*)det->At(ihit);
		    if (!hit) continue; // should not happen
			if (ihit>15) {
				cout<<"Error, more than 16 hits"<<endl;
				continue;
			}
            timeTofd[ihit] = hit->GetTime();
            xTofd[ihit] = hit->GetX();
            yTofd[ihit] = hit->GetY();
            qTofd[ihit] = hit->GetEloss();
            if(qTofd[ihit] > TofdQ) {
				TofdQ = qTofd[ihit];
				TofdX = xTofd[ihit];
			}
//		    cout<<"ToFD: "<<ihit<<" x: "<< xTofd[ihit] << " y: " << yTofd[ihit] 
//		        << " q: "<< qTofd[ihit] << " t: "<< timeTofd[ihit] << endl;

            fh_tofd_pos->Fill(xTofd[ihit],yTofd[ihit]); 
            fh_Cave_position->Fill(1000.,xTofd[ihit]);
            fh_tofd_charge->Fill(qTofd[ihit]); 
		
        }
        
        }

    }
    
    
    //fhChargeLosTofD->Fill(TofdQ,LosQ/2.);


    //----------------------------------------------------------------------
    // Fiber detectors
    //----------------------------------------------------------------------
    Double_t dtime = 0.0/0.0;
    Int_t test=0;
    Int_t counter=0;
    for(Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    { 
        Int_t iFib = 0;  
        FibMax[ifibcount]=0;
        auto detMapped = fMappedItems.at(DET_FI_FIRST + ifibcount);
        auto detHit = fHitItems.at(DET_FI_FIRST + ifibcount);

        const char* detName = fDetectorNames[DET_FI_FIRST + ifibcount];

        if(detMapped) 
        {
            Int_t nHits = detMapped->GetEntriesFast();
            std::vector<UInt_t> mapmt_num(512);
            std::vector<UInt_t> spmt_num(16);
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BBunchedFiberMappedData* hit = (R3BBunchedFiberMappedData*)detMapped->At(ihit);
                if (!hit) continue;

                // channel numbers are stored 1-based (1..n)
                Int_t iCha = hit->GetChannel();  // 1..

                if (hit->IsMAPMT() && hit->IsLeading()) 
                {
                    fh_channels_Fib[ifibcount]->Fill(iCha);    // Fill which clockTDC channel has events
                    ++mapmt_num.at(hit->GetChannel() - 1);     // multihit of a given clockTDC channel
                }

                if (!hit->IsMAPMT() && hit->IsLeading()) 
                {
                    ++spmt_num.at(hit->GetChannel() - 1);      // multihit of a given PADI channel
                }
            }

            for (int i = 0; i < 512; ++i) 
            { 
                auto m = mapmt_num.at(i);
                if(m > 0) fh_multihit_m_Fib[ifibcount]->Fill(i,m);  // multihit of a given clockTDC channel      
            }

            for (int i = 0; i < 16; ++i) 
            {
                auto s = spmt_num.at(i);

                if(s > 0) fh_multihit_s_Fib[ifibcount]->Fill(i,s); // multihit of a given PADI channel
            }
            
        }

        // begin HS_I
        /*
           ofstream CheckTxt;
           CheckTxt.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
           CheckTxt << "Direkt vor der if(detHit) loop! \n";
           CheckTxt.close();
           */


        //bool alreadyAvailable = false;  // HS_I
        bool doMyRoutine = false;
        const int maxAnodeNr = 256;                                                                                                                                             
        std::vector<double> x_Anode_y_Anode_Tot_MA (maxAnodeNr*2*9+1, -1000.);

        // end HS_I

        if(detHit) 
        {

          //printf("Checkmark One \n");
          ofstream CheckFile;
          CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
          CheckFile << "Checkmark One (detHit) \n";
          CheckFile.close();

          Int_t nHits = detHit->GetEntriesFast(); 
          //            cout<<"Fiber Hits: "<<nHits<<endl;

          std::vector<UInt_t> fiber_num(2049);

          Double_t xposfib = 0./0.;
          Double_t yposfib = 0./0.;
          Double_t totMax = 0.;     
          Double_t tfib=0./0., tof_fib_s= 0./0., tof_fib_m= 0./0.;
          Double_t randx;
          Int_t iFibMax = 0;
          Double_t totMax_MA=0.;
          Int_t iFibMax_MA=0;


          Double_t spmtMax;
          Double_t mapmtMax;
          Double_t tofMax;
          Double_t y;
          Double_t tof[14] = {0.,0.,0.,0.,60.,60.,0.,0.,0.,140.,140.,0.,140.,140.};
          Double_t z[14] = {-51.,-49.,0.,0.,45.,50.,0.,0.,0.,700.,670.,0.,630.,600.};
          Double_t x[14] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,-200.,-180.,0.,-150.,-130.};


          if(counter==0 && nHits>0){
            fFibEvents++;
            counter=1;
          }

          // begin HS_I
          
          if (nHits > 0) doMyRoutine = true;

          CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
          CheckFile << "nHits is " << nHits << endl;
          CheckFile.close();

          const int maxAnodesInRow = 16;   // maximum number of anodes in x direction per anode number in y direction (= row); number in both directions is equal


          //if (doMyRoutine == true && alreadyAvailable == false) {
          if (doMyRoutine == true) {
            //printf("detName: %s \n", detName);
            SetUpTheCoordinates(x_Anode_y_Anode_Tot_MA, maxAnodesInRow, maxAnodeNr, detName);
            //alreadyAvailable = true;
          }
          
          TH2F* hist2DAnodeHitMap = new TH2F(Form("2DAnodeHitMap_%d_%s", fNEvents, detName), " ", maxAnodesInRow, 0.5, maxAnodesInRow+0.5, maxAnodesInRow, 0.5, maxAnodesInRow+0.5);


          // end HS_I

          for (Int_t ihit = 0; ihit < nHits; ihit++) {

            //printf("Checkmark Two \n");
            CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
            CheckFile << "Checkmark Two (for ihit) \n";
            CheckFile.close();


            Double_t tMAPMT = 0./0.;
            Double_t tSPMT = 0./0.;

            R3BBunchedFiberHitData* hit = (R3BBunchedFiberHitData*)detHit->At(ihit);
            if (!hit) continue;

            iFib = hit->GetFiberId();  // 1..
            ++fiber_num.at(iFib);     // multihit of a given fiber
            auto m_fib = fiber_num.at(iFib);


            // times
            tMAPMT = hit->GetMAPMTTime_ns();
            tSPMT = hit->GetSPMTTime_ns();

            dtime = tMAPMT-tSPMT;       

            // "Push" two times in the same clock cycle:    
            while(dtime < -1024) 
            {
              tSPMT = tSPMT - 2048.; 
              dtime = tMAPMT - tSPMT;
            }
            while(dtime > 1024.) 
            {
              tSPMT = tSPMT + 2048.; 
              dtime = tMAPMT - tSPMT;
            }           

            // "Push" the Fib times in the same cycle with LOS:
            if(timeTofd[0]>0. && !(IS_NAN(timeTofd[0])))
            {
              while(tSPMT - timeTofd[0] < 4096.)
              {
                tMAPMT = tMAPMT + 2048.*4.; 
                tSPMT = tSPMT + 2048.*4.;         
              }       
              while(tSPMT - timeTofd[0] > 4096.)
              {
                tMAPMT = tMAPMT - 2048.*4.; 
                tSPMT = tSPMT - 2048.*4.;         
              }
            }



            // find fiber with most light
            //      Double_t ToT = sqrt(hit->GetSPMTToT_ns()*hit->GetMAPMTToT_ns());
            Double_t ToT = (3.*hit->GetMAPMTToT_ns()+hit->GetSPMTToT_ns())/4.;
            Double_t ToT_MA = hit->GetMAPMTToT_ns();

            // Not-calibrated ToF:   
            //      tfib = (tMAPMT + tSPMT) / 2.;
            tfib = tSPMT;
            if(tfib > 0. && !(IS_NAN(tfib)) && timeTofd[0]>0. && !(IS_NAN(timeTofd[0]))) tof_fib_s = tfib - timeTofd[0];   
            if(tMAPMT > 0. && !(IS_NAN(tMAPMT)) && timeTofd[0]>0. && !(IS_NAN(timeTofd[0]))) tof_fib_m = tMAPMT - timeTofd[0]; 


            //if(fNEvents<10000 && ifibcount == 13){
            if (fNEvents<10000 && iFib>900 && iFib<1020 && ifibcount==13){
              test=0;
              cout<<"Det: "<<ifibcount
                <<" Fiber: "<<iFib<<" MAToT: "<<hit->GetMAPMTToT_ns()
                <<" SToT: "<<hit->GetSPMTToT_ns()<<" sqrt: "<<ToT
                <<" tof: "<<tof_fib_s << " tM: "<<tof_fib_m << endl;
            }       

            //if(ToT>totMax && abs(tof_fib_s-tof[ifibcount])<20. && ToT<1000.) {
            if(ToT>totMax && ToT<1000.) {
              totMax=ToT;
              iFibMax=iFib;
              spmtMax=hit->GetSPMTToT_ns();
              mapmtMax=hit->GetMAPMTToT_ns();
              tofMax=tof_fib_s;
              // Not-calibrated position:  
              randx = (std::rand() / (float)RAND_MAX);
              xposfib = (-n_fiber[ifibcount]/2.+iFib+(0.5-randx)); 
              yposfib=tMAPMT-tSPMT;           
            }   

            //if(ToT_MA>totMax_MA && abs(tof_fib_s-tof[ifibcount])<20. && ToT_MA<1000.) {
            if(ToT_MA>totMax_MA) {
              totMax_MA=ToT_MA;
              iFibMax_MA=iFib;
            }   

            // My changes in this macro one // sign to search for all related lines in here: begin HS_I

            //PrepareDataForPlotting(x_Anode_y_Anode_Tot_MA, x_FibCoords, y_FibCoords, iFib, ToT_MA, maxAnodeNr);

            x_Anode_y_Anode_Tot_MA.at(iFib+maxAnodeNr*2*2) = ToT_MA;
            x_Anode_y_Anode_Tot_MA.at(iFib+maxAnodeNr*2*3) = tMAPMT;

            // end of this part from my changes one // end HS_I

          }  // end for(ihit)

          // begin HS_I
          if (doMyRoutine == true) {

            //bool stopDueToTimeDiff = false;
            PutHitsIntoHists(x_Anode_y_Anode_Tot_MA,/* fNEvents,*/ maxAnodesInRow, maxAnodeNr, hist2DAnodeHitMap, detName, ifibcount/*, stopDueToTimeDiff*/); 
            //if (stopDueToTimeDiff == true) printf("Das läuft immer noch !!!!!!!!!!!!!!!!!!!!!!!!!!\n");

            for (int FibDummy = 1; FibDummy <= maxAnodeNr * 2; FibDummy++){
              x_Anode_y_Anode_Tot_MA.at(maxAnodeNr*2*2+FibDummy) = -1000.; // reset the Tot_MA values in this vector to -1 before reusing the vector
              //int saveThisEntry = x_Anode_y_Anode_Tot_MA[300];
              //x_Anode_y_Anode_Tot_MA[300]=0;
              if (x_Anode_y_Anode_Tot_MA.at(FibDummy) < 1){
                //printf("problem with X (%f) at FibID %d! \n",x_Anode_y_Anode_Tot_MA[FibDummy], FibDummy);
                //for (int FibDummy2 = 1; FibDummy2 <= maxAnodeNr * 2; FibDummy2++) {
                //printf("FibId: %d, X: %f, Y: %f, ToT_MA: %f \n", FibDummy2, x_Anode_y_Anode_Tot_MA[FibDummy2], x_Anode_y_Anode_Tot_MA[FibDummy2+maxAnodeNr*2], x_Anode_y_Anode_Tot_MA[FibDummy2+maxAnodeNr*2*2]);
                //} // end for(int FibDummy2)
                CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
                CheckFile << "problem with X (" << x_Anode_y_Anode_Tot_MA.at(FibDummy) << ") at FibId " << FibDummy << endl;
                CheckFile.close();
              }
              //x_Anode_y_Anode_Tot_MA[300] = saveThisEntry;
              if (x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2) < 1) {
                //printf("problem with Y (%f) at FibID %d! \n",x_Anode_y_Anode_Tot_MA[FibDummy+maxAnodeNr*2], FibDummy);
                //for (int FibDummy2 = 1; FibDummy2 <= maxAnodeNr * 2; FibDummy2++) {
                //printf("FibId: %d, X: %f, Y: %f, ToT_MA: %f \n", FibDummy2, x_Anode_y_Anode_Tot_MA[FibDummy2], x_Anode_y_Anode_Tot_MA[FibDummy2+maxAnodeNr*2], x_Anode_y_Anode_Tot_MA[FibDummy2+maxAnodeNr*2*2]);
                //} // end for(int FibDummy2)
                CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
                CheckFile << "problem with Y (" << x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2) << ") at FibId " << FibDummy << endl;
                CheckFile.close();
              }
            } // end for (int FibDummy)

            /*
               for (int FibDummy = 1; FibDummy <= maxAnodeNr * 2; FibDummy++) {
               printf("FibId: %d, X: %f, Y: %f, ToT_MA: %f \n", FibDummy, x_Anode_y_Anode_Tot_MA[FibDummy], x_Anode_y_Anode_Tot_MA[FibDummy+maxAnodeNr*2], x_Anode_y_Anode_Tot_MA[FibDummy+maxAnodeNr*2*2]);
               } // end for(int FibDummy)
            */

          }

          doMyRoutine = false;
          if (hist2DAnodeHitMap->GetEntries() == 1) printf("Event with one Entry spotted! (%d) \n", fNEvents);
          x_Anode_y_Anode_Tot_MA.clear();
          delete hist2DAnodeHitMap;  

          // loop to cut out the following part of this function to save memory; make comment if everything beyond here is needed again
          
          /*
          if (ifibcount+1 >= NOF_FIB_DET) {                                          
            CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
            CheckFile << "Event bis zum Ende meiner MOD abgearbeitet! (Event #" << fNEvents << ")" << endl;
            CheckFile.close();
            printf("Return by me!\n");
            fNEvents++;
            return;
          }
          else if (ifibcount < NOF_FIB_DET-1){
            CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
            CheckFile << "Continue to next Fib from" << ifibcount << " (Fiber detector " << detName << ")" << endl;
            CheckFile.close(); 
            printf("Continue by me!\n");
            continue;
          }
          */

          // end HS_I


          //      if(abs(tof_fib_s-tof[ifibcount])<20.) 
          
          if(totMax>0) 
          {
            if(TofdQ>0. && TofdQ<100.) {
              fh_fibers_Fib[ifibcount]->Fill(iFibMax);  
              fh_fiber_Fib[ifibcount]->Fill(iFibMax);  
              fh_ToT_s_Fib[ifibcount]->Fill(iFibMax,spmtMax);
              fh_ToT_m_Fib[ifibcount]->Fill(iFibMax,mapmtMax);
              fh_time_Fib[ifibcount]->Fill(iFibMax,yposfib);
              fh_Fib_ToF[ifibcount]->Fill(iFibMax,tofMax);
              fh_xpos_Fib[ifibcount]->Fill(xposfib);  
              fh_ypos_Fib[ifibcount]->Fill(yposfib);  
              fh_Fibs_vs_Events[ifibcount]->Fill(fFibEvents,iFibMax); 
              fh_Fib_vs_Events[ifibcount]->Fill(fFibEvents,iFibMax_MA);   
              FibMax[ifibcount]=iFibMax;
            }
          }

          if (nHits>0) fh_mult_Fib[ifibcount]->Fill(nHits);      

          if(TofdQ>0. && TofdQ<100.) {
            fh_Fibs_vs_Tofd[ifibcount]->Fill(TofdX,iFibMax);
            //cout<<"test "<<TofdX<<"  "<<iFibMax<<endl;
          }
          fh_Cave_position->Fill(z[ifibcount],x[ifibcount]+xposfib);


          if(test) cout<<"fiber max: "<<iFibMax<<endl;



          }  // end if(aHit[ifibcount]) 

          if (ifibcount==13 && test){
            cout<<"   "<<endl;
            cout<<"new event! ******************************************************"<<endl;
          }

          //printf("ifibcount is %d, NOF_FIB_DET is %d \n", ifibcount, NOF_FIB_DET);

          } // end for(ifibcount)

          for(Int_t i = 0; i < NOF_FIB_DET; i++){
            for(Int_t j = i+1; j < NOF_FIB_DET; j++){
              if(fHitItems.at(i + DET_FI_FIRST) && fHitItems.at(j + DET_FI_FIRST))  {
                if(TofdQ>0. && TofdQ<100.) {
                  fh_Fib_vs_Fib[i][j]->Fill(FibMax[i],FibMax[j]);
                  if(FibMax[i]>0 && FibMax[j]>0) fh_Fib_dx[i][j]->Fill(FibMax[i],FibMax[j]-FibMax[i]);				
                }
              }
            }
          }

          // begin HS_I
          ofstream CheckFile;
          CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios:: app);
          CheckFile << "Ein Event abgearbeitet! (" << fNEvents << ")" << endl;
          CheckFile.close();
          // end HS_I

          fNEvents += 1;

}

void R3BGlobalAnalysis::FinishEvent()
{
  for(Int_t det = 0; det < DET_MAX; det++) {
    if(fMappedItems.at(det)) {
      fMappedItems.at(det)->Clear(); 
    }
    if(fCalItems.at(det)) {
      fCalItems.at(det)->Clear(); 
    }
    if(fHitItems.at(det)) {
      fHitItems.at(det)->Clear(); 
    }
  }       
}

void R3BGlobalAnalysis::FinishTask()
{    
  fh_Cave_position->Write();
  if(fMappedItems.at(DET_LOS)){
    fhChargeLosTofD->Write();
    fh_los_pos->Write();
  }
  if(fHitItems.at(DET_TOFD))
  {
    fh_tofd_pos->Write();
    fh_tofd_charge->Write();

  }    
  for(Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++) {  
    if(fMappedItems.at(ifibcount + DET_FI_FIRST)) 
    {
      fh_channels_Fib[ifibcount]->Write();
      fh_fibers_Fib[ifibcount]->Write();
      fh_fiber_Fib[ifibcount]->Write();
      fh_mult_Fib[ifibcount]->Write();
      fh_time_Fib[ifibcount]->Write();
      fh_xpos_Fib[ifibcount]->Write();
      fh_ypos_Fib[ifibcount]->Write();
      fh_multihit_m_Fib[ifibcount]->Write();   
      fh_multihit_s_Fib[ifibcount]->Write();
      fh_ToT_m_Fib[ifibcount]->Write();
      fh_ToT_s_Fib[ifibcount]->Write();
      fh_Fib_vs_Events[ifibcount]->Write();
      fh_Fibs_vs_Events[ifibcount]->Write();
      fh_Fibs_vs_Tofd[ifibcount]->Write();
      fh_Fib_ToF[ifibcount]->Write();
    }
  }   
    for(Int_t i = 0; i < NOF_FIB_DET; i++){
      for(Int_t j = i+1; j < NOF_FIB_DET; j++){
        if(fHitItems.at(i + DET_FI_FIRST) && fHitItems.at(j + DET_FI_FIRST)){ 
          fh_Fib_vs_Fib[i][j]->Write();
          fh_Fib_dx[i][j]->Write();
        }
      }
    }

}


// begin HS_I

void SetUpTheCoordinates(std::vector<double>& x_Anode_y_Anode_Tot_MA, int maxAnodesInRow, int maxAnodeNr, const char* detName)
{

  //printf("The function SetUpTheCoordinates is being executed! \n");
  ofstream CheckFile;
  CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
  CheckFile << "The function SetUpTheCoordinates is being executed! \n";
  CheckFile.close();

  const int maxFibsPerAnode = 2;
  int presentAnodeYCoord;
  int firstAnodeCoord;
  int nextAnodeNrInScheme;
  int yUsedIndex = 0;

  //std::cout << "detName is" <<  detName << endl;
   
  if (detName == "Fi3a" || detName == "Fi3b") {
    //printf("detName wurde erfolgreich als Fi3[] identifiziert. \n");
    firstAnodeCoord = 16; 
    nextAnodeNrInScheme = -4;
  }
  else {
    firstAnodeCoord = 1;
    nextAnodeNrInScheme = 4;
  }

  //printf("firstAnodeCoord is %d and nextAnodeNrInScheme is %d \n", firstAnodeCoord, nextAnodeNrInScheme);

  for (int dummyRow = 1; dummyRow <= maxAnodesInRow; dummyRow++) {
    presentAnodeYCoord = firstAnodeCoord;

    for (int dummyFiber = 1; dummyFiber <= maxAnodesInRow * 2; dummyFiber++) {

      x_Anode_y_Anode_Tot_MA.at((dummyRow - 1) * maxAnodesInRow * 2 + dummyFiber)                 = dummyRow;
      x_Anode_y_Anode_Tot_MA.at(maxAnodeNr * 2 +(dummyRow - 1) * maxAnodesInRow * 2 + dummyFiber) = presentAnodeYCoord;
      yUsedIndex++;

      if (yUsedIndex == maxFibsPerAnode) {
        presentAnodeYCoord += nextAnodeNrInScheme;
        yUsedIndex = 0;
      }

      if (presentAnodeYCoord > maxAnodesInRow || presentAnodeYCoord < 1) {

        if (detName == "Fi3a" || detName == "Fi3b") presentAnodeYCoord += (16-1);
        else presentAnodeYCoord += (-16+1); 
      }
                                                                                                                                                       
    } // end for (int dummyFiber)
  } // end for (int dummyRow)

  //for (int FibDummy = 1; FibDummy <= maxAnodeNr * 2; FibDummy++) {
    //printf("FibId: %d, X: %f, Y: %f, ToT_MA: %f \n", FibDummy, x_Anode_y_Anode_Tot_MA[FibDummy], x_Anode_y_Anode_Tot_MA[FibDummy+maxAnodeNr*2], x_Anode_y_Anode_Tot_MA[FibDummy+maxAnodeNr*2*2]);
  //} // end for(int FibDummy)
}

/*
void PrepareDataForPlotting(std::vector<double>& x_Anode_y_Anode_Tot_MA, std::vector<int>& x_FibCoords, std::vector<int>& y_FibCoords, int iFib, double ToT_MA, int maxAnodeNr) {

  //printf("The function PrepareDataForPlotting is being executed! \n");
  ofstream CheckFile;
  CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
  CheckFile << "The function PrepareDataForPlotting is being executed! \n"; 
  CheckFile.close();

  x_Anode_y_Anode_Tot_MA[iFib]                = x_FibCoords[iFib];
  x_Anode_y_Anode_Tot_MA[iFib+maxAnodeNr*2]   = y_FibCoords[iFib];
  x_Anode_y_Anode_Tot_MA[iFib+maxAnodeNr*2*2] = ToT_MA;

}
*/

void R3BGlobalAnalysis::PutHitsIntoHists(std::vector<double>& x_Anode_y_Anode_Tot_MA, /*const int fNEvents,*/ int maxAnodesInRow, int maxAnodeNr, TH2F* hist2DAnodeHitMap, const char* detName, const int ifibcount/*, bool stopDueToTimeDiff*/) {
//void PutHitsIntoHists(std::vector<double>& x_Anode_y_Anode_Tot_MA, const int fNEvents, int maxAnodesInRow, int maxAnodeNr, TH2F* hist2DAnodeHitMap, const char* detName, int& CrossTalkLabel) {

  //printf("The function PutHitsIntoHists is being executed! \n");
  ofstream CheckFile;
  CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
  CheckFile << "The function PutHitsIntoHists is being executed! \n"; 
  CheckFile.close();

  bool saveThisEvent = false;
  int binNumberOfMaxInX = -50;
  int binNumberOfMaxInY = -50;
  int binNumberOfMaxInZ = -50;
  int CrossTalkLabel = -1;     // 0: Single Peak | 1: cross talk on Anode (can also be between fibers near the Anode plane) | 2: fibrous cross talk | 3: 1 and 2 | 4: Fake case of Label 2

  TH1F* hist1DMAPMTTimeDistEvent = new TH1F(Form("MAPMTTimeDistEvent_%d_%s", fNEvents, detName), " ", 330000, 0, 16500);
  hist1DMAPMTTimeDistEvent->GetXaxis()->SetTitle("Time Stamp (tMAPMT) [ns]");
  hist1DMAPMTTimeDistEvent->GetYaxis()->SetTitle("Entries");

  //TH1F* hist1DMAPMTTimeDiffsEvent = new TH1F(Form("MAPMTTimeDiffsEvent_%d_%s", fNEvents, detName), " ", 40000, 0, 2000);
  TH1F* hist1DMAPMTTimeDiffsEvent = new TH1F(Form("MAPMTTimeDiffsEvent_%d_%s", fNEvents, detName), " ", 1000, 0, 100);
  //hist1DMAPMTTimeDiffsEvent->SetName(Form("MAPMTTimeDiffsEvent_%d_%s", fNEvents, detName));
  hist1DMAPMTTimeDiffsEvent->GetXaxis()->SetTitle("#Delta t (tMAPMT stamps to first tMAPMT of an event) [ns]");
  hist1DMAPMTTimeDiffsEvent->GetYaxis()->SetTitle("Entries");

  bool bouncer = false;
  double firstTimeValueOfEvent = 999999999999999.;
  std::vector<double> timeStampsOfEvent;

  for (int FibDummy=1; FibDummy <= maxAnodeNr*2; FibDummy++) {
    //CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    //CheckFile << "FibID " << FibDummy << ", x pos " << x_Anode_y_Anode_Tot_MA.at(FibDummy) << ", y pos " << x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2) << ", tot_Ma value " << x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*2) << ", MAPMT time stamp" << x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*3) << endl;
    //CheckFile.close();


    if(x_Anode_y_Anode_Tot_MA.at(FibDummy) == x_Anode_y_Anode_Tot_MA.at(FibDummy-1) && x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2) == x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2-1) /*&& x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*2) == x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*2-1)*/ && x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*3) == x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*3-1)) bouncer = true;
    else bouncer = false;

    if(bouncer == false) {
      if (x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*2) != -1000.) {
        hist2DAnodeHitMap->Fill(x_Anode_y_Anode_Tot_MA.at(FibDummy), x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2), x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*2));

        CheckFile.open("/u/schulte/testStation/CheckTimeDiffs_318.txt", ios::out | ios::app);   // FibTests082019_Run318_16oxygenNorm; to not overlook to change the run number here! 
        CheckFile << "Detector " << detName << " EventNr " << fNEvents << " FibID " << FibDummy << endl << "x = " << x_Anode_y_Anode_Tot_MA.at(FibDummy) << " y = " << x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2) << " ToT_MA = " << x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*2) << " tMAPMT = " << x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*3) << endl << "----------------------------------------------------------------------------------------------------------" << endl;
        CheckFile.close();

      }

      if (x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*3) != -1000.) {
        double TimeValueOfHit = x_Anode_y_Anode_Tot_MA.at(FibDummy+maxAnodeNr*2*3); 
        hist1DMAPMTTimeDistEvent         ->Fill(TimeValueOfHit);
        hist1DMAPMTTimeDistAll[ifibcount]->Fill(TimeValueOfHit);

        if (firstTimeValueOfEvent == 999999999999999.) firstTimeValueOfEvent = TimeValueOfHit;
        else if (TimeValueOfHit < firstTimeValueOfEvent) {
          timeStampsOfEvent.push_back(firstTimeValueOfEvent);
          firstTimeValueOfEvent = TimeValueOfHit;
        }
        else timeStampsOfEvent.push_back(TimeValueOfHit);
      }
    }
  } // end for (int FibDummy)

  for (int step = 0; step < timeStampsOfEvent.size(); step++) {
    const double deltaTime = timeStampsOfEvent.at(step) - firstTimeValueOfEvent;
    //if (deltaTime > 5) stopDueToTimeDiff = true;
    hist1DMAPMTTimeDiffsEvent    ->Fill(deltaTime);
    hist1DTimeDiffsAll[ifibcount]->Fill(deltaTime);
  } // end for (int step)
  timeStampsOfEvent.clear();
  //if (stopDueToTimeDiff == true){
    //printf("TimeDiff to big to be real CT!!!!!!!!!!!!!!!!!!!!\n");
    //fNEvents++;
    //return;
  //}
  //printf("Das kommt noch!!!!!!!!!!!!!!!!!!!!!!\n");

  hist2DAnodeHitMap->GetMaximumBin(binNumberOfMaxInX, binNumberOfMaxInY, binNumberOfMaxInZ);
  
  if (hist2DAnodeHitMap->GetEntries() >= 1 && 3 <= binNumberOfMaxInX && binNumberOfMaxInX  <= 13 && 3 <= binNumberOfMaxInY && binNumberOfMaxInY  <= 13) saveThisEvent = true;
  //saveThisEvent = true;

  //if (saveThisEvent == true) {

    CheckTheEventsCrossTalk(hist2DAnodeHitMap, CrossTalkLabel, ifibcount);

    if (/*CrossTalkLabel == 0 && */ saveThisEvent == true) {
    /*
      CheckFile.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
      CheckFile << "single peak in MAPMmap detected; return send" << endl;
      CheckFile.close();
      printf("Single Peak Event! Not Saved \n");
      return;
    */
      hist2DInvestigateMaxEntriesI[ifibcount] ->Fill(CrossTalkLabel, hist2DAnodeHitMap->GetBinContent(binNumberOfMaxInX,binNumberOfMaxInY)); 
      //hist3DInvestigateMaxEntriesII[ifibcount]->Fill(CrossTalkLabel, hist2DAnodeHitMap->GetBinContent(binNumberOfMaxInX,binNumberOfMaxInY), hist2DAnodeHitMap->GetEntries()); 
      if (CrossTalkLabel == 0) saveThisEvent = false;
    }


    if (saveThisEvent == true /*|| fNEvents == 10016589*/) {
    hist2DAnodeHitMap->GetXaxis()->SetTitle("x-Direction"); 
    hist2DAnodeHitMap->GetXaxis()->SetTitleOffset(1.4); 
    hist2DAnodeHitMap->GetYaxis()->SetTitle("y-Direction");
    hist2DAnodeHitMap->GetZaxis()->SetTitle("Time over threshold multi anode");
    if (CrossTalkLabel > 0) {
      hist2DAnodeHitMap->SetTitle(Form("Cross talk label %d", CrossTalkLabel));
      hist2DAnodeHitMap->SetName(Form("2DAnodeHitMap_%d_%s_CTLabel%d", fNEvents, detName, CrossTalkLabel));
      hist1DMAPMTTimeDiffsEvent->SetName(Form("MAPMTTimeDiffsEvent_%d_%s_%d", fNEvents, detName, CrossTalkLabel));
      hist1DMAPMTTimeDistEvent ->SetName(Form("MAPMTTimeDistEvent_%d_%s_%d", fNEvents, detName, CrossTalkLabel));
    }

    //for (int FibDummy=1; FibDummy <= maxAnodeNr*2; FibDummy++) {
    //} // end for (int FibDummy)

    TFile* outputII = TFile::Open(Form("/u/schulte/testStation/BeamData042019_Run318/OverviewHists%s.root", detName), "UPDATE");
    //TFile* outputII = TFile::Open(Form("/u/schulte/testStation/FibTests082019_Run15_16oxygenNorm/OverviewHists%s.root", detName), "UPDATE");
    outputII->cd();
    hist1DMAPMTTimeDistAll[ifibcount]->Write(hist1DMAPMTTimeDistAll[ifibcount]->GetName(), 2);
    hist1DTimeDiffsAll[ifibcount]->Write(hist1DTimeDiffsAll[ifibcount]->GetName(), 2);
    hist2DInvestigateMaxEntriesI[ifibcount]->Write(hist2DInvestigateMaxEntriesI[ifibcount]->GetName(), 2);
    if (hist1DMAPMTTimeDiffsEvent->GetBinContent(1) > 0) {
      hist1DMAPMTTimeDiffsEvent->Write();
      hist1DMAPMTTimeDistEvent->Write();
    }
    outputII->Close();
    //hist3DInvestigateMaxEntriesII[ifibcount]->Write(hist3DInvestigateMaxEntriesII[ifibcount]->GetName(), 2);
    delete outputII;
    if (CrossTalkLabel == 2 || CrossTalkLabel == 4) {
      TFile* output = TFile::Open(Form("/u/schulte/testStation/BeamData042019_Run318/FibHitInAnodeCoordsMaps_FiberDetector%s.root", detName), "UPDATE");
      //TFile* output = TFile::Open(Form("/u/schulte/testStation/FibTests082019_Run15_16oxygenNorm/FibHitInAnodeCoordsMaps_FiberDetector%s.root", detName), "UPDATE");
      output->cd();
      hist2DAnodeHitMap->Write();
      //hist1DMAPMTTimeDistEvent->Write();
      //hist1DMAPMTTimeDiffsEvent->GetXaxis()->SetRangeUser(hist1DMAPMTTimeDiffsEvent->FindFirstBinAbove(0,1)-5, hist1DMAPMTTimeDiffsEvent->FindLastBinAbove(0,1)+5);
      //hist1DMAPMTTimeDiffsEvent->Write();
      output->Close();
      delete output; 
    }

    if (CrossTalkLabel == 1) GetHistsCrossTalkLabel1(hist2DAnodeHitMap, ifibcount, detName);
    else if (CrossTalkLabel == 2) GetHistsCrossTalkLabel2(hist2DAnodeHitMap, ifibcount, detName, hist1DMAPMTTimeDistEvent, hist1DMAPMTTimeDiffsEvent);
    else if (CrossTalkLabel == 4) GetHistsCrossTalkLabel4(hist2DAnodeHitMap, ifibcount, detName, hist1DMAPMTTimeDistEvent, hist1DMAPMTTimeDiffsEvent);

    printf("Histograms Saved for Event #%d \n", fNEvents);
  }
  else {
    printf("xMax: %d, yMax: %d \n", binNumberOfMaxInX, binNumberOfMaxInY);
  }

  //if (CrossTalkLabel == 1) GetHistsCrossTalkLabel1(hist2DAnodeHitMap, ifibcount, detName);
  //else if (CrossTalkLabel == 2) GetHistsCrossTalkLabel2(hist2DAnodeHitMap, ifibcount, detName, hist1DMAPMTTimeDistEvent, hist1DMAPMTTimeDiffsEvent);
  //else if (CrossTalkLabel == 4) GetHistsCrossTalkLabel4(hist2DAnodeHitMap, ifibcount, detName, hist1DMAPMTTimeDistEvent, hist1DMAPMTTimeDiffsEvent);

  delete hist1DMAPMTTimeDistEvent;
  delete hist1DMAPMTTimeDiffsEvent;
}

void R3BGlobalAnalysis::CheckTheEventsCrossTalk(TH2F* hist2DAnodeHitMap, int& CrossTalkLabel, const int ifibcount){

  const int DiffFirstLastEntrX = hist2DAnodeHitMap->FindLastBinAbove(0,1) - hist2DAnodeHitMap->FindFirstBinAbove(0,1);
  const int DiffFirstLastEntrY = hist2DAnodeHitMap->FindLastBinAbove(0,2) - hist2DAnodeHitMap->FindFirstBinAbove(0,2);
  const int FirstBinWContX = hist2DAnodeHitMap->FindFirstBinAbove(0,1); 
  const int FirstBinWContY = hist2DAnodeHitMap->FindFirstBinAbove(0,2); 
  const int LastBinWContX  = hist2DAnodeHitMap->FindLastBinAbove(0,1); 
  const int LastBinWContY  = hist2DAnodeHitMap->FindLastBinAbove(0,2);
  int BetweenFirstLastBinsXY = -5;
  int OuterHitBins = 0;
  int xOfMaxBin = -1;
  int yOfMaxBin = -1;
  int zOfMaxBin = -1;
  hist2DAnodeHitMap->GetMaximumBin(xOfMaxBin, yOfMaxBin, zOfMaxBin);

  bool possibleLabel2Case = false;


  ofstream CheckLabel;

  if (FirstBinWContX+1 <= LastBinWContX-1 && FirstBinWContY+1 <= LastBinWContY-1) {

    BetweenFirstLastBinsXY = hist2DAnodeHitMap->Integral(FirstBinWContX+1, LastBinWContX-1, FirstBinWContY+1,LastBinWContY-1);
    possibleLabel2Case = true;
    
    CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckLabel << "possible Label 2 case detected!" << endl;
    CheckLabel.close();
  }

  if (FirstBinWContX+1 > LastBinWContX-1 && FirstBinWContY+1 <= LastBinWContY-1) {


    BetweenFirstLastBinsXY = hist2DAnodeHitMap->Integral(FirstBinWContX, LastBinWContX, FirstBinWContY+1, LastBinWContY-1);
    possibleLabel2Case = true;

    CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckLabel << "possible Label 2 case detected!" << endl;
    CheckLabel.close();
  }

  if (FirstBinWContX+1 <= LastBinWContX-1 && FirstBinWContY+1 > LastBinWContY-1) {


    BetweenFirstLastBinsXY = hist2DAnodeHitMap->Integral(FirstBinWContX+1, LastBinWContX-1, FirstBinWContY, LastBinWContY);
    possibleLabel2Case = true;

    CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckLabel << "possible Label 2 case detected!" << endl;
    CheckLabel.close();
  }

  if (possibleLabel2Case ==true) {
   for (int binX = FirstBinWContX; binX <= LastBinWContX; binX++) {
   if (hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(binX, FirstBinWContY)) > 0) OuterHitBins += 1;
   if (hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(binX, LastBinWContY)) > 0)  OuterHitBins += 1;
   } 
   if (DiffFirstLastEntrY >= 2) {
     for (int binY = FirstBinWContY+1; binY <= LastBinWContY-1; binY++) {
       if (hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX, binY)) > 0) OuterHitBins += 1;
       if (hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(LastBinWContX, binY)) > 0)  OuterHitBins += 1;
     } 
   }
   if (FirstBinWContY == LastBinWContY) OuterHitBins += -2;  //two Hits will be both counted twice by the logic above!
  }

  if (DiffFirstLastEntrX <= 2 && DiffFirstLastEntrY <= 2 && xOfMaxBin-1 <= FirstBinWContX && LastBinWContX <= xOfMaxBin+1 && yOfMaxBin-1 <= FirstBinWContY && LastBinWContY <= yOfMaxBin+1) {
    CrossTalkLabel = 1; 
    CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckLabel << "Label 1 given! (" << BetweenFirstLastBinsXY << ")" << endl;
    CheckLabel.close();
  }

  if (DiffFirstLastEntrX == 0 && DiffFirstLastEntrY == 0) {
    CrossTalkLabel = 0;
    CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckLabel << "Label 0 given! (" << BetweenFirstLastBinsXY << ")" << endl;
    CheckLabel.close();
  }

  // Label 3 conditions have to be worked over again /worked out newly if Label 3 needed!
  //if ((DiffFirstLastEntrX >= 3 && ((hist2DAnodeHitMap->Integral(FirstBinWContX+2, LastBinWContX-1) - hist2DAnodeHitMap->Integral(FirstBinWContX+2, FirstBinWContX+2, FirstBinWContY, LastBinWContY)) == 0 || (hist2DAnodeHitMap->Integral(FirstBinWContX+1, LastBinWContX-2) - hist2DAnodeHitMap->Integral(LastBinWContX-2, LastBinWContX-2, FirstBinWContY, LastBinWContY))  == 0)) || (DiffFirstLastEntrY >= 3 && ((hist2DAnodeHitMap->Integral(FirstBinWContY+2, LastBinWContY-1) - hist2DAnodeHitMap->Integral(FirstBinWContX, LastBinWContX, FirstBinWContY+2, FirstBinWContY+2)) == 0 || (hist2DAnodeHitMap->Integral(FirstBinWContY+1, LastBinWContY-2) - hist2DAnodeHitMap->Integral(FirstBinWContX, LastBinWContX, LastBinWContY-2, LastBinWContY-2)) == 0))) {

    //double XOne;double YOne; double XTwo; double YTwo; double XThree; double YThree; double XFour; double YFour;

    //if (DiffFirstLastEntrX == 3 && DiffFirstLastEntrY == 0){
      ////XOne   = hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX, FirstBinWContY)); 
      //XTwo   = hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX+1, FirstBinWContY)); 
      //XThree = hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX+2, FirstBinWContY)); 
      ////XFour  = hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX+3, FirstBinWContY)); 
    //}

    //if (DiffFirstLastEntrY == 3 && DiffFirstLastEntrX == 0){
      ////YOne   = hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX, FirstBinWContY)); 
      //YTwo   = hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX, FirstBinWContY+1)); 
      //YThree = hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX, FirstBinWContY+2)); 
      ////YFour  = hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(FirstBinWContX, FirstBinWContY+3)); 
    //}

    //if ((DiffFirstLastEntrX == 3 && DiffFirstLastEntrY == 0 && (XTwo != 0 && XThree != 0) || (XTwo == 0 && XThree == 0)) || (DiffFirstLastEntrY == 3 && DiffFirstLastEntrX == 0 && (YTwo != 0 && YThree != 0) || (YTwo == 0 && YThree == 0))){

    //}
    //else {
      //CrossTalkLabel = 3;
      //CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
      //CheckLabel << "Label 3 given! (" << BetweenFirstLastBinsXY << ")" << endl;
      //CheckLabel.close();
    //}
  //}
  
  if (BetweenFirstLastBinsXY == 0 && OuterHitBins == 2) {
    if ((DiffFirstLastEntrX == 0 && (DiffFirstLastEntrY == 4 || DiffFirstLastEntrY == 11) || DiffFirstLastEntrX == 1 && DiffFirstLastEntrY == 15)) {
    CrossTalkLabel = 2;
    CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckLabel << "Label 2 given! (" << BetweenFirstLastBinsXY << ")" << endl;
    CheckLabel.close();
    }
    else {
    CrossTalkLabel = 4;
    CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
    CheckLabel << "Label 4 given! (" << BetweenFirstLastBinsXY << ")" << endl;
    CheckLabel.close();
    }
  }

  CheckLabel.open("/u/schulte/testStation/Checks.txt", ios::out | ios::app);
  CheckLabel << "OuterHitBins is " << OuterHitBins << endl;
  CheckLabel.close();

  if (CrossTalkLabel == 1) hist2DMaxAnode_CTL1[ifibcount]->Fill(xOfMaxBin, yOfMaxBin);
  else if (CrossTalkLabel == 2) hist2DMaxAnode_CTL2[ifibcount]->Fill(xOfMaxBin, yOfMaxBin);
  else if (CrossTalkLabel == 4) hist2DMaxAnode_CTL4[ifibcount]->Fill(xOfMaxBin, yOfMaxBin);
    
}

void R3BGlobalAnalysis::GetHistsCrossTalkLabel1(TH2F* hist2DAnodeHitMap, const int ifibcount, const char* detName) {

  int neighbours = 0;
  double neighboursToTSum = 0;
  int xOfMaxBin = -1;
  int yOfMaxBin = -1;
  int zOfMaxBin = -1;
  hist2DAnodeHitMap->GetMaximumBin(xOfMaxBin, yOfMaxBin, zOfMaxBin);

  for(int x = xOfMaxBin-1; x <= xOfMaxBin+1; x++) {
    for(int y = yOfMaxBin-1; y <= yOfMaxBin+1; y++) {
      if(hist2DAnodeHitMap->GetBinContent(x,y) > 0) {
        neighboursToTSum += hist2DAnodeHitMap->GetBinContent(x,y);
        neighbours++;
      }
    }
  }
  neighboursToTSum += -1. * hist2DAnodeHitMap->GetBinContent(xOfMaxBin, yOfMaxBin);
  neighbours += -1;

  //printf("Number of neighbours is %d, their ToT_MA_Sum is %f \n", neighbours, neighboursToTSum);

  hist1DAnodeNeighbours[ifibcount]->Fill(neighbours);
  //hist1DAnodeNeighboursVsEventNr[ifibcount]->Fill(fNEvents, neighbours);
  hist1DAnodeNeighboursToTSum[ifibcount]->Fill(neighboursToTSum);
  hist2DAnodeNeighboursToTSumVsNeighbourNr[ifibcount]->Fill(neighbours, neighboursToTSum * 1./neighbours);

  //TFile* output = TFile::Open("/u/schulte/testStation/FibTests082019_Run15_16oxygenNorm/Label1CrossTalkHists.root", "UPDATE");
  TFile* output = TFile::Open("/u/schulte/testStation/BeamData042019_Run318/Label1CrossTalkHists.root", "UPDATE");
  output->cd();
  hist1DAnodeNeighbours[ifibcount]                   ->Write(hist1DAnodeNeighbours[ifibcount]->GetName(), 2);
  //hist1DAnodeNeighboursVsEventNr[ifibcount]          ->Write();
  hist1DAnodeNeighboursToTSum[ifibcount]             ->Write(hist1DAnodeNeighboursToTSum[ifibcount]->GetName(), 2);
  hist2DAnodeNeighboursToTSumVsNeighbourNr[ifibcount]->Write(hist2DAnodeNeighboursToTSumVsNeighbourNr[ifibcount]->GetName(), 2); 
  hist2DMaxAnode_CTL1[ifibcount]                     ->Write(hist2DMaxAnode_CTL1[ifibcount]->GetName(), 2);
  if (neighbours < 1 || 3 < neighbours) hist2DAnodeHitMap->Write();
  output->Close();
  delete output;

  //GetAverageToTMAPerNeighbour(fNEvents, neighbours, neighboursToTSum, detName);
}

void R3BGlobalAnalysis::GetHistsCrossTalkLabel2(TH2F* hist2DAnodeHitMap, const int ifibcount, const char* detName, TH1F* hist1DMAPMTTimeDistEvent, TH1F* hist1DMAPMTTimeDiffsEvent) {

  for (int x=0; x<hist2DAnodeHitMap->GetNbinsX(); x++) {
    for (int y=0; y<hist2DAnodeHitMap->GetNbinsY(); y++) {
      double binContent=hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(x,y));
      if (binContent != 0. && binContent < hist2DAnodeHitMap->GetMaximum()) hist2DSmallEntrVSBigEntr[ifibcount]->Fill(hist2DAnodeHitMap->GetMaximum(), binContent);
    }
  }
  //hist1DMAPMTTimeDist_CTL2[ifibcount]->Add(hist1DMAPMTTimeDistEvent);
  hist1DTimeDiffs_CTL2[ifibcount]->Add(hist1DMAPMTTimeDiffsEvent) ;
  TFile* output = TFile::Open("/u/schulte/testStation/BeamData042019_Run318/Label2CrossTalkHists.root", "UPDATE");
  //TFile* output = TFile::Open("/u/schulte/testStation/FibTests082019_Run15_16oxygenNorm/Label2CrossTalkHists.root", "UPDATE");
  output->cd();
  hist2DSmallEntrVSBigEntr[ifibcount]->Write(hist2DSmallEntrVSBigEntr[ifibcount]->GetName(), 2);
  hist2DMaxAnode_CTL2[ifibcount]->Write(hist2DMaxAnode_CTL2[ifibcount]->GetName(), 2);
  hist1DTimeDiffs_CTL2[ifibcount]->Write(hist1DTimeDiffs_CTL2[ifibcount]->GetName(), 2);
  output->Close();
  delete output;
}

void R3BGlobalAnalysis::GetHistsCrossTalkLabel4(TH2F* hist2DAnodeHitMap, const int ifibcount, const char* detName, TH1F* hist1DMAPMTTimeDistEvent, TH1F* hist1DMAPMTTimeDiffsEvent) {

  for (int x=0; x<hist2DAnodeHitMap->GetNbinsX(); x++) {
    for (int y=0; y<hist2DAnodeHitMap->GetNbinsY(); y++) {
      double binContent=hist2DAnodeHitMap->GetBinContent(hist2DAnodeHitMap->GetBin(x,y));
      if (binContent != 0. && binContent < hist2DAnodeHitMap->GetMaximum()) hist2DSmallEntrVSBigEntr_CTL4[ifibcount]->Fill(hist2DAnodeHitMap->GetMaximum(), binContent);
    }
  }
  //hist1DMAPMTTimeDist_CTL4[ifibcount]->Add(hist1DMAPMTTimeDistEvent);
  hist1DTimeDiffs_CTL4[ifibcount]->Add(hist1DMAPMTTimeDiffsEvent) ;
  TFile* output = TFile::Open("/u/schulte/testStation/BeamData042019_Run318/Label4CrossTalkHists.root", "UPDATE");
  //TFile* output = TFile::Open("/u/schulte/testStation/FibTests082019_Run15_16oxygenNorm/Label4CrossTalkHists.root", "UPDATE");
  output->cd();
  hist2DSmallEntrVSBigEntr_CTL4[ifibcount]->Write(hist2DSmallEntrVSBigEntr_CTL4[ifibcount]->GetName(), 2);
  hist2DMaxAnode_CTL4[ifibcount]->Write(hist2DMaxAnode_CTL4[ifibcount]->GetName(), 2);
  hist1DTimeDiffs_CTL4[ifibcount]->Write(hist1DTimeDiffs_CTL4[ifibcount]->GetName(), 2);
  output->Close();
  delete output;
}

void GetAverageToTMAPerNeighbour(int fNEvents, int neighbours, double neighboursToTSum, const char* detName) {
  TH1F* hist1DAnodeNeighboursToTSumVsNeighbourNr = new TH1F(Form("hist1DAnodeNeighboursToTSumVsNeighbourNr_Event%d_%s", fNEvents, detName), Form("On %s", detName), 9, 0, 9);
  hist1DAnodeNeighboursToTSumVsNeighbourNr->GetXaxis()->SetTitle("Number of neighbours around max");
  hist1DAnodeNeighboursToTSumVsNeighbourNr->GetYaxis()->SetTitle("Sum of ToT values of neighbours per number of neighbours");
  
  hist1DAnodeNeighboursToTSumVsNeighbourNr->Fill(neighbours, neighboursToTSum * 1./neighbours);
  TFile* output = TFile::Open("/u/schulte/testStation/BeamData042019_Run318/Label1CrossTalkHists.root", "UPDATE");
  //TFile* output = TFile::Open("/u/schulte/testStation/FibTests082019_Run15_16oxygenNorm/Label1CrossTalkHists.root", "UPDATE");
  output->cd();
  hist1DAnodeNeighboursToTSumVsNeighbourNr->Write();
  output->Close();
  delete output;


}

//end HS_I

ClassImp(R3BGlobalAnalysis)
