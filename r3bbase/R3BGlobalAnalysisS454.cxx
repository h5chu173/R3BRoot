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

// ------------------------------------------------------------
// -----                  R3BGlobalAnalysisS454                -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosMappedData.h"

#include "R3BBeamMonitorMappedData.h"

#include "R3BGlobalAnalysisS454.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BRoluCalData.h"
#include "R3BRoluMappedData.h"

#include "R3BPaddleCalData.h"

#include "R3BPspxCalData.h"
#include "R3BPspxMappedData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberMappedData.h"

#include "R3BMCTrack.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include "tracker_routines.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BGlobalAnalysisS454::R3BGlobalAnalysisS454()
    : R3BGlobalAnalysisS454("GlobalAnalysis", 1)
{
}

R3BGlobalAnalysisS454::R3BGlobalAnalysisS454(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fCuts(0)
    , fGhost(0)
    , fPairs(0)
    , fB(-1672)
    , fNEvents(0)
{
}

R3BGlobalAnalysisS454::~R3BGlobalAnalysisS454()
{
    for (int i = 0; i < NOF_FIB_DET; i++)
    {
        delete fh_xy_Fib[i];
        delete fh_mult_Fib[i];
        delete fh_ToT_Fib[i];
    }
}

InitStatus R3BGlobalAnalysisS454::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(INFO) << "R3BGlobalAnalysisS454::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    FairRunOnline* run = FairRunOnline::Instance();

    // Get objects for detectors on all levels
    fMCTrack = (TClonesArray*)mgr->GetObject("MCTrack");
    fTrack = (TClonesArray*)mgr->GetObject("Track");
    maxevent = mgr->CheckMaxEventNo();

    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    printf("Have %d fiber detectors.\n", NOF_FIB_DET);
    fMappedItems.push_back((TClonesArray*)mgr->GetObject(Form("%sMappedData", fDetectorNames[0])));
    if (NULL == fMappedItems.at(0))
    {
        printf("Could not find mapped data for '%s'.\n", fDetectorNames[0]);
    }
    fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCrystalCalData", fDetectorNames[0])));
    if (NULL == fCalItems.at(0))
    {
        printf("Could not find Cal data for '%s'.\n", fDetectorNames[0]);
    }
    fHitItems.push_back((TClonesArray*)mgr->GetObject(Form("%sHit", fDetectorNames[0])));
    if (NULL == fHitItems.at(0))
    {
        printf("Could not find hit data for '%s'.\n", fDetectorNames[0]);
    }
    for (int det = 1; det < DET_MAX; det++)
    {
        fMappedItems.push_back((TClonesArray*)mgr->GetObject(Form("%sMapped", fDetectorNames[det])));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCal", fDetectorNames[det])));
        if (NULL == fCalItems.at(det))
        {
            printf("Could not find Cal data for '%s'.\n", fDetectorNames[det]);
        }
        fHitItems.push_back((TClonesArray*)mgr->GetObject(Form("%sHit", fDetectorNames[det])));
        if (NULL == fHitItems.at(det))
        {
            printf("Could not find hit data for '%s'.\n", fDetectorNames[det]);
        }
    }

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    char strNameC[255];
    sprintf(strNameC, "Cave_C_position");
    fh_Cave_position = new TH2F(strNameC, "", 2100, -100., 2000., 1000, -500., 500.);

    //-----------------------------------------------------------------------
    // BeamMonitor

    // get the theoretical calib factors for SEETRAM
    Double_t fexp = float(fsens_SEE + 9);
    Double_t fpow = float(pow(10., fexp));
    calib_SEE = 135641.7786 * fpow;
    LOG(DEBUG) << fsens_SEE << ", " << fexp << ", " << fpow << ", " << calib_SEE << endl;

    fh_Tpat = new TH1F("Tpat", "Tpat", 20, 0, 20);
    fh_Tpat->GetXaxis()->SetTitle("Tpat value");

    fh_Trigger = new TH1F("Trigger", "Trigger all", 20, 0, 20);
    fh_Trigger->GetXaxis()->SetTitle("Trigger value");

    fh_IC = new TH1F("IC", "IC ", 1000, 0, 1000);
    fh_IC->GetXaxis()->SetTitle("spill number");
    fh_IC->GetYaxis()->SetTitle("IC counts");

    fh_SEE = new TH1F("SEETRAM", "SEETRAM ", 1000, 0, 1000);
    fh_SEE->GetXaxis()->SetTitle("spill number");
    fh_SEE->GetYaxis()->SetTitle("SEETRAM counts");

    fh_TOFDOR = new TH1F("TOFDOR", "TOFDOR ", 1000, 0, 1000);
    fh_TOFDOR->GetXaxis()->SetTitle("spill number");
    fh_TOFDOR->GetYaxis()->SetTitle("TOFDOR counts");

    //-----------------------------------------------------------------------
    // compare against MC Simulations

    fh_target_px = new TH1F("target_px", "target px ", 2000, -500., 500);
    fh_target_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_target_px->GetYaxis()->SetTitle("counts");

    fh_target_py = new TH1F("target_py", "target py ", 2000, -500., 500);
    fh_target_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_target_py->GetYaxis()->SetTitle("counts");

    fh_target_pz = new TH1F("target_pz", "target pz ", 30000, 0., 30000.);
    fh_target_pz->GetXaxis()->SetTitle("pz / MeV/c");
    fh_target_pz->GetYaxis()->SetTitle("counts");

    fh_target_p = new TH1F("target_p", "target p ", 30000, 0., 30000.);
    fh_target_p->GetXaxis()->SetTitle("pz / MeV/c");
    fh_target_p->GetYaxis()->SetTitle("counts");

    fh_px_He = new TH1F("px_He", " px He", 2000, -500., 500);
    fh_px_He->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_He->GetYaxis()->SetTitle("counts");

    fh_py_He = new TH1F("py_He", " py He", 2000, -500., 500);
    fh_py_He->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_He->GetYaxis()->SetTitle("counts");

    fh_pz_He = new TH1F("pz_He", " pz He", 30000, 0., 30000.);
    fh_pz_He->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_He->GetYaxis()->SetTitle("counts");

    fh_p_He = new TH1F("p_He", " p He", 30000, 0., 30000.);
    fh_p_He->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_He->GetYaxis()->SetTitle("counts");

    fh_px_C = new TH1F("px_C", " px C", 2000, -500., 500);
    fh_px_C->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_C->GetYaxis()->SetTitle("counts");

    fh_py_C = new TH1F("py_C", " py C", 2000, -500., 500);
    fh_py_C->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_C->GetYaxis()->SetTitle("counts");

    fh_pz_C = new TH1F("pz_C", " pz C", 30000, 0., 30000.);
    fh_pz_C->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_C->GetYaxis()->SetTitle("counts");

    fh_p_C = new TH1F("p_C", " p C", 30000, 0., 30000.);
    fh_p_C->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_C->GetYaxis()->SetTitle("counts");

    fh_target_xy = new TH2F("target_xy", "target xy ", 1000, -0.5, 0.5, 1000, -0.5, 0.5);
    fh_target_xy->GetXaxis()->SetTitle("x / cm");
    fh_target_xy->GetYaxis()->SetTitle("y / cm");

    fh_chi2 = new TH1F("chi2", "chi2 ", 1000, 0., 100);
    fh_chi2->GetXaxis()->SetTitle("Chi2");
    fh_chi2->GetYaxis()->SetTitle("counts");

    fh_dx = new TH1F("tracker_dx", "tracker dx ", 2000, -5., 5);
    fh_dx->GetXaxis()->SetTitle("dx / cm");
    fh_dx->GetYaxis()->SetTitle("counts");

    fh_dy = new TH1F("tracker_dy", "tracker dy ", 2000, -5., 5.);
    fh_dy->GetXaxis()->SetTitle("dy / cm");
    fh_dy->GetYaxis()->SetTitle("counts");

    fh_dz = new TH1F("tracker_dz", "tracker dz ", 2000, -5., 5.);
    fh_dz->GetXaxis()->SetTitle("dz / cm");
    fh_dz->GetYaxis()->SetTitle("counts");

    fh_dpx = new TH1F("tracker_dpx", "tracker dpx ", 2000, -10, 10);
    fh_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_dpx->GetYaxis()->SetTitle("counts");

    fh_dpy = new TH1F("tracker_dpy", "tracker dpy ", 2000, -10, 10);
    fh_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_dpy->GetYaxis()->SetTitle("counts");

    fh_dpz = new TH1F("tracker_dpz", "tracker dpz ", 2000, -10, 10);
    fh_dpz->GetXaxis()->SetTitle("dpz / percent");
    fh_dpz->GetYaxis()->SetTitle("counts");

    fh_dp = new TH1F("tracker_dp", "tracker dp ", 2000, -10, 10);
    fh_dp->GetXaxis()->SetTitle("dp / percent");
    fh_dp->GetYaxis()->SetTitle("counts");

    fh_thetax_dpx = new TH2F("tracker_thethax_dpx", "tracker theta_x vs dpx ", 1000, -50, 50, 200, -100., 100.);
    fh_thetax_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_thetax_dpx->GetYaxis()->SetTitle("theta_x / mrad");

    fh_dpy_dpx = new TH2F("tracker_dpy_dpx", "tracker dpy vs dpx ", 200, -100, 100, 200, -100., 100.);
    fh_dpy_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_dpy_dpx->GetYaxis()->SetTitle("dpy / percent");

    fh_thetay_dpy = new TH2F("tracker_thethay_dpy", "tracker theta_y vs dpy ", 1000, -50, 50, 200, -100., 100.);
    fh_thetay_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_thetay_dpy->GetYaxis()->SetTitle("theta_y / mrad");

    fh_thetax_dpy = new TH2F("tracker_thethax_dpy", "tracker theta_x vs dpy ", 1000, -50, 50, 200, -100., 100.);
    fh_thetax_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_thetax_dpy->GetYaxis()->SetTitle("theta_x / mrad");

    fh_thetay_dpx = new TH2F("tracker_thethay_dpx", "tracker theta_y vs dpx ", 1000, -50, 50, 200, -100., 100.);
    fh_thetay_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_thetay_dpx->GetYaxis()->SetTitle("theta_y / mrad");

    fh_thetax_dpx_abs =
        new TH2F("tracker_thethax_dpx_abs", "tracker theta_x vs dpx abs ", 1000, -50, 50, 200, -100., 100.);
    fh_thetax_dpx_abs->GetXaxis()->SetTitle("dpx");
    fh_thetax_dpx_abs->GetYaxis()->SetTitle("theta_x / mrad");

    fh_x_dpx = new TH2F("tracker_x_dpx", "tracker x vs dpx ", 1000, -50, 50, 200, -10., 10.);
    fh_x_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_x_dpx->GetYaxis()->SetTitle("x / cm");

    fh_y_dpy = new TH2F("tracker_y_dpy", "tracker y vs dpy ", 1000, -50, 50, 200, -10., 10.);
    fh_y_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_y_dpy->GetYaxis()->SetTitle("y / cm");

    //-----------------------------------------------------------------------
    // Fiber Detectors 1-NOF_FIB_DET

    char canvName[255];
    UInt_t Nmax = 1e7;
    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {

        if (fCalItems.at(DET_FI_FIRST + ifibcount) || fHitItems.at(DET_FI_FIRST + ifibcount))
        {

            const char* detName;
            const char* detName2;
            detName = fDetectorNames[DET_FI_FIRST + ifibcount];

            LOG(DEBUG) << "I am creating canvas " << detName << endl;

            // xy:
            fh_xy_Fib[ifibcount] =
                new TH2F(Form("%s_xy", detName), Form("%s xy", detName), 600, -30., 30., 2000, -1000., 1000.);
            fh_xy_Fib[ifibcount]->GetXaxis()->SetTitle("x / cm ");
            fh_xy_Fib[ifibcount]->GetYaxis()->SetTitle("y / cm");

            fh_xy_Fib_ac[ifibcount] = new TH2F(
                Form("%s_xy_ac", detName), Form("%s xy after cuts", detName), 600, -30., 30., 2000, -100., 100.);
            fh_xy_Fib_ac[ifibcount]->GetXaxis()->SetTitle("x / cm ");
            fh_xy_Fib_ac[ifibcount]->GetYaxis()->SetTitle("y / cm");

            // Multiplicity (number of hit fibers):
            fh_mult_Fib[ifibcount] = new TH1F(Form("%s_mult", detName), Form("%s # of fibers", detName), 500, 0., 500.);
            fh_mult_Fib[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            fh_mult_Fib_ac[ifibcount] =
                new TH1F(Form("%s_mult_ac", detName), Form("%s # of fibers after cuts", detName), 500, 0., 500.);
            fh_mult_Fib_ac[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib_ac[ifibcount]->GetYaxis()->SetTitle("Counts");

            // ToT MAPMT:
            fh_ToT_Fib[ifibcount] =
                new TH2F(Form("%s_tot_m", detName), Form("%s ToT of MAPMT", detName), 1200, -30., 30, 4000, 0., 400.);
            fh_ToT_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber x / cm");
            fh_ToT_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            fh_ToT_Fib_ac[ifibcount] = new TH2F(Form("%s_tot_m_ac", detName),
                                                Form("%s ToT of MAPMT after cuts", detName),
                                                1200,
                                                -30.,
                                                30,
                                                4000,
                                                0.,
                                                400.);
            fh_ToT_Fib_ac[ifibcount]->GetXaxis()->SetTitle("Fiber x / cm");
            fh_ToT_Fib_ac[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // ToF Tofd -> Fiber:
            fh_Fib_ToF[ifibcount] = new TH2F(
                Form("%s_tof", detName), Form("%s ToF Tofd to Fiber", detName), 1200, -30., 30, 10000, -1000., 1000.);
            fh_Fib_ToF[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF[ifibcount]->GetXaxis()->SetTitle("x / cm");

            fh_Fib_ToF_ac[ifibcount] = new TH2F(Form("%s_tof_ac", detName),
                                                Form("%s ToF Tofd to Fiber after cuts", detName),
                                                1200,
                                                -30.,
                                                30,
                                                10000,
                                                -1000.,
                                                1000.);
            fh_Fib_ToF_ac[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF_ac[ifibcount]->GetXaxis()->SetTitle("x / cm");

            // ToF Tofd -> Fiber vs. event number:
            fh_ToF_vs_Events[ifibcount] = new TH2F(Form("%s_tof_vs_events", detName),
                                                   Form("%s ToF Tofd to Fiber vs event number", detName),
                                                   10000,
                                                   0,
                                                   Nmax,
                                                   2200,
                                                   -5100,
                                                   5100);
            fh_ToF_vs_Events[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_ToF_vs_Events[ifibcount]->GetXaxis()->SetTitle("event number");

            fh_ToF_vs_Events_ac[ifibcount] = new TH2F(Form("%s_tof_vs_events_ac", detName),
                                                      Form("%s ToF Tofd to Fiber vs event number after cuts", detName),
                                                      10000,
                                                      0,
                                                      Nmax,
                                                      2200,
                                                      -5100,
                                                      5100);
            fh_ToF_vs_Events_ac[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_ToF_vs_Events_ac[ifibcount]->GetXaxis()->SetTitle("event number");

            // hit fiber number vs. event number:
            fh_Fib_vs_Events[ifibcount] = new TH2F(Form("%s_fib_vs_event", detName),
                                                   Form("%s Fiber # vs. Event #", detName),
                                                   10000,
                                                   0,
                                                   Nmax,
                                                   1200,
                                                   -30.,
                                                   30.);
            fh_Fib_vs_Events[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fib_vs_Events[ifibcount]->GetXaxis()->SetTitle("Event number");

            fh_Fib_vs_Events_ac[ifibcount] = new TH2F(Form("%s_fib_vs_event_ac", detName),
                                                      Form("%s Fiber # vs. Event # after cuts", detName),
                                                      10000,
                                                      0,
                                                      Nmax,
                                                      1200,
                                                      -30.,
                                                      30.);
            fh_Fib_vs_Events_ac[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fib_vs_Events_ac[ifibcount]->GetXaxis()->SetTitle("Event number");

            // hit fiber number vs. TofD position:
            fh_Fibs_vs_Tofd[ifibcount] = new TH2F(Form("%s_fib_vs_TofdX", detName),
                                                  Form("%s Fiber # vs. Tofd x-pos", detName),
                                                  200,
                                                  -100,
                                                  100,
                                                  200,
                                                  -100.,
                                                  100.);
            fh_Fibs_vs_Tofd[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fibs_vs_Tofd[ifibcount]->GetXaxis()->SetTitle("Tofd x / cm");

            fh_Fibs_vs_Tofd_ac[ifibcount] = new TH2F(Form("%s_fib_vs_TofdX_ac", detName),
                                                     Form("%s Fiber # vs. Tofd x-pos after cuts", detName),
                                                     200,
                                                     -100,
                                                     100,
                                                     200,
                                                     -100.,
                                                     100.);
            fh_Fibs_vs_Tofd_ac[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fibs_vs_Tofd_ac[ifibcount]->GetXaxis()->SetTitle("Tofd x / cm");

            // hit fiber vs. fiber position:

        } // end if(Mapped)

    } // end for(ifibcount)
    fh_Fib13_vs_Fib11 = new TH2F("fib13_vs_fib11", "Fiber 13 vs. Fiber 11", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib13_vs_Fib11->GetYaxis()->SetTitle("Fiber13");
    fh_Fib13_vs_Fib11->GetXaxis()->SetTitle("Fiber11");

    fh_Fib11_vs_Fib3a = new TH2F("fib11_vs_fib3a", "Fiber 11 vs. Fiber 3a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib11_vs_Fib3a->GetYaxis()->SetTitle("Fiber11");
    fh_Fib11_vs_Fib3a->GetXaxis()->SetTitle("Fiber3a");

    fh_Fib10_vs_Fib12 = new TH2F("fib10_vs_fib12", "Fiber 10 vs. Fiber 12", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib10_vs_Fib12->GetYaxis()->SetTitle("Fiber10");
    fh_Fib10_vs_Fib12->GetXaxis()->SetTitle("Fiber12");

    fh_Fib12_vs_Fib3b = new TH2F("fib12_vs_fib3b", "Fiber 12 vs. Fiber 3b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib12_vs_Fib3b->GetYaxis()->SetTitle("Fiber12");
    fh_Fib12_vs_Fib3b->GetXaxis()->SetTitle("Fiber3b");

    // dx between fibers vs x
    fh_Fib13_vs_Fib11_dx = new TH2F("fib13_fib11_dx", "dx of Fiber 13 and Fiber 11", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib13_vs_Fib11_dx->GetYaxis()->SetTitle("xFi13 - xFi11 / cm");
    fh_Fib13_vs_Fib11_dx->GetXaxis()->SetTitle("x Fi11 / cm");

    fh_Fib11_vs_Fib3a_dx = new TH2F("fib11_fib3a_dx", "dx of Fiber 11 and Fiber 3a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib11_vs_Fib3a_dx->GetYaxis()->SetTitle("xFi11 - xFi3a / cm");
    fh_Fib11_vs_Fib3a_dx->GetXaxis()->SetTitle("x Fi3a / cm");

    fh_Fib10_vs_Fib12_dx = new TH2F("fib10_fib12_dx", "dx of Fiber 10 and Fiber 12", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib10_vs_Fib12_dx->GetYaxis()->SetTitle("xFi10 - xFi12 / cm");
    fh_Fib10_vs_Fib12_dx->GetXaxis()->SetTitle("x Fi12 / cm");

    fh_Fib12_vs_Fib3b_dx = new TH2F("fib12_fib3b_dx", "dx of Fiber 12 and Fiber 3b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib12_vs_Fib3b_dx->GetYaxis()->SetTitle("xFi12 - xFi3b / cm");
    fh_Fib12_vs_Fib3b_dx->GetXaxis()->SetTitle("x Fi3b / cm");

    //---------------------------------------------------------------------------------------------------
    // TofD detector

    if (fHitItems.at(DET_TOFD) || fCalItems.at(DET_TOFD))
    {

        // xy:
        fh_xy_tofd = new TH2F("tofd_xy", "tofd xy", 2000, -100., 100., 2000, -100., 100.);
        fh_xy_tofd->GetXaxis()->SetTitle("x / cm ");
        fh_xy_tofd->GetYaxis()->SetTitle("y / cm");

        fh_xy_tofd_ac = new TH2F("tofd_xy_ac", "tofd xy after cuts", 2000, -100., 100., 2000, -100., 100.);
        fh_xy_tofd_ac->GetXaxis()->SetTitle("x / cm ");
        fh_xy_tofd_ac->GetYaxis()->SetTitle("y / cm");

        fh_tofd_charge = new TH1F("tofd_Q", "Charge of Tofd", 200, 0., 20.);
        fh_tofd_charge->GetXaxis()->SetTitle("x / cm ");
        fh_tofd_charge->GetYaxis()->SetTitle("y / cm");

        fh_tofd_charge_ac = new TH1F("tofd_Q_ac", "Charge of Tofd after cuts", 200, 0., 20.);
        fh_tofd_charge_ac->GetXaxis()->SetTitle("x / cm ");
        fh_tofd_charge_ac->GetYaxis()->SetTitle("y / cm");

        fh_tofd_mult = new TH1F("tofd_mult", "ToFD multiplicits ", 100, 0, 100);
        fh_tofd_mult->GetXaxis()->SetTitle("multiplicity");
        fh_tofd_mult->GetYaxis()->SetTitle("counts");

        fh_tofd_mult_ac = new TH1F("tofd_mult_ac", "ToFD multiplicits after cuts", 100, 0, 100);
        fh_tofd_mult_ac->GetXaxis()->SetTitle("multiplicity");
        fh_tofd_mult_ac->GetYaxis()->SetTitle("counts");

        fh_TimePreviousEvent = new TH1F("TimePreviousEvent", "Time between 2 particles ", 300000, 0, 3000);
        fh_TimePreviousEvent->GetXaxis()->SetTitle("time / ns");
        fh_TimePreviousEvent->GetYaxis()->SetTitle("counts");

        fh_tofd_q2_vs_q1 = new TH2F("tofd_q2_vs_q1", "tofd q2 vs. q1", 500, 0., 50., 500, 0., 50.);
        fh_tofd_q2_vs_q1->GetXaxis()->SetTitle("q1");
        fh_tofd_q2_vs_q1->GetYaxis()->SetTitle("q2");

        fh_tofd_q2_vs_q1_ac = new TH2F("tofd_q2_vs_q1_ac", "tofd q2 vs. q1 after cut", 500, 0., 50., 500, 0., 50.);
        fh_tofd_q2_vs_q1_ac->GetXaxis()->SetTitle("q1");
        fh_tofd_q2_vs_q1_ac->GetYaxis()->SetTitle("q2");
    }
    // Analysis
    fh_chiy_vs_chix = new TH2F("chiy_vs_chix", "chi y vs. q1", 500, 0., 500., 500, 0., 500.);
    fh_chiy_vs_chix->GetXaxis()->SetTitle("chi x");
    fh_chiy_vs_chix->GetYaxis()->SetTitle("chi y");

    fh_theta26 = new TH1F("theta26", "theta 26 ", 500, 0., 5);
    fh_theta26->GetXaxis()->SetTitle("angle / degree");
    fh_theta26->GetYaxis()->SetTitle("counts");

    fh_theta_16O = new TH1F("theta_16O", "theta 16O ", 500, 0., 5);
    fh_theta_16O->GetXaxis()->SetTitle("angle / degree");
    fh_theta_16O->GetYaxis()->SetTitle("counts");

    fh_theta26_cm = new TH1F("theta26_cm", "theta 26 cm ", 180, -90., 90.);
    fh_theta26_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_cm->GetYaxis()->SetTitle("counts");

    fh_phi26_cm = new TH1F("phi26_cm", "phi 26 cm", 360, 0., 360.);
    fh_phi26_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi26_cm->GetYaxis()->SetTitle("counts");

    fh_theta_4He_cm = new TH1F("theta_4He_cm", "Theta 4He cm", 180, 0., 180.);
    fh_theta_4He_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_4He_cm->GetYaxis()->SetTitle("counts");

    fh_phi_4He_cm = new TH1F("phi_4He_cm", "phi 4He cm", 360, 0., 360.);
    fh_phi_4He_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_4He_cm->GetYaxis()->SetTitle("counts");

    fh_theta_12C_cm = new TH1F("theta_12C_cm", "Theta 12C cm", 180, 0., 180.);
    fh_theta_12C_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_12C_cm->GetYaxis()->SetTitle("counts");

    fh_phi_12C_cm = new TH1F("phi_12C_cm", "phi 12C cm", 360, 0., 360.);
    fh_phi_12C_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_12C_cm->GetYaxis()->SetTitle("counts");

    fh_theta_16O_cm = new TH1F("theta_16O_cm", "Theta 16O cm", 180, 0., 180.);
    fh_theta_16O_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_16O_cm->GetYaxis()->SetTitle("counts");

    fh_phi_16O_cm = new TH1F("phi_16O_cm", "phi 16O cm", 360, 0., 360.);
    fh_phi_16O_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_16O_cm->GetYaxis()->SetTitle("counts");

    fh_theta26_simu = new TH1F("theta26_simu", "theta 26 simulation", 500, 0., 5);
    fh_theta26_simu->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_simu->GetYaxis()->SetTitle("counts");

    fh_Erel = new TH1F("Erel", "Erel ", 5000, 0., 100);
    fh_Erel->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel->GetYaxis()->SetTitle("counts");

    fh_ErelL = new TH1F("ErelL", "Erel links", 5000, 0., 100);
    fh_ErelL->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelL->GetYaxis()->SetTitle("counts");

    fh_ErelR = new TH1F("ErelR", "Erel rechts ", 5000, 0., 100);
    fh_ErelR->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelR->GetYaxis()->SetTitle("counts");

    fh_Erel_simu = new TH1F("Erel_simu", "Erel simulation", 6000, -10., 50);
    fh_Erel_simu->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel_simu->GetYaxis()->SetTitle("counts");

    fh_dErel_vs_x = new TH2F("dErel_vs_x", "delta Erel vs. x", 200, -100., 100., 100, -5., 5.);
    fh_dErel_vs_x->GetXaxis()->SetTitle("TofD x / cm");
    fh_dErel_vs_x->GetYaxis()->SetTitle("Delta Erel / MeV");

    fh_dErel_vs_y = new TH2F("dErel_vs_y", "delta Erel vs. y", 200, -100., 100., 100, -5., 5.);
    fh_dErel_vs_y->GetXaxis()->SetTitle("TofD y / cm");
    fh_dErel_vs_y->GetYaxis()->SetTitle("Delta Erel / MeV");

    return kSUCCESS;
}

void R3BGlobalAnalysisS454::Exec(Option_t* option)
{
    if (fNEvents / 100000. == (int)fNEvents / 100000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                  << " %) " << std::flush;

    // cout << "New event ******************************" << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    pHex = -100.;
    pHey = -100.;
    pHez = -100.;
    pCx = -100.;
    pCy = -100.;
    pCz = -100.;
    XHe = -100;
    YHe = -100;
    ZHe = -100;
    XC = -100;
    YC = -100;
    ZC = -100;
    Double_t px, py, pz;
    Double_t theta_16, theta_26;
    Double_t costh26;
    Double_t chiHex, chiHey, chiCx, chiCy;

    if (fTrack)
    {
        // read in Track parameter

        Int_t nHitsTrack = fTrack->GetEntriesFast();
        // cout << "Track hits: " << nHitsTrack << endl;

        for (Int_t l = 0; l < nHitsTrack; l++)
        {
            R3BTrack* aTrack = (R3BTrack*)fTrack->At(l);

            LOG(DEBUG) << "Charge " << aTrack->GetQ() << endl;

            if (aTrack->GetQ() == 2)
            {
                // 4He
                XHe = aTrack->GetX();
                YHe = aTrack->GetY();
                ZHe = aTrack->GetZ();

                pHex = aTrack->GetPx();
                pHey = aTrack->GetPy();
                pHez = aTrack->GetPz();

                chiHex = aTrack->GetChix();
                chiHey = aTrack->GetChiy();

                LOG(DEBUG) << "******************************************" << endl;
                LOG(DEBUG) << "Track In 4He"
                           << "x " << XHe << " y " << YHe << " z " << ZHe << endl;
                LOG(DEBUG) << "px " << pHex << " py " << pHey << " z " << pHez << endl;
            }
            if (aTrack->GetQ() == 6)
            {
                // 12C
                XC = aTrack->GetX();
                YC = aTrack->GetY();
                ZC = aTrack->GetZ();

                pCx = aTrack->GetPx();
                pCy = aTrack->GetPy();
                pCz = aTrack->GetPz();

                chiCx = aTrack->GetChix();
                chiCy = aTrack->GetChiy();

                LOG(DEBUG) << "******************************************" << endl;
                LOG(DEBUG) << "Track In 12C"
                           << "x " << XC << " y " << YC << " z " << ZC << endl;
                LOG(DEBUG) << "px " << pCx << " py " << pCy << " z " << pCz << endl;
            }
            if (aTrack->GetQ() == 8)
            {
                // 160
                Xf = aTrack->GetX();
                Yf = aTrack->GetY();
                Zf = aTrack->GetZ();

                Pxf = aTrack->GetPx();
                Pyf = aTrack->GetPy();
                Pzf = aTrack->GetPz();

                LOG(DEBUG) << "******************************************" << endl;
                LOG(DEBUG) << "Track In 16O"
                           << "x " << Xf << " y " << Yf << " z " << Zf << endl;
                LOG(DEBUG) << "px " << Pxf << " py " << Pyf << " z " << Pzf << endl;
            }
        }

        if (nHitsTrack > 1)
        {
            fh_chiy_vs_chix->Fill(chiHex, chiHey);
            fh_chiy_vs_chix->Fill(chiCx, chiCy);

            if (chiHex < 100 && chiHey < 100 && chiCx < 100 && chiCy < 100)
            {
                fh_target_xy->Fill(XHe * 100., YHe * 100.);

                fh_px_He->Fill(pHex);
                fh_py_He->Fill(pHey);
                fh_pz_He->Fill(pHez);

                fh_px_C->Fill(pCx);
                fh_py_C->Fill(pCy);
                fh_pz_C->Fill(pCz);

                // Calculate angle between alphs and C
                if (pCz == 0 || pHez == 0)
                    return;

                costh26 = (pHex * pCx + pHey * pCy + pHez * pCz) /
                          (sqrt(pHex * pHex + pHey * pHey + pHez * pHez) * sqrt(pCx * pCx + pCy * pCy + pCz * pCz));

                theta_26 = acos(costh26) * 180. / 3.14159; // opening angle

                fh_theta26->Fill(theta_26);

                Double_t pHe = sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2));
                Double_t pC = sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2));
                Double_t eHe = sqrt(pow(pHe, 2) + pow(mHe, 2)) - mHe;
                Double_t eC = sqrt(pow(pC, 2) + pow(mC, 2)) - mC;

                fh_p_He->Fill(pHe);
                fh_p_C->Fill(pC);

                Double_t Erela = sqrt(pow((mHe + mC + eHe + eC), 2) - pow(pHe, 2) - pow(pC, 2) -
                                      2 * pHe * pC * cos(theta_26 * 3.1415 / 180.)) -
                                 (mHe + mC); // Erel

                Double_t ptota = sqrt(pow(pHe, 2) + pow(pC, 2) + 2 * pHe * pC * cos(theta_26 * 3.1415 / 180.));

                Double_t m_inv = sqrt(mHe * mHe + mC * mC + 2. * sqrt(pC * pC + mC * mC) * sqrt(pHe * pHe + mHe * mHe) -
                                      2. * pHe * pC * cos(theta_26 * 3.1415 / 180.));
                Double_t Erelb = m_inv - mHe - mC;

                fh_Erel->Fill(Erela);
                if (pHex < 0. && pCx < 0.)
                    fh_ErelR->Fill(Erela);
                if (pHex > 0. && pCx > 0.)
                    fh_ErelL->Fill(Erela);

                //	fh_dErel_vs_x->Fill(xTest * 100., Erela - 1.);
                //	fh_dErel_vs_y->Fill(yTest * 100., Erela - 1.);
                // cout << "Theta 26: " << theta_26 << " Erel: " << Erela << " " << Erelb << endl;

                // transfer to cm system and make some rotations

                px = pHex + pCx;
                py = pHey + pCy;
                pz = pHez + pCz;

                if (pz != 0.)
                {
                    theta_16 = atan2(sqrt(pow(px, 2) + pow(py, 2)), pz) * 180. / 3.14159265;
                }

                fh_theta_16O->Fill(theta_16);

                Double_t n2[10], eaa, vaa, g2, v2[10], lt2[10][10];
                Double_t tmp1[10], tmp2[10], tmp3[10], tmp4[10];
                Double_t pxO_1, pyO_1, pzO_1, esum_3;
                Double_t psum_1, thcm, phcm, tanphicm;

                esum_3 = sqrt(pow((sqrt(pow(pHe, 2) + pow(pC, 2) + 2 * pHe * pC * cos(theta_26 * 3.14159 / 180.))), 2) +
                              pow(mO, 2)) -
                         mO;

                psum_1 = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
                thcm = acos(pz / psum_1);
                tanphicm = py / px;

                if (px > 0)
                {
                    phcm = atan(tanphicm) * 180. / 3.14159265;
                }

                if (px < 0. && py > 0.)
                {
                    phcm = (atan(tanphicm) * 180. / 3.14159265) + 180.;
                }

                if (px < 0. && py < 0.)
                {
                    phcm = (atan(tanphicm) * 180. / 3.14159265) - 180.;
                }

                if (phcm < 0)
                {
                    phcm = phcm + 360.;
                }
                phcm = phcm * 3.14159265 / 180.;

                fh_theta_16O_cm->Fill(thcm);
                fh_phi_16O_cm->Fill(phcm);

                n2[1] = px / ptota;
                n2[2] = py / ptota;
                n2[3] = pz / ptota;

                //       eaa   = esum_3;
                //     velocity of 6Li
                //       vaa=sqrt(eaa*eaa+2*mO*eaa)/(eaa+mO);

                //       cout<< " va1 "<< vaa<<endl;

                eaa = sqrt(mO * mO + psum_1 * psum_1);
                vaa = psum_1 / eaa;

                //       cout<< " va2 "<< vaa<<endl;

                g2 = 1.0 / sqrt(1.0 - vaa * vaa);

                v2[1] = vaa * n2[1];
                v2[2] = vaa * n2[2];
                v2[3] = vaa * n2[3];

                // Lorentz transformation
                lt2[1][1] = g2;
                lt2[1][2] = -1. * g2 * v2[1];
                lt2[1][3] = -1. * g2 * v2[2];
                lt2[1][4] = -1. * g2 * v2[3];
                lt2[2][1] = lt2[1][2];
                lt2[2][2] = 1.0 + (g2 - 1.0) * n2[1] * n2[1];
                lt2[2][3] = (g2 - 1.0) * n2[1] * n2[2];
                lt2[2][4] = (g2 - 1.0) * n2[1] * n2[3];
                lt2[3][1] = lt2[1][3];
                lt2[3][2] = lt2[2][3];
                lt2[3][3] = 1.0 + (g2 - 1.0) * n2[2] * n2[2];
                lt2[3][4] = (g2 - 1.0) * n2[2] * n2[3];
                lt2[4][1] = lt2[1][4];
                lt2[4][2] = lt2[2][4];
                lt2[4][3] = lt2[3][4];
                lt2[4][4] = 1.0 + (g2 - 1.0) * n2[3] * n2[3];

                //      four momentum of nucleus b
                tmp1[1] = sqrt(mHe * mHe + pHe * pHe);
                tmp1[2] = pHex;
                tmp1[3] = pHey;
                tmp1[4] = pHez;

                //     four momentum of nucleus c
                tmp2[1] = sqrt(mC * mC + pC * pC);
                tmp2[2] = pCx;
                tmp2[3] = pCy;
                tmp2[4] = pCz;

                //     four momentum of nucleus a for a check of the transformation
                tmp4[1] = sqrt(mO * mO + psum_1 * psum_1);
                tmp4[2] = px;
                tmp4[3] = py;
                tmp4[4] = pz;
                //       cout <<"compare " << tmp4[1]<<"  "<< esum_3+mO<<endl;

                /*
                                cout <<"in " << tmp4[1]<<endl;
                                cout <<"in " << tmp4[2]<<endl;
                                cout <<"in " << tmp4[3]<<endl;
                                cout <<"in " << tmp4[4]<<endl;
                */
                for (Int_t i1 = 1; i1 < 5; i1++)
                {
                    tmp3[i1] = 0.0;
                    for (Int_t i2 = 1; i2 < 5; i2++)
                    {
                        tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp1[i2];
                    }
                }

                Double_t pHex_1, pHey_1, pHez_1, pHe_1, pCx_1, pCy_1, pCz_1, pC_1;

                pHex_1 = tmp3[2];
                pHey_1 = tmp3[3];
                pHez_1 = tmp3[4];

                for (Int_t i1 = 1; i1 < 5; i1++)
                {
                    tmp3[i1] = 0.0;
                    for (Int_t i2 = 1; i2 < 5; i2++)
                    {
                        tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp2[i2];
                    }
                }

                pCx_1 = tmp3[2];
                pCy_1 = tmp3[3];
                pCz_1 = tmp3[4];

                for (Int_t i1 = 1; i1 < 5; i1++)
                {
                    tmp3[i1] = 0.0;
                    for (Int_t i2 = 1; i2 < 5; i2++)
                    {
                        tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp4[i2];
                    }
                }

                pxO_1 = tmp3[2];
                pyO_1 = tmp3[3];
                pzO_1 = tmp3[4];

                /*
                                cout <<"out " << tmp3[1]<<endl;
                                cout <<"out " << tmp3[2]<<endl;
                                cout <<"out " << tmp3[3]<<endl;
                                cout <<"out " << tmp3[4]<<endl;
                */

                Double_t pHex_cm, pHey_cm, pHez_cm, pHe_cm;
                Double_t pCx_cm, pCy_cm, pCz_cm, pC_cm;
                Double_t thet_cm_He, thet_cm_C, thet;
                Double_t theta_cm_He, theta_cm_C, theta;
                Double_t pxyd_cm, costhetacm_He, tanphi_d, phicm_He, pxyhe_cm;
                Double_t costhetacm_C, tanphi_he, phicm_C, theta_cm_26;
                Double_t phi_cm_26, pxO_cm, pyO_cm, pzO_cm;

                //     Vector for rotation
                n2[1] = py / sqrt(px * px + py * py);
                n2[2] = -px / sqrt(px * px + py * py);
                n2[3] = 0;

                //       cout <<"winkel " << thcm<<"  "<<theta6/180.*3.14159265<<endl;

                lt2[1][1] = cos(thcm) + n2[1] * n2[1] * (1. - cos(thcm));
                lt2[1][2] = n2[1] * n2[2] * (1. - cos(thcm)) - n2[3] * sin(thcm);
                lt2[1][3] = n2[1] * n2[3] * (1. - cos(thcm)) + n2[2] * sin(thcm);
                lt2[2][1] = n2[2] * n2[1] * (1. - cos(thcm)) + n2[3] * sin(thcm);
                lt2[2][2] = cos(thcm) + n2[2] * n2[2] * (1. - cos(thcm));
                lt2[2][3] = n2[2] * n2[3] * (1. - cos(thcm)) - n2[1] * sin(thcm);
                lt2[3][1] = n2[3] * n2[1] * (1. - cos(thcm)) - n2[2] * sin(thcm);
                lt2[3][2] = n2[3] * n2[2] * (1. - cos(thcm)) + n2[1] * sin(thcm);
                lt2[3][3] = cos(thcm) + n2[3] * n2[3] * (1. - cos(thcm));

                //      four momentum of nucleus b
                tmp1[1] = pHex_1;
                tmp1[2] = pHey_1;
                tmp1[3] = pHez_1;

                //     four momentum of nucleus c
                tmp2[1] = pCx_1;
                tmp2[2] = pCy_1;
                tmp2[3] = pCz_1;

                //     four momentum of nucleus a for a check of the transformation
                tmp4[1] = px;
                tmp4[2] = py;
                tmp4[3] = pz;
                //       cout <<"compare " << tmp4[1]<<"  "<< esum_3+mO<<endl;

                /*
                                cout <<"in1 " << tmp4[1]<<endl;
                                cout <<"in1 " << tmp4[2]<<endl;
                                cout <<"in1 " << tmp4[3]<<endl;
                */
                Double_t pHex_2, pHey_2, pHez_2, pHe_2, pCx_2, pCy_2, pCz_2, pC_2;
                Double_t pxO_2, pyO_2, pzO_2;

                for (Int_t i1 = 1; i1 < 4; i1++)
                {
                    tmp3[i1] = 0.0;
                    for (Int_t i2 = 1; i2 < 4; i2++)
                    {
                        tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp1[i2];
                    }
                }

                pHex_cm = tmp3[1];
                pHey_cm = tmp3[2];
                pHez_cm = tmp3[3];

                for (Int_t i1 = 1; i1 < 4; i1++)
                {
                    tmp3[i1] = 0.0;
                    for (Int_t i2 = 1; i2 < 4; i2++)
                    {
                        tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp2[i2];
                    }
                }

                pCx_cm = tmp3[1];
                pCy_cm = tmp3[2];
                pCz_cm = tmp3[3];

                for (Int_t i1 = 1; i1 < 4; i1++)
                {
                    tmp3[i1] = 0.0;
                    for (Int_t i2 = 1; i2 < 4; i2++)
                    {
                        tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp4[i2];
                    }
                }

                pxO_cm = tmp3[1];
                pyO_cm = tmp3[2];
                pzO_cm = tmp3[3];

                /*
                                cout <<"out1 " << tmp3[1]<<endl;
                                cout <<"out1 " << tmp3[2]<<endl;
                                cout <<"out1 " << tmp3[3]<<endl;
                */

                // 4He fragment

                pHe_cm = sqrt(pHex_cm * pHex_cm + pHey_cm * pHey_cm + pHez_cm * pHez_cm);

                costhetacm_He = pHez_cm / pHe_cm;
                theta_cm_He = acos(costhetacm_He) * 180. / 3.14159265;
                tanphi_d = pHey_cm / pHex_cm;

                if (pHex_cm > 0)
                {
                    phicm_He = atan(tanphi_d) * 180. / 3.14159265;
                }

                if (pHex_cm < 0. && pHey_cm > 0.)
                {
                    phicm_He = (atan(tanphi_d) * 180. / 3.14159265) + 180.;
                }

                if (pHex_cm < 0. && pHey_cm < 0.)
                {
                    phicm_He = (atan(tanphi_d) * 180. / 3.14159265) - 180.;
                }

                if (phicm_He < 0)
                {
                    phicm_He = phicm_He + 360.;
                }
                fh_theta_4He_cm->Fill(theta_cm_He);
                fh_phi_4He_cm->Fill(phicm_He);

                // 12C fragment

                pC_cm = sqrt(pCx_cm * pCx_cm + pCy_cm * pCy_cm + pCz_cm * pCz_cm);

                costhetacm_C = pCz_cm / pC_cm;
                theta_cm_C = acos(costhetacm_C) * 180. / 3.14159265;
                tanphi_he = pCy_cm / pCx_cm;

                if (pCx_cm > 0.)
                {
                    phicm_C = atan(tanphi_he) * 180. / 3.14159265;
                }

                if (pCx_cm < 0. && pCy_cm > 0.)
                {
                    phicm_C = (atan(tanphi_he) * 180. / 3.14159265) + 180.;
                }

                if (pCx_cm < 0. && pCy_cm < 0.)
                {
                    phicm_C = (atan(tanphi_he) * 180. / 3.14159265) - 180.;
                }

                if (phicm_C < 0.)
                {
                    phicm_C = phicm_C + 360.;
                }
                fh_theta_12C_cm->Fill(theta_cm_C);
                fh_phi_12C_cm->Fill(phicm_C);

                theta_cm_26 = theta_cm_He + theta_cm_C;
                phi_cm_26 = phicm_He - phicm_C;

                fh_theta26_cm->Fill(theta_cm_26);
                fh_phi26_cm->Fill(phi_cm_26);

            } // end if chi2
        }     // end if trackHits>1
    }         // end if fTrack

    fNEvents += 1;
    return;
    // --------------------------------------------------------------------------------
    if (header)
    {
        time = header->GetTimeStamp();

        if (time_start == 0 && time > 0)
        {
            time_start = time;
            fNEvents_start = fNEvents;
        }

        if (header->GetTrigger() == 12)
        {
            // spill start in nsec
            // cout << "spill start" << endl;
            num_spills++;
        }
        if (header->GetTrigger() == 13)
        {
            // spill end  in nsec
            // cout << "spill stop" << endl;
        }

        fh_Trigger->Fill(header->GetTrigger());
        //   check for requested trigger (Todo: should be done globablly / somewhere else)
        if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
            return;

        Int_t tpatbin;
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0)
                fh_Tpat->Fill(i + 1);
        }

        // fTpat = 1-16; fTpat_bit = 0-15
        Int_t fTpat_bit = fTpat - 1;
        Int_t itpat;
        Int_t tpatvalue;
        if (fTpat_bit >= 0)
        {
            itpat = header->GetTpat();
            tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
            if (tpatvalue == 0)
                return;
        }
    }

    // cout << "Spill: " << num_spills << endl;

    if (fMappedItems.at(DET_BMON))
    {
        unsigned long IC;
        unsigned long SEETRAM_raw;
        Double_t SEETRAM;
        unsigned long TOFDOR;

        auto detBmon = fMappedItems.at(DET_BMON);
        Int_t nHitsbm = detBmon->GetEntriesFast();
        // cout<<"Bmon hits: "<<nHitsbm<<endl;

        for (Int_t ihit = 0; ihit < nHitsbm; ihit++)
        {
            R3BBeamMonitorMappedData* hit = (R3BBeamMonitorMappedData*)detBmon->At(ihit);
            if (!hit)
                continue;

            IC = hit->GetIC(); // negative values if offset not high enough
            counts_IC += (double)IC;

            SEETRAM_raw = hit->GetSEETRAM();           // raw counts
            SEETRAM = (double)SEETRAM_raw * calib_SEE; // calibrated SEETRAM counts
            counts_SEE += SEETRAM;

            TOFDOR = hit->GetTOFDOR(); // only positive values possible
            counts_TofD += TOFDOR;

            if (fNEvents == fNEvents_start)
            {
                see_start = SEETRAM;
                ic_start = IC;
                tofdor_start = TOFDOR;
            }

            /*
            cout << "time " << time << endl;
            cout << "IC   " << IC << "  " << counts_IC << "  " << endl;
            cout << "SEE  " << SEETRAM_raw << "  " << counts_SEE << "  " << SEETRAM << endl;
            cout << "TofD " << TOFDOR << "  " << counts_TofD << "  " << endl;
            */

            // IC:
            Int_t yIC = IC - ic_start;
            fh_IC->Fill(num_spills, yIC);

            // SEETRAM:
            Int_t ySEE = SEETRAM - see_start;
            fh_SEE->Fill(num_spills, ySEE);
            // Double_t ySEE_part = (SEETRAM-see_mem)*fNorm*1.e+3-see_offset*calib_SEE;

            // TOFDOR:
            Int_t yTOFDOR = TOFDOR - tofdor_start;
            fh_TOFDOR->Fill(num_spills, yTOFDOR);
        }
    }

    Bool_t RoluCut = false;
    if (fMappedItems.at(DET_ROLU))
    {
        // rolu
        auto detRolu = fMappedItems.at(DET_ROLU);
        Int_t nHitsRolu = detRolu->GetEntriesFast();
        // cout<<"ROLU hits: "<<nHitsRolu<<endl;

        for (Int_t ihit = 0; ihit < nHitsRolu; ihit++)
        {
            R3BRoluMappedData* hitRolu = (R3BRoluMappedData*)detRolu->At(ihit);
            if (!hitRolu)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hitRolu->GetDetector(); // 1..
            Int_t iCha = hitRolu->GetChannel();  // 1..
            RoluCut = true;
        }
    }
    if (RoluCut)
    {
        // cout << "ROLU cut applied !!!" << endl;
        return;
    }

    Double_t xTest = 0.;
    Double_t yTest = 0.;

    Int_t max = 10000;
    Int_t detector[max];
    Double_t xdet[max];
    Double_t ydet[max];
    Double_t zdet[max];
    Int_t qdet[max];

    countdet = 0;

    Double_t track[12]; // x,y,z, px, py, pz
    Double_t chi[6];    // x,y,z, px, py, pz

    Int_t n_det = 8;
    if (fGhost)
        n_det = 9;

    Double_t x[n_det];
    Double_t y[n_det];
    Double_t z[n_det];
    Double_t q[n_det];
    Double_t t[n_det];
    Double_t x1[n_det];
    Double_t y1[n_det];
    Double_t z1[n_det];
    Double_t q1[n_det];
    Double_t t1[n_det];
    Double_t x2[n_det];
    Double_t y2[n_det];
    Double_t z2[n_det];
    Double_t q2[n_det];
    Double_t t2[n_det];
    Double_t xMax[n_det];
    Double_t yMax[n_det];
    Double_t zMax[n_det];
    Double_t qMax[n_det];
    Double_t tMax[n_det];

    Bool_t pat1[2 * n_det];
    Bool_t pat2[2 * n_det];

    Int_t id, id1, id2;

    Int_t det = 0;
    Int_t det1 = 0;
    Int_t det2 = 0;

    Double_t target[3];
    Double_t track1_best[6]; // x,y,z, px, py, pz
    Double_t track2_best[6]; // x,y,z, px, py, pz
    Double_t chi_best[6];    // chi2, chi2_red

    Double_t chi2;
    Double_t chi2_best;

    target[0] = 0.;
    target[1] = 0.;
    target[2] = 0.;

    Double_t track_parameter[6];        // x,y,z, px, py, pz
    Double_t chi_single_parameter[2];   // chi2, chi2_red
    Double_t track_parameter_best_l[6]; // x,y,z, px, py, pz
    Double_t track_parameter_best_r[6]; // x,y,z, px, py, pz
    Double_t chi2Minl[2];
    Double_t chi2Minr[2];

    for (int i = 0; i < n_det; i++)
    {
        x[i] = 0.;
        y[i] = 0.;
        z[i] = 0.;
        q[i] = 0.;
        t[i] = -1000.;

        x1[i] = 0.;
        y1[i] = 0.;
        z1[i] = 0.;
        q1[i] = 0.;
        t1[i] = -1000.;

        x2[i] = 0.;
        y2[i] = 0.;
        z2[i] = 0.;
        q2[i] = 0.;
        t2[i] = -1000.;

        xMax[i] = -1000.;
        yMax[i] = -1000.;
        zMax[i] = -1000.;
        qMax[i] = -1000.;
        tMax[i] = -1000.;
    }
    pHex = 0.;
    pHey = 0.;
    pHez = 0.;
    pCx = 0.;
    pCy = 0.;
    pCz = 0.;

    // is also number of ifibcount
    Int_t fi3a = 0;
    Int_t fi3b = 1;
    Int_t fi10 = 2;
    Int_t fi11 = 3;
    Int_t fi12 = 4;
    Int_t fi13 = 5;
    Int_t tofdr = 6;
    Int_t tofdl = 7;
    Int_t ghost = 8;

    Double_t tof = 0.;
    Bool_t pair = false;

    //    ofstream myfile;
    //    myfile.open ("events_v4.dat", ios::out | ios::app);
    for (Int_t i = 0; i < 10; i++)
    {
        tPrev[i] = -1000.;
        detPrev[i] = -1;
    }

    Bool_t coinc = false;

    auto detTofd = fHitItems.at(DET_TOFD);
    Int_t nHits = detTofd->GetEntriesFast();
    LOG(DEBUG) << "********************************" << endl;
    LOG(DEBUG) << "ToFD hits: " << nHits << endl;

    if (nHits > 0)
        fh_tofd_mult->Fill(nHits);

    if (nHits > 100)
        return;

    // loop over ToFD
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        R3BTofdHitData* hitTofd = (R3BTofdHitData*)detTofd->At(ihit);
        pair = false;
        if (fCuts && (hitTofd->GetTime() < 0. || hitTofd->GetTime() > 8192. * 5.))
            continue;
        if (IS_NAN(hitTofd->GetTime()))
            continue;
        if (hitTofd->GetX() / 100. <= -0.007)
        {
            // tof rechts
            det2 = tofdr;
            det = tofdr;
        }
        else
        {
            // tof links
            det2 = tofdl;
            det = tofdl;
        }
        x2[det2] = hitTofd->GetX() / 100.;
        y2[det2] = hitTofd->GetY() / 100.;
        z2[det2] = 0.;
        q2[det2] = hitTofd->GetEloss();
        t2[det2] = hitTofd->GetTime();
        id2 = hitTofd->GetDetId();

        x[det] = hitTofd->GetX() / 100.;
        y[det] = hitTofd->GetY() / 100.;
        z[det] = 0.;
        q[det] = hitTofd->GetEloss();
        t[det] = hitTofd->GetTime();
        id = hitTofd->GetDetId();

        LOG(DEBUG) << "t1 " << t1[det1] << endl;
        LOG(DEBUG) << "t2 " << t2[det2] << endl;

        // for 2 particle tracking, look for pairs
        if (abs(t2[det2] - t1[det1]) < 2. && fPairs)
        {
            // now we have two particle within a short time

            LOG(DEBUG) << "TofD ID " << id1 << "  " << id2 << endl;

            // check for good hits in tofd
            // if (fCuts && (id1 != 12 || id2 != 12))
            // continue;

            fh_tofd_q2_vs_q1->Fill(q1[det1], q2[det2]);
            LOG(DEBUG) << "q1 " << q1[det1] << endl;
            LOG(DEBUG) << "q2 " << q2[det2] << endl;

            if ((abs(q1[det1] - 2.) < 0.5 && abs(q2[det2] - 6.) < 0.5) ||
                (abs(q2[det2] - 2.) < 0.5 && abs(q1[det1] - 6.) < 0.5))
            {
                // Find pairs of 4He and 12C
                pair = true;

                countdet = 0;
                chi2_best = 1.E10;

                fh_tofd_q2_vs_q1_ac->Fill(q1[det1], q2[det2]);

                LOG(DEBUG) << endl;
                LOG(DEBUG) << "ToFD 1"
                           << " x: " << x1[det1] << " y: " << y1[det1] << " q: " << q1[det1] << " t: " << t1[det1]
                           << endl;

                LOG(DEBUG) << "ToFD 2"
                           << " x: " << x2[det2] << " y: " << y2[det2] << " q: " << q2[det2] << " t: " << t2[det2]
                           << endl;

                detector[countdet] = det1;
                xdet[countdet] = x1[det1];
                ydet[countdet] = y1[det1];
                zdet[countdet] = z1[det1];
                qdet[countdet] = (int)(q1[det1] + 0.5);
                countdet++;
                detector[countdet] = det2;
                xdet[countdet] = x2[det2];
                ydet[countdet] = y2[det2];
                zdet[countdet] = z2[det2];
                qdet[countdet] = (int)(q2[det2] + 0.5);
                countdet++;

                xTest = x1[det1];
                yTest = y1[det1];

                fh_xy_tofd->Fill(x1[det1] * 100., y1[det1] * 100.);
                fh_xy_tofd->Fill(x2[det2] * 100., y2[det2] * 100.);
                fh_tofd_charge->Fill(q1[det1]);
                fh_tofd_charge->Fill(q2[det2]);

                // fill histograms for ToFD
                fh_tofd_mult_ac->Fill(nHits);
                fh_xy_tofd_ac->Fill(x1[det1] * 100., y1[det1] * 100.);
                fh_xy_tofd_ac->Fill(x2[det2] * 100., y2[det2] * 100.);
                fh_tofd_charge_ac->Fill(q1[det1]);
                fh_tofd_charge_ac->Fill(q2[det2]);
            }
        }
        else
        {
            det1 = det2;
            x1[det1] = x2[det2];
            y1[det1] = y2[det2];
            z1[det1] = 0.;
            q1[det1] = q2[det2];
            t1[det1] = t2[det2];
            id1 = id2;
            pair = false;
        }
        if (!fPairs)
        {
            // single particle tracking
            fh_xy_tofd->Fill(x[det] * 100., y[det] * 100.);
            fh_tofd_charge->Fill(q[det]);

            if (fB == -1102)
            {
                if (q[det] < 10. || q[det] > 14.)
                    continue;
                // if(x[det] < 0.504 || x[det] > 0.518) continue;
            }
            if (fB == -1292)
            {
                if (q[det] < 8.5 || q[det] > 10.)
                    continue;
                // if(x[det] < 0.504 || x[det] > 0.518) continue;
            }
            if (fB == -1482)
            {
                if (q[det] < 7.5 || q[det] > 9.)
                    continue;
                // if(x[det] < 0.504 || x[det] > 0.518) continue;
            }
            if (fB == -1862)
            {
                if (q[det] < 6.5 || q[det] > 8.)
                    continue;
                // if(x[det] < 0.504 || x[det] > 0.518) continue;
            }
            if (fB == -2052)
            {
                if (q[det] < 6. || q[det] > 7.)
                    continue;
                // if(x[det] < 0.504 || x[det] > 0.518) continue;
            }
            if (fB == -2242)
            {
                if (q[det] < 5.8 || q[det] > 6.6)
                    continue;
                // if(x[det] < 0.504 || x[det] > 0.518) continue;
            }

            fh_tofd_mult_ac->Fill(nHits);
            fh_xy_tofd_ac->Fill(x[det] * 100., y[det] * 100.);
            fh_tofd_charge_ac->Fill(q[det]);
            LOG(DEBUG) << "ToFD"
                       << " x: " << x[det] << " y: " << y[det] << " q: " << q[det] << " t: " << t[det] << endl;
            q[det] = 8.;
            countdet = 0;
            detector[countdet] = det;
            xdet[countdet] = x[det];
            ydet[countdet] = y[det];
            zdet[countdet] = z[det];
            qdet[countdet] = (int)(q[det] + 0.5);
            countdet++;
        }

        if (!pair && fPairs)
            continue;

        for (int i = 0; i < n_det; i++)
        {
            xMax[i] = -1000.;
            yMax[i] = -1000.;
            zMax[i] = -1000.;
            qMax[i] = -1000.;
            tMax[i] = -1000.;
        }

        // loop over fiber 13
        auto detHit13 = fHitItems.at(DET_FI13);
        Int_t nHits13 = detHit13->GetEntriesFast();
        LOG(DEBUG) << "Fi13 hits: " << nHits13 << endl;

        Int_t mult13 = 0;
        for (Int_t ihit13 = 0; ihit13 < nHits13; ihit13++)
        {
            det = fi13;
            R3BBunchedFiberHitData* hit13 = (R3BBunchedFiberHitData*)detHit13->At(ihit13);
            x1[det] = hit13->GetX() / 100.;
            y1[det] = hit13->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit13->GetEloss();
            t1[det] = hit13->GetTime();

            // "Push" two times in the same clock cycle:
            while (t1[tofdr] - t1[det] < 2048. * 5. / 2.)
            {
                t1[det] -= 2048. * 5.;
            }
            while (t1[tofdr] - t1[det] > 2048. * 5. / 2.)
            {
                t1[det] += 2048. * 5.;
            }

            tof = t1[tofdr] - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits13);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofdr] * 100., x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);

            // Cuts on Fi13
            if (fCuts && y1[det] * 100. > 50.)
                continue;
            if (fCuts && y1[det] * 100. < -50.)
                continue;
            if (fCuts && q1[det] < 0)
                continue;
            if (fCuts && (tof < -280 || tof > 0))
                continue;
            if (fB == -1862)
            {
                if (x1[det] < -0.17 || x1[det] > -0.13)
                    continue;
            }
            if (fB == -2052)
            {
                if (x1[det] < -0.01 || x1[det] > 0.02)
                    continue;
            }
            if (fB == -2242)
            {
                if (x1[det] < 0.15 || x1[det] > 0.17)
                    continue;
            }

            mult13++;
            if (mult13 > 10)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult13);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofdr] * 100., x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);

            LOG(DEBUG) << "Fi13: " << ihit13 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                       << " t1: " << t1[det] << endl;

            detector[countdet] = det;
            xdet[countdet] = x1[det];
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = 0;
            countdet++;
        }
        // loop over fiber 11
        auto detHit11 = fHitItems.at(DET_FI11);
        Int_t nHits11 = detHit11->GetEntriesFast();
        LOG(DEBUG) << "Fi11 hits: " << nHits11 << endl;
        Int_t mult11 = 0;
        for (Int_t ihit11 = 0; ihit11 < nHits11; ihit11++)
        {
            det = fi11;
            R3BBunchedFiberHitData* hit11 = (R3BBunchedFiberHitData*)detHit11->At(ihit11);
            x1[det] = hit11->GetX() / 100.;
            y1[det] = hit11->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit11->GetEloss();
            t1[det] = hit11->GetTime();

            // "Push" two times in the same clock cycle:
            while (t1[tofdr] - t1[det] < 2048. * 5. / 2.)
            {
                t1[det] -= 2048. * 5.;
            }
            while (t1[tofdr] - t1[det] > 2048. * 5. / 2.)
            {
                t1[det] += 2048. * 5.;
            }

            tof = t1[tofdr] - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100);
            fh_mult_Fib[det]->Fill(nHits11);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofdr] * 100., x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);

            LOG(DEBUG) << "Fi11 bc: " << ihit11 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                       << " t1: " << t1[det] << endl;

            // Cuts on Fi11
            if (fCuts && y1[det] * 100. > 50.)
                continue;
            if (fCuts && y1[det] * 100. < -50.)
                continue;
            if (fCuts && q1[det] < 0)
                continue;
            if (fCuts && (tof < -350 || tof > 0))
                continue;
            if (fB == -1862)
            {
                if (x1[det] < -0.185 || x1[det] > -0.155)
                    continue;
            }
            if (fB == -2052)
            {
                if (x1[det] < -0.065 || x1[det] > -0.04)
                    continue;
            }
            if (fB == -2242)
            {
                if (x1[det] < 0.06 || x1[det] > 0.075)
                    continue;
            }

            mult11++;
            if (mult11 > 10)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100);
            fh_mult_Fib_ac[det]->Fill(mult11);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofdr] * 100., x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);

            LOG(DEBUG) << "Fi11: " << ihit11 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                       << " t1: " << t1[det] << endl;

            detector[countdet] = det;
            xdet[countdet] = x1[det];
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = 0;
            countdet++;
        }

        // loop over fiber 10
        auto detHit10 = fHitItems.at(DET_FI10);
        Int_t nHits10 = detHit10->GetEntriesFast();
        LOG(DEBUG) << "Fi10 hits: " << nHits10 << endl;
        Int_t mult10 = 0;
        for (Int_t ihit10 = 0; ihit10 < nHits10; ihit10++)
        {
            det = fi10;
            R3BBunchedFiberHitData* hit10 = (R3BBunchedFiberHitData*)detHit10->At(ihit10);
            x1[det] = hit10->GetX() / 100.;
            y1[det] = hit10->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit10->GetEloss();
            t1[det] = hit10->GetTime();

            // "Push" two times in the same clock cycle:
            while (t1[tofdl] - t1[det] < 2048. * 5. / 2.)
            {
                t1[det] -= 2048. * 5.;
            }
            while (t1[tofdl] - t1[det] > 2048. * 5. / 2.)
            {
                t1[det] += 2048. * 5.;
            }

            tof = t1[tofdl] - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits10);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofdl] * 100., x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);

            // Cuts on Fi10
            if (fCuts && q1[det] < 0)
                continue;
            if (fCuts && y1[det] * 100. > 50.)
                continue;
            if (fCuts && y1[det] * 100. < -50.)
                continue;
            if (fCuts && (tof < -180 || tof > -20))
                continue;
            // if(x1[det] > 0.144 || x1[det] < 0.14 ) continue;
            if (fB == -1102)
            {
                if (x1[det] < 0.135 || x1[det] > 0.155)
                    continue;
            }
            if (fB == -1292)
            {
                if (x1[det] < -0.01 || x1[det] > 0.015)
                    continue;
            }
            if (fB == -1482)
            {
                if (x1[det] < -0.155 || x1[det] > -0.115)
                    continue;
            }

            mult10++;
            if (mult10 > 10)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult10);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofdl] * 100., x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);

            LOG(DEBUG) << "Fi10: " << ihit10 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                       << " t1: " << t1[det] << endl;

            detector[countdet] = det;
            xdet[countdet] = x1[det];
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = 0;
            countdet++;
        }
        // loop over fiber 12
        auto detHit12 = fHitItems.at(DET_FI12);
        Int_t nHits12 = detHit12->GetEntriesFast();
        LOG(DEBUG) << "Fi12 hits: " << nHits12 << endl;
        Int_t mult12 = 0;
        for (Int_t ihit12 = 0; ihit12 < nHits12; ihit12++)
        {
            det = fi12;
            R3BBunchedFiberHitData* hit12 = (R3BBunchedFiberHitData*)detHit12->At(ihit12);
            x1[det] = hit12->GetX() / 100.;
            y1[det] = hit12->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit12->GetEloss();
            t1[det] = hit12->GetTime();

            // "Push" two times in the same clock cycle:
            while (t1[tofdl] - t1[det] < 2048. * 5. / 2.)
            {
                t1[det] -= 2048. * 5.;
            }
            while (t1[tofdl] - t1[det] > 2048. * 5. / 2.)
            {
                t1[det] += 2048. * 5.;
            }

            tof = t1[tofdl] - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits12);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofdl] * 100., x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);

            // Cuts on Fi12
            if (fCuts && q1[det] < 0)
                continue;
            if (fCuts && y1[det] * 100. > 50.)
                continue;
            if (fCuts && y1[det] * 100. < -50.)
                continue;
            if (fCuts && (tof < -180 || tof > -20))
                continue;
            // if(x1[det] > 0.050 || x1[det] < 0.048 ) continue;
            if (fB == -1102)
            {
                if (x1[det] < 0.04 || x1[det] > 0.06)
                    continue;
            }
            if (fB == -1292)
            {
                if (x1[det] < -0.07 || x1[det] > -0.046)
                    continue;
            }
            if (fB == -1482)
            {
                if (x1[det] < -0.175 || x1[det] > -0.16)
                    continue;
            }

            mult12++;
            if (mult12 > 10)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult12);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofdl] * 100., x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);

            LOG(DEBUG) << "Fi12: " << ihit12 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                       << " t1: " << t1[det] << endl;

            detector[countdet] = det;
            xdet[countdet] = x1[det];
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = 0;
            countdet++;
        }

        // loop over fiber 3a
        auto detHit3a = fHitItems.at(DET_FI3A);
        Int_t nHits3a = detHit3a->GetEntriesFast();
        LOG(DEBUG) << "Fi3a hits: " << nHits3a << endl;
        Int_t mult3a = 0;
        fh_mult_Fib[fi3a]->Fill(nHits3a);
        for (Int_t ihit3a = 0; ihit3a < nHits3a; ihit3a++)
        {
            det = fi3a;
            R3BBunchedFiberHitData* hit3a = (R3BBunchedFiberHitData*)detHit3a->At(ihit3a);
            x1[det] = hit3a->GetX() / 100.;
            y1[det] = hit3a->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit3a->GetEloss();
            t1[det] = hit3a->GetTime();

            // "Push" the Fib times in the same cycle with ToFD:
            // ToFD coarse counter has 8192 channels with 5 ns each
            // SPMT coarse counter has 2048 channels with 5 ns each

            while (t1[tofdr] - t1[det] < 2048. * 5. / 2.)
            {
                t1[det] -= 2048. * 5.;
            }
            while (t1[tofdr] - t1[det] > 2048. * 5. / 2.)
            {
                t1[det] += 2048. * 5.;
            }

            tof = t1[tofdr] - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits3a);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofdr] * 100., x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);

            // Cuts on Fi3a
            if (fCuts && y1[det] * 100. > 20.)
                continue;
            if (fCuts && y1[det] * 100. < -20.)
                continue;
            if (fCuts && q1[det] < 0)
                continue;
            if (fCuts && (tof < 10 || tof > 50))
                continue;

            mult3a++;
            if (mult3a > 10)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult3a);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofdr] * 100., x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);

            LOG(DEBUG) << "Fi3a " << ihit3a << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                       << " t1: " << tof << endl;

            detector[countdet] = det;
            xdet[countdet] = x1[det];
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = 0;
            countdet++;
        }

        // loop over fiber 3b
        auto detHit3b = fHitItems.at(DET_FI3B);
        Int_t nHits3b = detHit3b->GetEntriesFast();
        LOG(DEBUG) << "Fi3b hits: " << nHits3b << endl;
        Int_t mult3b = 0;
        fh_mult_Fib[fi3b]->Fill(nHits3b);

        for (Int_t ihit3b = 0; ihit3b < nHits3b; ihit3b++)
        {
            det = fi3b;
            R3BBunchedFiberHitData* hit3b = (R3BBunchedFiberHitData*)detHit3b->At(ihit3b);
            x1[det] = hit3b->GetX() / 100.;
            y1[det] = hit3b->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit3b->GetEloss();
            t1[det] = hit3b->GetTime();

            // "Push" the Fib times in the same cycle with ToFD:
            // ToFD coarse counter has 8192 channels with 5 ns each
            // SPMT coarse counter has 2048 channels with 5 ns each

            while (t1[tofdl] - t1[det] < 2048. * 5. / 2.)
            {
                t1[det] -= 2048. * 5.;
            }
            while (t1[tofdl] - t1[det] > 2048. * 5. / 2.)
            {
                t1[det] += 2048. * 5.;
            }

            tof = t1[tofdl] - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits3b);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofdl] * 100., x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);

            // Cuts on Fi3b
            if (fCuts && y1[det] * 100. > 20.)
                continue;
            if (fCuts && y1[det] * 100. < -20.)
                continue;
            if (fCuts && q1[det] < 0)
                continue;
            if (fCuts && (tof < 10 || tof > 50))
                continue;

            mult3b++;
            if (mult3b > 10)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult3b);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofdl] * 100., x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);

            LOG(DEBUG) << "Fi3b " << ihit3b << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                       << " t1: " << tof << endl;

            detector[countdet] = det;
            xdet[countdet] = x1[det];
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = 0;
            countdet++;
        }

        // Plots of correlations of Fiber detectors
        if (xMax[fi11] > -100 && xMax[fi13] > -100)
        {
            fh_Fib13_vs_Fib11->Fill(xMax[fi11] * 100., xMax[fi13] * 100.);
            fh_Fib13_vs_Fib11_dx->Fill(xMax[fi11] * 100., xMax[fi13] * 100. - xMax[fi11] * 100.);
        }
        if (xMax[fi3a] > -100 && xMax[fi11] > -100)
        {
            fh_Fib11_vs_Fib3a->Fill(xMax[fi3a] * 100., xMax[fi11] * 100.);
            fh_Fib11_vs_Fib3a_dx->Fill(xMax[fi3a] * 100., xMax[fi11] * 100. - xMax[fi3a] * 100.);
        }
        if (xMax[fi12] > -100 && xMax[fi10] > -100)
        {
            fh_Fib10_vs_Fib12->Fill(xMax[fi12] * 100., xMax[fi10] * 100.);
            fh_Fib10_vs_Fib12_dx->Fill(xMax[fi12] * 100., xMax[fi10] * 100. - xMax[fi12] * 100.);
        }
        if (xMax[fi3b] > -100 && xMax[fi12] > -100)
        {
            fh_Fib12_vs_Fib3b->Fill(xMax[fi3b] * 100., xMax[fi12] * 100.);
            fh_Fib12_vs_Fib3b_dx->Fill(xMax[fi3b] * 100., xMax[fi12] * 100. - xMax[fi3b] * 100.);
        }

        if (fGhost && (mult3a == 0 || mult3a == 0))
        {
            det = ghost;
            detector[countdet] = det;
            xdet[countdet] = 0.;
            ydet[countdet] = 0.;
            zdet[countdet] = 0.;
            qdet[countdet] = 0;
            countdet++;
        }

    } // end ToFD loop

    fNEvents += 1;
    //    myfile.close();
}
void R3BGlobalAnalysisS454::Output1(Double_t track[12], Double_t chi[6]) {}

void R3BGlobalAnalysisS454::Output2(Double_t track_parameter[12], Double_t chi_single_parameter[6])
{
    // compare

    LOG(DEBUG) << "******************************************" << endl;
    LOG(DEBUG) << "chi_sqingle_parameter " << chi_single_parameter[0] << "  " << chi_single_parameter[1] << endl;
    LOG(DEBUG) << "xyz: " << track_parameter[0] * 100. << "  " << track_parameter[1] * 100. << "  "
               << track_parameter[2] * 100. << endl;
    LOG(DEBUG) << "p: " << track_parameter[3] << "  " << track_parameter[4] << "  " << track_parameter[5] << endl;
    Double_t p_tot = sqrt(track_parameter[3] * track_parameter[3] + track_parameter[4] * track_parameter[4] +
                          track_parameter[5] * track_parameter[5]);
    Double_t thetax_simu = Pxf / Pzf * 1000.;
    Double_t thetay_simu = Pyf / Pzf * 1000.;

    fh_target_xy->Fill(track_parameter[0] * 100., track_parameter[1] * 100.);
    fh_target_px->Fill(track_parameter[3]);
    fh_target_py->Fill(track_parameter[4]);
    fh_target_pz->Fill(track_parameter[5]);
    fh_target_p->Fill(p_tot);
    fh_chi2->Fill(chi_single_parameter[0]);

    fh_thetax_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100., thetax_simu);
    fh_thetax_dpx_abs->Fill((Pxf - track_parameter[3]), thetax_simu);
    fh_thetay_dpy->Fill((Pyf - track_parameter[4]) / Pyf * 100., thetay_simu);

    fh_thetax_dpy->Fill((Pyf - track_parameter[4]) / Pyf * 100., thetax_simu);
    fh_thetay_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100., thetay_simu);
    fh_dpy_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100., (Pyf - track_parameter[4]) / Pyf * 100.);

    fh_x_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100., Xf);
    fh_y_dpy->Fill((Pyf - track_parameter[4]) / Pyf * 100., Yf);
    fh_dx->Fill(Xf - track_parameter[0] * 100.);
    fh_dy->Fill(Yf - track_parameter[1] * 100.);
    fh_dz->Fill(Zf - track_parameter[2] * 100.);
    fh_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100.);
    fh_dpy->Fill((Pyf - track_parameter[4]) / Pyf * 100.);
    fh_dpz->Fill((Pzf - track_parameter[5]) / Pzf * 100.);
    fh_dp->Fill((Pf_tot - p_tot) / Pf_tot * 100.);
}
void R3BGlobalAnalysisS454::FinishEvent()
{
    fTrack->Clear();
    for (Int_t det = 0; det < DET_MAX; det++)
    {
        if (fMappedItems.at(det))
        {
            fMappedItems.at(det)->Clear();
        }
        if (fCalItems.at(det))
        {
            fCalItems.at(det)->Clear();
        }
        if (fHitItems.at(det))
        {
            fHitItems.at(det)->Clear();
        }
    }
}

void R3BGlobalAnalysisS454::FinishTask()
{

    fh_Tpat->Write();
    fh_Trigger->Write();
    fh_Cave_position->Write();
    fh_TOFDOR->Write();
    fh_SEE->Write();
    fh_IC->Write();

    fh_target_xy->Write();
    fh_target_px->Write();
    fh_target_py->Write();
    fh_target_pz->Write();
    fh_target_p->Write();

    fh_px_He->Write();
    fh_py_He->Write();
    fh_pz_He->Write();
    fh_p_He->Write();
    fh_px_C->Write();
    fh_py_C->Write();
    fh_pz_C->Write();
    fh_p_C->Write();

    fh_chi2->Write();

    fh_dx->Write();
    fh_dy->Write();
    fh_dz->Write();
    fh_dpx->Write();
    fh_dpy->Write();
    fh_dpz->Write();
    fh_dp->Write();
    fh_thetax_dpx->Write();
    fh_thetax_dpx_abs->Write();
    fh_thetay_dpy->Write();
    fh_thetax_dpy->Write();
    fh_thetay_dpx->Write();
    fh_x_dpx->Write();
    fh_y_dpy->Write();
    fh_dpy_dpx->Write();

    if (fHitItems.at(DET_TOFD))
    {
        fh_xy_tofd->Write();
        fh_xy_tofd_ac->Write();
        fh_tofd_charge->Write();
        fh_tofd_charge_ac->Write();
        fh_TimePreviousEvent->Write();
        fh_tofd_mult->Write();
        fh_tofd_mult_ac->Write();
        fh_tofd_q2_vs_q1->Write();
        fh_tofd_q2_vs_q1_ac->Write();
    }

    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {
        if (fMappedItems.at(ifibcount + DET_FI_FIRST))
        {
            fh_xy_Fib[ifibcount]->Write();
            fh_xy_Fib_ac[ifibcount]->Write();
            fh_mult_Fib[ifibcount]->Write();
            fh_mult_Fib_ac[ifibcount]->Write();
            fh_ToT_Fib[ifibcount]->Write();
            fh_ToT_Fib_ac[ifibcount]->Write();
            fh_Fib_vs_Events[ifibcount]->Write();
            fh_Fib_vs_Events_ac[ifibcount]->Write();
            fh_Fibs_vs_Tofd[ifibcount]->Write();
            fh_Fibs_vs_Tofd_ac[ifibcount]->Write();
            fh_Fib_ToF[ifibcount]->Write();
            fh_Fib_ToF_ac[ifibcount]->Write();
            fh_ToF_vs_Events[ifibcount]->Write();
            fh_ToF_vs_Events_ac[ifibcount]->Write();
        }
    }

    fh_Fib13_vs_Fib11->Write();
    fh_Fib13_vs_Fib11_dx->Write();
    fh_Fib11_vs_Fib3a->Write();
    fh_Fib11_vs_Fib3a_dx->Write();
    fh_Fib10_vs_Fib12->Write();
    fh_Fib10_vs_Fib12_dx->Write();
    fh_Fib12_vs_Fib3b->Write();
    fh_Fib12_vs_Fib3b_dx->Write();

    fh_theta26->Write();
    fh_theta26_cm->Write();
    fh_theta_16O->Write();
    fh_phi26_cm->Write();
    fh_theta_4He_cm->Write();
    fh_phi_4He_cm->Write();
    fh_theta_12C_cm->Write();
    fh_phi_12C_cm->Write();
    fh_theta_16O_cm->Write();
    fh_phi_16O_cm->Write();
    fh_Erel->Write();
    fh_ErelR->Write();
    fh_ErelL->Write();
    fh_theta26_simu->Write();
    fh_Erel_simu->Write();
    fh_chiy_vs_chix->Write();
    fh_dErel_vs_x->Write();
    fh_dErel_vs_y->Write();
}
ClassImp(R3BGlobalAnalysisS454)
