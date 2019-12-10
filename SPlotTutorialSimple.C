/////////////////////////////////////////////////////////////////////////
//
// SPlot tutorial
// author: Kyle Cranmer
// date Dec. 2008
//
// This tutorial shows an example of using SPlot to unfold two distributions.
// The physics context for the example is that we want to know
// the isolation distribution for real electrons from Z events
// and fake electrons from QCD.  Isolation is our 'control' variable
// To unfold them, we need a model for an uncorrelated variable that
// discriminates between Z and QCD.  To do this, we use the invariant
// mass of two electrons.  We model the Z with a Gaussian and the QCD
// with a falling exponential.
//
// Note, since we don't have real data in this tutorial, we need to generate
// toy data.  To do that we need a model for the isolation variable for
// both Z and QCD.  This is only used to generate the toy data, and would
// not be needed if we had real data.
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"

// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;


// see below for implementation
void AddModel(RooWorkspace*);
void AddData(RooWorkspace*);
void DoSPlot(RooWorkspace*);
void MakePlots(RooWorkspace*);

void SPlotTutorialSimple()
{

  // Create a new workspace to manage the project.
  RooWorkspace* wspace = new RooWorkspace("myWS");

  AddModel(wspace);

  //AddData(wspace);

  DoSPlot(wspace);

  MakePlots(wspace);

  // cleanup
  delete wspace;

}


//____________________________________
void AddModel(RooWorkspace* ws){

  // set range of observable
  Double_t lowRange = 70, highRange = 110;

  // make a RooRealVar for the observables
  RooRealVar invMass("invMass", "M_{inv}", lowRange, highRange,"GeV");
  RooRealVar Mu_pT("Mu_pT", "Muon Pt", 0., 200., "GeV");
  //define dataset as args 
  RooArgSet args(invMass,Mu_pT);  
  RooDataSet *data =  new RooDataSet("data","data",args,0);

  /////////////////////////////////////////////
  // make 2-d model for Z including the invariant mass
  // distribution  and an muon Pt distribution from MC

  std::cout << "make z model" << std::endl;
  
  // mass model for Z
  RooRealVar mZ("mZ", "Z Mass", 91.2,"GeV");
  RooRealVar widthZ("widthZ", "Width of Voigtian", 2.5,"GeV");
  RooRealVar sigmaZVoig("sigmaZVoig", "Sigma of Voigtian", 1.173,"GeV");
  RooVoigtian mZModelVoig("mZModelVoig", "Z Model Voigtian", invMass, mZ, widthZ, sigmaZVoig);
  //add gaussian model
  RooRealVar sigmaZGauss("sigmaZGauss", "Sigma of Gaussian", 5.48,"GeV");
  RooGaussian mZModelGauss("mZModel", "Z 2 Model", invMass, mZ, sigmaZGauss);
  //add these two
  RooRealVar frac("frac", "fraction", 0.9017);
  RooAddPdf  mZModel("mZmodel","Gauss + Voig  Z models", RooArgList(mZModelVoig, mZModelGauss),frac);
    
  
  //MuonPt of Z from MC assuming it goes landau function
  RooRealVar m_pTMu("m_pTMu", "mean Mu Pt ", 30);
  RooRealVar sigmaMuPt("sigmaMuPt", "Sigma of Landau", 10,5,20,"GeV");
  RooLandau MuPtModel("MuPtModel", "Muon Pt Model", Mu_pT, m_pTMu, sigmaMuPt);
  
  // make the combined Z model (mass and muon Pt Product)
  RooProdPdf zModel("zModel", "4-d model for Z", RooArgSet(mZModel, MuPtModel));

  //////////////////////////////////////////////
  // make background model

  std::cout << "make background model" << std::endl;
  
  //Z Background assumed to be exponential 
  RooRealVar qcdMassDecayConst("qcdMassDecayConst", "Decay const for QCD mass spectrum", -0.05, -10, 10,"1/GeV");
  RooExponential qcdMassModel("qcdMassModel", "qcd Mass Model", invMass, qcdMassDecayConst);
  
  //Muon Pt Background is also assumed to me landau with lower mean  
  RooRealVar m_pTMuBak("m_pTMuBak", "mean Mu Pt Background", 5);
  RooRealVar sigmaMuPtBak("sigmaMuPtBak", "Sigma of Landau Background", 10,5,20,"GeV");
  RooLandau MuPtModelBak("MuPtModelBak", "Muon Pt Model Background", Mu_pT, m_pTMuBak, sigmaMuPt);
  
  // make the 2-d model of background of muon Pt and Z mass
  RooProdPdf qcdModel("qcdModel", "2-d model for QCD", RooArgSet(qcdMassModel, MuPtModelBak));

  //////////////////////////////////////////////
  // combined model

  // These variables represent the number of Z or QCD events
  // They will be fitted.
  RooRealVar zYield("zYield","fitted yield for Z",1000 ,0.,5000) ;
  RooRealVar qcdYield("qcdYield","fitted yield for QCD", 1000 ,0.,5000) ;

  // now make the combined model
  std::cout << "make full model" << std::endl;
  RooAddPdf model("model","z+qcd background models",
                  RooArgList(zModel, qcdModel),
                  RooArgList(zYield,qcdYield));


  // interesting for debugging and visualizing the model
  model.graphVizTree("fullModel.dot");

  std::cout << "import model" << std::endl;

  ws->import(model);
  std::cout << "finished importing model" << std::endl;
  std::cout << "Reading data from *.dat filel" << std::endl;
  data = RooDataSet::read("DataZMass_MuPT.dat",args);//check for random J from HZZ* MC             
  std::cout << "Now import dataset to workspace" << std::endl;
    
  // import data into workspace
  ws->import(*data, Rename("data"));
}

//____________________________________
void DoSPlot(RooWorkspace* ws){
  std::cout << "Calculate sWeights" << std::endl;
  
  // get what we need out of the workspace to do the fit
  RooAbsPdf* model = ws->pdf("model");
  RooRealVar* zYield = ws->var("zYield");
  RooRealVar* qcdYield = ws->var("qcdYield");
  RooDataSet* data = (RooDataSet*) ws->data("data");
  
  // fit the model to the data.
  model->fitTo(*data, Extended() );

  // The sPlot technique requires that we fix the parameters
  // of the model that are not yields after doing the fit.
  RooRealVar* sigmaZGauss = ws->var("sigmaZGauss");
  RooRealVar* sigmaZVoig = ws->var("sigmaZVoig");
  RooRealVar* widthZ = ws->var("widthZ");
  RooRealVar* qcdMassDecayConst = ws->var("qcdMassDecayConst");
  sigmaZGauss->setConstant();
  sigmaZVoig->setConstant();
  widthZ->setConstant();
  qcdMassDecayConst->setConstant();
    
  RooMsgService::instance().setSilentMode(true);
  
  // Now we use the SPlot class to add SWeights to our data set
  // based on our model and our yield variables
  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
                                               *data, model, RooArgList(*zYield,*qcdYield) );


  // Check that our weights have the desired properties
  
  std::cout << "Check SWeights:" << std::endl;

  
  std::cout << std::endl <<  "Yield of Z is "
            << zYield->getVal() << ".  From sWeights it is "
            << sData->GetYieldFromSWeight("zYield") << std::endl;
  
  
  std::cout << "Yield of QCD is "
            << qcdYield->getVal() << ".  From sWeights it is "
            << sData->GetYieldFromSWeight("qcdYield") << std::endl
            << std::endl;
  
  for(Int_t i=0; i < 10; i++)
    {
      std::cout << "z Weight   " << sData->GetSWeight(i,"zYield")
                << "   qcd Weight   " << sData->GetSWeight(i,"qcdYield")
                << "  Total Weight   " << sData->GetSumOfEventSWeight(i)
                << std::endl;
    }
  
  std::cout << std::endl;
  
  // import this new dataset with sWeights
  std::cout << "import new dataset with sWeights" << std::endl;
  ws->import(*data, Rename("dataWithSWeights"));
}

void MakePlots(RooWorkspace* ws){
  
  // Here we make plots of the discriminating variable (invMass) after the fit
  // and of the control variable (isolation) after unfolding with sPlot.
  std::cout << "make plots" << std::endl;

  // make our canvas
  TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
  cdata->Divide(1,3);
  
  // get what we need out of the workspace
  RooAbsPdf* model = ws->pdf("model");
  RooAbsPdf* zModel = ws->pdf("zModel");
  RooAbsPdf* qcdModel = ws->pdf("qcdModel");

  RooRealVar* Mu_pT = ws->var("Mu_pT");
  RooRealVar* invMass = ws->var("invMass");

  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

  // this shouldn't be necessary, need to fix something with workspace
  // do this to set parameters back to their fitted values.
  model->fitTo(*data, Extended() );

  //plot invMass for data with full model and individual componenets overlayed
  //  TCanvas* cdata = new TCanvas();
  cdata->cd(1);
  RooPlot* frame = invMass->frame(25) ;
  
  data->plotOn(frame ) ;
  model->plotOn(frame) ;
  model->plotOn(frame,Components(*zModel),LineStyle(kDashed), LineColor(kRed)) ;
  model->plotOn(frame,Components(*qcdModel),LineStyle(kDashed),LineColor(kGreen)) ;

  frame->SetTitle("Fit of model to discriminating variable");
  frame->Draw() ;


  // Now use the sWeights to show isolation distribution for Z and QCD.
  // The SPlot class can make this easier, but here we demonstrait in more
  // detail how the sWeights are used.  The SPlot class should make this
  // very easy and needs some more development.

  // Plot isolation for Z component.
  // Do this by plotting all events weighted by the sWeight for the Z component.
  // The SPlot class adds a new variable that has the name of the corresponding
  // yield + "_sw".
  cdata->cd(2);

  // create weightfed data set
  RooDataSet * dataw_z = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"zYield_sw") ;

  RooPlot* frame2 = Mu_pT->frame(25) ;
  dataw_z->plotOn(frame2, DataError(RooAbsData::SumW2) ) ;

  frame2->SetTitle("Muon Pt for Z Events");
  frame2->Draw() ;

  // Plot isolation for QCD component.
  // Eg. plot all events weighted by the sWeight for the QCD component.
  // The SPlot class adds a new variable that has the name of the corresponding
  // yield + "_sw".
  cdata->cd(3);
  RooDataSet * dataw_qcd = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"qcdYield_sw") ;
  RooPlot* frame3 = Mu_pT->frame(25) ;
  dataw_qcd->plotOn(frame3,DataError(RooAbsData::SumW2) ) ;

  frame3->SetTitle("Muon Pt For Background Events in Z");
  frame3->Draw() ;

  //  cdata->SaveAs("SPlot.gif");

}
