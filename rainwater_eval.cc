////////////////////////////////////////////////////////////////////////////
// Programm for evaluating ambient air data from the Picarro #insertModel //
// author: Sebastian Stezura // E-Mail: sebastian_stezura@gmx.net         //
// GitHub: CplusplusB                                                     //
// version: 1 // date: 03.01.2022                                         //
////////////////////////////////////////////////////////////////////////////

//read .csv files from month, event based rainsamples and flask
//read meteo data and correlated them to the rain samples
//plot some Graphs

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compile command:                                                                                                          //
// g++-10 rainwater_eval.cc tinyfiledialogs.c -ltbb `root-config --cflags --glibs --ldflags` -lMinuit -o ./Rainwater_eval.o  //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
// fancy cross-plattform file dialog, see also:                                 //
// /         \ hello.c v3.8.3 [Nov 1, 2020] zlib licence                        //
// |tiny file| Hello World file created [November 9, 2014]                      //
// | dialogs | Copyright (c) 2014 - 2020 Guillaume Vareille http://ysengrin.com //
// \____  ___/ http://tinyfiledialogs.sourceforge.net                           //
//      \|     git clone http://git.code.sf.net/p/tinyfiledialogs/code tinyfd   //
//////////////////////////////////////////////////////////////////////////////////
#include "tinyfiledialogs.h"

////////////////////
// C/C++ includes //
////////////////////
#include <iostream> //for Input/Output functions
#include <string> //for using strings
#include <vector> //for using vectors
#include <filesystem> //using c++-17 filesystem
#include <fstream> //for reading and writing to files
#include <sstream> //for reading files
#include <algorithm> //for different C/C++ functions
#include <ctime>  //for time
#include <execution> //for parallel stuff
#include <pthread.h> //multithreading
#include <thread> //multithreading

/////////////////////////////
// Root includes see also: //
// https://root.cern/      //
/////////////////////////////
#include <TMath.h>
#include <TApplication.h> //showing GUI
#include <TSystem.h> //System functions
#include <TCanvas.h> //Canvas for Grahps
#include <TFrame.h> //Frame for Graphs
#include <TGraph.h> //drawing Graphs
#include <TMultiGraph.h> //many Graphs in single Canvas
#include <TGraphErrors.h> //graphs with Error bars
#include <TF1.h> //for functions
#include <TH1F.h> //for Histogram with float
#include <TLatex.h> //for legends etc in Graphs
#include <TLine.h> //for drawing lines
#include <TMinuit.h> //for fittingTvirtualFitter
#include <TVirtualFitter.h> //for fitter
#include <TLinearFitter.h>
#include <TDatime.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TROOT.h>

using namespace std;
namespace fs = std::filesystem;


//class for data storage
class Data
{
public:
    vector<double> timed, timed_begin, timed_end, time_begin, time_end, O18, O18_sd, H2, H2_sd, H2O_m, flask_mass, month, D_excess, time_mean,time_mean_err, D_excess_sd, H2O_m_sd;
    vector<string> event_no, year, month_str;
    vector<int> flask_no;
    vector<double> windvel, contemp, rh1, rh2, grad, apress, o3g1, o3g3, no, ventemp, winddir, prec;
    vector<string> month_all = {"01","02","03","04","05","06","07","08","09","10","11","12"};
    Data(){};
    ~Data(){};
};

//getting Event Data from file
void getEvent(string datapath, Data &data, Data &data_flask)
{
    string year_r, date_r, time_r, flask_no_r, event_no_r, O18_r, O18_sd_r, H2_r, H2_sd_r, H2O_m_r;
    double time_beginr, time_endr, O18r, O18_sdr, H2r, H2_sdr, H2O_mr;
    int flask_nor;
    string date_begin, date_end;
    double month_mean;

    //Read file
    ifstream inFile(datapath);
    cout << "Reading file at " << datapath << " ..." << endl;

    if (inFile.is_open())
	{
        string line;

        while (getline(inFile, line))
		{
            remove(line.begin(), line.end(), ' ');
            stringstream stst(line);
            getline(stst,year_r,',');
            getline(stst,date_r,',');
            getline(stst,time_r,',');
            getline(stst,flask_no_r,',');
            getline(stst,event_no_r,',');
            getline(stst,O18_r,',');
            getline(stst,O18_sd_r,',');
            getline(stst,H2_r,',');
            getline(stst,H2_sd_r,',');
            getline(stst,H2O_m_r,',');

            //format time format YYYYMMHH00
            date_begin = year_r + date_r.substr(3,2) + date_r.substr(0,2) + time_r.substr(0,2) + "00";
            date_end = year_r + date_r.substr(9,2) + date_r.substr(6,2) + time_r.substr(3,2) + "00";

            try
            {
                time_beginr = stod(date_begin);
                time_endr = stod(date_end);
                flask_nor = stoi(flask_no_r);
                O18r = stod(O18_r);
                O18_sdr = stod(O18_sd_r);
                H2r = stod(H2_r);
                H2_sdr = stod(H2_sd_r);
                H2O_mr = stod(H2O_m_r);
            }
            catch (...)
            {
                cout << "Reading error: " << event_no_r << " at " << date_begin << "||" << date_end << endl;
            }
            month_mean = ((stod(date_end.substr(4,4)) - stod(date_begin.substr(4,4)))/2. + stod(date_begin.substr(4,4)))/100.;
            data.year.push_back(year_r);
            data.time_begin.push_back(time_beginr);
            data.time_end.push_back(time_endr);
            data.time_mean.push_back((time_beginr+time_endr)/2.);
            data.timed_begin.push_back(stod(date_r.substr(3,2)+date_r.substr(0,2)+time_r.substr(0,2)));
            data.timed_end.push_back(stod(date_r.substr(9,2)+date_r.substr(6,2)+time_r.substr(3,2)));
            data.flask_no.push_back(flask_nor);
            data.event_no.push_back(event_no_r);
            data.O18.push_back(O18r);
            data.O18_sd.push_back(O18_sdr);
            data.H2.push_back(H2r);
            data.H2_sd.push_back(H2_sdr);
            data.H2O_m.push_back(H2O_mr - data_flask.flask_mass[flask_nor-1]);
            data.H2O_m_sd.push_back(0.1);
            data.D_excess.push_back(H2r - O18r*8.);
            data.D_excess_sd.push_back(TMath::Sqrt((H2_sdr)*(H2_sdr)+(8.*O18_sdr)*(8.*O18_sdr)));
            data.month.push_back(month_mean);
            cout << date_begin << "||" << date_end << "=" << month_mean << endl;
        };
	};
    inFile.close();

};

//getting Month Data from file
void getMonth(string datapath, Data &data)
{
    string year_r, date_r, time_r, flask_no_r, event_no_r, O18_r, O18_sd_r, H2_r, H2_sd_r, H2O_m_no_r, H2O_m_r, month_r;
    double time_beginr, time_endr, O18r, O18_sdr, H2r, H2_sdr, H2O_mr, monthr;
    string date_begin, date_end;

    //Read file
    ifstream inFile(datapath);
    cout << "Reading file at " << datapath << " ..." << endl;

    if (inFile.is_open())
	{
        string line;
        int I = 1;
        while (getline(inFile, line))
		{
            remove(line.begin(), line.end(), ' ');
            stringstream stst(line);
            getline(stst,year_r,',');
            getline(stst,date_r,',');
            getline(stst,time_r,',');
            getline(stst,event_no_r,',');
            getline(stst,O18_r,',');
            getline(stst,O18_sd_r,',');
            getline(stst,H2_r,',');
            getline(stst,H2_sd_r,',');
            getline(stst,flask_no_r,',');
            getline(stst,H2O_m_no_r,',');
            getline(stst,H2O_m_r,',');
            getline(stst,month_r,',');


            //format time format
            date_begin = year_r + date_r.substr(0,2) + date_r.substr(3,2) + time_r.substr(0,2) + "00";
            date_end = year_r + date_r.substr(6,2) + date_r.substr(9,2) + time_r.substr(3,2) + "00";


            try
            {
                monthr = stod(month_r);
                time_beginr = stod(date_begin);
                time_endr = stod(date_end);
                O18r = stod(O18_r);
                O18_sdr = stod(O18_sd_r);
                H2r = stod(H2_r);
                H2_sdr = stod(H2_sd_r);
                H2O_mr = stod(H2O_m_r);
            }
            catch (...)
            {
                cout << "Reading error: " << event_no_r << " at " << date_begin << "||" << date_end << endl;
            }

            data.year.push_back(year_r);
            data.time_begin.push_back(time_beginr);
            data.time_end.push_back(time_endr);
            data.time_mean.push_back(I);
            data.O18.push_back(O18r);
            data.O18_sd.push_back(O18_sdr);
            data.H2.push_back(H2r);
            data.H2_sd.push_back(H2_sdr);
            data.H2O_m.push_back(H2O_mr);
            data.month.push_back(monthr);
            data.month_str.push_back(month_r);
            data.event_no.push_back(event_no_r);
            data.D_excess.push_back(H2r - O18r*8.);
            data.time_mean_err.push_back(1.);
            I++;

        };
	};
    inFile.close();

};

//getting Month Data from file
void getFlask(string datapath, Data &data)
{
    string flask_no_r, flask_mass_r;
    double flask_massr;
    int flask_nor;

    //Read file
    ifstream inFile(datapath);
    cout << "Reading file at " << datapath << " ..." << endl;

    if (inFile.is_open())
	{
        string line;
        int i = 0;

        while (getline(inFile, line))
		{
            if(i == 0){i++; continue;};
            remove(line.begin(), line.end(), ' ');
            stringstream stst(line);
            getline(stst,flask_no_r,',');
            getline(stst,flask_mass_r,',');

            try
            {
                flask_massr = stod(flask_mass_r);
                flask_nor = stoi(flask_no_r);
            }
            catch (...)
            {
                cout << "Reading error: " << flask_mass_r << " at " << flask_no_r << endl;
            }

            data.flask_no.push_back(flask_nor);
            data.flask_mass.push_back(flask_massr);
        };
	};
    inFile.close();

};

//getting Meteo Data from file
void getMeteo(string datapath, Data &data)
{
    string interval_r, windvel_r, contemp_r, rh1_r, rh2_r, grad_r, apress_r, o3g1_r, o3g3_r, no_r, ventemp_r, winddir_r, prec_r;
    //loop over alle files
    ifstream inFile (datapath);
    string name = datapath.substr(datapath.size()-27,27);
    cout << "Reading file " << name << " ..." << endl;

    // read signal values from file
	if (inFile.is_open())
	{
        int i = 0;
        string line;

        while (getline(inFile, line))
		{
            if (i < 1) {i++;continue;};
            remove(line.begin(), line.end(), ' ');
            stringstream stst(line);
            getline(stst,interval_r,';');
            getline(stst,windvel_r,';');
            getline(stst,contemp_r,';');
            getline(stst,rh1_r,';');
            getline(stst,rh2_r,';');
            getline(stst,grad_r,';');
            getline(stst,apress_r,';');
            getline(stst,o3g1_r,';');
            getline(stst,o3g3_r,';');
            getline(stst,no_r,';');
            getline(stst,ventemp_r,';');
            getline(stst,winddir_r,';');
            getline(stst,prec_r,';');
            try
            {
                remove(interval_r.begin(), interval_r.end(), '-');
                remove(interval_r.begin(), interval_r.end(), ':');

                data.timed.push_back(stod(interval_r));
                data.windvel.push_back(stod(windvel_r));
                data.contemp.push_back(stod(contemp_r));
                data.rh1.push_back(stod(rh1_r));
                data.rh2.push_back(stod(rh2_r));
                data.grad.push_back(stod(grad_r));
                data.apress.push_back(stod(apress_r));
                // data.o3g1.push_back(stod(o3g1_r));
                // data.o3g3.push_back(stod(o3g3_r));
                data.no.push_back(stod(no_r));
                data.ventemp.push_back(stod(ventemp_r));
                data.winddir.push_back(stod(winddir_r));
                data.prec.push_back(stod(prec_r));
                data.month_str.push_back(interval_r.substr(4,2));
                data.year.push_back(interval_r.substr(0,4));
                data.timed_begin.push_back(stod(interval_r.substr(4,6)));
                if (i == 20 || i == 40){cout << name << ": " << data.timed.back() << "||" << interval_r << "||" << windvel_r << "||" << contemp_r << "||" << rh1_r << "||" << rh2_r << "||" << grad_r << "||" << apress_r << "||" << o3g1_r << "||" << o3g3_r << "||" << no_r << "||" << ventemp_r << "||" << winddir_r << "||" << prec_r << endl;};
                if (i == 20){cout << "Year " << interval_r.substr(0,4) << endl;};
            }
            catch (...)
            {
                cout << "Problem reading Meteo at " << interval_r << endl;
            }
            i++;
        };
	};
    inFile.close();
};

// draw Graph
void drawEventGraph(Data &data, string evalpath)
{
    string name_file_all = evalpath + "/End/" + "rain_event_year.png";
    string name_file_lmwl = evalpath + "/End/" + "rain_event_LMWL.png";

    string O18title = "{}^{18}O Graph;Month of Year;#delta{}^{18}O";
    string H2title = "{}^{2}H Graph;Month of Year;#delta{}^{2}H";
    string Dextitle = "D-excess Graph;Month of Year;D_{excess}";
    string H2Otitle = "Precipitation Graph;Month of Year;H_{2}O[g]";

    vector<vector<double>> xO18,xH2,xDex,xH2O;
    vector<vector<double>> yO18,yH2,yDex,yH2O;

    vector<string> year;

    vector<int> color;

    for(int i = 0; i < 4; i++)
    {
        xO18.push_back(vector<double>());
        xH2.push_back(vector<double>());
        xDex.push_back(vector<double>());
        xH2O.push_back(vector<double>());
        yO18.push_back(vector<double>());
        yH2.push_back(vector<double>());
        yDex.push_back(vector<double>());
        yH2O.push_back(vector<double>());
        if(i == 0){year.push_back("2018");color.push_back(1);}; //black
        if(i == 1){year.push_back("2019");color.push_back(416);}; //gren
        if(i == 2){year.push_back("2020");color.push_back(600);}; //blue
        if(i == 3){year.push_back("2021");color.push_back(800);};  //orange
    };
    int width = 8000;
    int height = 3000;


    TCanvas *cAll = new TCanvas("All","All",0,0,width,height*4);
    cAll->Divide(1,4,0,0);
    cAll->GetFrame()->SetBorderSize(12);
    cAll->SetGrid();
    TMultiGraph *mgO18 = new TMultiGraph();
    mgO18->SetTitle(O18title.c_str());
    TMultiGraph *mgH2 = new TMultiGraph();
    mgH2->SetTitle(H2title.c_str());
    TMultiGraph *mgDex = new TMultiGraph();
    mgDex->SetTitle(Dextitle.c_str());
    TMultiGraph *mgH2O = new TMultiGraph();
    mgH2O->SetTitle(H2Otitle.c_str());
    vector<TGraph*> grO18, grH2, grDex, grH2O;
    for(int i = 0; i < data.O18.size();i++)
    {
        if(data.year[i] == "2018"){xO18[0].push_back(data.month[i]);yO18[0].push_back(data.O18[i]);};
        if(data.year[i] == "2019"){xO18[1].push_back(data.month[i]);yO18[1].push_back(data.O18[i]);};
        if(data.year[i] == "2020"){xO18[2].push_back(data.month[i]);yO18[2].push_back(data.O18[i]);};
        if(data.year[i] == "2021"){xO18[3].push_back(data.month[i]);yO18[3].push_back(data.O18[i]);};
        if(data.year[i] == "2018"){xH2[0].push_back(data.month[i]);yH2[0].push_back(data.H2[i]);};
        if(data.year[i] == "2019"){xH2[1].push_back(data.month[i]);yH2[1].push_back(data.H2[i]);};
        if(data.year[i] == "2020"){xH2[2].push_back(data.month[i]);yH2[2].push_back(data.H2[i]);};
        if(data.year[i] == "2021"){xH2[3].push_back(data.month[i]);yH2[3].push_back(data.H2[i]);};
        if(data.year[i] == "2018"){xDex[0].push_back(data.month[i]);yDex[0].push_back(data.D_excess[i]);};
        if(data.year[i] == "2019"){xDex[1].push_back(data.month[i]);yDex[1].push_back(data.D_excess[i]);};
        if(data.year[i] == "2020"){xDex[2].push_back(data.month[i]);yDex[2].push_back(data.D_excess[i]);};
        if(data.year[i] == "2021"){xDex[3].push_back(data.month[i]);yDex[3].push_back(data.D_excess[i]);};
        if(data.year[i] == "2018"){xH2O[0].push_back(data.month[i]);yH2O[0].push_back(data.H2O_m[i]);};
        if(data.year[i] == "2019"){xH2O[1].push_back(data.month[i]);yH2O[1].push_back(data.H2O_m[i]);};
        if(data.year[i] == "2020"){xH2O[2].push_back(data.month[i]);yH2O[2].push_back(data.H2O_m[i]);};
        if(data.year[i] == "2021"){xH2O[3].push_back(data.month[i]);yH2O[3].push_back(data.H2O_m[i]);};
    };
    for(int i = 0; i < 4; i++)
    {
        grO18.push_back(new TGraph(xO18[i].size(),xO18[i].data(),yO18[i].data()));
        grH2.push_back(new TGraph(xH2[i].size(),xH2[i].data(),yH2[i].data()));
        grDex.push_back(new TGraph(xDex[i].size(),xDex[i].data(),yDex[i].data()));
        grH2O.push_back(new TGraph(xH2O[i].size(),xH2O[i].data(),yH2O[i].data()));

        grO18[i]->SetTitle(year[i].c_str());grH2[i]->SetTitle(year[i].c_str());grDex[i]->SetTitle(year[i].c_str());grH2O[i]->SetTitle(year[i].c_str());
        grO18[i]->SetMarkerStyle(43);grH2[i]->SetMarkerStyle(43);grDex[i]->SetMarkerStyle(43);grH2O[i]->SetMarkerStyle(43);
        grO18[i]->SetMarkerColor(kBlack);grH2[i]->SetMarkerColor(kGreen);grDex[i]->SetMarkerColor(kOrange);grH2O[i]->SetMarkerColor(kBlue);
        grO18[i]->SetMarkerSize(5);grH2[i]->SetMarkerSize(5);grDex[i]->SetMarkerSize(5);grH2O[i]->SetMarkerSize(5);
        grO18[i]->SetLineWidth(8);grH2[i]->SetLineWidth(8);grDex[i]->SetLineWidth(8);grH2O[i]->SetLineWidth(8);
        grO18[i]->SetLineColor(color[i]);grH2[i]->SetLineColor(color[i]);grDex[i]->SetLineColor(color[i]);grH2O[i]->SetLineColor(color[i]);
        mgO18->Add(grO18[i]);mgH2->Add(grH2[i]);mgDex->Add(grDex[i]);mgH2O->Add(grH2O[i]);
    };


    cAll->cd(1);
    cout << "Drawing Month O18 ... " << endl;
    mgO18->Draw("ALP");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(2);
    cout << "Drawing Month H2 ... " << endl;
    mgH2->Draw("ALP");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(3);
    cout << "Drawing Month Dex ... " << endl;
    mgDex->Draw("ALP");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(4);
    cout << "Drawing Month H2O ... " << endl;
    mgH2O->Draw("ALP");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->Print(name_file_all.c_str());
    cAll->Close();


    vector<TGraph*> grLMWL;
    vector<TF1*> linFit;
    TCanvas *cLMWL = new TCanvas("LMWL","LMWL",0,0,2000,2000);
    cLMWL->Divide(2,2);
    vector<string> LMWLtitle;
    for(int i = 0; i < 4; i++)
    {
        cLMWL->cd(i+1);
        LMWLtitle.push_back("LMWL " + year[i] + ";#delta{}^{18}O;#delta{}^{2}H");
        linFit.push_back(new TF1(year[i].c_str(),"pol1",-100,100));
        grLMWL.push_back(new TGraph(yO18[i].size(),yO18[i].data(),yH2[i].data()));
        grLMWL[i]->SetTitle(LMWLtitle[i].c_str());
        grLMWL[i]->SetMarkerStyle(43);
        grLMWL[i]->SetMarkerColor(kBlue);
        grLMWL[i]->SetMarkerSize(3);
        grLMWL[i]->Draw("AP");
        linFit[i]->SetParameters(8.,10.);
        grLMWL[i]->Fit(linFit[i]);
        linFit[i] = grLMWL[i]->GetFunction(year[i].c_str());
        linFit[i]->SetLineColor(kRed);
        linFit[i]->SetLineStyle(4);
        linFit[i]->SetLineWidth(6);
        cLMWL->Update();
        gPad->Update();
        cLMWL->Update();
    };
    cLMWL->Print(name_file_lmwl.c_str());
    cLMWL->Close();
};

//draw Graph Event Meteo
void drawEventMeteoGraph(Data &data, vector<Data> &data_meteo, string evalpath)
{
    vector<double> temp, rh, windvel, winddir, grad;
    vector<string> month, year;

    vector<vector<double>> xTemp,xRh,xWindvel, xWinddir, xGrad;
    vector<vector<double>> yO18,yH2,yDex;
    vector<int> color;

    //preparing vectors
    for(int i = 0; i < 4; i++)
    {
        xTemp.push_back(vector<double>());
        xRh.push_back(vector<double>());
        xWindvel.push_back(vector<double>());
        xWinddir.push_back(vector<double>());
        xGrad.push_back(vector<double>());
        yO18.push_back(vector<double>());
        yH2.push_back(vector<double>());
        yDex.push_back(vector<double>());
        if(i == 0){year.push_back("2018");color.push_back(1);};
        if(i == 1){year.push_back("2019");color.push_back(416);};
        if(i == 2){year.push_back("2020");color.push_back(600);};
        if(i == 3){year.push_back("2021");color.push_back(800);};
    };

    for(int i = 0; i < data.timed_begin.size(); i++)
    {
        temp.clear();
        rh.clear();
        windvel.clear();
        winddir.clear();
        grad.clear();

        if(data.year[i] == "2018")
        {
            cout << "Year read and compare " << data.year[i] << endl;
            for(int j = 0; j < data_meteo[0].timed.size(); j++)
            {
                if(i == 0 && j == 0)
                {
                    cout << "Time compare " << data.timed_begin[i] << "||" << data_meteo[0].timed_begin[0] << "||" << data.timed_end[i] << endl;
                }
                if(data_meteo[0].timed_begin[j] >= data.timed_begin[i] &&  data_meteo[0].timed_begin[j] < data.timed_end[i])
                {
                    temp.push_back(data_meteo[0].ventemp[j]);
                    rh.push_back(data_meteo[0].rh2[j]);
                    windvel.push_back(data_meteo[0].windvel[j]);
                    winddir.push_back(data_meteo[0].winddir[j]);
                    grad.push_back(data_meteo[0].grad[j]);
                    if(i == 0)
                    {
                        cout << "Data meteo: " << data_meteo[0].ventemp[j] << "||" << data_meteo[0].rh2[j] << "||" << data_meteo[0].windvel[j] << "||" << data_meteo[0].winddir[j] << "||" << data_meteo[0].grad[j] << endl;
                    };
                };
            };
            cout << "Averaging in " << data.year[i] << "||" << TMath::Mean(temp.begin(),temp.end()) << endl;
            xTemp[0].push_back(TMath::Mean(temp.begin(),temp.end()));
            xRh[0].push_back(TMath::Mean(rh.begin(),rh.end()));
            xWindvel[0].push_back(TMath::Mean(windvel.begin(),windvel.end()));
            xWinddir[0].push_back(TMath::Mean(winddir.begin(),winddir.end()));
            xGrad[0].push_back(TMath::Mean(grad.begin(),grad.end()));
        };
        if(data.year[i] == "2019")
        {
            cout << "Year read and compare " << data.year[i] << endl;
            for(int j = 0; j < data_meteo[1].timed.size(); j++)
            {
                if(data_meteo[1].timed_begin[j] >= data.timed_begin[i] &&  data_meteo[1].timed_begin[j] < data.timed_end[i])
                {
                    temp.push_back(data_meteo[1].ventemp[j]);
                    rh.push_back(data_meteo[1].rh2[j]);
                    windvel.push_back(data_meteo[1].windvel[j]);
                    winddir.push_back(data_meteo[1].winddir[j]);
                    grad.push_back(data_meteo[1].grad[j]);
                };
            };
            cout << "Averaging in " << data.year[i] << "||" << TMath::Mean(temp.begin(),temp.end()) << endl;
            xTemp[1].push_back(TMath::Mean(temp.begin(),temp.end()));
            xRh[1].push_back(TMath::Mean(rh.begin(),rh.end()));
            xWindvel[1].push_back(TMath::Mean(windvel.begin(),windvel.end()));
            xWinddir[1].push_back(TMath::Mean(winddir.begin(),winddir.end()));
            xGrad[1].push_back(TMath::Mean(grad.begin(),grad.end()));
        };
        if(data.year[i] == "2020")
        {
            cout << "Year read and compare " << data.year[i] << endl;
            for(int j = 0; j < data_meteo[2].timed.size(); j++)
            {
                if(data_meteo[2].timed_begin[j] >= data.timed_begin[i] &&  data_meteo[2].timed_begin[j] < data.timed_end[i])
                {
                    temp.push_back(data_meteo[2].ventemp[j]);
                    rh.push_back(data_meteo[2].rh2[j]);
                    windvel.push_back(data_meteo[2].windvel[j]);
                    winddir.push_back(data_meteo[2].winddir[j]);
                    grad.push_back(data_meteo[2].grad[j]);
                };
            };
            cout << "Averaging in " << data.year[i] << "||" << TMath::Mean(temp.begin(),temp.end()) << endl;
            xTemp[2].push_back(TMath::Mean(temp.begin(),temp.end()));
            xRh[2].push_back(TMath::Mean(rh.begin(),rh.end()));
            xWindvel[2].push_back(TMath::Mean(windvel.begin(),windvel.end()));
            xWinddir[2].push_back(TMath::Mean(winddir.begin(),winddir.end()));
            xGrad[2].push_back(TMath::Mean(grad.begin(),grad.end()));
        };
        if(data.year[i] == "2021")
        {
            cout << "Year read and compare " << data.year[i] << endl;
            for(int j = 0; j < data_meteo[3].timed.size(); j++)
            {
                if(data_meteo[3].timed_begin[j] >= data.timed_begin[i] &&  data_meteo[3].timed_begin[j] < data.timed_end[i])
                {
                    temp.push_back(data_meteo[3].ventemp[j]);
                    rh.push_back(data_meteo[3].rh2[j]);
                    windvel.push_back(data_meteo[3].windvel[j]);
                    winddir.push_back(data_meteo[3].winddir[j]);
                    grad.push_back(data_meteo[3].grad[j]);
                };
            };
            cout << "Averaging in " << data.year[i] << "||" << TMath::Mean(temp.begin(),temp.end()) << endl;
            xTemp[3].push_back(TMath::Mean(temp.begin(),temp.end()));
            xRh[3].push_back(TMath::Mean(rh.begin(),rh.end()));
            xWindvel[3].push_back(TMath::Mean(windvel.begin(),windvel.end()));
            xWinddir[3].push_back(TMath::Mean(winddir.begin(),winddir.end()));
            xGrad[3].push_back(TMath::Mean(grad.begin(),grad.end()));
        };
    };

    for(int i = 0; i < xTemp.size(); i ++)
    {
        for(int j = 0; j < xTemp[i].size(); j++)
        {
            cout << "Vectors " << i << "-" << j << ":" << xTemp[i][j] << "||" << xRh[i][j] << "||" << xWindvel[i][j] << "||" << xWinddir[i][j] << "||" << xGrad[i][j] << endl;
        };
    };

    for(int i = 0; i < data.O18.size();i++)
    {
        if(data.year[i] == "2018"){yO18[0].push_back(data.O18[i]);};
        if(data.year[i] == "2019"){yO18[1].push_back(data.O18[i]);};
        if(data.year[i] == "2020"){yO18[2].push_back(data.O18[i]);};
        if(data.year[i] == "2021"){yO18[3].push_back(data.O18[i]);};
        if(data.year[i] == "2018"){yH2[0].push_back(data.H2[i]);};
        if(data.year[i] == "2019"){yH2[1].push_back(data.H2[i]);};
        if(data.year[i] == "2020"){yH2[2].push_back(data.H2[i]);};
        if(data.year[i] == "2021"){yH2[3].push_back(data.H2[i]);};
        if(data.year[i] == "2018"){yDex[0].push_back(data.D_excess[i]);};
        if(data.year[i] == "2019"){yDex[1].push_back(data.D_excess[i]);};
        if(data.year[i] == "2020"){yDex[2].push_back(data.D_excess[i]);};
        if(data.year[i] == "2021"){yDex[3].push_back(data.D_excess[i]);};
    };

    string name_file_temp = evalpath + "/End/Temperature_Event.png";
    string name_file_rh = evalpath + "/End/RH_Event.png";
    string name_file_windvel = evalpath + "/End/Windvel_Event.png";
    string name_file_winddir = evalpath + "/End/Winddir_Event.png";
    string name_file_grad = evalpath + "/End/Grad_Event.png";

    string temp_ax = "ventilated Temperature[C {}^{#circ}]";
    string rh_ax = "relative Humidity[%]";
    string windvel_ax = "Wind velocity[m/s]";
    string winddir_ax = "Wind direction[{}^{#circ}]";
    string grad_ax = "Global Radiation [W/m^2]";
    string O18_ax = "#delta {}^{18}O";
    string H2_ax = "#delta {}^{2}H";
    string Dex_ax = "Dexcess";

    //drawing Graphs
    int width = 2000;
    int height = 2000;

    TCanvas *cTemp = new TCanvas("Temp","Temp",0,0,width,height*3);
    cTemp->Divide(2,6);
    cTemp->GetFrame()->SetBorderSize(12);
    cTemp->SetGrid();
    vector<TGraph*> grTempO18, grTempH2, grTempDex;

    for(int i = 0; i < 4; i++)
    {
        grTempO18.push_back(new TGraph(yO18[i].size(),xTemp[i].data(),yO18[i].data()));
        grTempH2.push_back(new TGraph(yH2[i].size(),xTemp[i].data(),yH2[i].data()));
        grTempDex.push_back(new TGraph(yDex[i].size(),xTemp[i].data(),yDex[i].data()));

        grTempO18[i]->SetTitle((year[i] + ";" + temp_ax + ";" + O18_ax).c_str());grTempH2[i]->SetTitle((year[i] + ";" + temp_ax + ";" + H2_ax).c_str());grTempDex[i]->SetTitle((year[i] + ";" + temp_ax + ";" + Dex_ax).c_str());
        grTempO18[i]->SetMarkerStyle(43);grTempH2[i]->SetMarkerStyle(43);grTempDex[i]->SetMarkerStyle(43);
        grTempO18[i]->SetMarkerColor(kBlack);grTempH2[i]->SetMarkerColor(kGreen);grTempDex[i]->SetMarkerColor(kOrange);
        grTempO18[i]->SetMarkerSize(4);grTempH2[i]->SetMarkerSize(4);grTempDex[i]->SetMarkerSize(4);
        grTempO18[i]->SetLineWidth(8);grTempH2[i]->SetLineWidth(8);grTempDex[i]->SetLineWidth(8);
        grTempO18[i]->SetLineColor(color[i]);grTempH2[i]->SetLineColor(color[i]);grTempDex[i]->SetLineColor(color[i]);
        cTemp->cd(i+1);grTempO18[i]->Draw("AP");
        cTemp->cd(i+5);grTempH2[i]->Draw("AP");
        cTemp->cd(i+9);grTempDex[i]->Draw("AP");
    };
    cTemp->Update();
    cTemp->Update();
    cTemp->Update();

    cTemp->Print(name_file_temp.c_str());
    cTemp->Close();

    TCanvas *cRH = new TCanvas("RH","RH",0,0,width,height*3);
    cRH->Divide(2,6);
    cRH->GetFrame()->SetBorderSize(12);
    cRH->SetGrid();
    vector<TGraph*> grRHO18, grRHH2, grRHDex;

    for(int i = 0; i < 4; i++)
    {
        grRHO18.push_back(new TGraph(yO18[i].size(),xRh[i].data(),yO18[i].data()));
        grRHH2.push_back(new TGraph(yH2[i].size(),xRh[i].data(),yH2[i].data()));
        grRHDex.push_back(new TGraph(yDex[i].size(),xRh[i].data(),yDex[i].data()));

        grRHO18[i]->SetTitle((year[i] + ";" + rh_ax + ";" + O18_ax).c_str());grRHH2[i]->SetTitle((year[i] + ";" + rh_ax + ";" + H2_ax).c_str());grRHDex[i]->SetTitle((year[i] + ";" + rh_ax + ";" + Dex_ax).c_str());
        grRHO18[i]->SetMarkerStyle(43);grRHH2[i]->SetMarkerStyle(43);grRHDex[i]->SetMarkerStyle(43);
        grRHO18[i]->SetMarkerColor(kBlack);grRHH2[i]->SetMarkerColor(kGreen);grRHDex[i]->SetMarkerColor(kOrange);
        grRHO18[i]->SetMarkerSize(4);grRHH2[i]->SetMarkerSize(4);grRHDex[i]->SetMarkerSize(4);
        grRHO18[i]->SetLineWidth(8);grRHH2[i]->SetLineWidth(8);grRHDex[i]->SetLineWidth(8);
        grRHO18[i]->SetLineColor(color[i]);grRHH2[i]->SetLineColor(color[i]);grRHDex[i]->SetLineColor(color[i]);
        cRH->cd(i+1);grRHO18[i]->Draw("AP");
        cRH->cd(i+5);grRHH2[i]->Draw("AP");
        cRH->cd(i+9);grRHDex[i]->Draw("AP");
    };
    cRH->Update();
    cRH->Update();
    cRH->Update();

    cRH->Print(name_file_rh.c_str());
    cRH->Close();

    TCanvas *cWindvel = new TCanvas("Windvel","Windvel",0,0,width,height*3);
    cWindvel->Divide(2,6);
    cWindvel->GetFrame()->SetBorderSize(12);
    cWindvel->SetGrid();
    vector<TGraph*> grWVO18, grWVH2, grWVDex;

    for(int i = 0; i < 4; i++)
    {
        grWVO18.push_back(new TGraph(yO18[i].size(),xWindvel[i].data(),yO18[i].data()));
        grWVH2.push_back(new TGraph(yH2[i].size(),xWindvel[i].data(),yH2[i].data()));
        grWVDex.push_back(new TGraph(yDex[i].size(),xWindvel[i].data(),yDex[i].data()));

        grWVO18[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + O18_ax).c_str());grWVH2[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + H2_ax).c_str());grWVDex[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + Dex_ax).c_str());
        grWVO18[i]->SetMarkerStyle(43);grWVH2[i]->SetMarkerStyle(43);grWVDex[i]->SetMarkerStyle(43);
        grWVO18[i]->SetMarkerColor(kBlack);grWVH2[i]->SetMarkerColor(kGreen);grWVDex[i]->SetMarkerColor(kOrange);
        grWVO18[i]->SetMarkerSize(4);grWVH2[i]->SetMarkerSize(4);grWVDex[i]->SetMarkerSize(4);
        grWVO18[i]->SetLineWidth(8);grWVH2[i]->SetLineWidth(8);grWVDex[i]->SetLineWidth(8);
        grWVO18[i]->SetLineColor(color[i]);grWVH2[i]->SetLineColor(color[i]);grWVDex[i]->SetLineColor(color[i]);
        cWindvel->cd(i+1);grWVO18[i]->Draw("AP");
        cWindvel->cd(i+5);grWVH2[i]->Draw("AP");
        cWindvel->cd(i+9);grWVDex[i]->Draw("AP");
    };
    cWindvel->Update();
    cWindvel->Update();
    cWindvel->Update();

    cWindvel->Print(name_file_windvel.c_str());
    cWindvel->Close();

    TCanvas *cWinddir = new TCanvas("Winddir","Winddir",0,0,width,height*3);
    cWinddir->Divide(2,6);
    cWinddir->GetFrame()->SetBorderSize(12);
    cWinddir->SetGrid();
    vector<TGraph*> grWDO18, grWDH2, grWDDex;

    for(int i = 0; i < 4; i++)
    {
        grWDO18.push_back(new TGraph(yO18[i].size(),xWinddir[i].data(),yO18[i].data()));
        grWDH2.push_back(new TGraph(yH2[i].size(),xWinddir[i].data(),yH2[i].data()));
        grWDDex.push_back(new TGraph(yDex[i].size(),xWinddir[i].data(),yDex[i].data()));

        grWDO18[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + O18_ax).c_str());grWDH2[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + H2_ax).c_str());grWDDex[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + Dex_ax).c_str());
        grWDO18[i]->SetMarkerStyle(43);grWDH2[i]->SetMarkerStyle(43);grWDDex[i]->SetMarkerStyle(43);
        grWDO18[i]->SetMarkerColor(kBlack);grWDH2[i]->SetMarkerColor(kGreen);grWDDex[i]->SetMarkerColor(kOrange);
        grWDO18[i]->SetMarkerSize(4);grWDH2[i]->SetMarkerSize(4);grWDDex[i]->SetMarkerSize(4);
        grWDO18[i]->SetLineWidth(8);grWDH2[i]->SetLineWidth(8);grWDDex[i]->SetLineWidth(8);
        grWDO18[i]->SetLineColor(color[i]);grWDH2[i]->SetLineColor(color[i]);grWDDex[i]->SetLineColor(color[i]);
        cWinddir->cd(i+1);grWDO18[i]->Draw("AP");
        cWinddir->cd(i+5);grWDH2[i]->Draw("AP");
        cWinddir->cd(i+9);grWDDex[i]->Draw("AP");
    };
    cWinddir->Update();
    cWinddir->Update();
    cWinddir->Update();

    cWinddir->Print(name_file_winddir.c_str());
    cWinddir->Close();

    TCanvas *cGrad = new TCanvas("Grad","Grad",0,0,width,height*3);
    cGrad->Divide(2,6);
    cGrad->GetFrame()->SetBorderSize(12);
    cGrad->SetGrid();
    vector<TGraph*> grGRADO18, grGRADH2, grGRADDex;

    for(int i = 0; i < 4; i++)
    {
        grGRADO18.push_back(new TGraph(yO18[i].size(),xGrad[i].data(),yO18[i].data()));
        grGRADH2.push_back(new TGraph(yH2[i].size(),xGrad[i].data(),yH2[i].data()));
        grGRADDex.push_back(new TGraph(yDex[i].size(),xGrad[i].data(),yDex[i].data()));

        grGRADO18[i]->SetTitle((year[i] + ";" + grad_ax + ";" + O18_ax).c_str());grGRADH2[i]->SetTitle((year[i] + ";" + grad_ax + ";" + H2_ax).c_str());grGRADDex[i]->SetTitle((year[i] + ";" + grad_ax + ";" + Dex_ax).c_str());
        grGRADO18[i]->SetMarkerStyle(43);grGRADH2[i]->SetMarkerStyle(43);grGRADDex[i]->SetMarkerStyle(43);
        grGRADO18[i]->SetMarkerColor(kBlack);grGRADH2[i]->SetMarkerColor(kGreen);grGRADDex[i]->SetMarkerColor(kOrange);
        grGRADO18[i]->SetMarkerSize(4);grGRADH2[i]->SetMarkerSize(4);grGRADDex[i]->SetMarkerSize(4);
        grGRADO18[i]->SetLineWidth(8);grGRADH2[i]->SetLineWidth(8);grGRADDex[i]->SetLineWidth(8);
        grGRADO18[i]->SetLineColor(color[i]);grGRADH2[i]->SetLineColor(color[i]);grGRADDex[i]->SetLineColor(color[i]);
        cGrad->cd(i+1);grGRADO18[i]->Draw("AP");
        cGrad->cd(i+5);grGRADH2[i]->Draw("AP");
        cGrad->cd(i+9);grGRADDex[i]->Draw("AP");
    };
    cGrad->Update();
    cGrad->Update();
    cGrad->Update();

    cGrad->Print(name_file_grad.c_str());
    cGrad->Close();

    vector<double> xTemp_all, yO18_all, yH2_all;
    for(int i = 0; i < xTemp.size(); i++)
    {
        for(int j = 0; j < xTemp[i].size(); j++)
        {
            xTemp_all.push_back(xTemp[i][j]);
            yO18_all.push_back(yO18[i][j]);
            yH2_all.push_back(yH2[i][j]);
        };
    };
    string name_file_tempfit = evalpath + "/End/Temperature_Event_fit.png";

    TF1 *linFit1 = new TF1("Tempfit1","pol1",-300,100);
    TF1 *linFit2 = new TF1("Tempfit2","pol1",-300,100);
    TCanvas *cAll = new TCanvas("All","All",0,0,width,height);
    cAll->Divide(1,2,0,0);
    cAll->GetFrame()->SetBorderSize(12);
    cAll->SetGrid();
    cAll->cd(1);
    TGraph *grO18All = new TGraph(xTemp_all.size(),xTemp_all.data(),yO18_all.data());
    grO18All->SetTitle(("Temperature vs. #delta {}^{18}O;"+temp_ax+";"+O18_ax).c_str());
    grO18All->SetMarkerStyle(43);
    grO18All->SetMarkerColor(kBlack);
    grO18All->SetMarkerSize(4);
    grO18All->SetLineWidth(3);
    grO18All->SetLineColor(kRed);
    grO18All->Draw("AP");
    linFit1->SetParameters(8.,10.);
    grO18All->Fit(linFit1);
    linFit1 = grO18All->GetFunction("Tempfit1");
    linFit1->SetLineColor(kRed);
    linFit1->SetLineStyle(4);
    linFit1->SetLineWidth(3);
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(2);
    TGraph *grH2All = new TGraph(xTemp_all.size(),xTemp_all.data(),yH2_all.data());
    grH2All->SetTitle(("Temperature vs. #delta {}^{2}H;"+temp_ax+";"+H2_ax).c_str());
    grH2All->SetMarkerStyle(43);
    grH2All->SetMarkerColor(kBlack);
    grH2All->SetMarkerSize(4);
    grH2All->SetLineWidth(3);
    grH2All->SetLineColor(kRed);
    grH2All->Draw("AP");
    linFit2->SetParameters(8.,10.);
    grH2All->Fit(linFit2);
    linFit2 = grH2All->GetFunction("Tempfit2");
    linFit2->SetLineColor(kRed);
    linFit2->SetLineStyle(4);
    linFit2->SetLineWidth(3);
    cAll->Update();
    gPad->Update();
    cAll->Update();
    cAll->Print(name_file_tempfit.c_str());
    cAll->Close();
};

//draw Graph Month Graphs
void drawMonthGraph(Data &data, string evalpath)
{
    string name_file_all = evalpath + "/End/" + "rain_month_year.png";
    string name_file_lmwl = evalpath + "/End/" + "rain_month_LMWL.png";

    string O18title = "{}^{18}O Graph;Month of Year;#delta^{18}O";
    string H2title = "{}^{2}H Graph;Month of Year;#delta^{2}H";
    string Dextitle = "D-excess Graph;Month of Year;D_{excess}";
    string H2Otitle = "Precipitation Graph;Month of Year;H_{2}O[g]";

    vector<vector<double>> xO18,xH2,xDex,xH2O;
    vector<vector<double>> yO18,yH2,yDex,yH2O;

    vector<string> year;

    vector<int> color;

    for(int i = 0; i < 4; i++)
    {
        xO18.push_back(vector<double>());
        xH2.push_back(vector<double>());
        xDex.push_back(vector<double>());
        xH2O.push_back(vector<double>());
        yO18.push_back(vector<double>());
        yH2.push_back(vector<double>());
        yDex.push_back(vector<double>());
        yH2O.push_back(vector<double>());
        if(i == 0){year.push_back("2018");color.push_back(1);};
        if(i == 1){year.push_back("2019");color.push_back(416);};
        if(i == 2){year.push_back("2020");color.push_back(600);};
        if(i == 3){year.push_back("2021");color.push_back(800);};
    };
    int width = 8000;
    int height = 3000;


    TCanvas *cAll = new TCanvas("All","All",0,0,width,height*4);
    cAll->Divide(1,4,0,0);
    cAll->GetFrame()->SetBorderSize(12);
    cAll->SetGrid();
    TMultiGraph *mgO18 = new TMultiGraph();
    mgO18->SetTitle(O18title.c_str());
    TMultiGraph *mgH2 = new TMultiGraph();
    mgH2->SetTitle(H2title.c_str());
    TMultiGraph *mgDex = new TMultiGraph();
    mgDex->SetTitle(Dextitle.c_str());
    TMultiGraph *mgH2O = new TMultiGraph();
    mgH2O->SetTitle(H2Otitle.c_str());
    vector<TGraph*> grO18, grH2, grDex, grH2O;
    for(int i = 0; i < data.O18.size();i++)
    {
        if(data.year[i] == "2018"){xO18[0].push_back(data.month[i]);yO18[0].push_back(data.O18[i]);};
        if(data.year[i] == "2019"){xO18[1].push_back(data.month[i]);yO18[1].push_back(data.O18[i]);};
        if(data.year[i] == "2020"){xO18[2].push_back(data.month[i]);yO18[2].push_back(data.O18[i]);};
        if(data.year[i] == "2021"){xO18[3].push_back(data.month[i]);yO18[3].push_back(data.O18[i]);};
        if(data.year[i] == "2018"){xH2[0].push_back(data.month[i]);yH2[0].push_back(data.H2[i]);};
        if(data.year[i] == "2019"){xH2[1].push_back(data.month[i]);yH2[1].push_back(data.H2[i]);};
        if(data.year[i] == "2020"){xH2[2].push_back(data.month[i]);yH2[2].push_back(data.H2[i]);};
        if(data.year[i] == "2021"){xH2[3].push_back(data.month[i]);yH2[3].push_back(data.H2[i]);};
        if(data.year[i] == "2018"){xDex[0].push_back(data.month[i]);yDex[0].push_back(data.D_excess[i]);};
        if(data.year[i] == "2019"){xDex[1].push_back(data.month[i]);yDex[1].push_back(data.D_excess[i]);};
        if(data.year[i] == "2020"){xDex[2].push_back(data.month[i]);yDex[2].push_back(data.D_excess[i]);};
        if(data.year[i] == "2021"){xDex[3].push_back(data.month[i]);yDex[3].push_back(data.D_excess[i]);};
        if(data.year[i] == "2018"){xH2O[0].push_back(data.month[i]);yH2O[0].push_back(data.H2O_m[i]);};
        if(data.year[i] == "2019"){xH2O[1].push_back(data.month[i]);yH2O[1].push_back(data.H2O_m[i]);};
        if(data.year[i] == "2020"){xH2O[2].push_back(data.month[i]);yH2O[2].push_back(data.H2O_m[i]);};
        if(data.year[i] == "2021"){xH2O[3].push_back(data.month[i]);yH2O[3].push_back(data.H2O_m[i]);};
    };
    for(int i = 0; i < 4; i++)
    {
        grO18.push_back(new TGraph(xO18[i].size(),xO18[i].data(),yO18[i].data()));
        grH2.push_back(new TGraph(xH2[i].size(),xH2[i].data(),yH2[i].data()));
        grDex.push_back(new TGraph(xDex[i].size(),xDex[i].data(),yDex[i].data()));
        grH2O.push_back(new TGraph(xH2O[i].size(),xH2O[i].data(),yH2O[i].data()));

        grO18[i]->SetTitle(year[i].c_str());grH2[i]->SetTitle(year[i].c_str());grDex[i]->SetTitle(year[i].c_str());grH2O[i]->SetTitle(year[i].c_str());
        grO18[i]->SetMarkerStyle(43);grH2[i]->SetMarkerStyle(43);grDex[i]->SetMarkerStyle(43);grH2O[i]->SetMarkerStyle(43);
        grO18[i]->SetMarkerColor(kBlack);grH2[i]->SetMarkerColor(kGreen);grDex[i]->SetMarkerColor(kOrange);grH2O[i]->SetMarkerColor(kBlue);
        grO18[i]->SetMarkerSize(12);grH2[i]->SetMarkerSize(12);grDex[i]->SetMarkerSize(12);grH2O[i]->SetMarkerSize(12);
        grO18[i]->SetLineWidth(8);grH2[i]->SetLineWidth(8);grDex[i]->SetLineWidth(8);grH2O[i]->SetLineWidth(8);
        grO18[i]->SetLineColor(color[i]);grH2[i]->SetLineColor(color[i]);grDex[i]->SetLineColor(color[i]);grH2O[i]->SetLineColor(color[i]);
        mgO18->Add(grO18[i]);mgH2->Add(grH2[i]);mgDex->Add(grDex[i]);mgH2O->Add(grH2O[i]);
    };


    cAll->cd(1);
    cout << "Drawing Month O18 ... " << endl;
    mgO18->Draw("ACP");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(2);
    cout << "Drawing Month H2 ... " << endl;
    mgH2->Draw("ACP");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(3);
    cout << "Drawing Month Dex ... " << endl;
    mgDex->Draw("ACP");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(4);
    cout << "Drawing Month H2O ... " << endl;
    mgH2O->Draw("ACP");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->Print(name_file_all.c_str());
    cAll->Close();


    vector<TGraph*> grLMWL;
    vector<TF1*> linFit;
    TCanvas *cLMWL = new TCanvas("LMWL","LMWL",0,0,2000,2000);
    cLMWL->Divide(2,2);
    vector<string> LMWLtitle;
    for(int i = 0; i < 4; i++)
    {
        cLMWL->cd(i+1);
        LMWLtitle.push_back("LMWL " + year[i] + ";#delta^{18}O;#delta^{2}H");
        linFit.push_back(new TF1(year[i].c_str(),"pol1",-100,100));
        grLMWL.push_back(new TGraph(yO18[i].size(),yO18[i].data(),yH2[i].data()));
        grLMWL[i]->SetTitle(LMWLtitle[i].c_str());
        grLMWL[i]->SetMarkerStyle(43);
        grLMWL[i]->SetMarkerColor(kBlue);
        grLMWL[i]->SetMarkerSize(3);
        grLMWL[i]->Draw("AP");
        linFit[i]->SetParameters(8.,10.);
        grLMWL[i]->Fit(linFit[i]);
        linFit[i] = grLMWL[i]->GetFunction(year[i].c_str());
        linFit[i]->SetLineColor(kRed);
        linFit[i]->SetLineStyle(4);
        linFit[i]->SetLineWidth(6);
        cLMWL->Update();
        gPad->Update();
        cLMWL->Update();
    };
    cLMWL->Print(name_file_lmwl.c_str());
    cLMWL->Close();
};

//draw Graph Month Meteo
void drawMonthMeteoGraph(Data &data, vector<Data> &data_meteo, string evalpath)
{
    vector<double> temp, rh, windvel, winddir, grad;
    vector<string> month, year;

    vector<vector<double>> xTemp,xRh,xWindvel,xWinddir,xGrad;
    vector<vector<double>> yO18,yH2,yDex;

    vector<int> color;

    //preparing vectors
    for(int i = 0; i < 4; i++)
    {
        xTemp.push_back(vector<double>());
        xRh.push_back(vector<double>());
        xWindvel.push_back(vector<double>());
        xWinddir.push_back(vector<double>());
        xGrad.push_back(vector<double>());
        yO18.push_back(vector<double>());
        yH2.push_back(vector<double>());
        yDex.push_back(vector<double>());
        if(i == 0){year.push_back("2018");color.push_back(1);};
        if(i == 1){year.push_back("2019");color.push_back(416);};
        if(i == 2){year.push_back("2020");color.push_back(600);};
        if(i == 3){year.push_back("2021");color.push_back(800);};
    };
    string last_month_str;
    for(int i = 0; i < data_meteo.size(); i++)
    {
        last_month_str = "01";
        for(int j = 0; j < data_meteo[i].timed.size(); j++)
        {

            if(data_meteo[i].month_str[j] != last_month_str)
            {
                //do average
                xTemp[i].push_back(TMath::Mean(temp.begin(),temp.end()));
                xRh[i].push_back(TMath::Mean(rh.begin(),rh.end()));
                xWindvel[i].push_back(TMath::Mean(windvel.begin(),windvel.end()));
                xWinddir[i].push_back(TMath::Mean(winddir.begin(),winddir.end()));
                xGrad[i].push_back(TMath::Mean(grad.begin(),grad.end()));
                //clear vectors
                temp.clear();
                rh.clear();
                windvel.clear();
                winddir.clear();
                grad.clear();
            };
            temp.push_back(data_meteo[i].ventemp[j]);
            rh.push_back(data_meteo[i].rh2[j]);
            windvel.push_back(data_meteo[i].windvel[j]);
            winddir.push_back(data_meteo[i].winddir[j]);
            grad.push_back(data_meteo[i].grad[j]);
            last_month_str = data_meteo[i].month_str[j];
            if(j == data_meteo[i].timed.size()-1)
            {
                //do average
                xTemp[i].push_back(TMath::Mean(temp.begin(),temp.end()));
                xRh[i].push_back(TMath::Mean(rh.begin(),rh.end()));
                xWindvel[i].push_back(TMath::Mean(windvel.begin(),windvel.end()));
                xWinddir[i].push_back(TMath::Mean(winddir.begin(),winddir.end()));
                xGrad[i].push_back(TMath::Mean(grad.begin(),grad.end()));
                //clear vectors
                cout << data_meteo[i].month_str[j] << " Month evaluated at year " << data_meteo[i].year[j] << endl;
                temp.clear();
                rh.clear();
                windvel.clear();
                winddir.clear();
                grad.clear();
            };
        };
    };
    for(int i = 0; i < xTemp.size(); i ++)
    {
        for(int j = 0; j < xTemp[i].size(); j++)
        {
            cout << "Vectors " << i << "-" << j << ":" << xTemp[i][j] << "||" << xRh[i][j] << "||" << xWindvel[i][j] << "||" << xWinddir[i][j] << "||" << xGrad[i][j] << endl;
        };
    };

    for(int i = 0; i < data.O18.size();i++)
    {
        if(data.year[i] == "2018"){yO18[0].push_back(data.O18[i]);};
        if(data.year[i] == "2019"){yO18[1].push_back(data.O18[i]);};
        if(data.year[i] == "2020"){yO18[2].push_back(data.O18[i]);};
        if(data.year[i] == "2021"){yO18[3].push_back(data.O18[i]);};
        if(data.year[i] == "2018"){yH2[0].push_back(data.H2[i]);};
        if(data.year[i] == "2019"){yH2[1].push_back(data.H2[i]);};
        if(data.year[i] == "2020"){yH2[2].push_back(data.H2[i]);};
        if(data.year[i] == "2021"){yH2[3].push_back(data.H2[i]);};
        if(data.year[i] == "2018"){yDex[0].push_back(data.D_excess[i]);};
        if(data.year[i] == "2019"){yDex[1].push_back(data.D_excess[i]);};
        if(data.year[i] == "2020"){yDex[2].push_back(data.D_excess[i]);};
        if(data.year[i] == "2021"){yDex[3].push_back(data.D_excess[i]);};
    };

    string name_file_temp = evalpath + "/End/Temperature_Month.png";
    string name_file_rh = evalpath + "/End/RH_Month.png";
    string name_file_windvel = evalpath + "/End/Windvel_Month.png";
    string name_file_winddir = evalpath + "/End/Winddir_Month.png";
    string name_file_grad = evalpath + "/End/Grad_Month.png";

    string temp_ax = "ventilated Temperature[#circC]";
    string rh_ax = "relative Humidity[%]";
    string windvel_ax = "Wind velocity[m/s]";
    string winddir_ax = "Wind direction[#circ]";
    string grad_ax = "Global Radiation[W/m^2]";
    string O18_ax = "#delta^{18}O";
    string H2_ax = "#delta^{2}H";
    string Dex_ax = "Dexcess";

    //drawing Graphs
    int width = 2000;
    int height = 2000;
    //special case 2021
    // double tempmean,rhmean,windvelmean,winddirmean,gradmean;
    // tempmean = (xTemp[3][0]+xTemp[3][1])/2.;
    // rhmean = (xRh[3][0]+xRh[3][1])/2.;
    // windvelmean = (xWindvel[3][0]+xWindvel[3][1])/2.;
    // winddirmean = (xWinddir[3][0]+xWinddir[3][1])/2.;
    // gradmean = (xGrad[3][0]+xGrad[3][1])/2.;
    // xTemp[3].erase(xTemp[3].begin()+0);
    // xRh[3].erase(xRh[3].begin()+0);
    // xWindvel[3].erase(xWindvel[3].begin()+0);
    // xWinddir[3].erase(xWinddir[3].begin()+0);
    // xGrad[3].erase(xGrad[3].begin()+0);
    // xTemp[3][0] = tempmean;
    // xRh[3][0] = rhmean;
    // xWindvel[3][0] = windvelmean;
    // xWinddir[3][0] = winddirmean;
    // xGrad[3][0] = gradmean;
    // tempmean = (xTemp[3][6]+xTemp[3][7])/2.;
    // rhmean = (xRh[3][6]+xRh[3][7])/2.;
    // windvelmean = (xWindvel[3][6]+xWindvel[3][7])/2.;
    // winddirmean = (xWinddir[3][6]+xWinddir[3][7])/2.;
    // gradmean = (xGrad[3][6]+xGrad[3][7])/2.;
    // xTemp[3].erase(xTemp[3].begin()+6);
    // xRh[3].erase(xRh[3].begin()+6);
    // xWindvel[3].erase(xWindvel[3].begin()+6);
    // xWinddir[3].erase(xWinddir[3].begin()+6);
    // xGrad[3].erase(xGrad[3].begin()+6);
    // xTemp[3][6] = tempmean;
    // xRh[3][6] = rhmean;
    // xWindvel[3][6] = windvelmean;
    // xWinddir[3][6] = winddirmean;
    // xGrad[3][6] = gradmean;

    TCanvas *cTemp = new TCanvas("Temp","Temp",0,0,width,height*3);
    cTemp->Divide(2,6);
    cTemp->GetFrame()->SetBorderSize(12);
    cTemp->SetGrid();
    vector<TGraph*> grTempO18, grTempH2, grTempDex;

    for(int i = 0; i < 4; i++)
    {
        grTempO18.push_back(new TGraph(yO18[i].size(),xTemp[i].data(),yO18[i].data()));
        grTempH2.push_back(new TGraph(yH2[i].size(),xTemp[i].data(),yH2[i].data()));
        grTempDex.push_back(new TGraph(yDex[i].size(),xTemp[i].data(),yDex[i].data()));

        grTempO18[i]->SetTitle((year[i] + ";" + temp_ax + ";" + O18_ax).c_str());grTempH2[i]->SetTitle((year[i] + ";" + temp_ax + ";" + H2_ax).c_str());grTempDex[i]->SetTitle((year[i] + ";" + temp_ax + ";" + Dex_ax).c_str());
        grTempO18[i]->SetMarkerStyle(43);grTempH2[i]->SetMarkerStyle(43);grTempDex[i]->SetMarkerStyle(43);
        grTempO18[i]->SetMarkerColor(kBlack);grTempH2[i]->SetMarkerColor(kGreen);grTempDex[i]->SetMarkerColor(kOrange);
        grTempO18[i]->SetMarkerSize(4);grTempH2[i]->SetMarkerSize(4);grTempDex[i]->SetMarkerSize(4);
        grTempO18[i]->SetLineWidth(8);grTempH2[i]->SetLineWidth(8);grTempDex[i]->SetLineWidth(8);
        grTempO18[i]->SetLineColor(color[i]);grTempH2[i]->SetLineColor(color[i]);grTempDex[i]->SetLineColor(color[i]);
        cTemp->cd(i+1);grTempO18[i]->Draw("AP");
        cTemp->cd(i+5);grTempH2[i]->Draw("AP");
        cTemp->cd(i+9);grTempDex[i]->Draw("AP");
    };
    cTemp->Update();
    cTemp->Update();
    cTemp->Update();

    cTemp->Print(name_file_temp.c_str());
    cTemp->Close();

    TCanvas *cRH = new TCanvas("RH","RH",0,0,width,height*3);
    cRH->Divide(2,6);
    cRH->GetFrame()->SetBorderSize(12);
    cRH->SetGrid();
    vector<TGraph*> grRHO18, grRHH2, grRHDex;

    for(int i = 0; i < 4; i++)
    {
        grRHO18.push_back(new TGraph(yO18[i].size(),xRh[i].data(),yO18[i].data()));
        grRHH2.push_back(new TGraph(yH2[i].size(),xRh[i].data(),yH2[i].data()));
        grRHDex.push_back(new TGraph(yDex[i].size(),xRh[i].data(),yDex[i].data()));

        grRHO18[i]->SetTitle((year[i] + ";" + rh_ax + ";" + O18_ax).c_str());grRHH2[i]->SetTitle((year[i] + ";" + rh_ax + ";" + H2_ax).c_str());grRHDex[i]->SetTitle((year[i] + ";" + rh_ax + ";" + Dex_ax).c_str());
        grRHO18[i]->SetMarkerStyle(43);grRHH2[i]->SetMarkerStyle(43);grRHDex[i]->SetMarkerStyle(43);
        grRHO18[i]->SetMarkerColor(kBlack);grRHH2[i]->SetMarkerColor(kGreen);grRHDex[i]->SetMarkerColor(kOrange);
        grRHO18[i]->SetMarkerSize(4);grRHH2[i]->SetMarkerSize(4);grRHDex[i]->SetMarkerSize(4);
        grRHO18[i]->SetLineWidth(8);grRHH2[i]->SetLineWidth(8);grRHDex[i]->SetLineWidth(8);
        grRHO18[i]->SetLineColor(color[i]);grRHH2[i]->SetLineColor(color[i]);grRHDex[i]->SetLineColor(color[i]);
        cRH->cd(i+1);grRHO18[i]->Draw("AP");
        cRH->cd(i+5);grRHH2[i]->Draw("AP");
        cRH->cd(i+9);grRHDex[i]->Draw("AP");
    };
    cRH->Update();
    cRH->Update();
    cRH->Update();

    cRH->Print(name_file_rh.c_str());
    cRH->Close();

    TCanvas *cWindvel = new TCanvas("Windvel","Windvel",0,0,width,height*3);
    cWindvel->Divide(2,6);
    cWindvel->GetFrame()->SetBorderSize(12);
    cWindvel->SetGrid();
    vector<TGraph*> grWVO18, grWVH2, grWVDex;

    for(int i = 0; i < 4; i++)
    {
        grWVO18.push_back(new TGraph(yO18[i].size(),xWindvel[i].data(),yO18[i].data()));
        grWVH2.push_back(new TGraph(yH2[i].size(),xWindvel[i].data(),yH2[i].data()));
        grWVDex.push_back(new TGraph(yDex[i].size(),xWindvel[i].data(),yDex[i].data()));

        grWVO18[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + O18_ax).c_str());grWVH2[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + H2_ax).c_str());grWVDex[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + Dex_ax).c_str());
        grWVO18[i]->SetMarkerStyle(43);grWVH2[i]->SetMarkerStyle(43);grWVDex[i]->SetMarkerStyle(43);
        grWVO18[i]->SetMarkerColor(kBlack);grWVH2[i]->SetMarkerColor(kGreen);grWVDex[i]->SetMarkerColor(kOrange);
        grWVO18[i]->SetMarkerSize(4);grWVH2[i]->SetMarkerSize(4);grWVDex[i]->SetMarkerSize(4);
        grWVO18[i]->SetLineWidth(8);grWVH2[i]->SetLineWidth(8);grWVDex[i]->SetLineWidth(8);
        grWVO18[i]->SetLineColor(color[i]);grWVH2[i]->SetLineColor(color[i]);grWVDex[i]->SetLineColor(color[i]);
        cWindvel->cd(i+1);grWVO18[i]->Draw("AP");
        cWindvel->cd(i+5);grWVH2[i]->Draw("AP");
        cWindvel->cd(i+9);grWVDex[i]->Draw("AP");
    };
    cWindvel->Update();
    cWindvel->Update();
    cWindvel->Update();

    cWindvel->Print(name_file_windvel.c_str());
    cWindvel->Close();

    TCanvas *cWinddir = new TCanvas("Winddir","Winddir",0,0,width,height*3);
    cWinddir->Divide(2,6);
    cWinddir->GetFrame()->SetBorderSize(12);
    cWinddir->SetGrid();
    vector<TGraph*> grWDO18, grWDH2, grWDDex;

    for(int i = 0; i < 4; i++)
    {
        grWDO18.push_back(new TGraph(yO18[i].size(),xWinddir[i].data(),yO18[i].data()));
        grWDH2.push_back(new TGraph(yH2[i].size(),xWinddir[i].data(),yH2[i].data()));
        grWDDex.push_back(new TGraph(yDex[i].size(),xWinddir[i].data(),yDex[i].data()));

        grWDO18[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + O18_ax).c_str());grWDH2[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + H2_ax).c_str());grWDDex[i]->SetTitle((year[i] + ";" + windvel_ax + ";" + Dex_ax).c_str());
        grWDO18[i]->SetMarkerStyle(43);grWDH2[i]->SetMarkerStyle(43);grWDDex[i]->SetMarkerStyle(43);
        grWDO18[i]->SetMarkerColor(kBlack);grWDH2[i]->SetMarkerColor(kGreen);grWDDex[i]->SetMarkerColor(kOrange);
        grWDO18[i]->SetMarkerSize(4);grWDH2[i]->SetMarkerSize(4);grWDDex[i]->SetMarkerSize(4);
        grWDO18[i]->SetLineWidth(8);grWDH2[i]->SetLineWidth(8);grWDDex[i]->SetLineWidth(8);
        grWDO18[i]->SetLineColor(color[i]);grWDH2[i]->SetLineColor(color[i]);grWDDex[i]->SetLineColor(color[i]);
        cWinddir->cd(i+1);grWDO18[i]->Draw("AP");
        cWinddir->cd(i+5);grWDH2[i]->Draw("AP");
        cWinddir->cd(i+9);grWDDex[i]->Draw("AP");
    };
    cWinddir->Update();
    cWinddir->Update();
    cWinddir->Update();

    cWinddir->Print(name_file_winddir.c_str());
    cWinddir->Close();

    TCanvas *cGrad = new TCanvas("Grad","Grad",0,0,width,height*3);
    cGrad->Divide(2,6);
    cGrad->GetFrame()->SetBorderSize(12);
    cGrad->SetGrid();
    vector<TGraph*> grGRADO18, grGRADH2, grGRADDex;

    for(int i = 0; i < 4; i++)
    {
        grGRADO18.push_back(new TGraph(yO18[i].size(),xWinddir[i].data(),yO18[i].data()));
        grGRADH2.push_back(new TGraph(yH2[i].size(),xWinddir[i].data(),yH2[i].data()));
        grGRADDex.push_back(new TGraph(yDex[i].size(),xWinddir[i].data(),yDex[i].data()));

        grGRADO18[i]->SetTitle((year[i] + ";" + grad_ax + ";" + O18_ax).c_str());grGRADH2[i]->SetTitle((year[i] + ";" + grad_ax + ";" + H2_ax).c_str());grGRADDex[i]->SetTitle((year[i] + ";" + grad_ax + ";" + Dex_ax).c_str());
        grGRADO18[i]->SetMarkerStyle(43);grGRADH2[i]->SetMarkerStyle(43);grGRADDex[i]->SetMarkerStyle(43);
        grGRADO18[i]->SetMarkerColor(kBlack);grGRADH2[i]->SetMarkerColor(kGreen);grGRADDex[i]->SetMarkerColor(kOrange);
        grGRADO18[i]->SetMarkerSize(4);grGRADH2[i]->SetMarkerSize(4);grGRADDex[i]->SetMarkerSize(4);
        grGRADO18[i]->SetLineWidth(8);grGRADH2[i]->SetLineWidth(8);grGRADDex[i]->SetLineWidth(8);
        grGRADO18[i]->SetLineColor(color[i]);grGRADH2[i]->SetLineColor(color[i]);grGRADDex[i]->SetLineColor(color[i]);
        cGrad->cd(i+1);grGRADO18[i]->Draw("AP");
        cGrad->cd(i+5);grGRADH2[i]->Draw("AP");
        cGrad->cd(i+9);grGRADDex[i]->Draw("AP");
    };
    cGrad->Update();
    cGrad->Update();
    cGrad->Update();

    cGrad->Print(name_file_grad.c_str());
    cGrad->Close();

    vector<double> xTemp_all, yO18_all, yH2_all;
    for(int i = 0; i < xTemp.size(); i++)
    {
        for(int j = 0; j < xTemp[i].size(); j++)
        {
            xTemp_all.push_back(xTemp[i][j]);
            yO18_all.push_back(yO18[i][j]);
            yH2_all.push_back(yH2[i][j]);
        };
    };
    string name_file_tempfit = evalpath + "/End/Temperature_Month_fit.png";

    TF1 *linFit1 = new TF1("Tempfit1","pol1",-300,100);
    TF1 *linFit2 = new TF1("Tempfit2","pol1",-300,100);
    TCanvas *cAll = new TCanvas("All","All",0,0,width,height);
    cAll->Divide(1,2,0,0);
    cAll->GetFrame()->SetBorderSize(12);
    cAll->SetGrid();
    cAll->cd(1);
    TGraph *grO18All = new TGraph(xTemp_all.size(),xTemp_all.data(),yO18_all.data());
    grO18All->SetTitle(("Temperature vs. #delta {}^{18}O;"+temp_ax+";"+O18_ax).c_str());
    grO18All->SetMarkerStyle(43);
    grO18All->SetMarkerColor(kBlack);
    grO18All->SetMarkerSize(4);
    grO18All->SetLineWidth(3);
    grO18All->SetLineColor(kRed);
    grO18All->GetXaxis()->SetTitleSize(0.046);
    grO18All->GetYaxis()->SetTitleSize(0.046);
    grO18All->Draw("AP");
    linFit1->SetParameters(8.,10.);
    grO18All->Fit(linFit1);
    linFit1 = grO18All->GetFunction("Tempfit1");
    linFit1->SetLineColor(kRed);
    linFit1->SetLineStyle(4);
    linFit1->SetLineWidth(3);
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(2);
    TGraph *grH2All = new TGraph(xTemp_all.size(),xTemp_all.data(),yH2_all.data());
    grH2All->SetTitle(("Temperature vs. #delta {}^{2}H;"+temp_ax+";"+H2_ax).c_str());
    grH2All->SetMarkerStyle(43);
    grH2All->SetMarkerColor(kBlack);
    grH2All->SetMarkerSize(4);
    grH2All->SetLineWidth(3);
    grH2All->SetLineColor(kRed);
    grH2All->GetXaxis()->SetTitleSize(0.046);
    grH2All->GetYaxis()->SetTitleSize(0.046);
    grH2All->Draw("AP");
    linFit2->SetParameters(8.,10.);
    grH2All->Fit(linFit2);
    linFit2 = grH2All->GetFunction("Tempfit2");
    linFit2->SetLineColor(kRed);
    linFit2->SetLineStyle(4);
    linFit2->SetLineWidth(3);
    cAll->Update();
    gPad->Update();
    cAll->Update();
    cAll->Print(name_file_tempfit.c_str());
    cAll->Close();
};


int main(int argc, char* argv[])
{
    //////////////////////////////////
    // Name, date and path to files //
    //////////////////////////////////
    char const * lFilterPatterns[1]={"*.csv"};
    char const * lFilterPatternsm[1]={"*.dat"};
    string datapath_event = tinyfd_openFileDialog("Choose File with Event data", ".", 1, lFilterPatterns, NULL, 0);
    cout << "Data Event file: " << datapath_event << endl;
    string datapath_month = tinyfd_openFileDialog("Choose File with Month data", ".", 1, lFilterPatterns, NULL, 0);
    cout << "Data Month file: " << datapath_month << endl;
    string datapath_flask = tinyfd_openFileDialog("Choose File with Flask data", ".", 1, lFilterPatterns, NULL, 0);
    cout << "Data Flask file: " << datapath_flask << endl;
    vector<string> datapath_meteo;
    for(int i = 0; i < 4; i++)
    {
        datapath_meteo.push_back(tinyfd_openFileDialog("Choose File with Meteo data", ".", 1, lFilterPatternsm, NULL, 0));
        cout << "Data meteo file:" << datapath_meteo[i] << endl;
    };

    //choose directory for evaluation data
    string evalpath = tinyfd_selectFolderDialog("Choose Folder for Evaluation", ".");
    cout << "Evalpath " << evalpath << endl;

    //////////////////////////////////////////////////
    // Read files and store values //
    //////////////////////////////////////////////////
    Data data_event, data_month, data_flask;
    vector<Data> data_meteo;

    for(int i = 0; i < 4; i++)
    {
        data_meteo.push_back(Data());
        getMeteo(datapath_meteo[i],data_meteo[i]);
    };
    getFlask(datapath_flask, data_flask);
    getEvent(datapath_event, data_event, data_flask);
    getMonth(datapath_month, data_month);

    drawEventGraph(data_event, evalpath);
    drawEventMeteoGraph(data_event, data_meteo, evalpath);
    drawMonthGraph(data_month, evalpath);
    drawMonthMeteoGraph(data_month, data_meteo,evalpath);




    return 0;
}