////////////////////////////////////////////////////////////////////////////
// Programm for evaluating ambient air data from the Picarro #insertModel //
// author: Sebastian Stezura // E-Mail: sebastian_stezura@gmx.net         //
// GitHub: CplusplusB                                                     //
// version: 1 // date: 03.01.2022                                         //
////////////////////////////////////////////////////////////////////////////

//get corrected data from (eval_air_std.cc) and get meteo data
//draw Graphs
//write ambient data and correlated meteo data to one file

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compile command:                                                                                                                  //
// g++-10 ambient_eval_meteo.cc tinyfiledialogs.c -ltbb `root-config --cflags --glibs --ldflags` -lMinuit -o ./Ambient_eval_meteo.o  //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

using namespace std;
namespace fs = std::filesystem;


//class for data storage
class Data
{
public:
    vector<string> identifier, interval, hour, month_str;//, interval;
    vector<double> timed, O18, H2, Dexcess, H2O, temp_amb, rh1_amb, rh2_amb, windvel_amb, winddir_amb, prec_amb;
    vector<double> windvel, contemp, rh1, rh2, grad, apress, o3g1, o3g3, no, ventemp, winddir, prec;
    vector<TDatime> date;
    string file_name = "0";
    double temp_max, temp_min;
    vector<int> month_all{1,2,3,4,5,6,7,8,9,10,11,12};
    vector<int> month;
    vector<string> month_names{"space", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"};

    Data(){};
    ~Data(){};
};

//getting Meteo Data from file
void getData_meteo(string datapath, Data &data)
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
                data.date.push_back(TDatime());
                data.date.back().Set
                (
                    stoi(interval_r.substr(0,4)), // year
                    stoi(interval_r.substr(5,2)), //month
                    stoi(interval_r.substr(8,2)), //day
                    stoi(interval_r.substr(10,2)), //hour
                    stoi(interval_r.substr(13,2)), //minute
                    stoi(interval_r.substr(16,2)) //second
                );
                remove(interval_r.begin(), interval_r.end(), '-');
                remove(interval_r.begin(), interval_r.end(), ':');
                interval_r = interval_r.substr(0,interval_r.size()-4);
                data.timed.push_back(stod(interval_r));
                data.interval.push_back(interval_r);
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
                data.month.push_back(stod(interval_r.substr(4,2)));
                if (i == 20 || i == 40){cout << name << ": " << data.date.back().Convert() << "||" << interval_r << "||" << windvel_r << "||" << contemp_r << "||" << rh1_r << "||" << rh2_r << "||" << grad_r << "||" << apress_r << "||" << o3g1_r << "||" << o3g3_r << "||" << no_r << "||" << ventemp_r << "||" << winddir_r << "||" << prec_r << endl;};
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

//getting Ambient Data from file and Correct
void getData_amb(string datapath, Data &data, string &year)
{
    string time_r, O18_r, H2_r, Dexcess_r, H2O_r;
    double timer, O18r, H2r, Dexcessr, H2Or;
    TDatime date_code;

    //loop over alle files
    ifstream inFile(datapath);
    cout << "Reading file at " << datapath << " ..." << endl;

    // read signal values from file
	if (inFile.is_open())
	{
        int i = 0;
        string line;
        while (getline(inFile, line))
		{
            if (i < 1) {i++; continue;};

            stringstream stst(line);
            getline(stst,time_r,',');
            getline(stst,O18_r,',');
            getline(stst,H2_r,',');
            getline(stst,Dexcess_r, ',');
            getline(stst,H2O_r,',');

            try
            {
                date_code.Set
                (
                    stoi(time_r.substr(0,4)), //year
                    stoi(time_r.substr(4,2)), //month
                    stoi(time_r.substr(6,2)), //day
                    stoi(time_r.substr(8,2)), //hour
                    stoi(time_r.substr(10,2)),//minutes
                    stoi(time_r.substr(12,2)) //seconds
                );
                timer = stod(time_r);
            }
            catch (...)
            {
                cout << "Time Code: " << time_r << endl;
            }

            if (time_r.substr(0,4) != year){continue;};

            if(time_r.substr(0,4) == year)
            {
                H2Or = stod(H2O_r);
                O18r = stod(O18_r);
                H2r = stod(H2_r);
                Dexcessr = stod(Dexcess_r);
                data.timed.push_back(timer);
                data.date.push_back(date_code);
                data.H2O.push_back(H2Or);
                data.O18.push_back(O18r);
                data.H2.push_back(H2r);
                data.Dexcess.push_back(Dexcessr);
                data.hour.push_back(time_r.substr(8,2));
                data.month_str.push_back(time_r.substr(4,2));
            };
            i++;
		};
	};
    inFile.close();

};

//find temperature relation
int compare_temp(double date_amb, vector<string> &meteo_date, int start)
{
    string date_amb_str = to_string(date_amb);
    string meteo_date_sstr;
    string date_amb_sstr = date_amb_str.substr(4,6);
    double meteo_d = 0.;
    double diff = 10000.;
    double allowed_diff = 900.;
    for(int i = start; i < meteo_date.size(); i++)
    {
        try
        {
            meteo_date_sstr = meteo_date[i].substr(4,6);
        }
        catch (...)
        {
            // cout << "no data at " << i << " " << meteo_date[i] << " " << meteo_date.size() << endl;
        }
        if(date_amb_sstr == meteo_date_sstr)
        {
            meteo_d = stod(meteo_date[i]);
            diff = TMath::Abs(date_amb - meteo_d);
            if(diff <= allowed_diff){return i;};
        };
    };
    for(int i = start/2; i < meteo_date.size(); i++)
    {
        meteo_d = stod(meteo_date[i]);
        diff = TMath::Abs(date_amb - meteo_d);
        if(diff <= allowed_diff + 2000.){return i;};
    };
    cout << "Next value " << start+1 << " at " << meteo_date[start+1] << "||" << date_amb_str << endl;
    return start+1;
};

//correlate meteo and ambient data
void getMeteo_amb(Data &data_amb, Data &data_meteo, Data &data_meteo_amb)
{
    int start = 0;
    for(int i = 0; i < data_amb.timed.size(); i++)
    {
        start = compare_temp(data_amb.timed[i], data_meteo.interval, start);
        data_meteo_amb.timed.push_back(data_meteo.timed[start]);
        data_meteo_amb.interval.push_back(data_meteo.interval[start]);
        data_meteo_amb.windvel.push_back(data_meteo.windvel[start]);
        data_meteo_amb.contemp.push_back(data_meteo.contemp[start]);
        data_meteo_amb.rh1.push_back(data_meteo.rh1[start]);
        data_meteo_amb.rh2.push_back(data_meteo.rh2[start]);
        data_meteo_amb.ventemp.push_back(data_meteo.ventemp[start]);
        data_meteo_amb.winddir.push_back(data_meteo.winddir[start]);
        data_meteo_amb.prec.push_back(data_meteo.prec[start]);
        data_meteo_amb.grad.push_back(data_meteo.grad[start]);
    };
    cout << "Size of data Meteo: " << data_meteo_amb.timed.size() << endl;
    cout << "Size of data amb: " << data_amb.timed.size() << endl;
};

//draw Graphs
void drawGraphs(Data &data_amb, Data &data_meteo, string evalpath, string year)
{
    string name_file_O18 = evalpath + "/End/" + year + "_O18_overview" + ".png";
    string name_file_H2 = evalpath + "/End/" + year + "_H2_overview" + ".png";
    string name_file_Dex = evalpath + "/End/" + year + "_Dex_overview" + ".png";

    vector<string> O18Title;
    vector<string> H2Title;
    vector<string> DexTitle;
    vector<string> name_para;
    vector<string> unit;
    vector<string> name_amb_para;

    //name parameters
    name_para.push_back("Time");
    name_para.push_back("Interval");
    name_para.push_back("Wind velocity");
    name_para.push_back("Convective Temperature");
    name_para.push_back("Relative Humidity 1");
    name_para.push_back("Relative Humidity 2");
    name_para.push_back("Ventilated Temperature");
    name_para.push_back("Wind Direction");
    name_para.push_back("Global Radiation");

    //name ambient parameters
    name_amb_para.push_back("{}^{18}O");
    name_amb_para.push_back("{}^{2}H");
    name_amb_para.push_back("Dexcess");

    //units
    unit.push_back("[Month/Day/Hour/Minute/Second]");
    unit.push_back("[Month/Day/Hour/Minute/Second]");
    unit.push_back("[m/s]");
    unit.push_back("[#circC]");
    unit.push_back("[%]");
    unit.push_back("[%]");
    unit.push_back("[#circC]");
    unit.push_back("[#circ]");
    unit.push_back("[W/m^2]");

    //Graphs vectors
    vector<TGraph*> grO18;
    vector<TGraph*> grH2;
    vector<TGraph*> grDex;


    //name files begin with Wind velocity
    for(int i = 2; i < name_para.size(); i++)
    {
        O18Title.push_back(name_amb_para[0] + "-" + name_para[i] + " Graph Ambient Air " + year + ";#delta{}^{18}O;" + name_para[i] + unit[i]);
        H2Title.push_back(name_amb_para[1] + "-" + name_para[i] + " Graph Ambient Air " + year + ";#delta{}^{2}H;" + name_para[i] + unit[i]);
        DexTitle.push_back(name_amb_para[2] + "-" + name_para[i] + " Graph Ambient Air " + year + ";D_{excess};" + name_para[i] + unit[i]);
    };

    int width = 2500;
    int height = 1500;

    TCanvas *cO18 = new TCanvas("O18","O18",0,0,width,height*3);

    cO18->Divide(3,3);
    cO18->GetFrame()->SetBorderSize(12);
    cO18->SetGrid();
    cO18->cd();

    grO18.push_back(new TGraph(data_amb.O18.size(), data_amb.O18.data(), data_meteo.windvel.data()));
    grO18.push_back(new TGraph(data_amb.O18.size(), data_amb.O18.data(), data_meteo.contemp.data()));
    grO18.push_back(new TGraph(data_amb.O18.size(), data_amb.O18.data(), data_meteo.rh1.data()));
    grO18.push_back(new TGraph(data_amb.O18.size(), data_amb.O18.data(), data_meteo.rh2.data()));
    grO18.push_back(new TGraph(data_amb.O18.size(), data_amb.O18.data(), data_meteo.ventemp.data()));
    grO18.push_back(new TGraph(data_amb.O18.size(), data_amb.O18.data(), data_meteo.winddir.data()));
    grO18.push_back(new TGraph(data_amb.O18.size(), data_amb.O18.data(), data_meteo.grad.data()));

    for(int i = 0; i < grO18.size(); i++)
    {
        grO18[i]->SetTitle(O18Title[i].c_str());
        grO18[i]->SetMarkerStyle(43);
        grO18[i]->SetMarkerColor(kRed);
        grO18[i]->SetMarkerSize(2);
        cO18->cd(i+1);
        grO18[i]->Draw("AP");
        cO18->Update();
        gPad->Update();
        cO18->Update();
        cout << "Drawing Graph " << i+1 << " for O18" << endl;
    };
    cO18->Print(name_file_O18.c_str());
    cO18->Close();

    TCanvas *cH2 = new TCanvas("H2","H2",0,0,width,height*3);

    cH2->Divide(3,3);
    cH2->GetFrame()->SetBorderSize(12);
    cH2->SetGrid();
    cH2->cd();

    grH2.push_back(new TGraph(data_amb.H2.size(), data_amb.H2.data(), data_meteo.windvel.data()));
    grH2.push_back(new TGraph(data_amb.H2.size(), data_amb.H2.data(), data_meteo.contemp.data()));
    grH2.push_back(new TGraph(data_amb.H2.size(), data_amb.H2.data(), data_meteo.rh1.data()));
    grH2.push_back(new TGraph(data_amb.H2.size(), data_amb.H2.data(), data_meteo.rh2.data()));
    grH2.push_back(new TGraph(data_amb.H2.size(), data_amb.H2.data(), data_meteo.ventemp.data()));
    grH2.push_back(new TGraph(data_amb.H2.size(), data_amb.H2.data(), data_meteo.winddir.data()));
    grH2.push_back(new TGraph(data_amb.H2.size(), data_amb.H2.data(), data_meteo.grad.data()));

    for(int i = 0; i < grH2.size(); i++)
    {
        grH2[i]->SetTitle(H2Title[i].c_str());
        grH2[i]->SetMarkerStyle(43);
        grH2[i]->SetMarkerColor(kRed);
        grH2[i]->SetMarkerSize(2);
        cH2->cd(i+1);
        grH2[i]->Draw("AP");
        cH2->Update();
        gPad->Update();
        cH2->Update();
        cout << "Drawing Graph " << i+1 << " for H2" << endl;

    };

    cH2->Print(name_file_H2.c_str());
    cH2->Close();


    TCanvas *cDex = new TCanvas("Dex","Dex",0,0,width,height*3);

    cDex->Divide(3,3);
    cDex->GetFrame()->SetBorderSize(12);
    cDex->SetGrid();
    cDex->cd();

    grDex.push_back(new TGraph(data_amb.Dexcess.size(), data_amb.Dexcess.data(), data_meteo.windvel.data()));
    grDex.push_back(new TGraph(data_amb.Dexcess.size(), data_amb.Dexcess.data(), data_meteo.contemp.data()));
    grDex.push_back(new TGraph(data_amb.Dexcess.size(), data_amb.Dexcess.data(), data_meteo.rh1.data()));
    grDex.push_back(new TGraph(data_amb.Dexcess.size(), data_amb.Dexcess.data(), data_meteo.rh2.data()));
    grDex.push_back(new TGraph(data_amb.Dexcess.size(), data_amb.Dexcess.data(), data_meteo.ventemp.data()));
    grDex.push_back(new TGraph(data_amb.Dexcess.size(), data_amb.Dexcess.data(), data_meteo.winddir.data()));
    grDex.push_back(new TGraph(data_amb.Dexcess.size(), data_amb.Dexcess.data(), data_meteo.grad.data()));

    for(int i = 0; i < grDex.size(); i++)
    {
        grDex[i]->SetTitle(DexTitle[i].c_str());
        grDex[i]->SetMarkerStyle(43);
        grDex[i]->SetMarkerColor(kRed);
        grDex[i]->SetMarkerSize(2);
        cDex->cd(i+1);
        grDex[i]->Draw("AP");
        cDex->Update();
        gPad->Update();
        cDex->Update();
        cout << "Drawing Graph " << i+1 << " for Dex" << endl;
    };

    cDex->Print(name_file_Dex.c_str());
    cDex->Close();
};
//write data
void writeData(Data &data_amb, Data &data_meteo, string evalpath, string year)
{
    string OutputFileName;
    OutputFileName = "Ambient_data_meteo_" + year + "_corr.txt";
    ofstream outFile (evalpath + "/End/" + OutputFileName);
    outFile << "Time," << "O18," << "H2," << "Dexcess," << "H2O" << endl;
    cout << "Writing Data amb corrected to: " << OutputFileName << endl;

    for (int j = 0; j < data_amb.timed.size(); j++)
    {
        outFile << fixed << setprecision(3) << data_amb.timed[j] << "," << data_amb.O18[j] << "," << data_amb.H2[j] << "," << data_amb.Dexcess[j] << "," << data_amb.H2O[j] << endl;
    };
};

//draw Diurnal Graphs
void drawGraphDiurnal(Data &data_amb, string evalpath, string year)
{
    vector<vector<double>> O18,H2,Dex;
    vector<double> yO18,yH2,yDex,hour;
    vector<double> yO18_err,yH2_err,yDex_err,hour_err;
    for(int i = 0; i < 24; i++)
    {
        hour.push_back(i);
        hour_err.push_back(0.);
        O18.push_back(vector<double>());
        H2.push_back(vector<double>());
        Dex.push_back(vector<double>());
    }
    for(int i = 0; i < data_amb.timed.size(); i++)
    {
        if(data_amb.hour[i] == "00"){O18[0].push_back(data_amb.O18[i]);H2[0].push_back(data_amb.H2[i]);Dex[0].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "01"){O18[1].push_back(data_amb.O18[i]);H2[1].push_back(data_amb.H2[i]);Dex[1].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "02"){O18[2].push_back(data_amb.O18[i]);H2[2].push_back(data_amb.H2[i]);Dex[2].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "03"){O18[3].push_back(data_amb.O18[i]);H2[3].push_back(data_amb.H2[i]);Dex[3].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "04"){O18[4].push_back(data_amb.O18[i]);H2[4].push_back(data_amb.H2[i]);Dex[4].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "05"){O18[5].push_back(data_amb.O18[i]);H2[5].push_back(data_amb.H2[i]);Dex[5].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "06"){O18[6].push_back(data_amb.O18[i]);H2[6].push_back(data_amb.H2[i]);Dex[6].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "07"){O18[7].push_back(data_amb.O18[i]);H2[7].push_back(data_amb.H2[i]);Dex[7].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "08"){O18[8].push_back(data_amb.O18[i]);H2[8].push_back(data_amb.H2[i]);Dex[8].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "09"){O18[9].push_back(data_amb.O18[i]);H2[9].push_back(data_amb.H2[i]);Dex[9].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "10"){O18[10].push_back(data_amb.O18[i]);H2[10].push_back(data_amb.H2[i]);Dex[10].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "11"){O18[11].push_back(data_amb.O18[i]);H2[11].push_back(data_amb.H2[i]);Dex[11].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "12"){O18[12].push_back(data_amb.O18[i]);H2[12].push_back(data_amb.H2[i]);Dex[12].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "13"){O18[13].push_back(data_amb.O18[i]);H2[13].push_back(data_amb.H2[i]);Dex[13].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "14"){O18[14].push_back(data_amb.O18[i]);H2[14].push_back(data_amb.H2[i]);Dex[14].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "15"){O18[15].push_back(data_amb.O18[i]);H2[15].push_back(data_amb.H2[i]);Dex[15].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "16"){O18[16].push_back(data_amb.O18[i]);H2[16].push_back(data_amb.H2[i]);Dex[16].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "17"){O18[17].push_back(data_amb.O18[i]);H2[17].push_back(data_amb.H2[i]);Dex[17].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "18"){O18[18].push_back(data_amb.O18[i]);H2[18].push_back(data_amb.H2[i]);Dex[18].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "19"){O18[19].push_back(data_amb.O18[i]);H2[19].push_back(data_amb.H2[i]);Dex[19].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "20"){O18[20].push_back(data_amb.O18[i]);H2[20].push_back(data_amb.H2[i]);Dex[20].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "21"){O18[21].push_back(data_amb.O18[i]);H2[21].push_back(data_amb.H2[i]);Dex[21].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "22"){O18[22].push_back(data_amb.O18[i]);H2[22].push_back(data_amb.H2[i]);Dex[22].push_back(data_amb.Dexcess[i]);};
        if(data_amb.hour[i] == "23"){O18[23].push_back(data_amb.O18[i]);H2[23].push_back(data_amb.H2[i]);Dex[23].push_back(data_amb.Dexcess[i]);};
    };
    for(int i = 0; i < 24; i++)
    {
        yO18.push_back(TMath::Mean(O18[i].begin(),O18[i].end()));
        yH2.push_back(TMath::Mean(H2[i].begin(),H2[i].end()));
        yDex.push_back(TMath::Mean(Dex[i].begin(),Dex[i].end()));
        yO18_err.push_back(TMath::StdDev(O18[i].begin(),O18[i].end())/(O18[i].size()-1));
        yH2_err.push_back(TMath::StdDev(H2[i].begin(),H2[i].end())/(H2[i].size()-1));
        yDex_err.push_back(TMath::StdDev(Dex[i].begin(),Dex[i].end())/(Dex[i].size()-1));

        cout << "Data mean||STDEV for O18 at " << i << ":" << yO18[i] << "||" << yO18_err[i] << endl;
        cout << "Data mean||STDEV for H2 at " << i << ":" << yH2[i] << "||" << yH2_err[i] << endl;
        cout << "Data mean||STDEV for Dex at " << i << ":" << yDex[i] << "||" << yDex_err[i] << endl;
    };
    double mean_O18, mean_H2, mean_Dex;
    double mean_O18_err, mean_H2_err, mean_Dex_err;
    mean_O18 = TMath::Mean(yO18.begin(),yO18.end());
    mean_H2 = TMath::Mean(yH2.begin(),yH2.end());
    mean_Dex = TMath::Mean(yDex.begin(),yDex.end());
    mean_O18_err = TMath::StdDev(yO18.begin(),yO18.end());
    mean_H2_err = TMath::StdDev(yH2.begin(),yH2.end());
    mean_Dex_err = TMath::StdDev(yDex.begin(),yDex.end());

    double x,xerr;
    for(int i = 0; i < 24; i++)
    {
        //O18
        x = yO18[i]/mean_O18;
        xerr = TMath::Sqrt( (yO18_err[i]/mean_O18) * (yO18_err[i]/mean_O18) + (mean_O18_err * yO18[i] / mean_O18 / mean_O18) * (mean_O18_err * yO18[i] / mean_O18 / mean_O18) );
        yO18[i] = x;
        yO18_err[i] = xerr;
        //H2
        x = yH2[i]/mean_H2;
        xerr = TMath::Sqrt( (yH2_err[i]/mean_H2) * (yH2_err[i]/mean_H2) + (mean_H2_err * yH2[i] / mean_H2 / mean_H2) * (mean_H2_err * yH2[i] / mean_H2 / mean_H2) );
        yH2[i] = x;
        yH2_err[i] = xerr;
        //Dex
        x = yDex[i]/mean_Dex;
        xerr = TMath::Sqrt( (yDex_err[i]/mean_Dex) * (yDex_err[i]/mean_Dex) + (mean_Dex_err * yDex[i] / mean_Dex / mean_Dex) * (mean_Dex_err * yDex[i] / mean_Dex / mean_Dex) );
        yDex[i] = x;
        yDex_err[i] = xerr;

        cout << "Value O18: " << yO18[i] << "+-" << yO18_err[i] << endl;
        cout << "Value H2: " << yH2[i] << "+-" << yH2_err[i] << endl;
        cout << "Value Dex: " << yDex[i] << "+-" << yDex_err[i] << endl;
        cout << "========================" << endl;
    };

    string name_file = evalpath + "/End/" + year + "_diurnal_ambient.png";
    string Title = "Diurnal Graph " + year + ";Time[hour];normalised level";

    int width = 2000;
    int height = 2000;
    TCanvas *cAll = new TCanvas("All","All",0,0,width,height);
    cAll->GetFrame()->SetBorderSize(12);
    cAll->SetGrid();
    cAll->cd();

    TMultiGraph *mgr = new TMultiGraph();
    mgr->SetTitle(Title.c_str());

    TGraphErrors *grO18 = new TGraphErrors(hour.size(),hour.data(),yO18.data(),hour_err.data(),yO18_err.data());
    grO18->SetFillColor(kBlack);
    grO18->SetFillColorAlpha(kBlack,0.3);
    grO18->SetFillStyle(3003);
    grO18->SetLineWidth(4);
    grO18->SetLineColor(kBlack);
    TGraphErrors *grH2 = new TGraphErrors(hour.size(),hour.data(),yH2.data(),hour_err.data(),yH2_err.data());
    grH2->SetFillColor(kGreen);
    grH2->SetFillColorAlpha(kGreen,0.3);
    grH2->SetFillStyle(3003);
    grH2->SetLineWidth(4);
    grH2->SetLineColor(kGreen);
    TGraphErrors *grDex = new TGraphErrors(hour.size(),hour.data(),yDex.data(),hour_err.data(),yDex_err.data());
    grDex->SetFillColor(kOrange);
    grDex->SetFillColorAlpha(kOrange,0.3);
    grDex->SetFillStyle(3003);
    grDex->SetLineWidth(4);
    grDex->SetLineColor(kOrange);

    mgr->Add(grO18,"L");
    mgr->Add(grH2,"L");
    mgr->Add(grDex,"L");

    mgr->Draw("a3"); //a4 for smoothed, a3 for direct

    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->Print(name_file.c_str());
    cAll->Close();
};

//draw seasonal Graphs
void drawGraphSeason(Data &data_amb, Data &data_meteo, string evalpath, string year)
{
    //preparing vectors
    vector<vector<double>> xO18, yH2, temp, time;
    vector<string> Title, tempTitle;
    double max_O18, min_O18;
    double max_H2, min_H2;

    vector<TF1*> linFit;
    vector<string> fitname;

    // vector<TGraph*> grTemp;

    max_O18 = *std::max_element(data_amb.O18.begin(),data_amb.O18.end());
    min_O18 = *std::min_element(data_amb.O18.begin(),data_amb.O18.end());
    max_H2 = *std::max_element(data_amb.H2.begin(),data_amb.H2.end());
    min_H2 = *std::min_element(data_amb.H2.begin(),data_amb.H2.end());

    //seasons 0=Winter (Dec-Feb), 1=Spring (Ma-May), 2=Summer (Jun-Aug), 3=Autumn
    for(int i = 0; i < 4; i++)
    {
        xO18.push_back(vector<double>());
        yH2.push_back(vector<double>());
        temp.push_back(vector<double>());
        time.push_back(vector<double>());
        if(i == 0){Title.push_back(year + " Winter LMWL Ambient;#delta^{18}O;#delta^{2}H");};
        if(i == 1){Title.push_back(year + " Spring LMWL Ambient;#delta^{18}O;#delta^{2}H");};
        if(i == 2){Title.push_back(year + " Summer LMWL Ambient;#delta^{18}O;#delta^{2}H");};
        if(i == 3){Title.push_back(year + " Autumn LMWL Ambient;#delta^{18}O;#delta^{2}H");};
        if(i == 0){tempTitle.push_back(year + " Winter Temperature;Month;Temperature[C #circC]");};
        if(i == 1){tempTitle.push_back(year + " Spring Temperature;Month;Temperature[C #circC]");};
        if(i == 2){tempTitle.push_back(year + " Summer Temperature;Month;Temperature[C #circC]");};
        if(i == 3){tempTitle.push_back(year + " Autumn Temperature;Month;Temperature[C #circC]");};
        fitname.push_back("f" + i);
        linFit.push_back(new TF1(fitname[i].c_str(),"pol1",min_O18,max_O18));
    };
    for(int i = 0; i < data_amb.timed.size(); i++)
    {
        if(data_amb.month_str[i] == "01" || data_amb.month_str[i] == "02" || data_amb.month_str[i] == "12")
        {
            xO18[0].push_back(data_amb.O18[i]);
            yH2[0].push_back(data_amb.H2[i]);
        };

        if(data_amb.month_str[i] == "03" || data_amb.month_str[i] == "04" || data_amb.month_str[i] == "05")
        {
            xO18[1].push_back(data_amb.O18[i]);
            yH2[1].push_back(data_amb.H2[i]);
        };

        if(data_amb.month_str[i] == "06" || data_amb.month_str[i] == "07" || data_amb.month_str[i] == "08")
        {
            xO18[2].push_back(data_amb.O18[i]);
            yH2[2].push_back(data_amb.H2[i]);
        };

        if(data_amb.month_str[i] == "09" || data_amb.month_str[i] == "10" || data_amb.month_str[i] == "11")
        {
            xO18[3].push_back(data_amb.O18[i]);
            yH2[3].push_back(data_amb.H2[i]);
        };
    };
    for(int i = 0; i < data_meteo.ventemp.size(); i++)
    {
        if(data_meteo.month_str[i] == "01" || data_meteo.month_str[i] == "02" || data_meteo.month_str[i] == "12")
        {
            temp[0].push_back(data_meteo.ventemp[i]);
            time[0].push_back(data_meteo.month[i]);
        };

        if(data_meteo.month_str[i] == "03" || data_meteo.month_str[i] == "04" || data_meteo.month_str[i] == "05")
        {
            temp[1].push_back(data_meteo.ventemp[i]);
            time[1].push_back(data_meteo.month[i]);
        };

        if(data_meteo.month_str[i] == "06" || data_meteo.month_str[i] == "07" || data_meteo.month_str[i] == "08")
        {
            temp[2].push_back(data_meteo.ventemp[i]);
            time[2].push_back(data_meteo.month[i]);
        };

        if(data_meteo.month_str[i] == "09" || data_meteo.month_str[i] == "10" || data_meteo.month_str[i] == "11")
        {
            temp[3].push_back(data_meteo.ventemp[i]);
            time[3].push_back(data_meteo.month[i]);
        };
    };
    for(int i = 0; i < 4; i++)
    {
        cout << "Size of O18: " << i << "=" << xO18[i].size() << endl;
        cout << "Size of temp: " << i << "=" << temp[i].size() << endl;
    };

    //drawing
    string name_file = evalpath + "/End/" + year + "_LMWL_season.png";

    int width = 2000;
    int height = 2000;

    TCanvas *cAll = new TCanvas("All","All",0,0,width,height);
    cAll->Divide(2,2,0,0);
    cAll->GetFrame()->SetBorderSize(12);
    cAll->SetGrid();

    TF1 *lmwl = new TF1("LMWL","[0]+[1]*x",min_O18,max_O18);
    if(year == "2021"){lmwl->SetParameters(3.389,7.556);};
    if(year == "2020"){lmwl->SetParameters(3.850,7.572);};
    if(year == "2019"){lmwl->SetParameters(2.438,7.560);};
    if(year == "2018"){lmwl->SetParameters(1.947,7.617);};
    lmwl->SetLineColor(kBlack);
    lmwl->SetLineStyle(5);
    lmwl->SetLineWidth(3);

    cAll->cd(1);
    cout << "Drawing Graph Winter" << endl;
    TGraph *grWinter = new TGraph(xO18[0].size(),xO18[0].data(),yH2[0].data());
    grWinter->SetTitle(Title[0].c_str());
    grWinter->SetMarkerStyle(43);
    grWinter->SetMarkerColor(kBlue);
    grWinter->SetMarkerSize(2);
    grWinter->Draw("AP");
    linFit[0]->SetParameters(8.,10.);
    grWinter->Fit(linFit[0]);
    linFit[0] = grWinter->GetFunction(fitname[0].c_str());
    linFit[0]->SetLineColor(kRed);
    linFit[0]->SetLineStyle(4);
    linFit[0]->SetLineWidth(5);
    lmwl->Draw("SAME");
    grWinter->GetXaxis()->SetRangeUser(min_O18,max_O18);
    grWinter->GetYaxis()->SetRangeUser(min_H2,max_H2);
    cAll->Update();
    gPad->Update();
    cAll->Update();
    cout << "Drawing Graph Spring" << endl;
    cAll->cd(2);
    TGraph *grSpring = new TGraph(xO18[1].size(),xO18[1].data(),yH2[1].data());
    grSpring->SetTitle(Title[1].c_str());
    grSpring->SetMarkerStyle(43);
    grSpring->SetMarkerColor(kGreen);
    grSpring->SetMarkerSize(2);
    grSpring->Draw("AP");
    linFit[1]->SetParameters(8.,10.);
    grSpring->Fit(linFit[1]);
    linFit[1] = grSpring->GetFunction(fitname[1].c_str());
    linFit[1]->SetLineColor(kRed);
    linFit[1]->SetLineStyle(4);
    linFit[1]->SetLineWidth(5);
    lmwl->Draw("SAME");
    grSpring->GetXaxis()->SetRangeUser(min_O18,max_O18);
    grSpring->GetYaxis()->SetRangeUser(min_H2,max_H2);
    cAll->Update();
    gPad->Update();
    cAll->Update();
    cout << "Drawing Graph Summer" << endl;
    cAll->cd(3);
    TGraph *grSummer = new TGraph(xO18[2].size(),xO18[2].data(),yH2[2].data());
    grSummer->SetTitle(Title[2].c_str());
    grSummer->SetMarkerStyle(43);
    grSummer->SetMarkerColor(kOrange);
    grSummer->SetMarkerSize(2);
    grSummer->Draw("AP");
    linFit[2]->SetParameters(8.,10.);
    grSummer->Fit(linFit[2]);
    linFit[2] = grSummer->GetFunction(fitname[2].c_str());
    linFit[2]->SetLineColor(kRed);
    linFit[2]->SetLineStyle(4);
    linFit[2]->SetLineWidth(5);
    lmwl->Draw("SAME");
    grSummer->GetXaxis()->SetRangeUser(min_O18,max_O18);
    grSummer->GetYaxis()->SetRangeUser(min_H2,max_H2);
    cAll->Update();
    gPad->Update();
    cAll->Update();
    cout << "Drawing Graph Autumn" << endl;
    cAll->cd(4);
    TGraph *grAutumn = new TGraph(xO18[3].size(),xO18[3].data(),yH2[3].data());
    grAutumn->SetTitle(Title[3].c_str());
    grAutumn->SetMarkerStyle(43);
    grAutumn->SetMarkerColor(kOrange+4);
    grAutumn->SetMarkerSize(2);
    grAutumn->Draw("AP");
    linFit[3]->SetParameters(8.,10.);
    grAutumn->Fit(linFit[3]);
    linFit[3] = grAutumn->GetFunction(fitname[3].c_str());
    linFit[3]->SetLineColor(kRed);
    linFit[3]->SetLineStyle(4);
    linFit[3]->SetLineWidth(5);
    lmwl->Draw("SAME");
    grAutumn->GetXaxis()->SetRangeUser(min_O18,max_O18);
    grAutumn->GetYaxis()->SetRangeUser(min_H2,max_H2);
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->Print(name_file.c_str());
    cAll->Close();

    cout << "Mean Temperature w/s/s/a:";
    cout << TMath::Mean(temp[0].begin(),temp[0].end()) << "/";
    cout << TMath::Mean(temp[1].begin(),temp[1].end()) << "/";
    cout << TMath::Mean(temp[2].begin(),temp[2].end()) << "/";
    cout << TMath::Mean(temp[3].begin(),temp[3].end()) << endl;

};

//draw Temp
void drawGraphTemp(Data &data_amb, Data &data_meteo, string evalpath, string year)
{
    string name_file = evalpath + "/End/" + year + "_temp_corr" + ".png";
    string O18Title = "Graph Ambient vs. vent. Temperature"+ year+";#delta{}^{18}O;ventilated Temperature[#circC]";
    vector<vector<double>> xO18,yTemp;
    vector<TGraph*> gr;
    vector<int> color;
    for(int i = 0; i < 12; i++)
    {
        xO18.push_back(vector<double>());
        yTemp.push_back(vector<double>());
    };
    color.push_back(1); //black 1
    color.push_back(920); //gray 2
    color.push_back(632); //red 3
    color.push_back(416); //green 4
    color.push_back(600); //blue 5
    color.push_back(400); //yellow 6
    color.push_back(616); //Magenta 7
    color.push_back(432); //cyan 8
    color.push_back(800); //orange 9
    color.push_back(824); //spring+4 10
    color.push_back(845); //teal+5 11
    color.push_back(866); //azura+6 12
    for(int i = 0; i < data_amb.O18.size(); i++)
    {
        if(data_amb.month_str[i] == "01"){xO18[0].push_back(data_amb.O18[i]);yTemp[0].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "02"){xO18[1].push_back(data_amb.O18[i]);yTemp[1].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "03"){xO18[2].push_back(data_amb.O18[i]);yTemp[2].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "04"){xO18[3].push_back(data_amb.O18[i]);yTemp[3].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "05"){xO18[4].push_back(data_amb.O18[i]);yTemp[4].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "06"){xO18[5].push_back(data_amb.O18[i]);yTemp[5].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "07"){xO18[6].push_back(data_amb.O18[i]);yTemp[6].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "08"){xO18[7].push_back(data_amb.O18[i]);yTemp[7].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "09"){xO18[8].push_back(data_amb.O18[i]);yTemp[8].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "10"){xO18[9].push_back(data_amb.O18[i]);yTemp[9].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "11"){xO18[10].push_back(data_amb.O18[i]);yTemp[10].push_back(data_meteo.ventemp[i]);};
        if(data_amb.month_str[i] == "12"){xO18[11].push_back(data_amb.O18[i]);yTemp[11].push_back(data_meteo.ventemp[i]);};
    };

    TMultiGraph *mgO18 = new TMultiGraph();
    mgO18->SetTitle(O18Title.c_str());
    for(int i = 0; i < 12; i++)
    {
        gr.push_back(new TGraph(xO18[i].size(),xO18[i].data(),yTemp[i].data()));
        gr[i]->SetMarkerColor(color[i]);
        gr[i]->SetMarkerSize(11);
        mgO18->Add(gr[i]);
    };
    TCanvas *cO18 = new TCanvas("All","All",0,0,2000,2000);
    cO18->SetGrid();
    cO18->cd();
    mgO18->Draw("AP");
    cO18->Update();
    gPad->Update();
    cO18->Update();
    cO18->Print(name_file.c_str());
    cO18->Close();
};

//draw monthly ambient
void drawGraphMeanYear(Data &data_amb, Data &data_meteo, string evalpath, string year)
{
    string name_file = evalpath + "/End/" + year + "_mean_year.png";

    string O18Title = "Monthly mean Ambient Air Graph "+year+";Month;#delta{}^{18}O";
    string H2Title = "Monthly mean Ambient Air Graph "+year+";Month;#delta{}^{2}H";
    string TempTitle = "Monthly mean Ambient Air Graph "+year+";Month;vententilated Temperature[#circC]";

    TGraph *grO18, *grH2, *grTemp;
    vector<double> yO18, yH2, yTemp;
    vector<double> xMonth;
    vector<double> O18, H2, Temp;
    int month = 1;
    if(year == "2019"){month = 4;};
    if(year == "2018"){month = 4;};

    for(int i = 0; i < data_amb.timed.size(); i++)
    {
        if(stoi(data_amb.month_str[i]) > month)
        {
            yO18.push_back(TMath::Mean(O18.begin(),O18.end()));
            yH2.push_back(TMath::Mean(H2.begin(),H2.end()));
            yTemp.push_back(TMath::Mean(Temp.begin(),Temp.end()));
            cout << "Size of month:" << month << "=" << O18.size() << "|" << H2.size() << "|" << Temp.size() << endl;
            cout << "Size of month:" << month << "=" << yO18.size() << "|" << yH2.size() << "|" << yTemp.size() << endl;
            O18.clear();
            H2.clear();
            Temp.clear();
            xMonth.push_back(month);
            month = stoi(data_amb.month_str[i]);
        };
        if(stoi(data_amb.month_str[i]) == month)
        {
            O18.push_back(data_amb.O18[i]);
            H2.push_back(data_amb.H2[i]);
            Temp.push_back(data_meteo.ventemp[i]);
        };
    };
    yO18.push_back(TMath::Mean(O18.begin(),O18.end()));
    yH2.push_back(TMath::Mean(H2.begin(),H2.end()));
    yTemp.push_back(TMath::Mean(Temp.begin(),Temp.end()));
    cout << "Size of month:" << month << "=" << O18.size() << "|" << H2.size() << "|" << Temp.size() << endl;
    cout << "Size of month:" << month << "=" << yO18.size() << "|" << yH2.size() << "|" << yTemp.size() << endl;
    O18.clear();
    H2.clear();
    Temp.clear();
    xMonth.push_back(month);
    grO18 = new TGraph(yO18.size(),xMonth.data(),yO18.data());
    grH2 = new TGraph(yH2.size(),xMonth.data(),yH2.data());
    grTemp = new TGraph(yTemp.size(),xMonth.data(),yTemp.data());

    TCanvas *cAll = new TCanvas("c", "c", 0,0,2000,3000);
    cAll->SetGrid();
    cAll->Divide(1,3,0,0);
    cAll->cd(1);
    grO18->SetTitle(O18Title.c_str());
    grO18->SetMarkerStyle(43);
    grO18->SetMarkerSize(4);
    grO18->SetLineColor(kBlack);
    grO18->SetLineWidth(2);
    grO18->Draw("ALP");

    cAll->cd(2);
    grH2->SetTitle(H2Title.c_str());
    grH2->SetMarkerStyle(43);
    grH2->SetMarkerSize(4);
    grH2->SetLineColor(kGreen);
    grH2->SetLineWidth(2);
    grH2->Draw("ALP");

    cAll->cd(3);
    grTemp->SetTitle(TempTitle.c_str());
    grTemp->SetMarkerStyle(43);
    grTemp->SetMarkerSize(4);
    grTemp->SetLineColor(kRed);
    grTemp->SetLineWidth(2);
    grTemp->Draw("ALP");

    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->Print(name_file.c_str());
    cAll->Close();
};

int main(int argc, char* argv[])
{
    //////////////////////////////////
    // Name, date and path to files //
    //////////////////////////////////
    string year;
    cout << "Which year?" << endl;
    cin >> year;
    char const * lFilterPatterns[1]={"*.txt"};
    char const * lFilterPatterns1[1]={"*.dat"};
    string datapath_amb = tinyfd_openFileDialog("Choose File with Ambient data", ".", 1, lFilterPatterns, NULL, 0);
    cout << "Data Ambient file: " << datapath_amb << endl;

    string datapath_meteo = tinyfd_openFileDialog("Choose File with Meteo data", ".", 1, lFilterPatterns1, NULL, 0);
    cout << "Data Meteo file: " << datapath_meteo << endl;

    //choose directory for evaluation data
    string evalpath = tinyfd_selectFolderDialog("Choose Folder for Evaluation", ".");
    cout << "Evalpath " << evalpath << endl;

    //////////////////////////////////////////////////
    // Read files and store values //
    //////////////////////////////////////////////////
    Data data_amb, data_meteo, data_meteo_amb;

    cout << "################" << endl << "Reading Meteo data ..." << endl;
    getData_meteo(datapath_meteo, data_meteo);
    cout << "Finished." << endl << "################" << endl << "Reading Ambient Air data ..." << endl;
    getData_amb(datapath_amb, data_amb, year);
    cout << "Finished." << endl;
    cout << "################" << endl << "Comparing meteo and ambient ... " << endl;

    getMeteo_amb(data_amb, data_meteo, data_meteo_amb);
    cout << "Finished." << endl;
    cout << "################" << endl << "Drawing Graphs ... " << endl;
    drawGraphs(data_amb, data_meteo_amb, evalpath, year);
    cout << "Finished." << endl;
    cout << "################" << endl;
    writeData(data_amb, data_meteo_amb, evalpath, year);

    drawGraphDiurnal(data_amb, evalpath, year);
    drawGraphSeason(data_amb, data_meteo, evalpath, year);
    drawGraphTemp(data_amb, data_meteo_amb, evalpath, year);
    drawGraphMeanYear(data_amb, data_meteo_amb, evalpath, year);


    return 0;
}