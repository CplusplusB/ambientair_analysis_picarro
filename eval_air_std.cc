////////////////////////////////////////////////////////////////////////////
// Programm for evaluating ambient air data from the Picarro #insertModel //
// author: Sebastian Stezura // E-Mail: sebastian_stezura@gmx.net         //
// GitHub: CplusplusB                                                     //
// version: 1 // date: 03.01.2022                                         //
////////////////////////////////////////////////////////////////////////////

//read data from Standards_eval_end_data_YEAR.txt (standards_eval_corr.cc) and Ambient_data_YEAR.txt (ambient.cc)
//draw Graphs
//write corrected data to file Ambient_data_YEAR_corr.txt

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compile command:                                                                                                          //
// g++-10 eval_air_std.cc tinyfiledialogs.c -ltbb `root-config --cflags --glibs --ldflags` -lMinuit -o ./Eval_air_std.o  //
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

using namespace std;
namespace fs = std::filesystem;


//class for data storage
class Data
{
public:
    vector<string> identifier, interval;//, interval;
    string analysis_o;
    vector<double> timed, timed_err, timed_dated, H2O_mean, H2O_sd, O18, O18_sd, H2, H2_sd, H2O_sl, H2O_sl_sd;
    vector<double> H2O_mean_mean, O18_mean, H2_mean, timed_mean, temp, timed_mean_conv, D_excess;
    vector<double> windvel, contemp, rh1, rh2, grad, apress, o3g1, o3g3, no, ventemp, winddir, prec;
    vector<int> inj_nmb, first;
    vector<TDatime> date, date_mean;
    double O18_corr, O18_corr_sd, H2_corr, H2_corr_sd;
    double timed_corr, timed_corr_all, timed_mean_conv_corr;
    int corr = 0;
    string ID_name;
    string file_name = "0";
    double amb_O18_max, amb_O18_min, amb_H2_max, amb_H2_min, amb_H2O_max, amb_H2O_min;
    double temp_max, temp_min;
    vector<int> month_all{1,2,3,4,5,6,7,8,9,10,11,12};
    vector<int> month;
    vector<string> month_names{"space", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"};
    double min, max;

    Data(){};
    Data(string ID){ID_name = ID;};
    void Destroy(){
        timed.clear();
        date.clear();
        timed_err.clear();
        H2O_mean.clear();
        H2O_sd.clear();
        O18.clear();
        O18_sd.clear();
        H2.clear();
        H2_sd.clear();
        H2O_sl.clear();
        H2O_sl_sd.clear();
        inj_nmb.clear();
        first.clear();
        //windvel.clear();
        //contemp.clear();
        //rh1.clear();
        //rh2.clear();
        //grad.clear();
        //apress.clear();
        //o3g1.clear();
        //o3g3.clear();
        //no.clear();
        ventemp.clear();
        //winddir.clear();
        prec.clear();
        date.clear();
        H2O_mean_mean.clear();
        O18_mean.clear();
        H2_mean.clear();
        timed_mean.clear();
        temp.clear();
    };
    ~Data(){};
};

//correct O18 Injection
double correctO18_inj(int injnmb)
{
    if(injnmb >10){return 0.000;};
    vector<double> corr;
    corr.push_back(0.174);
    corr.push_back(0.099);
    corr.push_back(0.061);
    corr.push_back(0.025);
    corr.push_back(0.000);
    corr.push_back(0.000);
    corr.push_back(0.000);
    corr.push_back(0.000);
    corr.push_back(0.000);
    corr.push_back(0.000);
    return (corr[injnmb-1]);
};

//correct H2 Injection
double correctH2_inj(int injnmb)
{
    if(injnmb >10){return 0.000;};
    vector<double> corr;
    corr.push_back(3.915);
    corr.push_back(1.463);
    corr.push_back(0.819);
    corr.push_back(0.504);
    corr.push_back(0.293);
    corr.push_back(0.222);
    corr.push_back(0.112);
    corr.push_back(0.000);
    corr.push_back(0.000);
    corr.push_back(0.000);
    return (corr[injnmb-1]);
};

//correct O18 humidity
double correctO18_hum(double H2O)
{
    if(H2O >= 10000.){return 0;};
    if(H2O == 0){return 0;};
    double d_0 = -2.33817;
    double A = -2570.5;
    return (d_0-A/10000.)-(d_0-A/H2O);
};

//correct H2 humidity
double correctH2_hum(double H2O)
{
    if(H2O >= 10000.){return 0;};
    if(H2O == 0){return 0;};
    double d_0 = -42.251;
    double A = -10964.2;
    return (d_0-A/10000.)-(d_0-A/H2O);
};

//get first Datapoint in file
int getLineRawData(string datapath)
{
    int RawNum;
    ifstream inFile(datapath);
    if (inFile.is_open())
    {
        string line;
        int i = 1;
        while (getline(inFile, line))
        {
            if (line == "Raw Data")
            {
                RawNum = i+3;
                cout << "Line with first Data is: " << RawNum << endl;
            };
            i++;
        };
    };
    inFile.close();
    return RawNum;
};

//getting Data from file
void getData_std(string datapath, vector<Data> &data, string &year)
{
    string time_r, analysis_r, port_r, identifier_r, ignore_r, inj_nmb_r, H2O_mean_r, H2O_sd_r, O18_r, O18_sd_r, H2_r, H2_sd_r, temp_r, CH4_r, H2O_sl_r, first_r;
    double timer, inj_nmbr, H2O_meanr, H2O_sdr, O18r, O18_sdr, H2r, H2_sdr, CH4r, tempr, H2O_slr;
    TDatime date_code;
    // int ignore;
    //loop over alle files
    int RawNum = getLineRawData(datapath);
    ifstream inFile(datapath);
    cout << "Reading file at " << datapath << " ..." << endl;

    // read signal values from file
	if (inFile.is_open())
	{
        int i = 1;
        int j = 0;
        string line;
        while (getline(inFile, line))
		{
            if (i < RawNum) {i++; continue;};
            //if (i >= RawNum) {cout << "Here starts data, i = " << i << endl;};
            stringstream stst(line);
            getline(stst,time_r,',');
            getline(stst,analysis_r,',');
            getline(stst,port_r,',');
            getline(stst,identifier_r,',');
            getline(stst,ignore_r,',');
            getline(stst,inj_nmb_r,',');
            getline(stst,H2O_mean_r,',');
            getline(stst,H2O_sd_r,',');
            getline(stst,O18_r, ',');
            getline(stst,O18_sd_r,',');
            getline(stst,H2_r,',');
            getline(stst,H2_sd_r,',');
            getline(stst,temp_r,',');
            getline(stst,CH4_r,',');
            getline(stst,H2O_sl_r,',');
            getline(stst,first_r,',');
            if(i == RawNum)
            {
                cout << " | " << time_r << " | " <<  analysis_r << " | " <<  port_r << " | " <<  identifier_r << " | " <<  ignore_r << " | " << inj_nmb_r << " | " <<  H2O_mean_r << " | " <<  H2O_sd_r << " | " <<  O18_r << " | " <<  O18_sd_r << " | " <<  H2_r << " | " <<  H2_sd_r << " | " <<  temp_r << " | " <<  CH4_r << " | " <<  H2O_sl_r << " | " << endl;
            };

            try
            {
                date_code.Set
                (
                    stoi(time_r.substr(0,4)), //year
                    stoi(time_r.substr(4,2)), //month
                    stoi(time_r.substr(6,2)), //day
                    stoi(time_r.substr(8,2)), //hour
                    stoi(time_r.substr(10,2)), //minute
                    stoi(time_r.substr(12,2)) //second
                );
            }
            catch (...)
            {
                cout << "Time Code: " << time_r << " at " << analysis_r << endl;
            }

            try
            {
                timer = stod(time_r);
                H2O_meanr = stod(H2O_mean_r);
                H2O_sdr = stod(H2O_sd_r);
                O18r = stod(O18_r);
                O18_sdr = stod(O18_sd_r);
                H2r = stod(H2_r);
                H2_sdr = stod(H2_sd_r);
                CH4r = stod(CH4_r);
                tempr = stod(temp_r);
                H2O_slr = stod(H2O_sl_r);
                inj_nmbr = stoi(inj_nmb_r);
            }
            catch (...)
            {
                cout << "Time Code: " << time_r << " at " << analysis_r << endl;
            }


            if(i == RawNum)
            {
                cout << " | " << time_r << " | " <<  analysis_r << " | " <<  port_r << " | " <<  identifier_r << " | " <<  ignore_r << " | " << inj_nmb_r << " | " <<  H2O_meanr << " | " <<  H2O_sdr << " | " <<  O18r << " | " <<  O18_sdr << " | " <<  H2r << " | " <<  H2_sdr << " | " <<  tempr << " | " <<  CH4r << " | " <<  H2O_slr << " | " << endl;
            };
            if (time_r.substr(0,4) != year){continue;};
            if(H2O_slr > 1.8 || H2O_slr < 1.5)
            {
                cout << fixed << "Slope not in range at " << timer << endl;
                continue;
            };

            if (j == 0 && data.size() < 1){cout << "|j war 0|"; data.push_back(Data()); data[j].analysis_o = analysis_r;};
            if (analysis_r != data[j].analysis_o){j++; data.push_back(Data());};
            data[j].ID_name = identifier_r;
            data[j].analysis_o = analysis_r;
            data[j].timed.push_back(timer);
            data[j].timed_err.push_back(0.);
            data[j].timed_mean_conv.push_back(date_code.Convert());
            data[j].inj_nmb.push_back(inj_nmbr);
            data[j].H2O_mean.push_back(H2O_meanr);
            data[j].H2O_sd.push_back(H2O_sdr);
            data[j].O18.push_back(O18r);
            data[j].O18_sd.push_back(O18_sdr);
            data[j].H2.push_back(H2r);
            data[j].H2_sd.push_back(H2_sdr);
            data[j].H2O_sl.push_back(H2O_slr);
            data[j].first.push_back(stoi(first_r));
            data[j].month.push_back(stoi(time_r.substr(4,2)));



            if (0)
            {
                cout << "At i = " << i << " " << j << " with " << data[j].analysis_o << " at " << data[j].timed[0] << " O18 " << data[j].O18[0] << endl;
            };
            i++;
		};
	};
    cout << "Size of data is: " << data.size() << endl;
    inFile.close();

};

//Getting Average and Standard dev of vector for Standards
void average_stdev(vector<double> &x, double &mean, double &stdev)
{
    double mean_t = 0.;
    double sum = 0.;
    for(int i = 4; i < x.size(); i++)
    {
        mean_t = mean_t + x[i];
    };
    mean = mean_t / 6.;
    for(int i = 4; i < x.size(); i++)
    {
        sum = sum + (x[i] - mean)*(x[i] - mean);
    };
    stdev = std::sqrt(sum / 5.);
};

//Get all Standard dates for Memory correction of Ambient Air, Averaging Standards
void getStd_corr(vector<Data> &data)
{
    vector<double> O18, H2;

    for(int i = 0; i < data.size(); i++)
    {
        data[i].timed_corr_all = data[i].timed.back();
        data[i].timed_mean_conv_corr = data[i].timed_mean_conv.back();
        cout << fixed << setprecision(3) << "Standard all: " << data[i].timed_corr_all << endl;
    };
    for(int i = 0; i < data.size(); i++)
    {
        if(data[i].first.front() == 1){continue;};
        if(data[i].O18.size() < 10 || data[i].O18.size() > 10){continue;};
        data[i].timed_corr = data[i].timed.back();
        for(int j = 0; j < data[i].O18.size(); j++)
        {
            O18.push_back(data[i].O18[j]+correctO18_inj(data[i].inj_nmb[j]));
            H2.push_back(data[i].H2[j]+correctH2_inj(data[i].inj_nmb[j]));
        };
        average_stdev(O18, data[i].O18_corr, data[i].O18_corr_sd);
        average_stdev(H2, data[i].H2_corr, data[i].H2_corr_sd);
        data[i].corr = 1;
        O18.clear();
        H2.clear();
        cout << "Standard " << i << ":" << data[i].timed_corr << "||" << data[i].O18_corr << "||" << data[i].H2_corr << endl;
    };
};

//getting Ambient Data from file and Correct
void getData_amb(string datapath, Data &data, vector<Data> &data_std, string &year)
{
    string time_r, port_r, H2O_mean_r, O18_r, H2_r;
    double timer, H2O_meanr, O18r, H2r;
    TDatime date_code;

    //loop over alle files
    int RawNum = getLineRawData(datapath);
    ifstream inFile(datapath);
    cout << "Reading file at " << datapath << " ..." << endl;

    // read signal values from file
	if (inFile.is_open())
	{
        int i = 1;
        int j = 0;
        string line;
        while (getline(inFile, line))
		{
            if (i < RawNum) {i++; continue;};

            stringstream stst(line);
            getline(stst,time_r,',');
            getline(stst,port_r,',');
            getline(stst,H2O_mean_r,',');
            getline(stst,O18_r, ',');
            getline(stst,H2_r,',');


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
                H2O_meanr = stod(H2O_mean_r);
                O18r = stod(O18_r);
                H2r = stod(H2_r);
                data.timed.push_back(timer);
                data.date.push_back(date_code);
                data.H2O_mean.push_back(H2O_meanr);
                data.O18.push_back(O18r+correctO18_hum(H2O_meanr));
                data.H2.push_back(H2r+correctH2_hum(H2O_meanr));
            };
            i++;
		};
	};
    inFile.close();

};

//X seconds averaging Data
void meanXminData(Data &data, Data &data_new)
{
    int i = 0;
    double mean_H2O, mean_H2, mean_O18;
    TDatime date_code;
    double timed;
    int avetime = 60; //Averaging time for data
    //string date;
    cout << "size of data.timed: " << data.timed.size() << endl;
    while (i < data.timed.size()-avetime)
    {
        mean_H2O = 0;
        mean_H2 = 0;
        mean_O18 = 0;
        if (i == 0) {cout << "Data " << data.H2O_mean[0] << " " << data.H2[0] << " " << data.O18[0] << endl;};

        for (int j = 0; j < avetime; j++)
        {
            mean_H2O = mean_H2O + data.H2O_mean[i+j];
            mean_H2 = mean_H2 + data.H2[i+j];
            mean_O18 = mean_O18 + data.O18[i+j];
            //cout << data.H2O_mean_mean[i] << endl;
            i++;
        };
        // cout << data.date[i-avetime/2].GetDate() << " at " << i << " mean " << mean_H2O/avetime << "|" << mean_H2/avetime << "|" << mean_O18/avetime << endl;
        mean_H2O = mean_H2O / avetime;
        mean_H2 = mean_H2 / avetime;
        mean_O18 = mean_O18 / avetime;
        date_code = data.date[i-avetime/2];
        timed = data.timed[i-avetime/2];

        data_new.H2O_mean_mean.push_back(mean_H2O);
        data_new.H2_mean.push_back(mean_H2);
        data_new.O18_mean.push_back(mean_O18);
        data_new.timed_mean.push_back(timed);
        data_new.date_mean.push_back(date_code);
        //cout << data.H2O_mean_mean[data.H2O_mean_mean.size()-1] << endl;
        //cout << data.timed_mean[data.timed_mean.size()-1] << endl;
        //cout << "size: H2O|timed" << data.H2O_mean_mean.size() << "|" << data.timed_mean.size() << endl;


    };
    cout << "Size of data averaged: " << data_new.timed_mean.size() << endl;
};

//Compare times for Correction
double compare_dates(double date_amb, vector<double> date_comp, vector<double> corr)
{
    if(date_amb < date_comp[0])
    {
        return corr[0];
    };
    if(date_amb >= date_comp[date_comp.size()-1])
    {
        return corr[corr.size()-1];
    };
    if(date_amb >= date_comp.front())
    {
        for(int i = 0; i < date_comp.size()-1; i++)
        {
            if(date_amb >= date_comp[i])
            {
                if(date_amb < date_comp[i+1])
                {
                    return corr[i];
                };
            };
        };
    };
    return 0.;
};

//Compare times for Memory
bool compare_dates_mem(double date_amb, vector<double> date_comp, double skip)
{
    for(int i = 0; i < date_comp.size(); i++)
    {
        if(date_amb >= date_comp[i] && date_amb < date_comp[i]+skip)
        {
            return 1;
        };
    };
    return 0;
};

//Memory Correction and Calibration of Ambient Air
void memcorr_amb(Data &data_amb_mean, Data &data_amb_corr, vector<Data> &data)
{
    vector<double> std_date;
    vector<double> std_date_corr, O18_corr_r, H2_corr_r;
    double skip = 300.;
    double O18_corr = 0.; //value of Standard measured
    double H2_corr = 0.; //value of Standard measured
    double O18_true = -8.65; //O18 Value of Standard
    double H2_true = -61.55; //H2 value of Standard

    for(int i = 0; i < data.size(); i++)
    {
        std_date.push_back(data[i].timed_corr_all);
        if(data[i].corr == 1)
        {
            std_date_corr.push_back(data[i].timed_corr);
            O18_corr_r.push_back(data[i].O18_corr);
            H2_corr_r.push_back(data[i].H2_corr);
        };
    };
    int l = 0;
    for(int i = 0; i < std_date_corr.size(); i++)
    {
        cout << "Correction vector true: " << std_date_corr[i] << "||" << O18_corr_r[i] << "||" << H2_corr_r[i] << endl;
        if(l == 0)
        {
            cout << "Correction vector true last: " << std_date_corr[std_date_corr.size()-1] << "||" << O18_corr_r[O18_corr_r.size()-1] << "||" << H2_corr_r[H2_corr_r.size()-1] << endl;
            l++;
        };

    };
    cout << "Correction vectors count: " << std_date_corr.size() << endl;
    for(int i = 0; i < data_amb_mean.timed_mean.size(); i++)
    {

        if(compare_dates_mem(data_amb_mean.timed_mean[i], std_date,skip)){cout << "Skipped" << endl; continue;};
        //cout << endl;
        data_amb_corr.timed_mean.push_back(data_amb_mean.timed_mean[i]);
        data_amb_corr.date_mean.push_back(data_amb_mean.date_mean[i]);
        data_amb_corr.H2O_mean.push_back(data_amb_mean.H2O_mean_mean[i]);
        data_amb_corr.O18_mean.push_back(data_amb_mean.O18_mean[i]);
        data_amb_corr.H2_mean.push_back(data_amb_mean.H2_mean[i]);
        data_amb_corr.month.push_back(data_amb_mean.date_mean[i].GetMonth());
        data_amb_corr.timed_mean_conv.push_back(data_amb_mean.date_mean[i].Convert());

    };

    for(int i = 0; i < data_amb_corr.timed_mean.size(); i++)
    {
        data_amb_corr.O18_mean[i] = data_amb_corr.O18_mean[i] - O18_corr;
        data_amb_corr.H2_mean[i] = data_amb_corr.H2_mean[i] - H2_corr;
        data_amb_corr.D_excess.push_back(data_amb_corr.H2_mean[i] - 8. * data_amb_corr.O18_mean[i]);
    };
    // Standard per Standard correction
    for(int i = 0; i < data_amb_corr.timed_mean.size(); i++)
    {
        O18_corr = compare_dates(data_amb_corr.timed_mean[i], std_date_corr, O18_corr_r) - O18_true;
        H2_corr = compare_dates(data_amb_corr.timed_mean[i], std_date_corr, H2_corr_r) - H2_true;
        data_amb_corr.O18_mean[i] = data_amb_corr.O18_mean[i] - O18_corr;
        data_amb_corr.H2_mean[i] = data_amb_corr.H2_mean[i] - H2_corr;
        data_amb_corr.D_excess.push_back(data_amb_corr.H2_mean[i] - 8. * data_amb_corr.O18_mean[i]);
    };
};

//draw Graphs Yearplots
void drawGraphYear(Data &data_amb, string evalpath, string year)
{
    string name_file_year = evalpath + "/End/" + year + "_ambient_year" + ".png";
    string name_file_LMWL = evalpath + "/End/" + year + "_ambient_LMWL" + ".png";

    string O18title = "{}^{18}O Graph Ambient Air " + year + ";Date in GMT [mont/day];#delta^{18}O";
    string H2title = "{}^{2}H Graph Ambient Air " + year + ";Date in GMT [mont/day];#delta^{2}H";
    string H2Otitle = "H_{2}O Graph Ambient Air " + year + ";Date in GMT [mont/day];H_{2}O[ppm]";
    string Dextitle = "D-excess Graph Ambient Air " + year + ";Date in GMT [mont/day];D_{excess}";
    string LMWLtitle = "Local Meteoric Water Line Graph Ambient Air " + year + ";#delta^{18}O;#delta^{2}H";


    int width = 4000;
    int height = 1500;

    TDatime date_first;
    TDatime date_last;
    date_first.Set(
    stoi(year), //year
    1, //month
    1, //day
    0, //hour
    0, //minute
    0 //second
    );
    date_last.Set(
    stoi(year), //year
    12, //month
    31, //day
    23, //hour
    59, //minute
    59 //second
    );

    TCanvas *cYear = new TCanvas("Year","Year",0,0,width,height*4);

    cYear->Divide(1,4,0,0);
    cYear->cd(1);
    cYear->GetFrame()->SetBorderSize(12);
    cYear->SetGrid();

    TGraph *grO18 = new TGraph(data_amb.timed_mean.size(),data_amb.timed_mean_conv.data(),data_amb.O18_mean.data());
    grO18->SetTitle(O18title.c_str());
    grO18->SetMarkerStyle(43);
    grO18->SetMarkerColor(kBlack);
    grO18->SetMarkerSize(2);
    grO18->SetLineWidth(3);
    grO18->SetLineColor(kBlue);
    grO18->GetXaxis()->SetTimeFormat("%m/%d");
    grO18->GetXaxis()->SetTimeOffset(0,"gmt");
    grO18->GetXaxis()->SetNdivisions(-13, kFALSE);
    grO18->GetXaxis()->SetRangeUser(date_first.Convert(),date_last.Convert());
    grO18->Draw("AP");
    cYear->Update();
    gPad->Update();
    cYear->Update();

    cYear->cd(2);

    TGraph *grH2 = new TGraph(data_amb.timed_mean.size(),data_amb.timed_mean_conv.data(),data_amb.H2_mean.data());
    grH2->SetTitle(H2title.c_str());
    grH2->SetMarkerStyle(43);
    grH2->SetMarkerColor(kGreen);
    grH2->SetMarkerSize(2);
    grH2->SetLineWidth(3);
    grH2->SetLineColor(kBlue);
    grH2->GetXaxis()->SetTimeFormat("%m/%d");
    grH2->GetXaxis()->SetTimeOffset(0,"gmt");
    grH2->GetXaxis()->SetNdivisions(-13, kFALSE);
    grH2->GetXaxis()->SetRangeUser(date_first.Convert(),date_last.Convert());
    grH2->Draw("AP");
    cYear->Update();
    gPad->Update();
    cYear->Update();

    cYear->cd(3);


    TGraph *grDex = new TGraph(data_amb.timed_mean.size(),data_amb.timed_mean_conv.data(),data_amb.D_excess.data());
    grDex->SetTitle(Dextitle.c_str());
    grDex->SetMarkerStyle(43);
    grDex->SetMarkerColor(kOrange);
    grDex->SetMarkerSize(2);
    grDex->SetLineWidth(3);
    grDex->SetLineColor(kBlue);
    grDex->GetXaxis()->SetTimeFormat("%m/%d");
    grDex->GetXaxis()->SetTimeOffset(0,"gmt");
    grDex->GetXaxis()->SetNdivisions(13, kFALSE);
    grDex->GetYaxis()->SetRangeUser(-13., 28.);
    grDex->GetXaxis()->SetRangeUser(date_first.Convert(),date_last.Convert());
    grDex->Draw("AP");
    cYear->Update();
    gPad->Update();
    cYear->Update();

    cYear->cd(4);
    TGraph *grH2O = new TGraph(data_amb.timed_mean.size(),data_amb.timed_mean_conv.data(),data_amb.H2O_mean.data());
    grH2O->SetTitle(H2Otitle.c_str());
    grH2O->SetMarkerStyle(43);
    grH2O->SetMarkerColor(kBlue);
    grH2O->SetMarkerSize(2);
    grH2O->SetLineWidth(3);
    grH2O->SetLineColor(kBlue);
    grH2O->GetXaxis()->SetTimeFormat("%m/%d");
    grH2O->GetXaxis()->SetTimeOffset(0,"gmt");
    grH2O->GetXaxis()->SetNdivisions(-13, kFALSE);
    grH2O->GetXaxis()->SetRangeUser(date_first.Convert(),date_last.Convert());
    grH2O->Draw("AP");
    cYear->Update();
    gPad->Update();
    cYear->Update();

    cYear->Print(name_file_year.c_str());
    cYear->Close();

    TCanvas *cLMWL = new TCanvas("LMWL","LMWL",0,0,width,height);
    cLMWL->cd();
    cLMWL->GetFrame()->SetBorderSize(12);
    cLMWL->SetGrid();

    TGraph *grLMWL = new TGraph(data_amb.O18_mean.size(),data_amb.O18_mean.data(),data_amb.H2_mean.data());
    grLMWL->SetTitle(LMWLtitle.c_str());
    grLMWL->SetMarkerStyle(43);
    grLMWL->SetMarkerColor(kBlue);
    grLMWL->SetMarkerSize(2);
    grLMWL->Draw("AP");
    cLMWL->Update();
    gPad->Update();
    cLMWL->Update();

    TF1 *linFit = new TF1("f1", "pol1", -1000., 100.);

    linFit->SetParameters(8.,10.);
    grLMWL->Fit(linFit);
    linFit = grLMWL->GetFunction("f1");

    cout << "Parameter Intercept (0): " << linFit->GetParameter(0) << "+-" << linFit->GetParError(0) << endl;
    cout << "Parameter Slope (1): " << linFit->GetParameter(1) << "+-" << linFit->GetParError(1) << endl;

    linFit->SetLineColor(kRed);
    linFit->SetLineStyle(4);
    linFit->SetLineWidth(12);

    cLMWL->Update();
    gPad->Update();
    cLMWL->Update();

    cLMWL->Print(name_file_LMWL.c_str());
    cLMWL->Close();

    data_amb.amb_O18_max = *std::max_element(data_amb.O18_mean.begin(), data_amb.O18_mean.end());
    data_amb.amb_O18_min = *std::min_element(data_amb.O18_mean.begin(), data_amb.O18_mean.end());
    data_amb.amb_H2_max = *std::max_element(data_amb.H2_mean.begin(), data_amb.H2_mean.end());
    data_amb.amb_H2_min = *std::min_element(data_amb.H2_mean.begin(), data_amb.H2_mean.end());
    data_amb.amb_H2O_max = *std::max_element(data_amb.H2O_mean.begin(), data_amb.H2O_mean.end());
    data_amb.amb_H2O_min = *std::min_element(data_amb.H2O_mean.begin(), data_amb.H2O_mean.end());

    cout << "Max O18: " << data_amb.amb_O18_max << endl;
    cout << "Min O18: " << data_amb.amb_O18_min << endl;
    cout << "Max H2: " << data_amb.amb_H2_max << endl;
    cout << "Min H2: " << data_amb.amb_H2_min << endl;
    cout << "Max H2O: " << data_amb.amb_H2O_max << endl;
    cout << "Min H2O: " << data_amb.amb_H2O_min << endl;

    int N_d = 0; int N_u = 0; int N_dd = 0;
    for(int i = 0; i < data_amb.H2O_mean.size(); i++)
    {
        if(data_amb.H2O_mean[i] <= 10000.){N_d++;};
        if(data_amb.H2O_mean[i] >= 25000.){N_u++;};
        if(data_amb.H2O_mean[i] <= 7000.){N_dd++;};
    };
    cout << "Upper count: " << N_u << "||Lower count: " << N_d << "||Deeplow count: " << N_dd << endl;
};

//draw standard graphs and histograms
void drawStdGraph(vector<Data> &data_std, string evalpath, string year)
{
    string O18TitleHist = "{}^{18}O Standards VE Ambient " + year + ";#delta^{18}O;Counts";
    string H2TitleHist = "{}^{2}H Standards VE Ambient " + year + ";#delta^{2}H;Counts";
    string O18TitleAll = "{}^{18}O Standards VE Ambient " + year + ";Date in GMT [mont/day];{}^{18}O";
    string H2TitleAll = "{}^{2}H Standards VE Ambient " + year + ";Date in GMT [mont/day];{}^{2}H";

    string name_file_hist = evalpath + "/End/" + year + "_histo_std.png";
    string name_file_all = evalpath + "/End/" + year + "_all_std.png";

    vector<double> O18, H2, count;
    int N = 1;
    for(int i = 0; i < data_std.size(); i++)
    {
        if(data_std[i].corr == 0){continue;};
        O18.push_back(data_std[i].O18_corr);
        H2.push_back(data_std[i].H2_corr);
        count.push_back(data_std[i].timed_mean_conv_corr); N++;
    };

    int width = 1200;
    int height = 1200;

    TCanvas *cHist = new TCanvas("hist","Histograms",0,0,width,height);
    cHist->Divide(1,2);
    cHist->cd(1);

    TH1D *hO18 = new TH1D("Histogram statistics {}^{18}O",O18TitleHist.c_str(),150,-4,-2);
    hO18->FillN(O18.size(),O18.data(),NULL);
    hO18->GetXaxis()->SetTitleSize(0.044);
    hO18->GetYaxis()->SetTitleSize(0.044);
    hO18->Draw();

    cHist->cd(2);
    TH1D *hH2 = new TH1D("Histogram statistics {}^{2}H",H2TitleHist.c_str(),150,-44,-47);
    hH2->FillN(H2.size(),H2.data(),NULL);
    hH2->GetXaxis()->SetTitleSize(0.044);
    hH2->GetYaxis()->SetTitleSize(0.044);
    hH2->Draw();

    cHist->Update();
    gPad->Update();
    cHist->Update();
    cHist->Print(name_file_hist.c_str());
    cHist->Close();

    TCanvas *cAll = new TCanvas("All","All",0,0,width,height);

    cAll->Divide(1,2,0,0);
    cAll->GetFrame()->SetBorderSize(12);
    cAll->SetGrid();
    cAll->cd(1);
    TGraph *grO18 = new TGraph(O18.size(),count.data(),O18.data());
    grO18->SetTitle(O18TitleAll.c_str());
    grO18->SetMarkerStyle(43);
    grO18->SetMarkerColor(kBlack);
    grO18->SetMarkerSize(2);
    grO18->SetLineWidth(3);
    grO18->SetLineColor(kBlue);
    grO18->GetXaxis()->SetTimeFormat("%m/%d");
    grO18->GetXaxis()->SetTimeOffset(0,"gmt");
    grO18->GetXaxis()->SetNdivisions(-13, kFALSE);
    grO18->Draw("APL");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->cd(2);
    TGraph *grH2 = new TGraph(H2.size(),count.data(),H2.data());
    grH2->SetTitle(H2TitleAll.c_str());
    grH2->SetMarkerStyle(43);
    grH2->SetMarkerColor(kBlack);
    grH2->SetMarkerSize(2);
    grH2->SetLineWidth(3);
    grH2->SetLineColor(kBlue);
    grH2->GetXaxis()->SetTimeFormat("%m/%d");
    grH2->GetXaxis()->SetTimeOffset(0,"gmt");
    grH2->GetXaxis()->SetNdivisions(-13, kFALSE);
    grH2->Draw("APL");
    cAll->Update();
    gPad->Update();
    cAll->Update();

    cAll->Print(name_file_all.c_str());
    cAll->Close();


};

//Write amb Data to file
void writeData(Data &data_amb, string evalpath, string year)
{
    string OutputFileName;
    OutputFileName = "Ambient_data_" + year + "_corr.txt";
    ofstream outFile (evalpath + "/End/" + OutputFileName);
    outFile << "Time," << "O18," << "H2," << "Dexcess," << "H2O" << endl;
    cout << "Writing Data amb corrected to: " << OutputFileName << endl;

    for (int j = 0; j < data_amb.timed_mean.size(); j++)
    {
        outFile << fixed << setprecision(3) << data_amb.timed_mean[j] << "," << data_amb.O18_mean[j] << "," << data_amb.H2_mean[j] << "," << data_amb.D_excess[j] << "," << data_amb.H2O_mean[j] << endl;
    };
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
    string datapath_amb = tinyfd_openFileDialog("Choose File with Ambient data", ".", 1, lFilterPatterns, NULL, 0);
    cout << "Data Ambient file: " << datapath_amb << endl;
    string datapath_std = tinyfd_openFileDialog("Choose File with Standard data", ".", 1, lFilterPatterns, NULL, 0);
    cout << "Data Standard file: " << datapath_std << endl;

    //choose directory for evaluation data
    string evalpath = tinyfd_selectFolderDialog("Choose Folder for Evaluation", ".");
    cout << "Evalpath " << evalpath << endl;

    //////////////////////////////////////////////////
    // Read files and store values //
    //////////////////////////////////////////////////
    Data data_amb, data_amb_mean, data_amb_corr;
    vector<Data> data_std;
    cout << "################" << endl;
    cout << "Reading Standard data ...." << endl;
    getData_std(datapath_std, data_std, year);
    cout << "Finished." << endl;
    cout << "Averaging and Correcting Standard data ...." << endl;
    getStd_corr(data_std);
    cout << "Size of Std Data: " << data_std.size() << endl;
    cout << "Finished." << endl << "################" << endl << "Reading Ambient Air data ..." << endl;
    getData_amb(datapath_amb, data_amb, data_std, year);

    cout << "Finieshed." << endl << "################" << endl << "Averaging Ambient Air" << endl;
    meanXminData(data_amb, data_amb_mean);
    cout << "Finished." << endl;
    data_amb.Destroy();
    cout << "Old Ambient Air destroyed." << endl;
    cout << "################" << endl << "Clearing Ambient Air from Memory..." << endl;
    memcorr_amb(data_amb_mean, data_amb_corr, data_std);
    cout << "Finished." << endl;
    cout << "Size of Data corrected: " << data_amb_corr.timed_mean.size() << endl;
    data_amb_mean.Destroy();

    drawGraphYear(data_amb_corr, evalpath, year);
    drawStdGraph(data_std, evalpath, year);
    writeData(data_amb_corr, evalpath, year);




    return 0;
}