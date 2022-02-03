////////////////////////////////////////////////////////////////////////////
// Programm for evaluating ambient air data from the Picarro #insertModel //
// author: Sebastian Stezura // E-Mail: sebastian_stezura@gmx.net         //
// GitHub: CplusplusB                                                     //
// version: 1 // date: 03.01.2022                                         //
////////////////////////////////////////////////////////////////////////////

//reading csv files of the folder
//Draw Graphs to visualize them
//sort data to date
//write data to file Standards_eval_end_data_YEAR.txt

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compile command:                                                                                                                    //
// g++-10 standards_eval_corr.cc tinyfiledialogs.c -ltbb `root-config --cflags --glibs --ldflags` -lMinuit -o ./Standards_eval_corr.o  //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

using namespace std;
namespace fs = std::filesystem;


//class for data storage
class Data
{
public:
    vector<string> analysis, port, identifier, inj_nmb, timed, ignore, H2O_mean, H2O_sd, O18, O18_sd, H2, H2_sd, CH4, temp, H2O_sl, H2O_sl_sd, first;
    vector<TDatime> date;
    string ID_name;
    string file_name = "0";
    double min, max;

    Data(){};
    Data(string ID){ID_name = ID;};
    ~Data(){};
};

//getting files for evaluation
void getFiles(vector<string> &names, vector<string> &dates, vector<string> &adresses)
{
    string afilepath = tinyfd_selectFolderDialog("Choose Folder with csv Files", ".");

    string filesubstr;
    fs::path pafilepath(afilepath);
    fs::directory_iterator end_itr;
    for (fs::directory_iterator itr(pafilepath); itr != end_itr; ++itr)
    {
        if (fs::is_regular_file(itr->path()))
        {
            filesubstr = itr->path().filename().string();
            filesubstr = filesubstr.substr(filesubstr.size()-4,4);
            if (filesubstr == ".csv")
            {
                names.push_back(itr->path().filename().string());
                dates.push_back(itr->path().filename().string().substr(itr->path().filename().string().size()-19, 8));
                adresses.push_back(itr->path());
            };
        };
    };
    cout << "===================" << endl;
    cout << "Reading files" << endl;
    cout << "===================" << endl;

    for (int i = 0; i < names.size(); i++)
    {
        cout << i << ". " << names[i] << " from " << dates[i] << " at " << adresses[i] << endl;
    };
};

//check identifier
bool checkID_same(vector<string> identifier, string data_ID)
{
    bool ID_same = 0;
    for (int i = 0; i < identifier.size(); i++)
    {
        if (identifier[i] == data_ID){ID_same = 1;};
    };
    return ID_same;
};

//get identifers
void getID(vector<string> &identifiers)
{
    string afilepath = tinyfd_openFileDialog("Choose File with ID names", ".", 0, NULL, NULL, 0);
    ifstream inFile (afilepath);
    cout << "IDs:" << endl;
    if (inFile.is_open())
	{
        string line;
        while (getline(inFile, line))
        {
            remove(line.begin(), line.end(), ' ');
            identifiers.push_back(line);
            cout << line << endl;
        };
    };
    inFile.close();
};

//getting Data from files
void getData(string name, string files_adress, Data &data, vector<string> ID_names)
{
    string line_r, analysis_r, time_code_r, port_r, inj_nmb_r, O18w_r, H2w_r, H2Ow_mean_r, ignore_r, good_r, O18v_r, H2v_r, H2Ov_mean_r, identifier1_r, identifier2_r, gas_conf_r, time_mean_r, O18_sd_r, H2_sd_r, H2O_sd_r, O18_sl_r, H2_sl_r, H2O_sl_r, base_shift_r, slope_r, res_r, base_curv_r, interval_r, CH4_r, H2O_adj_r, H2O_shift_r, n2_r, temp_r, tray_r, sample_r, job_r, method_r, error_r;
    TDatime date_code;
    //loop over alle files
    ifstream inFile (files_adress);
    cout << "Reading file " << name << " ..." << endl;
    data.file_name = name;

    double O18, H2O_sl;
    // read signal values from file
	if (inFile.is_open())
	{
        int i = 0;
        string line;
        bool first_true = 0;
        string last_analysis = "0";
        string last_analysis_true = "0";

        while (getline(inFile, line))
		{
            remove(line.begin(), line.end(), ' ');
            stringstream stst(line);
            getline(stst,line_r,',');
            getline(stst,analysis_r,',');
            getline(stst,time_code_r,',');
            getline(stst,port_r,',');
            getline(stst,inj_nmb_r,',');
            getline(stst,O18w_r,',');
            getline(stst,H2w_r,',');
            getline(stst,H2Ow_mean_r,',');
            getline(stst,ignore_r,',');
            getline(stst,good_r,',');
            getline(stst,O18v_r,',');
            getline(stst,H2v_r,',');
            getline(stst,H2Ov_mean_r,',');
            getline(stst,identifier1_r,',');
            getline(stst,identifier2_r,',');
            getline(stst,gas_conf_r,',');
            getline(stst,time_mean_r,',');
            getline(stst,O18_sd_r,',');
            getline(stst,H2_sd_r,',');
            getline(stst,H2O_sd_r,',');
            getline(stst,O18_sl_r,',');
            getline(stst,H2_sl_r,',');
            getline(stst,H2O_sl_r,',');
            getline(stst,base_shift_r,',');
            getline(stst,slope_r,',');
            getline(stst,res_r,',');
            getline(stst,base_curv_r,',');
            getline(stst,interval_r,',');
            getline(stst,CH4_r,',');
            getline(stst,H2O_adj_r,',');
            getline(stst,H2O_shift_r,',');
            getline(stst,n2_r,',');
            getline(stst,temp_r,',');
            getline(stst,tray_r,',');
            getline(stst,sample_r,',');
            getline(stst,job_r,',');
            getline(stst,method_r,',');
            getline(stst,error_r,',');

            //format time format
            remove(time_code_r.begin(), time_code_r.end(), '/');
            remove(time_code_r.begin(), time_code_r.end(), ':');
            for (int j = 0; j < 4; j++) {if (!time_code_r.empty()) {time_code_r.pop_back();};};
            //time_code_r.insert(10,1,' ');

            if (i == 0) {i = 1; continue;};
            try
            {
                date_code.Set
                (
                    stoi(time_code_r.substr(0,4)),
                    stoi(time_code_r.substr(4,2)),
                    stoi(time_code_r.substr(6,2)),
                    stoi(time_code_r.substr(8,2)),
                    stoi(time_code_r.substr(10,2)),
                    stoi(time_code_r.substr(12,2))
                );
                if(O18w_r == ""){O18w_r = "0.";};
                O18 = stod(O18w_r);
                H2O_sl = stod(H2O_sl_r);
            }
            catch (...)
            {
                cout << "Time Code: " << time_code_r << " at " << time_mean_r << "||" << H2O_sl_r << endl;
            }

            if (time_code_r.substr(0,4) == "2018"){data.min = -3.9; data.max = -2.6;};
            if (time_code_r.substr(0,4) == "2019"){data.min = -3.9; data.max = -2.4;};
            if (time_code_r.substr(0,4) == "2020"){data.min = -3.2; data.max = -2.2;};
            if (time_code_r.substr(0,4) == "2021"){data.min = -3.2; data.max = -2.2;};

            if (checkID_same(ID_names, identifier2_r) && identifier2_r != "" && O18 >= data.min && O18 <= data.max)// && O18w_r != "")
            {
                //if(O18r < -20. && O18r > 100.){continue;};
                if(last_analysis_true == analysis_r)
                {
                    first_true = 1;
                }
                else
                {
                    first_true = 0;
                };
                if(line_r == "1" || line_r == "2")
                {
                    first_true = 1;
                    last_analysis_true = analysis_r;
                    cout << "First data " << line_r << " at " << time_code_r << endl;
                };

                if(first_true == 1)
                {
                    data.first.push_back("1");
                }
                else
                {
                    data.first.push_back("0");
                };
                data.port.push_back(port_r);
                data.analysis.push_back(analysis_r);
                data.timed.push_back(time_code_r);
                data.identifier.push_back(identifier2_r);
                //lines without data set to specific values to display them
                if(O18w_r == "")
                {
                    data.O18.push_back("0.");
                    data.H2O_mean.push_back("0");
                    data.ignore.push_back("-2");
                    data.H2O_sd.push_back("0.");
                    data.O18_sd.push_back("0.");
                    data.H2.push_back("0.");
                    data.H2_sd.push_back("0.");
                    data.temp.push_back("0.");
                    data.CH4.push_back("0.");
                    data.H2O_sl.push_back("0.");
                    cout << "Set 0 at " << analysis_r << " with ID " << identifier2_r << endl;

                }
                else
                {
                    data.ignore.push_back(ignore_r);
                    data.H2O_mean.push_back(H2Ow_mean_r);
                    data.H2O_sd.push_back(H2O_sd_r);
                    data.O18.push_back(O18w_r);
                    data.O18_sd.push_back(O18_sd_r);
                    data.H2.push_back(H2w_r);
                    data.H2_sd.push_back(H2_sd_r);
                    data.temp.push_back(temp_r);
                    data.CH4.push_back(CH4_r);
                    data.H2O_sl.push_back(H2O_sl_r);
                };

                data.inj_nmb.push_back(inj_nmb_r);
                data.date.push_back(date_code);
            };
		};
	};
    inFile.close();

};

//get unique identifiers
void getUnique(vector<string> &unique, Data &data)
{
    for (int i = 0; i < data.port.size(); i++)
    {
        if (!checkID_same(unique, data.identifier[i]))
        {
            unique.push_back(data.identifier[i]);
            cout << "Found dataname: " << data.identifier[i] << " in " << data.file_name << endl;
        };

    };

};

//check date
bool int_compare(int i, int j)
{
    return (i < j);
};

//write evaluated Data
void writeData(string year, string evalpath, vector<string> files_name, vector<Data> &data)
{
    time_t t = time(0);
    string c_time = ctime(&t);
    double value = 0;
    string curr_time = c_time.substr(8,2) + c_time.substr(4,3) + c_time.substr(20,4) + "_" + c_time.substr(11,2) + c_time.substr(14,2) + c_time.substr(17,2);
    cout << "current time: " << curr_time << endl;
    string OutputFileName;
    OutputFileName = "Standards_eval_end_data_corr" + year + ".txt";
    ofstream outFile (evalpath + "/" + OutputFileName);
    outFile << "Evaluation from: " << c_time << endl;
    outFile << "Type: Standards" << endl << endl;
    outFile << "Files evaluated: " << endl;
    for (int i = 0; i < files_name.size(); i++)
    {
        outFile << files_name[i] << endl;
    };
    outFile << endl << "========================" << endl << "Raw Data" << endl << "========================" << endl;
    outFile << "Time," << "Analaysis No," << "Port,"<< "Identifier," << "Ignore," << "Inj_nmb," << "H2O mean," << "H2O sd," << "O18," << "O18 sd," << "H2," << "H2 sd," << "temperature," << "CH4," << "H2O sl," << "first" << endl;

    for (int i = 0; i < data.size(); i++)
    {
        for (int j = 0; j < data[i].timed.size(); j++)
        {
            outFile << data[i].timed[j] << "," << data[i].analysis[j] << "," << data[i].port[j] << "," << data[i].identifier[j] << "," << data[i].ignore[j] << "," << data[i].inj_nmb[j] << "," << data[i].H2O_mean[j] << "," << data[i].H2O_sd[j] << "," << data[i].O18[j] << "," << data[i].O18_sd[j] << "," << data[i].H2[j] << "," << data[i].H2_sd[j] << "," << data[i].temp[j] << "," << data[i].CH4[j] << "," << data[i].H2O_sl[j] << "," << data[i].first[j] << endl;
        };
    };
};

//Draw Graph
void drawGraph(Data &data_std, string evalpath, string year)
{
    string name_file_o18 = evalpath + "/first/" + year + "_O18_" + data_std.analysis[0] + ".png";
    string name_file_h2 = evalpath + "/first/" + year + "_H2_" + data_std.analysis[0] + ".png";

    string O18title = "O^{18} Graph " + data_std.analysis[0] + " " + year + ";Time;#delta{}^{18}O";
    string H2title = "H^{2} Graph " + data_std.analysis[0] + " " + year + ";Time;#delta{}^{2}H";

    vector<vector<double>> x_o,y_o,x_o_err,y_o_err;
    vector<vector<double>> x_h,y_h,x_h_err,y_h_err;

    vector<double> x,y,x_err,y_err;

    int width = 8000;
    int height = 3000;

    string now_analysis;
    string last_analysis = "0";

    TCanvas *cO18 = new TCanvas("O^{18}","O^{18}",0,0,width,height);

    cO18->cd();
    cO18->GetFrame()->SetBorderSize(12);
    cO18->SetGrid();

    TMultiGraph *mgrO18 = new TMultiGraph();
    mgrO18->SetTitle(O18title.c_str());

    for (int i = 0; i < data_std.timed.size(); i++)
    {
        now_analysis = data_std.analysis[i];
        if(now_analysis != last_analysis)
        {
            if(data_std.first[i] == "1")
            {
                x.push_back(stod(data_std.inj_nmb[i]));
                y.push_back(stod(data_std.O18[i]));
                x_err.push_back(0.);
                // y_err.push_back(stod(data_std.O18_sd[i]));
                y_err.push_back(0.);
            }
            else
            {
                x_o.push_back(vector<double>());
                y_o.push_back(vector<double>());
                x_o_err.push_back(vector<double>());
                y_o_err.push_back(vector<double>());
                x_o.back().push_back(stod(data_std.inj_nmb[i])+x_o.size()/20.);
                y_o.back().push_back(stod(data_std.O18[i]));
                x_o_err.back().push_back(0.);
                y_o_err.back().push_back(0.);
                // y_o_err.back().push_back(stod(data_std.O18_sd[i]));
            };
        }
        else if(now_analysis == last_analysis)
        {
            if(data_std.first[i] == "1")
            {
                x.push_back(stod(data_std.inj_nmb[i]));
                y.push_back(stod(data_std.O18[i]));
                x_err.push_back(0.);
                //y_err.push_back(stod(data_std.O18_sd[i]));
                y_err.push_back(0.);
            }
            else
            {
                x_o.back().push_back(stod(data_std.inj_nmb[i])+x_o.size()/20.);
                y_o.back().push_back(stod(data_std.O18[i]));
                x_o_err.back().push_back(0.);
                y_o_err.back().push_back(0.);
                //y_o_err.back().push_back(stod(data_std.O18_sd[i]));
            };
        }
        else
        {
            cout << "Conditions broke at " << data_std.timed[i] << " ..." << endl;
        };
        last_analysis = now_analysis;
    };

    vector<TGraphErrors*> grO18;
    for(int i = 0; i < x_o.size(); i++)
    {
        grO18.push_back(new TGraphErrors(x_o[i].size(),x_o[i].data(),y_o[i].data(),x_o_err[i].data(),y_o_err[i].data()));
        grO18[i]->SetMarkerStyle(43);
        grO18[i]->SetMarkerColor(kRed);
        grO18[i]->SetMarkerSize(4);
        grO18[i]->SetLineWidth(3);
        grO18[i]->SetLineColor(kRed);
        mgrO18->Add(grO18[i],"LP");
    };

    TGraphErrors *grO18_first = new TGraphErrors(x.size(),x.data(),y.data(),x_err.data(),y_err.data());
    grO18_first->SetMarkerStyle(43);
    grO18_first->SetMarkerColor(kBlue);
    grO18_first->SetMarkerSize(5);
    grO18_first->SetLineWidth(3);
    grO18_first->SetLineColor(kBlue);
    mgrO18->Add(grO18_first,"LP");

    mgrO18->Draw("A");
    cO18->Update();
    gPad->Update();
    cO18->Update();

    cO18->Print(name_file_o18.c_str());
    cO18->Close();


    x.clear();
    y.clear();

    x_err.clear();
    y_err.clear();

    ///////////
    //H2 graphs
    ///////////
    last_analysis = "0";

    TCanvas *cH2 = new TCanvas("H^{2}","H^{2}",0,0,width,height);

    cH2->cd();
    cH2->GetFrame()->SetBorderSize(12);
    cH2->SetGrid();

    TMultiGraph *mgrH2 = new TMultiGraph();
    mgrH2->SetTitle(H2title.c_str());

    cout << "Sending Data to vector ..." << data_std.analysis[0] << endl;

    for (int i = 0; i < data_std.timed.size(); i++)
    {
        now_analysis = data_std.analysis[i];
        if(now_analysis != last_analysis)
        {
            if(data_std.first[i] == "1")
            {
                x.push_back(stod(data_std.inj_nmb[i]));
                y.push_back(stod(data_std.H2[i]));
                x_err.push_back(0.);
                y_err.push_back(0.);
                // y_err.push_back(stod(data_std.H2_sd[i]));
            }
            else
            {
                x_h.push_back(vector<double>());
                y_h.push_back(vector<double>());
                x_h_err.push_back(vector<double>());
                y_h_err.push_back(vector<double>());
                x_h.back().push_back(stod(data_std.inj_nmb[i])+x_h.size()/20.);
                y_h.back().push_back(stod(data_std.H2[i]));
                x_h_err.back().push_back(0.);
                // y_h_err.back().push_back(stod(data_std.H2_sd[i]));
                y_h_err.back().push_back(0.);
            };
        }
        else if(now_analysis == last_analysis)
        {
            if(data_std.first[i] == "1")
            {
                x.push_back(stod(data_std.inj_nmb[i]));
                y.push_back(stod(data_std.H2[i]));
                x_err.push_back(0.);
                y_err.push_back(0.);
                // y_err.push_back(stod(data_std.H2_sd[i]));
            }
            else
            {
                x_h.back().push_back(stod(data_std.inj_nmb[i])+x_h.size()/20.);
                y_h.back().push_back(stod(data_std.H2[i]));
                x_h_err.back().push_back(0.);
                // y_o_err.back().push_back(stod(data_std.H2_sd[i]));
                y_h_err.back().push_back(0.);
            };
        }
        else
        {
            cout << "Conditions broke at " << data_std.timed[i] << " ..." << endl;
        };
        last_analysis = now_analysis;
    };

    cout << "Drawing Graph H2 not first ..." << endl;
    vector<TGraphErrors*> grH2;
    for(int i = 0; i < x_o.size(); i++)
    {

        grH2.push_back(new TGraphErrors(x_h[i].size(),x_h[i].data(),y_h[i].data(),x_h_err[i].data(),y_h_err[i].data()));
        grH2[i]->SetMarkerStyle(43);
        grH2[i]->SetMarkerColor(kRed);
        grH2[i]->SetMarkerSize(4);
        grH2[i]->SetLineWidth(3);
        grH2[i]->SetLineColor(kRed);
        mgrH2->Add(grH2[i],"LP");
    };
    cout << "Drawing Graph H2 first ..." << endl;

    TGraphErrors *grH2_first = new TGraphErrors(x.size(),x.data(),y.data(),x_err.data(),y_err.data());
    grH2_first->SetMarkerStyle(43);
    grH2_first->SetMarkerColor(kBlue);
    grH2_first->SetMarkerSize(5);
    grH2_first->SetLineWidth(3);
    grH2_first->SetLineColor(kBlue);
    mgrH2->Add(grH2_first,"PL");

    mgrH2->Draw("A");
    cH2->Update();
    gPad->Update();
    cH2->Update();

    cH2->Print(name_file_h2.c_str());
    cH2->Close();
};

int main(int argc, char* argv[])
{
    //////////////////////////////////
    // Name, date and path to files //
    //////////////////////////////////
    vector<string> files_name;
    vector<string> files_date;
    vector<string> files_adress;
    string year;

    cout << "Which year?" << endl;
    cin >> year;

    getFiles(files_name, files_date, files_adress);
    //choose directory for evaluation data
    string evalpath = tinyfd_selectFolderDialog("Choose Folder for Evaluation", ".");
    cout << "Evalpath " << evalpath << endl;
    //cout << year << endl; //debug

    //////////////////////////////////////////////////
    // Read csv file(s) and store values in vectors //
    //////////////////////////////////////////////////
    vector<Data> data, data_sorted;
    vector<string> unique_ID;
    getID(unique_ID);

    vector<int> parItr;
    int reservedN = 200; //How much files are there
    for (int i = 0; i < reservedN; i++){ data.push_back(Data()); };

    for (int i = 0; i < files_adress.size(); i++){ parItr.push_back(i); };

    //loop over alle files and read signal
    auto files_loop = [files_name, files_adress, &data, &unique_ID](int i)
    {
        cout << i << ". ";
        getData(files_name[i], files_adress[i], data[i], unique_ID);
        data[i].file_name = files_name[i];
    };
    for_each(execution::par,parItr.begin(),parItr.end(), files_loop);

    //////////////
    //sort Data //
    //////////////

    string date_name;
    long int date_int;
    int year_akt, yearr;
    vector<long int> date;

    for (int i = 0; i < data.size(); i++)
    {
        if (data[i].file_name == "0"){continue;};
        date_name = data[i].file_name.substr(18,8) + data[i].file_name.substr(27,6);
        try
        {
            date_int = stol(date_name);
        }
        catch (...)
        {
            cout << "aborted because of '" << date_name << "'" << endl;
        }
        date.push_back(date_int);
    };

    cout << "sorting..." << endl;
    cout << "##########" << endl;
    std::sort(date.begin(), date.end(), int_compare);
    for (int i = date.size()-1; i != 0; i--)
    {
        date_name = to_string(date[i]);
        year_akt = stoi(date_name.substr(0,4));
        yearr = stoi(year);
        if (year_akt < yearr)
        {
            date.insert(date.begin(), date[i]);
            date.pop_back();
        };
    };

    for (int i = 0; i < date.size(); i++)
    {
        cout << date[i] << endl;
    };
    cout << endl;
    cout << "########" << endl;

    for (int i = 0; i < date.size(); i++)
    {
        for (int j = 0; j < data.size(); j++)
        {
            if (data[j].file_name == "0"){continue;};
            date_name = data[j].file_name.substr(18,8) + data[j].file_name.substr(27,6);
            //cout << "Date name: " << date_name << endl;
            date_int = stol(date_name);
            //cout << "Date in int: " << date_int << endl;
            if (date[i] == date_int) {data_sorted.push_back(data[j]); cout << "sorting in " << date_int << endl;};
        };
    };

    ///////////////
    //Write to file
    ///////////////
    cout << "Writing data to file ..." << endl;

    writeData(year, evalpath, files_name, data_sorted);

    ///////////////
    //Draw Graphs
    ///////////////

    //loop over alle files and read signal
    cout << "Drawing graphs ..." << endl;
    for(int i = 0; i < data.size(); i++)
    {
        if(data[i].O18.size() < 1){continue;};
        drawGraph(data[i],evalpath,year);
    };





    return 0;
}