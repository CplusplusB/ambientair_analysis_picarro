////////////////////////////////////////////////////////////////////////////
// Programm for evaluating ambient air data from the Picarro L-2130i      //
// author: Sebastian Stezura // E-Mail: sebastian_stezura@gmx.net         //
// GitHub: CplusplusB                                                     //
// version: 1 // date: 03.01.2022                                         //
////////////////////////////////////////////////////////////////////////////

// Get Data from .csv Files in Folder. Data determined by names.cc output file or choose own file
// write to file "Ambient_data_YEAR.txt"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compile command for linux:                                                                                            //
// g++-10 ambient.cc tinyfiledialogs.c -ltbb `root-config --cflags --glibs --ldflags` -lMinuit -o ./Ambient.o  //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
#include <TLinearFitter.h> //include linearFitter
#include <TDatime.h> //include Date storage

using namespace std;
namespace fs = std::filesystem;


//class for data storage
class Data
{
public:
    vector<string> port, timed, H2O_mean, O18, H2;
    vector<TDatime> date;
    string ID_name;
    string file_name = "0";

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

//check date
bool int_compare(int i, int j)
{
    return (i < j);
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
void getData(string name, string files_adress, Data &data)
{
    string line_r, analysis_r, time_code_r, port_r, inj_nmb_r, O18w_r, H2w_r, H2Ow_mean_r, ignore_r, good_r, O18v_r, H2v_r, H2Ov_mean_r, identifier1_r, identifier2_r, gas_conf_r, time_mean_r, O18_sd_r, H2_sd_r, H2O_sd_r, O18_sl_r, H2_sl_r, H2O_sl_r, base_shift_r, slope_r, res_r, base_curv_r, interval_r, CH4_r, H2O_adj_r, H2O_shift_r, n2_r, temp_r, tray_r, sample_r, job_r, method_r, error_r;
    TDatime date_code;
    double O18r,H2r;
    string last_analysis = "none";
    int memory = 0;
    //loop over alle files
    ifstream inFile (files_adress);
    cout << "Reading file " << name << " ..." << endl;
    data.file_name = name;
    // read signal values from file
	if (inFile.is_open())
	{
        int i = 0;
        string line;
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
                if(port_r == "Ambient")
                {
                    O18r = stod(O18v_r);
                    H2r = stod(H2v_r);
                };
            }
            catch (...)
            {
                cout << "Time Code: " << time_code_r << " at " << time_mean_r << endl;
            }
            //memory correction
            if (port_r == "Ambient" && last_analysis == "H2O")
            {
                memory++;
                if(memory <= 180){continue;};
            };
            if(memory >= 180)
            {
                memory = 0;
            };

            if (port_r == "Ambient")
            {
                data.port.push_back(port_r);
                data.timed.push_back(time_code_r);
                data.H2O_mean.push_back(H2Ov_mean_r);
                data.O18.push_back(O18v_r);
                data.H2.push_back(H2v_r);
            };
            last_analysis = gas_conf_r;
		};
	};
    //verbose identifiers
    // for (int i = 0; i < data.size(); i++)
    // {
    //     cout << data[i].ID_name << endl;
    // };
    inFile.close();

};

//write evaluated Data
void writeData(string year, string evalpath, vector<string> files_name, vector<Data> &data)
{
    time_t t = time(0);
    string c_time = ctime(&t);
    string curr_time = c_time.substr(8,2) + c_time.substr(4,3) + c_time.substr(20,4) + "_" + c_time.substr(11,2) + c_time.substr(14,2) + c_time.substr(17,2);
    cout << "current time: " << curr_time << endl;
    string OutputFileName;
    OutputFileName = "Ambient_data_" + year + ".txt";
    ofstream outFile (evalpath + "/" + OutputFileName);
    outFile << "Evaluation from: " << c_time << endl;
    outFile << "Type: Ambient" << endl << endl;
    outFile << "Files evaluated: " << endl;
    for (int i = 0; i < files_name.size(); i++)
    {
        outFile << files_name[i] << endl;
    };
    outFile << endl << "========================" << endl << "Raw Data" << endl << "========================" << endl;
    outFile << "Time," << "Port," << "H2O mean," << "O18," << "H2" << endl;

    for (int i = 0; i < data.size(); i++)
    {
        for (int j = 0; j < data[i].timed.size(); j++)
        {
            outFile << data[i].timed[j] << "," << data[i].port[j] << "," << data[i].H2O_mean[j] << "," << data[i].O18[j] << "," << data[i].H2[j] << endl;
        };
    };
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

    getFiles(files_name, files_date, files_adress);
    //choose directory for evaluation data
    string evalpath = tinyfd_selectFolderDialog("Choose Folder for Evaluation", ".");
    cout << "Evalpath " << evalpath << endl;
    cout << "############" << endl;
    cout << "Which year?" << endl;
    cout << "############" << endl;
    cin >> year;
    //cout << year << endl; //debug

    //////////////////////////////////////////////////
    // Read csv file(s) and store values in vectors //
    //////////////////////////////////////////////////
    vector<Data> data, data_sorted;

    vector<int> parItr, endposition;

    int reservedN = 200; //How much files are there
    for (int i = 0; i < reservedN; i++){ data.push_back(Data()); };

    for (int i = 0; i < files_adress.size(); i++){ parItr.push_back(i); };

    //loop over alle files and read signal
    auto files_loop = [files_name, files_adress, &data](int i)
    {
        cout << i << ". ";
        getData(files_name[i], files_adress[i], data[i]);
        data[i].file_name = files_name[i];
        //cout << data[i].file_name.substr(18,8) << data[i].file_name.substr(27,6) << endl;
    };
    for_each(execution::par,parItr.begin(),parItr.end(),files_loop);

    //sort Data to date
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
    data.clear();
    writeData(year, evalpath, files_name, data_sorted);

    return 0;
}

