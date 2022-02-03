////////////////////////////////////////////////////////////////////////////
// Programm for evaluating ambient air data from the Picarro #insertModel //
// author: Sebastian Stezura // E-Mail: sebastian_stezura@gmx.net         //
// GitHub: CplusplusB                                                     //
// version: 1 // date: 03.01.2022                                         //
////////////////////////////////////////////////////////////////////////////

// Get ID names in .csv Files in Folder. Write ID names in Standards_eval_names_YEAR_CURR_TIME.txt

////////////////////////////////////////////////////////////////////////////
// compile command:                                                       //
// g++-10 names.cc tinyfiledialogs.c -std=c++17 -o GetNames.o             //
////////////////////////////////////////////////////////////////////////////

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

using namespace std;
namespace fs = std::filesystem;


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

//getting Data from files
void getNames(string name, string files_adress, vector<string> &unique_ID)
{
    string line_r, analysis_r, time_code_r, port_r, inj_nmb_r, O18w_r, H2w_r, H2Ow_mean_r, ignore_r, good_r, O18v_r, H2v_r, H2Ov_mean_r, identifier1_r, identifier2_r, gas_conf_r, time_mean_r, O18_sd_r, H2_sd_r, H2O_sd_r, O18_sl_r, H2_sl_r, H2O_sl_r, base_shift_r, slope_r, res_r, base_curv_r, interval_r, CH4_r, H2O_adj_r, H2O_shift_r, n2_r, temp_r, tray_r, sample_r, job_r, method_r, error_r;
    //loop over alle files
    ifstream inFile (files_adress);
    cout << "Reading file " << name << " ..." << endl;
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
            if (!checkID_same(unique_ID, identifier2_r) && identifier2_r != "")
            {
                unique_ID.push_back(identifier2_r);
                cout << identifier2_r << endl;
            };
		};
	};
    inFile.close();

};

//write evaluated Data
void writeData(string year, string evalpath, vector<string> unique_ID)
{
    time_t t = time(0);
    string c_time = ctime(&t);
    string curr_time = c_time.substr(8,2) + c_time.substr(4,3) + c_time.substr(20,4) + "_" + c_time.substr(11,2) + c_time.substr(14,2) + c_time.substr(17,2);
    cout << "Current time: " << curr_time << endl;
    string OutputFileName;
    OutputFileName = "Standards_eval_names_" + year + "_" + curr_time + ".txt";
    cout << "Output File: " << evalpath + "/" + OutputFileName << endl;
    ofstream outFile (evalpath + OutputFileName);
    for (int i = 0; i < unique_ID.size(); i++)
    {
        outFile << unique_ID[i] << endl;
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
    year = files_date[0].substr(0,4);
    //cout << year << endl; //debug

    //////////////////////////////////////////////////
    // Read csv file(s) and store values in vectors //
    //////////////////////////////////////////////////
    vector<string> unique_ID;
    for (int i = 0; i < files_name.size(); i++)
    {
        cout << files_name[i] << " at " << files_adress[i] << endl;
        getNames(files_name[i], files_adress[i], unique_ID);
    };

    //write to file
    writeData(year, evalpath, unique_ID);



    return 0;
}