/* ----------------------------------------------------------------------------------------------------*/
//      FindKmer: Find all kmers from *.fasta files in a directory.
//      Created by Yangbo , Key Lab. for Network-based Intelligent Computing. Univ. of Jinan
//              Sunny, Stanford Medical School, Palo alto, CA, USA
//      Ver.1.00 Mar.31,2019.
//      Ver.1.10 April.3,2019.  -s -o command line arguments added.
//      Ver.2.00 April.4,2019, -s option deleted，output files added position information
//
//      Usage:  FindKmers  -d Directory(default to curent) -k Max_k(default to 4,max 8)
//      eg.:    FindKmers
//              FindKmers  -d ./InputData  -k 3 -o dog.txt
//              FindKmers  -k 6 -s
//
//  Copyright © 2019 Yangbo,Sunny. All rights reserved.
/*------------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <dirent.h>

using namespace std;

#define MAX_K   8

string PathStr;
int TotalFiles=0;
int Max_k=4;        //Default to 4
map <string,vector<int>> mapKmer;
bool bSingleFileOutput=true;            //Default to output each single file
string  OutFileName="Kmers.txt";        //Default to Kmers.txt

// Usage:  FindKmers  -d Directory -k Max_k
void DispHelp(bool bHelp)
{
    string Str1="\n"
    "FindKmer: Find all kmers from *.fasta files in a directory. (ver.2.00)\n"
    "   -Yangbo , Key Lab. for Network-based Intelligent Computing, Univ. of Jinan,China.\n"
    "   -Sunny,Stanford Medical School,Palo Alto,CA,USA.\n"
    "\n";
    string Str2=
    "Usage:  FindKmers  [options]"
    "\n"
    "Options:\n"
    "   -k INT      Max. length of kmers, consider all kmers of length 1,2,...,INT. default to 4, max to 8. \n"
    "   -d STR      Directory path where the FASTA files located, if not given, default to current directory.\n"
    "   -o STR      Output file name, the output file contains kmers for all input files. default to ""Kmers.txt""\n"
    "   if -o is not given, output individal file for each input file using name ended with kmers.txt.\n"
    "\n";
    cout<<Str1;
    if (bHelp) cout<<Str2;
};

int InsertKmer(map<string,vector<int>> &mKmer,string SeqLine)
{
    int K,Count=0;
    vector<int> vPosCount(1,0);          //Create a vector only contain 1 elements for total counter
    string Str;
    
    for (K=1;K<=Max_k;K++)
    {
        if (bSingleFileOutput)
            vPosCount.resize(SeqLine.length()-K+1,0);   //Resize the vector with position cunter elements
        for (int N=0; N<SeqLine.length()-K+1;N++)
        {
            Str=SeqLine.substr(N,K);
            pair<map<string,vector<int>>::iterator, bool> pair1 =mKmer.insert(pair<string,vector<int>>(Str,vPosCount));  //Try to insert the kmer, Number of occurrences =0
            if (bSingleFileOutput)
                pair1.first->second[N]++;           //If inserted and using single file output,increase the counter in this postion
            else
                pair1.first->second[0]++;           //If already exist and -o option, only increase vector[0] as total counter
            if (pair1.second) Count++;              //If inserted new kmer,Increase Total unmber of Kmers
        };
    }
    return Count;           //Return the kmers added in this line
};

int AddKmer(string FName,map<string,vector<int>> &mKmer)
{
    int iCount=0;
    ifstream fFile;                 //File name used for the sequence file
    fFile.open(FName.c_str());
    string Line;
    while(fFile)
    {
        getline(fFile,Line);
        if ((Line.length() == 0) || (Line[0]=='>')) continue;
        iCount+=InsertKmer(mKmer,Line);     //Process this line, find and add kmers,
    }
    fFile.close();
    return iCount;          //Total kmers added in this file
}

int CheckDirFiles(string PathStr)       //Check dir and get number of *.fasta files in the dir
{
    DIR *dir = opendir(PathStr.c_str());
    if (dir==NULL)
    {
       cout<<"Directory "+PathStr+" Not found!"<<endl;
       exit(1);
    }
    dirent* p = NULL;
    int FCount=0;
    while((p = readdir(dir)) != NULL)
    {
        if(p->d_name[0] == '.') continue;  //Skip the "." and ".." in directory
        if (strstr(p->d_name,".fasta"))
            FCount++;       //Get number of filres
    };
    closedir(dir);
    return FCount;
};

void OutputMapFile(map<string,vector<int>> &Map,string FileName)
{
    ofstream Fout;          //File used for output
    Fout.open(FileName);
    for( map<string,vector<int>>::iterator it = mapKmer.begin();it!=mapKmer.end();it++)
    {                   //Write N colums, first is kmer, followed with number of occurrences in each position
        Fout << it->first ;             //The kmer string
        for (int I=0;I<it->second.size();I++)       //The count, if -o option, only vector[0] used for counter
            Fout<<"\t"<< it->second[I];
        Fout<<endl;
    };
    cout<<"Kmers written to file: "+FileName<<endl;
    Fout.close();
    Map.clear();        //Clear the map, ready for next use
}

bool ProcessCMD_Line(int argc, const char *argv[])
{
    if (argc < 2)       // total number of arguments
    {
        DispHelp(true);     // if no argument in command  line, print help and exit
        return false;
    } else
        DispHelp(false);    //Only display program message
    
    string argStr=argv[1];
    PathStr = argv[0];      //Point to the executable FindKmers itself,
    PathStr.erase((int)PathStr.rfind('/'));     //Get the current dir name
    
    for (int i = 1; i < argc; i++)              //Process the  arguments
    {
        if (i != argc) {
            argStr=argv[i];
            if (argStr == "-k") {  // Max_k
                if (argv[i+1])
                    if ((Max_k = atoi(argv[i + 1]))>MAX_K) Max_k=MAX_K;  //Get Max_k
                i=i+1;
            } else if (argStr == "-d") {
                if (argv[i+1]) PathStr = argv[i + 1];          //Get input file path
                i=i+1;          //Skip next argument
            } else if (argStr=="-o") {
                if (argv[i+1])  OutFileName=argv[i+1];          //Get the total output file name
                bSingleFileOutput=false;        //Do not output single file
                i=i+1;      //Skip next argument
            }
        };
    };
    return true;
};
                        //The main procedure starts here
int main(int argc, const char * argv[])
{
    if (!ProcessCMD_Line(argc, argv))  exit(1);  //Exit if no enough command line arguments
    
    int FCount=CheckDirFiles(PathStr);        //Check the dir and get the total number of files
    if (FCount==0)
    {
        cout<<"No .fasta file found. Check the arguments entered in command line."<<endl;
        DispHelp(true);
        exit(1);            //Program ended
    }
    
    int KmersAdded=0;
    dirent* p = NULL;
    DIR *dir = opendir(PathStr.c_str());         //Re-open the dir
    while((p = readdir(dir)) != NULL)
    {
                    //only process *.fasta files, Skip the "." and ".." hiden file name and other files
        if (strstr(p->d_name,".fasta"))
        {
            string FName = PathStr +"/"+string(p->d_name);
            int Count=AddKmer(FName,mapKmer);
            TotalFiles++;           //Increase the counter
            cout<<"File "<<TotalFiles<< "/"<<FCount<<" k="<<Max_k<<": "<<FName<<" processed. "<<Count<<" kmers added."<< endl;
            KmersAdded+=Count;          //Total Count
            if (bSingleFileOutput)
                OutputMapFile(mapKmer, FName.replace(FName.find("fasta"),OutFileName.size(),OutFileName));     //Output the kmers for each input file
        }
    }
    cout<<endl<<TotalFiles<<" files processed. Totally "<<KmersAdded<<" kmers found."<< endl;
    closedir(dir);          //close the directory
   
    if (!bSingleFileOutput)
        OutputMapFile(mapKmer, PathStr+"/"+OutFileName);       //Output all kmers found in every input file into one output file.
}
