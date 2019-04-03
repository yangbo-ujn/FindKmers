/* ----------------------------------------------------------------------------------------------------*/
//      FindKmer: Find all kmers from *.fasta files in a directory.
//      Created by Yangbo , Key Lab. for Network-based Intelligent Computing. Univ. of Jinan
//              Sunny, Stanford Medical School, Palo alto, CA, USA
//      Ver.1.00 Mar.31,2019.
//      Ver.1.10 April.3,2019.  -s -o command line arguments added.
//
//      Usage:  FindKmers  -d Directory(default to curent) -k Max_k(default to 4,max 8)
//      eg.:    FindKmers
//              FindKmers  -d ./InputData  -k 3
//              FindKmers  -k 6
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
map <string,int> mapKmer;
bool bSingleFileOutput=false;           //Default to output a whole file
string  OutFileName="Kmers.txt";

// Usage:  FindKmers  -d Directory -k Max_k
void DispHelp(void)
{
    string Str="\n"
    "FindKmer: Find all kmers from *.fasta files in a directory. (ver.1.10)\n"
    "   -Yangbo , Key Lab. for Network-based Intelligent Computing, Univ. of Jinan,China.\n"
    "   -Sunny,Stanford Medical School,Palo Alto,CA,USA."
    "\n"
    "Usage:  FindKmers  [options]"
    "\n"
    "Options:\n"
    "   -k INT      Max. length of kmers, consider all kmers of length 1,2,...,INT. default to 4, max to 8. \n"
    "   -d STR      Directory path where the FASTA files located, if not given, default to current directory.\n"
    "   -o STR      Output file name, the output file contains kmers for all input files. default to ""Kmers.txt""\n"
    "   -s          -o is ignored, output individal file for each input file using name ended with kmers.txt.\n"
    "\n";
    cout<<Str;
};

int InsertKmer(map<string,int> &mKmer,string Line)
{
    int k,Count=0;
    string Str;
    for (k=1;k<=Max_k;k++)
    {
        for (int N=0; N<Line.length()-k+1;N++)
        {
            Str=Line.substr(N,k);
            pair<map<string,int>::iterator, bool> pair1 =mKmer.insert(pair<string,int>(Str,1));     //Try to insert the kmer, Number of occurrences =1
            if (pair1.second)
                Count++;            //Increase counter if inserted
            else
                pair1.first->second++;      //The same kmer found, only increase the counter of this kmer
        }
    }
    return Count;           //Return the kmers added in this line
};

int AddKmer(string FName,map<string,int> &mKmer)
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

void OutputMapFile(map<string,int> &Map,string FileName)
{
    ofstream Fout;          //File used for output
    Fout.open(FileName);
    for( map<string,int>::iterator it = mapKmer.begin();it!=mapKmer.end();it++)
        Fout << it->first << "\t" << it->second<< endl;         //Write two colums, first is kmer, second is the number of occurrences
    cout<<"Kmers written to file: "+FileName<<endl;
    Fout.close();
    Map.clear();        //Clear the map, ready for next use
}

bool ProcessCMD_Line(int argc, const char *argv[])
{
    if (argc < 2)       // total number of arguments
    {
        DispHelp();     // if no argument in command  line, print help and exit
        return false;
    };
    string argStr=argv[1];
    PathStr = argv[0];      //Point to the executable FindKmers itself,
    PathStr.erase((int)PathStr.rfind('/'));     //Get the current dir name
    
    for (int i = 1; i < argc; i++)              //Process the  arguments
    {
        if (i != argc) {
            argStr=argv[i];
            if (argStr == "-k") {  // Max_k
                if ((Max_k = atoi(argv[i + 1]))>MAX_K) Max_k=MAX_K;  //Get Max_k
                i=i+1;
            } else if (argStr == "-d") {
                PathStr = argv[i + 1];          //Get input file path
                i=i+1;          //Skip next argument
            } else if (argStr=="-s"){
                bSingleFileOutput=true;
            } else if (argStr=="-o") {
                OutFileName=argv[i+1];          //Get the total output file name
                i=i+1;      //Skip next argument
            }
        };
    };
    return true;
};
                        //The main procedure starts here
int main(int argc, const char * argv[])
{
    if (!ProcessCMD_Line(argc, argv))   exit(1);  //Exit if no enough command line arguments
    
    int FCount=CheckDirFiles(PathStr);        //Check the dir and get the total number of files
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
