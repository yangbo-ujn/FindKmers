/* ----------------------------------------------------------------------------------------------------*/
//      FindKmer: Find all kmers from *.fasta files in a directory. (ver.1.0)
//      Created by Yangbo , Key Lab. for Network-based Intelligent Computing. Univ. of Jinan
//              Sunny, Stanford Medical School, Palo alto, CA, USA
//      Ver.1.00 Mar.31,2019.
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

// Usage:  FindKmers  -d Directory -k Max_k
void DispHelp(void)
{
    string Str="\n"
    "FindKmer: Find all kmers from *.fasta files in a directory. (ver.1.0)\n"
    "   -by Yangbo , Key Lab. for Network-based Intelligent Computing. Univ. of Jinan\n"
    "Usage:  FindKmers  -d Directory(default to curent) -k Max_k(default to 4,max 8)\n"
    "eg.:    FindKmers\n"
    "        FindKmers  -d ./InputData  -k 3\n"
    "        FindKmers  -k 6\n\n";
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

bool ProcessCMD_Line(int argc, const char *argv[])
{
    if (argc < 2)       // total number of arguments
    {
        DispHelp();     // if no argument in command  line, print help and exit
        return false;
    };
    string argStr=argv[1];
    for (int i = 1; i < argc; i++)              //Process the  arguments
    {
        if (i != argc) {
            argStr=argv[i];
            if (argStr == "-k") {  // Max_k
                if ((Max_k = atoi(argv[i + 1]))>MAX_K) Max_k=MAX_K;  //Get Max_k
                i=i+1;
            } else if (argStr == "-d") {
                PathStr = argv[i + 1];          //Get input file path
                i=i+1;
            } else
            {
                PathStr = argv[0];      //Point to the executable FindKmers itself,
                PathStr.erase((int)PathStr.rfind('/'));     //Get the current dir name
            };
        };
    };
    return true;
};
                        //The main procedure starts here
int main(int argc, const char * argv[])
{
    if (!ProcessCMD_Line(argc, argv))   exit(1);  //Exit if no enough command line arguments
    
    int FCount=CheckDirFiles(PathStr);        //Check the dir and get the total number of files
    
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
        }
    }
    cout<<TotalFiles<<" files processed. Totally "<<mapKmer.size()<<" kmers found."<< endl;
    closedir(dir);          //close the directory
    
    ofstream Fout;          //File used for output
    Fout.open(PathStr+"/Kmers.txt");
    for( map<string,int>::iterator it = mapKmer.begin();it!=mapKmer.end();it++)
        Fout << it->first << "\t" << it->second<< endl;         //Write two colums, first is kmer, second is the number of occurrences
    cout<<"Kmers written to file: "+PathStr+"/Kmers.txt"<<endl;
    Fout.close();
}
