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
}
int InsertKmer(map<string,int> &mKmer,string Line)
{
    int k,Count=0;
    string Str;
    for (k=1;k<=Max_k;k++)
    {
        for (int N=0; N<Line.length()-k+1;N++)
        {
            Str=Line.substr(N,k);
            pair<map<string,int>::iterator, bool> pair1 =mKmer.insert(pair<string,int>(Str,k));
            if (pair1.second)
                Count++;            //Increase counter if inserted
        }
    }
    return Count;           //Return the kmers added in this line
}

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

int main(int argc, const char * argv[]) {
    if (argc < 2)   // total number of arguments
    {
        DispHelp(); // if only type FindKmer, print help and exit
        exit(1);
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
                PathStr = argv[0];      //The executable FindKmers itself
                PathStr.erase((int)PathStr.rfind('/'));     //Get the current dir
            };
        };
    };

    DIR* dir = opendir(PathStr.c_str());
    if (dir==NULL)
    {
        cout<<"Directory "+PathStr+" Not found!"<<endl;
        return 0;
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
    dir = opendir(PathStr.c_str());         //Re-open the dir
    while((p = readdir(dir)) != NULL)
    {
                    //Skip the "." and ".." hiden file name in UNIX
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
    
    ofstream Fout;
    Fout.open(PathStr+"/Kmers.txt");
    
    for( map<string,int>::iterator it = mapKmer.begin();it!=mapKmer.end();it++)
        Fout << it->first<< endl;
    cout<<"Kmers written to file: "+PathStr+"/Kmers.txt"<<endl;
    Fout.close();
}
