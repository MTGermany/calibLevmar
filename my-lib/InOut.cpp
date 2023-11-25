

//(jun 04) ACHTUNG: InOut gibt nun n_array als Zahl der Eintraege raus!
// max. Index also n_array-1  !!

//arne (may 05): habe <strstream.h> (VERALTET) ausgestauscht gegen sstream
//die klasse heisst nun istringstream anstatt istrstream

// c++ 
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <string>
#include <iomanip>
// #include <vector> in InOut.h because needed there as vector arg needed

 
using namespace std;

#include "InOut.h"


InOut::InOut(){;}

//######################################################
// determine column by the title line, e.g., title line=
// #t[s]   pos[m]  v[m/s]  lane    gap[m]  ID leadID  leadCl
// get_colnumber(filename, "gap[m]") results in 5
// returns 0 if no header found
// returns -1 if file is not accessible/does not exist

int InOut::get_colnumber(const char* fname, string header){ 
  ifstream  infile (fname, ios::in);
  if(!infile){
    cerr << "InOut.get_colnumber: Error opening file " 
	 << fname << " for reading" << endl;
    return -1;
  }

  if(false){
    cout<<"InOut.get_colnumber: fname="<<fname
	<<" header="<<header<<endl;
  }
  bool success=false;
  string line;

  do{
    getline(infile,line);
    prepare_line(line); // transforms ';' or ',' into spaces
    istringstream linestream (line);
    int icol=0;
    do{
      string word;
      linestream>>word;
      success=(word.compare(header)==0);
      if(false){
	cout<<"InOut.get_colnumber: line="<<line<<" word="<<word
	    <<" header="<<header<<" icol="<<icol
	    <<" success="<<success<<endl;
      }
      icol++;
    }
    while (linestream && (!success));

    if(success){return icol;}
  }
  while( (!infile.eof()) && (!(is_data_line(line))) && (!success));

  return 0;
  }





// string as arg does not work as filename. 
// Use string.c_str() if calling with string

bool InOut::fileExists(const char* fname){ 
  ifstream  infile (fname, ios::in);

  // formulate so cumbersom since nature of infile
  // (bool? file descriptor?) obscure

  if(!infile){return false;}
  else return true;
}


// get number of data lines (excluding empty line and comments)

int InOut::getNumberOfDataLines(const char* fname){
  ifstream  infile (fname, ios::in);
  if(!infile){
    cerr << "InOut.getNumberOfDataLines: Error opening file " 
	 << fname << " for reading" << endl;
    cout<<" setting number of data lines to -1"<<endl; // exit(-1);
    return -1;
  }
  
  int i=0;
  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  ){
    if( is_not_white(line,strlen(line)) && is_data_line(line) ){
      i++;
    }
  }
  infile.close();
  if(false){
    cout << "getNumberOfDataLines: file="<<fname
         <<" number of lines = " <<i<<endl;
  }
  return(i);
}

// get overall number of lines (including empty line and comments)

int InOut::getNumberOfLines(const char* fname){
  ifstream  infile (fname, ios::in);
  if(!infile) {
    cerr << "InOut.getNumberOfLines:Error opening file " 
	 << fname << " for reading" << endl;
    cerr<<" setting number of data lines to -1"<<endl; // exit(-1);
    return -1;
  }
  
  int i=0;
  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  ){
    i++;
  }
  infile.close();
  if(false){
    cout << "getNumberOfLines: file="<<fname
         <<" number of lines = " <<i<<endl;
  }
  return(i);
}


int InOut::getNumberOfCols(const char* fname){
  cout<<"in InOut.getNumberOfCols:: " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if(!infile) {
      cerr << "Error opening file " << fname << " for reading" << endl;
      exit(-1);
  }
  
  char line[LINEMAX];
  bool dataLineReached=false;

  //while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  //achtung darf nicht aufhoeren bei leerzeile...
  while ( (infile.getline(line,LINEMAX))&&(!dataLineReached)  ){
    if( is_not_white(line,strlen(line)) && is_data_line(line) ){
      dataLineReached=true;
    }
  }
  infile.close();
  int ncol=0;
  if(!dataLineReached){
    cerr<<"no data line found in "<<fname<<" return zero columns"<<endl;
    return(0);
  }
  else{
    istringstream linestream (line);
    double proforma;
    while(linestream){
      linestream >> proforma;
      ncol++;
    }
    ncol--;
    cout << "number of cols = " <<ncol<<endl;
  }
  return(ncol);
}

// get first line of file (e.g., title)

string InOut::getFirstLine(const char* fname){
  ifstream  infile (fname, ios::in);
  if( !infile){
    cerr << "Critical: InOut.getFirstDataLine: Error opening file " << fname 
	 << " for reading,"<<endl<<"          returning zero string" << endl;
    return "";
  }
  string line;
  getline(infile,line);
  return line;
}


// get first data line of file (e.g., to determine first date)

string InOut::getFirstDataLine(const char* fname){
  ifstream  infile (fname, ios::in);
  if( !infile){
    cerr << "Critical: InOut.getFirstDataLine: Error opening file " << fname 
	 << " for reading,"<<endl<<"          returning zero string" << endl;
    return "";
  }


  string line;
  string output="";
  bool success=false;
  while( (!success)&& getline(infile,line)){
    if(is_data_line(line)){
      output=line;
      success=true;
    }
  }
  if(!success){
    cerr<<"Critical: file " << fname << "has not a single data line"<<endl;
    return "";
  }
  else{
    return output;
  }
}

// get last data line of file (e.g., to determine last date)

string InOut::getLastDataLine(const char* fname){
  ifstream  infile (fname, ios::in);
  if( !infile){
    cerr << "Critical: InOut.getLastDataLine: Error opening file " << fname 
	 << " for reading,"<<endl<<"          returning zero string" << endl;
    return "";
  }


  string line;
  string output="";
  int nlines=0;
  while(getline(infile,line)){
    if(is_data_line(line)){
      output=line;
      nlines++;
    }
  }
  if(nlines==0){
    cerr<<"Critical: file " << fname << "has not a single data line"<<endl;
    return "";
  }
  else{
    return output;
  }
}







void  InOut::get_array(char* fname, int & n_array, double array1[])
{
  cout << "in InOut.get_array (double array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  //while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  while ( infile.getline(line,LINEMAX) )
  {
    if( is_data_line(line) && is_not_white(line,strlen(line)) ){
      prepare_line(line);
      istringstream linestream (line);
      linestream >> array1[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}

void  InOut::get_array (char* fname, int & n_array, int array1[]){
  cout << "in InOut.get_array (int array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      linestream >> array1[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}

//#######################################################################

// MT 2019-03 get only the last nlast data elements 
// of column col of file fname

void  InOut::get_col_last (const char* fname, int col, int nlast, 
			   int & n_array, int array1[]){

  ifstream  infile (fname, ios::in);
  if( !infile){
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }

  int ndata=getNumberOfDataLines(fname);

  if(ndata<=nlast){
    cerr<<"InOut.get_col_last: file "<<fname<<" has only "
	<<ndata<<" data lines which is <= nlast="<<nlast<<endl
	<<" providing the whole data column of "<<fname<<endl;
    get_col(fname, col, n_array, array1, false);
  }

  else{
    n_array=nlast;
    int i=0, count=0;
    int errcnt=0;
    char line[LINEMAX];
    char proforma[LINEMAX]; 
    while ( (infile.getline(line,LINEMAX))) {
      if(is_data_line(line)){i++;}
      if(i>ndata-nlast){
        prepare_line(line);
        istringstream linestream (line);
        for (int icol=1; icol<col; icol++){linestream >> proforma;}
        linestream >> array1[count];

	//error treatment
        if (!linestream){
	  errcnt++;
	  //cerr<<"error while parsing for double ... i="<<i<<" --> "
	  //   <<array1[i]<<" set to -999 errcnt="<<errcnt<<endl;
          array1[count]=-999;
	}
	count++;
      }
    }
    infile.close();
    n_array = count;
    if(errcnt>0){cerr<<"InOut.get_col: file="<<fname<<" col="<<col
		     <<": WARNING: not all data were int! Set to -999 in "
		     <<errcnt<<" cases!"<<endl;
    }
  }
}




//#######################################################################


// get column col from file fname 
// if log=true => log output (which file read etc)

void  InOut::get_col(const char* fname, int col,  int &n_array, int array1[]){
  get_col(fname, col,  n_array, array1, false);
}


void  InOut::get_col(const char* fname, int col,  int & n_array, int array1[], bool log){
  ifstream  infile (fname, ios::in);
  if( !infile){
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 

  int errcnt=0;
  while ( (infile.getline(line,LINEMAX))) {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      for (int icol=1; icol<col; icol++){linestream >> proforma;}
      linestream >> array1[i];

	//error treatment: MT dec18

      if (!linestream){
	  errcnt++;
	  //cerr<<"error while parsing for double ... i="<<i<<" --> "
	  //   <<array1[i]<<" set to -999 errcnt="<<errcnt<<endl;
          array1[i]=-999;
      }

      i++;
    }
  }
  infile.close();
  n_array = i;
  if(log){
    cout << "InOut.get_col(int): " << "fname = "<<fname 
         <<" col="<<col<<" n_array = " << n_array << endl;
  }
  if(errcnt>0){cerr<<"InOut.get_col: file="<<fname<<" col="<<col
		   <<": WARNING: not all data were int! Set to -999 in "
		   <<errcnt<<" cases!"<<endl;
  }

}

void  InOut::get_col(const char* fname, int col,  int &n_array, double array1[]){
  get_col(fname, col,  n_array, array1, false);
}


void  InOut::get_col(const char* fname, int col,  int &n_array, double array1[], bool log){
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  
  int errcnt=0;
  while ( (infile.getline(line,LINEMAX))) {
		
      if(is_data_line(line)){
	prepare_line(line);
	istringstream linestream (line);
	for (int icol=1; icol<col; icol++){linestream >> proforma;}
	linestream>>array1[i];

	//error treatment: arne feb 08

	if (!linestream){
	  errcnt++;
          array1[i]=-999;
	  if(false){
	    cerr<<"error while parsing for double ... i="<<i
	        <<" col="<<col<<" datum="<<array1[i]
	        <<" line=" <<line<<endl;
	  }
	} 
	i++;
      }
    }
	
  infile.close();
  n_array = i;
  if(log){
    cout << "get_col(double): " << "fname = "<<fname 
         <<" col="<<col<<" n_array = " << n_array << endl;
  }
  if(errcnt>0){cerr<<"InOut.get_col: file="<<fname<<" col="<<col
		   <<": WARNING: not all data were double! Set to -999 in "
		   <<errcnt<<" cases!"<<endl;
  }
}



void  InOut::getChar_col(const char* fname, int col, int & n_array, char array1[]){
  cout << "in InOut.getChar_col( char array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      for (int icol=1; icol<=col; icol++){linestream >> proforma;}
      array1[i]=proforma[0];
      i++;
    }
  }
  infile.close();
  n_array = i;
  //cout << "n_array = " << n_array << endl;
}



//arne 9-7-07
//MT2018-12: string as input parameter just does not work!! SUCK!!


void InOut::get_col(const char* fname,int col,int & n_array,string array1[]){ 
  get_col(fname, col,  n_array, array1, false);
}

void  InOut::get_col(const char* fname, int col,  int & n_array, 
		      string array1[], bool log){

  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  while ( (infile.getline(line,LINEMAX))) 
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      for (int icol=1; icol<col; icol++){linestream >> proforma;}
      linestream >> array1[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  if(log){
    cout << "in InOut.get_col(string array): " 
	 << "fname = "<<fname<<" col="<<col<<" n_array = "<< n_array<<endl;
  }
}

//###################################################################
//MT2023-02: get_col with vector argument

void InOut::get_col(const char* fname, int col,  int & n_array,
		    vector<double>& data){


  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }

  // start with an empty vector; otherwise push_back history dependent
  data.clear(); 

  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  double oneElement; // used for push_back
  
  int errcnt=0;
  while ( (infile.getline(line,LINEMAX))) {
    
		
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      for (int icol=1; icol<col; icol++){linestream >> proforma;}
      linestream>>oneElement;

	//error treatment for ints and doubles

      if (!linestream){
	  errcnt++;
          oneElement=-999;
	  if(false){
	    cerr<<"error while parsing for double ... i="<<i
	        <<" col="<<col<<" datum="<<oneElement
	        <<" line=" <<line<<endl;
	  }
      }
      
      data.push_back(oneElement);
      
      i++;
    }
  }
	
  infile.close();
  n_array = i;
  if(true){
    cout << "get_col(vector<double>&): " << "fname = "<<fname 
         <<" col="<<col<<" n_array = " << n_array << endl;
  }
  if(errcnt>0){cerr<<"InOut.get_col: file="<<fname<<" col="<<col
		   <<": WARNING: not all data were double! Set to -999 in "
		   <<errcnt<<" cases!"<<endl;
  }
}


void InOut::get_col(const char* fname, int col,  int & n_array,
		    vector<int>& data){


  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }

  // start with an empty vector; otherwise push_back history dependent
  data.clear(); 

  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  double oneElement; // used for push_back
  
  int errcnt=0;
  while ( (infile.getline(line,LINEMAX))) {
    
		
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      for (int icol=1; icol<col; icol++){linestream >> proforma;}
      linestream>>oneElement;

	//error treatment for ints and doubles

      if (!linestream){
	  errcnt++;
          oneElement=-999;
	  if(false){
	    cerr<<"error while parsing for double ... i="<<i
	        <<" col="<<col<<" datum="<<oneElement
	        <<" line=" <<line<<endl;
	  }
      }
      
      data.push_back(oneElement);
      
      i++;
    }
  }
	
  infile.close();
  n_array = i;
  if(true){
    cout << "get_col(vector<double>&): " << "fname = "<<fname 
         <<" col="<<col<<" n_array = " << n_array << endl;
  }
  if(errcnt>0){cerr<<"InOut.get_col: file="<<fname<<" col="<<col
		   <<": WARNING: not all data were double! Set to -999 in "
		   <<errcnt<<" cases!"<<endl;
  }
}








//###################################################################

void InOut::get_array (char* fname, int & n_array, 
      double array1[], double array2[])
{
  cout << "in InOut.get_array (double, double): " 
       << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;

  //  while (infile >> array1[i] >> array2[i]) i++;

  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      linestream >> array1[i] >> array2[i];
      cout << "array1[" << i <<"] = " << array1[i] 
             << " array2[" << i <<"] = " << array2[i] << endl;
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}


void  InOut::get_array (char* fname, int & n_array, 
      double array1[], double array2[], double array3[])
{
  cout << "in InOut.get_array(double,double,double): " 
       << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;

  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
        istringstream linestream (line);
        linestream >> array1[i] >> array2[i] >> array3[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}

void  InOut::get_array (char* fname, int & n_array, 
      int array1[], double array2[], double array3[])
{
  cout << "in InOut.get_array (int, double, double): " 
       << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;

  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
       linestream >> array1[i] >> array2[i] >> array3[i];
      cout <<"InOut::get_array(*,*,int[], ...):"
	   <<" line="<<line<<endl
	   <<" i="<<i<<" array1[i]="<<array1[i]<<endl;
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout <<"n_array = " << n_array << endl;
}






//################################################################
// determine grid numbers n1 and n2 of a data file representing a matrix
// each noncommented data line is organized to represent 
// data[0][0]
// data[0][1]
// ......
// data[0][n2-1]
// <empty line>
// data[1][0]
// .......
// data[n1-1][n2-1]
//################################################################


void InOut::getGridnumbers_n1n2(char* fname, int & n1, int & n2){
  cout << "in InOut.getGridnumber_n1n2: fname="<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }

  const int MAXSTR=500;
  char line[MAXSTR];
  
  int i1=0;
  int i2=0;
  while (!infile.eof()){
    i2=0;
    while ( (infile.getline(line,MAXSTR))
	    && (is_not_white(line,strlen(line)))  ){
      if( (line[0] != COMMENTCHAR) && 
          (line[0] != COMMENTCHAR2) &&
          (strlen(line)>0) ){
        i2++;
      }
    }
    if(i1==0){n2=i2;}
    i1++;
  }
  n1=i1-1;
  cout<<"InOut.getGridnumber_n1n2: n1="<<n1<<" n2="<<n2<<endl;
}



//################################################################
// grid data extraction, 
// if swap=false, then
// matrix[i1][i2] is read in from column column of a data file organized as
// icol=matrix
// matrix[0][0]
// matrix[0][1]
// ......
// matrix[0][n2-1]
// matrix[1][0]
// .......
// matrix[n1-1][n2-1]
//
// if swap=true, the transposed matrix matrix[i2][i1] is read
// empty line in input not needed here; 
// was used to get n1,n2 in previous function getGridnumbers_n1n2(...)
// use two function to be able to dimension the matrix before filling it
//################################################################

void InOut::get_matrix (char* fname, 
			int n1, int n2, bool swap,
			int column, double** matrix){

  cout << "in InOut.get_matrix: fname="<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  const int MAXSTR=500;
  char line[MAXSTR];

  int i1=0;
  int i2=0;
  int nrow=(swap) ? n2 : n1;
  int ncol=(swap) ? n1 : n2;
  double dummy;
  
  while ( (infile.getline(line,MAXSTR)) && (i1<n1) ){
    if( (line[0] != COMMENTCHAR) && 
        (line[0] != COMMENTCHAR2) &&
	is_not_white(line,strlen(line)) &&
        (strlen(line)>0) ){

      int irow=(swap) ? i2 : i1;
      int icol=(swap) ? i1 : i2;
      std::istringstream linestream (line);
      for(int ic=1; ic<column; ic++){linestream >> dummy;}
      linestream >>  matrix[irow][icol];

      if(false){
        cout<<" line="<<line<<endl;
        cout<<"i1="<<i1<<" i2="<<i2<<" irow="<<irow<<" icol="<<icol
	    <<" column="<<column<<endl;
        cout <<"  matrix[irow][icol]="<<matrix[irow][icol]<<endl<<endl;
        if(i1==3){exit(0);}
      }
      
      i2++;
      if(i2==n2){
	i1++;
	i2=0;
      }
    }
  }

  if(true){
    cout<<"Finished InOut.get_matrix:"<<endl;
    cout << "nrow="<<nrow<<" ncol="<<ncol<<" column="<<column<< endl;
    cout << "matrix["<<0<<"]["<<0<<"]="<< matrix[0][0]<<endl;
    cout << "matrix["<<0<<"]["<<1<<"]="<< matrix[0][1]<<endl;
    cout << "matrix["<<0<<"]["<<2<<"]="<< matrix[0][2]<<endl;
    cout << "matrix["<<0<<"]["<<3<<"]="<< matrix[0][3]<<endl;
    cout << "matrix["<<1<<"]["<<1<<"]="<< matrix[1][1]<<endl;
    cout << "matrix["<<2<<"]["<<1<<"]="<< matrix[2][1]<<endl;
     
    cout << "matrix["<<0<<"]["<<ncol-1<<"]="<< matrix[0][ncol-1]<<endl;
    cout << "matrix["<<nrow-1<<"]["<<0<<"]="<< matrix[nrow-1][0]<<endl;
    cout << "matrix["<<nrow-1<<"]["<<ncol-1<<"]="<< matrix[nrow-1][ncol-1]
	 <<endl;
  }
}







//################################################################
//  even more specialized function function  get_array2d 
//################################################################

void  InOut::get_array2d (char* fname, 
                  double unitTimeInDat_s, double unitSpaceInDat_m,
                  int & nt, int & nx, double & dt, double & dx, 
                  double array1[NXMAX][NTMAX+1],
                  double array2[NXMAX][NTMAX+1],
                  double array3[NXMAX][NTMAX+1])

{
  cout << "in get_array2d: " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  const int MAXSTR=500;
  char line[MAXSTR];
  double xvals[NXMAX];
  double tvals[NTMAX];

  int i=0;
  while (!(infile.eof())){

    int j=0;
    if (i>=NTMAX){
      cerr << "InOut.get_array2d: error: time index (second coordinate) i="<<i
	   <<" >=  NTMAX="<<NTMAX<<endl;
      exit(-1);
    }
    while (   (infile.getline(line,MAXSTR))  
           && (is_not_white(line,strlen(line))) 
          )

    {
      if( (line[0] != COMMENTCHAR) && 
          (line[0] != COMMENTCHAR2) &&
          (strlen(line)>0) )
      {
	//	cout <<"i="<<i<<" j="<<j<<endl;

        std::istringstream linestream (line);
        linestream >>  xvals[j] >>  tvals[i] 
                   >> array1[j][i] >> array2[j][i] >> array3[j][i];
        if(false)cout << "array1[" << i <<"][" << j <<"] = " << array1[j][i] << endl;
        j++;
      }
      nx = j-1; //!!
      nx = j;
    }
    i++;
  }

  infile.close();
  nt = i-2; //!!
  nt = i-1;
  dx = (nx>=2) ? unitSpaceInDat_m*(xvals[nx-1] - xvals[nx-2]) : xvals[nx-1];
  dt = (nt>=2) ? unitTimeInDat_s*(tvals[nt-1] - tvals[nt-2]) : xvals[nt-1];

  //  for (i=0; i<=nt; i++) cout << "tvals[i]="<<tvals[i]<<endl;
  cout << "nt = " << nt << "; nx = " << nx << endl;
  cout << "dt = " << dt << "; dx = " << dx << endl;
  cout << "array1[" << nx-1 <<"][" << nt-1 <<"] = " << array1[nx-1][nt-1] << endl
       << "array1[" << 0 <<"][" << 0 <<"] = " << array1[0][0] << endl
       << "array2[" << nx-1 <<"][" << nt-1 <<"] = " << array2[nx-1][nt-1] << endl
       << "array3[" << nx-1 <<"][" << nt-1 <<"] = " << array3[nx-1][nt-1] << endl;

  //exit(0);

}




void  InOut::write_array (char* fname, int nData, const double data_col1[],
			  char* name_col1){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", name_col1);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\n", data_col1[i]);
  }
  fclose(outfile);
} // InOut::write_array


void  InOut::write_array (char* fname, int nData, const double data_col1[],
const double data_col2[], char* name_col1, char* name_col2)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\n", name_col1, name_col2);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%f\n", data_col1[i], data_col2[i]);
  }
  fclose(outfile);
} // InOut::write_array



//<new apr19>

void InOut::write_array (char* fname, int nData, const double data_col1[],
			 const double data_col2[], 
			 char* titleString){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%1.4f\n", 
	    data_col1[i], data_col2[i]);
  }
  fclose(outfile);
} 

void InOut::write_array (char* fname, int nData, const int data_col1[],
			 const double data_col2[], 
			 char* titleString){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%i\t%1.3f\n", 
	    data_col1[i], data_col2[i]);
  }
  fclose(outfile);
} 

//</new>

void  InOut::write_array (char* fname, int nData, const int data_col1[],
const int data_col2[], char* name_col1, char* name_col2)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\n", name_col1, name_col2);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%i\t%i\n", data_col1[i], data_col2[i]);
  }
  fclose(outfile);
} // InOut::write_array

void  InOut::write_array (char* fname, int nData, const double data_col1[],
const double data_col2[], const double data_col3[], 
char* name_col1, char* name_col2, char* name_col3)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\n", name_col1, name_col2, name_col3);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%f\t%f\n", data_col1[i], data_col2[i], data_col3[i]);
  }
  fclose(outfile);

} // InOut::write_array


//<new aug16>

void InOut::write_array (char* fname, int nData, const double data_col1[],
			 const double data_col2[], const double data_col3[],
			 char* titleString){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%1.7f\t%1.7f\n", 
	    data_col1[i], data_col2[i], data_col3[i]);
  }
  fclose(outfile);
} 

//</new>



void  InOut::write_array (char* fname, int nData, const double data_col1[],
const double data_col2[], const double data_col3[], const double data_col4[], 
char* name_col1, char* name_col2, char* name_col3, char* name_col4)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\n", name_col1, name_col2, 
        name_col3, name_col4);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%6.8f\t%f\t%6.8f\n", data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i]);
  }
  fclose(outfile);

} // InOut::write_array

void  InOut::write_array (char* fname, int nData, const double data_col1[],
const int data_col2[], const double data_col3[], const double data_col4[], 
char* name_col1, char* name_col2, char* name_col3, char* name_col4)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\n", name_col1, name_col2, 
        name_col3, name_col4);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%i\t%f\t%6.8f\n", data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i]);
  }
  fclose(outfile);
} // InOut::write_array


//<new nov17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.4f\t%5.5f\t%5.5f\t%5.5f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i]);
  }
  fclose(outfile);
} 
//</new>


void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
		     const double data_col5[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%f\t%f\t%f\t%f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i]);
  }
  fclose(outfile);

} // InOut::write_array

void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const int data_col2[], 
                     const double data_col3[],const double data_col4[],
		     const double data_col5[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%i\t%f\t%f\t%f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i]);
  }
  fclose(outfile);

} // InOut::write_array




//<new nov17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
		     const double data_col5[],
		     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
      //fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
    fprintf(outfile, "%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i]);
  }
  fclose(outfile);
} 
//</new>



void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%f\t%f\t%f\t%f\t%f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i]);
  }
  fclose(outfile);
} 


//<new nov17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
		     const double data_col5[],const double data_col6[],
		     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i]);
  }
  fclose(outfile);
} 
//</new>




void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6, 
        name_col7);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i]);
  }
  fclose(outfile);

} // InOut::write_array


//<new mar17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],
                     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i]);
  }
  fclose(outfile);

} // InOut::write_array





void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6, 
        name_col7, name_col8);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i]);
  }
  fclose(outfile);

} // InOut::write_array

//<new mar17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],const double data_col8[],
                     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i]);
  }
  fclose(outfile);

} // InOut::write_array







void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8, 
                     char* name_col9)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6, 
        name_col7, name_col8, 
        name_col9);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i]);
  }
  fclose(outfile);

} // InOut::write_array

//<new mar17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],const double data_col8[],
                     const double data_col9[],
			  char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i]);
  }
  fclose(outfile);

} // InOut::write_array









void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],const double data_col10[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8, 
                     char* name_col9, char* name_col10)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6, 
        name_col7, name_col8, 
        name_col9, name_col10);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i], data_col10[i]);
  }
  fclose(outfile);

} // InOut::write_array

//<new mar17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],const double data_col8[],
                     const double data_col9[],const double data_col10[],
		     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i], data_col10[i]);
  }
  fclose(outfile);

} // InOut::write_array

void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],const double data_col8[],
                     const double data_col9[],const double data_col10[],
                     const double data_col11[],const double data_col12[],
		     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i], data_col10[i],
	    data_col11[i], data_col12[i]);

  }
  fclose(outfile);

} // InOut::write_array


//##########################################################
// general-purpose function  write_array2d 
// nata  => nLines
// nCols => generalizes above versions of write_array
//##########################################################

void  InOut::write_array_gen (char* fname, int nLines, int nCols,  
                     const double data[][LINEMAX], const char* headers){
  FILE  *outfile;                 
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }
  
  fprintf(outfile, "%s\n", headers);
  for(int i=0; i<nLines; i++){
    stringstream ss;
    ss<< setprecision(5); // with <<fixed<<setprecision(2) 2 decimal digits

    for(int icol=0; icol<nCols; icol++){
      ss<<data[icol][i]<<'\t';
    }
    string line=ss.str();
    fprintf(outfile, "%s\n", line.c_str());
    
    /*
    char line[MAXSTR];
    sprintf(line,""); // even if compiler complains: This line *is* necessary; otherwise overflow!
    for(int icol=0; icol<nCols; icol++){
      //cout <<"test: i="<<i<<" icol="<<icol<<" data[icol][i]="<<data[icol][i]<<endl;
      sprintf(line,"%s%5.4f\t", line, data[icol][i]);
    }
    //cout<<"test: line="<<line<<endl;
    fprintf(outfile, "%s\n", line);
    */
    
  }
}

void  InOut::write_array_gen (char* fname, int nLines, int nCols,  
                     const int data[][LINEMAX], const char* headers){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "%s\n", headers);
  for(int i=0; i<nLines; i++){

    stringstream ss;
    for(int icol=0; icol<nCols; icol++){
      ss<<data[icol][i]<<'\t';
    }
    string line=ss.str();
    fprintf(outfile, "%s\n", line.c_str());

    
    /*
    char line[MAXSTR];
    sprintf(line,""); // even if compiler complains: This line *is* necessary; otherwise overflow!
    for(int icol=0; icol<nCols; icol++){
      sprintf(line,"%s%i\t", line, data[icol][i]);
    }
    fprintf(outfile, "%s\n", line);
    */

    
  }
}



//########################################################################
//  function  get_array2d 
//########################################################################

/// reads in array[i][j] representing z(x_i,y_j) with
/// regular spacings 
/// x_i=xmin+i*dx, i=0 .. nx-1, 
/// y_j=ymin+j*dy, y=0 .. ny-1
/// from column col (>=3) of file fname (column count starts at 1); 
/// index i (first col) is the same in each block, 
/// index j (second col) varies every line; empty line between 
/// two blocks 
/// Output: 
/// xmin,dx,nx (calculated from first col)
/// ymin,dy,ny (calculated from second col)
/// array[i][j]

void  InOut::get_array2d_col(const char* fname, int col, 
		       double& xmin, double& dx, int& nx, 
		       double& ymin, double& dy, int& ny, 
		       double array[NXMAX][NYMAX]){
  cout << "in InOut.get_array2d: " 
       << "fname = "<<fname<<" col="<<col <<endl;
  if(col<3){
    cerr<<" InOut.get_array2d: Error: column index must be >=3"
	<<endl;
    exit(-1);
  }
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "InOut.get_array2d: Error opening file " << fname 
         << " for reading" << endl;
    exit(-1);
  }
  
  char line[LINEMAX];
  char proforma[LINEMAX]; 

  double xvals[NXMAX];
  double yvals[NYMAX];

  int i=0;
  int i_in_loop=0;
  while (!(infile.eof())){

    int j=0;
    while (   (infile.getline(line,LINEMAX))  
           && (is_not_white(line,strlen(line))) 
          ){
    if(is_data_line(line)){
      prepare_line(line);
        istringstream linestream (line);
        linestream >>  xvals[i] >>  yvals[j];
        for (int icol=3; icol<col; icol++){linestream >> proforma;}
        linestream >> array[i][j];
        cout << "array[" << i <<"][" << j <<"] = " << array[i][j] << endl;
        j++;
	i_in_loop=i;
      }
      ny = j;
    }
    i++;
    cout <<"after i++: i="<<i<<endl;
  }
  //nx = i-2; ???
  nx = i_in_loop+1;

  infile.close();

  xmin=xvals[0];
  ymin=yvals[0];
  dx = (nx>1) ? xvals[nx-1] - xvals[nx-2] : xvals[nx-1];
  dy = (ny>1) ? yvals[ny-1] - yvals[ny-2] : yvals[ny-1];

  if(true){
    cout <<"InOut.get_array2d_col successfully:\n";
    cout << "xmin = " << xmin<< "\tdx = " << dx << "\tnx = " << nx << endl;
    cout << "ymin = " << ymin<< "\tdy = " << dy << "\tny = " << ny << endl;
    cout << "array[" << 0  <<"][" << 0  <<"] = " << array[0][0]    << endl
	 << "array[" << 0  <<"][" << ny-1 <<"] = " << array[0][ny-1]  << endl
	 << "array[" << nx-1 <<"][" << 0  <<"] = " << array[nx-1][0]  << endl
	 << "array[" <<nx-1<<"]["<< ny-1 <<"] = " << array[nx-1][ny-1] << endl;
    }

}

//##########################################################
// general-purpose function  write_array2d 
// writes gnu-plottable file with col1=x value, col2=y value, col3=array=z(z,y) 
//##########################################################

// definition of arrays for calling write_array2d: see HINT in .h file
// see also write_array2d_fromarray1d

void  InOut::write_array2d (char* fname,
		           double xmin, double xmax, int nx,
		           double ymin, double ymax, int ny,
			    double array[][NYMAX],
			    char* titleString){

  if((nx>NXMAX)||(ny>NYMAX)){
    cerr <<"InOut::write_array2d:" 
                 <<" nx="<<nx<<" >= NXMAX="<<NXMAX<<" and/or "
	 <<" ny="<<ny<<"  >= NYMAX="<<NYMAX
	 <<" => Error!" << endl;
    exit(-1);
  }

  double dx=(xmax-xmin)/((double)(nx-1));
  double dy=(ymax-ymin)/((double)(ny-1));
  cout << "in write_array2d: "
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<" xmax="<<xmax <<" nx="<<nx<<" dx="<<dx
       <<" ymin="<<ymin<<" ymax="<<ymax <<" ny="<<ny<<" dy="<<dy <<endl;

  FILE  *outfile;                  

  outfile = fopen(fname,"w");

  fprintf(outfile, "#%s \n",titleString);

  for (int ix=0; ix<nx; ix++){
    for (int iy=0; iy<ny; iy++){
      double x=xmin+ix*dx;
      double y=ymin+iy*dy;
      fprintf(outfile,  "%.6f\t  %.6f\t %.6f\n", x,y,array[ix][iy]);
    }
    fprintf(outfile, "\n");        // newline mean new "sweep" for gnuplot
  }

  fclose(outfile);

} // end general-purpose write_array2d

//##########################################################
// general-purpose function  write_array2d_fromarray1d 
// writes gnu-plottable file with col1=x value, col2=y value, 
// col3=array=z(z,y) where z[i][j]=array1d[nx*i+j]
//##########################################################

// definition of arrays for calling write_array2d: see HINT in .h file

void  InOut::write_array2d_fromarray1d (char* fname,
		           double xmin, double xmax, int nx,
		           double ymin, double ymax, int ny,
			    double array1d[],
			    char* titleString){


  double dx=(xmax-xmin)/((double)(nx-1));
  double dy=(ymax-ymin)/((double)(ny-1));
  cout << "in write_array2d_fromarray1d: "
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<" xmax="<<xmax <<" nx="<<nx<<" dx="<<dx
       <<" ymin="<<ymin<<" ymax="<<ymax <<" ny="<<ny<<" dy="<<dy <<endl;

  FILE  *outfile;                  

  outfile = fopen(fname,"w");

  fprintf(outfile, "#%s \n",titleString);

  for (int ix=0; ix<nx; ix++){
    for (int iy=0; iy<ny; iy++){
      double x=xmin+ix*dx;
      double y=ymin+iy*dy;
      fprintf(outfile,  "%.6f\t  %.6f\t %.6f\n", x,y,array1d[nx*ix+iy]);
    }
    fprintf(outfile, "\n");        // newline mean new "sweep" for gnuplot
  }

  fclose(outfile);

} // end general-purpose write_array2d_fromarray1d

//#########################################################
//  special-purpose function  write_array2d (general purpose: above)
//#########################################################

  // definition of arrays for calling write_array2d: see HINT in .h file

void  InOut::write_array2d (char* fname, double dxout, double dtout, 
			    double xmin, double xmax, double tmin,
			    double tmax,
			    double array1[][NYMAX])
			     //double** array1,double** array2,double** array3)

{

  int nxout = (int)(0.5+(xmax-xmin)/dxout)  +1;
  int ntout = (int)(0.5+(tmax-tmin)/dtout)  +1;
  if((nxout>NXMAX)||(ntout>NYMAX)){
    cerr <<"InOut::write_array2d:" 
         <<" nxout="<<nxout<<">=NXMAX="<<NXMAX<<" and/or ntout="<<ntout
	 <<" >=NYMAX="<<NYMAX<<" => Error!" << endl;
    exit(-1);
  }

  cout << "in write_array2d: "<<endl
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<", xmax="<<xmax <<endl
       <<" tmin="<<tmin<<", tmax="<<tmax <<endl
       <<" dxout="<<dxout<<", dtout="<<dtout <<endl
       <<" nxout="<<nxout<<", ntout="<<ntout
       <<endl;
  FILE  *outfile;                  
  //  filecheck (outfile,fname);
  outfile = fopen(fname,"w");

  fprintf(outfile, 
    "# x(km) \t t(min) \t rho(1/km) \t v(km/h) \t Q(1/h) \t a(m^2/s^2) \n");

  for (int it=0; it<ntout; it++)
  {
  
    for (int j=0; j<nxout; j++)
    {
      double x_km          = (xmin+(xmax-xmin)*j/nxout)/1000.;
      double t_h           = (tmin+(tmax-tmin)*it/ntout)/3600.;
      double rho_invkm     = array1[j][it];
      //double v_kmh         = array2[j][it];
      //double Q_invh        = array3[j][it];
      double v_kmh = 0;
      double Q_invh = 0;
      fprintf(outfile,  "%.6f\t  %.6f\t %.6f\t  %.6f\t  %f\n",
                      x_km, t_h, rho_invkm, v_kmh, Q_invh);
    }
    fprintf(outfile, "\n");        // newline mean new t step for gnuplot
  }

  fclose(outfile);

}




//#########################################################
//  special-purpose function  write_array2d (general purpose: above)
//#########################################################

  // definition of arrays for calling write_array2d: see HINT in .h file

void  InOut::write_array2d (char* fname, double dxout, double dtout, 
    double xmin, double xmax, double tmin, double tmax,
    double array1[][NYMAX], double array2[][NYMAX],
    double array3[][NYMAX], double array4[][NYMAX])
		      //double** array1,double** array2,double** array3)

{

  int nxout = (int)((xmax-xmin)/dxout);
  int ntout = (int)((tmax-tmin)/dtout);
  if((nxout>NXMAX)||(ntout>NYMAX)){
    cerr <<"InOut::write_array2d:" 
         <<" nxout="<<nxout<<">=NXMAX="<<NXMAX<<" and/or ntout="<<ntout
	 <<" >=NYMAX="<<NYMAX<<" => Error!" << endl;
    exit(-1);
  }

  cout << "in write_array2d: "<<endl
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<", xmax="<<xmax <<endl
       <<" tmin="<<tmin<<", tmax="<<tmax <<endl
       <<" dxout="<<dxout<<", dtout="<<dtout <<endl
       <<" nxout="<<nxout<<", ntout="<<ntout
       <<endl;
  FILE  *outfile;                  
  //  filecheck (outfile,fname);
  outfile = fopen(fname,"w");

  fprintf(outfile, 
    "# x(km) \t t(min) \t rho(1/km) \t v(km/h) \t Q(1/h) \t wMatrix\n");

  for (int it=0; it<=ntout; it++)
  {
  
    for (int j=0; j<=nxout; j++)
    {
      double x_km          = (xmin+(xmax-xmin)*j/nxout)/1000.;
      double t_h         = (tmin+(tmax-tmin)*it/ntout)/3600.;
      double rho_invkm     = array1[j][it];
      double v_kmh         = array2[j][it];
      double Q_invh        = array3[j][it];
      fprintf(outfile,  "%.6f\t  %.6f\t %.6f\t  %.6f\t  %.6f\t  %f\n",
                      x_km, t_h, rho_invkm, v_kmh, Q_invh, array4[j][it]);
    }
    fprintf(outfile, "\n");        // newline mean new t step for gnuplot
  }

  fclose(outfile);

} // InOut::write_array2d





//############################################################
//  function getvar
//############################################################

// interprets lines beginning with COMMENTCHAR as comments;
// reads the first word of other lines;
// inteprets the remaining contents also of these lines as comments
//!!! geht nicht mit referenz-Werten (F... File lesen)

void InOut::getvar(FILE *fp, double *pdouble)
{
  char   str[LINEMAX];               // Dummy string for comments 
  //  const char COMMENTCHAR = '%';
  if (fgets(str,LINEMAX,fp))
  {
    while ( ( (str[0]==COMMENTCHAR)||(str[0]==COMMENTCHAR2)||(str[0]=='\\'))
        && fgets(str,LINEMAX,fp)){
    }
    sscanf(str,"%lf", pdouble);
    printf("%.3f\t %s",*pdouble,str);
  }
} // InOut::getvar


void InOut::getvar(FILE *fp, int *pint)
{
  char   str[LINEMAX];    
  const char COMMENTCHAR = '%';

  if (fgets(str,LINEMAX,fp))
  {
    while ( ( (str[0]==COMMENTCHAR)||(str[0]==COMMENTCHAR2)||(str[0]=='\\'))
        && fgets(str,LINEMAX,fp)){
    }
    sscanf(str,"%i",pint);
    printf("%i\t %s",*pint,str);
  }
}


void InOut::getvar(FILE *fp, long *plong)
{
  char   str[LINEMAX];    
  const char COMMENTCHAR = '%';

  if (fgets(str,LINEMAX,fp))
  {
    while ( ( (str[0]==COMMENTCHAR)||(str[0]==COMMENTCHAR2)||(str[0]=='\\'))
        && fgets(str,LINEMAX,fp)){
    }
    sscanf(str,"%li",plong);
    printf("%li\t %s",*plong,str);
  }
}


// transform "csv" into spaces
void InOut::prepare_line(char line[]){
  for (int j=0; j<LINEMAX; j++){
    if((line[j]==';')||(line[j]==',')){
      line[j]=' ';
    }
  }
}

void InOut::prepare_line(string & line){
  for (int j=0; j<int(line.length()); j++){
    if((line[j]==';')||(line[j]==',')){
      line[j]=' ';
    }
  }
}

int InOut::is_not_white(char line[], int nchar)
{
  int not_white=false;
  for (int i=0; i<=nchar-1; i++)  // last always '\n'
  {
     // cout << " line["<<i<<"]="<<line[i] << "nonwhite = "
      //    <<( ( (line[i]!=' ') && (line[i] !='\n') )? 1 : 0) << endl;
    if( (line[i]!=' ') && (line[i] !='\n') && (line[i] !='\t') ) not_white=true;
  }
  return(not_white);
}

bool InOut:: is_data_line(char line[]){
  return  (line[0] != COMMENTCHAR) && 
    (line[0] != COMMENTCHAR2) && 
    (strlen(line)>0) ;
}

bool InOut:: is_data_line(string line){
  return  (line[0] != COMMENTCHAR) && 
    (line[0] != COMMENTCHAR2) && 
    (!isalpha(line[0])) && // first character of line is not a letter [hack!]
    (line.length()>0) ;
}




void  InOut::write_array2d (char* fname, double dxout, double dtout, 
    double xmin, double xmax, double tmin, double tmax,
    double array1[][NYMAX+1], double array2[][NYMAX+1],
			    double array3[][NYMAX+1], double array4[][NYMAX+1],
			    bool data3d)
			   

{

  int nxout = (int)((xmax-xmin)/dxout)+1; // write_array2d (multivariate data)
  int ntout = (int)((tmax-tmin)/dtout)+1; // write_array2d (multivariate data)
  if((nxout>NXMAX)||(ntout>NYMAX)){
    cerr <<"InOut::write_array2d:" 
         <<" nxout="<<nxout<<">=NXMAX="<<NXMAX<<" and/or ntout="<<ntout
	 <<" >=NYMAX="<<NYMAX<<" => Error!" << endl;
    exit(-1);
  }

  cout << "in write_array2d: "<<endl
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<", xmax="<<xmax <<endl
       <<" tmin="<<tmin<<", tmax="<<tmax <<endl
       <<" dxout="<<dxout<<", dtout="<<dtout <<endl
       <<" nxout="<<nxout<<", ntout="<<ntout
       <<endl;
  FILE  *outfile;                  
  //  filecheck (outfile,fname);
  outfile = fopen(fname,"w");

  if(data3d){
    fprintf(outfile, "# x(inp units) \t t(inp units) \t rho \t v(inp units) \t Q\t wMatrix\n");
  }
  else{
      fprintf(outfile, "# x(km) \t t(h) \t\t rho(1/km) \t v(km/h) \t Q(1/h) \t wMatrix\n");
  }
  
  // (martin nov07) no longer case distinction necessary
  
  for (int it=0; it<ntout; it++){
     for (int j=0; j<nxout; j++) {
       double x          = (xmin+(xmax-xmin)*j/(nxout-1));
       double t        = (tmin+(tmax-tmin)*it/(ntout-1));
       double rho    = array1[j][it];
       double v         = array2[j][it];
       double Q       = array3[j][it];
       fprintf(outfile,  "%.6f\t  %.6f\t %.6f\t  %.6f\t  %.6f\t  %f\n",
                      x, t, rho, v, Q, array4[j][it]);
     }

     fprintf(outfile, "\n");        // newline mean new t step for gnuplot
  }

  fclose(outfile);

}


