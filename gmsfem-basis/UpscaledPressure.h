#ifndef UPSCALEDPRESSURE_H
#define UPSCALEDPRESSURE_H

#include <dolfin.h>
#include <fstream>

using namespace dolfin;

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

double readDoubleNumbers(const std::string &s, std::vector<double> &v){
    std::istringstream is(s);
    double n;
    while( is >> n ) 
        v.push_back( n );
    return v.size();
}

int readIntNumbers(const std::string &s, std::vector<int> &v){
    std::istringstream is(s);
    int n;
    while( is >> n ) 
        v.push_back( n );
    return v.size();
} 

void loadVector(std::string fileVec, std::vector<double> &vals){
    std::ifstream infile(fileVec.data());
    std::string line;
    while( getline(infile,line) ){
        std::vector<std::string> vec = split(line, ' ');
        int dof = std::atoi(vec[0].data());
        double val = std::atof(vec[1].data());
        vals.push_back(val);
    }
    infile.close();
}

class MyMatrix {
  public:
    std::vector<std::vector<int> > rowconn;
    MyMatrix(std::string fileMat, int sizeA){
        info(" ====== LOAD MATRIX %d ====== ", sizeA);
        rowconn.resize(sizeA);
        MatCreateSeqAIJ(PETSC_COMM_SELF, sizeA, sizeA, sizeA, 0, &Ap);
        std::string line;
        std::ifstream inMat(fileMat.data());
        while( getline(inMat, line) ){
            std::vector<std::string> vec = split(line, ' ');
            int i = std::atoi(vec[0].data());
            int j = std::atoi(vec[1].data());
            double val = std::atof(vec[2].data());
            // set value
            MatSetValues(Ap, 1, &i, 1, &j ,&val, INSERT_VALUES);
            // connection by matrix
            rowconn[i].push_back(j);
        }
        inMat.close();
        MatAssemblyBegin(Ap, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Ap, MAT_FINAL_ASSEMBLY);
    }

    double get(int i, int j){
        double val;
        MatGetValues(Ap,  1, &i, 1, &j, &val);
        return val;
    }

    Mat& getMat(){
        return Ap;
    }

 private:
    Mat Ap;
};

#endif