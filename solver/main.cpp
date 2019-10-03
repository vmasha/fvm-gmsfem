#include <dolfin.h>
#include <fstream>
#include "quadmesh.h"
#include "Pressure.h"

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

double kr(double p){
    return std::exp(-0.1 * std::fabs(p));
}

// ./solver F 50 1.0e-4 80 80 ../data/modelF/ ../data/out/ 1 ./ err.txt
// ./solver C 50 1.0e-4 80 80 ../data/modelF/ ../data/out/ 1936 ../data/modelMs/R100 err.txt
int main(int argc, char** argv){
    parameters["reorder_dofs_serial"] = false;
    parameters["allow_extrapolation"] = true;

    Timer timer("solverFC");
    timer.start();

    /*
    * types:
    *   F - fine
    *   C - coarse
    */
    std::string type = argv[1];

    int timeIt = std::atoi(argv[2]); // 200
    double ta = std::atof(argv[3]);
    double Nx = std::atoi(argv[4]);
    int Ny = std::atoi(argv[5]);
    std::string indir = argv[6];
    std::string outdir = argv[7];

    int Nc = std::atoi(argv[8]); // 1024
    std::string infileR = argv[9];
    std::string errfilename = argv[10];

    // system
    std::string infileT = indir + "mat-K.txt";
    std::string infileRhs = indir + "rhs.txt";
    
    // time
    double timeMax = ta*timeIt; 
    double dt = timeMax/timeIt;
    int tdel = timeIt/50;
    info("Time = %f and dt = %f", timeMax, dt);
    
    // fine grid parameters
    double Lx = 1.0, Ly = 1.0;
    double hh = Lx/Nx; 
    double volK = hh*hh;
    int n = Nx*Ny;
    info("FINE SIZE %d x %d =  %d", Nx, Ny, n);

    // system parameters
    double ap = 1.0;
    double uinit = 0.0;
    
    // fine mesh and space for saving
    Mesh mesh2; 
    build(mesh2, Nx, Ny, 0, Lx, 0, Ly);
    auto mesh = std::make_shared<Mesh>(mesh2);
    mesh->init(); info(*mesh);
    auto W = std::make_shared<Pressure::FunctionSpace>(mesh);
    Function up(W), up2(W);
    File filep(outdir + "p.pvd"), filep2(outdir + "p2.pvd");
    
    // load T
    Mat PT;
    MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 5, 0, &PT);
    std::string line;
    info("load T from %s", infileT.data());
    std::ifstream inMat(infileT.data());
    while( getline(inMat, line) ){
        std::vector<std::string> vec = split(line, ' ');
        int i = std::atoi(vec[0].data());
        int j = std::atoi(vec[1].data());
        double val = std::atof(vec[2].data());
        MatSetValues(PT, 1, &i, 1, &j ,&val, INSERT_VALUES);
    }
    inMat.close();
    MatAssemblyBegin(PT, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(PT, MAT_FINAL_ASSEMBLY);
    info("loaded T %d-%d", n, n);

    // load Rhs
    PETScVector q(PETSC_COMM_SELF, n);
    info("load Rhs from %s", infileRhs.data());
    std::ifstream inRhs(infileRhs.data());
    while( getline(inRhs, line) ){
        std::vector<std::string> vec = split(line, ' ');
        int i = std::atoi(vec[0].data());
        double val = std::atof(vec[1].data());
        q.setitem(i, val); 
    }
    q.apply("insert");
    inRhs.close();
    info("loaded Rhs q in (%g, %g)", q.min(), q.max());

    // load R
    Mat PR, PRT;
    if(type == "C"){
        MatCreateSeqAIJ(PETSC_COMM_SELF, Nc, n, n, 0, &PR);
        info("load R from %s", infileR.data());
        std::ifstream inMatR(infileR.data());
        while( getline(inMatR, line) ){
            std::vector<std::string> vec = split(line, ' ');
            int i = std::atoi(vec[0].data());
            int j = std::atoi(vec[1].data());
            double val = std::atof(vec[2].data());
            MatSetValues(PR, 1, &i, 1, &j ,&val, INSERT_VALUES);
        }
        inMatR.close();
        MatAssemblyBegin(PR, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(PR, MAT_FINAL_ASSEMBLY);
        MatTranspose(PR, MAT_INITIAL_MATRIX, &PRT);
        info("COARSE: loaded R %d-%d", Nc, n);
    }

    // init solver
    LinearSolver solverF("default", "default");  
    LinearSolver solverC("default", "default");  
    PETScVector x(PETSC_COMM_SELF, n);
    PETScVector xms(PETSC_COMM_SELF, n);
    PETScVector b(PETSC_COMM_SELF, n);
    PETScVector b2(PETSC_COMM_SELF, n);
    Mat PA, PA2;
    // Mat PRA, PRART;
    MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 5, 0, &PA);
    MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 5, 0, &PA2);
    Vec pxms, pbc;
    VecCreateSeq(PETSC_COMM_SELF, n, &pxms);
    VecCreateSeq(PETSC_COMM_SELF, Nc, &pbc);
    info("System init %d", n);
    PETScVector bc(PETSC_COMM_SELF, Nc);
    PETScVector xc(PETSC_COMM_SELF, Nc);
    info("System init %d", Nc);

    // initial conditions
    for(int I = 0; I < n; I++){
        x.setitem(I, uinit);
        xms.setitem(I, uinit);
    }
    x.apply("insert"); 
    xms.apply("insert");
    
    // info("TIME Pinit %g", timer1.stop());
    // timer1.start();

    remove(errfilename.data()); 
    std::ofstream errfile(errfilename, std::ios_base::app);

    PetscInt ncols;
    const PetscInt *cols;
    const PetscReal *vals;
    double val, fval;     

    double t = 0;
    int tcounter = 0;
    while(t <= timeMax && tcounter <= timeIt){
        info("\n Time[%d] %g, dt = %g", tcounter, t, dt);
        
        // ----- PRESSURE FINE -----   
        for(int I = 0; I < n; I++){
            MatGetRow(PT, I, &ncols, &cols, &vals);
            double diagval = 0, diagval2 = 0;
            for(int ej = 0; ej < ncols; ej++){
                int J = cols[ej];
                if(I != J){ // off diag
                    // f
                    val = vals[ej]*( kr(x.getitem(I)) + kr(x.getitem(J)) )/2;
                    MatSetValues(PA, 1, &I, 1, &J ,&val, INSERT_VALUES);
                    diagval += val;
                    // ms 
                    val = vals[ej]*( kr(xms.getitem(I)) + kr(xms.getitem(J)) )/2;
                    MatSetValues(PA2, 1, &I, 1, &J ,&val, INSERT_VALUES);
                    diagval2 += val;
                }
            }
            // f
            val = -diagval + ap/dt*volK;
            MatSetValues(PA, 1, &I, 1, &I ,&val, INSERT_VALUES);
            fval = q.getitem(I) + ap/dt*volK*x.getitem(I);
            b.setitem(I, fval);
            // ms
            val = -diagval2 + ap/dt*volK;
            MatSetValues(PA2, 1, &I, 1, &I ,&val, INSERT_VALUES);
            fval = q.getitem(I) + ap/dt*volK*xms.getitem(I);
            b2.setitem(I, fval);
        }
        MatAssemblyBegin(PA, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(PA, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(PA2, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(PA2, MAT_FINAL_ASSEMBLY);
        b.apply("insert"); 
        b2.apply("insert"); 
        info("A [%d, %d] update", n, n);
        // solve
        PETScMatrix *A = new PETScMatrix(PA);  
        solverF.solve(*A, x, b);
        info("solve: p in (%g, %g)", x.min(), x.max());
        delete A;

        // ----- PRESSURE COARSE -----   
        if(type == "C"){
            Mat PRA, PRART;
            MatMatMult(PR, PA2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &PRA);
            MatMatMult(PRA, PRT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &PRART);
            MatMult(PR, b2.vec(), pbc);
            // solve
            PETScMatrix *A2 = new PETScMatrix(PRART);  
            for(int I = 0; I < Nc; I++){
                VecGetValues(pbc, 1, &I, &val);
                bc.setitem(I, val);
            }
            bc.apply("insert");
            info("rhsC in (%g, %g)", bc.min(), bc.max());
            solverC.solve(*A2, xc, bc);
            delete A2;
            MatDestroy(&PRA);
            MatDestroy(&PRART);
            // project 
            MatMult(PRT, xc.vec(), pxms);
            for(int I = 0; I < n; I++){
                VecGetValues(pxms, 1, &I, &val);
                xms.setitem(I, val);
            }
            xms.apply("insert");
            info("COARSE: solve p in (%g, %g)", xms.min(), xms.max());
        }
        
        // error and save
        double errvt = 0.0, sumvt = 0.0;
        double errp = 0.0, sump = 0.0;
        double x1, x2;
        for(int I = 0; I < n; I++){
            // f
            x1 = x.getitem(I);
            up.vector().get()->setitem(I, x1);
            // ms
            if(type == "C"){
                x2 = xms.getitem(I);
                up2.vector().get()->setitem(I, x2);
                // error
                errp += (x1 - x2)* (x1 - x2); 
                sump += x1*x1;
                // velocity error
                MatGetRow(PT, I, &ncols, &cols, &vals);
                for(int ej = 0; ej < ncols; ej++){
                    int I2 = cols[ej];
                    if(I != I2){ // off diag
                        if(type == "C"){
                            // err vel
                            double ft1 = vals[ej]*( kr(x.getitem(I)) + kr(x.getitem(I2)) )/2;
                            double ft2 = vals[ej]*( kr(xms.getitem(I)) + kr(xms.getitem(I2)) )/2;
                            double ut1 = -ft1*(x.getitem(I) - x.getitem(I2));
                            double ut2 = -ft2*(xms.getitem(I) - xms.getitem(I2));
                            errvt += (ut1 - ut2)*(ut1 - ut2);    
                            sumvt += ut1*ut1;
                        }
                    }
                }
                //
            }
        }
        up.vector().get()->apply("insert");  
        if(tcounter%tdel == 0){
            filep << up;
            info("saved");
        }
        if(type == "C"){
            up2.vector().get()->apply("insert"); 
            if(tcounter%tdel == 0){
                filep2 << up2;
                info("saved");
            }
            errp = std::sqrt(errp/sump);
            errvt = std::sqrt(errvt/sumvt);
            info("ERROR p %f, vel %f", errp*100, errvt*100);
            errfile << tcounter << " " << errp*100 << " " << errvt*100 << " " << "\n";
        }
        
        t += dt;
        tcounter++; 
    }
    errfile.close();
    info("SUCCESS %d time steps", tcounter);
    // PetscFinalize();
    return 0;
}
