#include "quadmesh.h"
#include "Pressure.h"
#include "UpscaledPressure.h"

class ShapeFunction : public Expression{
    double vx[8];
    double vy[8];
    double vz[8];
    int mainind;
    double hx, hy, hz;
    int _dim;
 public:
    
    ShapeFunction(int dim, std::string filename):_dim(dim){
        // load cell borders 
        std::ifstream inlines(filename.data());
        std::string line;
        getline(inlines, line);
        std::vector<double> borders;
        readDoubleNumbers(line, borders);
        inlines.close();

        double minX = borders[0];
        double maxX = borders[1];
        hx = maxX - minX;
        double minY = borders[2];
        double maxY = borders[3];
        hy = maxY - minY;
        double minZ, maxZ;
        if(_dim == 3){
            minZ = borders[4];
            maxZ = borders[5];      
            hz = maxZ - minZ;
        }
        // z
        vx[0] = minX;  vy[0] = minY;  vz[0] = minZ;
        vx[1] = maxX;  vy[1] = minY;  vz[1] = minZ;
        vx[2] = minX;  vy[2] = maxY;  vz[2] = minZ;
        vx[3] = maxX;  vy[3] = maxY;  vz[3] = minZ;
        // z
        vx[4] = minX;  vy[4] = minY;  vz[4] = maxZ;
        vx[5] = maxX;  vy[5] = minY;  vz[5] = maxZ;
        vx[6] = minX;  vy[6] = maxY;  vz[6] = maxZ;
        vx[7] = maxX;  vy[7] = maxY;  vz[7] = maxZ;
    }

    void setMain(std::vector<double> vcoords, double meps){
        if(_dim == 2){
            for(int i = 0; i < 4; i++){
                // info("%d: %g, %g", i, vx[i], vy[i]);
                if( std::abs(vx[i] - vcoords[0]) < meps && std::abs(vy[i] - vcoords[1]) < meps )
                    mainind = i;
            }
            // info("%d: %g, %g", mainind, vcoords[0], vcoords[1]);
        }else{
            for(int i = 0; i < 8; i++){
                // info("%d: %g, %g, %g", i, vx[i], vy[i], vz[i]);
                if( std::abs(vx[i] - vcoords[0]) < meps && std::abs(vy[i] - vcoords[1]) < meps && std::abs(vz[i] - vcoords[2]) < meps )
                    mainind = i;
            }
            // info("%d: %g, %g, %g", mainind, vcoords[0], vcoords[1], vcoords[2]);
        }
        info("MAIN VERTEX FOR SHAPE FUNCTIOB %d (%g, %g, %g)", mainind, vx[mainind], vy[mainind], vz[mainind]);       
    }

    void eval(Array<double>& values, const Array<double>& x) const{
        values[0] = 0.0;
        bool inCell;
        if(_dim == 2)
            inCell = (x[0] >= vx[0] && x[0] <= vx[1]) && (x[1] >= vy[0] && x[1] <= vy[2]);
        else
            inCell = (x[0] >= vx[0] && x[0] <= vx[1]) && (x[1] >= vy[0] && x[1] <= vy[2]) && (x[2] >= vz[0] && x[2] <= vz[4]);
        if(!inCell)
            return;
        if(mainind == 0 || mainind == 4){
            values[0] =  (x[0] - vx[1])/hx * (x[1] - vy[3])/hy;
            if(_dim == 3 && mainind == 0)
                values[0] = -(x[0] - vx[1])/hx * (x[1] - vy[3])/hy*(x[2] - vz[4])/hz;
            if(_dim == 3 && mainind == 4)
                values[0] =  (x[0] - vx[1])/hx * (x[1] - vy[3])/hy*(x[2] - vz[0])/hz;
        }
        if(mainind == 1 || mainind == 5){
            values[0] = -(x[0] - vx[0])/hx * (x[1] - vy[2])/hy;
            if(_dim == 3 && mainind == 1)
                values[0] =  (x[0] - vx[0])/hx * (x[1] - vy[2])/hy*(x[2] - vz[5])/hz;
            if(_dim == 3 && mainind == 5)
                values[0] = -(x[0] - vx[0])/hx * (x[1] - vy[2])/hy*(x[2] - vz[1])/hz;
        }
        if(mainind == 2 || mainind == 6){
            values[0] = -(x[0] - vx[3])/hx * (x[1] - vy[1])/hy;
            if(_dim == 3 && mainind == 2)
                values[0] =  (x[0] - vx[3])/hx * (x[1] - vy[1])/hy*(x[2] - vz[6])/hz;
            if(_dim == 3 && mainind == 6)
                values[0] = -(x[0] - vx[3])/hx * (x[1] - vy[1])/hy*(x[2] - vz[2])/hz;
        }
        if(mainind == 3 || mainind == 7){
            values[0] =  (x[0] - vx[2])/hx * (x[1] - vy[0])/hy;
            if(_dim == 3 && mainind == 3)
                values[0] = -(x[0] - vx[2])/hx * (x[1] - vy[0])/hy*(x[2] - vz[7])/hz;
            if(_dim == 3 && mainind == 7)
                values[0] =  (x[0] - vx[2])/hx * (x[1] - vy[0])/hy*(x[2] - vz[3])/hz;
        }
    }
};

class QuadCell : public SubDomain {
 public:
    QuadCell(double minX, double minY, double maxX, double maxY, double eps): _eps(eps){
        _minX = minX;
        _maxX = maxX;
        _minY = minY;
        _maxY = maxY;        
        // info("Quad cell [%g, %g] - [%g, %g]", _minX, _maxX, _minY, _maxY);
    }

    bool inside(const Array<double> &x, bool on_boundary) const{
        if (x[0] < (_minX - _eps) || x[0] > (_maxX + _eps) || 
            x[1] < (_minY - _eps) || x[1] > (_maxY + _eps))
                return false;
        return true;
    }

 private:
    double _minX, _maxX, _minY, _maxY, _minZ, _maxZ;
    double _eps;
};

// ./demo_poisson 2 ../../data/3level-frac/out-10u/ ../../data/3level-frac/omega-5/ 34 1 0.025 40 8 ../../data/3level-frac/model-efm/
// ./demo_poisson 2 ../../data/3level-frac/out-10u/ ../../data/3level-frac/omega-5/ 34 1 0.025 40 8 ../../data/3level-frac/model-efm-40/
int main(int argc, char** argv){
    PetscInitialize(&argc, &argv, (char*)0, "");
    parameters["reorder_dofs_serial"] = false;
    parameters["allow_extrapolation"] = true;
    
    /*
    * type: 
    *   2 EFM
    *   3 NLMC
    */
    int type = std::atoi(argv[1]);
    std::string outDir = argv[2];
    std::string omDir = argv[3];
    int cind = std::atoi(argv[4]);

    int save0 = atoi(argv[5]);
    bool saveit = (save0 == 1)?true:false;

    std::string resultdir = outDir + "results/";

    std::string dofFile   = outDir + "dof/dof-" + std::to_string(cind);// + ".txt";
    std::string eigenFile = outDir + "eigen/" + std::to_string(cind) + "-m";
    std::string countFile = outDir + "eigen/basisCount-" + std::to_string(cind) + ".txt";
    std::string volFile   = outDir + "eigen/basisVol-" + std::to_string(cind) + ".txt";
    remove(dofFile.data());
    remove(eigenFile.data());
    remove(countFile.data());
    remove(volFile.data());

    // int NS = atoi(argv[6]);
    // double hh = atof(argv[6]);
    int NL = atoi(argv[7]);

    int ecount =  std::atoi(argv[8]);

    std::string modelDir = argv[9]; 

    // global domain
    double lx = atof(argv[6]);//1.0;
    double ly = lx;

     double hh = lx/NL;
    
    // -----------------
    // ----- OMEGA -----
    // -----------------
    std::string omfile = omDir + "w" + std::to_string(cind) + ".txt";
    std::ifstream inlines(omfile.data());
    std::string line;
    getline(inlines, line);
    std::vector<int> omcells;
    readIntNumbers(line, omcells);
    inlines.close();
    // get coords
    std::string omfile2 = omDir + "v" + std::to_string(cind) + ".txt";
    std::ifstream inlines2(omfile2.data());
    std::string line2;
    getline(inlines2, line2);
    std::vector<double> vcoords;
    readDoubleNumbers(line2, vcoords);
    inlines2.close();
    // all cells
    double minX = vcoords[0], minY = vcoords[1];
    double maxX = vcoords[0], maxY = vcoords[1];
    for (int ii = 0; ii < omcells.size(); ii++){
        int ci = omcells[ii];
        std::string cfile = omDir + "c" + std::to_string(ci) + ".txt";
        // load cell coords
        std::ifstream inlines(cfile.data());
        std::string line;
        getline(inlines, line);
        std::vector<double> borders;
        readDoubleNumbers(line, borders);
        inlines.close();
        // find local domain coord for QUAD cells
        double minXi = borders[0], minYi = borders[2];
        double maxXi = borders[1], maxYi = borders[3];
        minX = (minXi < minX)?minXi:minX;
        minY = (minYi < minY)?minYi:minY;
        maxX = (maxXi > maxX)?maxXi:maxX;
        maxY = (maxYi > maxY)?maxYi:maxY;
    }

    int NSx = (maxX - minX + 1.0e-10)/hh;
    int NSy = (maxY - minY+ 1.0e-10)/hh;
    info("MESH %g and %g", (maxX - minX)/hh,  (maxY - minY)/hh);
    info("MESH %d, %d", NSx, NSy);

    // -----------------    
    // ----- MESH ------
    // -----------------
    int facetCountInCell = 4;
    Mesh mesh2; 
    build(mesh2, NSx, NSy, minX, maxX, minY, maxY);
    auto mesh = std::make_shared<Mesh>(mesh2);
    mesh->init(); info(*mesh);
    if(saveit){
        File filePvD1(resultdir + "mesh-" + std::to_string(cind) + ".pvd");
        filePvD1 << *mesh;
    }
    // function
    auto W = std::make_shared<Pressure::FunctionSpace>(mesh);
    Function p(W);
    
    double myeps = mesh->hmax()*1.0e-5;
    int mdim = mesh->topology().dim();
    int Nc = mesh->num_cells();

    // -----------------------
    // --------- POU ---------
    // -----------------------
    std::vector<double> puvec(Nc);
    Function pu(W);
    Function pu2(W);
    // all cells
    for (int ii = 0; ii < omcells.size(); ii++){
        int ci = omcells[ii];
        std::string cfile = omDir + "c" + std::to_string(ci) + ".txt";
        ShapeFunction fpou(mdim, cfile);    
        fpou.setMain(vcoords, myeps);

        pu2.vector().get()->zero();
        pu2 = fpou;
        for(int di = 0; di < pu2.vector().get()->size(); di++){
            double val = pu2.vector().get()->getitem(di);
            if(val > 1.0e-10){
                pu.vector().get()->setitem(di, val);
                puvec[di] = val;
            }
        }
        pu.vector().get()->apply("insert");
    }
    info("pou (%g, %g)", pu.vector().get()->min(), pu.vector().get()->max());
    if(saveit){
        File fileMsbU(resultdir + "pou-" + std::to_string(cind) + ".pvd");
        fileMsbU << pu;
    }

    // -----------------    
    // --- CELL MAP ----
    // -----------------
    Mesh fmesh2; 
    build(fmesh2, NL, NL, 0, lx, 0, ly);
    auto finemesh = std::make_shared<Mesh>(fmesh2);
    finemesh->init();  info(*finemesh);
    // File filePvD2(resultdir + "meshf.pvd");
    // filePvD2 << *finemesh;
    MeshFunction<std::size_t> fineSundomForOmega(finemesh, mdim);
    fineSundomForOmega = 0;
    QuadCell qcell(minX, minY, maxX, maxY, myeps);
    qcell.mark(fineSundomForOmega, 1);
    int fNc = finemesh->num_cells();
    std::vector<int> cellmap(Nc);
    for(int i = 0; i < Nc; i++){
        Cell lcell(*mesh, i);
        Point lpoint = lcell.midpoint();
        for(int j = 0; j < fNc; j++){
            if(fineSundomForOmega[j] == 1){
                Cell gcell(*finemesh, j);
                double dd = lpoint.distance(gcell.midpoint());
                // info("dist %g, g%d (%g, %g) and l%d (%g, %g)", dd, j, gcell.midpoint().x(), gcell.midpoint().y(), i, lpoint.x(), lpoint.y());
                if(dd < myeps){
                    cellmap[i] = j;
                    // info("cell map %d - %d", i, j);
                }
            }
        }
    }
   
    // ------------------------------
    // ------ LOAD GLOBAL A, M ------
    // ------------------------------
    std::string filenameDOF = modelDir + "dof100";
    std::string filenameM   = modelDir + "mat-S.txt";
    std::string filenameK   = modelDir + "mat-K.txt";
    std::string filenameMdiag = modelDir + "mat-Sdiag.txt";
    info("LOAD DOF %s", filenameDOF.c_str());
    // DOF
    std::vector<double> rowvol;
    std::vector<int> rowcell;
    std::vector<int> rowcont;
    std::ifstream inDOF(filenameDOF.data());
    while( getline(inDOF,line) ){
        std::vector<std::string> vec = split(line, ' ');
        int dof = std::atoi(vec[0].data());
        int ci = std::atoi(vec[1].data());
        int mc = std::atoi(vec[2].data());
        double pi = std::atof(vec[3].data());
        double vol = std::atof(vec[4].data());
        rowvol.push_back(vol);
        rowcell.push_back(ci);
        rowcont.push_back(mc);
        // info("dof %d, ci %d, mc %d, p %g, vol %g", dof, ci, mc, pi, vol);
    }
    inDOF.close();
    int Ngl = rowvol.size();
    // MATRICES
    MyMatrix matK(filenameK, Ngl);
    MyMatrix matM(filenameM, Ngl);
    std::vector<double> rowmass;
    loadVector(filenameMdiag, rowmass);
    // CONNECTIONS
    std::vector<std::vector<int> > rowconnK = matK.rowconn;
    std::vector<std::vector<int> > rowconnM = matM.rowconn;
    info("LOAD GLOBAL DOF %d", Ngl);

    // ---------------------------
    // ------ EXTRACT LOCAL ------
    // ---------------------------
    std::vector<int> dofmap;
    std::vector<int> dofcell;
    std::vector<int> dofcont;
    std::vector<double> dofvol;
    for(int i = 0; i < Ngl; i++){
        int gci = rowcell[i];
        int mci = rowcont[i];
        double voli = rowvol[i];
        // info("global %d", gci);
        bool containcell = (std::find(cellmap.begin(), cellmap.end(), gci) != cellmap.end());
        if(containcell){
            // info("global222 %d", gci);
            dofmap.push_back(i);
            int lci = find(cellmap.begin(), cellmap.end(), gci) - cellmap.begin();
            dofcell.push_back(lci);
            dofcont.push_back(mci);
            dofvol.push_back(voli);
        }
    }
    int sizeA = dofmap.size();
    info("LOCAL DOF %d", sizeA);

    // ------------------------
    // ------ LOCAL A, M ------
    // ------------------------
    Mat Ap, Mp;
    int maxNZ = sizeA;
    MatCreateSeqAIJ(PETSC_COMM_SELF, sizeA, sizeA, maxNZ, 0, &Ap);
    MatCreateSeqAIJ(PETSC_COMM_SELF, sizeA, sizeA, 1, 0, &Mp);
    double val;        
    for(int li = 0; li < sizeA; li++){
        int gi = dofmap[li];
        double sumT = 0;
        std::vector<int> iconn = rowconnK[gi];
        for(int lj = 0; lj < sizeA; lj++){
            int gj = dofmap[lj];
            if( std::find(iconn.begin(), iconn.end(), gj) != iconn.end() ){
                if(gi != gj){
                    val = matK.get(gi, gj);
                    MatSetValues(Ap, 1, &li, 1, &lj ,&val, ADD_VALUES);
                    sumT += -val;
                }
            }
        } 
        val = sumT;
        MatSetValues(Ap, 1, &li, 1, &li ,&val, INSERT_VALUES);
        // mass
        val =  rowmass[gi];
        MatSetValues(Mp, 1, &li, 1, &li ,&val, INSERT_VALUES);
    }
    MatAssemblyBegin(Ap, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Ap, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Mp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Mp, MAT_FINAL_ASSEMBLY);
    info("create A, M[%d, %d]", sizeA, sizeA);

    // ----------------------------
    // --------- SPECTRAL BASIS ---
    // ----------------------------
    auto As = std::make_shared<PETScMatrix>(Ap);
    auto Ms = std::make_shared<PETScMatrix>(Mp);
    Mat PR, PRT;
    PETScVector uv(PETSC_COMM_SELF, sizeA);

    // Eigen SOLVER
    SLEPcEigenSolver *esolver = new SLEPcEigenSolver(MPI_COMM_SELF, As, Ms);
    esolver->parameters["solver"] = "krylov-schur";
    // esolver->parameters["problem_type"] = "gen_hermitian";
    // if(haveFrac){
        esolver->parameters["spectral_transform"] = "shift-and-invert";
        esolver->parameters["spectral_shift"] = 1.0e-5;//1.0e-8;//1.0e-5;
    // }
    esolver->parameters["spectrum"] = "smallest magnitude";// magnitude real
    esolver->solve(ecount+2);
    int ccount = esolver->get_number_converged();
    ecount = (ecount > ccount)?ccount:ecount;
    // info("ecount %d", ecount);

    double lmd, v;
    PETScVector rx, vx;   
    int eigcounter = 0;
    for(int i = 0; i < ecount; i++){
        esolver->get_eigenpair(lmd, v, rx, vx, i);
        info("lambda[%d] = %g, z in (%g, %g)", i, lmd, rx.min(), rx.max());
        // if(lmd > -1.0e-3)
        // if(lmd < 1.0e-1 && lmd > -1.0e-3)
            eigcounter++;
    }
    info("ECOUNT %d", ecount);
    ecount = eigcounter;
    info("ECOUNT %d", ecount);
    
    // save
    PetscInt arraydof[sizeA];
    for(int i = 0; i < sizeA; i++)
        arraydof[i] = dofmap[i];
    int fdd;
    PetscViewer viewerDof;
    PetscViewerBinaryOpen(PETSC_COMM_SELF, dofFile.data(), FILE_MODE_WRITE, &viewerDof);
    PetscViewerBinaryGetDescriptor(viewerDof, &fdd);
    PetscBinaryWrite(fdd, &arraydof, sizeA, PETSC_INT, PETSC_FALSE);
    PetscViewerDestroy(&viewerDof);
    // save bin
    int fd1;
    PetscViewer viewerSol1;
    PetscViewerBinaryOpen(PETSC_COMM_SELF,  eigenFile.data(), FILE_MODE_WRITE, &viewerSol1);
    PetscViewerBinaryGetDescriptor(viewerSol1, &fd1);
    // vol
    std::ofstream ofileVol(volFile, std::ios_base::app);   

    PetscScalar *array;
    File filep0(resultdir + "ev-m.pvd");
    int ccounter = 0;
    for(int i = 0; i < ecount; i++){
        double vol = 1.0;
        esolver->get_eigenpair(lmd, v, rx, vx, i);
        // if(lmd < 1.0e-1 && lmd > -1.0e-3){
        if(lmd > -1.0e-3){
            ofileVol << ccounter << " " << vol << "\n";
            info("lambda[%d] = %g, z in (%g, %g)", i, lmd, rx.min(), rx.max());

            p.vector().get()->zero();
            for(int jj = 0; jj < sizeA; jj++){
                int cjj = dofcell[jj];
                double mval = rx.getitem(jj);
                if(dofcont[jj] == 0)
                    p.vector().get()->setitem(cjj, mval);
                double rval =  puvec[cjj] * mval;
                rx.setitem(jj, rval);
                // ebasisFile << rval << std::endl;
            }
            rx.apply("insert");
            p.vector().get()->apply("insert");

            VecGetArray(rx.vec(), &array);
            PetscBinaryWrite(fd1, array, sizeA, PETSC_DOUBLE, PETSC_FALSE);

            if(saveit)
                filep0 << p;
            ccounter++;
        }
    }
    // ebasisFile.close();
    PetscViewerDestroy(&viewerSol1);    
    ofileVol.close();

    // count
    std::ofstream ofile(countFile, std::ios_base::app);   
    ofile << ccounter << std::endl << sizeA << "\n";
    ofile.close();
    
 
    info("SUCCESS %d", ccounter);
//    PetscFinalize();
    return 0;
}
