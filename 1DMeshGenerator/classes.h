/*Classes for AMR & Cartesian Mesh Generator*/



  class c_Node {

     public:
     double nodePos[DIM];
     int nodeID;

  };



  class c_Cell_AMR: public c_Node {

     public:
     c_Node* node;
     int _myID;
     bool _active;
     int _myLevel;
     c_Cell_AMR* _parent;
     c_Cell_AMR* _children;

     c_Cell_AMR() {
       _active = true;
       _parent = NULL;
       _children = NULL;
       node = new c_Node[(int)pow(2,DIM)];
     }

     ~c_Cell_AMR() {
     //   delete [] node;
       _parent = NULL;
       if (_children != NULL) {
         delete [] _children;
         _children = NULL;
       }
     }

     bool refinementFlag();
     void evaluateCellBasedOnRefinementCri();
     void refineCell();
     void createChild(c_Cell_AMR&, int);
     int countOwnActiveCells(int&);

  };
    


  class c_Cell_Root: public c_Cell_AMR {

     public:
     c_Cell_Root() {
       _myLevel =  0;
       _active =  true;
       _parent =  NULL;
       _children =  NULL;
       node = new c_Node[(int)pow(2,DIM)];
     }

     ~c_Cell_Root(){
       delete [] node;
       if (_children != NULL){
          delete [] _children;
          _children = NULL;
       }
     }

     void createAMRCells();

  };



  class c_GridStoringAMR{

     protected:
     int _numberOfCells;
     int _numberOfCellCenterValues;

     public:
     float *x;
     float *y;
     int *cellCon;
     float **cellData;
     class c_Mesh* _Mesh;


     #ifdef ThreeD
     float *z;
     #endif

     c_GridStoringAMR(int loc_numberOfCells, int numberOfCellCenterValues){

       _numberOfCells = loc_numberOfCells;
       _numberOfCellCenterValues = numberOfCellCenterValues;
       x = NULL;
       y = NULL;
       x = new float[(int)pow(2,DIM)*_numberOfCells];
       y = new float[(int)pow(2,DIM)*_numberOfCells];
     #ifdef ThreeD
       z = NULL;
       z = new float[(int)pow(2,DIM)*_numberOfCells];
     #endif
       cellCon = NULL;
       cellCon = new int[(int)pow(2,DIM)*_numberOfCells];
       cellData = NULL;
       cellData = new float *[_numberOfCellCenterValues];
       for (int i=0; i<_numberOfCellCenterValues; i++) {
           cellData[i] = new float [_numberOfCells];
       }

     }

     ~c_GridStoringAMR() {

       delete [] x;
       delete [] y;
       x = NULL;
       y = NULL;
     #ifdef ThreeD
       delete [] z;
       z = NULL;
     #endif
       delete [] cellCon;
       cellCon = NULL;

       for (int i=0;i<_numberOfCellCenterValues;i++){
           delete [] cellData[i];
           cellData[i]=NULL;
       }

       delete [] cellData;
       cellData=NULL;

     }

     void setNodes();
     void recursiveSetNodes(c_Cell_AMR&, int &);
     void setValues();
     void recursiveSetValues(c_Cell_AMR&, int &);

  };



  class c_Mesh: public c_Cell_Root {

     public:
     int meshNodesInDir[DIM];
     double meshMin[DIM];
     double meshMax[DIM];
     int nRoots;
     int totalActiveCellsInTheMesh;
     #ifdef ThreeD
     c_Cell_Root ***root = NULL;
     #elif TwoD
     c_Cell_Root **root = NULL;
     #endif

     void allocateMesh();
     void initializeRootCells();
     void performAMRProcedure();
     void getIJKfromRootCellID(int, int*);
     void outputCartesianMeshInTecplot();
     void outputCartesianConnectivityInTecplot();
     void countActiveCellsInTheMesh();
     void outputAMRMeshInTecplot();
     void deallocateMesh();

  };

