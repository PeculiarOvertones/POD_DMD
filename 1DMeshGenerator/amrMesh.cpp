#include "main.h"

  int main() {
         
     c_Mesh mesh;
     
     mesh.allocateMesh();
     mesh.initializeRootCells();
//     mesh.performAMRProcedure();
     mesh.outputCartesianMeshInTecplot();
     mesh.outputCartesianConnectivityInTecplot();
//     mesh.countActiveCellsInTheMesh();
//     mesh.outputAMRMeshInTecplot();
     mesh.deallocateMesh();
     
     return 0;

  }


  
  void c_Mesh:: allocateMesh() {
  
     for (int i=0; i<DIM; i++) {
         meshMin[i] = GLO_MeshMin[i];
         meshMax[i] = GLO_MeshMax[i];
     }
     
     nRoots = 1;

     for (int i=0; i<DIM; i++) {
         meshNodesInDir[i] = GLO_Nodes[i];
         nRoots *= (meshNodesInDir[i]-1);
     }

     cout << "nRoots: " << nRoots << endl;
     
     #ifdef ThreeD
     root = new c_Cell_Root **[meshNodesInDir[0]-1];
     
     for (int i=0; i<meshNodesInDir[0]-1; i++) {
         root[i] = new c_Cell_Root *[meshNodesInDir[1]-1];
     }

     for (int i=0; i<meshNodesInDir[0]-1; i++) {
         for (int j=0; j<meshNodesInDir[1]-1; j++) {
             root[i][j] = new c_Cell_Root [meshNodesInDir[2]-1];
         }
     }

     #elif TwoD
     root = new c_Cell_Root *[meshNodesInDir[0]-1];
     
     for (int i=0; i<meshNodesInDir[0]-1; i++) {
         root[i] = new c_Cell_Root [meshNodesInDir[1]-1];
     }
     #endif
     
     cout << "Roots are allocated." << '\n';
   
  } 
  
  
  
  void c_Mesh:: deallocateMesh() {
  
     #ifdef ThreeD

     for (int i=0; i<meshNodesInDir[0]-1; i++){
         for (int j=0; j<meshNodesInDir[1]-1; j++){
             delete [] root[i][j];
             root[i][j]=NULL;
         }
     }
     
     for (int i=0; i<meshNodesInDir[0]-1; i++){
         delete [] root[i];
         root[i]=NULL;
     }
     #elif TwoD

     for (int i=0; i<meshNodesInDir[0]-1; i++){
         for (int j=0; j<meshNodesInDir[1]-1; j++){
             delete [] root[i];
             root[i]=NULL;
         }
     }

     #endif
     delete [] root;
     root=NULL;
     cout << "Roots are deallocated." << '\n';
   
  }


  
  void c_Mesh:: initializeRootCells() {
  
     double dx[DIM];
     for (int i=0; i<DIM; i++){
     dx[i] = (meshMax[i] - meshMin[i]) / double(meshNodesInDir[i]-1);
     }
     
     int dum_cellID = 0;
     int dum_nodeID = 0;
     
     #ifdef ThreeD
     for (int i=0; i<meshNodesInDir[0]-1; i++){
         for (int j=0; j<meshNodesInDir[1]-1; j++){
             for (int k=0; k<meshNodesInDir[2]-1; k++){

                 root[i][j][k]._myID=  dum_cellID;
                 dum_cellID++;
     
     
                 root[i][j][k].node[0].nodeID = dum_nodeID;
                 dum_nodeID++;
                 root[i][j][k].node[0].nodePos[0] = meshMin[0] + dx[0] * double(i);
                 root[i][j][k].node[0].nodePos[1] = meshMin[1] + dx[1] * double(j);
                 root[i][j][k].node[0].nodePos[2] = meshMin[2] + dx[2] * double(k);
       
                 root[i][j][k].node[1].nodeID = dum_nodeID;
                 dum_nodeID++;
                 root[i][j][k].node[1].nodePos[0] = meshMin[0] + dx[0] * double(i);
                 root[i][j][k].node[1].nodePos[1] = meshMin[1] + dx[1] * double(j);
                 root[i][j][k].node[1].nodePos[2] = meshMin[2] + dx[2] * double(k+1);
       
                 root[i][j][k].node[2].nodeID = dum_nodeID;
                 dum_nodeID++;
                 root[i][j][k].node[2].nodePos[0] = meshMin[0] + dx[0] * double(i);
                 root[i][j][k].node[2].nodePos[1] = meshMin[1] + dx[1] * double(j+1);
                 root[i][j][k].node[2].nodePos[2] = meshMin[2] + dx[2] * double(k);
       
                 root[i][j][k].node[3].nodeID = dum_nodeID;
                 dum_nodeID++;
                 root[i][j][k].node[3].nodePos[0] = meshMin[0] + dx[0] * double(i);
                 root[i][j][k].node[3].nodePos[1] = meshMin[1] + dx[1] * double(j+1);
                 root[i][j][k].node[3].nodePos[2] = meshMin[2] + dx[2] * double(k+1);
       
                 root[i][j][k].node[4].nodeID = dum_nodeID;
                 dum_nodeID++;
                 root[i][j][k].node[4].nodePos[0] = meshMin[0] + dx[0] * double(i+1);
                 root[i][j][k].node[4].nodePos[1] = meshMin[1] + dx[1] * double(j);
                 root[i][j][k].node[4].nodePos[2] = meshMin[2] + dx[2] * double(k);
       
                 root[i][j][k].node[5].nodeID = dum_nodeID;
                 dum_nodeID++;
                 root[i][j][k].node[5].nodePos[0] = meshMin[0] + dx[0] * double(i+1);
                 root[i][j][k].node[5].nodePos[1] = meshMin[1] + dx[1] * double(j);
                 root[i][j][k].node[5].nodePos[2] = meshMin[2] + dx[2] * double(k+1);
       
                 root[i][j][k].node[6].nodeID = dum_nodeID;
                 dum_nodeID++;
                 root[i][j][k].node[6].nodePos[0] = meshMin[0] + dx[0] * double(i+1);
                 root[i][j][k].node[6].nodePos[1] = meshMin[1] + dx[1] * double(j+1);
                 root[i][j][k].node[6].nodePos[2] = meshMin[2] + dx[2] * double(k);
       
                 root[i][j][k].node[7].nodeID = dum_nodeID;
                 dum_nodeID++;
                 root[i][j][k].node[7].nodePos[0] = meshMin[0] + dx[0] * double(i+1);
                 root[i][j][k].node[7].nodePos[1] = meshMin[1] + dx[1] * double(j+1);
                 root[i][j][k].node[7].nodePos[2] = meshMin[2] + dx[2] * double(k+1);
       
             }  
         }    
     }
     #elif TwoD

     for (int i=0; i<meshNodesInDir[0]-1; i++){
         for (int j=0; j<meshNodesInDir[1]-1; j++){

             root[i][j]._myID=  dum_cellID;
             dum_cellID++;
             root[i][j]._myLevel =  0;
             root[i][j]._active =  true;
             root[i][j]._parent =  NULL;
      
             root[i][j].node = new c_Node[(int)pow(2,DIM)];
      
             root[i][j].node[0].nodeID = dum_nodeID;
             dum_nodeID++;
             root[i][j].node[0].nodePos[0] = meshMin[0] + dx[0] * double(i);
             root[i][j].node[0].nodePos[1] = meshMin[1] + dx[1] * double(j);
      
             root[i][j].node[1].nodeID = dum_nodeID;
             dum_nodeID++;
             root[i][j].node[1].nodePos[0] = meshMin[0] + dx[0] * double(i);
             root[i][j].node[1].nodePos[1] = meshMin[1] + dx[1] * double(j+1);
      
             root[i][j].node[2].nodeID = dum_nodeID;
             dum_nodeID++;
             root[i][j].node[2].nodePos[0] = meshMin[0] + dx[0] * double(i+1);
             root[i][j].node[2].nodePos[1] = meshMin[1] + dx[1] * double(j);
      
             root[i][j].node[3].nodeID = dum_nodeID;
             dum_nodeID++;
             root[i][j].node[3].nodePos[0] = meshMin[0] + dx[0] * double(i+1);
             root[i][j].node[3].nodePos[1] = meshMin[1] + dx[1] * double(j+1);
      
        }    
     }
     
     #endif
     cout << "Roots are initialized." << '\n';
  
  }

  
  
  void c_Mesh:: getIJKfromRootCellID(int id, int* index){

     #ifdef ThreeD
     index[1] = int(id / ((meshNodesInDir[0]-1)*(meshNodesInDir[2]-1) ));
    
     index[0] = id % (meshNodesInDir[0]-1);
    
     index[2] = (id - index[1] * (meshNodesInDir[2]-1)*(meshNodesInDir[0]-1) - index[0])/
              (meshNodesInDir[2]-1);

     #elif TwoD
     
     index[0] = int(id / (meshNodesInDir[1]-1));
     index[1] = id % (meshNodesInDir[1]-1);
     
     #endif
   
  }


  
  
  void c_Mesh:: performAMRProcedure() {

     #ifdef ThreeD
     for (int i=0;i<meshNodesInDir[0]-1;i++){
        for (int j=0;j<meshNodesInDir[1]-1;j++){
           for (int k=0;k<meshNodesInDir[1]-1;k++){
               if (root[i][j][k]._myID == 0) {
                  root[i][j][k].createAMRCells();
               }
           }
        }
     } 
     #elif TwoD
     for (int i=0;i<meshNodesInDir[0]-1;i++){
        for (int j=0;j<meshNodesInDir[1]-1;j++){
            if (root[i][j]._myID == 0) {
               root[i][j].createAMRCells();
            }
        }
     } 
     #endif
     cout << "AMR Procedure is performed." << endl;
  }
  
  
  void c_Cell_Root:: createAMRCells() {

     evaluateCellBasedOnRefinementCri();

  }
  
  
  bool c_Cell_AMR:: refinementFlag(){

     if (_myLevel < 2 ) {
     return true;
     }
  
     else if (_myID == 0 && _myLevel < 4) {
     return true;
     }
     
     return false;

  }

  
  
  void c_Cell_AMR:: evaluateCellBasedOnRefinementCri(){

     if (refinementFlag()){
        refineCell();
     }
     }
     
     
     void c_Cell_AMR:: refineCell(){
     _active = false;
     _children = new c_Cell_AMR[(int)pow(2,DIM)];
     
     for (int i=0; i<(int)pow(2,DIM); i++){
         createChild(_children[i], i);
     }
     
     for (int i=0; i<(int)pow(2,DIM); i++){
         _children[i].evaluateCellBasedOnRefinementCri();
     }
     
  }

    
  
  void c_Cell_AMR:: createChild(c_Cell_AMR& child, int child_id){

     child._myLevel = _myLevel + 1;
     child._myID = child_id;
     child._parent = this;
     #ifdef ThreeD
     double len_scale[3];
     
     len_scale[0] = double((child_id)/4) * 0.5;
     
     if((child_id==0)||(child_id==1)||(child_id==4)||(child_id==5)){
       len_scale[1] = 0.0;
     }
     else{
       len_scale[1] = 0.5;
     }
     
     len_scale[2] = (child_id)%2 * 0.5;
     
     child.node[0].nodePos[0] = node[0].nodePos[0] + len_scale[0] *
     (node[4].nodePos[0]-node[0].nodePos[0]);
     child.node[0].nodePos[1] = node[0].nodePos[1] + len_scale[1] *
     (node[2].nodePos[1]-node[0].nodePos[1]);
     child.node[0].nodePos[2] = node[0].nodePos[2] + len_scale[2] *
     (node[1].nodePos[2]-node[0].nodePos[2]);
     
     child.node[1].nodePos[0] = node[1].nodePos[0] + len_scale[0] *
     (node[5].nodePos[0]-node[1].nodePos[0]);
     child.node[1].nodePos[1] = node[1].nodePos[1] + len_scale[1] *
     (node[3].nodePos[1]-node[1].nodePos[1]);
     child.node[1].nodePos[2] = node[0].nodePos[2] + (len_scale[2]+0.5) *
     (node[1].nodePos[2]-node[0].nodePos[2]);
     
     child.node[2].nodePos[0] = node[2].nodePos[0] + len_scale[0] *
     (node[6].nodePos[0]-node[2].nodePos[0]);
     child.node[2].nodePos[1] = node[0].nodePos[1] + (len_scale[1]+0.5) *
     (node[2].nodePos[1]-node[0].nodePos[1]);
     child.node[2].nodePos[2] = node[2].nodePos[2] + len_scale[2] *
     (node[3].nodePos[2]-node[2].nodePos[2]);
     
     child.node[3].nodePos[0] = node[3].nodePos[0] + len_scale[0] *
     (node[7].nodePos[0]-node[3].nodePos[0]);
     child.node[3].nodePos[1] = node[1].nodePos[1] + (len_scale[1]+0.5) *
     (node[3].nodePos[1]-node[1].nodePos[1]);
     child.node[3].nodePos[2] = node[2].nodePos[2] + (len_scale[2]+0.5) *
     (node[3].nodePos[2]-node[2].nodePos[2]);
     
     child.node[4].nodePos[0] = node[0].nodePos[0] + (len_scale[0]+0.5) *
     (node[4].nodePos[0]-node[0].nodePos[0]);
     child.node[4].nodePos[1] = node[4].nodePos[1] + len_scale[1] *
     (node[6].nodePos[1]-node[4].nodePos[1]);
     child.node[4].nodePos[2] = node[4].nodePos[2] + len_scale[2] *
     (node[5].nodePos[2]-node[4].nodePos[2]);
     
     child.node[5].nodePos[0] = node[1].nodePos[0] + (len_scale[0]+0.5) *
     (node[5].nodePos[0]-node[1].nodePos[0]);
     child.node[5].nodePos[1] = node[5].nodePos[1] + len_scale[1] *
     (node[7].nodePos[1]-node[5].nodePos[1]);
     child.node[5].nodePos[2] = node[4].nodePos[2] + (len_scale[2]+0.5) *
     (node[5].nodePos[2]-node[4].nodePos[2]);
     
     child.node[6].nodePos[0] = node[2].nodePos[0] + (len_scale[0]+0.5) *
     (node[6].nodePos[0]-node[2].nodePos[0]);
     child.node[6].nodePos[1] = node[4].nodePos[1] + (len_scale[1]+0.5) *
     (node[6].nodePos[1]-node[4].nodePos[1]);
     child.node[6].nodePos[2] = node[6].nodePos[2] + len_scale[2] *
     (node[7].nodePos[2]-node[6].nodePos[2]);
     
     child.node[7].nodePos[0] = node[3].nodePos[0] + (len_scale[0]+0.5) *
     (node[7].nodePos[0]-node[3].nodePos[0]);
     child.node[7].nodePos[1] = node[5].nodePos[1] + (len_scale[1]+0.5) *
     (node[7].nodePos[1]-node[5].nodePos[1]);
     child.node[7].nodePos[2] = node[6].nodePos[2] + (len_scale[2]+0.5) *
     (node[7].nodePos[2]-node[6].nodePos[2]);
     
     #elif TwoD
     double len_scale[2];
     
     if((child_id==0)||(child_id==1)){
       len_scale[0] = 0.0;
     }
     else{
       len_scale[0] = 0.5;
     }
     
     if((child_id==0)||(child_id==2)){
       len_scale[1] = 0.0;
     }
     else{
       len_scale[1] = 0.5;
     }
     
     
     child.node[0].nodePos[0] = node[0].nodePos[0] + len_scale[0] *
     (node[2].nodePos[0]-node[0].nodePos[0]);
     child.node[0].nodePos[1] = node[0].nodePos[1] + len_scale[1]*
     (node[1].nodePos[1]-node[0].nodePos[1]);
     
     child.node[1].nodePos[0] = node[1].nodePos[0] + len_scale[0] *
     (node[3].nodePos[0]-node[1].nodePos[0]);
     child.node[1].nodePos[1] = node[0].nodePos[1] + (len_scale[1]+0.5)*
     (node[1].nodePos[1]-node[0].nodePos[1]);
     
     child.node[2].nodePos[0] = node[0].nodePos[0] + (len_scale[0]+0.5) *
     (node[2].nodePos[0]-node[0].nodePos[0]);
     child.node[2].nodePos[1] = node[2].nodePos[1] + len_scale[1] *
     (node[3].nodePos[1]-node[2].nodePos[1]);
     
     child.node[3].nodePos[0] = node[1].nodePos[0] + (len_scale[0]+0.5)*
     (node[3].nodePos[0]-node[1].nodePos[0]);
     child.node[3].nodePos[1] = node[2].nodePos[1] + (len_scale[1]+0.5)*
     (node[3].nodePos[1]-node[2].nodePos[1]);
     
     #endif

  }


  
  void c_Mesh:: countActiveCellsInTheMesh() {

     totalActiveCellsInTheMesh = 0;
     int eachNumberOfCells;
     #ifdef ThreeD
     for (int i=0;i<meshNodesInDir[0]-1;i++){
        for (int j=0;j<meshNodesInDir[1]-1;j++){
           for (int k=0;k<meshNodesInDir[1]-1;k++){
               eachNumberOfCells = 0;
               totalActiveCellsInTheMesh += root[i][j][k].countOwnActiveCells(eachNumberOfCells);     
               cout << "ActiveCells: " << totalActiveCellsInTheMesh << endl;  
           }
        }
     }
     #elif TwoD
     for (int i=0;i<meshNodesInDir[0]-1;i++){
        for (int j=0;j<meshNodesInDir[1]-1;j++){
               eachNumberOfCells = 0;
               totalActiveCellsInTheMesh += root[i][j].countOwnActiveCells(eachNumberOfCells);     
               cout << "ActiveCells: " << totalActiveCellsInTheMesh << endl;  
        }
     }
     #endif  

  }
  
  
  
  int c_Cell_AMR:: countOwnActiveCells(int &cellCounter) {
  
     if (_active == true){
             cellCounter++;
     }
     
     else {
            for(int i=0; i<(int)pow(2,DIM); i++){
               _children[i].countOwnActiveCells(cellCounter);
            }
     }
     
     return cellCounter;
   
  }


  void c_Mesh:: outputAMRMeshInTecplot() {
  
     int cellNum = totalActiveCellsInTheMesh;
     int nodePerCell = (int)pow(2,DIM);//Hexahedron
     int nodeNum = nodePerCell * cellNum;
     int varSize = 3;
     
     
     ofstream outmyfile;
     #ifdef ThreeD
     outmyfile.open("my3DAMRMesh.tec");
     outmyfile << "title = \"SUGAR TECPLOT\"" << '\n';
     outmyfile << "variables = \"x(m)\", \"y(m)\", \"z(m)\", \"Active Flag\", \"Id\", \"Depth\"" << '\n';
     outmyfile << " zone T = \"Mesh\", n=" << nodeNum << " ,e=" << cellNum << " ,DATAPACKING=BLOCK, ZONETYPE=FEBRICK" << '\n';
     outmyfile << " VARLOCATION = ([4-6] = CELLCENTERED)" << '\n';
     outmyfile << '\n';
     #elif TwoD
     outmyfile.open("my2DAMRMesh.tec");
     outmyfile << "title = \"SUGAR TECPLOT\"" << '\n';
     outmyfile << "variables = \"x(m)\", \"y(m)\",  \"Active Flag\", \"Id\", \"Depth\"" << '\n';
     outmyfile << " zone T = \"Mesh\", n=" << nodeNum << " ,e=" << cellNum << " ,DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << '\n';
     outmyfile << " VARLOCATION = ([3-5] = CELLCENTERED)" << '\n';
     outmyfile << '\n';
     #endif
     
        c_GridStoringAMR *AMR = new c_GridStoringAMR(totalActiveCellsInTheMesh, varSize);
        AMR->_Mesh = this;
        AMR->setNodes();
        AMR->setValues();
     
        int **cells=NULL;
        cells = new int *[cellNum];
        for(int i =0; i<cellNum; i++){
          cells[i] = new int [nodePerCell];
        }
     
        for(int i=0;i<cellNum;i++){
          for(int j=0;j<nodePerCell;j++){
            cells[i][j]=AMR->cellCon[nodePerCell*i+j];
          }
        }
     
        outmyfile << '\n';
        for(int i =0; i<nodeNum; i++){
          outmyfile << AMR->x[i] << " ";
          if((i%1000 == 0) && (i != 0)){
            outmyfile << '\n';
          }
        }
        outmyfile << '\n';
        for(int i =0; i<nodeNum; i++){
          outmyfile << AMR->y[i] << " ";
          if((i%1000 == 0) && (i != 0)){
            outmyfile << '\n';
          }
        }
        outmyfile << '\n';
     #ifdef ThreeD
        for(int i =0; i<nodeNum; i++){
          outmyfile << AMR->z[i] << " ";
          if((i%1000 == 0) && (i != 0)){
            outmyfile << '\n';
          }
        }
        outmyfile << '\n';
     #endif
        for(int i =0; i<varSize; i++){
          for(int j =0; j<cellNum; j++){
            outmyfile << AMR->cellData[i][j] << " ";
            if((j%1000 == 0) && (j != 0)){
            outmyfile << '\n';
            }
          }
          outmyfile << '\n';
        }
     
        outmyfile << '\n';
        for(int i =0; i<cellNum; i++){
     #ifdef ThreeD
          outmyfile << cells[i][0]+1 << " ";
          outmyfile << cells[i][1]+1 << " ";
          outmyfile << cells[i][3]+1 << " ";
          outmyfile << cells[i][2]+1 << " ";
          outmyfile << cells[i][4]+1 << " ";
          outmyfile << cells[i][5]+1 << " ";
          outmyfile << cells[i][7]+1 << " ";
          outmyfile << cells[i][6]+1 << " ";
          outmyfile << '\n';
     #elif TwoD
          outmyfile << cells[i][0]+1 << " ";
          outmyfile << cells[i][1]+1 << " ";
          outmyfile << cells[i][3]+1 << " ";
          outmyfile << cells[i][2]+1 << " ";
          outmyfile << '\n';
     #endif
     
        }
       for (int i=0;i<cellNum;i++){
         delete [] cells[i];
         cells[i]=NULL;
       }
       delete [] cells;
       cells=NULL;
     
       delete AMR;
       AMR = NULL;
     
     outmyfile.close();
     cout << "AMR mesh is outputted for the Mesh." << '\n';
  
  }

  void c_Mesh:: outputCartesianConnectivityInTecplot(){
     int cellNum = nRoots;
     int nodePerCell = (int)pow(2,DIM);
     int nodeNum = nodePerCell * cellNum;
     
     ofstream outmyfile;
     outmyfile.open("connectivity.dat");
     int **cells=NULL;
     cells = new int *[cellNum];
     for (int i =0; i<cellNum; i++){
         cells[i] = new int [nodePerCell];
     }
     
     int dumNodeCounterConn=0;
     for (int i=0;i<cellNum;i++){
         int index[DIM];
         getIJKfromRootCellID(i, index);
         for (int j=0;j<nodePerCell;j++){
             cells[i][j]=dumNodeCounterConn;
             dumNodeCounterConn++;
         }
     }
     
      outmyfile << '\n';
     for (int i =0; i<cellNum; i++){
     #ifdef ThreeD
          outmyfile << cells[i][0]+1 << " ";
          outmyfile << cells[i][1]+1 << " ";
          outmyfile << cells[i][3]+1 << " ";
          outmyfile << cells[i][2]+1 << " ";
          outmyfile << cells[i][4]+1 << " ";
          outmyfile << cells[i][5]+1 << " ";
          outmyfile << cells[i][7]+1 << " ";
          outmyfile << cells[i][6]+1 << " ";
          outmyfile << '\n';
     #elif TwoD
          outmyfile << cells[i][0]+1 << " ";
          outmyfile << cells[i][1]+1 << " ";
          outmyfile << cells[i][3]+1 << " ";
          outmyfile << cells[i][2]+1 << " ";
          outmyfile << '\n';
     #endif
     }
     

     for (int i=0;i<cellNum;i++){
         delete [] cells[i];
         cells[i]=NULL;
     }
     delete [] cells;
     cells=NULL;
     cout << "connectivity is outputted! " << endl;
     outmyfile.close();

  }


  void c_Mesh:: outputCartesianMeshInTecplot() {
  
     int cellNum = nRoots;
     int nodePerCell = (int)pow(2,DIM);
     int nodeNum = nodePerCell * cellNum;

     double dx[DIM];
     for (int i=0; i<DIM; i++){
     dx[i] = (meshMax[i] - meshMin[i]) / double(meshNodesInDir[i]-1);
     }


	 int ICells = GLO_Nodes[0]-1;
	 int JCells = GLO_Nodes[2]-1;
	 int KCells = GLO_Nodes[1]-1;
     
     ofstream outmyfile;
     #ifdef ThreeD
     outmyfile.open("my1DMesh.dat");
     outmyfile << "title = \"SUGAR TECPLOT\"" << '\n';
     outmyfile << "variables = \"X\", \"Y\", \"Z\", \"Var\"" << '\n';
//     outmyfile << " zone T = \"Mesh\", n=" << nodeNum << " ,e=" << cellNum << " ,DATAPACKING=BLOCK, ZONETYPE=FEBRICK" << '\n';
     outmyfile << " zone T = \"Mesh\", I= "<< ICells << ", J=" << JCells << ", K="<< KCells << ", DATAPACKING=BLOCK" << '\n';
     outmyfile << " VARLOCATION = ([4] = CELLCENTERED)" << '\n';
     outmyfile << '\n';
     #elif TwoD
     outmyfile.open("my2DMesh.dat");
     outmyfile << "title = \"SUGAR TECPLOT\"" << '\n';
     outmyfile << "variables = \"X\", \"Z\", \"ID\", \"Var\"" << '\n';
//     outmyfile << " zone T = \"Mesh\", n=" << nodeNum << " ,e=" << cellNum << " ,DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << '\n';
     outmyfile << " zone T = \"Mesh\", I= "<< ICells << ", J=" << JCells << ", DATAPACKING=BLOCK" << '\n';
     outmyfile << " VARLOCATION = ([3-4] = CELLCENTERED)" << '\n';
     outmyfile << '\n';
     #endif
//unordered data FEQUADRILATERAL     
//     for (int p=0; p<DIM; p++) {   
//         for (int i=0; i<cellNum; i++){
//             int index[DIM];
//             getIJKfromRootCellID(i, index);
//
//             for (int j=0; j<nodePerCell; j++ ){
//     #ifdef ThreeD
//                 outmyfile << setw(15) << root[index[0]][index[1]][index[2]].node[j].nodePos[p];
//     #elif TwoD
//                 outmyfile << setw(15) << root[index[0]][index[1]].node[j].nodePos[p];
//     #endif
//             }
//          if((i%100 == 0) && (i != 0)){
//            outmyfile << '\n';
//          }
//         }
//         outmyfile << '\n';
//         outmyfile << '\n';
//     }
//Ordered data
     for (int p=0; p<DIM; p++) {   
		int counter =0;
		int index;
		for (int k=0; k<KCells; k++){
			for (int j=0; j<JCells; j++){
				for (int i=0; i<ICells; i++){
					if(p==0){
						index=i;
					}
					else if (p==1){
						index=k;
					}
					else if(p==2){
						index=j;
					}
					outmyfile << setw(15) << dx[p]*index;
          			if((counter%100 == 0) && (counter != 0)){
          			  outmyfile << '\n';
          			}
					counter++;
				}
			}
		}
         outmyfile << '\n';
         outmyfile << '\n';
     }
     
//     for (int i=0; i<cellNum; i++){
//             int index[DIM];
//             getIJKfromRootCellID(i, index);
//     #ifdef ThreeD
//         outmyfile << setw(15) << root[index[0]][index[1]][index[2]]._myID;
//     #elif TwoD
//         outmyfile << setw(15) << root[index[0]][index[1]]._myID;
//     #endif
//          if((i%100 == 0) && (i != 0)){
//            outmyfile << '\n';
//          }
//     }
//     outmyfile << '\n';
//     outmyfile << '\n';
//     int **cells=NULL;
//     cells = new int *[cellNum];
//     for (int i =0; i<cellNum; i++){
//         cells[i] = new int [nodePerCell];
//     }
//     
//     int dumNodeCounterConn=0;
//     for (int i=0;i<cellNum;i++){
//         int index[DIM];
//         getIJKfromRootCellID(i, index);
//         for (int j=0;j<nodePerCell;j++){
//             cells[i][j]=dumNodeCounterConn;
//             dumNodeCounterConn++;
//         }
//     }
//     
//      outmyfile << '\n';
//     for (int i =0; i<cellNum; i++){
//     #ifdef ThreeD
//          outmyfile << cells[i][0]+1 << " ";
//          outmyfile << cells[i][1]+1 << " ";
//          outmyfile << cells[i][3]+1 << " ";
//          outmyfile << cells[i][2]+1 << " ";
//          outmyfile << cells[i][4]+1 << " ";
//          outmyfile << cells[i][5]+1 << " ";
//          outmyfile << cells[i][7]+1 << " ";
//          outmyfile << cells[i][6]+1 << " ";
//          outmyfile << '\n';
//     #elif TwoD
//          outmyfile << cells[i][0]+1 << " ";
//          outmyfile << cells[i][1]+1 << " ";
//          outmyfile << cells[i][3]+1 << " ";
//          outmyfile << cells[i][2]+1 << " ";
//          outmyfile << '\n';
//     #endif
//     }
//     
//
//     for (int i=0;i<cellNum;i++){
//         delete [] cells[i];
//         cells[i]=NULL;
//     }
//     delete [] cells;
//     cells=NULL;
//     
//     
     cout << "Tecplot files are outputted for the Mesh." << '\n';
     
     outmyfile.close();

  }



  
  void c_GridStoringAMR::setNodes(){
  
     #ifdef ThreeD
          int counter_nodes = 0;
          for(int i=0; i<_Mesh->nRoots; i++){
             int index[DIM];
             _Mesh->getIJKfromRootCellID(i, index);
             if(_Mesh->root[index[0]][index[1]][index[2]]._active == true){
               for(int j=0; j<8; j++){
                 x[counter_nodes] = _Mesh->root[index[0]][index[1]][index[2]].node[j].nodePos[0];
                 y[counter_nodes] = _Mesh->root[index[0]][index[1]][index[2]].node[j].nodePos[1];
                 z[counter_nodes] = _Mesh->root[index[0]][index[1]][index[2]].node[j].nodePos[2];
                 cellCon[counter_nodes] = counter_nodes;
                 counter_nodes++;
               }
             }
            else{
             for(int j=0; j<8; j++){
                recursiveSetNodes(_Mesh->root[index[0]][index[1]][index[2]]._children[j], counter_nodes);
             }
            }
          }
     
     #elif TwoD
          int counter_nodes = 0;
          for(int i=0; i<_Mesh->nRoots; i++){
             int index[DIM];
             _Mesh->getIJKfromRootCellID(i, index);
             if(_Mesh->root[index[0]][index[1]]._active == true){
               for(int j=0; j<(int)pow(2,DIM); j++){
                 x[counter_nodes] = _Mesh->root[index[0]][index[1]].node[j].nodePos[0];
                 y[counter_nodes] = _Mesh->root[index[0]][index[1]].node[j].nodePos[1];
                 cellCon[counter_nodes] = counter_nodes;
                 counter_nodes++;
               }
             }
            else{
             for(int j=0; j<(int)pow(2,DIM); j++){
                recursiveSetNodes(_Mesh->root[index[0]][index[1]]._children[j], counter_nodes);
             }
            }
          }
     
     #endif

  }


  
  void c_GridStoringAMR::recursiveSetNodes(c_Cell_AMR &cell, int &counter_nodes){
  
          if(cell._active == true){
            for(int j=0; j<(int)pow(2,DIM); j++){
              x[counter_nodes] = cell.node[j].nodePos[0];
              y[counter_nodes] = cell.node[j].nodePos[1];
     #ifdef ThreeD
              z[counter_nodes] = cell.node[j].nodePos[2];
     #endif
              cellCon[counter_nodes] = counter_nodes;
              counter_nodes++;
            }
          }
          else{
            for(int j=0; j<(int)pow(2,DIM); j++){
              recursiveSetNodes(cell._children[j], counter_nodes);
            }
          }
   
  }

  
  
  void c_GridStoringAMR::setValues(){

     #ifdef ThreeD
          int counter_cells = 0;
          for(int i=0; i<_Mesh->nRoots; i++){
             int index[DIM];
             _Mesh->getIJKfromRootCellID(i, index);
             if(_Mesh->root[index[0]][index[1]][index[2]]._active == true){
              cellData[0][counter_cells] =_Mesh->root[index[0]][index[1]][index[2]]._active;
              cellData[1][counter_cells] =_Mesh->root[index[0]][index[1]][index[2]]._myID;
              cellData[2][counter_cells] =_Mesh->root[index[0]][index[1]][index[2]]._myLevel;
              counter_cells++;
             }
            else{
             for(int j=0; j<(int)pow(2,DIM); j++){
                recursiveSetValues(_Mesh->root[index[0]][index[1]][index[2]]._children[j], counter_cells);
             }
            }
          }
     #elif TwoD
          int counter_cells = 0;
          for(int i=0; i<_Mesh->nRoots; i++){
             int index[DIM];
             _Mesh->getIJKfromRootCellID(i, index);
             if(_Mesh->root[index[0]][index[1]]._active == true){
              cellData[0][counter_cells] =_Mesh->root[index[0]][index[1]]._active;
              cellData[1][counter_cells] =_Mesh->root[index[0]][index[1]]._myID;
              cellData[2][counter_cells] =_Mesh->root[index[0]][index[1]]._myLevel;
              counter_cells++;
             }
            else{
             for(int j=0; j<(int)pow(2,DIM); j++){
                recursiveSetValues(_Mesh->root[index[0]][index[1]]._children[j], counter_cells);
             }
            }
          }
     
     #endif

  }


  
  void c_GridStoringAMR::recursiveSetValues(c_Cell_AMR &cell, int &counter_cells){
  
      if (cell._active == true){
           cellData[0][counter_cells] = cell._active;
           cellData[1][counter_cells] = cell._myID;
           cellData[2][counter_cells] = cell._myLevel;
         counter_cells++;
  
      }
      else{
         for(int j=0; j<(int)pow(2,DIM); j++){
           recursiveSetValues(cell._children[j], counter_cells);
         }
      }

  }
  
