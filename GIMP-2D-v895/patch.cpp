// Philip Wallstedt 2004-2009

#include "patch.h"
using namespace std;

// Cycle through each managed array and resize it
void patch::resizeParts(int s){
   for(unsigned v=0;v<parts.size();v+=1)parts[v]->resize(s);
   numberOfParts=s;
}

// Cycle through each managed array and resize it
void patch::resizeNodes(int s){
   for(unsigned v=0;v<nodes.size();v+=1)nodes[v]->resize(s);
   numberOfNodes=s;
}

// Copy a particle; not sure if I use this anywhere
int patch::copyPart(const patch& q,int i){
   for(unsigned v=0;v<parts.size();v+=1)parts[v]->push_back(q.parts[v],i);
   numberOfParts+=1;
   return Npart()-1;
}

patch::patch(const int Nx,
             const int Ny,
             const double bx,
             const double by,
             const double ex,
             const double ey,
             const int Ng,
             const double th):
               regionBegin(bx,by),
               regionEnd(ex,ey),
               dx((ex-bx)/double(Nx)),
               dy((ey-by)/double(Ny)),
               thick(th),
               Nghost(Ng),
               I(Nx+2*Ng+1),
               J(Ny+2*Ng+1){
   // taking control of the managed arrays
   parts=commonPartsPtr;commonPartsPtr.clear();numberOfParts=0;
   nodes=commonNodesPtr;commonNodesPtr.clear();
   resizeNodes(I*J); // lots of memory for nodes allocated here
   incCount=0;
   // laying out nodes in a grid
   for(int j=0;j<J;j+=1){
      const double y=(j-Ng)*dy+regionBegin.y;
      for(int i=0;i<I;i+=1){
         const double x=(i-Ng)*dx+regionBegin.x;
         const int n=j*I+i;
         gx[n]=Vector2(x,y);
      }
   }
}








