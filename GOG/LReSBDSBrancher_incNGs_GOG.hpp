
//***********************************************************************//

#include <gecode/int/branch.hh>
#include <gecode/int.hh>
#include <gecode/GOG/GOG.hpp>

namespace Gecode {
 


	//the following are used for cs_sbds

	int s_index;
	IntArgs MVH;
	IntArgs Aa,Bb,Cc,Dd,E,F,G,H,reF;
	IntArgs WatNum,VarWNum,FailNum;
	int ChoH;
	int nNum, nSym,_myc,_psize,nSymparcom;
	int dosize;
	IntArgs PartA; 
	IntArgs MyA,MyB;
	int size1;
	int size2;
	IntArgs ReA;
	int e_f_myb;
	int e_p_myb;
	int ct;

	Index_class _symmetries (int id, int index, int val){
	   
		 if(reF[id]==0)
		 {	if(val==Aa[id]) val=Bb[id];
			else if(val==Bb[id]) val=Aa[id];}
		 else  
		 if(reF[id]==1)
		 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
			else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;}
		 else if(reF[id]==2)
		 {	if(index%nNum==Aa[id]-1) index=index+Bb[id]-Aa[id];
			else if(index%nNum==Bb[id]-1) index=index-Bb[id]+Aa[id];}
		 else if (reF[id]==3)
		 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
			else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
			if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
			else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
		 }
		 else if (reF[id]==4)
		 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
			else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
			if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
			else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
			if(index>=nNum*(E[id]-1)&&index<nNum*E[id]) index=index+(F[id]-E[id])*nNum;
			else if(index>=nNum*(F[id]-1)&&index<nNum*F[id]) index=index-(F[id]-E[id])*nNum;
		 }
		 else if (reF[id]==5)
		 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
			else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
			if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
			else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
			if(index%nNum==G[id]-1) index=index+H[id]-G[id];
			else if(index%nNum==H[id]-1) index=index-H[id]+G[id];
		 }
		 else if (reF[id]==6)
		 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
			else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
			if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
			else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
			if(index>=nNum*(E[id]-1)&&index<nNum*E[id]) index=index+(F[id]-E[id])*nNum;
			else if(index>=nNum*(F[id]-1)&&index<nNum*F[id]) index=index-(F[id]-E[id])*nNum;
			if(index%nNum==G[id]-1) index=index+H[id]-G[id];
			else if(index%nNum==H[id]-1) index=index-H[id]+G[id];
		 }
		 return Index_class(index,val); 
	}  
	  


//**************************************************//
//**************************************************//
//**************************************************//
//**************************************************//
//**************************************************//
//*************   branch ***************************//



 template<class View, int n, class Val, unsigned int a>
class LReSBDSBrancher : public ViewValBrancher<View,n,Val,a> {
    typedef typename ViewBrancher<View,n>::BranchFilter BranchFilter;
  public:
     
	int start;
 
	int p_size;
	
	int p_myb;
 
  protected:
    /// Constructor for cloning \a b
    LReSBDSBrancher(Space& home, bool share, LReSBDSBrancher& b);
    /// Constructor for creation
	
	
    LReSBDSBrancher(Home home, 
                 ViewArray<View>& x, 
                 ViewSel<View>* vs[n], 
                 ValSelCommitBase<View,Val>* vsc,
                 BranchFilter bf,IntVarValPrint vvp);
  public:
    /// Return choice
    virtual const Choice* choice(Space& home);
    /// Return choice
    virtual const Choice* choice(const Space& home, Archive& e);
    /// Perform commit for choice \a c and alternative \a b
    virtual ExecStatus commit(Space& home, const Choice& c, unsigned int b);
    /// Perform cloning
    virtual Actor* copy(Space& home, bool share);
    /// Perform dispose
    virtual size_t dispose(Space& home);
    /// Delete brancher and return its size
    static BrancherHandle post(Home home,
                               ViewArray<View>& x, 
                               ViewSel<View>* vs[n],
                               ValSelCommitBase<View,Val>* vsc,
                               BranchFilter bf,IntVarValPrint vvp);
  };

  template<class View, int n, class Val, unsigned int a> 
  LReSBDSBrancher<View,n,Val,a>
  ::LReSBDSBrancher(Home home, ViewArray<View>& x, 
                 ViewSel<View>* vs[n],
                 ValSelCommitBase<View,Val>* vsc,
                 BranchFilter bf,IntVarValPrint vvp)
    : ViewValBrancher<View,n,Val,a>(home, x, vs, vsc, bf,vvp) 
  {		
		
		PartA=IntArgs(2*x.size());
		p_size=0;
		start=0;
		p_myb=0;
		home.notice(*this, AP_DISPOSE);
		 
	  
  }

  template<class View, int n, class Val, unsigned int a>
  forceinline BrancherHandle
  LReSBDSBrancher<View,n,Val,a>::
  post(Home home, ViewArray<View>& x, 
       ViewSel<View>* vs[n], ValSelCommitBase<View,Val>* vsc,
       BranchFilter bf,IntVarValPrint vvp) {
		return *new (home) LReSBDSBrancher<View,n,Val,a>(home,x,vs,vsc,bf,vvp);
  }

  template<class View, int n, class Val, unsigned int a>
  forceinline
  LReSBDSBrancher<View,n,Val,a>::
  LReSBDSBrancher(Space& home, bool shared, LReSBDSBrancher<View,n,Val,a>& b)
    : ViewValBrancher<View,n,Val,a>(home,shared,b) {
	  
		//std::cout<<"the mystack is "<<b.SymObject.mystack<<"\n";
		 
		 
		p_size=b.p_size;
		start=b.start;
		p_myb=b.p_myb;
		p_myb=e_p_myb;
		//y.update(home,shared,b.y); 
		//std::cout<<"copy over\n";
		

  }
  
  template<class View, int n, class Val, unsigned int a>
  Actor*
  LReSBDSBrancher<View,n,Val,a>::copy(Space& home, bool shared) {
    //std::cout<<"copy now! \n";
    return new (home) LReSBDSBrancher<View,n,Val,a>(home,shared,*this);
  }


  // Compute choice
  template<class View, int n, class Val, unsigned int a>
  const Choice*
  LReSBDSBrancher<View,n,Val,a>::choice(Space& home) { //std::cout<<"choice \n";
	return ViewValBrancher<View,n,Val,a>::choice(home);
 
  }

 template<class View, int n, class Val, unsigned int a>
  const Choice*
  LReSBDSBrancher<View,n,Val,a>::choice(const Space& home, Archive& e) {

  
    return ViewValBrancher<View,n,Val,a>::choice(home,e);
  } 
template<class View, int n, class Val, unsigned int a>
  size_t
  LReSBDSBrancher<View,n,Val,a>::dispose(Space& home) {
    home.ignore(*this,AP_DISPOSE);
    (void) ViewValBrancher<View,n,Val,a>::dispose(home);
    return sizeof(LReSBDSBrancher<View,n,Val,a>);
  }
  template<class View, int n, class Val, unsigned int a>
  ExecStatus
  LReSBDSBrancher<View,n,Val,a>
  ::commit(Space& home, const Choice& c, unsigned int b) {
  
		const PosValChoice<Val>& pvc
		  = static_cast<const PosValChoice<Val>&>(c);
		int pos = pvc.pos().pos;
		int val = pvc.val();
	   
		_myc=1-b;
 
		if(b==0)
		{      	
			p_myb=e_p_myb;
			 
			 //ajust symmetries
			 PartA[p_size*2]=pos;
			 PartA[p_size*2+1]=val;
			 p_size++;
			 _psize=p_size;
			 e_p_myb=p_myb;
			 
			 GECODE_ME_CHECK(this->x[pos].eq(home,val));
			 StatusStatistics c;
			 home.status(c); 
			 
			   
        }
		else
		{		
			//symmetry breaking constraint

			MyB[p_myb++]=pos; MyB[p_myb++]=val; MyB[p_myb++]=p_size-1;
			 
			 
			e_p_myb=p_myb;
			 
			_psize=p_size;
			GECODE_ME_CHECK(this->x[pos].nq(home,val));
			StatusStatistics c;
			home.status(c);  
 
        }
		
  

		return ES_OK;
	}
  
//**************************************************//
//**************************************************//
//**************************************************//
//**************************************************//
//**************************************************//
//*************   post ***************************//

  
 BrancherHandle
  branch(Home home, const IntVarArgs& x,
         IntVarBranch vars, IntValBranch vals,
         int n,int m,int dosize0,IntBranchFilter bf=NULL) {
    using namespace Int;
    if (home.failed()) return BrancherHandle();
    ViewArray<IntView> xv(home,x);
    ViewSel<IntView>* vs[1] = { 
      Branch::viewselint(home,vars) 
    };
	int a=n*(n-1)/2;
	int b=m*(m-1)/2;

	nSym=n*m*n*m;
	nSymparcom=n-1+m-1;
	std::cout<<nSym<<"\t"<<n<<","<<m<<"\n";
	nNum=m;
	 Aa=IntArgs(nSym);
	 Bb=IntArgs(nSym);
	 Cc=IntArgs(nSym);
	 Dd=IntArgs(nSym);
	 reF=IntArgs(nSym);
	 
	 
 
	IntArgs _A(a+b);
	IntArgs _B(a+b);
	
 
	int r=0;
	int r1=0;
	for(int i=0;i<n-1;i++)
	{
	  int j=i+1;
	  for(int j=i+1;j<n;j++)
      {
        _A[r1]=i+1;
        _B[r1++]=j+1;

	  }
	  Aa[r]=i+1;
	  Bb[r]=i+2;
	 
	  reF[r++]=1;
    }
	
	int re0=r;

    for(int i=0;i<m-1;i++)
    {  
	   int j=i+1;
	   for(int j=i+1;j<m;j++)
       {
         _A[r1]=i+1;
         _B[r1++]=j+1;
       } 
	     Aa[r]=i+1;
         Bb[r]=i+2;
         reF[r++]=2;
	}
	  
	for(int i=0;i<a;i++)
		for(int j=0;j<b;j++)
		{ 	Aa[r]=_A[i];
			Bb[r]=_B[i];
			Cc[r]=_A[a+j];
			Dd[r]=_B[a+j];
		 
			reF[r++]=3;
		 
		}
 
 
    nSym=r;  
 
	std::cout<<"------------nSym is "<<nSym<<"\n";
	s_index=x.size()*dosize0;
	dosize=dosize0;

	//post global constraints
	ReA=IntArgs(x.size()*dosize0*2);
	size1=6;
	//create MyA,MyB
	MyA=IntArgs((nSym*(x.size()+1))*size1);
	MyB=IntArgs(s_index*3);
	 
	for(int i=0;i<nSym;i++)
	{	  
		MyA[i*size1+0]=1;//whether the symmetry is broken or subsumed or not
		MyA[i*size1+1]=i;//the index of the symmetry
		MyA[i*size1+2]=-1;//first pointer.pos
		MyA[i*size1+3]=-1;//first pointer.val
		MyA[i*size1+4]=-1;//the position of first pointer
		MyA[i*size1+5]=0;//the position of rhs assignment
	}
		
	gog(home,xv,r); 
 
		
	//build a columwise heuristic
	MVH=IntArgs(x.size());
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
		 	MVH[i*n+j]=j*m+i;
 
    	 
    return LReSBDSBrancher<IntView,1,int,2>::post
        (home,xv,vs,Branch::valselcommitint(home,x.size(),vals),bf,NULL);
    }
}


  

