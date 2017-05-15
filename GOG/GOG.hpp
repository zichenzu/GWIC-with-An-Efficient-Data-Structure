#include <gecode/int/branch.hh>
#include <gecode/int.hh>
#include <gecode/GOG/head.hh>

namespace Gecode{  
 
  
	//extern IntArgs Aa,Bb,Cc,Dd,reF;

	extern int nNum, nSym, _myc, _psize,nSymparcom;
	extern IntArgs PartA; 
	extern IntArgs MyA,MyB;
	extern int size1,size2;
	extern IntArgs ReA;
	extern int e_p_myb;
	extern IntArgs WatNum,VarWNum,FailNum;
	extern int ChoH;
	extern int dosize;	
	template<class View>
	 class GOG : public Propagator {
	 protected: 
		ViewArray<View> x;

		int stack;
		int _nSym;//number of unbroken symmetries
		int _csize;//number of current size of MyA;
		int re;
		 
	 public: 
		//post
		GOG(Space& home,  ViewArray<View>& x0, int r);
		//copy
		GOG(Space& home, bool share, GOG<View>& p);
		virtual Propagator* copy(Space& home, bool share);
		//cost
		virtual PropCost cost(const Space&, const ModEventDelta&) const;
		//propagation
		virtual ExecStatus propagate(Space& home, const ModEventDelta&);
		//post
		static ExecStatus post(Space& home, ViewArray<View>& x0, int r);
		//dispose
		virtual size_t dispose(Space& home);
 
	};
 
	// posting
	template<class View>
	inline
	GOG<View>::GOG(Space& home,   ViewArray<View>& x0, int r) : Propagator(home), x(x0){
		//initialize 
		_nSym=_csize=r;
		stack=0;
		re=0;
		 
		x.subscribe(home,*this,Int::PC_INT_DOM); 
		 
	} 
	template<class View>
	ExecStatus
	GOG<View>::post(Space& home, ViewArray<View>& x0, int r) {
		(void) new (home) GOG(home,x0,r); 
		return ES_OK; 
	} 
	// disposal 
	template<class View>
	forceinline size_t
	GOG<View>::dispose(Space& home) {
		 
		x.cancel(home,*this,Int::PC_INT_DOM);
		 
		ChoH=-1;
		(void) Propagator::dispose(home);
		return sizeof(*this); 
	} 
	// copying
	template<class View>
	forceinline
	GOG<View>::GOG(Space& home, bool share, GOG<View>& p) : Propagator(home,share,p),
														_nSym(p._nSym),_csize(p._csize),stack(p.stack),re(p.re){
												 
		p.stack++;
		x.update(home,share,p.x); 
 
		int r=0;
		 
		for(int i=0;i<_csize;i++)
			if(MyA[(stack*nSym+i)*size1+0]!=-2)//not subsumed yet but false 
			{	
				 
				for(int j=0;j<size1;j++)
				{	
					MyA[((stack+1)*nSym+r)*size1+j]=MyA[(stack*nSym+i)*size1+j];
				}
				
				r++;
			}
	  
		 
		p._csize=r;
		 
	} 
	template<class View> Propagator*  GOG<View>::copy(Space& home, bool share) {
		return new (home) GOG(home,share,*this); 
	} 
	// cost computation 
	template<class View>
	PropCost
	GOG<View>::cost(const Space&, const ModEventDelta&) const {
		return PropCost::linear(PropCost::LO,int(100000000000));
	} 
	
	 
	// propagation 
	template<class View>
	ExecStatus
	GOG<View>::propagate(Space& home, const ModEventDelta&) { 
	 
	 
		int mre=re;
		re=0;
		if(_myc==1&&_nSym>0)//new assignments are added
		{	//std::cout<<"@ "<<_csize<<","<<stack*nSym+0<<":"<<MyA[stack*nSym+0][0][0]<<"\n";
			for(int i=0;i<_csize;i++)
				if(MyA[(stack*nSym+i)*size1+0]==1&& MyA[(stack*nSym+i)*size1+2]==-1 &&	MyA[(stack*nSym+i)*size1+3]==-1)//need to find the first pointer
				{	 
					int s;
					for(s=MyA[(stack*nSym+i)*size1+4]+1;s<_psize;s++)
					{
						 
						Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],PartA[s*2],PartA[1+s*2]);
						if(!x[g.index].in(g.val))//now the symmetry is broken
						{	MyA[(stack*nSym+i)*size1+0]=-2;//subsumed now
							 
							_nSym--;
							break;
						}
						if(x[g.index].in(g.val)&&!x[g.index].assigned())//now find the first pointer
						{
						
							MyA[(stack*nSym+i)*size1+2]=g.index;
							MyA[(stack*nSym+i)*size1+3]=g.val;
							MyA[(stack*nSym+i)*size1+4]=s; 
							 
							//add extra checks in order to prune the 
							int ns=0;
							int t;
							int p=MyA[(stack*nSym+i)*size1+5];
							for(t=s+1;t<_psize;t++) 
							{
								
								 	
								for(;p<e_p_myb&& MyB[p+2]<t;p+=3)
								{	
									Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],MyB[p],MyB[p+1]);
									if(x[g.index].in(g.val))
									{	
										ns=1;
										break;
									}
								}
								
								if(ns==1) break;
								 
								Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],PartA[t*2],PartA[1+t*2]);
								if(!x[g.index].in(g.val))//now the symmetry is broken
								{	
									ns=-1; 
									break;
								}
								 
								
				 
							}
							if(ns==-1) 
							{	 
								MyA[(stack*nSym+i)*size1+0]=-2;//subsumed now
								 
								_nSym--;
							}
							 
							break;
						}else{
							int p=MyA[(stack*nSym+i)*size1+5];
							for(;p<e_p_myb&& MyB[p+2]<s+1;p+=3)
							{	 
								Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],MyB[p],MyB[p+1]);
								if(x[g.index].in(g.val))
								{	
									GECODE_ME_CHECK(x[g.index].nq(home,g.val));
									ReA[re++]=g.index;
									ReA[re++]=g.val;
								}
							}
							MyA[(stack*nSym+i)*size1+5]=p;
						}
							
					}
					
					if(s==_psize) MyA[(stack*nSym+i)*size1+4]=s-1;
					
				}
		}

		while(re!=-1)
		{	if(_myc!=1)
				re=0;
			_myc=0; 
		 
			//find new watched literals
			for(int i=0;i<_csize;i++)
			{	 
				if(MyA[(stack*nSym+i)*size1+0]!=-2 && MyA[(stack*nSym+i)*size1+2]!=-1 &&	MyA[(stack*nSym+i)*size1+3]!=-1)//not subsumed yet and have constraints
				{	 
					if(!x[MyA[(stack*nSym+i)*size1+2]].in(MyA[(stack*nSym+i)*size1+3]))//the symmetry is broken 
					{	
						 
						MyA[(stack*nSym+i)*size1+0]=-2;
						_nSym--;
						 
						continue;
					}	
					if(x[MyA[(stack*nSym+i)*size1+2]].in(MyA[(stack*nSym+i)*size1+3])&&x[MyA[(stack*nSym+i)*size1+2]].assigned())//need to find a new first watched literal
					{	
						int p=MyA[(stack*nSym+i)*size1+5];
						for(;p<e_p_myb&& MyB[p+2]<MyA[(stack*nSym+i)*size1+4]+1;p+=3)
						{	 
							Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],MyB[p],MyB[p+1]);
							if(x[g.index].in(g.val))
							{	
								//std::cout<<"prune2 "<<MyA[(stack*nSym+i)*size1+1]<<":"<<MyB[p]<<","<<MyB[p+1]<<"\t"<<g.index<<","<<g.val<<"\n";
								GECODE_ME_CHECK(x[g.index].nq(home,g.val));
								ReA[re++]=g.index;
								ReA[re++]=g.val;
							}
						}
						MyA[(stack*nSym+i)*size1+5]=p;
						int s;
						//effect the RHSs
						 
						for(s=MyA[(stack*nSym+i)*size1+4]+1;s<_psize;s++)
						{
							 
							Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],PartA[2*s],PartA[1+2*s]);
							if(!x[g.index].in(g.val))//now the symmetry is broken
							{	
								
								MyA[(stack*nSym+i)*size1+0]=-2;
								 
								_nSym--;
								break;
							}
							if(x[g.index].in(g.val)&&!x[g.index].assigned())//now find the first pointer
							{
								MyA[(stack*nSym+i)*size1+2]=g.index;
								MyA[(stack*nSym+i)*size1+3]=g.val;
								MyA[(stack*nSym+i)*size1+4]=s; 
									int ns=0;
									int t;
									int p=MyA[(stack*nSym+i)*size1+5];
									for(t=s+1;t<_psize;t++) 
									{
							 
										 
										for(;p<e_p_myb&& MyB[p+2]<t;p+=3)
										{	Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],MyB[p],MyB[p+1]);
											if(x[g.index].in(g.val))
											{	
												ns=1;
												break;
											}
										}
										
										if(ns==1) break;
									 
										Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],PartA[t*2],PartA[1+t*2]);
										if(!x[g.index].in(g.val))//now the symmetry is broken
										{	
											ns=-1; 
											break;
										}
										 
										
						 
									}
									if(ns==-1) 
									{
										//std::cout<<"prune\t";
										MyA[(stack*nSym+i)*size1+0]=-2;//subsumed now
										//std::cout<<"sub "<<stack*nSym+i<<"\n";
										_nSym--;
									}
								 
									break;
							}
							else{
								int p=MyA[(stack*nSym+i)*size1+5];
								for(;p<e_p_myb&& MyB[p+2]<s+1;p+=3)
								{	Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],MyB[p],MyB[p+1]);
									if(x[g.index].in(g.val))
									{	
										//std::cout<<"prune2 "<<MyA[(stack*nSym+i)*size1+1]<<":"<<MyB[p]<<","<<MyB[p+1]<<"\t"<<g.index<<","<<g.val<<"\n";
										GECODE_ME_CHECK(x[g.index].nq(home,g.val));
										ReA[re++]=g.index;
										ReA[re++]=g.val;
									}
								}
								MyA[(stack*nSym+i)*size1+5]=p;
							}
						 
						}
						 
						//-----------------------------------------------------------------effect all nogoods whose LHSs are before s
						
						
						if(MyA[(stack*nSym+i)*size1+0]==-2) continue;
						 
						if(s==_psize)
						{
							//std::cout<<"set the first pointer\n";
							MyA[(stack*nSym+i)*size1+2]=-1;
							MyA[(stack*nSym+i)*size1+3]=-1;
							MyA[(stack*nSym+i)*size1+4]=s-1;
						}
						
					}
					
				
					
				}
				else if(MyA[(stack*nSym+i)*size1+0]!=-2)
				{	int p=MyA[(stack*nSym+i)*size1+5];
					for(;p<e_p_myb;p+=3)
						{	Index_class g=(*_symmetries)(MyA[(stack*nSym+i)*size1+1],MyB[p],MyB[p+1]);
							if(x[g.index].in(g.val))
							{	
								//std::cout<<"prune3 "<<MyA[(stack*nSym+i)*size1+1]<<":"<<MyB[p]<<","<<MyB[p+1]<<"\t"<<g.index<<","<<g.val<<"\n";
								GECODE_ME_CHECK(x[g.index].nq(home,g.val));
								 
								ReA[re++]=g.index;
								ReA[re++]=g.val;
							}
						}
					MyA[(stack*nSym+i)*size1+5]=p;
				}
				
				 
				
				
			}	
			//std::cout<<"\n"; 
			 
			for(int i=0;i<re;i+=2)
			{	 
				MyB[e_p_myb++]=ReA[i]; MyB[e_p_myb++]=ReA[i+1]; MyB[e_p_myb++]=_psize-1;
			} 
			if(re==0) re=-1;
			
		}
		
		if(_nSym==0) 
		{	 
			return home.ES_SUBSUMED(*this);
		}
		return ES_FIX;
		
				
			
			
	
	
	} 


	void gog(Space& home, ViewArray<Int::IntView> xv, int r) { 
		// constraint post function 
 
		if (GOG<Int::IntView>::post(home,xv,r) != ES_OK)
			home.fail(); 
	} 
}