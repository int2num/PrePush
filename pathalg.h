//
// Created by root on 17-5-9.
//
#ifndef CSPALLOC_PATHALG_H
#include"limits.h"
#define CSPALLOC_PATHALG_H
#define INFCOST INT_MAX/2
#include<bits/stdc++.h>
#include"edge.h"
#define ML 50
#define BS 5
#define WD 15
#define inf INT_MAX/2
using namespace std;
class pairless {
    public:
        bool operator()(pair<int,int>a,pair<int,int>b){
            return a.second>b.second;
        }
};
class algbase {
    protected:
        vector<int> getrout(int &s, int &t, vector<edge> &edges, vector<int> &pre) {
            vector<int> rout;
            int pp = pre[t];
            while (pp >= 0) {
                rout.push_back(pp);
                pp = pre[edges[pp].s];
            }
            reverse(rout.begin(), rout.end());
            return rout;
        }
    public:
        algbase(){};
        virtual bool cutcake(int)=0;
        virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo)=0;
	 	virtual vector<int> routalg(int s,int t,int bw)=0;
};

class dijkstor:public algbase{
    public:
		vector<vector<int>>neibour;
		vector<int>ancestor;
		int edgesize,nodenum,pesize,pnodesize,maxbw;
		vector<vector<vector<int>>>mask;
		vector<edge>edges;
		vector<int>dist;
		vector<int>pre;
		vector<int>leveln;
		vector<int>exn2n;
		vector<vector<int>>rela;
		vector<int>order;
        dijkstor(){};
        vector<int> topsort()
        {
        	cout<<" in top sort "<<endl;
        	queue<int>zero;
        	vector<int>orz(nodenum,-1);
        	order=orz;
        	for(int i=0;i<pnodesize;i++)
        		zero.push(i);
        	int biao=0;
			while(!zero.empty())
			{
				int node=zero.front();
				zero.pop();
				order[node]=biao++;
				for(int i=0;i<neibour[node].size();i++)
				{
					if((--ancestor[edges[neibour[node][i]].t])==0)
							zero.push(edges[neibour[node][i]].t);
				}
			}
			return order;
        }
        virtual bool cutcake(int index){
        	cout<<"cut "<<index<<endl;
        	if(maxbw-(index+1)*10>=0)
        			maxbw-=(index+1)*10;
			else
				{
					cout<<"failure"<<endl;
					return false;
				}
        	vector<int>tmp;
        	for(int i=0;i<edgesize;i++)
        		tmp.push_back(i);
        	mask[index].push_back(tmp);
        	leveln[index]++;
        	return true;
        }
        virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf){
        	maxbw=500;
        	rela=relate;
        	edgesize=extenedges.size(),nodenum=ginf.enodesize,pesize=ginf.pesize,pnodesize=ginf.pnodesize;
        	exn2n=ginf.exn2n;
        	edges=extenedges;
        	vector<vector<int>>nd(nodenum,vector<int>());
        	neibour=nd;
        	vector<int>ad(nodenum,0);
        	ancestor=ad;
        	for(int i=0;i<edgesize;i++)
        		{
        			neibour[extenedges[i].s].push_back(extenedges[i].t);
        			ancestor[extenedges[i].t]++;
        		}
    		vector<int>bl(BS,0);
    		leveln=bl;
    		vector<vector<vector<int>>>nm(BS,vector<vector<int>>());
    		mask=nm;
    		vector<int>dd(nodenum*ML,inf);
    		dist=dd;
    		vector<int>pp(nodenum*ML,-1);
    		pre=pp;
    		topsort();
        }
        virtual vector<int> routalg(int s,int t,int bw){
        	int index=bw/10-1;
        	if(leveln[index]==0)if(!cutcake(index))return vector<int>();
        	int tnode=-1;
        	while(true)
        	{
				vector<int>visited(nodenum,0);
				vector<int>pre(nodenum,-1);
				int vflag=1;
				queue<int>que;
				que.push(s);
	        	visited[s]=1;
				while(!que.empty()&&vflag)
				{
					int node=que.front();
					que.pop();
					for(int i=0;i<neibour[node].size();i++)
					{
						int to=neibour[node][i];
						if(visited[to]==0){pre[to]=node;que.push(to);visited[to]=1;}
						else{continue;}
						if(exn2n[to]==t){tnode=to;vflag=0;break;}
					}

				}
				vector<int>rout;
				int pp=pre[tnode];
				cout<<"tnode is "<<tnode<<endl;
				while(pp>=0)
				{
					rout.push_back(pp);
					pp=pre[pp];
				}
				cout<<endl;
				if(rout.size()>0)
				{
					for(int i=0;i<rout.size();i++)
						cout<<rout[i]<<" ";
					cout<<endl;
					return rout;
				}
				else
				{
					if(!cutcake(index))return vector<int>();
				}
        	}
	 	}
};
class parallelor:public algbase
{
	private:
		edge *dev_edges,*aedges;
		int*dev_m,*m,*dev_pre,*pre,*pred,*dev_pred,*dev_d,*d,*dev_mask,*mask,*dev_leveln,*leveln;
		int*dev_rela,*rela;
		int *dev_chan,*chan;
		int presize,dsize,masksize,levelnsize;
		int edgesize,nodenum,pesize,pnodesize;
		int neisize,duansize;
		int *choosel,*dev_choosel;
		int *rout,*dev_rout;
		int *routn,*dev_routn,*order,*dev_order;
		vector<int>hleveln,ancestor;
		int *ordernode,*dev_ordernode;
		int maxbw;
		int *dev_qian,*qian,*dev_qbeg,*qbeg,*dev_qsize,*qsize;
		epair *dev_nei,*nei;
		int *dev_duan,*duan;
		int *dev_beg,*beg;
		vector<vector<int>>neibour;
		void allocate(int maxn,int maxedges);
		void copydata(int s,vector<edge>&edges,int nodenum);
		void dellocate();
	public:
	 	 parallelor();
	 	 void topsort();
	 	 virtual bool cutcake(int index);
	 	 virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo);
	 	 virtual vector<int> routalg(int s,int t,int bw);
	 	 virtual ~parallelor(){dellocate();}
};
class paircmp{
    operator ()(pair<int,int>a,pair<int,int>b)
    {
        return a.second<b.second;
    }
};
class PreflowPush:public algbase{
public:
	vector<int>order;
    vector<edges>alledges;
	virtual bool cutcake(int){};
	virtual vector<int> routalg(int s,int t,int bw){};
	virtual void init(vector<edge> &extenedges, vector<vector<int>> &relate, ginfo gg) {
        dijkstor plo;
        alledges=extenedges;
        int n=gg.pnodesize;
        plo.init(extenedges,relate,gg);
        order=plo.topsort();
        cout<<gg.pnodesize<<endl;
        vector<pair<int,int>>st;
        vector<vector<int>>neibour;
        for(int i=0;i<1000;i++)
        {
            int s=rand()%n;
            int t=s;
            while(t==s)t=rand()%n;
            t+WD*n;
            st.push_back(make_pair(s,t));
        }
        vector<pair<int,int>>vv;
        for(int i=0;i<order.size();i++)
            vv.push_back(make_pair(i,order[i]));
        sort(vv.begin(),vv.end());
        start=gg.pnodesize*(WD+1);
        end=start+1;
        vector<vector<int>>(end+1,vector<int>);
        for(int i=0;i<st.size();i++)
        {
            alledges.push_back(edge(start,st.first));
            alledges.push_back(edge(st.second,end));
        }
        for(int i=0;i<alledges.size();i++)
            neibour[alledges[i].s].push_back(alledges[i].t);

    };
};
#endif //CSPALLOC_PATHALG_H
