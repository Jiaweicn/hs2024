//producing energy loss correlation for E599 experiment
#include "/home/cai/prjs/include/cpp/outColor.h"
#include "mass.h"
#include "input.h"
#include <boost/filesystem.hpp>

void setBrsAdd(TTree *ti,double *ene){
	ti->Branch("eTarget", &ene[0],"eTarget/D");
	ti->Branch("eGR", &ene[1],"eGR/D");
	ti->Branch("eAir0", &ene[2],"eAir0/D");
	ti->Branch("eVDC1Window0", &ene[3],"eVDC1Window0/D");
	ti->Branch("eVDC1Cathode", &ene[4],"eVDC1Cathode/D");
	ti->Branch("eVDC1", &ene[5],"eVDC1/D");
	ti->Branch("eVDC1Window1", &ene[6],"eVDC1Window1/D");
	ti->Branch("eHeBoxWindow0", &ene[7],"eHeBoxWindow0/D");
	ti->Branch("eHeBox", &ene[8],"eHeBox/D");
	ti->Branch("eHeBoxWindow1", &ene[9],"eHeBoxWindow1/D");
	ti->Branch("eVDC2Window0", &ene[10],"eVDC2Window0/D");
	ti->Branch("eVDC2Cathode", &ene[11],"eVDC2Cathode/D");
	ti->Branch("eVDC2", &ene[12],"eVDC2/D");
	ti->Branch("eVDC2Window1", &ene[13],"eVDC2Window1/D");
	ti->Branch("eAir1", &ene[14],"eAir1/D");
	ti->Branch("eBoxShading0", &ene[15],"eBoxShading0/D");
	ti->Branch("ePs1Shading0", &ene[16],"ePs1Shading0/D");
	ti->Branch("ePs1", &ene[17],"ePs1/D");
	ti->Branch("ePs1Shading1", &ene[18],"ePs1Shading1/D");
	ti->Branch("eAir2", &ene[19],"eAir2/D");
	ti->Branch("ePs2Shading0", &ene[20],"ePs2Shading0/D");
	ti->Branch("ePs2", &ene[21],"ePs2/D");
	ti->Branch("ePs2Shading1", &ene[22],"ePs2Shading1/D");
	ti->Branch("eAir3", &ene[23],"eAir3/D");
}

void readSRIM(char file[50], TGraph *grdedx){
	ifstream fin(file);
	//cout<<" file: "<<file<<endl;
	string line, unit;
	double e, dedx0,dedx1;//keV, keV/um
	while(getline(fin, line)){
		stringstream ss;
		ss.str(line);
		ss>>e>>unit>>dedx0>>dedx1;//reading SRIM file
		if(grdedx->GetN()==0 && unit!="keV") continue;
		if(unit=="keV") e=e/1000./*MeV*/;
		else if(unit=="MeV") NULL;
		else break;

		grdedx->SetPoint(grdedx->GetN(),e, (dedx0+dedx1)/1000./*MeV/um*/);
	}
	fin.close();
}
double getdE(const double &e0,TGraph *grdedx,const double &thickness){
	if(e0<0.) return 0.;
	double e=e0, x=0.;//um
	const double dx=0.1;//um
	while(x<thickness){
		x=x+dx;
		e=e-dx*grdedx->Eval(e);
		if(e<0.) break;
	}
	return e0-e;
}

const string layers[]={"target","GRWindow","Air","VDC1Window0","VDC1Cathode","VDC1","VDC1Window1","HeBoxWindow0","HeBox","HeBoxWindow1","VDC2Window0","VDC2Cathode","VDC2","VDC2Window1","Air","BoxShading0","ps1Shading0","ps1","ps1Shading1","Air","ps2Shading0","ps2","ps2Shading1","Air"};//names of layers
const short nlayers=sizeof(layers)/sizeof(layers[0]);

void fillEloss(const string &ejectile,const double &tkelow,const double &tkehigh,double *ene,TTree *ti,TGraph *gr1,TGraph *gr2,TGraph *gr3){
	//aramid: aramid
	const string mats[]={"c12","aramid","air","aramid","aramid","c1h4he5","aramid","aramid","he","aramid","aramid","aramid","c1h4he5","aramid","air","mylar","mylar","c9h10","mylar","air","mylar","c9h10","mylar","air"};//names of materials
	const double matThick[]={2,50,60000,12.5,18,60000,12.5,12,190000,12,12.5,18,60000,12.5,55000,24,10,1000,10,70000,10,3000,10,70000};//thicknesses of materials in um
	const double kinstep=.1;//MeV
	const double angle=45.*TMath::DegToRad();//angle of layers
	double kin0=tkelow;
	while(kin0<tkehigh){
		double kin=kin0;
		for(short i=0;i<nlayers;i++){
			char filename[50];
			sprintf(filename,"../include/dedxSRIM/%s-%s.txt",ejectile.c_str(),mats[i].c_str());
			if(!boost::filesystem::exists(filename)){
				cout<<filename<<" does not exit!\n";
				exit(1);
			}
			auto *grdedx=new TGraph();
			readSRIM(filename,grdedx);
			if(i==0) ene[i]=getdE(kin,grdedx,matThick[i]);
			else if(i>0) ene[i]=getdE(kin,grdedx,matThick[i]/cos(angle));
			kin=kin-ene[i];
		}
		if(ene[12]>0.01) gr1->SetPoint(gr1->GetN(),ene[12],ene[5]);
		if(ene[17]>0.01) gr2->SetPoint(gr2->GetN(),ene[17],ene[12]);
		if(ene[21]>0.01) gr3->SetPoint(gr3->GetN(),ene[21],ene[17]);
		ti->Fill();
		kin0=kin0+kinstep;
	}
}

void rcnpPID(double ex=4.3,double theta=2.5){
	string target,ejectile;
	double mp=0., mr=0., mt, qp0, qp;//beam energy, qp, charge of ejectile
	cout<<"  ex energy of residual"<<": "<<GREEN<<ex<<RESET<<"MeV\n";
	cout<<"  outgoing angle(lab) of "<<ejectile<<": "<<GREEN<<theta<<RESET<<"deg\n";
	theta=theta*TMath::DegToRad();
	inputTarget(target);
	inputEjectile(ejectile);
	inputMass(target,ejectile,mp,mt,mr,qp0);
	const double eb=mb+tb;
	const double pb=sqrt(eb*eb-mb*mb);
	double cost=(mt*mt+mb*mb+mp*mp+2.*mt*eb-TMath::Sq(ex+mr))/2.;
	double ca=TMath::Sq(mt+eb)-TMath::Sq(pb*cos(theta));
	double cb=-2.*cost*pb*cos(theta);
	double cc=TMath::Sq((mt+eb)*mp)-cost*cost;
	double pc=(-cb+sqrt(cb*cb-4*ca*cc))/(2*ca);
	double pclow=(1.-mmtacp/2)*pc;
	double pchigh=(1.+mmtacp/2)*pc;
	cout<<" BRho setting: "<<GREEN<<3.3356e-3*pc/qp0<<RESET<<" Tm\n";
	cout<<" momemntum: "<<GREEN<<pc<<" ("<<YELLOW<<pclow<<", "<<pchigh<<RESET<<")MeV/c\n\n";
	char rootout[30];
	sprintf(rootout,"%s_%s.root",ejectile.c_str(),target.c_str());
	auto *fo = TFile::Open(rootout,"recreate");
	TTree *ts[nejectiles];
	TGraph *gr1s[nejectiles],*gr2s[nejectiles],*gr3s[nejectiles];
	auto *mgr1=new TMultiGraph();
	auto *mgr2=new TMultiGraph();
	auto *mgr3=new TMultiGraph();
	mgr1->SetName("mgr1");
	mgr2->SetName("mgr2");
	mgr3->SetName("mgr3");
	mgr1->SetTitle("PID corr.;VDC2/MeV;VDC1/MeV");
	mgr2->SetTitle("PID corr.;PS1/MeV;VDC2/MeV");
	mgr3->SetTitle("PID corr.;PS2/MeV;PS1/MeV");
	double ene[nlayers]={0.};
	auto *lg1=new TLegend();
	lg1->SetName("lg1");
	lg1->SetFillStyle(0);
	lg1->SetBorderSize(0);
	auto *lg2=new TLegend();
	lg2->SetName("lg2");
	lg2->SetFillStyle(0);
	lg2->SetBorderSize(0);
	auto *lg3=new TLegend();
	lg3->SetName("lg3");
	lg3->SetFillStyle(0);
	lg3->SetBorderSize(0);
	short colorNO=1;
	for(short i=0;i<nejectiles;i++){
		gr1s[i]=new TGraph();
		gr2s[i]=new TGraph();
		gr3s[i]=new TGraph();
		gr1s[i]->SetMarkerStyle(8);
		gr2s[i]->SetMarkerStyle(8);
		gr3s[i]->SetMarkerStyle(8);
		gr1s[i]->SetMarkerSize(2);
		gr2s[i]->SetMarkerSize(2);
		gr3s[i]->SetMarkerSize(2);
		gr1s[i]->SetMarkerColor(colorNO);
		gr2s[i]->SetMarkerColor(colorNO);
		gr3s[i]->SetMarkerColor(colorNO);
		gr1s[i]->SetLineColor(colorNO);
		gr2s[i]->SetLineColor(colorNO);
		gr3s[i]->SetLineColor(colorNO);
		if(colorNO>=10){
			gr1s[i]->SetMarkerColor(colorNO+1);
			gr2s[i]->SetMarkerColor(colorNO+1);
			gr3s[i]->SetMarkerColor(colorNO+1);
			gr1s[i]->SetLineColor(colorNO+1);
			gr2s[i]->SetLineColor(colorNO+1);
			gr3s[i]->SetLineColor(colorNO+1);
		}
		colorNO++;
		gr1s[i]->SetLineWidth(2);
		gr2s[i]->SetLineWidth(2);
		gr3s[i]->SetLineWidth(2);
		inputMass(target,ejectiles[i],mp,mt,mr,qp);
		ts[i]=new TTree();
		ts[i]->SetName(ejectiles[i].c_str());
		setBrsAdd(ts[i],ene);
		double excenter=pc2ex(pc,eb,pb,qp,qp0,theta,mp,mt,mr);
		double exlow=pc2ex(pchigh,eb,pb,qp,qp0,theta,mp,mt,mr);
		double exhigh=pc2ex(pclow,eb,pb,qp,qp0,theta,mp,mt,mr);
		double tkecenter=pc2tke(pc,qp,qp0,theta,mp,mt,mr);
		double tkelow=pc2tke(pclow,qp,qp0,theta,mp,mt,mr);
		double tkehigh=pc2tke(pchigh,qp,qp0,theta,mp,mt,mr);
		cout<<"Calculating "<<GREEN<<ejectiles[i]<<RESET<<", ex "<<excenter<<" tke "<<tkecenter<<endl;
		if(exhigh<0.){
			cout<<" max excitation < 0MeV\n";
			continue;
		}
		if(exlow<0.){
			cout<<" min excitation < 0MeV, Corrected to be 0\n";
			double p0MeV=ex2pc(0.0,eb,pb,qp,qp0,theta,mp,mt,mr);
			exlow=pc2ex(p0MeV,eb,pb,qp,qp0,theta,mp,mt,mr);
			tkehigh=pc2tke(p0MeV,qp,qp0,theta,mp,mt,mr);
		}
		cout<<" filling eloss tree for "<<GREEN<<ejectiles[i]<<RESET<<endl;
		fillEloss(ejectiles[i],tkelow,tkehigh,ene,ts[i],gr1s[i],gr2s[i],gr3s[i]);
		if(gr1s[i]->GetN()>0){
			mgr1->Add(gr1s[i]);
			lg1->AddEntry(gr1s[i],ejectiles[i].c_str());
		}
		if(gr2s[i]->GetN()>0){
			mgr2->Add(gr2s[i]);
			lg2->AddEntry(gr2s[i],ejectiles[i].c_str());
		}
		if(gr3s[i]->GetN()>0){
			mgr3->Add(gr3s[i]);
			lg3->AddEntry(gr3s[i],ejectiles[i].c_str());
		}
	}
	fo->cd();
	for(short i=0;i<nejectiles;i++) ts[i]->Write();
	mgr1->Write();
	lg1->Write();
	mgr2->Write();
	lg2->Write();
	mgr3->Write();
	lg3->Write();
	fo->Close();
	cout<<"Saved in: "<<GREEN<<rootout<<RESET<<endl;
}
