//return the excitation energy of the residuals for the ejectiles with same magnetic rigidity
//need to set the excitation energy and outgoing particles
//reaction: a+b->c+d

#include "/home/cai/prjs/include/cpp/outColor.h"
#include "mass.h"
#include "input.h"

void evtRteEst(double ex=2.,double theta=2.){//from the excitation energy of residue to the momentum of ejectile
	string target="c12",ejectile="li6";
	double mp=0., mr=0., mt, qp0, qp;//qp, charge of ejectile
	cout<<"incident beam is: "<<BOLDGREEN<<beam<<" ("<<tb<<" MeV)"<<RESET<<endl;

	cout<<"  ex of residual"<<": "<<BOLDGREEN<<ex<<RESET<<"MeV\n";
	cout<<"  outgoing angle(lab)"<<ejectile<<": "<<BOLDGREEN<<theta<<RESET<<"deg\n";
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
	cout<<" momemntum: "<<GREEN<<pc<<YELLOW<<" ("<<pclow<<", "<<pchigh<<RESET<<")MeV/c\n\n";
	for(short i=0;i<nejectiles;i++){
		ejectile=ejectiles[i];
		cout<<"ejectile: "<<ejectile<<endl;
		inputMass(target,ejectile,mp,mt,mr,qp);
		double excenter=pc2ex(pc,eb,pb,qp,qp0,theta,mp,mt,mr);
		double exlow=pc2ex(pchigh,eb,pb,qp,qp0,theta,mp,mt,mr);
		double exhigh=pc2ex(pclow,eb,pb,qp,qp0,theta,mp,mt,mr);
		double tkecenter=pc2tke(pc,qp,qp0,theta,mp,mt,mr);
		double tkelow=pc2tke(pclow,qp,qp0,theta,mp,mt,mr);
		double tkehigh=pc2tke(pchigh,qp,qp0,theta,mp,mt,mr);
		if(exhigh<0.){
			cout<<" max ex < 0 MeV\n";
			continue;
		}
		if(exlow<0.){
			cout<<" min ex "<<exlow<<" < 0 MeV, changed to be 0\n";
			double p0MeV=ex2pc(0.,eb,pb,qp,qp0,theta,mp,mt,mr);
			exlow=pc2ex(p0MeV,eb,pb,qp,qp0,theta,mp,mt,mr);
			tkehigh=pc2tke(p0MeV,qp,qp0,theta,mp,mt,mr);
		}
		cout<<setprecision(3)<<fixed;
		cout<<" ex: "<<GREEN<<excenter<<YELLOW<<" ("<<exlow<<", "<<exhigh<<RESET<<")MeV,";
		cout<<" tke: "<<GREEN<<tkecenter<<YELLOW<<" ("<<tkelow<<", "<<tkehigh<<RESET<<")MeV\n";
	}
}
