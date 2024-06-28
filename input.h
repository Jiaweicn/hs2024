const double mmtacp=0.05;//momentum acceptance, 5%, GR
//const double mmtacp=0.3;//momentum acceptance, 30%, LAS

const string beam="a";
const double mb=mhe4;//rest mass of beam particle
const double tb=100.;//beam energy
string targets[]={"sn120","c12"};
string ejectiles[]={"h1","h2","h3","he3","he4","he6","he8","li6","li7","li8","be7","be8","be9","be10"};
//string ejectiles[]={"h1","h2","h3","he3","he4","he6","he8","li6","li7","li8","be7","be8","be9","be10","b8","b10","b11","b12"};
const short nejectiles=sizeof(ejectiles)/sizeof(ejectiles[0]);

void inputMass(const string &target,const string &ejectile,double &mp,double &mt,double &mr,double &qp){
	if(target==targets[0]){//rest mass of ejectile and residue
		mt=msn120;//mass of target 120Sn
		if(ejectile==ejectiles[0]){mp=mh1;mr=msb123;qp=1.;}
		else if(ejectile==ejectiles[1]){mp=mh2;mr=msb122;qp=1.;}
		else if(ejectile==ejectiles[2]){mp=mh3;mr=msb121;qp=1.;}
		else if(ejectile==ejectiles[3]){mp=mhe3;mr=msn121;qp=2.;}
		else if(ejectile==ejectiles[4]){mp=mhe4;mr=msn120;qp=2.;}
		else if(ejectile==ejectiles[5]){mp=mhe6;mr=msn118;qp=2.;}
		else if(ejectile==ejectiles[6]){mp=mhe8;mr=msn116;qp=2.;}
		else if(ejectile==ejectiles[7]){mp=mli6;mr=min118;qp=3.;}
		else if(ejectile==ejectiles[8]){mp=mli7;mr=min117;qp=3.;}
		else if(ejectile==ejectiles[9]){mp=mli8;mr=min116;qp=3.;}
		else if(ejectile==ejectiles[10]){mp=mbe7;mr=mcd117;qp=4.;}
		else if(ejectile==ejectiles[11]){mp=mbe8;mr=mcd116;qp=4.;}
		else if(ejectile==ejectiles[12]){mp=mbe9;mr=mcd115;qp=4.;}
		else if(ejectile==ejectiles[13]){mp=mbe10;mr=mcd114;qp=4.;}
	}
	else if(target==targets[1]){//12C target
		mt=mc12;//mass of target
		if(ejectile==ejectiles[0]){mp=mh1;mr=mn15;qp=1.;}
		else if(ejectile==ejectiles[1]){mp=mh2;mr=mn14;qp=1.;}
		else if(ejectile==ejectiles[2]){mp=mh3;mr=mn13;qp=1.;}
		else if(ejectile==ejectiles[3]){mp=mhe3;mr=mc13;qp=2.;}
		else if(ejectile==ejectiles[4]){mp=mhe4;mr=mc12;qp=2.;}
		else if(ejectile==ejectiles[5]){mp=mhe6;mr=mc10;qp=2.;}
		else if(ejectile==ejectiles[6]){mp=mhe8;mr=mc8;qp=2.;}
		else if(ejectile==ejectiles[7]){mp=mli6;mr=mb10;qp=3.;}
		else if(ejectile==ejectiles[8]){mp=mli7;mr=mb9;qp=3.;}
		else if(ejectile==ejectiles[9]){mp=mli8;mr=mb8;qp=3.;}
		else if(ejectile==ejectiles[10]){mp=mbe7;mr=mbe9;qp=4.;}
		else if(ejectile==ejectiles[11]){mp=mbe8;mr=mbe8;qp=4.;}
		else if(ejectile==ejectiles[12]){mp=mbe9;mr=mbe7;qp=4.;}
		else if(ejectile==ejectiles[13]){mp=mbe10;mr=mbe6;qp=4.;}
	}
}
void inputTarget(string &target){//input target correctly
	bool inputloop=true;
	string target0="c12";
	const short nts=sizeof(targets)/sizeof(targets[0]);
	while(inputloop){//target input loop
		cout<<"input target: ";for(short i=0;i<nts;i++)cout<<targets[i]<<" ";cout<<", default "<<target0<<"): ";
		getline(cin,target);
		if(target==""){
			target=target0;
			inputloop=false;
		}
		else{
			for(short i=0;i<nts;i++){
				if(target==targets[i]){
					//cout<<"target is: "<<BOLDGREEN<<target<<RESET<<endl;
					inputloop=false;
					break;
				}
				if(i==nts-1 && inputloop) cout<<"no such target, try again!\n";
			}
		}
		if(!inputloop) break;
	}
}
bool firstshow=true;
void inputEjectile(string &ejectile){//input ejectile correctly
	bool inputloop=true;
	string ejectile0="li6";
	const short nps=sizeof(ejectiles)/sizeof(ejectiles[0]);
	while(inputloop){//ejectile input loop
		if(firstshow){
			firstshow=false;
			cout<<"input ejectile(";for(short i=0;i<nps;i++)cout<<ejectiles[i]<<" ";cout<<"\n  default "<<ejectile0<<", q to quit: ";
		}
		else cout<<"input ejectile: ";
		getline(cin,ejectile);
		if(ejectile==""){
			ejectile=ejectile0;
			inputloop=false;
		}
		else if(ejectile=="q") exit(1);
		else{
			for(short i=0;i<nps;i++){
				if(ejectile==ejectiles[i]){
					//cout<<"ejectile is: "<<BOLDGREEN<<ejectile<<RESET<<endl;
					inputloop=false;
					break;
				}
				if(i==nps-1 && inputloop) cout<<"no such ejectile,try again!\n";
			}
		}
		if(!inputloop) break;
	}
}
double pc2tke(const double &pc0,double &qp,double &qp0,double &theta,double &mp,double &mt,double &mr){
	double pc=pc0*qp/qp0;
	double ep=sqrt(pc*pc+mp*mp);
	double tke=ep-mp;
	return tke;
}
double ex2pc(const double &ex,const double &eb,const double &pb,double &qp,double &qp0,double &theta,double &mp,double &mt,double &mr){
	double M=pow(ex+mr,2)-(mb*mb+mt*mt+mp*mp+2.*eb*mt);
	double pc=(-4.*M*pb*cos(theta)+sqrt(pow(4.*M*pb*cos(theta),2)-4.*(pow(2.*(mt+eb),2)-pow(2.*pb*cos(theta),2))*(pow(2.*(mt+eb)*mp,2)-M*M)))/(2.*(pow(2.*(mt+eb),2)-pow(2.*pb*cos(theta),2)));

	return pc*qp0/qp;
}
double pc2ex(const double &pc0,const double &eb,const double &pb,double &qp,double &qp0,double &theta,double &mp,double &mt,double &mr){
	double pc=pc0*qp/qp0;
	double ep=sqrt(pc*pc+mp*mp);
	double ex=sqrt(mb*mb+mt*mt+mp*mp+2*mt*(eb-ep)-2*eb*ep+2*pb*pc*cos(theta))-mr;
	return ex;
}
