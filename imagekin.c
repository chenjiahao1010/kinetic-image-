//double a,b,c,d; //a+b->c+d
/*class particle
{	
public:
		int M;//ion mass
		void init(double m ,double v,double e)
		double V,E;//velocity energy
}

void particle::init(double m,double v,double e)
{
		M=m
		V=v;
		E=e;
}*/
#include <iostream>
using namespace std;
void imagekin()
{double ekin;//incident energy

//particle a,b,c,d;
double am,bm,cm,dm;

ekin=96.0;//束流能量6*16MeV
double u=931.4940;
double Exeme=15;//emerge核激发的能量
double Eth_emerge=14.437;//emerge核的激发阈值
double Eth_reco=20;//反冲核的激发阈值,用不到时，default什么值都可以，这里写成20，不影响


am=16.00*u ;//入射粒子
bm=12.01*u;//靶
cm=16.00*u+Exeme;//emerge粒子（包含了激发能）
dm=12.01*u;//反冲粒子（也可以包含激发能）



//cout<<"请输入入射能量Ekin(MeV)： "<<endl;
//cin>>ekin;
//cout<<"请输入反应前后的（a+b->c+d）粒子质量(MeV) "<<endl;
//cin>> am >> bm >> cm >> dm;
//cin>>am;
//cin>>bm;
//cin>>cm;
//cin>>dm;
double av=sqrt(ekin*2/am);
double bv=0,be=0;
double cv=0,ce=0,dv=0,de=0;

		double		Q=-(cm+dm-am-bm);
	    double	Ec=bm/(am+bm)*ekin;
		double mccenter,cvcenter,dvcenter;//velocity of mass & c& d in CM
		cvcenter=sqrt(2*(Ec+Q)/(cm+cm*cm/dm));//emerge核质心系速度
		dvcenter=cm*cvcenter/dm;//反冲核质心系速度
		mccenter=am/(am+bm)*av;//质心在lab中的速度
//对速度归一来画图
double CEN=0.3;//质心画图的速度矢量长度
double CVCEN=cvcenter/mccenter*CEN;//emerge画图的矢量长度
double DVCEN=dvcenter/mccenter*CEN;

//画一级圆圈和箭头
TCanvas *c1=new TCanvas("c1");
c1->Range(0,0,1,1);
c1->SetWindowSize(800,800);
TPaveLabel * par= new TPaveLabel(0.1,0.9,0.9,0.99,"kinetic image of reaction by homepro");
par->SetFillColor(42);
par->Draw();
stringstream exstring;
exstring<<Exeme;
TString exstrings="Ex of emerge nuclei: "+exstring.str()+"MeV";
TPaveLabel * parex= new TPaveLabel(0.1,0.8,0.3,0.85,exstrings);

parex->SetFillColor(31);
parex->Draw();

double ax=0.2;//arrow ax of origin
double ay=0.5;//arrow ay of origin
TArrow *masscenter=new TArrow(ax,ay,ax+CEN,ay,0.05,">");
masscenter->Draw();
TEllipse *circleC=new TEllipse(ax+CEN,ay,CVCEN,CVCEN);
circleC->SetLineColor(2);
circleC->SetFillColorAlpha(2,0.5);
circleC->SetLineWidth(3);
circleC->Draw();

   auto tex = new TLatex(ax+CEN+CVCEN,ay,"16O");
 //tex->new TLatex()
   tex->SetTextColorAlpha(2, 0.6);
   tex->SetTextSize(0.04);	 
   tex->SetTextAngle(0);
   tex->Draw();


tex=new TLatex(ax+CEN,ay-DVCEN,"12C");
   tex->SetTextColorAlpha(4, 0.6);
   tex->SetTextSize(0.04);	 
   tex->SetTextAngle(0);
   tex->Draw();


double theta1=TMath::ASin(CVCEN/CEN);
double theta2=TMath::ASin(DVCEN/CEN);

TArrow *tangentline_c=new TArrow(ax,ay,ax+CEN-CVCEN*TMath::Sin(theta1),CVCEN*TMath::Cos(theta1)+ay,0.05,">");
TArrow *tangentline_d=new TArrow(ax,ay,ax+CEN-DVCEN*TMath::Sin(theta2),-DVCEN*TMath::Cos(theta2)+ay,0.05,">");


//std::string theta11=std::to_string(theta1);

double theta1an=theta1/(TMath::Pi())*180;
double theta2an=theta2/(TMath::Pi())*180;
//std::string theta22=std::to_string(theta2);
stringstream theta10,theta20;
theta10<<theta1an;
theta20<<theta2an;
TString theta11=theta10.str()+"#circ";
TString theta22=theta20.str()+"#circ";
tex=new TLatex(ax+CEN-CVCEN*TMath::Sin(theta1),ay+CVCEN*TMath::Cos(theta1),theta11);
   tex->SetTextColorAlpha(2, 0.6);
   tex->SetTextSize(0.03);	 
   tex->SetTextAngle(0);
   tex->Draw();


tex=new TLatex(ax+CEN-DVCEN*TMath::Sin(theta2),ay-DVCEN*TMath::Cos(theta2),theta22);
   tex->SetTextColorAlpha(4, 0.6);
   tex->SetTextSize(0.03);	 
   tex->SetTextAngle(0);
   tex->Draw();

   

tangentline_c->Draw();
tangentline_d->Draw();

TEllipse *circleD=new TEllipse(ax+CEN,ay,DVCEN,DVCEN);
circleD->SetFillColorAlpha(4,0.5);
circleD->SetLineColorAlpha(4,0.5);
circleD->SetLineWidth(3);
circleD->Draw();

//画破裂核的圆圈（二级）有特异性，需要看怎么破裂
   double emeangdeg=15.0;
   double emeang=emeangdeg/180*TMath::Pi();//碎片角度
   //针对16O->4alpha
   double malpha=4.003*u;
   double ekalpha=(Exeme-Eth_emerge)/4;
   double valpha=sqrt(2*ekalpha/malpha);
//这里4α认为各项同性，所以画一个圈圈就行了
    TEllipse *circleC1=new TEllipse(ax+CEN+CVCEN*TMath::Cos(emeang),ay+CVCEN*TMath::Sin(emeang),valpha,valpha);
    circleC1->SetFillColorAlpha(3,0.5);
    circleC1->SetLineColorAlpha(3,0.5);
    circleC1->SetLineWidth(3);
    circleC1->Draw();

//计算质心origin和alpha的外边的切线的箭头
double alphax=ax+CEN+CVCEN*TMath::Cos(emeang);
double alphay=ay+CVCEN*TMath::Sin(emeang);
double stringalpha=sqrt((alphax-ax)*(alphax-ax)+(alphay-ay)*(alphay-ay));
double phi1=TMath::ASin(valpha/stringalpha);//原点与alpha圆心连线与外切线之间的夹角
double phi2=TMath::ATan((alphay-ay)/(alphax-ax));//原点与alpha圆心连线与水平线之间的夹角
double phi=phi1+phi2;
double tangentstring=stringalpha*TMath::Cos(phi1);//原点到alpha外切点的距离
double tangentx=ax+tangentstring*TMath::Cos(phi);
double tangenty=ay+tangentstring*TMath::Sin(phi);

TArrow *tangentline_c1=new TArrow(ax,ay,tangentx,tangenty,0.05,">");
tangentline_c1->SetLineColor(3);
tangentline_c1->Draw();

   tex=new TLatex(tangentx,tangenty,"#alpha");
   tex->SetTextColorAlpha(3, 0.6);
   tex->SetTextSize(0.04);	 
   tex->SetTextAngle(0);
   tex->Draw();

stringstream thetaalpha;
double tandeg=phi/TMath::Pi()*180.0;
thetaalpha<<tandeg;
stringstream emeangdegs,apeangs;
emeangdegs<<emeangdeg;
double phi1deg=phi1/TMath::Pi()*180.0;
double apeang=2.0*phi1deg;
apeangs<<apeang;
TString thetaalphas="#splitline{"+thetaalpha.str()+"#circ (arbitrary)}{"+"aperture angle :  "+apeangs.str()+"#circ  In CM: "+emeangdegs.str()+"#circ}";
tex=new TLatex(tangentx-0.2,tangenty-0.05,thetaalphas);
   tex->SetTextColorAlpha(3, 0.6);
   tex->SetTextSize(0.02);	 
   tex->SetTextAngle(0);
   tex->Draw();


}



