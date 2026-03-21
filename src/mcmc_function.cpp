#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
using namespace Rcpp ;
using namespace RcppParallel;


//########################################################################################################################################
// function to generate normal random variable

double rnorm(double a, double b) { // a: mean, b: s.d.
	double c=a+b*sum(rnorm(1));
    return c;
}

//########################################################################################################################################
// function to generate normal random variable

double gen_gamma(const double a, const double b) { // a: mean, b: s.d.
    return R::rgamma(a,b);
}

// directly call gamma(double) for gamma function


//########################################################################################################################################
// function to generate normal random variable

double d_gamma(const double x, const double a, const double b) { // a: mean, b: s.d.
    return R::dgamma(x,a,b,0);
}


//########################################################################################################################################
// function to generate binomial random number

int gen_binom(const double p){
double cut=R::runif(0,1);
int out=0;
if (cut<p){
out=1;
}
return out;
}


//########################################################################################################################################
// function to general multiviariable normal given sigma

NumericVector rmnorm(arma::mat sigma) {
int ncols=sigma.n_cols;
arma::rowvec c=arma::randn(1,ncols);
arma::rowvec a=c*arma::chol(sigma);   
NumericVector b=NumericVector(a.begin(),a.end());   
return b;
}

//########################################################################################################################################
//function to compute the prior likelihood 
// [[Rcpp::export]]
double prior_loglik(NumericVector para, int n_seasons){
// check if the para are within their possible range
double out=1;
int b1;
int hai_start = 42 + 3 * n_seasons;
int para_end = hai_start + 2 * n_seasons - 1;
out*=sum(dunif(NumericVector::create(para(0)),0.00000001,0.9999999));
for (b1=hai_start-1;b1>=1;--b1){
out*=sum(dunif(NumericVector::create(para(b1)),0.01,100.0));
}
for (b1=para_end;b1>=hai_start;--b1){
out*=sum(dnorm(NumericVector::create(para(b1)),0.0,1.0));
}
// if the prior is outside the parameter space
if (out==0){
out=-9999999;
}
else{
out=log(out);	
}
return out;
}


//########################################################################################################################################
//function to do simulation
// [[Rcpp::export]]
List sim_data(NumericMatrix data1,
NumericMatrix ILI,
NumericVector para,
NumericVector para2,
int hai_start){
// clone the data first
NumericMatrix data11(clone(data1)); // record the true titer
NumericMatrix data111(clone(data1)); // record the output with measurement error
// the augmented data to store information
NumericMatrix data21(data11.nrow(),10);
NumericMatrix record(data11.nrow(),30);
NumericMatrix record2(1,1);
// 1-3 hhID, member, age
// 4, infection or not 
// 5, infection time
// 6, boosting
// 7 waning rate 

int b1;
int b2;
int b3;
//int b4;

for (b1=data11.nrow()-1;b1>=0;--b1){
// first impute the basic information here
data21(b1,0)=data11(b1,0);	
data21(b1,1)=data11(b1,1);
data21(b1,2)=data11(b1,2);
// boosting, cross-boosting and waning parameter
data21(b1,5)=R::rgamma(para[18+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);
data21(b1,6)=R::rgamma(para[19+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);

// here generate the baseline titer level
// first select 0-9
double titerlevel1=R::runif(0,1);
for (b3=9;b3>=0;--b3){
titerlevel1-=para2[b3+10*(data21(b1,2)>=1)+20*data111(b1,11)];
if (titerlevel1<0){
data11(b1,8)=b3;
// break out from the loop 
break;
}
}

data11(b1,8)+=R::runif(0,1);


// here start to generate the infection status
// 499 is the change point
for (b2=data11(b1,3)+1;b2<=data11(b1,4);++b2){
if (data21(b1,3)==0){ // at risk 
double hazard1=ILI(b2-1,data111(b1,12))*para[42+data21(b1,2)+3*data111(b1,11)]*(data21(b1,3)==0); //*pow(3,b2>509); //+3*(b2>509)
int inf=gen_binom(1-exp(-hazard1*exp(data11(b1,8)*para[hai_start+2*data11(b1,11)])));
if (inf==1){
data21(b1,3)=1;
data21(b1,4)=b2;
}
}
}

// based on the infection status, generate the antibody respose
for (b2=0;b2<=1;++b2){
// first check if imputing is necessary
if (data1(b1,6+b2)!=-1){
// first check if there is infection with possible boosting in that interval
int cond1=(data1(b1,5+b2)<data21(b1,4))&&(data21(b1,4)<=data1(b1,6+b2));
if (cond1==0){
// non-infection	
data11(b1,9+b2)=data11(b1,8+b2)*exp(-(data1(b1,6+b2)-data1(b1,5+b2))*data21(b1,6)/365.0);
}
else {
// infection
double beforeinf=data11(b1,8+b2)*exp(-(data21(b1,4)-data1(b1,5+b2))*data21(b1,6)/365.0); 
data11(b1,9+b2)=(beforeinf+data21(b1,5))*exp(-(data1(b1,6+b2)-data21(b1,4)-1)*data21(b1,6)/365.0);
}
data11(b1,9+b2)=data11(b1,9+b2)*(data11(b1,9+b2)<17)*(data11(b1,9+b2)>0)+17*(data11(b1,9+b2)>17);
}
}

// here to add the measurement error 
for (b2=0;b2<=2;++b2){
if (data1(b1,5+b2)!=-1){	

double probtemp[20];
double probtemptotal=0;

for (b3=9;b3>=0;--b3){
// the random error
probtemp[b3]=para[0];
probtemp[10+b3]=0.000001;
}

double measured=data11(b1,8+b2); //+para[0]*data11(b1,12);

if (measured<1){
// two-fold error at 0 
probtemp[0]+=1.0-measured/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[1]+=measured/(2.0*exp(para[1+1*(data11(b1,12))]));
}
else{
// two-fold error at other step	
int index=floor(measured);
probtemp[index]+=1.0-(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]))-(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index-1]+=(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index+1]+=(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
}

for (b3=19;b3>=0;--b3){	
probtemptotal+=probtemp[b3];
//if (b2==0){
//record(b1,b3)=probtemp[b3];
//}
}

double titerlevel=R::runif(0,1);
titerlevel*=probtemptotal;
//record2(b1,b2+b4*8)=titerlevel;
//record2(b1,b2+b4*8+4)=probtemptotal;
for (b3=19;b3>=0;--b3){
titerlevel-=probtemp[b3];
if (titerlevel<0){
data111(b1,8+b2)=b3;
// break out from the loop 
break;
}
}

}
}

// missing baseline
if (data1(b1,8)==-1){
data111(b1,8)=-1;
}

}





return List::create(_[""]=data11,
_[""]=data111,	
_[""]=data21,
_[""]=record,
_[""]=record2);
} 




//########################################################################################################################################
//########################################################################################################################################
// the main body of the parallel function
struct LogLik:public Worker{
// source vector
RMatrix<double> out1;
RMatrix<double> out2;
RMatrix<double> out3;
RMatrix<double> data11;
RMatrix<double> data111;
RMatrix<double> data21;
RMatrix<double> ILI;
RVector<double> para;
RVector<double> para2;
int level1;
int level2;
int level3;
int season;
int hai_start;
// destination vector
// initialize with source and destination
LogLik(NumericMatrix out1,
NumericMatrix out2,
NumericMatrix out3,
NumericMatrix data11,
NumericMatrix data111,
NumericMatrix data21,
NumericMatrix ILI,
NumericVector para,
NumericVector para2,
int level1,
int level2,
int level3,
int season,
int hai_start)
:out1(out1),out2(out2),out3(out3),data11(data11),data111(data111),data21(data21),ILI(ILI),para(para),para2(para2),level1(level1),level2(level2),level3(level3),season(season),hai_start(hai_start){}
void operator()(std::size_t begin, std::size_t end) {

// section to write the parallel version
// functor (pass input and output matrixes)
for (unsigned int b1=begin;b1<end;++b1){
int b2;
int b3;
//int b4;
//int b5;

if (level1){
// first level to record the measurment error likelihood
for (b2=2;b2>=0;--b2){
if (data111(b1,8+b2)!=-1){

double probtemp[20];
double probtemptotal=0;


for (b3=9;b3>=0;--b3){
// the random error
probtemp[b3]=para[0];
probtemp[10+b3]=0.000001;
}

double measured=data11(b1,8+b2); //+para[0]*data11(b1,12);

if (measured<1){
// two-fold error at 0 
probtemp[0]+=1.0-measured/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[1]+=measured/(2.0*exp(para[1+1*(data11(b1,12))]));
}
else{
// two-fold error at other step	
int index=floor(measured);
probtemp[index]+=1.0-(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]))-(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index-1]+=(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index+1]+=(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
}


for (b3=19;b3>=0;--b3){	
probtemptotal+=probtemp[b3];
//if (b2==0){
//record(b1,b3)=probtemp[b3];
//}
}


out1(b1,b2)=log(probtemp[int(data111(b1,8+b2))])-log(probtemptotal);

}
}
}





if (level2){
if (data21(b1,3)==1){
out2(b1,0)=R::dgamma(data21(b1,5),para[18+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0,1);
out2(b1,1)=R::dgamma(data21(b1,6),para[19+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0,1);
}
if (data111(b1,8)!=-1){
out2(b1,2)=log(para2[10*(data21(b1,2)>=1)+floor(data11(b1,8))+20*data111(b1,11)]);
}
}


if (level3){
if ((season==0)||(data11(b1,11)==season-1)){
out3(b1,0)=0;
out3(b1,1)=0;
double h1;
double hisus=exp(data11(b1,8)*para[hai_start+2*data11(b1,11)]);
for (b2=data11(b1,3)+1;b2<=data11(b1,4);++b2){
if ((data21(b1,3)==0)||(b2<=data21(b1,4))){
h1=hisus*ILI(b2-1,data111(b1,12))*para[42+data21(b1,2)+3*data111(b1,11)]*((data21(b1,3)==0)||(b2<=data21(b1,4))); //*pow(3,b2>509); //+3*(b2>509)
if (b2==data21(b1,4)){
out3(b1,1)=log(1-exp(-h1));
}
else{
out3(b1,0)-=h1;
}
}	
}
}
}

}
}
};

//########################################################################################################################################
// function to likelihood
// given the full data set, compute the likelihood
// [[Rcpp::export]]
List loglik(NumericMatrix data11,
NumericMatrix data111,
NumericMatrix data21,
NumericMatrix ILI,
NumericVector para,
NumericVector para2,
int level1,
int level2,
int level3,
int season,
NumericMatrix blankmatrix,
int hai_start){
// check if the para are within their possible range
NumericMatrix out1(data11.nrow(),3);
NumericMatrix out2(data11.nrow(),3);
NumericMatrix out3=clone(blankmatrix);
// call parallel program
LogLik loglik(out1,out2,out3,data11,data111,data21,ILI,para,para2,level1,level2,level3,season,hai_start);
// call parallelFor to do the work
parallelFor(0,data11.nrow(),loglik);

return List::create(_[""]=out1,
	_[""]=out2,
	_[""]=out3);
}






//########################################################################################################################################
//########################################################################################################################################
// the main body of the parallel function
struct AllUpdate:public Worker{
// source vector
RMatrix<double> data11out;
RMatrix<double> data21out;
RMatrix<double> data11pro;
RMatrix<double> data21pro;
RMatrix<double> loglik1out;
RMatrix<double> loglik2out;
RMatrix<double> loglik3out;
RMatrix<double> loglik1pro;
RMatrix<double> loglik2pro;
RMatrix<double> loglik3pro;
RMatrix<double> mcmcrecord;
RMatrix<double> data11;
RMatrix<double> data111;
RMatrix<double> data21;
RMatrix<double> ILI;
RVector<double> para;
RVector<double> para2;
RMatrix<double> loglik1;
RMatrix<double> loglik2;
RMatrix<double> loglik3;
RMatrix<double> temprecord;
int hai_start;
// destination vector
// initialize with source and destination
AllUpdate(NumericMatrix data11out,
NumericMatrix data21out,
NumericMatrix data11pro,
NumericMatrix data21pro,
NumericMatrix loglik1out,
NumericMatrix loglik2out,
NumericMatrix loglik3out,
NumericMatrix loglik1pro,
NumericMatrix loglik2pro,
NumericMatrix loglik3pro,
NumericMatrix mcmcrecord,
NumericMatrix data11,
NumericMatrix data111,
NumericMatrix data21,
NumericMatrix ILI,
NumericVector para,
NumericVector para2,
NumericMatrix loglik1,
NumericMatrix loglik2,
NumericMatrix loglik3,
NumericMatrix temprecord,
int hai_start)
:data11out(data11out),data21out(data21out),data11pro(data11pro),data21pro(data21pro),loglik1out(loglik1out),loglik2out(loglik2out),loglik3out(loglik3out),loglik1pro(loglik1pro),loglik2pro(loglik2pro),loglik3pro(loglik3pro),mcmcrecord(mcmcrecord),data11(data11),data111(data111),data21(data21),ILI(ILI),para(para),para2(para2),loglik1(loglik1),loglik2(loglik2),loglik3(loglik3),temprecord(temprecord),hai_start(hai_start){}

void operator()(std::size_t begin, std::size_t end) {

// section to write the parallel version
// functor (pass input and output matrixes)
for (unsigned int b1=begin;b1<end;++b1){
int b2;
int b3;
//int b4;
//int b5;

// here need to design how many update in the step
// 1-4: baseline AT, inf time, boosting and waning for sh1
int ifupdate[4];
int numberifupdate=2;
ifupdate[0]=1;

if (data21(b1,3)==1){
ifupdate[1]=1;
ifupdate[2]=1;		
numberifupdate+=2;
}
else{
ifupdate[1]=0;
ifupdate[2]=0;		
}
ifupdate[3]=1;	


//ifupdate[1]=0;
//ifupdate[2]=0;
//ifupdate[3]=0;


int totalupdate=0;
int getnumber=1+floor((pow(2,numberifupdate)-1)*R::runif(0,1));   
for (b2=3;b2>=0;--b2){
if (ifupdate[b2]==1){
ifupdate[b2]=getnumber%2;
getnumber/=2;
++totalupdate;   
}
}
mcmcrecord(b1,4)=totalupdate;


double proratio=0;

// the get the propose infection time
// first declare the length
if (ifupdate[1]==1){	
int problength=data111(b1,4);

double prob[problength];
double totalprob=0;

for (b2=problength-1;b2>=data111(b1,3);--b2){
prob[b2]=ILI(b2,data111(b1,12));
totalprob+=prob[b2];
}


// here propose the infection time

double gen=R::runif(0,1);
gen*=totalprob;
for (b2=problength-1;b2>=data111(b1,3);--b2){
gen-=prob[b2];
if (gen<0){
data21pro(b1,4)=b2+1;
// break out from the loop 
break;
}
}

proratio+=log(ILI(data21pro(b1,4)-1,data111(b1,12)))-log(ILI(data21(b1,4)-1,data111(b1,12)));

}



// baseline AT update		
// here the proposal distribution based on the P(baseline) to increase the acceptance 
// then the proposal and likelihood will cancel out	
if (ifupdate[0]==1){	
double titerlevel1=R::runif(0,1);
for (b3=9;b3>=0;--b3){
titerlevel1-=para2[b3+10*(data21(b1,2)>=1)+20*data111(b1,11)];
if (titerlevel1<0){
data11pro(b1,8)=b3;
// break out from the loop 
break;
}
}
data11pro(b1,8)+=R::runif(0,1);
}

// boosting/cross boosting update	
if (ifupdate[2]==0){	
data21pro(b1,5)=data21(b1,5); 
}
// waning
if (ifupdate[3]==0){
data21pro(b1,6)=data21(b1,6); 
}


// here based on the new inf time, recompute the AT dynamic
// based on the infection status, generate the antibody respose
for (b2=0;b2<=1;++b2){
// first check if imputing is necessary
if (data11pro(b1,6+b2)!=-1){
// first check if there is infection with possible boosting in that interval
int cond1=(data11pro(b1,5+b2)<data21pro(b1,4))&&(data21pro(b1,4)<=data11pro(b1,6+b2));

if (cond1==0){
data11pro(b1,9+b2)=data11pro(b1,8+b2)*exp(-(data11pro(b1,6+b2)-data11pro(b1,5+b2))*data21pro(b1,6)/365.0);
}
else{
// it is infection 
double beforeinf=data11pro(b1,8+b2)*exp(-(data21pro(b1,4)-data11pro(b1,5+b2))*data21pro(b1,6)/365.0);	
data11pro(b1,9+b2)=(beforeinf+data21pro(b1,5))*exp(-(data11pro(b1,6+b2)-data21pro(b1,4)-1)*data21pro(b1,6)/365.0);
}
data11pro(b1,9+b2)=data11pro(b1,9+b2)*(data11pro(b1,9+b2)<17)*(data11pro(b1,9+b2)>0)+17*(data11pro(b1,9+b2)>17);
}
}


// here compute the likelihood for propose

// first level to record the measurment error likelihood
for (b2=2;b2>=0;--b2){
if (data111(b1,8+b2)!=-1){

double probtemp[20];
double probtemptotal=0;

for (b3=9;b3>=0;--b3){
// the random error
probtemp[b3]=para[0];
probtemp[10+b3]=0.000001;
}

double measured=data11pro(b1,8+b2); //+para[0]*data11(b1,12);

if (measured<1){
// two-fold error at 0 
probtemp[0]+=1.0-measured/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[1]+=measured/(2.0*exp(para[1+1*(data11(b1,12))]));
}
else{
// two-fold error at other step	
int index=floor(measured);
probtemp[index]+=1.0-(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]))-(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index-1]+=(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index+1]+=(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
}

for (b3=19;b3>=0;--b3){	
probtemptotal+=probtemp[b3];
//if (b2==0){
//record(b1,b3)=probtemp[b3];
//}
}

loglik1pro(b1,b2)=log(probtemp[int(data111(b1,8+b2))])-log(probtemptotal);

}

}

// second level
// everything is gibbs sampler, not need this level

// third level
loglik3pro(b1,0)=0;
loglik3pro(b1,1)=0;
double h1;
double hisus=exp(data11pro(b1,8)*para[hai_start+2*data11(b1,11)]);
for (b2=data11pro(b1,3)+1;b2<=data11pro(b1,4);++b2){
if ((data21pro(b1,3)==0)||(b2<=data21pro(b1,4))){
h1=hisus*ILI(b2-1,data111(b1,12))*para[42+data21pro(b1,2)+3*data111(b1,11)]*((data21pro(b1,3)==0)||(b2<=data21pro(b1,4))); //*pow(3,b2>509); //+3*(b2>509)
if (b2==data21pro(b1,4)){
loglik3pro(b1,1)=log(1-exp(-h1));
}
else{
loglik3pro(b1,0)-=h1;
}
}
}

// here do the metropolis hasting update
double liknew=0;
double likold=0;
for (b2=2;b2>=0;--b2){
liknew+=loglik1pro(b1,b2);
likold+=loglik1(b1,b2);	
}
// don't add level2, beacuse using gibbs sampler
for (b2=1;b2>=0;--b2){
liknew+=loglik3pro(b1,b2);
likold+=loglik3(b1,b2);	
}


double loglikratio=liknew-likold;
double accept_pro=pow(exp(1),loglikratio-proratio);
if (data21pro(b1,6)<pow(10,-150)){
accept_pro=0;
}
mcmcrecord(b1,1)=accept_pro;
mcmcrecord(b1,2)=loglikratio;	
mcmcrecord(b1,3)=proratio;
if (gen_binom(accept_pro)){
mcmcrecord(b1,0)=1;	
// if accept, make the out to the same as the proposal
for (b2=10;b2>=8;--b2){
data11out(b1,b2)=data11pro(b1,b2);
}
for (b2=6;b2>=4;--b2){
data21out(b1,b2)=data21pro(b1,b2);
}
for (b2=2;b2>=0;--b2){
loglik1out(b1,b2)=loglik1pro(b1,b2);		
}
for (b2=1;b2>=0;--b2){
loglik3out(b1,b2)=loglik3pro(b1,b2);	
}
}
else{
mcmcrecord(b1,0)=-1;	
}


// here need to modify the output 
// because using gibbs sampler, so the P(AT) is ignored
// here need to put back the baseline AT titer
if (data111(b1,8)!=-1){
loglik2out(b1,2)=log(para2[floor(data11out(b1,8))+10*(data21(b1,2)>=1)+20*data111(b1,11)]);	
}
if (data21out(b1,3)==1){
loglik2out(b1,0)=R::dgamma(data21out(b1,5),para[18+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0,1);
loglik2out(b1,1)=R::dgamma(data21out(b1,6),para[19+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0,1);
}
else{
loglik2out(b1,0)=0;
loglik2out(b1,1)=0;	
}


}
}
};


//########################################################################################################################################
// function to update infection time, boosting/waning parameter, baseline AT titer
// but not the infection status
// [[Rcpp::export]]
List all_update(NumericMatrix data11,
NumericMatrix data111,
NumericMatrix data21,
NumericMatrix ILI,
NumericVector para,
NumericVector para2,
NumericMatrix loglik1,
NumericMatrix loglik2,
NumericMatrix loglik3,
int hai_start){

NumericMatrix data11pro(clone(data11));
NumericMatrix data21pro(clone(data21));
NumericMatrix data11out(clone(data11));
NumericMatrix data21out(clone(data21));
NumericMatrix loglik1pro(clone(loglik1));
NumericMatrix loglik2pro(clone(loglik2));
NumericMatrix loglik3pro(clone(loglik3));
NumericMatrix loglik1out(clone(loglik1));
NumericMatrix loglik2out(clone(loglik2));
NumericMatrix loglik3out(clone(loglik3));
// 1.accept/reject
NumericMatrix mcmcrecord(data11.nrow(),10);
NumericMatrix temprecord(data11.nrow(),10);

int b1;



// here put the boosting waning etc
for (b1=data11.nrow()-1;b1>=0;--b1){  	
// first impute the basic information here
data21pro(b1,5)=R::rgamma(para[18+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);
data21pro(b1,6)=R::rgamma(para[19+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);
}

// call parallel program
AllUpdate allupdate(data11out,data21out,data11pro,data21pro,loglik1out,loglik2out,loglik3out,loglik1pro,loglik2pro,loglik3pro,mcmcrecord,data11,data111,data21,ILI,para,para2,loglik1,loglik2,loglik3,temprecord,hai_start);
// call parallelFor to do the work
parallelFor(0,data21.nrow(),allupdate);


double a1=0;
double a2=0;
for (b1=mcmcrecord.nrow()-1;b1>=0;--b1){
if (mcmcrecord(b1,0)==1){
++a1;
}
if (mcmcrecord(b1,0)!=0){
++a2;	
}
}


return List::create(_[""]=data11out,
_[""]=data21out,
_[""]=loglik1out,
_[""]=loglik2out,
_[""]=loglik3out,
_[""]=a1/a2,
_[""]=mcmcrecord,
_[""]=data11pro,
_[""]=data21pro,
_[""]=loglik1pro,
_[""]=loglik2pro,
_[""]=loglik3pro,
_[""]=temprecord);
}



//########################################################################################################################################
// function to add infection for reversible jump MCMC
// [[Rcpp::export]]
List add_remove_infection(NumericMatrix data11,
NumericMatrix data111,
NumericMatrix data21,
NumericMatrix ILI,
NumericVector para,
NumericVector para2,
NumericMatrix loglik1,
NumericMatrix loglik2,
NumericMatrix loglik3,
int hai_start){

NumericMatrix data11pro(clone(data11));
NumericMatrix data21pro(clone(data21));
NumericMatrix data11out(clone(data11));
NumericMatrix data21out(clone(data21));
NumericMatrix loglik1pro(clone(loglik1));
NumericMatrix loglik2pro(clone(loglik2));
NumericMatrix loglik3pro(clone(loglik3));
NumericMatrix loglik1out(clone(loglik1));
NumericMatrix loglik2out(clone(loglik2));
NumericMatrix loglik3out(clone(loglik3));
// 1.accept/reject
NumericMatrix mcmcrecord(data11.nrow(),10);
NumericMatrix temprecord(1,1);

int b1;
int b2;
int b3;
//int b4;   

// first compute the number to add or delete
int numberofinfect=0;
int numberofall=data11.nrow();
for (b1=data11.nrow()-1;b1>=0;--b1){
if (data21(b1,3)==1){
++numberofinfect;
}
}

// here to do the add/delete step
int add=gen_binom(0.5);
int select=0;
double proratio=0;
int finalaccept=0;
//######################################################################################################################################################################################################
// try to add infection
if (add){
if (numberofall-numberofinfect>0){
// first get the random number
select=(int)floor(R::runif(0, numberofall-numberofinfect));
// here get the subject and strain
// use b1 to indicate the subject
for (b3=data11.nrow()-1;b3>=0;--b3){
if (data21(b3,3)==0){
--select;
if (select<0){
b1=b3;
break;
}
}
}



// change from 0 to 1, then need to add 
data21pro(b1,3)=1;
data21pro(b1,5)=R::rgamma(para[18+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);
data21pro(b1,6)=R::rgamma(para[19+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);	
//##################################################################
// infection time
int problength=data111(b1,4);

double prob[problength];
double totalprob=0;

for (b2=problength-1;b2>=data111(b1,3);--b2){
prob[b2]=ILI(b2,data111(b1,12));
totalprob+=prob[b2];
}



//for (b2=problength-1;b2>=0;--b2){
//temprecord(b1,b2)=prob[b2];	
//}

// here propose the infection time

double gen=R::runif(0,1);
gen*=totalprob;
for (b2=problength-1;b2>=data111(b1,3);--b2){
gen-=prob[b2];
if (gen<0){
data21pro(b1,4)=b2+1;
// break out from the loop 
break;
}
}

proratio+=(log(ILI(data21pro(b1,4)-1,data111(b1,12)))-log(totalprob));


//##################################################################
// baseline AT titer
double titerlevel1=R::runif(0,1);
for (b3=9;b3>=0;--b3){
titerlevel1-=para2[b3+10*(data21(b1,2)>=1)+20*data111(b1,11)];
if (titerlevel1<0){
data11pro(b1,8)=b3;
// break out from the loop 
break;
}
}
data11pro(b1,8)+=R::runif(0,1);


// here based on the new inf time, recompute the AT dynamic
// based on the infection status, generate the antibody respose
for (b2=0;b2<=1;++b2){
// first check if imputing is necessary
if (data11pro(b1,6+b2)!=-1){
// first check if there is infection with possible boosting in that interval
int cond1=(data11pro(b1,5+b2)<data21pro(b1,4))&&(data21pro(b1,4)<=data11pro(b1,6+b2));
if (cond1==0){
data11pro(b1,9+b2)=data11pro(b1,8+b2)*exp(-(data11pro(b1,6+b2)-data11pro(b1,5+b2))*data21pro(b1,6)/365.0);
}
else {
// it is infection 
double beforeinf=data11pro(b1,8+b2)*exp(-(data21pro(b1,4)-data11pro(b1,5+b2))*data21pro(b1,6)/365.0);
data11pro(b1,9+b2)=(beforeinf+data21pro(b1,5))*exp(-(data11pro(b1,6+b2)-data21pro(b1,4)-1)*data21pro(b1,6)/365.0);
}
data11pro(b1,9+b2)=data11pro(b1,9+b2)*(data11pro(b1,9+b2)<17)*(data11pro(b1,9+b2)>0)+17*(data11pro(b1,9+b2)>17);
}
}


// here compute the likelihood for propose
// first level to record the measurment error likelihood
for (b2=2;b2>=0;--b2){
if (data111(b1,8+b2)!=-1){

double probtemp[20];
double probtemptotal=0;

for (b3=9;b3>=0;--b3){
// the random error
probtemp[b3]=para[0];
probtemp[10+b3]=0.000001;
}

double measured=data11pro(b1,8+b2); //+para[0]*data11(b1,12);

if (measured<1){
// two-fold error at 0 
probtemp[0]+=1.0-measured/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[1]+=measured/(2.0*exp(para[1+1*(data11(b1,12))]));
}
else{
// two-fold error at other step	
int index=floor(measured);
probtemp[index]+=1.0-(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]))-(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index-1]+=(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index+1]+=(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
}

for (b3=19;b3>=0;--b3){	
probtemptotal+=probtemp[b3];
//if (b2==0){
//record(b1,b3)=probtemp[b3];
//}
}
loglik1pro(b1,b2)=log(probtemp[int(data111(b1,8+b2))])-log(probtemptotal);

}

}

// second level 
// in this update nothing from level 2 need to be added
// because they are all using gibbs sampler


// third level
loglik3pro(b1,0)=0;
loglik3pro(b1,1)=0;
double h1;
double hisus=exp(data11pro(b1,8)*para[hai_start+2*data11(b1,11)]);

for (b2=data11pro(b1,3)+1;b2<=data11pro(b1,4);++b2){
if ((data21pro(b1,3)==0)||(b2<=data21pro(b1,4))){
h1=hisus*ILI(b2-1,data111(b1,12))*para[42+data11pro(b1,2)+3*data111(b1,11)]*((data21pro(b1,3)==0)||(b2<=data21pro(b1,4))); //*pow(3,b2>509); //+3*(b2>509)
if (b2==data21pro(b1,4)){
loglik3pro(b1,1)=log(1-exp(-h1));
}
else{
loglik3pro(b1,0)-=h1;
}
}
}

// here do the metropolis hasting update
double liknew=0;
double likold=0;
for (b2=2;b2>=0;--b2){
liknew+=loglik1pro(b1,b2);
likold+=loglik1(b1,b2);
}
// don't add level 2
for (b2=1;b2>=0;--b2){
liknew+=loglik3pro(b1,b2);
likold+=loglik3(b1,b2);
}


double loglikratio=liknew-likold;
double numberratio=(numberofall-numberofinfect)/(numberofinfect+1.0);
double accept_pro=pow(exp(1),loglikratio-proratio)*numberratio;
if (data21pro(b1,6)<pow(10,-150)){
accept_pro=0;	
}
mcmcrecord(b1,1)=accept_pro;
mcmcrecord(b1,2)=loglikratio;	
mcmcrecord(b1,3)=proratio;
mcmcrecord(b1,6)=numberratio;
if (gen_binom(accept_pro)){
finalaccept=1;	
mcmcrecord(b1,0)=1;	
// if accept, make the out to the same as the proposal
for (b2=6;b2>=3;--b2){
data21out(b1,b2)=data21pro(b1,b2);
}
for (b2=10;b2>=8;--b2){
data11out(b1,b2)=data11pro(b1,b2);
}
for (b2=2;b2>=0;--b2){
loglik1out(b1,b2)=loglik1pro(b1,b2);		
}
for (b2=1;b2>=0;--b2){
loglik3out(b1,b2)=loglik3pro(b1,b2);	
}
 
loglik2out(b1,0)=R::dgamma(data21out(b1,5),para[18+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0,1);
loglik2out(b1,1)=R::dgamma(data21out(b1,6),para[19+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0,1);



}
else{
mcmcrecord(b1,0)=-1;	
}
// because using gibbs sampler, so the P(AT) is ignored
// here need to put back the baseline AT titer
if (data111(b1,8)!=-1){
loglik2out(b1,2)=log(para2[floor(data11out(b1,8))+10*(data21(b1,2)>=1)+20*data111(b1,11)]);	
}

mcmcrecord(b1,5)=data21pro(b1,4);



}
}

//######################################################################################################################################################################################################

// try to remove infection
else{
if (numberofinfect>0){
// first get the random number
select=(int)floor(R::runif(0, numberofinfect));
// here get the subject and strain
// use b1 to indicate the subject
for (b3=data11.nrow()-1;b3>=0;--b3){
if (data21(b3,3)==1){
--select;
if (select<0){
b1=b3;
break;
}
}
}	




data21pro(b1,3)=0;	
data21pro(b1,4)=0;
data21pro(b1,6)=R::rgamma(para[19+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);	

//##################################################################
// infection time
int problength=data111(b1,4);

double prob[problength];
double totalprob=0;

for (b2=problength-1;b2>=data111(b1,3);--b2){
prob[b2]=ILI(b2,data111(b1,12));
totalprob+=prob[b2];
}

// here need this to compute the proposal ratio


proratio-=(log(ILI(data21(b1,4)-1,data111(b1,12)))-log(totalprob));

//##################################################################
// baseline AT titer
double titerlevel1=R::runif(0,1);
for (b3=9;b3>=0;--b3){
titerlevel1-=para2[b3+10*(data21(b1,2)>=1)+20*data111(b1,11)];
if (titerlevel1<0){
data11pro(b1,8)=b3;
// break out from the loop 
break;
}
}
data11pro(b1,8)+=R::runif(0,1);

// here based on the new inf time, recompute the AT dynamic
// based on the infection status, generate the antibody respose
for (b2=0;b2<=1;++b2){
// first check if imputing is necessary
if (data11pro(b1,6+b2)!=-1){
// first check if there is infection with possible boosting in that interval
int cond1=(data11pro(b1,5+b2)<data21pro(b1,4))&&(data21pro(b1,4)<=data11pro(b1,6+b2));
if (cond1==0){
data11pro(b1,9+b2)=data11pro(b1,8+b2)*exp(-(data11pro(b1,6+b2)-data11pro(b1,5+b2))*data21pro(b1,6)/365.0);
}
else {
// it is infection 
double beforeinf=data11pro(b1,8+b2)*exp(-(data21pro(b1,4)-data11pro(b1,5+b2))*data21pro(b1,6)/365.0);
data11pro(b1,9+b2)=(beforeinf+data21pro(b1,5))*exp(-(data11pro(b1,6+b2)-data21pro(b1,4)-1)*data21pro(b1,6)/365.0);
}
data11pro(b1,9+b2)=data11pro(b1,9+b2)*(data11pro(b1,9+b2)<17)*(data11pro(b1,9+b2)>0)+17*(data11pro(b1,9+b2)>17);
}
}


// here compute the likelihood for propose
// first level to record the measurment error likelihood
for (b2=2;b2>=0;--b2){

if (data111(b1,8+b2)!=-1){

double probtemp[20];
double probtemptotal=0;

for (b3=9;b3>=0;--b3){
// the random error
probtemp[b3]=para[0];
}



for (b3=9;b3>=0;--b3){
// the random error
probtemp[b3]=para[0];
probtemp[10+b3]=0.000001;
}

double measured=data11pro(b1,8+b2); //+para[0]*data11(b1,12);

if (measured<1){
// two-fold error at 0 
probtemp[0]+=1.0-measured/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[1]+=measured/(2.0*exp(para[1+1*(data11(b1,12))]));
}
else{
// two-fold error at other step	
int index=floor(measured);
probtemp[index]+=1.0-(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]))-(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index-1]+=(1.0-(measured-floor(measured)))/(2.0*exp(para[1+1*(data11(b1,12))]));
probtemp[index+1]+=(measured-floor(measured))/(2.0*exp(para[1+1*(data11(b1,12))]));
}

for (b3=19;b3>=0;--b3){	
probtemptotal+=probtemp[b3];
//if (b2==0){
//record(b1,b3)=probtemp[b3];
//}
}
loglik1pro(b1,b2)=log(probtemp[int(data111(b1,8+b2))])-log(probtemptotal);

}
}

// second level
// in this update nothing from level 2 need to be added
// because they are all using gibbs sampler

// third level
loglik3pro(b1,0)=0;
loglik3pro(b1,1)=0;
double h1;
double hisus=exp(data11pro(b1,8)*para[hai_start+2*data11(b1,11)]);

for (b2=data11pro(b1,3)+1;b2<=data11pro(b1,4);++b2){
if ((data21pro(b1,3)==0)||(b2<=data21pro(b1,4))){
h1=hisus*ILI(b2-1,data111(b1,12))*para[42+data11pro(b1,2)+3*data111(b1,11)]*((data21pro(b1,3)==0)||(b2<=data21pro(b1,4))); //*pow(3,b2>509); //+3*(b2>509)
if (b2==data21pro(b1,4)){
loglik3pro(b1,1)=log(1-exp(-h1));
}
else{
loglik3pro(b1,0)-=h1;
}
}
}

// here do the metropolis hasting update
double liknew=0;
double likold=0;
for (b2=2;b2>=0;--b2){
liknew+=loglik1pro(b1,b2);
likold+=loglik1(b1,b2);
}
// don't add level 2
for (b2=1;b2>=0;--b2){
liknew+=loglik3pro(b1,b2);
likold+=loglik3(b1,b2);
}

double loglikratio=liknew-likold;
double numberratio=numberofinfect/(numberofall-numberofinfect+1.0);
double accept_pro=pow(exp(1),loglikratio-proratio)*numberratio;
if (data21pro(b1,6)<pow(10,-150)){
accept_pro=0;	
}
mcmcrecord(b1,1)=accept_pro;
mcmcrecord(b1,2)=loglikratio;	
mcmcrecord(b1,3)=proratio;
mcmcrecord(b1,6)=numberratio;
if (gen_binom(accept_pro)){
finalaccept=1;	
mcmcrecord(b1,0)=1;	
// if accept, make the out to the same as the proposal
for (b2=10;b2>=8;--b2){
data11out(b1,b2)=data11pro(b1,b2);
}
for (b2=6;b2>=3;--b2){
data21out(b1,b2)=data21pro(b1,b2);
}
for (b2=2;b2>=0;--b2){
loglik1out(b1,b2)=loglik1pro(b1,b2);		
}
for (b2=1;b2>=0;--b2){
loglik3out(b1,b2)=loglik3pro(b1,b2);	
}

// no longer contribute to boosting or waning
loglik2out(b1,0)=0;
loglik2out(b1,1)=0;
  
}
else{
mcmcrecord(b1,0)=-1;	
}
// because using gibbs sampler, so the P(AT) is ignored
// here need to put back the baseline AT titer
if (data111(b1,8)!=-1){
loglik2out(b1,2)=log(para2[floor(data11out(b1,8))+10*(data21(b1,2)>=1)+20*data111(b1,11)]);	
}


}
}






return List::create(_[""]=data11out,
_[""]=data21out,
_[""]=loglik1out,
_[""]=loglik2out,
_[""]=loglik3out,
_[""]=finalaccept,
_[""]=mcmcrecord,
_[""]=data11pro,
_[""]=data21pro,
_[""]=loglik1pro,
_[""]=loglik2pro,
_[""]=loglik3pro,
_[""]=temprecord,
_[""]=numberofinfect,
_[""]=numberofall);
}











//##############################################################################################################################################
//##############################################################################################################################################
// function for mcmc
// [[Rcpp::export]]
List mcmc(NumericMatrix input1,
NumericMatrix input2,
NumericMatrix input3,
NumericMatrix ILI,
int mcmc_n,             // length of mcmc stain
NumericVector int_para, // initial parameter
NumericVector int_para2,
NumericVector int_para3,
NumericVector paraseason,
NumericVector move,     // which one should move in the model
NumericVector sigma,
NumericVector sigma3,
int burnin,
int thinning,
int n_seasons){

// compute hai_start once
int hai_start = 42 + 3 * n_seasons;

// create the vector for use
int b0;
int b1;
int b2;
//int b3;
//int b4;
int moveindex;

// here data111 is the input that can not be changed
// data21 is those underlying variables
// data11 is computed from data21
NumericMatrix data11(clone(input1));
NumericMatrix data111(clone(input2));
NumericMatrix data21(clone(input3));


// matrix to record LL
// need to set number of parameter here
NumericMatrix p_para(mcmc_n,int_para.length());
NumericMatrix p_para_r(mcmc_n,sum(move));
p_para(0,_)=int_para;
moveindex=sum(move)-1;
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){
p_para_r(0,moveindex)=p_para(0,b1);
--moveindex;
}	
}

// for the second parameter vectors
NumericMatrix baselinecount(mcmc_n,int_para2.length());
NumericMatrix p_para2(mcmc_n,int_para2.length());
p_para2(0,_)=int_para2;

// for the third parameter vectors
NumericMatrix p_para3(mcmc_n,int_para3.length());
p_para3(0,_)=int_para3;


// here set the initial point for missing
// first assume there is no missing infection
// i.e. all -1 mean no infection

for (b1=data21.nrow()-1;b1>=0;--b1){

// the initial infection status is put in the input file
data21(b1,3)=gen_binom(0.2);	
data21(b1,4)=0;	
// here set the initial point for the infection time
if (data21(b1,3)==1){
data21(b1,4)=floor((data111(b1,3)+data111(b1,4))/2);
}

// here set the initial point for the baseline titer
if (data111(b1,8)>=0){
data11(b1,8)=data111(b1,8)+R::runif(0,1);
}
if (data111(b1,8)==-1){
data11(b1,8)=data111(b1,9)+R::runif(0,1);
}


// update boosting and waning
data21(b1,5)=R::rgamma(int_para[18+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);
data21(b1,6)=R::rgamma(int_para[19+4*(data11(b1,13))+2*(data11(b1,2)>0)],1.0);
// here base on the input, update the data1
// based on the infection status, generate the antibody respose
for (b2=0;b2<=1;++b2){
// first check if imputing is necessary
if (data11(b1,6+b2)!=-1){
// first check if there is infection with possible boosting in that interval
int cond1=(data11(b1,5+b2)<data21(b1,4))&&(data21(b1,4)<=data11(b1,6+b2));
if (cond1==0){
data11(b1,9+b2)=data11(b1,8+b2)*exp(-(data11(b1,6+b2)-data11(b1,5+b2))*data21(b1,6)/365.0);
}
else {
// it is infection 
double beforeinf=data11(b1,8+b2)*exp(-(data21(b1,4)-data11(b1,5+b2))*data21(b1,6)/365.0);
data11(b1,9+b2)=(beforeinf+data21(b1,5))*exp(-(data11(b1,6+b2)-data21(b1,4)-1)*data21(b1,6)/365.0);
}
data11(b1,9+b2)=data11(b1,9+b2)*(data11(b1,9+b2)<10)*(data11(b1,9+b2)>0)+10*(data11(b1,9+b2)>10);
}
}

}



// first row is the overall matrix, other three row is the indiviudal likelihood
NumericVector acceptrate(int_para.length());
NumericVector acceptrate3(int_para3.length());
NumericMatrix LL1(mcmc_n,4); // the original likelihood
NumericMatrix LL2(mcmc_n,int_para3.length()); // the likelihood for hyperparameter
NumericMatrix LL3(mcmc_n,4); // for the update except the infeciton status
NumericMatrix LL4(mcmc_n,4);
NumericMatrix LL5(mcmc_n,4);
NumericMatrix totalinfectionnumber(mcmc_n,1);
NumericMatrix updateacceptrate(mcmc_n,2);
// here select the first 20 people to record
int rownumber=floor((mcmc_n-burnin)/thinning);
NumericMatrix impute_record_titer1(rownumber,data11.nrow());
IntegerMatrix impute_record_infstatus1(rownumber,data11.nrow());
IntegerMatrix impute_record_inftime1(rownumber,data11.nrow());
NumericMatrix impute_record_inftime1pro(1,1);
NumericMatrix impute_record_waning1(rownumber,data11.nrow());
NumericMatrix impute_record_boosting1(rownumber,data11.nrow());
//####################################################################################################################################
// initial step

//####################################################################################################################################
// compute likelihood

NumericMatrix blankmatrix(data11.nrow(),2);

List loglikall=loglik(data11,data111,data21,ILI,p_para(0,_),p_para2(0,_),1,1,1,0,blankmatrix,hai_start);
List loglikallpro;
NumericMatrix loglik1=loglikall(0);
NumericMatrix loglik2=loglikall(1);
NumericMatrix loglik3=loglikall(2);
NumericMatrix loglik1pro;
NumericMatrix loglik2pro;
NumericMatrix loglik3pro;

LL1(0,1)=sum(loglik1);
LL1(0,2)=sum(loglik2);
LL1(0,3)=sum(loglik3);
LL1(0,0)=LL1(0,1)+LL1(0,2)+LL1(0,3);

NumericVector temploglik(4);
NumericVector newloglik(4);
temploglik(0)=LL1(0,0)+prior_loglik(p_para(0,_),n_seasons);
temploglik(1)=LL1(0,1);
temploglik(2)=LL1(0,2);
temploglik(3)=LL1(0,3);

// here also need to fill the LL3
for (b1=int_para3.length()-1;b1>=0;--b1){	
LL2(0,b1)=log(tgamma(10*p_para3(0,b1)))-10*log(tgamma(p_para3(0,b1)));
for (b2=9;b2>=0;--b2){
LL2(0,b1)+=(p_para3(0,b1)-1)*log(p_para2(0,10*b1+b2));
}
}



double loglikeratio;
double accept_pro;
NumericVector pro_para(int_para.length());

totalinfectionnumber(0,0)=sum(data21(_,3));

int recordcurrentindex=0;
//####################################################################################################################################
// main mcmc step

//####################################################################################################################################
for (b0=1;b0<mcmc_n;++b0){


// after 500 step, then set the sigma to be the empirical sigma
if ((b0>500)&&(b0%200==0)){
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){	
NumericVector temp1(b0-1);
for (b2=b0-2;b2>=0;--b2){
temp1(b2)=p_para(b2,b1);	
}
sigma(b1)=sd(temp1);
// tuning
if (acceptrate(b1)<0.1){
sigma(b1)*=0.5;
}	
if ((acceptrate(b1)<0.15)&&(acceptrate(b1)>0.1)){
sigma(b1)*=0.8;
}
if ((acceptrate(b1)<0.2)&&(acceptrate(b1)>0.15)){
sigma(b1)*=0.95;
}
if ((acceptrate(b1)<0.4)&&(acceptrate(b1)>0.3)){
sigma(b1)*=1.05;
}
if ((acceptrate(b1)<0.9)&&(acceptrate(b1)>0.4)){
sigma(b1)*=1.2;
}
if (acceptrate(b1)>0.9){
sigma(b1)*=2;
}
}
}
// for sigma3

for (b1=int_para3.length()-1;b1>=0;--b1){
NumericVector temp1(b0-1);
for (b2=b0-2;b2>=0;--b2){
temp1(b2)=p_para3(b2,b1);	
}
sigma3(b1)=sd(temp1);
// tuning
if (acceptrate3(b1)<0.1){
sigma3(b1)*=0.5;
}	
if ((acceptrate3(b1)<0.15)&&(acceptrate3(b1)>0.1)){
sigma3(b1)*=0.8;
}
if ((acceptrate3(b1)<0.2)&&(acceptrate3(b1)>0.15)){
sigma3(b1)*=0.95;
}
if ((acceptrate3(b1)<0.4)&&(acceptrate3(b1)>0.3)){
sigma3(b1)*=1.05;
}
if ((acceptrate3(b1)<0.9)&&(acceptrate3(b1)>0.4)){
sigma3(b1)*=1.2;
}
if (acceptrate3(b1)>0.9){
sigma3(b1)*=2;
}
}

}


// metorpolis-hasing update on parameter
for (b1=0;b1<int_para.length();++b1){
if (move(b1)){
int level=2;
if (b1>=42){
level=3;	
}
if (b1<=17){ 
level=1;
}
pro_para=p_para(b0-1,_);
for (b2=b1-1;b2>=0;--b2){
pro_para(b2)=p_para(b0,b2);	
}
pro_para(b1)+=rnorm(0.0,sigma(b1));
newloglik(0)=prior_loglik(pro_para,n_seasons);
if (newloglik(0)> -9999999){
// level 3 
if (level==3){
loglikallpro=loglik(data11,data111,data21,ILI,pro_para,p_para2(b0-1,_),0,0,1,paraseason(b1),loglik3,hai_start);
NumericMatrix tempoutput=loglikallpro(2);
loglik3pro=clone(tempoutput);
newloglik(1)=temploglik(1);
newloglik(2)=temploglik(2);
newloglik(3)=sum(loglik3pro);
}
if (level==1){ // level1
loglikallpro=loglik(data11,data111,data21,ILI,pro_para,p_para2(b0-1,_),1,0,0,0,loglik3,hai_start);
NumericMatrix tempoutput=loglikallpro(0);
loglik1pro=clone(tempoutput);
newloglik(1)=sum(loglik1pro);
newloglik(2)=temploglik(2);
newloglik(3)=temploglik(3);
}
if (level==2){  //level2
loglikallpro=loglik(data11,data111,data21,ILI,pro_para,p_para2(b0-1,_),0,1,0,0,loglik3,hai_start);
NumericMatrix tempoutput=loglikallpro(1);
loglik2pro=clone(tempoutput);
newloglik(1)=temploglik(1);
newloglik(2)=sum(loglik2pro);
newloglik(3)=temploglik(3);
}
newloglik(0)+=newloglik(1)+newloglik(2)+newloglik(3);
loglikeratio=newloglik(0)-temploglik(0);
accept_pro=pow(exp(1),loglikeratio);
//if (b1==0){
//record(b0,0)=loglikeratio;
//record(b0,1)=newloglik(0);	
//record(b0,2)=newloglik(1);	
//record(b0,3)=newloglik(2);	
//record(b0,4)=newloglik(3);	
//record(b0,5)=temploglik(0);	
//record(b0,6)=temploglik(1);	
//record(b0,7)=temploglik(2);	
//record(b0,8)=temploglik(3);		
//}
}
else{
accept_pro=0;	
}
if(gen_binom(accept_pro)){
if (level==1){
loglik1=clone(loglik1pro);	
temploglik(1)=newloglik(1);
}
if (level==2){
loglik2=clone(loglik2pro);		
temploglik(2)=newloglik(2);
}
if (level==3){
loglik3=clone(loglik3pro);		
temploglik(3)=newloglik(3);	
}
p_para(b0,b1)=pro_para(b1);
temploglik(0)=newloglik(0);
acceptrate(b1)*=(b0-1);
acceptrate(b1)+=1;
acceptrate(b1)/=b0;
}
else{
p_para(b0,b1)=p_para(b0-1,b1);
acceptrate(b1)*=(b0-1);
acceptrate(b1)/=b0;
}
}
else {
p_para(b0,b1)=p_para(b0-1,b1);
}
}

LL1(b0,0)=temploglik(0)-prior_loglik(p_para(b0,_),n_seasons);
LL1(b0,1)=temploglik(1);
LL1(b0,2)=temploglik(2);
LL1(b0,3)=temploglik(3);



// here to update p_para2
for (b1=baselinecount.ncol()-1;b1>=0;--b1){
baselinecount(b0,b1)=0;	
} 

for (b1=data11.nrow()-1;b1>=0;--b1){
if (data111(b1,8)!=-1){
++baselinecount(b0,10*(data21(b1,2)>=1)+floor(data11(b1,8))+20*data111(b1,11));
}
}

NumericVector tempcount(int_para3.length());
for (b1=int_para3.length()-1;b1>=0;--b1){
tempcount(b1)=0;	
}
for (b1=baselinecount.ncol()-1;b1>=0;--b1){
p_para2(b0,b1)=R::rgamma(baselinecount(b0,b1)+p_para3(b0-1,b1/10),1.0);
tempcount(b1/10)+=p_para2(b0,b1);
}
for (b1=baselinecount.ncol()-1;b1>=0;--b1){
p_para2(b0,b1)/=tempcount(b1/10);
}

// here need to update the likelihood
loglikallpro=loglik(data11,data111,data21,ILI,p_para(b0,_),p_para2(b0,_),0,1,0,0,loglik3,hai_start);
NumericMatrix tempoutput=loglikallpro(1);
loglik2pro=clone(tempoutput);
temploglik(2)=sum(loglik2pro);


// here to update p_para3
// need to write the dirichlet likelihood here
// here conduct metropolis hasting
//p_para3(b0,_)=int_para3;

double hyperoldlik=0;
double hypernewlik=0;

for (b1=int_para3.length()-1;b1>=0;--b1){
p_para3(b0,b1)=p_para3(b0-1,b1)+rnorm(0.0,sigma3(b1));
hyperoldlik=log(tgamma(10*p_para3(b0-1,b1)))-10*log(tgamma(p_para3(b0-1,b1)));
for (b2=9;b2>=0;--b2){
hyperoldlik+=(p_para3(b0-1,b1)-1)*log(p_para2(b0,10*b1+b2));
}
if (p_para3(b0,b1)<0){
accept_pro=0;
}
else{
hypernewlik=log(tgamma(10*p_para3(b0,b1)))-10*log(tgamma(p_para3(b0,b1)));
for (b2=9;b2>=0;--b2){
hypernewlik+=(p_para3(b0,b1)-1)*log(p_para2(b0,10*b1+b2));
}
loglikeratio=hypernewlik-hyperoldlik;
accept_pro=pow(exp(1),loglikeratio);
}
// here do the metropolis hasting
// accept
if(gen_binom(accept_pro)){
acceptrate3(b1)*=(b0-1);
acceptrate3(b1)+=1;
acceptrate3(b1)/=b0;
LL2(b0,b1)=hypernewlik;
}
else{
LL2(b0,b1)=hyperoldlik;
p_para3(b0,b1)=p_para3(b0-1,b1);
acceptrate3(b1)*=(b0-1);
acceptrate3(b1)/=b0;
}
} 



// update the boosting waning parameter without changing infection status
List allupdate=all_update(data11,data111,data21,ILI,p_para(b0,_),p_para2(b0,_),loglik1,loglik2,loglik3,hai_start);
NumericMatrix tempoutput1=allupdate(0);
data11=clone(tempoutput1);
NumericMatrix tempoutput2=allupdate(1);
data21=clone(tempoutput2);
NumericMatrix tempoutput3=allupdate(2);
loglik1=clone(tempoutput3);
NumericMatrix tempoutput4=allupdate(3);
loglik2=clone(tempoutput4);
NumericMatrix tempoutput5=allupdate(4);
loglik3=clone(tempoutput5);
updateacceptrate(b0,0)=allupdate(5);

// update the impute likelihood
LL3(b0,1)=sum(loglik1);
LL3(b0,2)=sum(loglik2);
LL3(b0,3)=sum(loglik3);
LL3(b0,0)=LL3(b0,1)+LL3(b0,2)+LL3(b0,3);

temploglik(0)=LL3(b0,0)+prior_loglik(p_para(b0,_),n_seasons);
temploglik(1)=LL3(b0,1);
temploglik(2)=LL3(b0,2);
temploglik(3)=LL3(b0,3);



// add/remove infection step
for (b1=18;b1>=0;--b1){
List addremoveinfection=add_remove_infection(data11,data111,data21,ILI,p_para(b0,_),p_para2(b0,_),loglik1,loglik2,loglik3,hai_start);
NumericMatrix tempoutput11=addremoveinfection(0);
data11=clone(tempoutput11);
NumericMatrix tempoutput12=addremoveinfection(1);
data21=clone(tempoutput12);
NumericMatrix tempoutput13=addremoveinfection(2);
loglik1=clone(tempoutput13);
NumericMatrix tempoutput14=addremoveinfection(3);
loglik2=clone(tempoutput14);
NumericMatrix tempoutput15=addremoveinfection(4);
loglik3=clone(tempoutput15);
}

List addremoveinfection=add_remove_infection(data11,data111,data21,ILI,p_para(b0,_),p_para2(b0,_),loglik1,loglik2,loglik3,hai_start);
NumericMatrix tempoutput11=addremoveinfection(0);
data11=clone(tempoutput11);
NumericMatrix tempoutput12=addremoveinfection(1);
data21=clone(tempoutput12);
NumericMatrix tempoutput13=addremoveinfection(2);
loglik1=clone(tempoutput13);
NumericMatrix tempoutput14=addremoveinfection(3);
loglik2=clone(tempoutput14);
NumericMatrix tempoutput15=addremoveinfection(4);
loglik3=clone(tempoutput15);

updateacceptrate(b0,1)=addremoveinfection(5);


// update the impute likelihood
LL4(b0,1)=sum(loglik1);
LL4(b0,2)=sum(loglik2);
LL4(b0,3)=sum(loglik3);
LL4(b0,0)=LL4(b0,1)+LL4(b0,2)+LL4(b0,3);

temploglik(0)=LL4(b0,0)+prior_loglik(p_para(b0,_),n_seasons);
temploglik(1)=LL4(b0,1);
temploglik(2)=LL4(b0,2);
temploglik(3)=LL4(b0,3);


if (b0>=burnin){
if (b0%(thinning)==thinning-1){	
for (b1=data11.nrow()-1;b1>=0;--b1){
// record infection status	
impute_record_infstatus1(recordcurrentindex,b1)=data21(b1,3);
// record infection time	
impute_record_inftime1(recordcurrentindex,b1)=data21(b1,4);
// here record the waning rate
impute_record_waning1(recordcurrentindex,b1)=data21(b1,6);
// here record the boosting
impute_record_boosting1(recordcurrentindex,b1)=data21(b1,5);
// here record the baseline AT
impute_record_titer1(recordcurrentindex,b1)=data11(b1,8);
}
++recordcurrentindex;
}
}

totalinfectionnumber(b0,0)=sum(data21(_,3));


// move the matirx to another matrix to store the parameter and compute the correlation matrix
moveindex=sum(move)-1;
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){
p_para_r(b0,moveindex)=p_para(b0,b1);
--moveindex;
}	
}


if (b0%1000==0){
Rcout << b0 << std::endl;
}

}


return List::create(_[""]=p_para,
_[""]=p_para2,
_[""]=p_para3,
_[""]=LL1,
_[""]=LL2,
_[""]=LL3,
_[""]=LL4,
_[""]=updateacceptrate,
_[""]=totalinfectionnumber,
_[""]=data11,
_[""]=data111,
_[""]=data21,
_[""]=impute_record_infstatus1,
_[""]=impute_record_inftime1,
_[""]=impute_record_waning1,
_[""]=impute_record_boosting1,
_[""]=impute_record_titer1);
} 

/*
_[""]=impute_record_boosting21,
_[""]=impute_record_boosting22,
_[""]=impute_record_waning1,
_[""]=impute_record_waning2,
_[""]=impute_record_titer1,
_[""]=impute_record_titer2
*/







