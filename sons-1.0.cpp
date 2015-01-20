#include <iostream>

using namespace std;

#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <ctime>
#include <cmath>
#include <map>

/**************************************************************************************************/

struct Similarity {
	int x, y, count;
	float J, Jhat_se, L, Lhat_se, Uhat, Uhat_se, Vhat, Vhat_se, Uobs, Vobs, Aotu_shared, Botu_shared, UVhat, UVhat_se, V, V_se, Chao, jaccard, sorenson, thetaN, thetaYC, thetaYC_se;
};

/**************************************************************************************************/

struct Accession {
	int order;
	int group;
	int otu;
};

/**************************************************************************************************/

void zero_out(Similarity& info)
{	
	info.J = 0;				info.Jhat_se = 0;
	info.L = 0;				info.Lhat_se = 0;
	info.Uhat = 0;			info.Uhat_se = 0;
	info.Vhat = 0;			info.Vhat_se = 0;
	info.UVhat = 0;			info.UVhat_se = 0;
	info.V = 0;				info.V_se = 0;
	info.Chao = 0;
	info.Uobs = 0;			info.Vobs = 0;
	info.Aotu_shared = 0;	info.Botu_shared = 0;
	info.jaccard = 0;		info.sorenson = 0;
	info.thetaN = 0;		info.thetaYC = 0;		info.thetaYC_se = 0;
}

/**************************************************************************************************/

template<typename T>
string tostring(const T&x){
	ostringstream out;
	out << x;
	return out.str();
}

/**************************************************************************************************/

void usageError(char *name)
{
	cerr << "Usage: " << name << "" << endl;
	cerr << "Options:\n";
	exit(1);
}

/**************************************************************************************************/

class Sims {
public:
	Sims(vector<vector<int> >, int, int, map<int,int>, vector<vector<int> >&, string, vector<string>, ostream&);
private:
//	Functions and variables for calculating similarity values	
	void get_fs(vector<int>, vector<int>);
	void calc_UV(vector<int>, vector<int>);
	void calc_thetaYC(vector<int>, vector<int>, float&, float&);
	void calc_thetaN(vector<int>, vector<int>, float&);
	float calc_Jhat(void);
	float calc_Lhat(void);
	void JL_bootstrap(vector<int>, vector<int>, int, float&, float&, float&, float&, float&);
	int f1p, f2p, fp1, fp2, f11;
	float Uhat, Vhat;

//	Functions and variables for calculating shared richness
	float calc_Chao(void);
	void get_rarefs(vector<int>, vector<int>, vector<int>&, vector<int>&);
	float calc_c12(vector<int>, vector<int>);
	void calc_T(vector<int>, vector<int>, vector<vector<int> >&);
	void calc_G(vector<vector<int> >, float, vector<float>&);
	float calc_V(vector<float>, float);
	void V_bootstrap(vector<int>, vector<int>, float, int, float&);
	int rf1p, rfp1, rf11;
	int Vabund, Vrare;
	int n, m;

//	General purpose functions and variables
	void read(string);
	void get_nonzero(int, int, vector<int>&, vector<int>&);
	void get_XYdata(vector<int>, vector<int>, vector<int>&, vector<int>&);
	int assemblages;
	int tot_species;
	int Xindiv, Xspec;
	int Xshared_sum, Yshared_sum;
	int Yindiv, Yspec;
	int shared;
	int NXsing, NYsing;

	vector<vector<int> > list_data;
};

/**************************************************************************************************/
	
Sims::Sims(vector<vector<int> > input_data, int iterations, int group, map<int,int> seqs_in_group, 
					vector<vector<int> >& indices, string distance, vector<string> names, ostream& f)
{
	list_data = input_data;
	tot_species = list_data.size();
	assemblages = list_data[0].size();
	
	Similarity data;
	zero_out(data);
	int index = 0;
	
	
	for(int i=0;i<assemblages;i++){
		
		if(i != group){
			string key, key_names;
		
			int i_gt_group = 0;
			
			if(i<group)	{
				key = tostring(i) + "\t" + tostring(group) + "\t" + tostring(indices[i][group]++);
				key_names = names[i] + "\t" + names[group];
				
			}
			else		{
				key = tostring(group) + "\t" + tostring(i) + "\t" + tostring(indices[group][i]++);	
				key_names = names[group] + "\t" + names[i];
				i_gt_group = 1;
			}
	
			vector<int> X;			vector<int> Y;
			vector<int> Xshared;	vector<int> Yshared;
			vector<int> Xrare;		vector<int> Yrare;

			vector<vector<int> > T;
			vector<float> G(3,0.0);
			
			get_nonzero(i, group, X, Y);
			get_XYdata(X, Y, Xshared, Yshared);
			
					int flag = 0;
			if(Xindiv + Yindiv == seqs_in_group[i] + seqs_in_group[group]){	flag = 1;}
				
			if(Xindiv > 0){			data.Uobs = (float)Xshared_sum / (float)Xindiv;			}
			if(Yindiv > 0){			data.Vobs = (float)Yshared_sum / (float)Yindiv;			}

			float thetaN=0, thetaYC=0, thetaYC_se=0;
			
			if(Xindiv > 1 && Yindiv > 1){
				calc_thetaYC(X, Y, thetaYC, thetaYC_se);
			}
			
			
//			cout << i << '\t' << group << '\t' << Xindiv << '\t' << Yindiv << endl;
			
			if(Xspec != 0){			data.Aotu_shared = (float)shared / (float)Xspec;		}
			if(Yspec != 0){			data.Botu_shared = (float)shared / (float)Yspec;		}
			
			data.jaccard = (float)shared / (float)(Xspec + Yspec - shared);
			data.sorenson = 2.0*(float)shared / (float)(Xspec + Yspec);	

			float J=0, L=0, Jhat_se=0, Lhat_se=0, Uhat_se=0, Vhat_se=0, UVhat_se=0, V=0, V_se=0, Chao=0;
			fp1=0,fp2=0,f1p=0,f2p=0;
			
			int NXsing_save, NYsing_save;
			
			if(shared != 0){
				
				thetaN = data.Uobs * data.Vobs / (data.Uobs+data.Vobs-data.Uobs*data.Vobs);	

				get_fs(Xshared, Yshared);
				calc_UV(Xshared, Yshared);	
				
				J = calc_Jhat();				L = calc_Lhat();				Chao = calc_Chao();
				
				data.Uhat = Uhat;				data.Vhat = Vhat;				data.UVhat = Uhat * Vhat;

				NXsing_save = NXsing;			NYsing_save = NYsing;
				NXsing = NXsing_save;			NYsing = NYsing_save;

				get_rarefs(Xshared, Yshared, Xrare, Yrare);
				float c12 = calc_c12(Xrare, Yrare);	

				if(c12>0){
					calc_T(Xrare, Yrare, T);		
					calc_G(T, c12, G);					
					V = calc_V(G,c12);
				}
				
				if(flag == 1){
					JL_bootstrap(X, Y, iterations, Lhat_se, Jhat_se, Uhat_se, Vhat_se, UVhat_se);
					if(V > 1){
						V_bootstrap(Xshared, Yshared, V, iterations, V_se);
					}
				}
			}
			else{
				J = 0;				Jhat_se = 0;
				L = 0;				Lhat_se = 0;
				data.Uhat = 0;		data.Uhat_se = 0;
				data.Vhat = 0;		data.Vhat_se = 0;
				data.UVhat = 0;		data.UVhat_se = 0;
				V = 0;				V_se = 0;				Chao = 0;
			}
			
			data.J = J;						data.Jhat_se = Jhat_se;
			data.L = L;						data.Lhat_se = Lhat_se;
			
			data.Uhat_se = Uhat_se;			data.Vhat_se = Vhat_se;			data.UVhat_se = UVhat_se;
			data.V = V;						data.V_se = V_se;				data.Chao = Chao;
			data.thetaYC = thetaYC;			data.thetaYC_se = thetaYC_se;	data.thetaN = thetaN;
			
			if(flag == 1){
				cout << distance << "\t" << key_names << "\t" << setprecision(4);

				if(i_gt_group == 0){
					cout << data.Uhat << "\t" << data.Uhat_se << "\t";
					cout << data.Vhat << "\t" << data.Vhat_se <<"\t";
				}
				else{
					cout << data.Vhat << "\t" << data.Vhat_se <<"\t";
					cout << data.Uhat << "\t" << data.Uhat_se << "\t";
				}
				
				cout << data.UVhat << "\t" << data.UVhat_se <<"\t";
				cout << data.J << "\t" << data.Jhat_se << "\t";
				cout << data.L << "\t" << data.Lhat_se << "\t";
				cout << setprecision(1) << data.V << "\t" << data.V_se << "\t";
				cout << data.Chao << "\t";
			
				if(i_gt_group == 0)	{	cout << setprecision(4) << data.Uobs << "\t" << data.Vobs << "\t";}
				else				{	cout << setprecision(4) << data.Vobs << "\t" << data.Uobs << "\t";}

				if(i_gt_group == 0)	{	cout << setprecision(4) << data.Aotu_shared << "\t" << data.Botu_shared << "\t";}
				else				{	cout << setprecision(4) << data.Botu_shared << "\t" << data.Aotu_shared << "\t";}
				
				cout << data.jaccard << "\t" << data.sorenson << '\t' << thetaYC << '\t' << thetaYC_se << '\t' << thetaN << endl;
			}

			f << distance << '\t' << key << "\t" << setprecision(4);	
		
			if(i_gt_group == 0){
				f << data.Uhat << "\t" << data.Uhat_se << "\t";
				f << data.Vhat << "\t" << data.Vhat_se <<"\t";
			}
			else{
				f << data.Vhat << "\t" << data.Vhat_se <<"\t";
				f << data.Uhat << "\t" << data.Uhat_se << "\t";
			}
			f << data.UVhat << "\t" << data.UVhat_se <<"\t";
			f << data.J << "\t" << data.Jhat_se << "\t";
			f << data.L << "\t" << data.Lhat_se << "\t";
			f << data.V << "\t" << data.V_se << "\t";
			f << data.Chao << "\t";

			if(i_gt_group == 0)	{	f << data.Uobs << "\t" << data.Vobs << "\t";}
			else				{	f << data.Vobs << "\t" << data.Uobs << "\t";}

			if(i_gt_group == 0)	{	f << data.Aotu_shared << "\t" << data.Botu_shared << "\t";}
			else				{	f << data.Botu_shared << "\t" << data.Aotu_shared << "\t";}
			
			f << data.jaccard << "\t" << data.sorenson << '\t' << data.thetaYC << '\t' << data.thetaYC_se << '\t' << data.thetaN << endl;
		}
	}
}

/**************************************************************************************************/

void Sims::get_nonzero(int a1, int a2, vector<int>& X, vector<int>& Y)
{
	for(int i=0;i<tot_species;i++){
		if(list_data[i][a1] != 0 || list_data[i][a2] != 0){
			X.push_back(list_data[i][a1]);
			Y.push_back(list_data[i][a2]);
		}
	}
}

/**************************************************************************************************/

void Sims::get_XYdata(vector<int> X, vector<int> Y, vector<int>& Xshared, vector<int>& Yshared)
{
	Xindiv=0, Yindiv=0;
	Xspec=0, Yspec=0;
	NXsing=0, NYsing=0;
	
	Xshared_sum = 0;
	Yshared_sum = 0;
	
	for(int i=0;i<X.size();i++){
		Xindiv += X[i];
		Yindiv += Y[i];

		if(X[i] <= 10)		{	NXsing += X[i];	}
		if(Y[i] <= 10)		{	NYsing += Y[i];	}

		if(X[i] != 0)		{	Xspec++;		}
		if(Y[i] != 0)		{	Yspec++;		}

		if(X[i] != 0 && Y[i] != 0){
			Xshared.push_back(X[i]);
			Yshared.push_back(Y[i]);
			
			Xshared_sum += X[i];
			Yshared_sum += Y[i];

		}
	}

	//cout << Xindiv << '\t' << Yindiv << '\t' << shared << endl;
			//cout << Xshared_sum << '\t' << Yshared_sum << endl;

	shared = Xshared.size();
}

/**************************************************************************************************/

void Sims::calc_thetaYC(vector<int> X, vector<int> Y, float& thetaYC, float& thetaYC_se){

	float a = 0, b = 0, d = 0, pminq_sq = 0, q_sq = 0;
	float pi3 = 0, qi3 = 0, pq2 = 0, p2q = 0;
	
//	float pi, qi;
		
	for(int z=0;z<X.size();z++){
		float pi = (float)X[z]/(float)Xindiv;
		float qi = (float)Y[z]/(float)Yindiv;
//		float pi_1 = (float)(X[z]-1)/(float)(Xindiv-1);  //bias-corrected theta gives values that are negative and >1 
//		float qi_1 = (float)(Y[z]-1)/(float)(Yindiv-1);  //bias-corrected theta gives values that are negative and >1

		a += pi*pi;	
		b += qi*qi;
		d += pi*qi;

		pi3 += pow(pi,3);
		qi3 += pow(qi,3);
		pq2 += pi * pow(qi,2);
		p2q += qi * pow(pi,2);

	}
	
	double den = a+b-d;
	
	if(den > 0){
		
		thetaYC = d / den;
		
		float varA = (4.0/(float)Xindiv) * (pi3 - a*a);
		float varB = (4.0/(float)Yindiv) * (qi3 - b*b);
		float varD =  pq2/(float)Xindiv + p2q/(float)Yindiv - (1/(float)Xindiv+1/(float)Yindiv)*pow(d,2);
		float covAD = (2.0/(float)Xindiv) * (p2q - a * d);
		float covBD = (2.0/(float)Yindiv) * (pq2 - b * d);	
		
		double thetaYC_var = (pow(d,2)*(varA + varB) + pow(a+b,2)*varD-2*(a+b)*d*(covAD+covBD))/pow(den,4);

		if(thetaYC_var>0)	{	thetaYC_se = pow(thetaYC_var, 0.5);		}
		else				{	thetaYC_se = 0;							}
	}
	else{
		thetaYC = 0;
		thetaYC_se = 0;
	}

//	cout << setprecision(4);
//	cout << a << '\t' << b << '\t' << d << '\t';
//	cout << pi3 << '\t' << qi3 << '\t' << pq2 << '\t' << p2q << '\t';
//	cout << varA << '\t' << varB << '\t' << varD << '\t' << covAD << '\t' << covBD << '\t';
//	cout << pq2/(float)Xindiv << '\t' << p2q/(float)Yindiv << '\t' << (1./(float)Xindiv+1./(float)Yindiv)*pow(d,2) << '\t';
//	cout << setprecision(6) << thetaYC << '\t' << '\t' << var << endl;
}

/**************************************************************************************************/

void Sims::calc_thetaN(vector<int> Xshared, vector<int> Yshared, float& thetaN){

	float pi, qi, pi_sum, qi_sum;
	
	for(int z=0;z<shared;z++){
		pi = (float)Xshared[z]/(float)Xindiv;
		qi = (float)Yshared[z]/(float)Yindiv;
		pi_sum += pi;
		qi_sum += qi;
	}
	
	thetaN = pi_sum * qi_sum / (pi_sum + qi_sum - pi_sum * qi_sum);
	cout << pi_sum << '\t' << qi_sum << '\t' << thetaN << endl;

//	cout << setprecision(4);
//	cout << a << '\t' << b << '\t' << d << '\t';
//	cout << pi3 << '\t' << qi3 << '\t' << pq2 << '\t' << p2q << '\t';
//	cout << varA << '\t' << varB << '\t' << varD << '\t' << covAD << '\t' << covBD << '\t';
//	cout << pq2/(float)Xindiv << '\t' << p2q/(float)Yindiv << '\t' << (1./(float)Xindiv+1./(float)Yindiv)*pow(d,2) << '\t';
//	cout << setprecision(6) << thetaYC << '\t' << '\t' << var << endl;
}

/**************************************************************************************************/

void Sims::get_fs(vector<int> Xshared, vector<int> Yshared)
{
	f1p = 0, f2p = 0;
	fp1 = 0, fp2 = 0;
	f11 = 0;

	for(int i=0;i<shared;i++){
		if(Xshared[i] == 1)		{	f1p++;		}
		else if(Xshared[i] == 2){	f2p++;		}

		if(Yshared[i] == 1)		{	fp1++;		}
		else if(Yshared[i] == 2){	fp2++;		}

		if(Xshared[i] == 1 && Yshared[i] == 1){	f11++;	}
	}
	
	if(fp2 == 0){	fp2 = 1;	}
	if(f2p == 0){	f2p = 1;	}
}

/**************************************************************************************************/

void Sims::get_rarefs(vector<int> Xshared, vector<int> Yshared, vector<int>& Xrare, vector<int>& Yrare)
{
	rf1p = 0, rfp1 = 0, rf11 = 0, Vabund = 0, Vrare = 0;
	int N1p = 0, Np1 = 0, rare_n = 0, rare_m = 0;

	int l_shared = Xshared.size();
	
	for(int i=0;i<l_shared;i++){
		if(Xshared[i] > 10 || Yshared[i] > 10){
			Vabund++;
		}
		else{
			Vrare++;
			Xrare.push_back(Xshared[i]);
			Yrare.push_back(Yshared[i]);
			
			rare_n += Xshared[i];
			rare_m += Yshared[i];
			
			if(Xshared[i] == 1)						{	rf1p++;		}
			if(Yshared[i] == 1)						{	rfp1++;		}
			if(Xshared[i] == 1 && Yshared[i] == 1)	{	rf11++;		}
		}
	}

	n = NXsing;
	m = NYsing;
	
}

/**************************************************************************************************/

void Sims::calc_UV(vector<int> Xshared, vector<int> Yshared)
{
	float U_first_term = 0.0;
	float V_first_term = 0.0;

	float U_sec_term = 0.0;
	float V_sec_term = 0.0;
	
	for(int i=0;i<shared;i++){
		U_first_term += (float)Xshared[i];
		V_first_term += (float)Yshared[i];

		if(Yshared[i] == 1){	U_sec_term += (float)Xshared[i];	}
		if(Xshared[i] == 1){	V_sec_term += (float)Yshared[i];	}
	}
//	cout << Xindiv << '\t' << Yindiv << endl;
	U_first_term /= (float)Xindiv;
	V_first_term /= (float)Yindiv;

	U_sec_term = U_sec_term * ((float)(Yindiv-1)/(float)Yindiv) * (0.5*(float)fp1/(float)fp2) / (float)Xindiv;
	V_sec_term = V_sec_term * ((float)(Xindiv-1)/(float)Xindiv) * (0.5*(float)f1p/(float)f2p) / (float)Yindiv;

	Uhat = U_first_term + U_sec_term;
	Vhat = V_first_term + V_sec_term;
	
	if(Uhat > 1){	Uhat = 1;	}
	if(Vhat > 1){	Vhat = 1;	}
}

/**************************************************************************************************/

float Sims::calc_Chao(void)
{

	if(fp2 == 0 || f2p == 0){
		return (float)shared + (float)f11*(float)(f1p*fp1)/float(4*(f2p+1)*(fp2+1))
								+ float(f1p*(f1p-1))/float(2*f2p+2)
								+ float(fp1*(fp1-1))/float(2*fp2+2);
	}
	else{
		return (float)shared + (float)f11*(float)(f1p*fp1)/float(4*f2p*fp2)
								+ float(f1p*f1p)/float(2*f2p)
								+ float(fp1*fp1)/float(2*fp2);
	}
		
}

/**************************************************************************************************/

float Sims::calc_Jhat(void)
{
	return Uhat * Vhat / (Uhat + Vhat - Uhat * Vhat);	
}

/**************************************************************************************************/

float Sims::calc_Lhat(void)
{
	return 2.0 * Uhat * Vhat / (Uhat + Vhat);	
}

/**************************************************************************************************/

void Sims::JL_bootstrap(vector<int> X, vector<int> Y, int iter, float& Lhat_se, float& Jhat_se, float& Uhat_se, float& Vhat_se, float& UVhat_se)
{
	int D = X.size();
	
	float Jhatsq = 0.0;
	float Lhatsq = 0.0;
	float Uhatsq = 0.0;
	float Vhatsq = 0.0;
	float UVhatsq = 0.0;
	
	float Jhatsum = 0.0;
	float Lhatsum = 0.0;
	float Uhatsum = 0.0;
	float Vhatsum = 0.0;
	float UVhatsum = 0.0;
	
	for(int i=0;i<iter;i++){
		vector<int> Xrand(D, 0);		vector<int> Yrand(D, 0);
		vector<int> Xshared;			vector<int> Yshared;
	
		for(int j=0;j<D;j++){
			int z = int((float)(D) * (float)(rand()) / ((float)RAND_MAX+1.0));
			Xrand[j] = X[z];
			Yrand[j] = Y[z];
		}

		get_XYdata(Xrand, Yrand, Xshared, Yshared);
	
		float Jhat, Lhat, UVhat;
		if(shared != 0){
			get_fs(Xshared, Yshared);
			calc_UV(Xshared, Yshared);
			Jhat = calc_Jhat();
			Lhat = calc_Lhat();
			UVhat = Uhat * Vhat;
		}
		else{
			Jhat = 0;
			Lhat = 0;
			UVhat = 0;
			Uhat = 0;
			Vhat = 0;
		}
		Jhatsq += Jhat*Jhat;
		Lhatsq += Lhat*Lhat;
		Uhatsq += Uhat*Uhat;
		Vhatsq += Vhat*Vhat;
		UVhatsq += UVhat*UVhat;
		
		Jhatsum += Jhat;
		Lhatsum += Lhat;
		Uhatsum += Uhat;
		Vhatsum += Vhat;
		UVhatsum += UVhat;
	}

	Lhat_se = (float)pow(double(Lhatsq - (Lhatsum*Lhatsum)/(float)iter)/(float)(iter - 1),0.5);
	Jhat_se = (float)pow(double(Jhatsq - (Jhatsum*Jhatsum)/(float)iter)/(float)(iter - 1),0.5);
	
	Uhat_se = (float)pow(double(Uhatsq - (Uhatsum*Uhatsum)/(float)iter)/(float)(iter - 1),0.5);
	Vhat_se = (float)pow(double(Vhatsq - (Vhatsum*Vhatsum)/(float)iter)/(float)(iter - 1),0.5);
	UVhat_se = (float)pow(double(UVhatsq - (UVhatsum*UVhatsum)/(float)iter)/(float)(iter - 1),0.5);

	//if(tostring(Lhat_se) == "nan")	{	Lhat_se = 0;	}
	//if(tostring(Jhat_se) == "nan")	{	Jhat_se = 0;	}
	//if(tostring(Uhat_se) == "nan")	{	Uhat_se = 0;	}
	//if(tostring(Vhat_se) == "nan")	{	Vhat_se = 0;	}
	//if(tostring(UVhat_se) == "nan")	{	UVhat_se = 0;	}

}

/**************************************************************************************************/

float Sims::calc_c12(vector<int> Xrare, vector<int> Yrare){

	int numerator = 0;
	int denomenator = 0;
	int num_rare = Xrare.size();
	
	for(int i=0;i<num_rare;i++){
		if(Xrare[i] == 1)					{	numerator += Yrare[i];	}
		if(Yrare[i] == 1)					{	numerator += Xrare[i];	}
		if(Xrare[i] == 1 && Yrare[i] == 1)	{	numerator--;			}
		
		denomenator += (Xrare[i] * Yrare[i]);
	}
	
	return (1-(float)numerator/(float)denomenator);
}

/**************************************************************************************************/

void Sims::calc_T(vector<int> Xrare, vector<int> Yrare, vector<vector<int> >& T)
{
	T.resize(3);
	T[0].resize(3,0);
	T[1].resize(3,0);
	T[2].resize(3,0);
	
	for(int i=0;i<Vrare;i++){
		T[1][0] += Xrare[i];
		T[0][1] += Yrare[i];
		T[1][1] += (Xrare[i] * Yrare[i]);
		T[2][1] += (Xrare[i] * (Xrare[i] - 1) * Yrare[i]);
		T[1][2] += (Xrare[i] * (Yrare[i] - 1) * Yrare[i]);
		T[2][2] += (Xrare[i] * Yrare[i] * (Xrare[i] - 1) * (Yrare[i] - 1));
	}
}

/**************************************************************************************************/

void Sims::calc_G(vector<vector<int> > T, float C, vector<float>& G)
{
	float NO = (float)Vrare/C;
	
	G[0] = NO * float(n)/float(n-1) * (float(T[2][1])/float(T[1][0]))/float(T[1][1]) - 1.0;
	
	G[1] = float(m * NO * T[1][2]) / float((m-1)*T[0][1]*T[1][1]) - 1.0;
	
	G[2] = (NO * NO) * (((float(T[2][2])/(float(T[0][1]))/float(T[1][0]))/float(T[1][1])));
	G[2] = G[2] * ( float(n)/float(n-1) ) * ( float(m)/float(m-1) );
	G[2] = G[2] - float(NO*T[1][1])/float(T[0][1]*T[1][0]);
	G[2] = G[2] - G[1] - G[0];

	if(n == 1 || m == 1 || T[0][1] == 0 || T[1][0] == 0 || T[1][1]){
		G[0] = 0;
		G[1] = 0;
		G[2] = 0;
	}
}

/**************************************************************************************************/

float Sims::calc_V(vector<float> G, float C)
{
	float V = Vabund + (Vrare + rf1p*G[0] + rfp1*G[1] + rf11*G[2]) / C;
	
//	if(tostring(V) == "nan"){//	{	return V;	}
//		cout << endl << G[0] << '\t' << G[1] << '\t' << G[2] << '\t' << C << endl << endl;
//	}

	if (V>0)	{	return V;	}
	else		{	return 0;	}
}

/**************************************************************************************************/

void Sims::V_bootstrap(vector<int> Xshared, vector<int> Yshared, float V_saved, int iter, float& V_se)
{
	int N = int(V_saved);
	int M = Xshared.size();
		
	float V_sq = 0.0;
	float V_sum = 0.0;
//	vector<float> V(iter, 0.0);

	for(int i=0;i<iter;i++){
		vector<int> XsharedR(N,0);		vector<int> YsharedR(N,0);
		vector<int> Xrare;				vector<int> Yrare;
		
		vector<vector<int> > T;			vector<float> G(3,0);
	
		for(int j=0;j<N;j++){
			int z = int((float)(M) * (float)(rand()) / ((float)RAND_MAX+1.0));
			XsharedR[j] = Xshared[z];
			YsharedR[j] = Yshared[z];
		}

		get_rarefs(XsharedR, YsharedR, Xrare, Yrare);	

		float c12 = calc_c12(Xrare, Yrare);
		if( c12 > 0 ){
			calc_T(Xrare, Yrare, T);
			calc_G(T, c12, G);
//			V[i] = calc_V(G, c12);
			float V = calc_V(G, c12);
			V_sq += V * V;
			V_sum += V;
		}
		else{	i--;		}
	}
//	sort(V.begin(),V.end());
//	V_se = V[int(iter*0.025)];
	
	V_se = (float)pow(double(V_sq - (V_sum*V_sum)/(float)iter)/(float)(iter - 1),0.5);
}

/**************************************************************************************************/

class GetData {
	public:
		GetData(string, string, int, int);
	private:
		void get_order_groups(string, int);
		void get_otus(string);
		void collect(vector<Accession>, int, ostream&);
		void format_sim_output(string);
		void format_otu_output(string);

		string current_distance;
		int num_groups;
		int num_sequences;
		int iterations;
		int row_count;
		vector<int> row_length;
		vector<string> distances;
		vector<string> group_id;
		map<int,int> seqs_per_group;
		
		map<string,Accession> accession_data;
};

/**************************************************************************************************/

GetData::GetData(string list_file, string name_file, int iter, int jumble)
{
	map<string,int> names_ids;
	iterations = iter;
	row_count = 0;
	
	get_order_groups(name_file, jumble);
	get_otus(list_file);
	format_sim_output(list_file);
	format_otu_output(list_file);
}


/**************************************************************************************************/

void GetData::get_order_groups(string fname, int jumble)
{
	ifstream f(fname.c_str());
	if(!f) {
		cerr << "Error: Could not open " << fname << endl;
		exit(1);
	}
	
	map<string,int> group_lookup;
	
	int nindex = 0, gindex = 0;
	
	while(f.peek() != EOF){
		string name, group;
		f >> name >> group;								//	read in name and group id

		if(name == "" && group == ""){break;}
		f.get();	

		accession_data[name].order = nindex;			//	give each sequence name an index
		nindex++;										//	increment index for next go around
		
		if(group_lookup.count(group) == 0){				//  if this group has never been seen before create it
			group_lookup[group] = gindex;				//	assign a number to each group
			group_id.push_back(group);
			gindex++;									//	increment group index for next go around
		}

		seqs_per_group[group_lookup[group]]++;			//	increment the number of sequences in that group

		accession_data[name].group = group_lookup[group];//	assign group number to accession
		
		num_groups  = group_lookup.size();				//	determine the number of groups
		num_sequences = accession_data.size();			//	determine the number of sequences
	}

	if(jumble == 1){
		vector<int> order(num_sequences);
		for(int i=0;i<num_sequences;i++){			order[i] = i;			}

		for(int j=num_sequences-1;j>=0;j--){
			int z = int((float)(j+1) * (float)(rand()) / ((float)RAND_MAX+1.0));	
			int t = order[z]; order[z]=order[j];order[j]=t;
		}
	
		map<string,Accession>::iterator pos;
		for(pos=accession_data.begin();pos!=accession_data.end();pos++){
			pos->second.order = order[pos->second.order];
		
		}
	}
}

/**************************************************************************************************/

void GetData::get_otus(string fname)
{
	ifstream f(fname.c_str());
	if(!f) {
		cerr << "Error: Could not open " << fname << endl;
		exit(1);
	}

	string c_outname = "collect_data";
	ofstream coll(c_outname.c_str(), ios::trunc);
	

	cout << "Dist\tA\tB\tUhat\tUhat_se\tVhat\tVhat_se\tUVhat\tUV_se\tJhat\tJhat_se\tLhat\tLhat_se\tVest\tVest_se\tChao\tUobs\tVobs\tAotu_shared\tBotu_shared\tJclas\tLclas\tthetaYC\tthetaYC_se\tthetaN\n";

	while(f.peek() != EOF){
		vector<Accession> order_group_otu(num_sequences);
		
		f >> current_distance;
		distances.push_back(current_distance);
		
		int otus;			f >> otus;
		char space;			space = f.get();
		int index = 0;
	
		char d = f.get();
		
		while(d != '\n'){
			string name;
			if(d != '\t' && d != ',' && d != '\n'){
				
				while(d != ',' && d!='\t' && d!='\n'){						
					name += d;
					d = f.get();
				}
				
				if(accession_data.count(name) == 0){
					cout << name << " not found, check *.names file\n";
					exit(1);
				}
				order_group_otu[accession_data[name].order].group = accession_data[name].group;
				order_group_otu[accession_data[name].order].otu = index;
			}
			else if(d == '\t'){
				d = f.get();
				index++;
			}
			else if(d == ','){
				d = f.get();
			}	
		}
		collect(order_group_otu, otus, coll);
	}
}

/**************************************************************************************************/

void GetData::collect(vector<Accession> order_data, int num_otus, ostream& f)
{
	vector<vector<int> > list_data(num_otus);
	for(int i=0;i<num_otus;i++){		list_data[i].resize(num_groups, 0);			}
	
	vector<vector<int> > indices(num_groups);
	for(int i=0;i<num_groups;i++){		indices[i].resize(num_groups,0);			}

	for(int i=0;i<num_sequences;i++){
		list_data[order_data[i].otu][order_data[i].group]++;
		Sims coll(list_data, iterations, order_data[i].group, seqs_per_group, indices, current_distance, group_id, f);
	}

	string count_file = "count_data";
	ofstream otu_counts(count_file.c_str(), ios::app);
	
	for(int i=0;i<num_groups;i++){
		otu_counts << current_distance << '\t' << group_id[i] << '\t' << num_otus << '\t';
	
		for(int j=0;j<num_otus;j++){
			otu_counts << list_data[j][i] << '\t';
		}
		otu_counts << endl;
	
		row_count++;
		row_length.push_back(num_otus+3);
	}
	f << ";\n";
}

/**************************************************************************************************/

void GetData::format_sim_output(string list_file)
{
	cout << "Cleaning up...\n";
	string input = "collect_data";
	ifstream f(input.c_str());
	if(!f) {
		cerr << "Error: Could not open " << input << endl;
		exit(1);
	}

	if(list_file.find_last_of(".")!=string::npos){
		list_file.erase(list_file.find_last_of(".")+1);
	}
	else{
		list_file += ".";
	}
	string output_collect = list_file + "sons";
	string output_ltt = output_collect + ".ltt";
	
	ofstream sort_data(output_collect.c_str(), ios::trunc);
	sort_data.setf(ios::fixed, ios::floatfield);
	sort_data.setf(ios::showpoint);

	ofstream ltt (output_ltt.c_str(), ios::trunc);
	ltt.setf(ios::fixed, ios::floatfield);
	ltt.setf(ios::showpoint);

	if(!sort_data || !ltt) {
		cerr << "Error: Could not open " << output_collect << endl;
		exit(1);
	}

	sort_data << "Distance\tA\tB\tSSamp\tUhat\tVhat\tUVhat\tJhat\tLhat\tVest\tSharedChao1\tUobs\tVobs\tAotu_shared\tBotu_shared\tJclass\tSclass\tthetaYC\tthetaYC_se\tthetaN\n";
	ltt << "Distance\tA\tB\tUhat\tUhat_se\tVhat\tVhat_se\tUVhat\tUVhat_se\tJhat\tJhat_se\tLhat\tLhat_se\tVest\tVest_se\tSharedChao1\tUobs\tVobs\tAotu_shared\tBotu_shared\tJclass\tSclass\tthetaYC\tthetaYC_se\tthetaN\n";

	string distance;
	f >> distance;
	
	while(f.peek() != EOF){
		string prev_distance = distance;

		vector<vector<vector<Similarity> > > output(num_groups);
		for(int i=0;i<num_groups;i++){
			output[i].resize(num_groups);
		}

		while(distance != ";"){
			Similarity line_data;
			f >> line_data.x		>> line_data.y			>> line_data.count;
			f >> line_data.Uhat		>> line_data.Uhat_se;
			f >> line_data.Vhat		>> line_data.Vhat_se;
			f >> line_data.UVhat	>> line_data.UVhat_se;
			f >> line_data.J		>> line_data.Jhat_se;
			f >> line_data.L		>> line_data.Lhat_se;
			f >> line_data.V		>> line_data.V_se;
			f >> line_data.Chao;
			f >> line_data.Uobs		>> line_data.Vobs;
			f >> line_data.Aotu_shared		>> line_data.Botu_shared;
			f >> line_data.jaccard	>> line_data.sorenson;
			f >> line_data.thetaYC	>> line_data.thetaYC_se >> line_data.thetaN;
			
			output[line_data.x][line_data.y].push_back(line_data);

			prev_distance = distance;
			f >> distance;
		}
		
		f >> distance;
		
		for(int i=0;i<num_groups;i++){
			for(int j=i+1;j<num_groups;j++){
				int sequences = output[i][j].size();
				for(int z=0;z<sequences;z++){
					sort_data << prev_distance << '\t';
					sort_data << group_id[output[i][j][z].x] << '\t' << group_id[output[i][j][z].y] << '\t';
					sort_data << output[i][j][z].count+1 << '\t' << setw(7) << setprecision(5);
					sort_data << output[i][j][z].Uhat << '\t' << output[i][j][z].Vhat << '\t' << output[i][j][z].UVhat << '\t';
					sort_data << output[i][j][z].J << '\t' << output[i][j][z].L << '\t';
					sort_data << output[i][j][z].V << '\t' << output[i][j][z].Chao << '\t';
					sort_data << output[i][j][z].Uobs << '\t' << output[i][j][z].Vobs << '\t';
					sort_data << output[i][j][z].Aotu_shared << '\t' << output[i][j][z].Botu_shared << '\t';
					sort_data << output[i][j][z].jaccard << '\t' << output[i][j][z].sorenson << '\t';
					sort_data << output[i][j][z].thetaYC << '\t' << output[i][j][z].thetaYC_se << '\t';
					sort_data << output[i][j][z].thetaN << endl;
					sort_data.flush();
				}
				ltt << prev_distance << '\t';
				ltt << group_id[output[i][j][sequences-1].x] << '\t' << group_id[output[i][j][sequences-1].y] << '\t';
				ltt << setw(7) << setprecision(5);
				ltt << output[i][j][sequences-1].Uhat		<< '\t' << output[i][j][sequences-1].Uhat_se	<< '\t';
				ltt << output[i][j][sequences-1].Vhat		<< '\t' << output[i][j][sequences-1].Vhat_se	<< '\t';
				ltt << output[i][j][sequences-1].UVhat		<< '\t' << output[i][j][sequences-1].UVhat_se	<< '\t';
				ltt << output[i][j][sequences-1].J			<< '\t' << output[i][j][sequences-1].Jhat_se	<< '\t';
				ltt << output[i][j][sequences-1].L			<< '\t' << output[i][j][sequences-1].Lhat_se	<< '\t';
				ltt << output[i][j][sequences-1].V			<< '\t' << output[i][j][sequences-1].V_se		<< '\t';
				ltt << output[i][j][sequences-1].Chao		<< '\t';
				ltt << output[i][j][sequences-1].Uobs		<< '\t' << output[i][j][sequences-1].Vobs	<< '\t';
				ltt << output[i][j][sequences-1].Aotu_shared		<< '\t' << output[i][j][sequences-1].Botu_shared	<< '\t';
				ltt << output[i][j][sequences-1].jaccard	<< '\t' << output[i][j][sequences-1].sorenson	<< '\t';
				ltt << output[i][j][sequences-1].thetaYC	<< '\t' << output[i][j][sequences-1].thetaYC_se	<< '\t';
				ltt << output[i][j][sequences-1].thetaN	<< '\n';
				ltt.flush();
			}
		}
	}
	remove("collect_data");
}

/**************************************************************************************************/

void GetData::format_otu_output(string list_file)
{
	string input = "count_data";
	
	ifstream f(input.c_str());
	if(!f) {
		cerr << "Error: Could not open " << input << endl;
		exit(1);
	}

	vector<vector<string> > otu_data(row_count);
	for(int i=0;i<row_count;i++){
	
		otu_data[i].resize(row_length[0]);
		
		for(int j=0;j<row_length[i];j++){
			f >> otu_data[i][j];
		}
		if(i>0){
			for(int j=row_length[i];j<row_length[0];j++){
				otu_data[i][j] = " ";
			}
		}
	}
		
	if(list_file.find_last_of(".")!=string::npos){
		list_file.erase(list_file.find_last_of(".")+1);
	}
	else{
		list_file += ".";
	}
	
	string output_otu_data = list_file + "sons.otu";
	ofstream otu_data_file(output_otu_data.c_str(), ios::trunc);
	otu_data_file.setf(ios::fixed, ios::floatfield);
	otu_data_file.setf(ios::showpoint);

	if(!otu_data_file) {
		cerr << "Error: Could not open " << output_otu_data << endl;
		exit(1);
	}
	
	for(int j=0;j<otu_data[0].size();j++){
		for(int i=0;i<row_count;i++){
			if(otu_data[i][j]== " "){	break;	}
			otu_data_file << otu_data[i][j] << '\t';
		}
		otu_data_file << endl;
	}
	remove("count_data");
}

/**************************************************************************************************/

int main(int argc, char *argv[])
{
	srand( (unsigned)time( NULL ) );
	cout.setf(ios::fixed, ios::floatfield);
	cout.setf(ios::showpoint);
	cerr.setf(ios::fixed, ios::floatfield);
	cerr.setf(ios::showpoint);

	string filename, list_file, names_file;
	int iter = 1000;
	int jumble = 0;
	
	char **p;
	if(argc>1){	
		for(p=argv+1;p<argv+argc;p++){
			if(strcmp(*p, "-i")==0){
				if(++p>=argv+argc)
					usageError(argv[0]);
				istringstream f(*p);
				if(!(f >> iter))
					usageError(argv[0]);
			}
			if(strcmp(*p, "-list")==0){
				if(++p>=argv+argc)
					usageError(argv[0]);
				istringstream f(*p);
				if(!(f >> list_file))
					usageError(argv[0]);
			}
			if(strcmp(*p, "-names")==0){
				if(++p>=argv+argc)
					usageError(argv[0]);
				istringstream f(*p);
				if(!(f >> names_file))
					usageError(argv[0]);
			}
			if(strcmp(*p, "-jumble")==0){
				if(p>=argv+argc)
					usageError(argv[0]);
				jumble=1; //sample with replacement
			}
			else{
				istringstream f(*p);
				if(!(f >> filename))
					usageError(argv[0]);
			}
		}
	}
	else
		usageError(argv[0]);

	remove("collect_data");
	remove("count_data");

	GetData list_data(list_file, names_file, iter, jumble);
	
	return 0;
}

