#ifndef WAVES_H_INCLUDED
#define WAVES_H_INCLUDED

#include <cmath>

double squareWave(double amplitude, double frequency, double sampleTime){

	double tmp;
	double posInWave=std::modf(sampleTime * frequency,&tmp);
	
	double rval =posInWave >.5?-amplitude : amplitude;

	return rval;
	
}





double sineWave(double amplitude, double frequency, double sampleTime){

	const double pi2 =2.0* 3.14159265;
	double sinVal = sin( (pi2 *(frequency) )*sampleTime);
	
	//std::cout<<sinVal<<std::endl;

	return sinVal * amplitude;

}
 
double whistleWave(double amplitude, double frequency, double sampleTime){

	return sineWave (.8*amplitude, frequency, sampleTime) + 
		sineWave (.2*amplitude, 20.0*frequency, sampleTime); 

		
}

double sawtoothWave(double amplitude, double frequency, double sampleTime){
	 
	double tmp = (sampleTime * frequency);

	return amplitude * ((tmp - floor(tmp)) *2.0) - 1.0;
	 
}


double tanWave(double amplitude, double frequency, double sampleTime){
	const double pi= 3.14159265;

	double tanVal = tan(pi *frequency*sampleTime);
	
	return tanVal *amplitude;
}


  
void initTable(std::vector<double>& table, double frequency ){
	
	table.resize(frequency * 2);///\todo fractional frequency a problem

	//std::cout<<table.size()<<std::endl;

	//system("pause");
	for(int i =0; i < table.size();i+=2){
	try{
		
		if(i+1>=table.size()){
		
			throw 'x';
		}

		table[i]   =  double(rand()) / RAND_MAX;
		table[i+1] = -double(rand()) / RAND_MAX;
		//std::cout<<i+1<<std::endl;
	}catch(...){
	
		std::cout<<"i is "<<i<<'\n'<<"size is "<<table.size()<<std::endl;

	}
	}

	/*const uint passes=1;
	//boxFilter<std::vector<double>::iterator>(table.begin(), table.end());
	for(int i = 0; i < passes; i++){
	for(int i =0; i < table.size()-2; i++){
	
		table[i] = .5*(table[i] + table[i+2]);

	}

	}*/
}
 


void initTable(std::vector<double>& table){
	
	table.resize(1024);///\todo fractional frequency a problem

	//std::cout<<table.size()<<std::endl;

	//system("pause");
	for(int i =0; i < table.size();i++){
	 
		
	
		table[i]   =  ((double(rand()) / RAND_MAX) - .5)*2.0;
	
	  
	 
	}

	 
	  
}




///\todo seems to run at half pitch too
///\todo nonsquare periodic noise.  random, smooth noise, and weighted
double squareNoise(double amplitude, double frequency, double sampleTime){

	static std::vector<double> table;	

	if(!table.size()){
		initTable(table);
	}

	
	uint idx = sampleTime * frequency;
	idx<<=1;//times 2 because we're indexing a table of values, and the wave should change direction twice per iteration
	idx %= table.size();

	double randVal = table[idx];

	return amplitude*randVal;//6(val>0.0?-1.0:1.0);

}



#endif