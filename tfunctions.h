#ifndef TFUNCTIONS_H_INCLUDED
#define TFUNCTIONS_H_INCLUDED

#include <boost/bind.hpp>
#include <vector>
#include <cstdint>


typedef unsigned int uint;

typedef std::vector<std::uint16_t> AudioSample16;
typedef std::vector<std::uint8_t> AudioSample8;

#include <cmath>

//typedef   std::uint8_t _8bit;
//typedef 1
//typedef    _8bit bufferType; //we'll just default everything to 16 bit sound right now
#include "waves.h"

#include <boost/function.hpp>

double round(double in){

	double tmp=0.0;
	double frac = modf(in, &tmp);

	return frac > .5? ceil(in):tmp;

}

double sampleToTime(uint sample, uint samplingRate){

	return double(sample) / double(samplingRate); 
	
}


double returnConstant(double time, double constant){
	 
	return constant;

}
 

double lerp(double begin, double end, double t){
	return (1.0 - t) * begin + t * end;
}

///\todo this is linear: should use different interpolation functions: pass in as an argument
double ASDRAmplitude(
	double time, 
	double fullAmplitude, 
	double sustainAmplitude, 
	double releasePoint,
	double attackTime, 
	double decayTime,
	double sustainTime, 
	double releaseTime){

		double vol;
		//go from nothing to full volume
		double attackT = std::min<double>(  time/attackTime, 1.0);///\todo max here makes cool harmonic
		vol = lerp(0.0, fullAmplitude,attackT );
		
		//now go from full volume to sustain
		double abc = (time - attackTime) / decayTime;

		abc = std::max(0.0, abc);
		abc = std::min(1.0, abc);

		vol = lerp(vol, sustainAmplitude, abc);
	
		//now go from sustain to 0
		double def = (time - releasePoint) / releaseTime;

		def = std::max(0.0, def);
		def = std::min(1.0, def);

		vol = lerp(vol, 0.0, def);



		//double linT =time/releasePoint;
		//vol = lerp(1.0, 0.0, std::min(1.0, linT ));


		//std::cout<<"at time "<<time<<" val is "<<vol<<std::endl;

		return vol;
}





struct Instrument{

	boost::function<double(uint,uint)> tOfSample;
	boost::function<double (double, double, double)> wave;
	boost::function<double (double)> fOfT;
	boost::function<double (double)> aOfT;

	//default initialize everything to the default instrument
	Instrument()
	
		:tOfSample(boost::bind(sampleToTime,_1,_2)),
		wave(&sineWave),
		fOfT(boost::bind(returnConstant, _1, 440.0)),
		aOfT(boost::bind(returnConstant, _1, 1.0))

	{
	
	}

	//this is the 'pipeline' - each of the subfunctions called is a shader when we port this to opencl
	double getValueAtSample(uint i, uint samplingRate){
	
		double sampleTime = tOfSample(i, samplingRate);	

		double amplitude = aOfT(sampleTime);

		double frequency = fOfT(sampleTime);

		double val =  wave(amplitude, frequency, sampleTime);

		val = std::min<double>(val, 1.0);
		val = std::max<double>(val,-1.0);

		return val;


	}



};



template <class BufferType>
class SoundBuffer{
 
	std::vector<BufferType> buffer;
	double length;
	uint samplingRate;

public:

	std::vector<BufferType>& getBuffer(){
		return buffer;
	}

	uint getSamplingRate(){
		return samplingRate;
	}

	uint totalSamples(){
		return buffer.size();
	}

	double getDuration(){
		return length;
	} 

	SoundBuffer(double lengthSeconds, uint samplesPerSecond =44100 )
		:buffer(lengthSeconds * samplesPerSecond, 0),
		samplingRate(samplesPerSecond), length(lengthSeconds)
	{
	
	}

	///\todo bounds checking
	void updateRange(double beginTime, double endTime, Instrument& instrument){
		
		uint start = beginTime * samplingRate;
		uint end = endTime * samplingRate;
		
		updateRange(start, end, instrument);
		
	}
	 
	

//	template <class T> 
//	T toOutputFormat(double in){
//		T tmp;
//		return tmp;
//	}
	
	BufferType toOutputFormat(double d);


	

	void updateRange(uint beginSample, uint endSample, Instrument& instrument){
	
		int sample = 0;
		for(int i = beginSample; i < endSample; i++, sample++){
		
			buffer[i] = toOutputFormat (instrument.getValueAtSample(sample, samplingRate));

		}

	}


};


//converts doubles to integer types
	 
	std::uint8_t SoundBuffer<std::uint8_t>::toOutputFormat (double in){
	 
		std::uint8_t tmp=std::numeric_limits<std::uint8_t>::max()>>1;	//openal sound library expects 0 to map to -1, MAX/2 to 0, and MAX to 1.0
		auto rval = static_cast<std::uint8_t>(round((in * tmp) + tmp));
		return rval;

	}
	 

	//converts doubles to integer types
	
	std::int16_t SoundBuffer<std::int16_t>::toOutputFormat(double in){
		std::int16_t tmp=std::numeric_limits<std::int16_t>::max();	 
		auto rval = static_cast<std::int16_t>(round((in * tmp)));
		return rval;

	}


typedef SoundBuffer<std::uint8_t> SoundBuffer8Bit;
typedef SoundBuffer<std::int16_t> SoundBuffer16Bit;
 

template <class T>
void initializeBuffer(T& buffer, double seconds, uint samplesPerSecond){

	uint bSize =static_cast<uint>( samplesPerSecond * seconds );

	buffer.resize(bSize);

}


double sawtoothWave(unsigned int pos, double maxAmplitude, double frequency, uint samplingRate){

	double secondsIn = (double(pos) / double(samplingRate));
	double tmp = (secondsIn * frequency);

	//return  maxAmplitude*((tmp - floor(tmp)));
	return maxAmplitude* ((tmp - floor(tmp)) *2.0) - 1.0;
	 
}


double sineWave(unsigned int pos, double maxAmplitude, double frequency, uint samplingRate){

	const double pi2 =2.0* 3.14159265;
	double sinVal = sin( ((pi2 *(frequency) /  double(samplingRate)))*double(pos));
	
	return sinVal *maxAmplitude;

}


double tanWave(unsigned int pos, double maxAmplitude, double frequency, uint samplingRate){
	const double pi2 =2.0* 3.14159265;
	double tanVal = tan( ((pi2 *(frequency) /  double(samplingRate)))*double(pos));
	
	return tanVal *maxAmplitude;
}




/*
template<class IT>
void boxFilter(IT begin, IT end, int passes=1){

	for(int p = 0; p<passes;p++){
	
		for(IT it = begin; it != end;){
		

			IT old = it;

			it++;

			*old=*it;

			//IT old = it;
		//	it++;



			//double tmp= .5 * (*one + *(it));
			//*old=tmp;

		}

	}
	

}*/


///\todo inv frequency precalculate
double squareWave(unsigned int pos, double maxAmplitude, double frequency, uint samplingRate){

	//find the position within the wave as a fract.  if < .5, return maxAmplitude.  if >.5, return -maxAmplitude

	double secondsIn = (double(pos) / samplingRate);
	double tmp;
	double posInWave=std::modf(secondsIn * frequency,&tmp);
	
	double rval =posInWave >.5?-maxAmplitude : maxAmplitude;

	//std::cout<<rval<<std::endl;

	return rval;
	//double val=sineWave(  pos,   maxAmplitude,   .10*frequency,   samplingRate);
//	return val>0.50?0.0:1.0;

	 
	//return (double(pos % samplingRate) / (double)samplingRate) *frequency

}




double whiteNoise(){

	return ((double(rand()) / RAND_MAX) - .5) * 2.0;


}

 

const double majorScaleCoeff[] = {1.0, 9.0/8.0,5.0/4.0, 4.0/3.0,3.0/2.0,5.0/3.0,15.0/8.0,2.0};
const uint octave = 7;
const uint fifth = 4;
const uint fourth = 3;
const uint rootNote = 0;
const uint third = 2;



///\todo seems to run at half pitch too
///\todo nonsquare periodic noise.  random, smooth noise, and weighted
double squareNoise(unsigned int pos, double maxAmplitude, double frequency, uint samplingRate){

	double val=sineWave(  pos,   maxAmplitude,   frequency,   samplingRate);
	
	static std::vector<double> table;	

	if(!table.size()){
		initTable(table,frequency);
	}

//	double randVal = double(rand())  / RAND_MAX;

	double secondsIn = (double(pos) / samplingRate);

	double tmp;
	double posInWave=std::modf(secondsIn * frequency,&tmp);
	
	uint idx =tmp*2;

	idx+= posInWave>.5?1:0;

	//double rval =posInWave >.5?-maxAmplitude : maxAmplitude;
	idx %= table.size();

	double randVal = table[idx];

	return maxAmplitude*randVal;//6(val>0.0?-1.0:1.0);

}


template<class IT>
void updateRange(IT begin, IT end, double freq)
{ 
	for(IT it = begin; it != end;it++){
	
		const uint samplingRate = 44100;
		 
		//using full amplitude

		auto pos = it - begin;

		double val;
		
		if(  pos == (end-begin)/2){
		
			freq *= majorScaleCoeff[1];

		}

		/*
		val= sineWave (pos, .250, freq, samplingRate);
		val+= squareWave (pos, .250, ( freq*majorScaleCoeff[third])  , samplingRate);
		val+= sawtoothWave (pos, .250, ( freq*majorScaleCoeff[fifth])  , samplingRate);
		val+= sineWave (pos, .250, ( freq*majorScaleCoeff[octave])  , samplingRate);
		*/
		//val=   sawtoothWave (pos, 1.0, freq, samplingRate);
		//val *= 1 + .335*whiteNoise();

		if(true){
		if(  pos < (end-begin)/2){
			//val= tanWave (pos, 1.0, freq, samplingRate);
			
			//val=   squareNoise (pos, .80, freq, samplingRate)
			//	+(squareNoise (pos, .20, 20.0*freq, samplingRate))
				//;
			//val*=val;
			val=   squareNoise (pos, 1.0, freq, samplingRate)//*
				//abs(sineWave (pos, 1.0, freq, samplingRate))
				;
			///that's an interesting instrument
		}else{
			//val= sineWave (pos, .75, freq, samplingRate);

			//val= squareWave(pos, .15, freq * majorScaleCoeff[third],samplingRate);
			
			//val+= sawtoothWave(pos, .15, freq * majorScaleCoeff[octave],samplingRate);
			
			//double tfreq =freq;
			double fcoeff = 1.0;
			
			double amp = 1.0/4.0;

			
			//val = amp*sineWave (pos , 1.0, freq * fcoeff, samplingRate);
			//val += .05*amp * sineWave(pos, 1.0, freq * majorScaleCoeff[2], samplingRate);
			
			
			//val = (.99*val + .01*(double(rand()) / RAND_MAX )) ;

			//val = (.025 * whiteNoise()) + (.95 * val) ;
			val = 0;
			 
			
			for(int i = 0; i <4; i++){
				
				val += squareNoise (pos , amp, freq * fcoeff, samplingRate);

				//amp *=.5;
				fcoeff+=.50;

			}
		 
		}
		}

		val = std::min<double>(val, 1.0);
		val = std::max<double>(val,-1.0);
		 

	//	std::cout<<val<<std::endl;

		///\todo:: to get what i had before, double the max val, not half
		double maxVal = std::numeric_limits<std::uint8_t>::max() >>  1 ; //half, since neutral is in the middle of the range

		//std::cout<<uint(maxVal)<<std::endl;
		//system("pause");

		*it =    round((maxVal * (val) )+ maxVal); 

	//	std::cout<<*it<<std::endl;
		  

	}
}
#endif