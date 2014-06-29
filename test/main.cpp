
#include <iostream>
#include <al.h>
#include <alc.h>
#include <vector>
#include <cstdint>

#include "tfunctions.h"

template <class T> 
void playSound(SoundBuffer<T>& buff){

	ALCdevice* device=alcOpenDevice(0); //takes string identifying device
	//todo: use enumeration extension to query which sound device user wants to use.
	//should be on menu, and should have a test

	if(!device){
	
		std::cout<<"Unable to open device!"<<std::endl;
		system("pause");
		exit(0);
	}

	ALCcontext* context=alcCreateContext(device,NULL); 
	alcMakeContextCurrent(context);  

	if(alGetError()){
		std::cout<<"Problem during setup"<<std::endl;
		system("pause");
		exit(0);
	}

	ALuint buffer;
	alGenBuffers(1,&buffer);

	if(alGetError()){///\todo have a proper error handling method that has a switch statement, logs exceptions
		std::cout<<"Problem creating buffer"<<std::endl;
		system("pause");
		exit(0);
	}



	auto enumVal = sizeof(T) == 1 ? AL_FORMAT_MONO8: (AL_FORMAT_MONO16);///\todo set this to 8: WHAT THE FUCK!

	std::cout<<"Size is "<<sizeof(T)<<std::endl;
	system("pause");

	///\todo this enum is hard coded to 16 when buffer type could be 8 or something else
	alBufferData(buffer, enumVal, &buff.getBuffer()[0], buff.getBuffer().size() * sizeof(T),buff.getSamplingRate()); 

	ALuint source;
	alGenSources(1, &source);

	alSourcei(source, AL_BUFFER,buffer);

	alSourcePlay(source);
	
	system("pause");

}



int main(void){

	double tones[15] = {261.63,261.63, 233.08, 261.63, 293.66, 311.13, 293.66, 233.08, 261.63, 261.63, 233.08, 261.63, 349.23, 311.13, 293.66 };

	for(int i = 0; i < 15; i++){
	
		//tones[i]*=.50;

	}

	SoundBuffer16Bit buff(15.0); //create a sound buffer 10 seconds long

	Instrument inst; //will default to an instrument that always plays a 440 hz A at full volume

	inst.aOfT = boost::bind(ASDRAmplitude, 
		_1, //time
		.50, //full volume
		.25, //sustain volume
		.80, //release point
		.250, //attack time
		.25,	//.25, //decay time
		.25, //sustain time
		.25);//release time

	for(int i = 0 ; i < 15; i++){
	
		inst.fOfT = boost::bind(returnConstant, _1, tones[i]);	

		buff.updateRange(float(i), float(i+1), inst);

	}

	//buff.updateRange(0.0, 10.0, inst);

	playSound(buff);

	return 0;

}



int main2(void){

	const double root=100;//25;//261.63;
	const double tone = majorScaleCoeff[fifth]*root;

	AudioSample8 test;
	initializeBuffer(test, 4.0,44100);
	updateRange(test.begin(), test.end(),root);

	/*
	int val = 0;
	int its =0;
	for(int i = 0; i < test.size(); i++){
	
		int tval = test[i];

		if(tval != val){
		
			val = tval;
			its++;
		}

	}

	std::cout<<its<<std::endl;
	system("pause");
	*/

	AudioSample8 test2;
	initializeBuffer(test2, 2.0,44100);
	updateRange(test2.begin(), test2.end(),tone);
	
	//need offset

	ALCdevice* device=alcOpenDevice(0); //takes string identifying device
	//todo: use enumeration extension to query which sound device user wants to use.
	//should be on menu, and should have a test
	

	if(!device){
	
		std::cout<<"Unable to open device!"<<std::endl;
		system("pause");
		return 0;

	}

	ALCcontext* context=alcCreateContext(device,NULL); 
	alcMakeContextCurrent(context);  

	if(alGetError()){
		std::cout<<"Problem during setup"<<std::endl;
		system("pause");
		return 0;
	}

	ALuint buffer;
	alGenBuffers(1,&buffer);

	ALuint buffer2;
	alGenBuffers(1,&buffer2);



	if(alGetError()){///\todo have a proper error handling method that has a switch statement, logs exceptions
		std::cout<<"Problem creating buffer"<<std::endl;
		system("pause");
		return 0;
	}

	alBufferData(buffer, AL_FORMAT_MONO8, &test[0], test.size() * sizeof(std::uint8_t),44100); 
	
	alBufferData(buffer2, AL_FORMAT_MONO8, &test2[0], test2.size() * sizeof(std::uint8_t),44100); 


	ALuint source;
	alGenSources(1, &source);

	ALuint source2;
	alGenSources(1, &source2);

	alSourcei(source, AL_BUFFER,buffer);

	alSourcei(source2, AL_BUFFER,buffer2);

	alSourcePlay(source);
	//alSourcePlay(source2);
	
	
	system("pause");
}
