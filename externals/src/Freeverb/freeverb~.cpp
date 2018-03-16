#include <m_pd.h>
#include "Accelerate/Accelerate.h"
#include <string.h>

#include "tuning.h"

class comb
{
public:
    comb();
	~comb();
    void	setbuffer(float *buf, int size);
	inline  float process(float inp);
	void	setBlockSize(int);
	inline  void processAndSumBlock(float*, float*);
    void	mute();
    void	setdamp(float val);
    float	getdamp();
    void	setfeedback(float val);
    float	getfeedback();

private:
	inline  void processAndSum(float*, float*, int);
	
	float	feedback;
	float*	filterstore;
	float   lastfilterstore;
	float*	tempbuffer;
	float	*buffer;
	int		bufsize;
	int		bufidx;
	float	damp1;
	float	damp2;
	float   dc;
	int		blocksize;
};


// Big to inline - but crucial for speed

inline float comb::process(float input)
{
	float output;
    
	output = buffer[bufidx];
	lastfilterstore = (output * damp2) + (lastfilterstore * damp1) + dc;
	buffer[bufidx] = input + (lastfilterstore * feedback);
    
	if(++bufidx >= bufsize)
		bufidx = 0;
    
	return output;
}

inline void comb::processAndSumBlock(float* in, float* out)
{
	for(int samplesleft = blocksize; samplesleft > 0;)
	{
		int numsamples = bufsize - bufidx >= blocksize ? blocksize : bufsize - bufidx;
		if(numsamples > samplesleft)
			numsamples = samplesleft;
		processAndSum(in, out, numsamples);
		in += numsamples;
		out += numsamples;
		samplesleft -= numsamples;
	}
}

inline void comb::processAndSum(float* in, float* out, int numsamples)
{
	// out[0 to numsamples-1] += buffer[bufidx to bufidx+numsamples-1]
	vDSP_vadd(buffer + bufidx, 1, out, 1, out, 1, numsamples);
	
	// filterstore[n] = buffer[bufidx + n] * damp2 + filterstore[n - 1] * damp1
	vDSP_vsmsa(buffer + bufidx, 1, &damp2, &dc, tempbuffer, 1, numsamples);
	float* f = filterstore;
	float* g = filterstore;
	float* buf = tempbuffer;
	*f++ = (*buf++) + lastfilterstore * damp1;
	int count = numsamples - 1;
	while(count--)
		*f++ = (*buf++) + (*g++) * damp1;
	lastfilterstore = *g;
	
	// buffer[bufidx to bufidx+numsamples-1] = in[0 to numsamples-1] + feedback * filterstore[0 to numsamples-1]
	vDSP_vsma(filterstore, 1, &feedback, in, 1, buffer + bufidx, 1, numsamples);
	
	bufidx += numsamples;
	if(bufidx >= bufsize)
		bufidx = 0;
}

comb::comb()
: filterstore(0), lastfilterstore(0), tempbuffer(0), buffer(0), bufsize(0), bufidx(0), dc(0)
{
	setBlockSize(1);
}

comb::~comb()
{
	delete[] filterstore;
	delete[] tempbuffer;
}

void comb::setbuffer(float *buf, int size)
{
	buffer = buf; 
	bufsize = size;
}

void comb::setBlockSize( int size )
{
	blocksize = size;
	delete[] filterstore;
	delete[] tempbuffer;
	filterstore = new float[size];
	tempbuffer = new float[size];
	for(int i = 0; i < size; ++i)
		filterstore[i] = 0;
}

void comb::mute()
{
	for(int i = 0; i < bufsize; ++i)
		buffer[i] = 0;
	for(int i = 0; i < blocksize; ++i)
		filterstore[i] = 0;
}

void comb::setdamp(float val) 
{
	damp1 = val; 
	damp2 = 1-val;
}

float comb::getdamp() 
{
	return damp1;
}

void comb::setfeedback(float val) 
{
	feedback = val;
	dc = val < 1 ? 1e-100 : 0;
}

float comb::getfeedback() 
{
	return feedback;
}

class allpass
{
public:
    allpass();
    void	setbuffer(float *buf, int size);
    void setBlockSize(int size) { blocksize = size; }
	inline  float	process(float inp);
    inline void processBlock( float* input, float* output );
    inline void process( float* input, float* output, int samples );
    void	mute();
    void	setfeedback(float val);
    float	getfeedback();
    // private:
	float	feedback;
	float	*buffer;
	int		bufsize;
	int		bufidx;
    int blocksize;
	float	dc;
};


// Big to inline - but crucial for speed

inline float allpass::process(float input)
{
	float output;
	float bufout;
	
	bufout = buffer[bufidx] + dc;
	output = -input + bufout;
	buffer[bufidx] = input + (bufout*feedback);
    
	if(++bufidx>=bufsize) bufidx = 0;
    
	return output;
}

inline void allpass::processBlock(float* in, float* out)
{
	for(int samplesleft = blocksize; samplesleft > 0;)
	{
		int numsamples = bufsize - bufidx >= blocksize ? blocksize : bufsize - bufidx;
		if(numsamples > samplesleft)
			numsamples = samplesleft;
		process(in, out, numsamples);
		in += numsamples;
		out += numsamples;
		samplesleft -= numsamples;
	}
}

inline void allpass::process(float* in, float* out, int numsamples)
{
	//output = -input + buffer[bufidx];
    vDSP_vsub(in, 1, buffer + bufidx, 1, out, 1, numsamples);
    
	//buffer[bufidx] = input + buffer[bufidx] * feedback;
	vDSP_vsma(buffer + bufidx, 1, &feedback, in, 1, buffer + bufidx, 1, numsamples);
    
	bufidx += numsamples;
	if(bufidx >= bufsize)
		bufidx = 0;
}

allpass::allpass()
{
	bufidx = 0;
    blocksize = 1;
}

void allpass::setbuffer(float *buf, int size) 
{
	buffer = buf; 
	bufsize = size;
}

void allpass::mute()
{
	for (int i=0; i<bufsize; i++)
		buffer[i]=0;
}

void allpass::setfeedback(float val) 
{
	feedback = val;
	dc = val < 1 ? 1e-100 : 0;
}

float allpass::getfeedback() 
{
	return feedback;
}

class revmodel
{
public:
    revmodel();
    void	mute();
    void	processreplace(float *input, float *outputL, float *outputR, long numsamples);
    void	setroomsize(float value);
    float	getroomsize();
    void	setdamp(float value);
    float	getdamp();
    void	setwet(float value);
    float	getwet();
    void	setdry(float value);
    float	getdry();
    void	setwidth(float value);
    float	getwidth();
    void	setmode(float value);
    float	getmode();
    
    int blockSize() { return blocksize; }
    void setBlockSize( int );
    
private:
    void	update();
    
private:
	float	gain;
	float	roomsize,roomsize1;
	float	damp,damp1;
	float	mode;
	int   blocksize;
	float* swapBuffer;
    
	// The following are all declared inline 
	// to remove the need for dynamic allocation
	// with its subsequent error-checking messiness
    
	// Comb filters
	comb	combL[numcombs];
	comb	combR[numcombs];
    
	// Allpass filters
	allpass	allpassL[numallpasses];
	allpass	allpassR[numallpasses];
    
	// Buffers for the combs
	float	bufcombL1[combtuningL1];
	float	bufcombR1[combtuningR1];
	float	bufcombL2[combtuningL2];
	float	bufcombR2[combtuningR2];
	float	bufcombL3[combtuningL3];
	float	bufcombR3[combtuningR3];
	float	bufcombL4[combtuningL4];
	float	bufcombR4[combtuningR4];
	float	bufcombL5[combtuningL5];
	float	bufcombR5[combtuningR5];
	float	bufcombL6[combtuningL6];
	float	bufcombR6[combtuningR6];
	float	bufcombL7[combtuningL7];
	float	bufcombR7[combtuningR7];
	float	bufcombL8[combtuningL8];
	float	bufcombR8[combtuningR8];
    
	// Buffers for the allpasses
	float	bufallpassL1[allpasstuningL1];
	float	bufallpassR1[allpasstuningR1];
	float	bufallpassL2[allpasstuningL2];
	float	bufallpassR2[allpasstuningR2];
	float	bufallpassL3[allpasstuningL3];
	float	bufallpassR3[allpasstuningR3];
	float	bufallpassL4[allpasstuningL4];
	float	bufallpassR4[allpasstuningR4];
};

revmodel::revmodel()
: swapBuffer(0)
{
	// Tie the components to their buffers
	combL[0].setbuffer(bufcombL1,combtuningL1);
	combR[0].setbuffer(bufcombR1,combtuningR1);
	combL[1].setbuffer(bufcombL2,combtuningL2);
	combR[1].setbuffer(bufcombR2,combtuningR2);
	combL[2].setbuffer(bufcombL3,combtuningL3);
	combR[2].setbuffer(bufcombR3,combtuningR3);
	combL[3].setbuffer(bufcombL4,combtuningL4);
	combR[3].setbuffer(bufcombR4,combtuningR4);
	combL[4].setbuffer(bufcombL5,combtuningL5);
	combR[4].setbuffer(bufcombR5,combtuningR5);
	combL[5].setbuffer(bufcombL6,combtuningL6);
	combR[5].setbuffer(bufcombR6,combtuningR6);
	combL[6].setbuffer(bufcombL7,combtuningL7);
	combR[6].setbuffer(bufcombR7,combtuningR7);
	combL[7].setbuffer(bufcombL8,combtuningL8);
	combR[7].setbuffer(bufcombR8,combtuningR8);
	allpassL[0].setbuffer(bufallpassL1,allpasstuningL1);
	allpassR[0].setbuffer(bufallpassR1,allpasstuningR1);
	allpassL[1].setbuffer(bufallpassL2,allpasstuningL2);
	allpassR[1].setbuffer(bufallpassR2,allpasstuningR2);
	allpassL[2].setbuffer(bufallpassL3,allpasstuningL3);
	allpassR[2].setbuffer(bufallpassR3,allpasstuningR3);
	allpassL[3].setbuffer(bufallpassL4,allpasstuningL4);
	allpassR[3].setbuffer(bufallpassR4,allpasstuningR4);
    
	// Set default values
	allpassL[0].setfeedback(0.5f);
	allpassR[0].setfeedback(0.5f);
	allpassL[1].setfeedback(0.5f);
	allpassR[1].setfeedback(0.5f);
	allpassL[2].setfeedback(0.5f);
	allpassR[2].setfeedback(0.5f);
	allpassL[3].setfeedback(0.5f);
	allpassR[3].setfeedback(0.5f);
	setroomsize(initialroom);
	setdamp(initialdamp);
	setmode(initialmode);
    
    setBlockSize(64);
    
	// Buffer will be full of rubbish - so we MUST mute them
	mute();
}

void revmodel::setBlockSize( int size )
{
	for(int i = 0; i < numcombs; ++i)
	{
		combL[i].setBlockSize(size);
		combR[i].setBlockSize(size);
	}
	for(int i = 0; i < numallpasses; ++i)
	{
		allpassL[i].setBlockSize(size);
		allpassR[i].setBlockSize(size);
	}
	delete[] swapBuffer;
	swapBuffer = new float[size];
	blocksize = size;
}

void revmodel::mute()
{
	if (getmode() >= freezemode)
		return;
    
	for (int i=0;i<numcombs;i++)
	{
		combL[i].mute();
		combR[i].mute();
	}
	for (int i=0;i<numallpasses;i++)
	{
		allpassL[i].mute();
		allpassR[i].mute();
	}
}

void revmodel::processreplace(float *input, float *outputL, float *outputR, long numsamples)
{
	if(numsamples != blockSize())
	{
		post("blueroom~: Block size changed to %d!", numsamples);
		setBlockSize(numsamples);
	}

	vDSP_vsmul(input, 1, &gain, swapBuffer, 1, numsamples);
	
	vDSP_vclr(outputL, 1, numsamples);
	vDSP_vclr(outputR, 1, numsamples);

	combL[0].processAndSumBlock(swapBuffer, outputL);
	combR[0].processAndSumBlock(swapBuffer, outputR);
	combL[1].processAndSumBlock(swapBuffer, outputL);
	combR[1].processAndSumBlock(swapBuffer, outputR);
	combL[2].processAndSumBlock(swapBuffer, outputL);
	combR[2].processAndSumBlock(swapBuffer, outputR);
	combL[3].processAndSumBlock(swapBuffer, outputL);
	combR[3].processAndSumBlock(swapBuffer, outputR);
	combL[4].processAndSumBlock(swapBuffer, outputL);
	combR[4].processAndSumBlock(swapBuffer, outputR);
	combL[5].processAndSumBlock(swapBuffer, outputL);
	combR[5].processAndSumBlock(swapBuffer, outputR);
	combL[6].processAndSumBlock(swapBuffer, outputL);
	combR[6].processAndSumBlock(swapBuffer, outputR);
	combL[7].processAndSumBlock(swapBuffer, outputL);
	combR[7].processAndSumBlock(swapBuffer, outputR);

    allpassL[0].processBlock(outputL, swapBuffer);
    allpassL[1].processBlock(swapBuffer, outputL);
    allpassL[2].processBlock(outputL, swapBuffer);
    allpassL[3].processBlock(swapBuffer, outputL);
    allpassR[0].processBlock(outputR, swapBuffer);
    allpassR[1].processBlock(swapBuffer, outputR);
    allpassR[2].processBlock(outputR, swapBuffer);
    allpassR[3].processBlock(swapBuffer, outputR);
}

void revmodel::update()
{
    // Recalculate internal values after parameter change
    
	int i;
    
	if (mode >= freezemode)
	{
		roomsize1 = 1;
		damp1 = 0;
		gain = muted;
	}
	else
	{
		roomsize1 = roomsize;
		damp1 = damp;
		gain = fixedgain;
	}
    
	for(i=0; i<numcombs; i++)
	{
		combL[i].setfeedback(roomsize1);
		combR[i].setfeedback(roomsize1);
	}
    
	for(i=0; i<numcombs; i++)
	{
		combL[i].setdamp(damp1);
		combR[i].setdamp(damp1);
	}
}

// The following get/set functions are not inlined, because
// speed is never an issue when calling them, and also
// because as you develop the reverb model, you may
// wish to take dynamic action when they are called.

void revmodel::setroomsize(float value)
{
	roomsize = (value*scaleroom) + offsetroom;
	update();
}

float revmodel::getroomsize()
{
	return (roomsize-offsetroom)/scaleroom;
}

void revmodel::setdamp(float value)
{
	damp = value*scaledamp;
	update();
}

float revmodel::getdamp()
{
	return damp/scaledamp;
}

void revmodel::setmode(float value)
{
	mode = value;
	update();
}

float revmodel::getmode()
{
	if (mode >= freezemode)
		return 1;
	else
		return 0;
}

static t_class* freeverb_tilde_class;	// pointer to the freeverb~ class

static const char* VERSION_STRING = "blueroom~ reverb v1.0";

typedef struct _freeverb_tilde
{
    t_object mPdObject;	// pd object structure

	revmodel* mModel;
} t_freeverb_tilde;

extern "C"
{
	void fverb_tilde_setup();
	void* freeverb_tilde_new();
	void freeverb_tilde_free( t_freeverb_tilde* );
	void freeverb_tilde_setroomsize( t_freeverb_tilde*, t_floatarg );
	void freeverb_tilde_setdamping( t_freeverb_tilde*, t_floatarg );
	void freeverb_tilde_setfreeze( t_freeverb_tilde*, t_floatarg );
	void freeverb_tilde_clear( t_freeverb_tilde* );
	void freeverb_tilde_dsp( t_freeverb_tilde*, t_signal** );
	t_int* freeverb_tilde_perform( t_int* );
}

void* freeverb_tilde_new()
{
    t_freeverb_tilde* obj = (t_freeverb_tilde*)pd_new(freeverb_tilde_class);

    outlet_new(&obj->mPdObject, gensym("signal"));
    outlet_new(&obj->mPdObject, gensym("signal"));

	obj->mModel = new revmodel;
	
    return (void*)obj;
}

void freeverb_tilde_free( t_freeverb_tilde* obj )
{
	delete obj->mModel;
}

void freeverb_tilde_setroomsize( t_freeverb_tilde* obj, t_floatarg roomsize )
{
	obj->mModel->setroomsize(roomsize);
}

void freeverb_tilde_setdamping( t_freeverb_tilde* obj, t_floatarg damping )
{
	obj->mModel->setdamp(damping);
}

void freeverb_tilde_setfreeze( t_freeverb_tilde* obj, t_floatarg freeze )
{
	obj->mModel->setmode(freeze);
}

void freeverb_tilde_clear( t_freeverb_tilde* obj )
{
	obj->mModel->mute();
}

t_int* freeverb_tilde_perform( t_int* w )
{
    t_freeverb_tilde* obj = (t_freeverb_tilde*)(w[1]);

	obj->mModel->processreplace((t_sample*)(w[3]), (t_sample*)(w[4]), (t_sample*)(w[5]), (int)(w[2]));

    return (w + 6);
}

void freeverb_tilde_dsp( t_freeverb_tilde* obj, t_signal** sp )
{
    dsp_add(freeverb_tilde_perform, 5, obj, sp[0]->s_n, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec);
}

void fverb_tilde_setup()
{
    freeverb_tilde_class = class_new(gensym("fverb~"), (t_newmethod)freeverb_tilde_new, (t_method)freeverb_tilde_free, sizeof(t_freeverb_tilde), CLASS_DEFAULT, (t_atomtype)0);

    class_addmethod(freeverb_tilde_class, nullfn, gensym("signal"), (t_atomtype)0);
    class_addmethod(freeverb_tilde_class, (t_method)freeverb_tilde_dsp, gensym("dsp"), (t_atomtype)0);
    class_addmethod(freeverb_tilde_class, (t_method)freeverb_tilde_setroomsize, gensym("roomsize"), A_FLOAT, (t_atomtype)0);
    class_addmethod(freeverb_tilde_class, (t_method)freeverb_tilde_setdamping, gensym("damping"), A_FLOAT, (t_atomtype)0);
    class_addmethod(freeverb_tilde_class, (t_method)freeverb_tilde_setfreeze, gensym("freeze"), A_FLOAT, (t_atomtype)0);
    class_addmethod(freeverb_tilde_class, (t_method)freeverb_tilde_clear, gensym("clear"), (t_atomtype)0);
    
    post(VERSION_STRING);
}
