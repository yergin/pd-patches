#include <m_pd.h>
#include <string.h>
#include <math.h>
#include "Accelerate/Accelerate.h"

#define PI 3.14159265359f

static t_class* mu_tilde_class;	// pointer to the mu~ class

#define DF2
#define TABLE_LOOKUP

class Allpass2
{
public:
	Allpass2();
	
	inline float process( float in, float c, float k );
	void clear();
	
private:
#ifdef DF2
	float mW1;
	float mW2;
#else
	float mX1;
	float mX2;
	float mY1;
	float mY2;
#endif
};

Allpass2::Allpass2()
{
	clear();
}

inline float Allpass2::process( float in, float c, float k )
{
	// a1 = -k, a2 = c, b0 = -c, b1 = k, b2 = 1
#ifdef DF2
	float W0 = in - k * mW1 + c * mW2;
	float out = -c * W0 + k * mW1 + mW2;
	mW2 = mW1;
	mW1 = W0;
	return out;
#else
	float out = -c * in + k * mX1 + mX2 - k * mY1 + c * mY2;
	mX2 = mX1;
	mX1 = in;
	mY2 = mY1;
	mY1 = out;
	return out;
#endif
}

void Allpass2::clear()
{
#ifdef DF2
	mW1 = 0;
	mW2 = 0;
#else
	mX1 = 0;
	mX2 = 0;
	mY1 = 0;
	mY2 = 0;
#endif
}


class CircularBuffer
{
public:
	CircularBuffer( int length );
	~CircularBuffer();
	
	void write( float );
	float advanceAndRead();
	void clear();
	
private:
	int mLength;
	float* mBuffer;
	int mPosition;
};

CircularBuffer::CircularBuffer( int length )
	: mLength(length), mPosition(0)
{
	mBuffer = new float[length];
	clear();
}

CircularBuffer::~CircularBuffer()
{
	delete[] mBuffer;
}

void CircularBuffer::write( float sample )
{
	mBuffer[mPosition] = sample;
}

float CircularBuffer::advanceAndRead()
{
	mPosition = (mPosition + 1) % mLength;
	return mBuffer[mPosition];
}

void CircularBuffer::clear()
{
	for(int i = 0; i < mLength; ++i)
		mBuffer[i] = 0;
}


class TableLookup
{
public:
	TableLookup()
		: mSize(0), mTable(0)
	{ }
	
	void generate( float min, float max, int size );
	inline float lookup( float ) const;
	inline void lookup( float* in, float* out, int length );
	
	void debug();
	
protected:
	virtual float function( float ) const = 0;
	
private:
	int mSize;
	float* mTable;
	float mOffset;
	float mScale;
	float mMin;
	float mMax;
};

void TableLookup::generate( float min, float max, int size )
{
	mSize = size < 2 ? 2 : size;
	delete mTable;
	mMin = min;
	mMax = max;
	mTable = new float[mSize];
	mScale = static_cast<float>(mSize - 1) / (mMax - mMin);
	mOffset = -mMin * mScale;
	
	float invScale = 1.0f / mScale;
	for(int i = 0; i < mSize; ++i)
		mTable[i] = function((static_cast<float>(i) - mOffset) * invScale);
}

inline float TableLookup::lookup( float value ) const
{
	if(value <= mMin)
		return mTable[0];
	else if(value >= mMax)
		return mTable[mSize-1];
		
	float f = value * mScale + mOffset;
	int index = (static_cast<int>(f));
	
	return mTable[index] + (mTable[index + 1] - mTable[index]) * (f - static_cast<float>(index));
}

inline void TableLookup::lookup( float* in, float* out, int length )
{
#if 0
	vDSP_vtabi(in, 1, &mScale, &mOffset, mTable, mSize, out, 1, length);
#else
	while(length--)
		*out++ = lookup(*in++);
#endif
}

void TableLookup::debug()
{
	post("TableLookup table:");
	
	float invScale = 1.0f / mScale;
	for(int i = 0; i < mSize; ++i)
	{
		float val = (static_cast<float>(i) - mOffset) * invScale;
		post("%d:\t%f\t-> %f\t~ %f", i, val, function(val), lookup(val));
	}
}


class DampLookup : public TableLookup
{
public:
	DampLookup() { generate(0, 2.3, 64); }
	
protected:
	// TableLookup
	float function( float value ) const { return (tanf(value) - 1) / (tanf(value) + 1); }
};

class KLookup : public TableLookup
{
public:
	KLookup() { setSampleRate(44100); }
	void setSampleRate( float );
	
protected:
	// TableLookup
	float function( float value ) const { return -cosf(value * 2 * PI / mSampleRate); }
	
private:
	float mSampleRate;
};

void KLookup::setSampleRate( float samplerate )
{
	mSampleRate = samplerate;
	generate(0, mSampleRate / 2, 128);
}


typedef struct _mu_tilde
{
    t_object mPdObject;	// pd object structure

	Allpass2** mAllpass2Array;
	int mAllpass2Count;
	CircularBuffer* mBuffer;
	float* mDampBuffer;
	float* mKBuffer;
	int mBlockSize;
	float mFeedback;
	float mDampFactor;
	float mFreqFactor;
} t_mu_tilde;

static const char* version = "mu~ resonator v1.0";

static DampLookup dampLookup;
static KLookup kLookup;

extern "C"
{
	void mu_tilde_setup();
	void* mu_tilde_new( t_floatarg, t_floatarg, t_floatarg, t_floatarg );
	void mu_tilde_free( t_mu_tilde* );
	void mu_tilde_setfeedback( t_mu_tilde*, t_floatarg );
	void mu_tilde_setdensity( t_mu_tilde*, t_floatarg );
	void mu_tilde_clear( t_mu_tilde* );
	void mu_tilde_dsp( t_mu_tilde*, t_signal** );
	t_int* mu_tilde_perform( t_int* );
}

void* mu_tilde_new( t_floatarg allpassCount, t_floatarg loopDelay, t_floatarg density, t_floatarg feedback )
{
    t_mu_tilde* obj = (t_mu_tilde*)pd_new(mu_tilde_class);

	inlet_new(&obj->mPdObject, &obj->mPdObject.ob_pd, &s_signal, &s_signal);

    outlet_new(&obj->mPdObject, gensym("signal"));
    outlet_new(&obj->mPdObject, gensym("signal"));

	obj->mAllpass2Count = (int)allpassCount;
	obj->mAllpass2Array = new Allpass2*[obj->mAllpass2Count];
	for(int ap = 0; ap < obj->mAllpass2Count; ++ap)
		obj->mAllpass2Array[ap] = new Allpass2;
	obj->mBuffer = new CircularBuffer(loopDelay < 1 ? 1 : loopDelay);
	mu_tilde_setfeedback(obj, feedback);
	mu_tilde_setdensity(obj, (density == 0) ? 1 : density);
	obj->mFreqFactor = 2 * PI / sys_getsr();
	
	obj->mBlockSize = 64;
	obj->mDampBuffer = new float[obj->mBlockSize];
	obj->mKBuffer = new float[obj->mBlockSize];
	
    return (void*)obj;
}

void mu_tilde_free( t_mu_tilde* obj )
{
	for(int ap = 0; ap < obj->mAllpass2Count; ++ap)
		delete obj->mAllpass2Array[ap];
	delete[] obj->mAllpass2Array;
	delete obj->mBuffer;
}

void mu_tilde_setfeedback( t_mu_tilde* obj, t_floatarg feedback )
{
	obj->mFeedback = feedback;
}

void mu_tilde_setdensity( t_mu_tilde* obj, t_floatarg density )
{
	obj->mDampFactor = PI / sys_getsr() * density;
}

void mu_tilde_clear( t_mu_tilde* obj )
{
	for(int i = 0; i < obj->mAllpass2Count; ++i)
		obj->mAllpass2Array[i]->clear();
	obj->mBuffer->clear();	
}

t_int* mu_tilde_perform( t_int* w )
{
	static const float DC = 1e-100;
	float one = 1.0f;
	
    t_mu_tilde* obj = (t_mu_tilde*)(w[1]);
	int n = (int)(w[2]);
	t_sample* audioIn = (t_sample*)(w[3]);
	t_sample* freqIn = (t_sample*)(w[4]);
	t_sample* outLeft = (t_sample*)(w[5]);
	t_sample* outRight = (t_sample*)(w[6]);
	
	if(n != obj->mBlockSize)
	{
		delete[] obj->mDampBuffer;
		delete[] obj->mKBuffer;
		obj->mBlockSize = n;
		obj->mDampBuffer = new float[obj->mBlockSize];
		obj->mKBuffer = new float[obj->mBlockSize];
	}
	
	/*
		c = tanf(*freqIn * obj->mDampFactor);
		c = (c - 1) / (c + 1);
		k = -cosf(*freqIn++ * obj->mFreqFactor) * (1 - c);
			OR
		c = dampLookup.lookup(*freqIn * obj->mDampFactor);
		k = kLookup.lookup(*freqIn++) * (1 - c);
	*/
	// use mKBuffer as temp for dampLookup
	vDSP_vsmul(freqIn, 1, &obj->mDampFactor, obj->mKBuffer, 1, n);
	dampLookup.lookup(obj->mKBuffer, obj->mDampBuffer, n);
	
	kLookup.lookup(freqIn, obj->mKBuffer, n);
	vDSP_vsbm(&one, 0, obj->mDampBuffer, 1, obj->mKBuffer, 1, obj->mKBuffer, 1, n);
	
	float in;
	float* c = obj->mDampBuffer;
	float* k = obj->mKBuffer;
	while(n--)
	{
		in = *audioIn++ + obj->mFeedback * obj->mBuffer->advanceAndRead() + DC;
		int ap = 0;
		while(ap != obj->mAllpass2Count - 1)
			in = obj->mAllpass2Array[ap++]->process(in, *c, *k);
		*outRight++ = in;
		*outLeft = obj->mAllpass2Array[ap]->process(in, *c++, *k++);
		obj->mBuffer->write(*outLeft++);
	}

    return (w + 7);
}

void mu_tilde_dsp( t_mu_tilde* obj, t_signal** sp )
{
    dsp_add(mu_tilde_perform, 6, obj, sp[0]->s_n, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec);
}

void mu_tilde_setup()
{
	kLookup.setSampleRate(sys_getsr());
	
    mu_tilde_class = class_new(gensym("mu~"), (t_newmethod)mu_tilde_new, (t_method)mu_tilde_free, sizeof(t_mu_tilde), CLASS_DEFAULT, A_FLOAT, A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT, (t_atomtype)0);

    class_addmethod(mu_tilde_class, nullfn, gensym("signal"), (t_atomtype)0);
    class_addmethod(mu_tilde_class, (t_method)mu_tilde_dsp, gensym("dsp"), (t_atomtype)0);
    class_addmethod(mu_tilde_class, (t_method)mu_tilde_setfeedback, gensym("feedback"), A_FLOAT, (t_atomtype)0);
    class_addmethod(mu_tilde_class, (t_method)mu_tilde_setdensity, gensym("density"), A_FLOAT, (t_atomtype)0);
    class_addmethod(mu_tilde_class, (t_method)mu_tilde_clear, gensym("clear"), (t_atomtype)0);

	post(version);
#if 0
	dampLookup.debug();
	kLookup.debug();
#endif
}
