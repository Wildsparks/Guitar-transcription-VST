/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include <string> 

//==============================================================================
Jacode_iiiAudioProcessorEditor::Jacode_iiiAudioProcessorEditor (Jacode_iiiAudioProcessor& p) //constructor
    : AudioProcessorEditor (&p), processor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.

	//size of the windows
    setSize (1600, 700);

	//FPS
	startTimerHz(30);  //interval en miliseconde divisé par 8 pour faire des 60000 / BPM / 8

	//background
	background = ImageCache::getFromMemory(BinaryData::backgroundImagev3_png, BinaryData::backgroundImage_pngSize);

	//--------------------------------//
	//--------onset sliders-----------//
	//--------------------------------//

	//--------------------------------//
	addAndMakeVisible(thresholdSlider);
	thresholdSlider.setRange(0, 100);            
	thresholdSlider.setTextValueSuffix(" %");    
	thresholdSlider.setValue(25.0);
	thresholdSlider.addListener(this);
	thresholdSlider.setColour(Slider::textBoxTextColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(thresholdLabel);
	thresholdLabel.setText("Threshold :", dontSendNotification);
	thresholdLabel.attachToComponent(&thresholdSlider, true); 
	thresholdLabel.setColour(Label::textColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(betaresSlider);
	betaresSlider.setRange(1, 8, 1);
	betaresSlider.setTextValueSuffix("");
	betaresSlider.setValue(5.0);
	betaresSlider.addListener(this);
	betaresSlider.setColour(Slider::textBoxTextColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(betaresLabel);
	betaresLabel.setText("betares : 1e-", dontSendNotification);
	betaresLabel.attachToComponent(&betaresSlider, true);
	betaresLabel.setColour(Label::textColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(nfftSlider);
	nfftSlider.setRange(10, 22, 1);
	nfftSlider.setTextValueSuffix("");
	nfftSlider.setValue(17.0);
	nfftSlider.addListener(this);
	nfftSlider.setColour(Slider::textBoxTextColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(nfftLabel);
	nfftLabel.setText("nfft : 2^", dontSendNotification);
	nfftLabel.attachToComponent(&nfftSlider, true);
	nfftLabel.setColour(Label::textColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(distanceMaxPluck);
	distanceMaxPluck.setRange(6.000001, 64.3);
	distanceMaxPluck.setTextValueSuffix(" cm");
	distanceMaxPluck.setValue(64.3);
	distanceMaxPluck.addListener(this);
	distanceMaxPluck.setColour(Slider::textBoxTextColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(distanceMaxPluckLabel);
	distanceMaxPluckLabel.setText("distance :", dontSendNotification);
	distanceMaxPluckLabel.attachToComponent(&distanceMaxPluck, true);
	distanceMaxPluckLabel.setColour(Label::textColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(resPluck);
	resPluck.setRange(1, 3, 1);
	resPluck.setTextValueSuffix(" cm");
	resPluck.setValue(1.0);
	resPluck.addListener(this);
	resPluck.setColour(Slider::textBoxTextColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(resPluckLabel);
	resPluckLabel.setText("resolution : 1e-", dontSendNotification);
	resPluckLabel.attachToComponent(&resPluck, true);
	resPluckLabel.setColour(Label::textColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(nbHarmo);
	nbHarmo.setRange(10, 30, 1);
	nbHarmo.setTextValueSuffix("");
	nbHarmo.setValue(20.0);
	nbHarmo.addListener(this);
	nbHarmo.setColour(Slider::textBoxTextColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(nbHarmoLabel);
	nbHarmoLabel.setText("Harmonics Betta :", dontSendNotification);
	nbHarmoLabel.attachToComponent(&nbHarmo, true);
	nbHarmoLabel.setColour(Label::textColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(nbHarmoF0);
	nbHarmoF0.setRange(2, 8, 1);
	nbHarmoF0.setTextValueSuffix("");
	nbHarmoF0.setValue(5.0);
	nbHarmoF0.addListener(this);
	nbHarmoF0.setColour(Slider::textBoxTextColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(nbHarmoF0Label);
	nbHarmoF0Label.setText("Harmonics F0 :", dontSendNotification);
	nbHarmoF0Label.attachToComponent(&nbHarmoF0, true);
	nbHarmoF0Label.setColour(Label::textColourId, Colours::black);
	//--------------------------------//

	addAndMakeVisible(string_1);
	string_1.addListener(this);
	string_1.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_2);
	string_2.addListener(this);
	string_2.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_3);
	string_3.addListener(this);
	string_3.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_4);
	string_4.addListener(this);
	string_4.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_5);
	string_5.addListener(this);
	string_5.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_6);
	string_6.addListener(this);
	string_6.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_7);
	string_7.addListener(this);
	string_7.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_8);
	string_8.addListener(this);
	string_8.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_9);
	string_9.addListener(this);
	string_9.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_10);
	string_10.addListener(this);
	string_10.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_11);
	string_11.addListener(this);
	string_11.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_12);
	string_12.addListener(this);
	string_12.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_13);
	string_13.addListener(this);
	string_13.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_14);
	string_14.addListener(this);
	string_14.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_15);
	string_15.addListener(this);
	string_15.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_16);
	string_16.addListener(this);
	string_16.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_17);
	string_17.addListener(this);
	string_17.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_18);
	string_18.addListener(this);
	string_18.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_19);
	string_19.addListener(this);
	string_19.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_20);
	string_20.addListener(this);
	string_20.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_21);
	string_21.addListener(this);
	string_21.setSize(40, 15);
	//---..---//
	addAndMakeVisible(string_22);
	string_22.addListener(this);
	string_22.setSize(40, 15);

	//--------------------------------//
	addAndMakeVisible(DfullLabel);
	DfullLabel.setText("D_full :", dontSendNotification);
	DfullLabel.attachToComponent(&string_1, true);
	DfullLabel.setColour(Label::textColourId, Colours::black);

	addAndMakeVisible(DcoreLabel);
	DcoreLabel.setText("D_core :", dontSendNotification);
	DcoreLabel.attachToComponent(&string_7, true);
	DcoreLabel.setColour(Label::textColourId, Colours::black);

	addAndMakeVisible(TensionLabel);
	TensionLabel.setText("Tension :", dontSendNotification);
	TensionLabel.attachToComponent(&string_13, true);
	TensionLabel.setColour(Label::textColourId, Colours::black);

	addAndMakeVisible(rhoCore);
	rhoCore.setText("RhoCore :", dontSendNotification);
	rhoCore.attachToComponent(&string_19, true);
	rhoCore.setColour(Label::textColourId, Colours::black);

	addAndMakeVisible(rhoWrapping);
	rhoWrapping.setText("RhoWrap:", dontSendNotification);
	rhoWrapping.attachToComponent(&string_20, true);
	rhoWrapping.setColour(Label::textColourId, Colours::black);

	addAndMakeVisible(steelShearModulus);
	steelShearModulus.setText("ShearMod :", dontSendNotification);
	steelShearModulus.attachToComponent(&string_21, true);
	steelShearModulus.setColour(Label::textColourId, Colours::black);

	addAndMakeVisible(steelYoungModulus);
	steelYoungModulus.setText("YoungMod :", dontSendNotification);
	steelYoungModulus.attachToComponent(&string_22, true);
	steelYoungModulus.setColour(Label::textColourId, Colours::black);

	//--------------------------------//
	addAndMakeVisible(reloadButton);
	reloadButton.setSize(200, 25);
	reloadButton.setColour(ToggleButton::textColourId, Colours::black);
	reloadButton.setButtonText("Reload model");
	reloadButton.addListener(this);

	addAndMakeVisible(switchPluck);
	switchPluck.setButtonText("activate plucking position");
	switchPluck.setSize(150, 50);
	switchPluck.setColour(ToggleButton::textColourId, Colours::black);
	switchPluck.setColour(ToggleButton::tickColourId, Colours::black); 
	switchPluck.setColour(ToggleButton::tickDisabledColourId, Colours::black);
	switchPluck.addListener(this);

	string_1.setText("0.010");
	string_2.setText("0.013");
	string_3.setText("0.0172");
	string_4.setText("0.026");
	string_5.setText("0.036");
	string_6.setText("0.046");

	string_7.setText("0.010");
	string_8.setText("0.013");
	string_9.setText("0.0172");
	string_10.setText("0.0146");
	string_11.setText("0.016");
	string_12.setText("0.018");

	string_13.setText("16.5");
	string_14.setText("16");
	string_15.setText("17.5");
	string_16.setText("18.5");
	string_17.setText("20");
	string_18.setText("17");

	string_19.setText("7950");

	string_20.setText("6000");

	string_21.setText("79.3e9");

	string_22.setText("2.27e11");

}

Jacode_iiiAudioProcessorEditor::~Jacode_iiiAudioProcessorEditor()                            //destructor 
{

}

//=================================METHODE======================================
void Jacode_iiiAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
    g.setColour (Colours::darkgrey);
    g.setFont (15.0f);	
	g.drawImage(background,
		0,
		0,
		1600,
		700,
		0,
		0,
		1641,
		700
	);
	processor.drawFrame(g);
}

void Jacode_iiiAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
	auto sliderLeft = 80;
	thresholdSlider.setBounds(sliderLeft, 650, getWidth() - sliderLeft, 40);

	nfftSlider.setBounds(sliderLeft+20,0, getWidth() - sliderLeft-1100, 15);
	betaresSlider.setBounds(sliderLeft+20, 15, getWidth() - sliderLeft-1100, 15);
	nbHarmoF0.setBounds(sliderLeft+20, 30, getWidth() - sliderLeft - 1100, 15);
	nbHarmo.setBounds(sliderLeft+20, 45, getWidth() - sliderLeft - 1100, 15);

	switchPluck.setBounds(sliderLeft+1000, 25, getWidth() - sliderLeft - 1100, 15);
	resPluck.setBounds(sliderLeft + 1050, 0, getWidth() - sliderLeft - 1100, 15);
	distanceMaxPluck.setBounds(sliderLeft + 1050, 15, getWidth() - sliderLeft - 1100, 15);

	reloadButton.setBounds(1280, 35, 1600, 40);

	string_1.setBounds(sliderLeft + 120, 420, sliderLeft + 120 + 40, 15);
	string_2.setBounds(sliderLeft + 160, 420, sliderLeft + 160 + 40, 15);
	string_3.setBounds(sliderLeft + 200, 420, sliderLeft + 200 + 40, 15);
	string_4.setBounds(sliderLeft + 240, 420, sliderLeft + 240 + 40, 15);
	string_5.setBounds(sliderLeft + 280, 420, sliderLeft + 280 + 40, 15);
	string_6.setBounds(sliderLeft + 320, 420, sliderLeft + 320 + 40, 15);

	string_7.setBounds(sliderLeft + 440, 420, sliderLeft + 440 + 40, 15);
	string_8.setBounds(sliderLeft + 480, 420, sliderLeft + 480 + 40, 15);
	string_9.setBounds(sliderLeft + 520, 420, sliderLeft + 520 + 40, 15);
	string_10.setBounds(sliderLeft + 560, 420, sliderLeft + 560 + 40, 15);
	string_11.setBounds(sliderLeft + 600, 420, sliderLeft + 600 + 40, 15);
	string_12.setBounds(sliderLeft + 640, 420, sliderLeft + 640 + 80, 15);

	string_13.setBounds(sliderLeft + 760, 420, sliderLeft + 760 + 40, 15);
	string_14.setBounds(sliderLeft + 800, 420, sliderLeft + 800 + 40, 15);
	string_15.setBounds(sliderLeft + 840, 420, sliderLeft + 840 + 40, 15);
	string_16.setBounds(sliderLeft + 880, 420, sliderLeft + 880 + 40, 15);
	string_17.setBounds(sliderLeft + 920, 420, sliderLeft + 920 + 40, 15);
	string_18.setBounds(sliderLeft + 960, 420, sliderLeft + 960 + 80, 15);

	string_19.setBounds(sliderLeft + 1080, 420, sliderLeft + 1080 + 40, 15);

	string_20.setBounds(sliderLeft + 1200, 420, sliderLeft + 1200 + 40, 15);

	string_21.setBounds(sliderLeft + 1320, 420, sliderLeft + 1320 + 40, 15);

	string_22.setBounds(sliderLeft + 1440, 420, sliderLeft + 1440 + 80, 15);
}

void Jacode_iiiAudioProcessorEditor::timerCallback()
{
	if (processor.getNextBlockReady())
	{
		processor.drawNextFrame();
		repaint();
	}
}

void Jacode_iiiAudioProcessorEditor::sliderValueChanged(Slider* slider)
{
	if (slider == &thresholdSlider)
	{
		thresholdSlider.setValue(thresholdSlider.getValue(), dontSendNotification);
		processor.setThresholdValue(thresholdSlider.getValue());
	}
	else 	if (slider == &betaresSlider)
	{
		betaresSlider.setValue(betaresSlider.getValue(), dontSendNotification);
		processor.setbetaresValue(betaresSlider.getValue());
	}
	else 	if (slider == &nfftSlider)
	{
		nfftSlider.setValue(nfftSlider.getValue(), dontSendNotification);
		processor.setnfftValue(nfftSlider.getValue());
	}
	else 	if (slider == &nbHarmo)
	{
		nbHarmo.setValue(nbHarmo.getValue(), dontSendNotification);
		processor.setNBHarmonics(nbHarmo.getValue());
	}
	else 	if (slider == &resPluck)
	{
		resPluck.setValue(resPluck.getValue(), dontSendNotification);
		processor.setResPluck(resPluck.getValue());
	}
	else 	if (slider == &distanceMaxPluck)
	{
		distanceMaxPluck.setValue(distanceMaxPluck.getValue(), dontSendNotification);
		processor.setDistanceMaxPluck(distanceMaxPluck.getValue());
	}
	else 	if (slider == &nbHarmoF0)
	{
		nbHarmoF0.setValue(nbHarmoF0.getValue(), dontSendNotification);
		processor.setNBHarmonicsF0(nbHarmoF0.getValue());
	}
}

void Jacode_iiiAudioProcessorEditor::buttonClicked(Button* button)
{
	if (button == &switchPluck)                                                     
	{
		processor.setPluckOnOff(switchPluck.getToggleState());
	}
	else 	if (button == &reloadButton)
	{
		nfftSlider.setValue(nfftSlider.getValue()+1.0, dontSendNotification);
		processor.reloadModel(string_1.getText().getDoubleValue(), string_2.getText().getDoubleValue(), string_3.getText().getDoubleValue(), string_4.getText().getDoubleValue(), string_5.getText().getDoubleValue(),
			string_6.getText().getDoubleValue(), string_7.getText().getDoubleValue(), string_8.getText().getDoubleValue(), string_9.getText().getDoubleValue(), string_10.getText().getDoubleValue(),
			string_11.getText().getDoubleValue(), string_12.getText().getDoubleValue(), string_13.getText().getDoubleValue(), string_14.getText().getDoubleValue(), string_15.getText().getDoubleValue(),
			string_16.getText().getDoubleValue(), string_17.getText().getDoubleValue(), string_18.getText().getDoubleValue(), string_19.getText().getDoubleValue(), string_20.getText().getDoubleValue(),
			std::stod(string_21.getText().toStdString()), std::stod(string_22.getText().toStdString()));
	}
}
