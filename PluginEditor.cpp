/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include <string> 

#define BPM 110

//==============================================================================
Jacode_iiiAudioProcessorEditor::Jacode_iiiAudioProcessorEditor (Jacode_iiiAudioProcessor& p) //constructor
    : AudioProcessorEditor (&p), processor (p), value(0)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (1200, 1000);

	labelvalue.setBounds(100, 100, 200, 200);
	labelvalue.setFont(60.0f);
	labelvalue.setText("load...", sendNotification);
	addAndMakeVisible(&labelvalue);
	startTimerHz(60);  //interval en miliseconde divisé par 8 pour faire des 60000 / BPM / 8

	
}

Jacode_iiiAudioProcessorEditor::~Jacode_iiiAudioProcessorEditor()                            //destructor 
{

}

//=================================METHODE======================================
void Jacode_iiiAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
    g.setColour (Colours::white);
    g.setFont (15.0f);
	processor.drawFrame(g);

}

void Jacode_iiiAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
}

void Jacode_iiiAudioProcessorEditor::timerCallback()
{

	if (processor.getNextFFTBlockReady())
	{
		labelvalue.setText(std::to_string(processor.getAfficheValue()), sendNotification);
		processor.drawNextFrameOfSpectrum();
		processor.setNextFFTBlockReady(false);
		repaint();
	}
}