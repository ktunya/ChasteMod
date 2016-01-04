#ifndef _STEPPPERCHOICE_HPP_
#define _STEPPPERCHOICE_HPP_

class StepperChoice {
	
	public:

		/** Enum representing available timesteppers **/
    	enum Steppers {EULER, RK4, BACKWARDEULER, ADAMSMOULTON, DOP853, LAST}; 
};

#endif /*_STEPPPERCHOICE_HPP_*/