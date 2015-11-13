class StepSizeException : public std::exception {
	
public:

	double 		displacement;
	double 		suggestedNewStep;
	const char* message;

	StepSizeException(double errDisplacement, double errNewStep, const char* errMessage):
		displacement(errDisplacement),
		suggestedNewStep(errNewStep),
		message(errMessage),
		std::exception()
	{
	}

	const char* what(){
		return message;
	}
};