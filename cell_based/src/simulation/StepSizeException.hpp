class StepSizeException : public std::exception {
	
public:

	double 	displacement;
	double 	suggestedNewStep;
	const std::string message;
	bool endSimulation;

	StepSizeException(double errDisplacement, double errNewStep, const std::string errMessage, bool errEndSim):
		displacement(errDisplacement),
		suggestedNewStep(errNewStep),
		message(errMessage),
		endSimulation(errEndSim),
		std::exception()
	{
	}

	virtual const char* what() const throw() {
		return message.c_str();
	}

	~StepSizeException() throw() {} 
};