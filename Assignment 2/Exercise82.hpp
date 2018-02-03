#ifndef SUBMISSION_EXERCISE82_HPP_
#define SUBMISSION_EXERCISE82_HPP_

template<typename T>
T CalcAbs(T val) {
	// write a single function to calculate absolute value of an int or double
	T result;
	if (val < 0)
    {
        result = -val;
    }
    else
    {
        result = val;
    }
    return result;
}



#endif /* SUBMISSION_EXERCISE82_HPP_ */
