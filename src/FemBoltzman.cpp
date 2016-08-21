
#include <iostream>
#include <armadillo>
#include "FiniteElemNL.h"
#include "FemBoltzman.h"

using namespace std;
using namespace arma;

/**
 * This is a Java JNI interface file. It runs calls the libraries present
 * and passes the parameters.
 *
 *
 * input double vector has 5 fields:
 *
 * h - index 0.
 * x0 - index 1.
 * x1 - index 2.
 * y0 - index 3.
 * y1 - index 4.
 *
 *
 * @param *env - Java environment (provided)
 * @param jobj - Java object (provided)
 * @param jarray - Java double array to extract parameters from
 * @param femeqn - String representing the name of the FEM equation to solve.
 * @return returns a double array representing the solution to the PDE.
 */
JNIEXPORT jdoubleArray JNICALL Java_net_tedkwan_javafemjni_FemBoltzman_cFEMnLin(
		JNIEnv *env, jobject jobj, jdoubleArray jarray, jstring femeqn, jint levels) {

	// Extract Java inputs and convert to C++ objects.
	jboolean isCopy1;
	jboolean isCopy2;
	const char *s = env->GetStringUTFChars(femeqn, &isCopy2);
	int J=levels;
	string femeqnc(s);
	jdouble* srcArrayElems = env->GetDoubleArrayElements(jarray, &isCopy1);

	// Set parameters.
	vector<double> pars(5);
	pars[0] = srcArrayElems[0];
	pars[1] = srcArrayElems[1];
	pars[2] = srcArrayElems[2];
	pars[3] = srcArrayElems[3];
	pars[4] = srcArrayElems[4];

	// Run finite element method.
	FiniteElemNL femp(pars, femeqnc,J);

	// Release memory.
	if (isCopy2 == JNI_TRUE) {
		env->ReleaseStringUTFChars(femeqn, s);
	}

	// Get output information.
	vec u = femp.u;
	mat node=femp.mesh.node;

	// Setup array elements for output.
	vector<double> res;
	vec outvec=join_vert(node.col(0),join_vert(node.col(1),u));
	for (uword j = 0; j < outvec.n_elem; j++) {
		// Add all values to output vector.
		res.push_back(outvec(j));
	}

	// Release memory of input array.
	if (isCopy1 == JNI_TRUE) {
		env->ReleaseDoubleArrayElements(jarray, srcArrayElems, JNI_ABORT);
	}
	// Create output array.
	int lenres = res.size();
	jdoubleArray result = env->NewDoubleArray(lenres);
	env->SetDoubleArrayRegion(result, 0, lenres, &res[0]);

	return result;
}
