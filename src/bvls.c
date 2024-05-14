//
//    This source file is part of the software to infer antigenic trees.
//    Copyright (C) 2012  Lars Steinbrueck
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <jni.h>
#include <stdio.h>
#include <stdlib.h>

// pointer to FORTRAN90 subroutine
extern void bvls_ ();
//extern void bvls(float A[][2], float B[] , float BND[][2], float X[], float *RNORM, int *NSETP, float W[], int INDEX[], int *IERR, int *DIMN, int *DIMM);

JNIEXPORT jfloatArray JNICALL Java_NNLSsolver_solveBVLS (JNIEnv * env, jobject jobj, jfloatArray A_java, jfloatArray b_java, jfloatArray bnd_java, jint dim_n, jint dim_m) {
	// variable declaration
	int	i, j, nstp, ierr;
	float	rnorm;
	// local pointers to arrays
	float	*A_tmp = (*env)->GetFloatArrayElements (env, A_java, NULL);
	float	*b_tmp = (*env)->GetFloatArrayElements (env, b_java, NULL);
	float	*bnd_tmp = (*env)->GetFloatArrayElements (env, bnd_java, NULL);
	float	**A;
	float	**bnd;
	float	*x = (float*) malloc (sizeof(float) * dim_m);
	// return vector
	jfloatArray	resultVec = (*env)->NewFloatArray(env, dim_m);
	float	*W = (float*) malloc (sizeof(float) * dim_n);
	float	*index = (float*) malloc (sizeof(float) * dim_n);
	
	// compute BVLS
	bvls_ (A_tmp, b_tmp, bnd_tmp, x, &rnorm, &nstp, W, index, &ierr, &dim_n, &dim_m);
	
	// copy results to java readable variable
	(*env)->SetFloatArrayRegion(env, resultVec, 0, dim_m, x);
	
	// free variables
	(*env)->ReleaseFloatArrayElements (env, A_java, A_tmp, 0);
	(*env)->ReleaseFloatArrayElements (env, b_java, b_tmp, 0);
	(*env)->ReleaseFloatArrayElements (env, bnd_java, bnd_tmp, 0);
	
	free (x);
	free (W);
	free (index);
	
	return (resultVec);
}
