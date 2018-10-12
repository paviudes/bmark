/*
 * weights.h
 *
 *  Created on: Aug 22, 2015
 *      Author: ndelfosse
 */

#ifndef BENCHMARKING_C_WEIGHTS_H_
#define BENCHMARKING_C_WEIGHTS_H_

extern int MaxWeight(struct tiling G, struct tiling H);
extern int *FaceMeasurementWeights(struct tiling G, int maxWeight);

#endif /* BENCHMARKING_C_WEIGHTS_H_ */
