#include <string.h>

#define KEWB_FORCE_INLINE	inline __attribute__((__always_inline__))
#define KERNEL_SIZE 3
#define KERNEL_SCALE 9

// according to valgrind the cache miss rates are very low:
// L1 -> 0.01%
// L2 -> 0.02%
// L3 -> 0.9%
typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} pixel;

typedef struct {
	// maybe short is faster
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {
	const int imageScale = n*n;
	const int dim = n;
	const int bytes = imageScale * 3;			// sizeof(pixel) = 3
	/*
	* [1, 1, 1]
	* [1, 1, 1]
	* [1, 1, 1]
	*/
	//int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
	/*
	* [-1, -1, -1]
	* [-1, 9, -1]
	* [-1, -1, -1]
	*/
	//int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};
	// pointer to last pixel == Image[N][N]
	pixel* limit = (pixel *)image->data + imageScale - 1;
	pixel* dst = (pixel *)image->data;	// applying convolution directly on the input
	pixel_sum pixel_sum;				// accumlator
	pixel* endLine;
	register pixel* iter1;						// pointer to image[i-1][j-1]
	register pixel* iter2;						// pointer to image[i][j-1]
	register pixel* iter3;						// pointer to image[i+1][j-1]
	pixel* src = malloc(bytes);			// allocate space for source image
	pixel* firstDst = dst + dim + 1;
	pixel* firstSrc = src + dim + 1;
	pixel* temp = src;
	if (flag == '1') {
		/***BLUR CONVOLUTION***/
		memcpy(src, dst, bytes);		// src = image.copy()
		// src --> src[kernelsize/2][kernelsize/2] = src[1][1] = src[1*dim + 1] = src[dim+1]
		dst = firstDst;
		src = firstSrc; 				// point to first valid pixel image[1][1]
		endLine = dst + dim - 2;		// the address of the last pixel in the current output line
		iter1 = temp;					// pointer to image[i-1][j-1]
		iter2 = src - 1;				// pointer to image[i][j-1]
		iter3 = iter2 + dim;			// pointer to image[i+1][j-1]
		do {
			do {
				/***Applying the kernel***/
				// todo -> apply SIMD
				// all reds:
				pixel_sum.red = iter1->red + (iter1 + 1)->red + (iter1 + 2)->red
				+ iter2->red + (iter2 + 1)->red + (iter2 + 2)->red
				+ iter3->red + (iter3 + 1)->red + (iter3 + 2)->red;
				// all greens
				pixel_sum.green = iter1->green + (iter1 + 1)->green + (iter1 + 2)->green 
				+ iter2->green + (iter2 + 1)->green + (iter2 + 2)->green 
				+ iter3->green + (iter3 + 1)->green + (iter3 + 2)->green;
				// all blues
				pixel_sum.blue = iter1->blue + (iter1 + 1)->blue + (iter1 + 2)->blue 
				+ iter2->blue + (iter2 + 1)->blue + (iter2 + 2)->blue
				+ iter3->blue + (iter3 + 1)->blue + (iter3 + 2)->blue;
				/**Assign sum to pixel at [i,j]**/
				/* 32 bits signed integer i, can be clamped into 8 bits unsigned char
				with this branchless code:
				clampIntToRange(int i, 0, 255):
					i &= -(!(i >> 31));
					i |= ((255 - i) >> 31); */
				// clamp values to match the range [0,255]
				pixel_sum.red = pixel_sum.red / 9; 						// divide by kernel's weight, todo optimize
				pixel_sum.green = pixel_sum.green / 9; 					// divide by kernel's weight, todo optimize
				pixel_sum.blue = pixel_sum.blue / 9; 					// divide by kernel's weight, todo optimize
				pixel_sum.red &= -(!(pixel_sum.red >> 31));
				pixel_sum.green &= -(!(pixel_sum.green >> 31));
				pixel_sum.blue &= -(!(pixel_sum.blue >> 31));
				// copy answers to pixel
				dst->red = pixel_sum.red | ((255 - pixel_sum.red) >> 31);
				dst->green = pixel_sum.green | ((255 - pixel_sum.green) >> 31);
				dst->blue = pixel_sum.blue | ((255 - pixel_sum.blue) >> 31);
				/***Moving to next pixel***/
				++src;
				++dst;
				// adjusting the kernel pixels iterators
				++iter1;
				++iter2;
				++iter3;
				// pixel sums iterator
			} while (dst != endLine);
			endLine += dim;					// calculate next line's end address
			// move src and dst to next rows, add kernelsize which is 3 in our case
			// but we add kernelsize - 1 because inner loop already added 1
			src += 2;
			dst += 2;
			// adjusting the kernel pixels iterators
			iter1 += 2;
			iter2 += 2;
			iter3 += 2;
		} while (endLine < limit);
		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);

		/***SHARP CONVOLUTION***/
		memcpy(temp, image->data, bytes);	// src = image.copy()
		src = firstSrc; 					// point to first valid pixel image[1][1]
		dst = firstDst;
		endLine = dst + dim - 2;			// the address of the last pixel in the current output line
		iter1 = temp;						// pointer to image[i-1][j-1]
		iter2 = src - 1;					// pointer to image[i][j-1]
		iter3 = iter2 + dim;				// pointer to image[i+1][j-1]
		do {
			do {
				/***Applying the kernel***/
				// todo -> apply SIMD
				// all reds:
				pixel_sum.red = -(iter1->red + (iter1 + 1)->red + (iter1 + 2)->red
				+ iter2->red + (iter2 + 2)->red
				+ iter3->red + (iter3 + 1)->red + (iter3 + 2)->red) + (9*((iter2 + 1)->red));
				// all greens
				pixel_sum.green = -(iter1->green + (iter1 + 1)->green + (iter1 + 2)->green
				+ iter2->green + (iter2 + 2)->green
				+ iter3->green + (iter3 + 1)->green + (iter3 + 2)->green) + (9*((iter2 + 1)->green));
				// all blues
				pixel_sum.blue = -(iter1->blue + (iter1 + 1)->blue + (iter1 + 2)->blue
				+ iter2->blue + (iter2 + 2)->blue
				+ iter3->blue + (iter3 + 1)->blue + (iter3 + 2)->blue) + (9*((iter2 + 1)->blue));

				/**Assign sum to pixel at [i,j]**/
				/* 32 bits signed integer i, can be clamped into 8 bits unsigned char
				with this branchless code:
				clampIntToRange(int i, 0, 255):
					i &= -(!(i >> 31));
					i |= ((255 - i) >> 31); */
				// clamp values to match the range [0,255]
				pixel_sum.red &= -(!(pixel_sum.red >> 31));
				pixel_sum.green &= -(!(pixel_sum.green >> 31));
				pixel_sum.blue &= -(!(pixel_sum.blue >> 31));
				// copy answers to pixel
				dst->red = pixel_sum.red | ((255 - pixel_sum.red) >> 31);
				dst->green = pixel_sum.green | ((255 - pixel_sum.green) >> 31);
				dst->blue = pixel_sum.blue | ((255 - pixel_sum.blue) >> 31);
				/***Moving to next pixel***/
				++src;
				++dst;
				// adjusting the kernel pixels iterators
				++iter1;
				++iter2;
				++iter3;
			} while (dst != endLine);
			endLine += dim;					// calculate next line's end address
			// move src and dst to next rows, add kernelsize which is 3 in our case
			// but we add kernelsize - 1 because inner loop already added 1
			src += 2;
			dst += 2;
			// adjusting the kernel pixels iterators
			iter1 += 2;
			iter2 += 2;
			iter3 += 2;

		} while (endLine < limit);
		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);
	} else {
		/***Calculating the sum of each pixel in advance***/
		// allocate 4 bytes integer for each pixel
		int* pixel_sums = malloc(imageScale << 2);
		pixel* pixel_iter = dst;					// pixels iterator
		pixel* last = limit;						// last pixel address + 1
		int* sum_iter = pixel_sums;					// pixel_sums iterator
		for ( ;pixel_iter < last; pixel_iter += 8, sum_iter += 8) {
			// pixel_sums[i] = sum_of_pixel(image_data[i]);
			// write the sum of pixel at pixel_iter to sum_iter
			*sum_iter = ((int) pixel_iter->red) + ((int) pixel_iter->green) + ((int) pixel_iter->blue);
			*(sum_iter  + 1) = ((int) (pixel_iter + 1)->red) + ((int) (pixel_iter + 1)->green) + ((int) (pixel_iter + 1)->blue);
			*(sum_iter  + 2) = ((int) (pixel_iter + 2)->red) + ((int) (pixel_iter + 2)->green) + ((int) (pixel_iter + 2)->blue);
			*(sum_iter  + 3) = ((int) (pixel_iter + 3)->red) + ((int) (pixel_iter + 3)->green) + ((int) (pixel_iter + 3)->blue);
			*(sum_iter  + 4) = ((int) (pixel_iter + 4)->red) + ((int) (pixel_iter + 4)->green) + ((int) (pixel_iter + 4)->blue);
			*(sum_iter  + 5) = ((int) (pixel_iter + 5)->red) + ((int) (pixel_iter + 5)->green) + ((int) (pixel_iter + 5)->blue);
			*(sum_iter  + 6) = ((int) (pixel_iter + 6)->red) + ((int) (pixel_iter + 6)->green) + ((int) (pixel_iter + 6)->blue);
			*(sum_iter  + 7) = ((int) (pixel_iter + 7)->red) + ((int) (pixel_iter + 7)->green) + ((int) (pixel_iter + 7)->blue);
		}
		for ( ;pixel_iter < last; ++pixel_iter, ++sum_iter) {
			// pixel_sums[i] = sum_of_pixel(image_data[i]);
			// write the sum of pixel at pixel_iter to sum_iter
			*sum_iter = ((int) pixel_iter->red) + ((int) pixel_iter->green) + ((int) pixel_iter->blue);
		}
		int* todel = pixel_sums;

		/***FILTERED BLUR***/
		memcpy(src, dst, bytes);			// src = image.copy()
		src = firstSrc; 					// point to first valid pixel image[1][1]
		dst = firstDst;
		pixel* endLine = dst + dim - 2;		// the address of the last pixel in the current output line
		iter1 = temp;				// pointer to image[i-1][j-1]
		iter2 = src - 1;				// pointer to image[i][j-1]
		iter3 = iter2 + dim;			// pointer to image[i+1][j-1]
		// pointer to sumArray[i-1][j-1]
		int newlineOffset = dim - 2;
		register pixel* minPos = iter1;
		register pixel* maxPos = iter1;
		register int min_value;
		register int max_value;
		register int sum;	// accumlator
		int *sums;			// pixel sums
		do {
			do {
				/***Applying the kernel***/
				// todo -> apply SIMD
				// all reds:
				pixel_sum.red = iter1->red + (iter1 + 1)->red + (iter1 + 2)->red
				+ iter2->red + (iter2 + 1)->red + (iter2 + 2)->red
				+ iter3->red + (iter3 + 1)->red + (iter3 + 2)->red;
				// all greens
				pixel_sum.green = iter1->green + (iter1 + 1)->green + (iter1 + 2)->green 
				+ iter2->green + (iter2 + 1)->green + (iter2 + 2)->green 
				+ iter3->green + (iter3 + 1)->green + (iter3 + 2)->green;
				// all blues
				pixel_sum.blue = iter1->blue + (iter1 + 1)->blue + (iter1 + 2)->blue 
				+ iter2->blue + (iter2 + 1)->blue + (iter2 + 2)->blue 
				+ iter3->blue + (iter3 + 1)->blue + (iter3 + 2)->blue;
				/***Filter***/
				// for each pixel in kernel -> calculate it's sum
				sums = pixel_sums;	// pointer to sums array
				minPos = iter1;
				maxPos = iter1;
				// buttom left pixel, iter1 points to buttom left pixel::
				sum = *sums;
				// no branch needed here, because the first sum is also the current minimum and maximum
				min_value = sum;
				max_value = sum;

				// buttom middle pixel:
				++iter1;				// advance pointer
				++sums;
				sum = *sums;
				if (min_value < sum) {
					// if sum > min_value, and we know that
					// 		min_value == max_value
					// we conclude that ==> sum >= max_value
					if (max_value != sum) {
						max_value = sum;
						++maxPos;			// maxPos points to iter1
					}
				} else {
					// sum <= min_value
					min_value = sum;
					++minPos; 				// minPos points to iter1
				}

				// buttom right pixel:
				++iter1;
				++sums;
				sum = *sums;
				if (min_value < sum) {
					// check if sum is new maximun
					if (sum > max_value) {
						// set sum as new maximum
						max_value = sum;
						maxPos = iter1;		// minPos points to current pixel
					}
				} else {
					// set sum as new minimum
					min_value = sum;
					minPos = iter1;			// maxPos points to current pixel
				}

				// center left pixel, iter2 points to that place
				sums += newlineOffset;
				sum = *sums;
				if (min_value < sum) {
					// check if sum is new maximun
					if (sum <= max_value) {} else {
						// set sum as new maximum
						max_value = sum;
						maxPos = iter2;
					}
				} else {
					// set sum as new minimum
					min_value = sum;
					minPos = iter2;
				}

				// center middle pixel
				++iter2;
				++sums;
				sum = *sums;
				if (min_value < sum) {
					// check if sum is new maximun
					if (sum <= max_value) {} else {
						// set sum as new maximum
						max_value = sum;
						maxPos = iter2;
					}
				} else {
					// set sum as new minimum
					min_value = sum;
					minPos = iter2;
				}

				// center right pixel:
				++iter2;
				++sums;
				sum = *sums;
				if (min_value < sum) {
					// check if sum is new maximun
					if (sum <= max_value) {} else {
						// set sum as new maximum
						max_value = sum;
						maxPos = iter2;
					}
				} else {
					// set sum as new minimum
					min_value = sum;
					minPos = iter2;
				}

				// top left pixel, iter3 points to that pixel
				sums += newlineOffset;
				sum = *sums;
				if (min_value < sum) {
					// check if sum is new maximun
					if (sum <= max_value) {} else {
						// set sum as new max
						// set sum as new minimum
						max_value = sum;
						maxPos = iter3;
					}
				} else {
					// set sum as new minimum
					min_value = sum;
					minPos = iter3;
				}

				// top middle pixel:
				++iter3;
				++sums;
				sum = *sums;
				if (min_value < sum) {
					// check if sum is new maximun
					if (sum <= max_value) {} else {
						// set sum as new maximum
						max_value = sum;
						maxPos = iter3;
					}
				} else {
					// set sum as new minimum
					min_value = sum;
					minPos = iter3;
				}

				// top right pixel:
				++iter3;
				++sums;
				sum = *sums;
				if (min_value < sum) {
					// check if sum is new maximun
					if (sum <= max_value) {} else {
						// set sum as new maximum
						max_value = sum;
						maxPos = iter3;
					}
				} else {
					// set sum as new minimum
					min_value = sum;
					minPos = iter3;
				}

				// filter out the neibour pixels with min and max sums
				pixel pMin = *minPos;
				pixel pMax = *maxPos;
				iter1 -= 2;						// adjusting the iterators because we moved them 
				pixel_sum.red -= pMin.red + pMax.red;
				iter2 -= 2;						// adjusting the iterators because we moved them 
				pixel_sum.green -= pMin.green + pMax.green;
				iter3 -= 2;						// adjusting the iterators because we moved them 
				pixel_sum.blue -= pMin.blue + pMax.blue;
				/**Assign sum to pixel at [i,j]**/
				/* 32 bits signed integer i, can be clamped into 8 bits unsigned char
				with this branchless code:
				clampIntToRange(int i, 0, 255):
					i &= -(!(i >> 31));
					i |= ((255 - i) >> 31); */
				// clamp values to match the range [0,255]
				pixel_sum.red = pixel_sum.red / 7; 					// divide by kernel's weight, todo optimize
				pixel_sum.green = pixel_sum.green / 7; 				// divide by kernel's weight, todo optimize
				pixel_sum.blue = pixel_sum.blue / 7; 					// divide by kernel's weight, todo optimize
				pixel_sum.red &= -(!(pixel_sum.red >> 31));
				pixel_sum.green &= -(!(pixel_sum.green >> 31));
				pixel_sum.blue &= -(!(pixel_sum.blue >> 31));
				// copy answers to pixel
				dst->red = pixel_sum.red | ((255 - pixel_sum.red) >> 31);
				dst->green = pixel_sum.green | ((255 - pixel_sum.green) >> 31);
				dst->blue = pixel_sum.blue | ((255 - pixel_sum.blue) >> 31);
				/***Moving to next pixel***/
				++src;
				++dst;
				// adjusting the kernel pixels iterators
				++iter1;
				++iter2;
				++iter3;
				// pixel sums iterator
				++pixel_sums;
			} while (dst != endLine);
			endLine += dim;					// calculate next line's end address
			// move src and dst to next rows, add kernelsize which is 3 in our case
			// but we add kernelsize - 1 because inner loop already added 1
			src += 2;
			dst += 2;
			// adjusting the kernel pixels iterators
			iter1 += 2;
			iter2 += 2;
			iter3 += 2;
			pixel_sums += 2;
		} while (endLine < limit);
		free(todel);
		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);
		/***SHARP CONVOLUTION***/
		memcpy(temp, image->data, bytes);	// src = image.copy()
		src = firstSrc; 					// point to first valid pixel image[1][1]
		dst = firstDst;
		endLine = dst + dim - 2;			// the address of the last pixel in the current output line
		iter1 = temp;						// pointer to image[i-1][j-1]
		iter2 = src - 1;					// pointer to image[i][j-1]
		iter3 = iter2 + dim;				// pointer to image[i+1][j-1]
		do {
			do {
				/***Applying the kernel***/
				// todo -> apply SIMD
				// all reds:
				pixel_sum.red = -(iter1->red + (iter1 + 1)->red + (iter1 + 2)->red
				+ iter2->red + (iter2 + 2)->red
				+ iter3->red + (iter3 + 1)->red + (iter3 + 2)->red) + (9*((iter2 + 1)->red));
				// all greens
				pixel_sum.green = -(iter1->green + (iter1 + 1)->green + (iter1 + 2)->green
				+ iter2->green + (iter2 + 2)->green
				+ iter3->green + (iter3 + 1)->green + (iter3 + 2)->green) + (9*((iter2 + 1)->green));
				// all blues
				pixel_sum.blue = -(iter1->blue + (iter1 + 1)->blue + (iter1 + 2)->blue
				+ iter2->blue + (iter2 + 2)->blue
				+ iter3->blue + (iter3 + 1)->blue + (iter3 + 2)->blue) + (9*((iter2 + 1)->blue));

				/**Assign sum to pixel at [i,j]**/
				/* 32 bits signed integer i, can be clamped into 8 bits unsigned char
				with this branchless code:
				clampIntToRange(int i, 0, 255):
					i &= -(!(i >> 31));
					i |= ((255 - i) >> 31); */
				// clamp values to match the range [0,255]
				pixel_sum.red &= -(!(pixel_sum.red >> 31));
				pixel_sum.green &= -(!(pixel_sum.green >> 31));
				pixel_sum.blue &= -(!(pixel_sum.blue >> 31));
				// copy answers to pixel
				dst->red = pixel_sum.red | ((255 - pixel_sum.red) >> 31);
				dst->green = pixel_sum.green | ((255 - pixel_sum.green) >> 31);
				dst->blue = pixel_sum.blue | ((255 - pixel_sum.blue) >> 31);
				/***Moving to next pixel***/
				++src;
				++dst;
				// adjusting the kernel pixels iterators
				++iter1;
				++iter2;
				++iter3;
			} while (dst != endLine);
			endLine += dim;					// calculate next line's end address
			// move src and dst to next rows, add kernelsize which is 3 in our case
			// but we add kernelsize - 1 because inner loop already added 1
			src += 2;
			dst += 2;
			// adjusting the kernel pixels iterators
			iter1 += 2;
			iter2 += 2;
			iter3 += 2;

		} while (endLine < limit);
		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);	
	}
	//***free allocated memory***//
	free(temp);
}

/*
int* calculate_sums(pixel* imageData, pixel* limit, int imageScale) {
	// allocate 4 bytes integer for each pixel
	int* pixel_sums = malloc(imageScale << 2);
	pixel* pixel_iter = imageData;				// pixels iterator
	pixel* last = limit;						// last pixel address + 1
	int* sum_iter = pixel_sums;					// pixel_sums iterator
	for ( ;pixel_iter < last; pixel_iter += 8, sum_iter += 8) {
		// pixel_sums[i] = sum_of_pixel(image_data[i]);
		// write the sum of pixel at pixel_iter to sum_iter
		*sum_iter = ((int) pixel_iter->red) + ((int) pixel_iter->green) + ((int) pixel_iter->blue);
		*(sum_iter  + 1) = ((int) (pixel_iter + 1)->red) + ((int) (pixel_iter + 1)->green) + ((int) (pixel_iter + 1)->blue);
		*(sum_iter  + 2) = ((int) (pixel_iter + 2)->red) + ((int) (pixel_iter + 2)->green) + ((int) (pixel_iter + 2)->blue);
		*(sum_iter  + 3) = ((int) (pixel_iter + 3)->red) + ((int) (pixel_iter + 3)->green) + ((int) (pixel_iter + 3)->blue);
		*(sum_iter  + 4) = ((int) (pixel_iter + 4)->red) + ((int) (pixel_iter + 4)->green) + ((int) (pixel_iter + 4)->blue);
		*(sum_iter  + 5) = ((int) (pixel_iter + 5)->red) + ((int) (pixel_iter + 5)->green) + ((int) (pixel_iter + 5)->blue);
		*(sum_iter  + 6) = ((int) (pixel_iter + 6)->red) + ((int) (pixel_iter + 6)->green) + ((int) (pixel_iter + 6)->blue);
		*(sum_iter  + 7) = ((int) (pixel_iter + 7)->red) + ((int) (pixel_iter + 7)->green) + ((int) (pixel_iter + 7)->blue);
	}
	for ( ;pixel_iter < last; ++pixel_iter, ++sum_iter) {
		// pixel_sums[i] = sum_of_pixel(image_data[i]);
		// write the sum of pixel at pixel_iter to sum_iter
		*sum_iter = ((int) pixel_iter->red) + ((int) pixel_iter->green) + ((int) pixel_iter->blue);
	}
	return pixel_sums;
}
*/


/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
/*
void smooth(int dim, pixel *src, pixel *dst, int *kernel, int kernelScale, bool filter, pixel* limit, int* pixel_sums) {
	// src -> src[kernelsize/2][kernelsize/2] = src[1][1] = src[1*dim + 1] = src[dim+1]
	src += dim + 1; // point to first valid pixel
	dst += dim + 1;
	pixel* endLine = dst + dim - 2;
	do {
		do {
			applyKernel(dim, src, dst, kernel, kernelScale, filter, pixel_sums);
			++src;
			++dst;
		} while (dst != endLine);
		// calculate next line's limit
		endLine += dim;
		// move src and dst to next rows, add kernelsize which is 3 in our case
		// but we add kernelsize - 1 because inner loop already added 1
		src += 2;
		dst += 2;
	} while (endLine < limit);
}
*/
/*
void charsToPixels(Image *charsImg, pixel* pixels, int bytes) {
	// each pixel is 3 bytes, so we have len*3 bytes to copy
	// the compiler will calculate len*3 efficently
	memcpy(pixels, image->data, bytes);
}

void pixelsToChars(pixel* pixels, Image *charsImg, int bytes) {
	// each pixel is 3 bytes, so we have len*3 bytes to copy
	// the compiler will calculate len*3 efficently
	memcpy(image->data, pixels, bytes);
}
// todo copy memcpy() soruce to here
void copyPixels(pixel* src, pixel* dst, int bytes) {
	// each pixel is 3 bytes, so we have len*3 bytes to copy
	// the compiler will calculate len*3 efficently
	memcpy(dst, src, bytes);
}
*/
//smooth(image->sizeX, src, dst, (int *) kernel, kernelScale, filter, limit, pixel_sums);
/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
/*
inline static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {
	/* 32 bits signed integer i, can be clamped into 8 bits unsigned char
		with this branchless code:
	 clampIntToRange(int i, 0, 255):
	 	i &= -(!(i >> 31));
     	i |= ((255 - i) >> 31);
	*/
	// clamp values to match the range [0,255]
	/*
	sum.red = sum.red / kernelScale; 		// divide by kernel's weight, todo optimize
	sum.green = sum.green / kernelScale; 	// divide by kernel's weight, todo optimize
	sum.blue = sum.blue / kernelScale; 		// divide by kernel's weight, todo optimize
	sum.red &= -(!(sum.red >> 31));
	sum.green &= -(!(sum.green >> 31));
	sum.blue &= -(!(sum.blue >> 31));
	// copy answers to pixel
	current_pixel->red = sum.red | ((255 - sum.red) >> 31);
	current_pixel->green = sum.green | ((255 - sum.green) >> 31);
	current_pixel->blue = sum.blue | ((255 - sum.blue) >> 31);
	return;
}
*/
/*
 *  Applies kernel for pixel at (i,j)
 */
/*
inline static void applyKernel(int dim, pixel *src, pixel* dst, int *kernel, int kernelScale, bool filter, int* pixel_sums) {
	pixel_sum pixel_sum;
	pixel* iter1 = src - dim - 1;	// pointer to image[i-1][j-1]
	pixel* iter2 = src - 1;			// pointer to image[i][j-1]
	pixel* iter3 = iter2 + dim;		// pointer to image[i+1][j-1]
	// todo -> optimize these defines
	int w0 = kernel[0];
	int w1 = kernel[1];
	int w2 = kernel[2];
	int w3 = kernel[3];
	int w4 = kernel[4];
	int w5 = kernel[5];
	int w6 = kernel[6];
	int w7 = kernel[7];
	int w8 = kernel[8];
	// todo -> apply SIMD
	// all reds:
	pixel_sum.red = w0*iter1->red + w1*(iter1 + 1)->red + w2*(iter1 + 2)->red
	+ w3*iter2->red + w4*(iter2 + 1)->red + w5*(iter2 + 2)->red
	+ w6*iter3->red + w7*(iter3 + 1)->red + w8*(iter3 + 2)->red;
	// all greens
	pixel_sum.green = w0*iter1->green + w1*(iter1 + 1)->green + w2*(iter1 + 2)->green 
	+ w3*iter2->green + w4*(iter2 + 1)->green + w5*(iter2 + 2)->green 
	+ w6*iter3->green + w7*(iter3 + 1)->green + w8*(iter3 + 2)->green;
	// all blues
	pixel_sum.blue = w0*iter1->blue + w1*(iter1 + 1)->blue + w2*(iter1 + 2)->blue 
	+ w3*iter2->blue + w4*(iter2 + 1)->blue + w5*(iter2 + 2)->blue 
	+ w6*iter3->blue + w7*(iter3 + 1)->blue + w8*(iter3 + 2)->blue;
	
	if (filter) {
		// todo: cache sums
		// for each pixel in kernel -> calculate it's sum
		register pixel* minPos = iter1;
		register pixel* maxPos = iter1;
		register int min_value;
		register int max_value;
		register int sum;

		// the offset of the top left pixel
		//int delta = dst - (pixel *)image->data - dim - 1;
		int newlineOffset = dim - 2;
		int* sum_array = pixel_sums + (dst - (pixel *)image->data - dim - 1);

		// buttom left pixel, iter1 points to buttom left pixel::
		sum = *sum_array;
		// no branch needed here, because the first sum is also the current minimum and maximum
		min_value = sum;
		max_value = sum;

		// buttom middle pixel:
		++iter1;				// advance pointer
		++sum_array;
		sum = *sum_array;
		if (min_value < sum) {
			// if sum > min_value, and we know that
			// 		min_value == max_value
			// we conclude that ==> sum >= max_value
			if (max_value != sum) {
				max_value = sum;
				++maxPos;			// maxPos points to iter1
			}
		} else {
			// sum <= min_value
			min_value = sum;
			++minPos; 				// minPos points to iter1
		}

		// buttom right pixel:
		++iter1;
		++sum_array;
		sum = *sum_array;
		if (min_value < sum) {
			// check if sum is new maximun
			if (sum > max_value) {
				// set sum as new maximum
				max_value = sum;
				maxPos = iter1;		// minPos points to current pixel
			}
		} else {
			// set sum as new minimum
			min_value = sum;
			minPos = iter1;			// maxPos points to current pixel
		}

		// center left pixel, iter2 points to that place
		sum_array += newlineOffset;
		sum = *sum_array;
		if (min_value < sum) {
			// check if sum is new maximun
			if (sum <= max_value) {} else {
				// set sum as new maximum
				max_value = sum;
				maxPos = iter2;
			}
		} else {
			// set sum as new minimum
			min_value = sum;
			minPos = iter2;
		}

		// center middle pixel
		++iter2;
		++sum_array;
		sum = *sum_array;
		if (min_value < sum) {
			// check if sum is new maximun
			if (sum <= max_value) {} else {
				// set sum as new maximum
				max_value = sum;
				maxPos = iter2;
			}
		} else {
			// set sum as new minimum
			min_value = sum;
			minPos = iter2;
		}

		// center right pixel:
		++iter2;
		++sum_array;
		sum = *sum_array;
		if (min_value < sum) {
			// check if sum is new maximun
			if (sum <= max_value) {} else {
				// set sum as new maximum
				max_value = sum;
				maxPos = iter2;
			}
		} else {
			// set sum as new minimum
			min_value = sum;
			minPos = iter2;
		}

		// top left pixel, iter3 points to that pixel
		sum_array += newlineOffset;
		sum = *sum_array;
		if (min_value < sum) {
			// check if sum is new maximun
			if (sum <= max_value) {} else {
				// set sum as new max
			// set sum as new minimum
			max_value = sum;
			maxPos = iter3;
			}
		} else {
			// set sum as new minimum
			min_value = sum;
			minPos = iter3;
		}

		// top middle pixel:
		++iter3;
		++sum_array;
		sum = *sum_array;
		if (min_value < sum) {
			// check if sum is new maximun
			if (sum <= max_value) {} else {
				// set sum as new maximum
				max_value = sum;
				maxPos = iter3;
			}
		} else {
			// set sum as new minimum
			min_value = sum;
			minPos = iter3;
		}

		// top right pixel:
		++iter3;
		++sum_array;
		sum = *sum_array;
		if (min_value < sum) {
			// check if sum is new maximun
			if (sum <= max_value) {} else {
				// set sum as new maximum
				max_value = sum;
				maxPos = iter3;
			}
		} else {
			// set sum as new minimum
			min_value = sum;
			minPos = iter3;
		}

		// filter out the neibour pixels with min and max sums
		pixel pMin = *minPos;
		pixel pMax = *maxPos;
		pixel_sum.red -= pMin.red + pMax.red;
		pixel_sum.green -= pMin.green + pMax.green;
		pixel_sum.blue -= pMin.blue + pMax.blue;
	}
	// assign kernel's result to pixel at [i,j]
	assign_sum_to_pixel(dst, pixel_sum, kernelScale);
}
*/