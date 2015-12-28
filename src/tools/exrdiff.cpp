#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <cmath>


#define TINYEXR_IMPLEMENTATION
#include "3rdparty/tinyexr.h"


static bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha);
static void WriteEXR(const char *name, float *pixels, int xRes, int yRes);


static void usage()
{
    fprintf(stderr, "usage: exrdiff [-o difffile.exr] [-d diff tolerance %%] <foo.exr> <bar.exr>\n");
    exit(1);
}

int main(int argc, char *argv[]) 
{
    const char *outfile = NULL;
    const char *imageFile1 = NULL, *imageFile2 = NULL;
    float tol = 0.f;

    if (argc == 1)
        usage();
    
    for (int i = 1; i < argc; ++i)
    {
        if (!strcmp(argv[i], "-o"))
        {
            if (!argv[i+1])
                usage();
            outfile = argv[i+1];
            ++i;
        }
        else if (!strcmp(argv[i], "-d"))
        {
            if (!argv[i+1])
                usage();
            tol = atof(argv[i+1]);
            ++i;
        }
        else if (!imageFile1)
        {
            imageFile1 = argv[i];
        }
        else if (!imageFile2)
        {
            imageFile2 = argv[i];
        }
        else
        {
            usage();
        }
    }

    float *im1, *im2;
    int r1[2], r2[2];
    bool hasAlpha;
    if (!ReadEXR(imageFile1, im1, r1[0], r1[1], hasAlpha))
    {
        printf("couldn't read image %s\n", imageFile1);
        return 1;
    }
    assert(hasAlpha);
    
    if (!ReadEXR(imageFile2, im2, r2[0], r2[1], hasAlpha))
    {
        printf("couldn't read image %s\n", imageFile2);
        return 1;
    }
    assert(hasAlpha);
    
    if (r1[0] != r2[0] || r1[1] != r2[1])
    {
        printf("%s/%s:\n\tresolutions don't match! (%d,%d) vs (%d,%d)\n",
               imageFile1, imageFile2, r1[0], r1[1], r2[0], r2[1]);
        return 1;
    }

    float *diffImage = NULL;
    if (outfile != NULL)
        diffImage = new float[4 * r1[0] * r1[1]];

    double sum1 = 0.f, sum2 = 0.f;
    int smallDiff = 0, bigDiff = 0;
    double mse = 0.f;
    
    for (int i = 0; i < 4*r1[0]*r1[1]; ++i)
    {
        if (diffImage)
            diffImage[i] = fabsf(im1[i] - im2[i]);
        
        if (im1[i] == 0 && im2[i] == 0) 
            continue;
        if ((i % 4) == 3) // alpha channel
            continue;

        sum1 += im1[i];
        sum2 += im2[i];
        float d = fabsf(im1[i] - im2[i]) / im1[i];
        mse += (im1[i] - im2[i]) * (im1[i] - im2[i]);
        
        if (d > .005)
            ++smallDiff;
        if (d > .05)
            ++bigDiff;
    }
    
    double avg1 = sum1 / (3. * r1[0] * r1[1]);
    double avg2 = sum2 / (3. * r1[0] * r1[1]);
    double avgDelta = (avg1-avg2) / std::min(avg1, avg2);
    
    if ((tol == 0. && (bigDiff > 0 || smallDiff > 0)) ||
        (tol > 0. && 100.f * fabs(avgDelta) > tol))
    {
        printf("%s %s\n\tImages differ: %d big (%.2f%%), %d small (%.2f%%)\n"
               "\tavg 1 = %g, avg2 = %g (%f%% delta)\n"
               "\tMSE = %g, RMS = %.3f%%\n",
               imageFile1, imageFile2,
               bigDiff, 100.f * float(bigDiff) / (3 * r1[0] * r1[1]),
               smallDiff, 100.f * float(smallDiff) / (3 * r1[0] * r1[1]),
               avg1, avg2, 100. * avgDelta,
               mse / (3. * r1[0] * r1[1]),
               100. * sqrt(mse / (3. * r1[0] * r1[1])));
        
        if (outfile)
            WriteEXR(outfile, diffImage, r1[0], r1[1]);
        
        return 1;
    }

    return 0;
}

static bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha)
{
    hasAlpha = true;
    const char* err;
    
    int ret = LoadEXR(&rgba, &xRes, &yRes, name, &err);
    if (ret != 0)
    {
        printf("Unable to read image file \"%s\": %s", name, err);
        return false;
    }
    
    printf("Read EXR image %s (%d x %d)", name, xRes, yRes);
    
    return true;
}

static void WriteEXR(const char *name, float *pixels, int xRes, int yRes)
{
    // TODO support scaling, offsetting of output image.
    EXRImage image;
    InitEXRImage(&image);

    int numPixels = xRes * yRes;

    std::vector<float> rBuf, gBuf, bBuf;
    rBuf.resize(numPixels);
    gBuf.resize(numPixels);
    bBuf.resize(numPixels);

    for (int i = 0; i < numPixels; i++) {
        rBuf[i] = pixels[3*i+0];
        gBuf[i] = pixels[3*i+1];
        bBuf[i] = pixels[3*i+2];
    }

    float* imagePtr[3];
    imagePtr[0] = &(bBuf.at(0));
    imagePtr[1] = &(gBuf.at(0));
    imagePtr[2] = &(rBuf.at(0));

    // Must be BGR(A) order, since most of EXR viewers expect this channel order.
    const char* channel_names[] = {"B", "G", "R"}; // "B", "G", "R", "A" for RGBA image
    image.channel_names = channel_names;
    
    image.images = (unsigned char**)imagePtr;
    image.width = xRes;
    image.height = yRes;
    image.compression = TINYEXR_COMPRESSIONTYPE_ZIP;
    image.num_channels = 3;

    image.pixel_types = new int[image.num_channels];
    image.requested_pixel_types = new int[image.num_channels];
    
    for (int i = 0; i < image.num_channels; i++) {
      image.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
      image.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // pixel type of output image to be stored in .EXR
    }

    const char* err;
    int ret = SaveMultiChannelEXRToFile(&image, name, &err);
    if (ret != 0) {
        printf("Unable to write image file \"%s\": %s", name, err);
    }
    printf("Saved EXR image %s (%d x %d)\n", name, image.width, image.height);

    delete[] image.pixel_types;
    delete[] image.requested_pixel_types;
}
