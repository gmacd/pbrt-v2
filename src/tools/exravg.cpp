#include <stdio.h>
#include <stdlib.h>


#define TINYEXR_IMPLEMENTATION
#include "3rdparty/tinyexr.h"


static bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha);

int main(int argc, char *argv[]) 
{
    if (argc < 2) {
        fprintf(stderr, "usage: exravg [file1.exr] <file2.exr> ...\n");
        return 1;
    }

    float *rgba, *orig_rgba = NULL;
    int xRes = 0, yRes = 0;
    bool hasAlpha = false;
    float a = 0;
    int file;

    for (file = 1 ; file < argc ; file++)
    {
        if (ReadEXR(argv[file], rgba, xRes, yRes, hasAlpha))
        {
            orig_rgba = rgba;
            a = 0;
            for (int i = 0; i < xRes*yRes; ++i)
            {
                for (int j = 0; j < 3; ++j)
                    a += rgba[j];
                rgba += hasAlpha ? 4 : 3;
            }
        }
        
        printf("%s: Average value %f\n", argv[file], a / (3.f * xRes * yRes));
        delete[] orig_rgba;
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
