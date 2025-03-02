#include "stdio.h"
#include "stdlib.h"

void read_file(const char *fn, size_t *sze, char **ret) {
    FILE *f     = fopen(fn, "rb");

    fseek(f, 0, SEEK_END);
    *sze = ftell(f);
    fseek(f, 0, SEEK_SET);

    *ret = (char*) malloc( (*sze)* sizeof(*ret) );

    printf("reading %lu byte from %s\n", *sze, fn);
    if ( *sze != fread(*ret, sizeof(char), *sze, f) ) {
        printf("error reading file %s\n", fn);
        exit(0);
    }
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: read_header <filename>\n");
        return 0;
    }

    char *dat;
    size_t sze;
    read_file(argv[1], &sze, &dat);

    unsigned char *datep = (unsigned char*) dat + 176;

    printf("%x %x %x %x\n", datep[0], datep[1], datep[2], datep[3]);

    unsigned int *epoch = (unsigned int*) datep;
    printf("%u", *epoch);

    return 0;
}
