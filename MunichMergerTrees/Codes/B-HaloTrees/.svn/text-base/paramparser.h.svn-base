#ifndef SIMPOSTPROC_PARSER_H
#define SIMPOSTPROC_PARSER_H

#ifdef __cplusplus
extern "C" {
#endif

void *openParams (const char *file);
void closeParams (void *vpar);
int isVar(void *vpar, const char *key);
int getVarInt (void *vpar, const char *key, int *result);
int getVarLongLong (void *vpar, const char *key, long long *result);
int getVarDouble (void *vpar, const char *key, double *result);
int getVarString (void *vpar, const char *key, char **result);
int getVarStringMem (void *vpar, const char *key, char **result);
int getVarStringProt (void *vpar, const char *key, char *result, int capacity);
int getVarStringUnProt (void *vpar, const char *key, char *result);

#ifdef __cplusplus
}
#endif

#endif
