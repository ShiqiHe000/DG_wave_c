#ifndef DG_INDEX_LOCAL_GLOBAL_H
#define DG_INDEX_LOCAL_GLOBAL_H

int Index_local_to_global(int rank, int l_index);

int Index_global_to_local(int rank, int g_index);

#endif
