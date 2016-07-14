/*
* This file is part of the Anthony Lomax C Library.
*
* Copyright (C) 2006-2010 Anthony Lomax <anthony@alomax.net www.alomax.net>
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 0

#define MAX_NZEROS 8
#define MAX_NPOLES 8
double xv[MAX_NZEROS+1], yv[MAX_NPOLES+1];

#include "timedomain_memory.h"
#include "timedomain_filter.h"

#include "timedomain_filter_1.c"
#include "timedomain_filter_2.c"
#include "timedomain_filter_3.c"
#include "timedomain_filter_4.c"
#include "timedomain_filter_6_VEL.c"
#include "timedomain_filter_6_DISP.c"

