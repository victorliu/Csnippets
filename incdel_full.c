#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Incremental delaunay triangulation using halfedge datastructure.
 *
 * A halfedge is an oriented edge. Each interior edge in the triangulation
 * is associated with two halfedges pointing in opposite directions.
 * Each triangle is represented as a ring of three halfedges pointing
 * to each other in CCW order.
 */

/* Vertex data structure */
typedef struct{
	double r[2]; /* x and y coordinates */
	int h; /* index of an outgoing halfedge */
	/* If a point lies on the boundary of the convex hull, then the halfedge
	 * is guaranteed to be "rewound" in the sense that it is on the boundary.
	 */
} Vert;

/* Halfedge data structure */
typedef struct{
	int next; /* index of the next halfedge in the ring */
	int flip; /* index of the opposite halfedge across an edge (-1 if none) */
	int from; /* index of origin vertex */
	int face; /* index of the parent face */
	int edge; /* index of the parent edge */
} Half;

typedef struct{
	int h;
} Edge;

typedef struct{
	int h;
} Face;

typedef struct{
	/* vector of vertices */
	int nv, nv_alloc;
	Vert *v;
	
	/* vector of halfedges */
	int nh, nh_alloc;
	Half *h;
	
	/* vector of edges */
	int ne, ne_alloc;
	Edge *e;
	
	/* vector of faces */
	int nf, nf_alloc;
	Face *f;
	/* The halfedge of a face is the one (of three) with the smallest index */
	
	int boundary; /* pointer to a boundary halfedge */
} DT;

#define NEXT(H) (T->h[H].next)
#define FLIP(H) (T->h[H].flip)
#define FROM(H) (T->h[H].from)
#define FACE(H) (T->h[H].face)
#define EDGE(H) (T->h[H].edge)
#define SETHALF(H,NXT,OPP,ORG,TRI,EDG) do{ \
	T->h[H].next = NXT; \
	T->h[H].flip = OPP; \
	T->h[H].from = ORG; \
	T->h[H].face = TRI; \
	T->h[H].edge = EDG; \
}while(0)

void DT_dump(const DT *T){
	int i;
	printf("{\nnxt : opp : org : tri\n");
	for(i = 0; i < T->nh; ++i){
		printf("%d\t%3d : %3d : %3d : %3d\n",
			i, NEXT(i), FLIP(i), FROM(i), FACE(i)
		);
	}
	printf("}\n");
}

/* Create a new initial DT with the following layout:
 *                c
 *                +
 *              .'|
 *            .'  |
 *       e2 .'    |
 *        .'h2  h1|e1
 *      .'    f0  |
 *    .'    h0    |
 *   +------------+
 *  a      e0      b
 */
DT* DT_new(const Vert *a, const Vert *b, const Vert *c){
	DT *T = (DT*)malloc(sizeof(DT));
	/* Initialize vertex list */
	T->nv = 3;
	T->nv_alloc = 4;
	T->v = (Vert*)malloc(sizeof(Vert) * T->nv_alloc);
	/* Initialize halfedge list */
	T->nh = 3;
	T->nh_alloc = 4;
	T->h = (Half*)malloc(sizeof(Half) * T->nh_alloc);
	/* Initialize edge list */
	T->ne = 3;
	T->ne_alloc = 4;
	T->e = (Edge*)malloc(sizeof(Edge) * T->ne_alloc);
	/* Initialize face list */
	T->nf = 1;
	T->nf_alloc = 4;
	T->f = (Face*)malloc(sizeof(Face) * T->nf_alloc);
	/* Add the vertices */
	memcpy(&T->v[0], a, sizeof(Vert));
	memcpy(&T->v[1], b, sizeof(Vert));
	memcpy(&T->v[2], c, sizeof(Vert));
	T->v[0].h = 0;
	T->v[1].h = 1;
	T->v[2].h = 2;
	/* Create the new halfedges */
	SETHALF(0, 1, -1, 0, 0, 0);
	SETHALF(1, 2, -1, 1, 0, 1);
	SETHALF(2, 0, -1, 2, 0, 2);
	T->f[0].h = 0;
	T->boundary = 0;
	return T;
}
void DT_destroy(DT *T){
	if(NULL == T){ return; }
	free(T->h);
	free(T->f);
	free(T->e);
	free(T->v);
	free(T);
}

inline double orient2d(
	const double a[2], const double b[2], const double c[2]
){
	return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
}
double incircle(
	const double a[2], const double b[2], const double c[2], const double d[2]
){ /* Returns positive if the point d is inside the circle abc */
	return
		(a[0]*a[0] + a[1]*a[1]) * orient2d(b, c, d) -
		(b[0]*b[0] + b[1]*b[1]) * orient2d(a, c, d) +
		(c[0]*c[0] + c[1]*c[1]) * orient2d(a, b, d) -
		(d[0]*d[0] + d[1]*d[1]) * orient2d(a, b, c);
}

/* Returns a halfedge of a face in ih (smallest index ih)
 * Returns 0 if inside a triangle, 1 if on exterior
 */
int DT_locate(const DT *T, const double r[2], int *ih){
	while(1){
		const int h0 = *ih;
		const int h1 = NEXT(h0);
		const int h2 = NEXT(h1);
		int nih;
		if(orient2d(T->v[FROM(h0)].r, T->v[FROM(h1)].r, r) < 0){
			nih = FLIP(*ih);
			if(-1 == nih){ return 1; }else{ *ih = nih; }
		}else if(orient2d(T->v[FROM(h1)].r, T->v[FROM(h2)].r, r) < 0){
			*ih = NEXT(*ih);
			nih = FLIP(*ih);
			if(-1 == nih){ return 1; }else{ *ih = nih; }
		}else if(orient2d(T->v[FROM(h2)].r, T->v[FROM(h0)].r, r) < 0){
			*ih = NEXT(*ih);
			*ih = NEXT(*ih);
			nih = FLIP(*ih);
			if(-1 == nih){ return 1; }else{ *ih = nih; }
		}else{
			if(h1 < *ih){ *ih = h1; }
			if(h2 < *ih){ *ih = h2; }
			return 0;
		}
	}
}

void DT_flip(DT *T, int h){
	/* assert(0 <= h && h < T->nh && -1 != FLIP(h)); */
	const int h0 = h;
	const int h1 = FLIP(h0), h2 = NEXT(h0), h3 = NEXT(h2);
	const int h4 = NEXT(h1), h5 = NEXT(h4);
	const int ia = FROM(h0), ib = FROM(h2);
	const int ic = FROM(h3), id = FROM(h5);
	const int f0 = FACE(h0), f1 = FACE(h1);
	/* Flips an edge in its parent quadrilateral
	 *
	 *     c                  b          c                  b
	 *      +----------------+            +----------------+
	 *      |       h2     .'|            |`.     h2       |
	 *      |   f0       .'  |            |  `.      f1    |
	 *      |          .'    |            |    `.          |
	 *      |h3    h0.'      |            |h3    `.h1    h5|
	 *      |      .' h1   h5|    ===>    |      h0`.      |
	 *      |    .'          |            |          `.    |
	 *      |  .'      f1    |            |   f0       `.  |
	 *      |.'     h4       |            |       h4     `.|
	 *      +----------------+            +----------------+
	 *     a                  d          a                  d
	 */
	NEXT(h0) = h3; FROM(h0) = id;
	NEXT(h1) = h5; FROM(h1) = ic;
	NEXT(h2) = h1; FACE(h2) = f1;
	NEXT(h3) = h4;
	NEXT(h4) = h0; FACE(h4) = f0;
	NEXT(h5) = h2;
	T->v[ia].h = h4; T->v[ib].h = h2;
	T->v[ic].h = h3; T->v[id].h = h5;
	
	/* Re-set the face pointers */
	T->f[f0].h = h0;
	if(h3 < h0){ T->f[f0].h = h3; }
	if(h4 < T->f[f0].h){ T->f[f0].h = h4; }
	
	T->f[f1].h = h1;
	if(h2 < h1){ T->f[f1].h = h2; }
	if(h5 < T->f[f1].h){ T->f[f1].h = h5; }
}

/* h should be a halfedge emanating from the newly added vertex.
 * We assume that h is "rewound" so that if the new vertex is on
 * the boundary, h is on the boundary.
 */
void DT_fixup(DT *T, int h){
	const int h0 = h;
	const int ip = FROM(h);
	h = NEXT(h);
	do{
		const int ht = FLIP(h);
		if(-1 != ht){
			const int ia = FROM(ht);
			const int ib = FROM(NEXT(ht));
			const int ic = FROM(NEXT(NEXT(ht)));
			if(/*orient2d(T->v[ia].r, T->v[ic].r, T->v[ip].r) < 0 &&*/
				incircle(T->v[ia].r, T->v[ib].r, T->v[ic].r, T->v[ip].r) > 0
			){
				DT_flip(T, h);
				h = NEXT(NEXT(h));
			}else{ /* advance */
				h = FLIP(NEXT(h));
				if(-1 == h || h == h0){ break; }
				h = NEXT(h);
			}
		}else{
			h = FLIP(NEXT(h));
			if(-1 == h || h == h0){ break; }
			h = NEXT(h);
		}
	}while(1);
}

/* Allocate more halfedges */
static void DT_addhalfs(DT *T, int n){ /* Only performs allocation */
	if(T->nh+n >= T->nh_alloc){
		while(T->nh+n >= T->nh_alloc){ T->nh_alloc *= 2; }
		T->h = (Half*)realloc(T->h, sizeof(Half) * T->nh_alloc);
	}
}
static void DT_addedges(DT *T, int n){ /* Only performs allocation */
	if(T->ne+n >= T->ne_alloc){
		while(T->ne+n >= T->ne_alloc){ T->ne_alloc *= 2; }
		T->e = (Edge*)realloc(T->e, sizeof(Edge) * T->ne_alloc);
	}
}
static void DT_addface(DT *T, int h){
	if(T->nf >= T->nf_alloc){
		T->nf_alloc *= 2;
		T->f = (Face*)realloc(T->f, sizeof(Face) * T->nf_alloc);
	}
	T->f[T->nf].h = h;
	T->nf++;
}

/* Add a vertex to the triangulation.
 *
 * The incremental algorithm works as follows:
 *   Locate the triangle that contains the new vertex p.
 *   If p is not in an triangle (outside convex hull) then
 *     Join p to all boundary vertices visible from p.
 *   Else
 *     Connect p to the 3 vertices of the triangle containing it.
 *   While an edge in the boundary of the 1-ring of p is not Delaunay,
 *     Flip it.
 *
 * The point location is performed by walking.
 * The visibility testing is performed by walking along the boundary.
 * The edge flipping is determined by checking the triangles adjacent
 * to the 1-ring of p, to see if p is contained in the circumcircle.
 */
void DT_add(DT *T, const Vert *p){
	int h0 = 0;
	int outside;
	const int ip = T->nv;
	
	/* Add the vertex*/
	if(T->nv >= T->nv_alloc){
		T->nv_alloc *= 2;
		T->v = (Vert*)realloc(T->v, sizeof(Vert) * T->nv_alloc);
	}
	memcpy(&T->v[T->nv], p, sizeof(Vert));
	T->nv++;
	
	outside = DT_locate(T, p->r, &h0);
	if(outside){
		int hp = T->nh+1; /* Pointer to current prev boundary edge wrt p */
		int hn = T->nh+2; /* Pointer to current next boundary edge wrt p */
		int ia, ib;
		/* Find all visible edges, crawl out from h0 */
		int ht; /* temporary halfedge pointer that sits on the boundary */
		
		/* Add first triangle */
		ia = FROM(h0);
		ib = FROM(NEXT(h0));
		FLIP(h0) = T->nh+0;
		DT_addhalfs(T, 3);
		DT_addedges(T, 2);
		SETHALF(T->nh+0, hp   , h0, ib, T->nf, EDGE(h0));
		SETHALF(T->nh+1, hn   , -1, ia, T->nf, T->ne+0 );
		SETHALF(T->nh+2, T->nh, -1, ip, T->nf, T->ne+1 );
		T->e[T->ne+0].h = T->nh+1;
		T->e[T->ne+1].h = T->nh+2;
		DT_addface(T, T->nh);
		T->v[T->nv-1].h = hn;
		T->nh += 3; T->ne += 2;
		
		ht = h0;
		do{ /* loop until not visible */
			/* advance ht forwards */
			ht = NEXT(ht);
			while(-1 != FLIP(ht)){
				ht = NEXT(FLIP(ht));
			}
			ia = FROM(ht);
			ib = FROM(NEXT(ht));
			if(orient2d(T->v[ia].r, T->v[ib].r, p->r) < 0){ /* if visible */
				FLIP(ht) = T->nh+0;
				FLIP(hn) = T->nh+1;
				DT_addhalfs(T, 3);
				DT_addedges(T, 1);
				SETHALF(T->nh+0, T->nh+1, ht, ib, T->nf, EDGE(ht));
				SETHALF(T->nh+1, T->nh+2, hn, ia, T->nf, EDGE(hn));
				SETHALF(T->nh+2, T->nh+0, -1, ip, T->nf, T->ne   );
				T->e[T->ne].h = T->nh+2;
				DT_addface(T, T->nh);
				hn = T->nh+2;
				
				T->v[T->nv-1].h = hn;
				T->boundary = hn;
				T->nh += 3; T->ne++;
			}else{ break; }
		}while(1);
		
		ht = h0;
		do{ /* loop until not visible */
			/* advance ht backwards */
			ht = NEXT(NEXT(ht));
			while(-1 != FLIP(ht)){
				ht = NEXT(NEXT(FLIP(ht)));
			}
			ia = FROM(ht);
			ib = FROM(NEXT(ht));
			if(orient2d(T->v[ia].r, T->v[ib].r, p->r) < 0){ /* if visible */
				FLIP(ht) = T->nh+0;
				FLIP(hp) = T->nh+2;
				DT_addhalfs(T, 3);
				DT_addedges(T, 1);
				SETHALF(T->nh+0, T->nh+1, ht, ib, T->nf, EDGE(ht));
				SETHALF(T->nh+1, T->nh+2, -1, ia, T->nf, T->ne   );
				SETHALF(T->nh+2, T->nh+0, hp, ip, T->nf, EDGE(hp));
				T->e[T->ne].h = T->nh+1;
				DT_addface(T, T->nh);
				
				hp = T->nh+1;
				T->nh += 3; T->ne++;
			}else{ break; }
		}while(1);
	}else{ /* inside convex hull */
		const int h1 = NEXT(h0), h2 = NEXT(h1);
		const int ia = FROM(h0), ib = FROM(h1), ic = FROM(h2);
		const int f0 = FACE(h0);
		/*                    c
		 *                    +
		 *                   /|\
		 *                  / | \
		 *                 /  |  \
		 *                /  4|5  \
		 *               / fl1|fl0 \
		 *              /h2   |   h1\
		 *             /      +p     \
		 *            /  1 .-' `-. 2  \
		 *           /  .-'0 f0  3`-.  \
		 *          /.-'     h0      `-.\
		 *        a+---------------------+b
		 */
		/* connect p to vertices */
		DT_addhalfs(T, 6);
		SETHALF(T->nh+0, h0     , T->nh+1, ip, f0     , T->ne+0);
		SETHALF(T->nh+1, T->nh+4, T->nh+0, ia, T->nf+1, T->ne+0);
		SETHALF(T->nh+2, h1     , T->nh+3, ip, T->nf+0, T->ne+1);
		SETHALF(T->nh+3, T->nh+0, T->nh+2, ib, f0     , T->ne+1);
		SETHALF(T->nh+4, h2     , T->nh+5, ip, T->nf+1, T->ne+2);
		SETHALF(T->nh+5, T->nh+2, T->nh+4, ic, T->nf+0, T->ne+2);
		DT_addedges(T, 3);
		T->e[T->ne+0].h = T->nh+0;
		T->e[T->ne+1].h = T->nh+2;
		T->e[T->ne+2].h = T->nh+4;
		NEXT(h0) = T->nh+3;
		NEXT(h1) = T->nh+5; FACE(h1) = T->nf+0;
		NEXT(h2) = T->nh+1; FACE(h2) = T->nf+1;
		T->v[T->nv-1].h = T->nh;
		T->f[f0].h = h0;
		DT_addface(T, h1);
		DT_addface(T, h2);
		T->nh += 6; T->ne += 3;
	}
	DT_fixup(T, T->v[T->nv-1].h);
}

double frand(){ return (double)rand()/(double)RAND_MAX; }

int main(int argc, char *argv[]){
	int i;
	int n = atoi(argv[1]);
	Vert a, b, c;
	DT *T;
	a.r[0] = 0;
	a.r[1] = 0;
	b.r[0] = 1;
	b.r[1] = 0;
	c.r[0] = 1;
	c.r[1] = 1;
	T = DT_new(&a, &b, &c);
	srand(0);
	for(i = 0; i < n; ++i){
		a.r[1] = frand();
		a.r[0] = frand();
		/* printf("%% adding %d: %g, %g\n", (int)(i+3), a.r[0], a.r[1]); */
		DT_add(T, &a);
	}
	/*DT_dump(T);*/
	
	printf(
		"%%!\n"
		"72 72 scale\n"
		"1 4 translate\n"
		"0.001 setlinewidth\n"
		"/Helvetica findfont 0.03 scalefont setfont\n"
	);
	for(i = 0; i < T->nh; ++i){
		int ih = i;
		int ia = T->h[ih].from;
		int ib = T->h[T->h[ih].next].from;
		int ic = T->h[T->h[T->h[ih].next].next].from;
		const double a[2] = { T->v[ia].r[0], T->v[ia].r[1] };
		const double b[2] = { T->v[ib].r[0], T->v[ib].r[1] };
		const double c[2] = { T->v[ic].r[0], T->v[ic].r[1] };
		const double cc[2] = { (a[0]+b[0]+c[0])/3., (a[1]+b[1]+c[1])/3. };
		printf("newpath %g %g moveto %g %g lineto %g %g lineto stroke\n",
			0.94*a[0]+0.02*b[0]+0.04*cc[0],
			0.94*a[1]+0.02*b[1]+0.04*cc[1],
			0.94*b[0]+0.02*a[0]+0.04*cc[0],
			0.94*b[1]+0.02*a[1]+0.04*cc[1],
			0.84*b[0]+0.08*a[0]+0.08*cc[0],
			0.84*b[1]+0.08*a[1]+0.08*cc[1]
		);
	}
	
	for(i = 0; i < T->nf; ++i){
		int ih = T->f[i].h;
		int ia = T->h[ih].from;
		int ib = T->h[T->h[ih].next].from;
		int ic = T->h[T->h[T->h[ih].next].next].from;
		const double a[2] = { T->v[ia].r[0], T->v[ia].r[1] };
		const double b[2] = { T->v[ib].r[0], T->v[ib].r[1] };
		const double c[2] = { T->v[ic].r[0], T->v[ic].r[1] };
		printf("%g %g moveto %g %g lineto %g %g lineto closepath stroke\n",
			a[0], a[1], b[0], b[1], c[0], c[1]
		);
		/*printf("newpath %g %g moveto (%d) show\n", cc[0], cc[1], (int)i);*/
	}
	for(i = 0; i < T->nv; ++i){
		printf("newpath %g %g moveto (%d) show\n",
			T->v[i].r[0], T->v[i].r[1], (int)i
		);
	}
	printf("showpage\n");
	DT_destroy(T);
	return 0;
}

