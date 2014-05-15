#include <stdio.h>
#include <math.h>

/*     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have */
/*       the same meanings as the corresponding arguments of BOBYQB. */
/*     KOPT is the index of the optimal interpolation point. */
/*     KNEW is the index of the interpolation point that is going to be moved. */
/*     ADELT is the current trust region bound. */
/*     XNEW will be set to a suitable new position for the interpolation point */
/*       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region */
/*       bounds and it should provide a large denominator in the next call of */
/*       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the */
/*       straight lines through XOPT and another interpolation point. */
/*     XALT also provides a large value of the modulus of the KNEW-th Lagrange */
/*       function subject to the constraints that have been mentioned, its main */
/*       difference from XNEW being that XALT-XOPT is a constrained version of */
/*       the Cauchy step within the trust region. An exception is that XALT is */
/*       not calculated if all components of GLAG (see below) are zero. */
/*     ALPHA will be set to the KNEW-th diagonal element of the H matrix. */
/*     CAUCHY will be set to the square of the KNEW-th Lagrange function at */
/*       the step XALT-XOPT from XOPT for the vector XALT that is returned, */
/*       except that CAUCHY is set to zero if XALT is not calculated. */
/*     GLAG is a working space vector of length N for the gradient of the */
/*       KNEW-th Lagrange function at XOPT. */
/*     HCOL is a working space vector of length NPT for the second derivative */
/*       coefficients of the KNEW-th Lagrange function. */
/*     W is a working space vector of length 2N that is going to hold the */
/*       constrained Cauchy step from XOPT of the Lagrange function, followed */
/*       by the downhill version of XALT when the uphill step is calculated. */

void altmov_(
	int n, int npt, double *xpt, 
	double *xopt, double *bmat, double *zmat, int ndim,
	double *sl, double *su,
	int kopt, int knew, double adelt,
	double *xnew, double *xalt,
	double *alpha,
	double *cauchy,
	/* Workspaces */
	double *glag, /* length n */
	double *hcol, /* length npt */
	double *w     /* length 2*n */
){
	/* System generated locals */
	int xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1, 
		zmat_offset;
	double d__1, d__2;

	/* Local variables */
	int i, j, k;
	double ha, gw, one, diff, half;
	int ilbd, isbd;
	double slbd;
	int iubd;
	double vlag, subd, temp;
	int ksav;
	double step, zero, curv;
	int iflag;
	double scale, csave, tempa, tempb, tempd, const__, sumin, 
		ggfree;
	int ibdsav;
	double dderiv, bigstp, predsq, presav, distsq, stpsav, wfixsq, 
		wsqsav;




	/* Parameter adjustments */
	zmat_dim1 = npt;
	zmat_offset = 1 + zmat_dim1;
	zmat -= zmat_offset;
	xpt_dim1 = npt;
	xpt_offset = 1 + xpt_dim1;
	xpt -= xpt_offset;
	--xopt;
	bmat_dim1 = ndim;
	bmat_offset = 1 + bmat_dim1;
	bmat -= bmat_offset;
	--sl;
	--su;
	--xnew;
	--xalt;
	--glag;
	--hcol;
	--w;





	/* Set the first NPT components of W to the leading elements of the */
	/* KNEW-th column of the H matrix. */
	half = .5;
	one = 1.;
	zero = 0.;
	const__ = one + sqrt(2.);
	for(k = 1; k <= npt; ++k){
		hcol[k] = 0;
	}
	for(j = 1; j <= npt - n - 1; ++j){
		temp = zmat[knew + j * zmat_dim1];
		for(k = 1; k <= npt; ++k){
			hcol[k] += temp * zmat[k + j * zmat_dim1];
		}
	}
	*alpha = hcol[knew];
	ha = half * *alpha;

	/* Calculate the gradient of the KNEW-th Lagrange function at XOPT. */

	for(i = 1; i <= n; ++i){
		glag[i] = bmat[knew + i * bmat_dim1];
	}
	for(k = 1; k <= npt; ++k){
		temp = zero;
		for(j = 1; j <= n; ++j){
			temp += xpt[k + j * xpt_dim1] * xopt[j];
		}
		temp = hcol[k] * temp;
		for(i = 1; i <= n; ++i){
			glag[i] += temp * xpt[k + i * xpt_dim1];
		}
	}

	/* Search for a large denominator along the straight lines through XOPT */
	/* and another interpolation point. SLBD and SUBD will be lower and upper */
	/* bounds on the step along each of these lines in turn. PREDSQ will be */
	/* set to the square of the predicted denominator for each line. PRESAV */
	/* will be set to the largest admissible value of PREDSQ that occurs. */

	presav = zero;
	for(k = 1; k <= npt; ++k){
		if(k == kopt){
			continue;
		}
		dderiv = zero;
		distsq = zero;
		for(i = 1; i <= n; ++i){
			temp = xpt[k + i * xpt_dim1] - xopt[i];
			dderiv += glag[i] * temp;
			distsq += temp * temp;
		}
		subd = adelt / sqrt(distsq);
		slbd = -subd;
		ilbd = 0;
		iubd = 0;
		sumin = subd;
		if(1 < sumin){ sumin = 1; }

		/* Revise SLBD and SUBD if necessary because of the bounds in SL and SU. */

		for(i = 1; i <= n; ++i){
			temp = xpt[k + i * xpt_dim1] - xopt[i];
			if(temp > zero){
			if(slbd * temp < sl[i] - xopt[i]){
				slbd = (sl[i] - xopt[i]) / temp;
				ilbd = -i;
			}
			if(subd * temp > su[i] - xopt[i]){
				subd = (su[i] - xopt[i]) / temp;
				if(sumin > subd){ subd = sumin; }
				iubd = i;
			}
			}else if(temp < zero){
			if(slbd * temp > su[i] - xopt[i]){
				slbd = (su[i] - xopt[i]) / temp;
				ilbd = i;
			}
			if(subd * temp < sl[i] - xopt[i]){
				subd = (sl[i] - xopt[i]) / temp;
				if(sumin > subd){ subd = sumin; }
				iubd = -i;
			}
			}
		}

		/* Seek a large modulus of the KNEW-th Lagrange function when the index */
		/* of the other interpolation point on the line through XOPT is KNEW. */

		if(k == knew){
			diff = dderiv - one;
			step = slbd;
			vlag = slbd * (dderiv - slbd * diff);
			isbd = ilbd;
			temp = subd * (dderiv - subd * diff);
			if(fabs(temp) > fabs(vlag)){
			step = subd;
			vlag = temp;
			isbd = iubd;
			}
			tempd = half * dderiv;
			tempa = tempd - diff * slbd;
			tempb = tempd - diff * subd;
			if(tempa * tempb < zero){
			temp = tempd * tempd / diff;
			if(fabs(temp) > fabs(vlag)){
				step = tempd / diff;
				vlag = temp;
				isbd = 0;
			}
			}

		/* Search along each of the other lines through XOPT and another point. */

		}else{
			step = slbd;
			vlag = slbd * (one - slbd);
			isbd = ilbd;
			temp = subd * (one - subd);
			if(fabs(temp) > fabs(vlag)){
			step = subd;
			vlag = temp;
			isbd = iubd;
			}
			if(subd > half){
				if(fabs(vlag) < .25){
					step = half;
					vlag = .25;
					isbd = 0;
				}
			}
			vlag *= dderiv;
		}

		/*     Calculate PREDSQ for the current line search and maintain PRESAV. */

		temp = step * (one - step) * distsq;
		predsq = vlag * vlag * (vlag * vlag + ha * temp * temp);
		if(predsq > presav){
			presav = predsq;
			ksav = k;
			stpsav = step;
			ibdsav = isbd;
		}
	}

	/*     Construct XNEW in a way that satisfies the bound constraints exactly. */

	for(i = 1; i <= n; ++i){
		xnew[i] = xopt[i] + stpsav * (xpt[ksav + i * xpt_dim1] - xopt[i]);
		if(su[i] < xnew[i]){ xnew[i] = su[i]; }
		if(sl[i] > xnew[i]){ xnew[i] = sl[i]; }
	}
	if(ibdsav < 0){
		xnew[-ibdsav] = sl[-ibdsav];
	}
	if(ibdsav > 0){
		xnew[ibdsav] = su[ibdsav];
	}

	/*     Prepare for the iterative method that assembles the constrained Cauchy */
	/*     step in W. The sum of squares of the fixed components of W is formed in */
	/*     WFIXSQ, and the free components of W are set to BIGSTP. */

	bigstp = adelt + adelt;
	iflag = 0;
L100:
	wfixsq = zero;
	ggfree = zero;
	for(i = 1; i <= n; ++i){
		w[i] = zero;
		tempa = xopt[i] - sl[i];
		if(glag[i] < tempa){ tempa = glag[i]; }
		tempb = xopt[i] - su[i];
		if(glag[i] > tempb){ tempb = glag[i]; }
		if(tempa > zero || tempb < zero){
			w[i] = bigstp;
			ggfree += glag[i] * glag[i];
		}
	}
	if(ggfree == zero){
		*cauchy = zero;
		return;
	}

	/*     Investigate whether more components of W can be fixed. */

L120:
	do{
		temp = adelt * adelt - wfixsq;
		if(temp > zero){
			wsqsav = wfixsq;
			step = sqrt(temp / ggfree);
			ggfree = zero;
			for(i = 1; i <= n; ++i){
				if(w[i] == bigstp){
					temp = xopt[i] - step * glag[i];
					if(temp <= sl[i]){
						w[i] = sl[i] - xopt[i];
						wfixsq += w[i] * w[i];
					}else if(temp >= su[i]){
						w[i] = su[i] - xopt[i];
						wfixsq += w[i] * w[i];
					}else{
						ggfree += glag[i] * glag[i];
					}
				}
			}
			if(wfixsq > wsqsav && ggfree > zero){
				goto L120;
			}
		}
	}while(0);

	/*     Set the remaining free components of W and all components of XALT, */
	/*     except that W may be scaled later. */

	gw = zero;
	for(i = 1; i <= n; ++i){
		if(w[i] == bigstp){
			w[i] = -step * glag[i];
			xalt[i] = xopt[i] + w[i];
			if(su[i] < xalt[i]){ xalt[i] = su[i]; }
			if(sl[i] > xalt[i]){ xalt[i] = sl[i]; }
		}else if(w[i] == zero){
			xalt[i] = xopt[i];
		}else if(glag[i] > zero){
			xalt[i] = sl[i];
		}else{
			xalt[i] = su[i];
		}
		gw += glag[i] * w[i];
	}

	/*     Set CURV to the curvature of the KNEW-th Lagrange function along W. */
	/*     Scale W by a factor less than one if that can reduce the modulus of */
	/*     the Lagrange function at XOPT+W. Set CAUCHY to the final value of */
	/*     the square of this function. */

	curv = zero;
	for(k = 1; k <= npt; ++k){
		temp = zero;
		for(j = 1; j <= n; ++j){
			temp += xpt[k + j * xpt_dim1] * w[j];
		}
		curv += hcol[k] * temp * temp;
	}
	if(iflag == 1){
		curv = -curv;
	}
	if(curv > -gw && curv < -const__ * gw){
		scale = -gw / curv;
		for(i = 1; i <= n; ++i){
			xalt[i] = xopt[i] + scale * w[i];
			if(su[i] < xalt[i]){ xalt[i] = su[i]; }
			if(sl[i] > xalt[i]){ xalt[i] = sl[i]; }
		}
		/* Computing 2nd power */
		d__1 = half * gw * scale;
		*cauchy = d__1 * d__1;
	}else{
	/* Computing 2nd power */
		d__1 = gw + half * curv;
		*cauchy = d__1 * d__1;
	}

	/* If IFLAG is zero, then XALT is calculated as before after reversing */
	/* the sign of GLAG. Thus two XALT vectors become available. The one that */
	/* is chosen is the one that gives the larger value of CAUCHY. */

	if(iflag == 0){
		for(i = 1; i <= n; ++i){
			glag[i] = -glag[i];
			w[n + i] = xalt[i];
		}
		csave = *cauchy;
		iflag = 1;
		goto L100;
	}
	if(csave > *cauchy){
		for(i = 1; i <= n; ++i){
			xalt[i] = w[n + i];
		}
		*cauchy = csave;
	}
}

/*     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same */
/*       meanings as the corresponding arguments of BOBYQB. */
/*     DELTA is the trust region radius for the present calculation, which */
/*       seeks a small value of the quadratic model within distance DELTA of */
/*       XOPT subject to the bounds on the variables. */
/*     XNEW will be set to a new vector of variables that is approximately */
/*       the one that minimizes the quadratic model within the trust region */
/*       subject to the SL and SU constraints on the variables. It satisfies */
/*       as equations the bounds that become active during the calculation. */
/*     D is the calculated trial step from XOPT, generated iteratively from an */
/*       initial value of zero. Thus XNEW is XOPT+D after the final iteration. */
/*     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated */
/*       when D is updated. */
/*     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is */
/*       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the */
/*       I-th variable has become fixed at a bound, the bound being SL(I) or */
/*       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This */
/*       information is accumulated during the construction of XNEW. */
/*     The arrays S, HS and HRED are also used for working space. They hold the */
/*       current search direction, and the changes in the gradient of Q along S */
/*       and the reduced D, respectively, where the reduced D is the same as D, */
/*       except that the components of the fixed variables are zero. */
/*     DSQ will be set to the square of the length of XNEW-XOPT. */
/*     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise */
/*       it is set to the least curvature of H that occurs in the conjugate */
/*       gradient searches that are not restricted by any constraints. The */
/*       value CRVMIN=-1.0D0 is set, however, if all of these searches are */
/*       constrained. */

/*     A version of the truncated conjugate gradient is applied. If a line */
/*     search is restricted by a constraint, then the procedure is restarted, */
/*     the values of the variables that are at their bounds being fixed. If */
/*     the trust region boundary is reached, then further changes may be made */
/*     to D, each one being in the two dimensional space that is spanned */
/*     by the current D and the gradient of Q at XOPT+D, staying on the trust */
/*     region boundary. Termination occurs when the reduction in Q seems to */
/*     be close to the greatest reduction that can be achieved. */

void trsbox_(int n, int npt, double *xpt, 
	double *xopt, double *gopt, double *hq, double *pq, 
	double *sl, double *su, double delta, double *xnew, 
	double *d, double *gnew, double *xbdi, double *s, 
	double *hs, double *hred, double *dsq, double *crvmin)
{
    /* System generated locals */
    int xpt_offset;

    /* Local variables */
    int i, j, k, ih;
    double ds;
    int iu;
    double dhd, dhs, cth, shs, sth, ssq, beta, sdec, blen;
    int iact, nact;
    double angt, qred;
    int isav;
    double temp, xsav, xsum, angbd, dredg, sredg;
    int iterc;
    double resid, delsq, ggsav, tempa, tempb, ratio, sqstp, redmax, dredsq, redsav, gredsq, rednew;
    int itcsav;
    double rdprev, rdnext, stplen, stepsq;
    int itermax;

	const double half = .5;
	const double one = 1.;
	const double onemin = -1.;
	const double zero = 0.;

	/* The sign of GOPT(I) gives the sign of the change to the I-th variable */
	/* that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether */
	/* or not to fix the I-th variable at one of its bounds initially, with */
	/* NACT being set to the number of fixed variables. D and GNEW are also */
	/* set for the first iteration. DELSQ is the upper bound on the sum of */
	/* squares of the free variables. QRED is the reduction in Q so far. */

    iterc = 0;
    nact = 0;
    sqstp = 0;
    for(i = 0; i < n; ++i){
		xbdi[i] = 0;
		if(xopt[i] <= sl[i]){
			if(gopt[i] >= 0){
				xbdi[i] = onemin;
			}
		}else if(xopt[i] >= su[i]){
			if(gopt[i] <= 0){
				xbdi[i] = one;
			}
		}
		if(xbdi[i] != 0){
			++nact;
		}
		d[i] = 0;
		gnew[i] = gopt[i];
    }
    delsq = delta * delta;
    qred = zero;
    *crvmin = onemin;


    /* Parameter adjustments */
    xpt_offset = 1 + npt;
    xpt -= xpt_offset;
    --xopt;
    --gopt;
    --hq;
    --pq;
    --sl;
    --su;
    --xnew;
    --d;
    --gnew;
    --xbdi;
    --s;
    --hs;
    --hred;
    
	/* Set the next search direction of the conjugate gradient method. It is */
	/* the steepest descent direction initially and when the iterations are */
	/* restarted because a variable has just been fixed by a bound, and of */
	/* course the components of the fixed variables are zero. ITERMAX is an */
	/* upper bound on the indices of the conjugate gradient iterations. */

L20:
    beta = zero;
L30:
    stepsq = zero;
    for(i = 1; i <= n; ++i){
		if(xbdi[i] != zero){
			s[i] = zero;
		}else if(beta == zero){
			s[i] = -gnew[i];
		}else{
			s[i] = beta * s[i] - gnew[i];
		}
		stepsq += s[i] * s[i];
    }
    if(stepsq == zero){
		goto L190;
    }
    if(beta == zero){
		gredsq = stepsq;
		itermax = iterc + n - nact;
    }
    if(gredsq * delsq <= qred * 1e-4 * qred){
		goto L190;
    }

	/* Multiply the search direction by the second derivative matrix of Q and */
	/* calculate some scalars for the choice of steplength. Then set BLEN to */
	/* the length of the the step to the trust region boundary and STPLEN to */
	/* the steplength, ignoring the simple bounds. */

    goto L210;
L50:
    resid = delsq;
    ds = zero;
    shs = zero;
    for(i = 1; i <= n; ++i){
		if(xbdi[i] == zero){
			resid -= d[i] * d[i];
			ds    += s[i] * d[i];
			shs   += s[i] * hs[i];
		}
    }
    if(resid <= zero){
		goto L90;
    }
    temp = sqrt(stepsq * resid + ds * ds);
    if(ds < zero){
		blen = (temp - ds) / stepsq;
    }else{
		blen = resid / (temp + ds);
    }
    stplen = blen;
    if(shs > zero){
		stplen = gredsq / shs;
		if(blen < stplen){ stplen = blen; }
    }

	/* Reduce STPLEN if necessary in order to preserve the simple bounds, */
	/* letting IACT be the index of the new constrained variable. */

    iact = 0;
    for(i = 1; i <= n; ++i){
		if(s[i] != 0){
			xsum = xopt[i] + d[i];
			if(s[i] > zero){
			temp = (su[i] - xsum) / s[i];
			}else{
			temp = (sl[i] - xsum) / s[i];
			}
			if(temp < stplen){
				stplen = temp;
				iact = i;
			}
		}
    }

	/* Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q. */

    sdec = zero;
    if(stplen > zero){
		++iterc;
		temp = shs / stepsq;
		if(iact == 0 && temp > zero){
			if(temp < *crvmin){ *crvmin = temp; }
			if(*crvmin == onemin){
				*crvmin = temp;
			}
		}
		ggsav = gredsq;
		gredsq = zero;
		for(i = 1; i <= n; ++i){
			gnew[i] += stplen * hs[i];
			if(xbdi[i] == zero){
				gredsq += gnew[i] * gnew[i];
			}
			d[i] += stplen * s[i];
		}
		sdec = stplen * (ggsav - half * stplen * shs);
		if(sdec < 0){ sdec = 0; }
		qred += sdec;
    }

	/* Restart the conjugate gradient method if it has hit a new bound. */

    if(iact > 0){
		++nact;
		xbdi[iact] = one;
		if(s[iact] < zero){
			xbdi[iact] = onemin;
		}
		delsq -= d[iact] * d[iact];
		if(delsq <= zero){
			goto L90;
		}
		goto L20;
    }

	/* If STPLEN is less than BLEN, then either apply another conjugate */
	/* gradient iteration or RETURN. */

    if(stplen < blen){
		if(iterc == itermax){
			goto L190;
		}
		if(sdec <= qred * .01){
			goto L190;
		}
		beta = gredsq / ggsav;
		goto L30;
    }
L90:
    *crvmin = zero;

/*     Prepare for the alternative iteration by calculating some scalars and */
/*     by multiplying the reduced D by the second derivative matrix of Q. */

L100:
    if(nact >= n - 1){
		goto L190;
    }
    dredsq = 0;
    dredg = 0;
    gredsq = 0;
    for(i = 1; i <= n; ++i){
		if(xbdi[i] == 0){
			dredsq += d[i] * d[i];
			dredg += d[i] * gnew[i];
			gredsq += gnew[i] * gnew[i];
			s[i] = d[i];
		}else{
			s[i] = zero;
		}
    }
    itcsav = iterc;
    goto L210;

/*     Let the search direction S be a linear combination of the reduced D */
/*     and the reduced G that is orthogonal to the reduced D. */

L120:
    ++iterc;
    temp = gredsq * dredsq - dredg * dredg;
    if(temp <= qred * 1e-4 * qred){
		goto L190;
    }
    temp = sqrt(temp);
    for(i = 1; i <= n; ++i){
		if(xbdi[i] == zero){
			s[i] = (dredg * d[i] - dredsq * gnew[i]) / temp;
		}else{
			s[i] = zero;
		}
    }
    sredg = -temp;

/*     By considering the simple bounds on the variables, calculate an upper */
/*     bound on the tangent of half the angle of the alternative iteration, */
/*     namely ANGBD, except that, if already a free variable has reached a */
/*     bound, there is a branch back to label 100 after fixing that variable. */

    angbd = one;
    iact = 0;
    for(i = 1; i <= n; ++i){
		if(xbdi[i] == zero){
			tempa = xopt[i] + d[i] - sl[i];
			tempb = su[i] - xopt[i] - d[i];
			if(tempa <= zero){
				++nact;
				xbdi[i] = onemin;
				goto L100;
			}else if(tempb <= zero){
				++nact;
				xbdi[i] = one;
				goto L100;
			}
			ratio = one;
			ssq = d[i] * d[i] + s[i] * s[i];
			temp = ssq - (xopt[i] - sl[i]) * (xopt[i] - sl[i]);
			if(temp > zero){
				temp = sqrt(temp) - s[i];
				if(angbd * temp > tempa){
					angbd = tempa / temp;
					iact = i;
					xsav = onemin;
				}
			}
			temp = ssq - (su[i] - xopt[i]) * (su[i] - xopt[i]);
			if(temp > zero){
				temp = sqrt(temp) + s[i];
				if(angbd * temp > tempb){
					angbd = tempb / temp;
					iact = i;
					xsav = one;
				}
			}
		}
    }

	/* Calculate HHD and some curvatures for the alternative iteration. */

    goto L210;
L150:
    shs = zero;
    dhs = zero;
    dhd = zero;
    for(i = 1; i <= n; ++i){
		if(xbdi[i] == zero){
			shs += s[i] * hs[i];
			dhs += d[i] * hs[i];
			dhd += d[i] * hred[i];
		}
    }

	/* Seek the greatest reduction in Q for a range of equally spaced values */
	/* of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of */
	/* the alternative iteration. */

    redmax = zero;
    isav = 0;
    redsav = zero;
    iu = (int) (angbd * 17. + 3.1);
    for(i = 1; i <= iu; ++i){
		angt = angbd * (double) i / (double) iu;
		sth = (angt + angt) / (one + angt * angt);
		temp = shs + angt * (angt * dhd - dhs - dhs);
		rednew = sth * (angt * dredg - sredg - half * sth * temp);
		if(rednew > redmax){
			redmax = rednew;
			isav = i;
			rdprev = redsav;
		}else if(i == isav + 1){
			rdnext = rednew;
		}
		redsav = rednew;
    }

	/* Return if the reduction is zero. Otherwise, set the sine and cosine */
	/* of the angle of the alternative iteration, and calculate SDEC. */

    if(isav == 0){
		goto L190;
    }
    if(isav < iu){
		temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
		angt = angbd * ((double) isav + half * temp) / (double) iu;
    }
    cth = (one - angt * angt) / (one + angt * angt);
    sth = (angt + angt) / (one + angt * angt);
    temp = shs + angt * (angt * dhd - dhs - dhs);
    sdec = sth * (angt * dredg - sredg - half * sth * temp);
    if(sdec <= zero){
		goto L190;
    }

	/* Update GNEW, D and HRED. If the angle of the alternative iteration */
	/* is restricted by a bound on a free variable, that variable is fixed */
	/* at the bound. */

    dredg = zero;
    gredsq = zero;
    for(i = 1; i <= n; ++i){
		gnew[i] = gnew[i] + (cth - one) * hred[i] + sth * hs[i];
		if(xbdi[i] == zero){
			d[i] = cth * d[i] + sth * s[i];
			dredg += d[i] * gnew[i];
			gredsq += gnew[i] * gnew[i];
		}
		hred[i] = cth * hred[i] + sth * hs[i];
    }
    qred += sdec;
    if(iact > 0 && isav == iu){
		++nact;
		xbdi[iact] = xsav;
		goto L100;
    }

	/* If SDEC is sufficiently small, then RETURN after setting XNEW to */
	/* XOPT+D, giving careful attention to the bounds. */

    if(sdec > qred * .01){
		goto L120;
    }
L190:
    *dsq = zero;
    for(i = 1; i <= n; ++i){
		xnew[i] = xopt[i] + d[i];
		if(su[i] < xnew[i]){ xnew[i] = su[i]; }
		if(sl[i] > xnew[i]){ xnew[i] = sl[i]; }
		if(xbdi[i] == onemin){
			xnew[i] = sl[i];
		}
		if(xbdi[i] == one){
			xnew[i] = su[i];
		}
		d[i] = xnew[i] - xopt[i];
		*dsq += d[i] * d[i];
    }
    return;
    
	/* The following instructions multiply the current S-vector by the second */
	/* derivative matrix of the quadratic model, putting the product in HS. */
	/* They are reached from three different parts of the software above and */
	/* they can be regarded as an external subroutine. */

L210:
    ih = 0;
    for(j = 1; j <= n; ++j){
		hs[j] = zero;
		for(i = 1; i <= j; ++i){
			++ih;
			if(i < j){
				hs[j] += hq[ih] * s[i];
			}
			hs[i] += hq[ih] * s[j];
		}
    }
    for(k = 1; k <= npt; ++k){
		if(pq[k] != zero){
			temp = zero;
			for(j = 1; j <= n; ++j){
				temp += xpt[k + j * npt] * s[j];
			}
			temp *= pq[k];
			for(i = 1; i <= n; ++i){
				hs[i] += temp * xpt[k + i * npt];
			}
		}
    }
    if(*crvmin != zero){
		goto L50;
    }
    if(iterc > itcsav){
		goto L150;
    }
    for(i = 1; i <= n; ++i){
		hred[i] = hs[i];
    }
    goto L120;
}

/*
The arrays BMAT and ZMAT are updated, as required by the new position
of the interpolation point that has the index KNEW. The vector VLAG has
N+NPT components, set on entry to the first NPT and last N components
of the product Hw in equation (4.11) of the Powell (2006) paper on
NEWUOA. Further, BETA is set on entry to the value of the parameter
with that name, and DENOM is set to the denominator of the updating
formula. Elements of ZMAT may be treated as zero if their moduli are
at most ZTEST. The first NDIM elements of W are used for working space.
*/

void update_(
	int n, int npt,
	double *bmat, /* size ndim by ? */
	double *zmat, /* size npt by ? */
	int ndim,
	double *vlag, /* length n+npt */
	double beta, double denom,
	int knew, /* 0 based index */
	double *w /* size ndim workspace */
){
    int i, j;
    double tau;
    double alpha, tempa, tempb, ztest;


	const int nptm = npt - n - 1;
	ztest = 0;
	for(j = 0; j < nptm; ++j){
		int k;
		for(k = 0; k < npt; ++k) {
			double tmp = fabs(zmat[k+j*npt]);
			if(tmp > ztest){ ztest = tmp; }
		}
	}
	ztest *= 1e-20;
	
	/* Apply the rotations that put zeros in the KNEW-th row of ZMAT. */
	for(j = 1; j < nptm; ++j){
		if(fabs(zmat[knew+j*npt]) > ztest){
			double temp = hypot(zmat[knew+0*npt], zmat[knew+j*npt]);
			tempa = zmat[knew+0*npt] / temp;
			tempb = zmat[knew+j*npt] / temp;
			for(i = 0; i < npt; ++i){
				temp = tempa * zmat[i+0*npt] + tempb * zmat[i+j*npt];
				zmat[i+j*npt] = tempa * zmat[i+j*npt] - tempb * zmat[i+0*npt];
				zmat[i+0*npt] = temp;
			}
		}
		zmat[knew+j*npt] = 0;
	}

	/* Put the first NPT components of the KNEW-th column of HLAG into W, */
	/* and calculate the parameters of the updating formula. */
	for(i = 0; i < npt; ++i) {
		w[i] = zmat[knew+0*npt] * zmat[i+0*npt];
	}
	alpha = w[knew];
	tau = vlag[knew];
	vlag[knew] -= 1.;

	/* Complete the updating of ZMAT. */
	{
		double temp = sqrt(denom);
		tempb = zmat[knew+0*npt] / temp;
		tempa = tau / temp;
	}
	for(i = 0; i < npt; ++i){
		zmat[i+0*npt] = tempa * zmat[i+0*npt] - tempb * vlag[i];
	}

	/* Finally, update the matrix BMAT. */
	for(j = 0; j < n; ++j){
		const int jp = npt + j;
		w[jp] = bmat[knew+j*ndim];
		tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
		tempb = (-beta * w[jp] - tau * vlag[jp]) / denom;
		for(i = 0; i <= jp; ++i){
			bmat[i+j*ndim] = bmat[i+j*ndim] + tempa * vlag[i] + tempb * w[i];
			if(i >= npt){
				bmat[jp+(i-npt)*ndim] = bmat[i+j*ndim];
			}
		}
	}
}

static inline int iabs(int i){ return (i < 0) ? -i : i; }
static inline int imax(int i, int j){ return (i > j) ? i : j; }

/*     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT, */
/*       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as */
/*       the corresponding arguments of BOBYQB on the entry to RESCUE. */
/*     NF is maintained as the number of calls of CALFUN so far, except that */
/*       NF is set to -1 if the value of MAXFUN prevents further progress. */
/*     KOPT is maintained so that FVAL(KOPT) is the least calculated function */
/*       value. Its correct value must be given on entry. It is updated if a */
/*       new least function value is found, but the corresponding changes to */
/*       XOPT and GOPT have to be made later by the calling program. */
/*     DELTA is the current trust region radius. */
/*     VLAG is a working space vector that will be used for the values of the */
/*       provisional Lagrange functions at each of the interpolation points. */
/*       They are part of a product that requires VLAG to be of length NDIM. */
/*     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and */
/*       PTSAUX(2,J) specify the two positions of provisional interpolation */
/*       points when a nonzero step is taken along e_J (the J-th coordinate */
/*       direction) through XBASE+XOPT, as specified below. Usually these */
/*       steps have length DELTA, but other lengths are chosen if necessary */
/*       in order to satisfy the given bounds on the variables. */
/*     PTSID is also a working space array. It has NPT components that denote */
/*       provisional new positions of the original interpolation points, in */
/*       case changes are needed to restore the linear independence of the */
/*       interpolation conditions. The K-th point is a candidate for change */
/*       if and only if PTSID(K) is nonzero. In this case let p and q be the */
/*       int parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p */
/*       and q are both positive, the step from XBASE+XOPT to the new K-th */
/*       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise */
/*       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or */
/*       p=0, respectively. */
/*     The first NDIM+NPT elements of the array W are used for working space. */
/*     The final elements of BMAT and ZMAT are set in a well-conditioned way */
/*       to the values that are appropriate for the new interpolation points. */
/*     The elements of GOPT, HQ and PQ are also revised to the values that are */
/*       appropriate to the final quadratic model. */
void rescue_(
	int n, int npt,
	double *xl, double *xu,
	double (*func)(void *data, int n, const double *x), void *data,
	int iprint, int maxfun, double *xbase, 
	double *xpt, double *fval, double *xopt, double *gopt,
	 double *hq, double *pq, double *bmat, double *zmat, 
	int ndim, double *sl, double *su, int *nf, 
	double delta, int *kopt, double *vlag, double *
	ptsaux, double *ptsid, double *w
){
    /* System generated locals */
    int xpt_offset, bmat_offset, zmat_offset, i__1, i__2, i__3;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    double f;
    int i__, j, k, ih, jp, ip, iq, np, iw;
    double xp, xq, den;
    int ihp;
    double one;
    int ihq, jpn, kpt;
    double sum, diff, half, beta;
    int kold;
    double winc;
    int nrem, knew;
    double temp, bsum;
    int nptm;
    double zero, hdiag, fbase, sfrac, denom, vquad, sumpq;
	    
	extern void update_(int, int, double *, double *, 
	    int, double *, double, double, int, double *);
    double dsqmin, distsq, vlmxsq;


    /* Parameter adjustments */
    zmat_offset = 1 + npt;
    zmat -= zmat_offset;
    xpt_offset = 1 + npt;
    xpt -= xpt_offset;
    --xl;
    --xu;
    --xbase;
    --fval;
    --xopt;
    --gopt;
    --hq;
    --pq;
    bmat_offset = 1 + ndim;
    bmat -= bmat_offset;
    --sl;
    --su;
    --vlag;
    ptsaux -= 3;
    --ptsid;
    --w;

    /* Function Body */
    half = .5;
    one = 1.;
    zero = 0.;
    np = n + 1;
    sfrac = half / (double) np;
    nptm = npt - np;

/*     Shift the interpolation points so that XOPT becomes the origin, and set */
/*     the elements of ZMAT to zero. The value of SUMPQ is required in the */
/*     updating of HQ below. The squares of the distances from XOPT to the */
/*     other interpolation points are set at the end of W. Increments of WINC */
/*     may be added later to these squares to balance the consideration of */
/*     the choice of point that is going to become current. */

    sumpq = zero;
    winc = zero;
    for (k = 1; k <= npt; ++k) {
		distsq = zero;
		for (j = 1; j <= n; ++j) {
			xpt[k + j * npt] -= xopt[j];
			distsq += xpt[k + j * npt] * xpt[k + j * npt];
		}
		sumpq += pq[k];
		w[ndim + k] = distsq;
		if(distsq > winc){ winc = distsq; }
		for (j = 1; j <= nptm; ++j) {
			zmat[k + j * npt] = zero;
		}
    }

/*     Update HQ so that HQ and PQ define the second derivatives of the model */
/*     after XBASE has been shifted to the trust region centre. */

    ih = 0;
    for (j = 1; j <= n; ++j) {
		w[j] = half * sumpq * xopt[j];
		for (k = 1; k <= npt; ++k) {
			w[j] += pq[k] * xpt[k + j * npt];
		}
		for (i__ = 1; i__ <= j; ++i__) {
			++ih;
			hq[ih] = hq[ih] + w[i__] * xopt[j] + w[j] * xopt[i__];
		}
    }

/*     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and */
/*     also set the elements of PTSAUX. */

    for (j = 1; j <= n; ++j) {
		xbase[j] += xopt[j];
		sl[j] -= xopt[j];
		su[j] -= xopt[j];
		xopt[j] = zero;
		ptsaux[2*j+1] = delta;
		if(su[j] < ptsaux[2*j+1]){ ptsaux[2*j+1] = su[j]; }
		ptsaux[2*j+2] = -delta;
		if(sl[j] > ptsaux[2*j+2]){ ptsaux[2*j+2] = sl[j]; }
		if (ptsaux[(2*j) + 1] + ptsaux[(2*j) + 2] < zero) {
			temp = ptsaux[(2*j) + 1];
			ptsaux[(2*j) + 1] = ptsaux[(2*j) + 2];
			ptsaux[(2*j) + 2] = temp;
		}
		if(fabs(ptsaux[2*j+2]) < half * fabs(ptsaux[2*j+1])){
			ptsaux[2*j+2] = half * ptsaux[2*j+1];
		}
		for (i__ = 1; i__ <= ndim; ++i__) {
			bmat[i__ + j * ndim] = zero;
		}
    }
    fbase = fval[*kopt];

/*     Set the identifiers of the artificial interpolation points that are */
/*     along a coordinate direction from XOPT, and set the corresponding */
/*     nonzero elements of BMAT and ZMAT. */

    ptsid[1] = sfrac;
    for (j = 1; j <= n; ++j) {
		jp = j + 1;
		jpn = jp + n;
		ptsid[jp] = (double) j + sfrac;
		if (jpn <= npt) {
			ptsid[jpn] = (double) j / (double) np + sfrac;
			temp = one / (ptsaux[2*j+1] - ptsaux[2*j+2]);
			bmat[jp + j * ndim] = -temp + one / ptsaux[2*j+1];
			bmat[jpn + j * ndim] = temp + one / ptsaux[2*j+2];
			bmat[j * ndim + 1] = -bmat[jp + j * ndim] - bmat[jpn + j * ndim];
			zmat[j * npt + 1] = M_SQRT2 / fabs(ptsaux[2*j+1] * ptsaux[2*j+2]);
			zmat[jp + j * npt] = zmat[j * npt + 1] * ptsaux[2*j+2] * temp;
			zmat[jpn + j * npt] = -zmat[j * npt + 1] * ptsaux[2*j+1] * temp;
		} else {
			bmat[j * ndim + 1] = -one / ptsaux[(2*j) + 1];
			bmat[jp + j * ndim] = one / ptsaux[(2*j) + 1];
			bmat[j + npt + j * ndim] = -half * (ptsaux[2*j+1] * ptsaux[2*j+1]);
		}
    }

/* Set any remaining identifiers with their nonzero elements of ZMAT. */

    if (npt >= n + np) {
		for (k = 2*np; k <= npt; ++k) {
			iw = (int) (((double) (k - np) - half) / (double) (n));
			ip = k - np - iw * n;
			iq = ip + iw;
			if (iq > n) {
				iq -= n;
			}
			ptsid[k] = (double) ip + (double) iq / (double) np + 
				sfrac;
			temp = one / (ptsaux[(ip << 1) + 1] * ptsaux[(iq << 1) + 1]);
			zmat[(k - np) * npt + 1] = temp;
			zmat[ip + 1 + (k - np) * npt] = -temp;
			zmat[iq + 1 + (k - np) * npt] = -temp;
			zmat[k + (k - np) * npt] = temp;
		}
    }
    nrem = npt;
    kold = 1;
    knew = *kopt;

/* Reorder the provisional points in the way that exchanges PTSID(KOLD) */
/* with PTSID(KNEW). */

L80:
    for (j = 1; j <= n; ++j) {
		temp = bmat[kold + j * ndim];
		bmat[kold + j * ndim] = bmat[knew + j * ndim];
		bmat[knew + j * ndim] = temp;
    }
    for (j = 1; j <= nptm; ++j) {
		temp = zmat[kold + j * npt];
		zmat[kold + j * npt] = zmat[knew + j * npt];
		zmat[knew + j * npt] = temp;
    }
    ptsid[kold] = ptsid[knew];
    ptsid[knew] = zero;
    w[ndim + knew] = zero;
    --nrem;
    if (knew != *kopt) {
		temp = vlag[kold];
		vlag[kold] = vlag[knew];
		vlag[knew] = temp;

		/* Update the BMAT and ZMAT matrices so that the status of the KNEW-th */
		/* interpolation point can be changed from provisional to original. The */
		/* branch to label 350 occurs if all the original points are reinstated. */
		/* The nonnegative values of W(NDIM+K) are required in the search below. */

		update_(n, npt, &bmat[bmat_offset], &zmat[zmat_offset], ndim, &vlag[1], beta, denom, knew-1, &w[1]);
		if (nrem == 0) {
			return;
		}
		for (k = 1; k <= npt; ++k) {
			w[ndim + k] = fabs(w[ndim + k]);
		}
    }

/*     Pick the index KNEW of an original interpolation point that has not */
/*     yet replaced one of the provisional interpolation points, giving */
/*     attention to the closeness to XOPT and to previous tries with KNEW. */

L120:
    dsqmin = zero;
    for (k = 1; k <= npt; ++k) {
		if (w[ndim + k] > zero) {
			if (dsqmin == zero || w[ndim + k] < dsqmin) {
			knew = k;
			dsqmin = w[ndim + k];
			}
		}
    }
    if (dsqmin == zero) {
		goto L260;
    }

	/* Form the W-vector of the chosen original interpolation point. */

    for (j = 1; j <= n; ++j) {
		w[npt + j] = xpt[knew + j * npt];
    }
    for (k = 1; k <= npt; ++k) {
		sum = zero;
		if (k == *kopt) {
		} else if (ptsid[k] == zero) {
			for (j = 1; j <= n; ++j) {
				sum += w[npt + j] * xpt[k + j * npt];
			}
		} else {
			ip = (int) ptsid[k];
			if (ip > 0) {
				sum = w[npt + ip] * ptsaux[(ip << 1) + 1];
			}
			iq = (int) ((double) np * ptsid[k] - (double) (ip * np));
			if (iq > 0) {
				iw = 1;
				if (ip == 0) {
					iw = 2;
				}
				sum += w[npt + iq] * ptsaux[iw + (iq << 1)];
			}
		}
		w[k] = half * sum * sum;
    }

/* Calculate VLAG and BETA for the required updating of the H matrix if */
/* XPT(KNEW,.) is reinstated in the set of interpolation points. */

    for (k = 1; k <= npt; ++k) {
		sum = zero;
		for (j = 1; j <= n; ++j) {
			sum += bmat[k + j * ndim] * w[npt + j];
		}
		vlag[k] = sum;
    }
    beta = zero;
    for (j = 1; j <= nptm; ++j) {
		sum = zero;
		for (k = 1; k <= npt; ++k) {
			sum += zmat[k + j * npt] * w[k];
		}
		beta -= sum * sum;
		for (k = 1; k <= npt; ++k) {
			vlag[k] += sum * zmat[k + j * npt];
		}
    }
    bsum = zero;
    distsq = zero;
    for (j = 1; j <= n; ++j) {
		sum = zero;
		for (k = 1; k <= npt; ++k) {
			sum += bmat[k + j * ndim] * w[k];
		}
		jp = j + npt;
		bsum += sum * w[jp];
		i__2 = ndim;
		for (ip = npt + 1; ip <= i__2; ++ip) {
			sum += bmat[ip + j * ndim] * w[ip];
		}
		bsum += sum * w[jp];
		vlag[jp] = sum;
		distsq += xpt[knew + j * npt] * xpt[knew + j * npt];
    }
    beta = half * distsq * distsq + beta - bsum;
    vlag[*kopt] += one;

	/* KOLD is set to the index of the provisional interpolation point that is */
	/* going to be deleted to make way for the KNEW-th original interpolation */
	/* point. The choice of KOLD is governed by the avoidance of a small value */
	/* of the denominator in the updating calculation of UPDATE. */

    denom = zero;
    vlmxsq = zero;
    for (k = 1; k <= npt; ++k) {
		if (ptsid[k] != zero) {
			hdiag = zero;
			for (j = 1; j <= nptm; ++j) {
				hdiag += zmat[k + j * npt] * zmat[k + j * npt];
			}
			den = beta * hdiag + vlag[k] * vlag[k];
			if (den > denom) {
				kold = k;
				denom = den;
			}
		}
		if(vlag[k] * vlag[k] > vlmxsq){ vlmxsq = vlag[k] * vlag[k]; }
    }
    if (denom <= vlmxsq * .01) {
		w[ndim + knew] = -w[ndim + knew] - winc;
		goto L120;
    }
    goto L80;

	/* When label 260 is reached, all the final positions of the interpolation */
	/* points have been chosen although any changes have not been included yet */
	/* in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart */
	/* from the shift of XBASE, the updating of the quadratic model remains to */
	/* be done. The following cycle through the new interpolation points begins */
	/* by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero, */
	/* except that a RETURN occurs if MAXFUN prohibits another value of F. */

L260:
	for (kpt = 1; kpt <= npt; ++kpt) {
		if (ptsid[kpt] == zero) { continue; }
		if (*nf >= maxfun) {
			*nf = -1;
			break;
		}
		ih = 0;
		for (j = 1; j <= n; ++j) {
			w[j] = xpt[kpt + j * npt];
			xpt[kpt + j * npt] = zero;
			temp = pq[kpt] * w[j];
			for (i__ = 1; i__ <= j; ++i__) {
			++ih;
			hq[ih] += temp * w[i__];
			}
		}
		pq[kpt] = zero;
		ip = (int) ptsid[kpt];
		iq = (int) ((double) np * ptsid[kpt] - (double) (ip * np));
		if (ip > 0) {
			xp = ptsaux[2*ip+1];
			xpt[kpt + ip * npt] = xp;
		}
		if (iq > 0) {
			xq = ptsaux[2*iq+1];
			if (ip == 0) {
			xq = ptsaux[2*iq+2];
			}
			xpt[kpt + iq * npt] = xq;
		}

		/* Set VQUAD to the value of the current model at the new point. */

		vquad = fbase;
		if (ip > 0) {
			ihp = (ip + ip * ip) / 2;
			vquad += xp * (gopt[ip] + half * xp * hq[ihp]);
		}
		if (iq > 0) {
			ihq = (iq + iq * iq) / 2;
			vquad += xq * (gopt[iq] + half * xq * hq[ihq]);
			if (ip > 0) {
			iw = imax(ihp,ihq) - iabs(ip - iq);
			vquad += xp * xq * hq[iw];
			}
		}
		for (k = 1; k <= npt; ++k) {
			temp = zero;
			if (ip > 0) {
				temp += xp * xpt[k + ip * npt];
			}
			if (iq > 0) {
				temp += xq * xpt[k + iq * npt];
			}
			vquad += half * pq[k] * temp * temp;
		}

		/* Calculate F at the new interpolation point, and set DIFF to the factor */
		/* that is going to multiply the KPT-th Lagrange function when the model */
		/* is updated to provide interpolation to the new function value. */

		for (i__ = 1; i__ <= n; ++i__) {
			w[i__] = xbase[i__] + xpt[kpt + i__ * npt];
			if(xl[i__] > w[i__]){ w[i__] = xl[i__]; }
			if(xu[i__] < w[i__]){ w[i__] = xu[i__]; }
			if (xpt[kpt + i__ * npt] == sl[i__]) {
				w[i__] = xl[i__];
			}
			if (xpt[kpt + i__ * npt] == su[i__]) {
				w[i__] = xu[i__];
			}
		}
		++(*nf);
		f = func(data, n, &w[1]);
		if(iprint == 3) {
			printf("Function number %d    F = %g    The corresponding X is:\n", *nf, f);
			for (i__ = 1; i__ <= n; ++i__) {
				printf("\t%g\n", w[i__]);
			}
		}
		fval[kpt] = f;
		if (f < fval[*kopt]) {
			*kopt = kpt;
		}
		diff = f - vquad;

		/* Update the quadratic model. The RETURN from the subroutine occurs when */
		/* all the new interpolation points are included in the model. */

		for (i__ = 1; i__ <= n; ++i__) {
			gopt[i__] += diff * bmat[kpt + i__ * ndim];
		}
		for (k = 1; k <= npt; ++k) {
			sum = zero;
			for (j = 1; j <= nptm; ++j) {
				sum += zmat[k + j * npt] * zmat[kpt + j * npt];
			}
			temp = diff * sum;
			if (ptsid[k] == zero) {
				pq[k] += temp;
			} else {
				ip = (int) ptsid[k];
				iq = (int) ((double) np * ptsid[k] - (double) (ip 
					* np));
				ihq = (iq * iq + iq) / 2;
				if (ip == 0) {
					hq[ihq] += temp * (ptsaux[2*iq+2] * ptsaux[2*iq+2]);
				} else {
					ihp = (ip * ip + ip) / 2;
					hq[ihp] += temp * (ptsaux[2*ip+1] * ptsaux[2*ip+1]);
					if (iq > 0) {
						hq[ihq] += temp * (ptsaux[2*iq+1] * ptsaux[2*iq+1]);
						iw = imax(ihp,ihq) - iabs(iq - ip);
						hq[iw] += temp * ptsaux[2*ip+1] * ptsaux[2*iq+1];
					}
				}
			}
		}
		ptsid[kpt] = zero;
    }
}

/* The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the */
/*   same as the corresponding arguments in SUBROUTINE BOBYQA. */
/* The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU */
/*   are the same as the corresponding arguments in BOBYQB, the elements */
/*   of SL and SU being set in BOBYQA. */
/* GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but */
/*   it is set by PRELIM to the gradient of the quadratic model at XBASE. */
/*   If XOPT is nonzero, BOBYQB will change it to its usual value later. */
/* NF is maintaned as the number of calls of CALFUN so far. */
/* KOPT will be such that the least calculated value of F so far is at */
/*   the point XPT(KOPT,.)+XBASE in the space of the variables. */

/* SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
/* BMAT and ZMAT for the first iteration, and it maintains the values of */
/* NF and KOPT. The vector X is also changed by PRELIM. */

void prelim_(
	int n, int npt,
	double *x, double *xl, double *xu,
	double (*func)(void *data, int n, const double *x), void *data,
	double rhobeg, int iprint, 
	int maxfun, double *xbase, double *xpt, double *fval,
	 double *gopt, double *hq, double *pq, double *bmat, 
	double *zmat, int ndim, double *sl, double *su, 
	int *nf, int *kopt
){
	/* System generated locals */
	int xpt_offset, bmat_offset, zmat_offset;
	double d__1, d__2, d__3, d__4;

	/* Local variables */
	double f;
	int i, j, k, ih, nfm;
	int nfx, ipt, jpt;
	double fbeg, diff, temp, stepa, stepb;
	int itemp;

	static const double half = .5;
	static const double one = 1.;
	static const double two = 2.;
	static const double zero = 0.;
	const double rhosq = rhobeg * rhobeg;
	const double recip = one / rhosq;
	const int np = n + 1;

/*     Set XBASE to the initial vector of variables, and set the initial */
/*     elements of XPT, BMAT, HQ, PQ and ZMAT to zero. */

    for (j = 0; j < n; ++j) {
		xbase[j] = x[j];
		for(k = 0; k < npt; ++k) {
			xpt[k + j * npt] = zero;
		}
		for (i = 0; i < ndim; ++i) {
			bmat[i + j * ndim] = zero;
		}
    }
    for (ih = 0; ih < n * np / 2; ++ih) {
		hq[ih] = zero;
    }
    for (k = 0; k < npt; ++k) {
		pq[k] = zero;
	}
	for (j = 0; j < npt - np; ++j) {
		for (k = 0; k < npt; ++k) {
			zmat[k + j * npt] = zero;
		}
    }

    /* Parameter adjustments */
    zmat_offset = 1 + npt;
    zmat -= zmat_offset;
    xpt_offset = 1 + npt;
    xpt -= xpt_offset;
    --x;
    --xl;
    --xu;
    --xbase;
    --fval;
    --gopt;
    --hq;
    --pq;
    bmat_offset = 1 + ndim;
    bmat -= bmat_offset;
    --sl;
    --su;
    
/*     Begin the initialization procedure. NF becomes one more than the number */
/*     of function values so far. The coordinates of the displacement of the */
/*     next initial interpolation point from XBASE are set in XPT(NF+1,.). */

    *nf = 0;
L50:
    nfm = *nf;
    nfx = *nf - n;
    ++(*nf);
    if (nfm <= 2*n) {
		if (nfm >= 1 && nfm <= n) {
			stepa = rhobeg;
			if (su[nfm] == zero) {
				stepa = -stepa;
			}
			xpt[*nf + nfm * npt] = stepa;
		} else if (nfm > n) {
			stepa = xpt[*nf - n + nfx * npt];
			stepb = -rhobeg;
			if (sl[nfx] == zero) {
				stepb = two * rhobeg;
				if(su[nfx] < stepb){ stepb = su[nfx]; }
			}
			if (su[nfx] == zero) {
				stepb = -two * rhobeg;
				if(sl[nfx] > stepb){ stepb = sl[nfx]; }
			}
			xpt[*nf + nfx * npt] = stepb;
		}
    } else {
		itemp = (nfm - np) / n;
		jpt = nfm - itemp * n - n;
		ipt = jpt + itemp;
		if (ipt > n) {
			itemp = jpt;
			jpt = ipt - n;
			ipt = itemp;
		}
		xpt[*nf + ipt * npt] = xpt[ipt + 1 + ipt * npt];
		xpt[*nf + jpt * npt] = xpt[jpt + 1 + jpt * npt];
    }

/* Calculate the next value of F. The least function value so far and */
/* its index are required. */

    for (j = 1; j <= n; ++j) {
		x[j] = xbase[j] + xpt[*nf + j * npt];
		if(xl[j] > x[j]){ x[j] = xl[j]; }
		if(xu[j] < x[j]){ x[j] = xu[j]; }
		if(xpt[*nf + j * npt] == sl[j]) {
			x[j] = xl[j];
		}
		if (xpt[*nf + j * npt] == su[j]) {
			x[j] = xu[j];
		}
    }
    f = func(data, n, &x[1]);
    if (iprint == 3) {
		printf("Function number %d    F = %g    The corresponding X is:\n", *nf, f);
		for (i = 1; i <= n; ++i) {
			printf("\t%g\n", x[i]);
		}
    }
    fval[*nf] = f;
    if (*nf == 1) {
		fbeg = f;
		*kopt = 1;
    } else if (f < fval[*kopt]) {
		*kopt = *nf;
    }

	/* Set the nonzero initial elements of BMAT and the quadratic model in the */
	/* cases when NF is at most 2n+1. If NF exceeds N+1, then the positions */
	/* of the NF-th and (NF-N)-th interpolation points may be switched, in */
	/* order that the function value at the first of them contributes to the */
	/* off-diagonal second derivative terms of the initial quadratic model. */

    if (*nf <= 2*n + 1) {
		if (*nf >= 2 && *nf <= n + 1) {
			gopt[nfm] = (f - fbeg) / stepa;
			if (npt < *nf + n) {
				bmat[nfm * ndim + 1] = -one / stepa;
				bmat[*nf + nfm * ndim] = one / stepa;
				bmat[npt + nfm + nfm * ndim] = -half * rhosq;
			}
		} else if (*nf >= n + 2) {
			ih = nfx * (nfx + 1) / 2;
			temp = (f - fbeg) / stepb;
			diff = stepb - stepa;
			hq[ih] = two * (temp - gopt[nfx]) / diff;
			gopt[nfx] = (gopt[nfx] * stepb - temp * stepa) / diff;
			if (stepa * stepb < zero) {
				if (f < fval[*nf - n]) {
					fval[*nf] = fval[*nf - n];
					fval[*nf - n] = f;
					if (*kopt == *nf) {
						*kopt = *nf - n;
					}
					xpt[*nf - n + nfx * npt] = stepb;
					xpt[*nf + nfx * npt] = stepa;
				}
			}
			bmat[nfx * ndim + 1] = -(stepa + stepb) / (stepa * stepb);
			bmat[*nf + nfx * ndim] = -half / xpt[*nf - n + nfx * npt];
			bmat[*nf - n + nfx * ndim] = -bmat[nfx * ndim + 1] - bmat[*nf + nfx * ndim];
			zmat[nfx * npt + 1] = sqrt(two) / (stepa * stepb);
			zmat[*nf + nfx * npt] = sqrt(half) / rhosq;
			zmat[*nf - n + nfx * npt] = -zmat[nfx * npt + 1] - zmat[*nf + nfx * npt];
		}

		/* Set the off-diagonal second derivatives of the Lagrange functions and */
		/* the initial quadratic model. */

    } else {
		ih = ipt * (ipt - 1) / 2 + jpt;
		zmat[nfx * npt + 1] = recip;
		zmat[*nf + nfx * npt] = recip;
		zmat[ipt + 1 + nfx * npt] = -recip;
		zmat[jpt + 1 + nfx * npt] = -recip;
		temp = xpt[*nf + ipt * npt] * xpt[*nf + jpt * npt];
		hq[ih] = (fbeg - fval[ipt + 1] - fval[jpt + 1] + f) / temp;
    }
    if (*nf < npt && *nf < maxfun) {
		goto L50;
    }
}

void bobyqb_(
	int n, int npt,
	double *x, double *xl, double *xu,
	double (*func)(void *data, int n, const double *x), void *data,
	double rhobeg, double rhoend, int iprint, int maxfun, double *xbase, 
	double *xpt, double *fval, double *xopt, double *gopt,
	double *hq, double *pq, double *bmat, double *zmat, 
	int ndim, double *sl, double *su, double *xnew, 
	double *xalt, double *d__, double *vlag, double *w
){
    /* System generated locals */
    int xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1, 
	    zmat_offset, i__1, i__2, i__3;
    double d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    double f;
    int i__, j, k, ih, nf, jj, nh, ip, jp;
    double dx;
    int np;
    double den, one, ten, dsq, rho, sum, two, diff, half, beta, 
	    gisq;
    int knew;
    double temp, suma, sumb, bsum, fopt;
    int kopt, nptm;
    double zero, curv;
    int ksav;
    double gqsq, dist, sumw, sumz, diffa, diffb, diffc, hdiag;
    int kbase;
    double alpha, delta, adelt, denom, fsave, bdtol, delsq;
    int nresc, nfsav;
    double ratio, dnorm, vquad, pqold, tenth;
    int itest;
    double sumpq, scaden;
    double errbig, cauchy, fracsq, biglsq, densav;
    extern void update_(int, int, double *, 
	    double *, int, double *, double, double,
	     int, double *);
    double bdtest;
    extern void rescue_(
		int n, int npt,
		double *xl, double *xu,
		double (*func)(void *data, int n, const double *x), void *data,
		int iprint, int maxfun, double *xbase, 
		double *xpt, double *fval, double *xopt, double *gopt,
		 double *hq, double *pq, double *bmat, double *zmat, 
		int ndim, double *sl, double *su, int *nf, 
		double delta, int *kopt, double *vlag, double *
		ptsaux, double *ptsid, double *w
	);
	extern void prelim_(
		int n, int npt,
		double *x, double *xl, double *xu,
		double (*func)(void *data, int n, const double *x), void *data,
		double rhobeg, int iprint, 
		int maxfun, double *xbase, double *xpt, double *fval,
		 double *gopt, double *hq, double *pq, double *bmat, 
		double *zmat, int ndim, double *sl, double *su, 
		int *nf, int *kopt
	);
    double crvmin, frhosq;
    extern void altmov_(int, int, double *, 
	    double *, double *, double *, int, double *, double *,
		int, int, double,
		double *, double *, double *, double *, double *, double *, double *);
    double distsq;
    extern void trsbox_(int, int, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *);
    int ntrits;
    double xoptsq;


/*     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN */
/*       are identical to the corresponding arguments in SUBROUTINE BOBYQA. */
/*     XBASE holds a shift of origin that should reduce the contributions */
/*       from rounding errors to values of the model and Lagrange functions. */
/*     XPT is a two-dimensional array that holds the coordinates of the */
/*       interpolation points relative to XBASE. */
/*     FVAL holds the values of F at the interpolation points. */
/*     XOPT is set to the displacement from XBASE of the trust region centre. */
/*     GOPT holds the gradient of the quadratic model at XBASE+XOPT. */
/*     HQ holds the explicit second derivatives of the quadratic model. */
/*     PQ contains the parameters of the implicit second derivatives of the */
/*       quadratic model. */
/*     BMAT holds the last N columns of H. */
/*     ZMAT holds the factorization of the leading NPT by NPT submatrix of H, */
/*       this factorization being ZMAT times ZMAT^T, which provides both the */
/*       correct rank and positive semi-definiteness. */
/*     NDIM is the first dimension of BMAT and has the value NPT+N. */
/*     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively. */
/*       All the components of every XOPT are going to satisfy the bounds */
/*       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when */
/*       XOPT is on a constraint boundary. */
/*     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the */
/*       vector of variables for the next call of CALFUN. XNEW also satisfies */
/*       the SL and SU constraints in the way that has just been mentioned. */
/*     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW */
/*       in order to increase the denominator in the updating of UPDATE. */
/*     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT. */
/*     VLAG contains the values of the Lagrange functions at a new point X. */
/*       They are part of a product that requires VLAG to be of length NDIM. */
/*     W is a one-dimensional array that is used for working space. Its length */
/*       must be at least 3*NDIM = 3*(NPT+N). */

/*     Set some constants. */

    /* Parameter adjustments */
    zmat_dim1 = npt;
    zmat_offset = 1 + zmat_dim1;
    zmat -= zmat_offset;
    xpt_dim1 = npt;
    xpt_offset = 1 + xpt_dim1;
    xpt -= xpt_offset;
    --x;
    --xl;
    --xu;
    --xbase;
    --fval;
    --xopt;
    --gopt;
    --hq;
    --pq;
    bmat_dim1 = ndim;
    bmat_offset = 1 + bmat_dim1;
    bmat -= bmat_offset;
    --sl;
    --su;
    --xnew;
    --xalt;
    --d__;
    --vlag;
    --w;

    /* Function Body */
    half = .5;
    one = 1.;
    ten = 10.;
    tenth = .1;
    two = 2.;
    zero = 0.;
    np = n + 1;
    nptm = npt - np;
    nh = n * np / 2;

/*     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
/*     BMAT and ZMAT for the first iteration, with the corresponding values of */
/*     of NF and KOPT, which are the number of calls of CALFUN so far and the */
/*     index of the interpolation point at the trust region centre. Then the */
/*     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is */
/*     less than NPT. GOPT will be updated if KOPT is different from KBASE. */

    prelim_(n, npt, &x[1], &xl[1], &xu[1], func, data, rhobeg, iprint, maxfun, &xbase[1], 
	    &xpt[xpt_offset], &fval[1], &gopt[1], &hq[1], &pq[1], &bmat[bmat_offset], &zmat[zmat_offset], ndim, &sl[1], &su[1], &nf, &kopt);
    xoptsq = zero;
    for (i__ = 1; i__ <= n; ++i__) {
		xopt[i__] = xpt[kopt + i__ * xpt_dim1];
	/* L10: */
	/* Computing 2nd power */
		d__1 = xopt[i__];
		xoptsq += d__1 * d__1;
    }
    fsave = fval[1];
    if (nf < npt) {
		if (iprint > 0) {
			printf("Return from BOBYQA because CALFUN has been called MAXFUN times.\n");
		}
		goto L720;
    }
    kbase = 1;

/*     Complete the settings that are required for the iterative procedure. */

    rho = rhobeg;
    delta = rho;
    nresc = nf;
    ntrits = 0;
    diffa = zero;
    diffb = zero;
    itest = 0;
    nfsav = nf;

/*     Update GOPT if necessary before the first iteration and after each */
/*     call of RESCUE that makes a call of CALFUN. */

L20:
    if (kopt != kbase) {
	ih = 0;
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++ih;
		if (i__ < j) {
		    gopt[j] += hq[ih] * xopt[i__];
		}
/* L30: */
		gopt[i__] += hq[ih] * xopt[j];
	    }
	}
	if (nf > npt) {
	    i__2 = npt;
	    for (k = 1; k <= i__2; ++k) {
		temp = zero;
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
/* L40: */
		    temp += xpt[k + j * xpt_dim1] * xopt[j];
		}
		temp = pq[k] * temp;
		i__1 = n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
		    gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
		}
	    }
	}
    }

/*     Generate the next point in the trust region that provides a small value */
/*     of the quadratic model subject to the constraints on the variables. */
/*     The int NTRITS is set to the number "trust region" iterations that */
/*     have occurred since the last "alternative" iteration. If the length */
/*     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to */
/*     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW. */

L60:
    trsbox_(n, npt, &xpt[xpt_offset], &xopt[1], &gopt[1], &hq[1], &pq[1], &sl[
	    1], &su[1], delta, &xnew[1], &d__[1], &w[1], &w[np], &w[np + n],
	     &w[np + (n << 1)], &w[np + n * 3], &dsq, &crvmin);
    dnorm = sqrt(dsq);
    if(delta < dnorm){ dnorm = delta; }
    if (dnorm < half * rho) {
	ntrits = -1;
	distsq = (ten * rho) * (ten * rho);
	if (nf <= nfsav + 2) {
	    goto L650;
	}

/*     The following choice between labels 650 and 680 depends on whether or */
/*     not our work with the current RHO seems to be complete. Either RHO is */
/*     decreased or termination occurs if the errors in the quadratic model at */
/*     the last three interpolation points compare favourably with predictions */
/*     of likely improvements to the model within distance HALF*RHO of XOPT. */

	errbig = diffa;
	if(diffb > errbig){ errbig = diffb; }
	if(diffc > errbig){ errbig = diffc; }
	frhosq = rho * .125 * rho;
	if (crvmin > zero && errbig > frhosq * crvmin) {
	    goto L650;
	}
	bdtol = errbig / rho;
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    bdtest = bdtol;
	    if (xnew[j] == sl[j]) {
		bdtest = w[j];
	    }
	    if (xnew[j] == su[j]) {
		bdtest = -w[j];
	    }
	    if (bdtest < bdtol) {
		curv = hq[(j + j * j) / 2];
		i__2 = npt;
		for (k = 1; k <= i__2; ++k) {
		    curv += pq[k] * (xpt[k + j * xpt_dim1] * xpt[k + j * xpt_dim1]);
		}
		bdtest += half * curv * rho;
		if (bdtest < bdtol) {
		    goto L650;
		}
	    }
	}
	goto L680;
    }
    ++ntrits;

/*     Severe cancellation is likely to occur if XOPT is too far from XBASE. */
/*     If the following test holds, then XBASE is shifted so that XOPT becomes */
/*     zero. The appropriate changes are made to BMAT and to the second */
/*     derivatives of the current model, beginning with the changes to BMAT */
/*     that do not depend on ZMAT. VLAG is used temporarily for working space. */

L90:
    if (dsq <= xoptsq * .001) {
		fracsq = xoptsq * .25;
		sumpq = zero;
		for (k = 1; k <= npt; ++k) {
			sumpq += pq[k];
			sum = -half * xoptsq;
			for (i__ = 1; i__ <= n; ++i__) {
				sum += xpt[k + i__ * xpt_dim1] * xopt[i__];
			}
			w[npt + k] = sum;
			temp = fracsq - half * sum;
			for (i__ = 1; i__ <= n; ++i__) {
				w[i__] = bmat[k + i__ * bmat_dim1];
				vlag[i__] = sum * xpt[k + i__ * xpt_dim1] + temp * xopt[i__];
				ip = npt + i__;
				for (j = 1; j <= i__; ++j) {
					bmat[ip + j * bmat_dim1] = bmat[ip + j * bmat_dim1] + w[i__] * vlag[j] + vlag[i__] * w[j];
				}
			}
		}

	/* Then the revisions of BMAT that depend on ZMAT are calculated. */

		for (jj = 1; jj <= nptm; ++jj) {
			sumz = zero;
			sumw = zero;
			i__2 = npt;
			for (k = 1; k <= i__2; ++k) {
				sumz += zmat[k + jj * zmat_dim1];
				vlag[k] = w[npt + k] * zmat[k + jj * zmat_dim1];
				sumw += vlag[k];
			}
			for (j = 1; j <= n; ++j) {
				sum = (fracsq * sumz - half * sumw) * xopt[j];
				for (k = 1; k <= npt; ++k) {
					sum += vlag[k] * xpt[k + j * xpt_dim1];
				}
				w[j] = sum;
				for (k = 1; k <= npt; ++k) {
					bmat[k + j * bmat_dim1] += sum * zmat[k + jj * zmat_dim1];
				}
			}
			for (i__ = 1; i__ <= n; ++i__) {
				ip = i__ + npt;
				temp = w[i__];
				for (j = 1; j <= i__; ++j) {
					bmat[ip + j * bmat_dim1] += temp * w[j];
				}
			}
		}

	/*     The following instructions complete the shift, including the changes */
	/*     to the second derivative parameters of the quadratic model. */

		ih = 0;
		for (j = 1; j <= n; ++j) {
			w[j] = -half * sumpq * xopt[j];
			i__1 = npt;
			for (k = 1; k <= i__1; ++k) {
				w[j] += pq[k] * xpt[k + j * xpt_dim1];
				xpt[k + j * xpt_dim1] -= xopt[j];
			}
			for (i__ = 1; i__ <= j; ++i__) {
				++ih;
				hq[ih] = hq[ih] + w[i__] * xopt[j] + xopt[i__] * w[j];
				bmat[npt + i__ + j * bmat_dim1] = bmat[npt + j + i__ * bmat_dim1];
			}
		}
		for (i__ = 1; i__ <= n; ++i__) {
			xbase[i__] += xopt[i__];
			xnew[i__] -= xopt[i__];
			sl[i__] -= xopt[i__];
			su[i__] -= xopt[i__];
			xopt[i__] = zero;
		}
		xoptsq = zero;
    }
    if (ntrits == 0) {
		goto L210;
    }
    goto L230;

/*     XBASE is also moved to XOPT by a call of RESCUE. This calculation is */
/*     more expensive than the previous shift, because new matrices BMAT and */
/*     ZMAT are generated from scratch, which may include the replacement of */
/*     interpolation points whose positions seem to be causing near linear */
/*     dependence in the interpolation conditions. Therefore RESCUE is called */
/*     only if rounding errors have reduced by at least a factor of two the */
/*     denominator of the formula for updating the H matrix. It provides a */
/*     useful safeguard, but is not invoked in most applications of BOBYQA. */

L190:
    nfsav = nf;
    kbase = kopt;
    rescue_(n, npt, &xl[1], &xu[1], func, data, iprint, maxfun, &xbase[1], &xpt[
	    xpt_offset], &fval[1], &xopt[1], &gopt[1], &hq[1], &pq[1], &bmat[
	    bmat_offset], &zmat[zmat_offset], ndim, &sl[1], &su[1], &nf,
	    delta, &kopt, &vlag[1], &w[1], &w[n + np], &w[ndim + np]);

/*     XOPT is updated now in case the branch below to label 720 is taken. */
/*     Any updating of GOPT occurs after the branch below to label 20, which */
/*     leads to a trust region iteration as does the branch to label 60. */

    xoptsq = zero;
    if (kopt != kbase) {
		for (i__ = 1; i__ <= n; ++i__) {
			xopt[i__] = xpt[kopt + i__ * xpt_dim1];
			xoptsq += xopt[i__] * xopt[i__];
		}
    }
    if (nf < 0) {
		nf = maxfun;
		if (iprint > 0) {
			printf("Return from BOBYQA because CALFUN has been called MAXFUN times.\n");
		}
		goto L720;
    }
    nresc = nf;
    if (nfsav < nf) {
		nfsav = nf;
		goto L20;
    }
    if (ntrits > 0) {
		goto L60;
    }

	/* Pick two alternative vectors of variables, relative to XBASE, that */
	/* are suitable as new positions of the KNEW-th interpolation point. */
	/* Firstly, XNEW is set to the point on a line through XOPT and another */
	/* interpolation point that minimizes the predicted value of the next */
	/* denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL */
	/* and SU bounds. Secondly, XALT is set to the best feasible point on */
	/* a constrained version of the Cauchy step of the KNEW-th Lagrange */
	/* function, the corresponding value of the square of this function */
	/* being returned in CAUCHY. The choice between these alternatives is */
	/* going to be made when the denominator is calculated. */

L210:
    altmov_(n, npt, &xpt[xpt_offset], &xopt[1], &bmat[bmat_offset], &zmat[
	    zmat_offset], ndim, &sl[1], &su[1], kopt, knew, adelt, &xnew[1]
	    , &xalt[1], &alpha, &cauchy, &w[1], &w[np], &w[ndim + 1]);
    for (i__ = 1; i__ <= n; ++i__) {
		d__[i__] = xnew[i__] - xopt[i__];
    }

/*     Calculate VLAG and BETA for the current choice of D. The scalar */
/*     product of D with XPT(K,.) is going to be held in W(NPT+K) for */
/*     use when VQUAD is calculated. */

L230:
    for (k = 1; k <= npt; ++k) {
		suma = zero;
		sumb = zero;
		sum = zero;
		i__2 = n;
		for (j = 1; j <= i__2; ++j) {
			suma += xpt[k + j * xpt_dim1] * d__[j];
			sumb += xpt[k + j * xpt_dim1] * xopt[j];
			sum += bmat[k + j * bmat_dim1] * d__[j];
		}
		w[k] = suma * (half * suma + sumb);
		vlag[k] = sum;
		w[npt + k] = suma;
    }
    beta = zero;
    for (jj = 1; jj <= nptm; ++jj) {
		sum = zero;
		for (k = 1; k <= npt; ++k) {
			sum += zmat[k + jj * zmat_dim1] * w[k];
		}
		beta -= sum * sum;
		for (k = 1; k <= npt; ++k) {
			vlag[k] += sum * zmat[k + jj * zmat_dim1];
		}
    }
    dsq = zero;
    bsum = zero;
    dx = zero;
    for (j = 1; j <= n; ++j) {
		d__1 = d__[j];
		dsq += d__1 * d__1;
		sum = zero;
		for (k = 1; k <= npt; ++k) {
			sum += w[k] * bmat[k + j * bmat_dim1];
		}
		bsum += sum * d__[j];
		jp = npt + j;
		for (i__ = 1; i__ <= n; ++i__) {
			sum += bmat[jp + i__ * bmat_dim1] * d__[i__];
		}
		vlag[jp] = sum;
		bsum += sum * d__[j];
		dx += d__[j] * xopt[j];
    }
    beta = dx * dx + dsq * (xoptsq + dx + dx + half * dsq) + beta - bsum;
    vlag[kopt] += one;

/*     If NTRITS is zero, the denominator may be increased by replacing */
/*     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if */
/*     rounding errors have damaged the chosen denominator. */

	if (ntrits == 0) {
		denom = vlag[knew] * vlag[knew] + alpha * beta;
		if (denom < cauchy && cauchy > zero) {
			i__2 = n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			xnew[i__] = xalt[i__];
			d__[i__] = xnew[i__] - xopt[i__];
			}
			cauchy = zero;
			goto L230;
		}
		if (denom <= half * (vlag[knew] * vlag[knew])) {
			if (nf > nresc) {
				goto L190;
			}
			if (iprint > 0) {
				printf("Return from BOBYQA because of much cancellation in a denominator.\n");
			}
			goto L720;
		}

		/* Alternatively, if NTRITS is positive, then set KNEW to the index of */
		/* the next interpolation point to be deleted to make room for a trust */
		/* region step. Again RESCUE may be called if rounding errors have damaged */
		/* the chosen denominator, which is the reason for attempting to select */
		/* KNEW before calculating the next value of the objective function. */
    } else {
		delsq = delta * delta;
		scaden = zero;
		biglsq = zero;
		knew = 0;
		for (k = 1; k <= npt; ++k) {
			if (k == kopt) {
				goto L350;
			}
			hdiag = zero;
			for (jj = 1; jj <= nptm; ++jj) {
				hdiag += zmat[k + jj * zmat_dim1] * zmat[k + jj * zmat_dim1];
			}
			den = beta * hdiag + vlag[k] * vlag[k];
			distsq = zero;
			for (j = 1; j <= n; ++j) {
				distsq += (xpt[k + j * xpt_dim1] - xopt[j]) * (xpt[k + j * xpt_dim1] - xopt[j]);
			}
			{
				double d3 = distsq / delsq;
				temp = d3*d3;
				if(one > temp){ temp = one; }
			}
			if (temp * den > scaden) {
			scaden = temp * den;
			knew = k;
			denom = den;
			}
			{
				double d2 = temp * (vlag[k] * vlag[k]);
				if(d2 > biglsq){ biglsq = d2; }
			}
L350:
			;
		}
		if (scaden <= half * biglsq) {
			if (nf > nresc) {
				goto L190;
			}
			if (iprint > 0) {
				printf("Return from BOBYQA because of much cancellation in a denominator.\n");
			}
			goto L720;
		}
    }

	/*     Put the variables for the next calculation of the objective function */
	/*       in XNEW, with any adjustments for the bounds. */


	/*     Calculate the value of the objective function at XBASE+XNEW, unless */
	/*       the limit on the number of calculations of F has been reached. */

L360:
    for (i__ = 1; i__ <= n; ++i__) {
		x[i__] = xbase[i__] + xnew[i__];
		if(xl[i__] > x[i__]){ x[i__] = xl[i__]; }
		if(xu[i__] < x[i__]){ x[i__] = xu[i__]; }
		if (xnew[i__] == sl[i__]) {
			x[i__] = xl[i__];
		}
		if (xnew[i__] == su[i__]) {
			x[i__] = xu[i__];
		}
    }
    if (nf >= maxfun) {
		if (iprint > 0) {
			printf("Return from BOBYQA because CALFUN has been called MAXFUN times.\n");
		}
		goto L720;
    }
    ++nf;
    f = func(data, n, &x[1]);
    if (iprint == 3) {
		printf("Function number %d    F = %g\n    The corresponding X is:\n", nf, f);
		for (i__ = 1; i__ <= n; ++i__) {
			printf("\t%g\n", x[i__]);
		}
    }
    if (ntrits == -1) {
		fsave = f;
		goto L720;
    }

/*     Use the quadratic model to predict the change in F due to the step D, */
/*       and set DIFF to the error of this prediction. */

    fopt = fval[kopt];
    vquad = zero;
    ih = 0;
    for (j = 1; j <= n; ++j) {
		vquad += d__[j] * gopt[j];
		for (i__ = 1; i__ <= j; ++i__) {
			++ih;
			temp = d__[i__] * d__[j];
			if (i__ == j) {
			temp = half * temp;
			}
			vquad += hq[ih] * temp;
		}
    }
    for (k = 1; k <= npt; ++k) {
		vquad += half * pq[k] * (w[npt + k] * w[npt + k]);
    }
    diff = f - fopt - vquad;
    diffc = diffb;
    diffb = diffa;
    diffa = fabs(diff);
    if (dnorm > rho) {
		nfsav = nf;
    }

/*     Pick the next value of DELTA after a trust region step. */

	if (ntrits > 0) {
		if (vquad >= zero) {
			if (iprint > 0) {
				printf("Return from BOBYQA because a trust region step has failed to reduce Q.\n");
			}
			goto L720;
		}
		ratio = (f - fopt) / vquad;
		if (ratio <= tenth) {
			delta = half * delta;
			if(dnorm < delta){ delta = dnorm; }
		} else if (ratio <= .7) {
			delta = half * delta;
			if(dnorm > delta){ delta = dnorm; }
		} else {
			delta = half * delta;
			if(2*dnorm > delta){ delta = 2*dnorm; }
		}
		if (delta <= rho * 1.5) {
			delta = rho;
		}

	/*     Recalculate KNEW and DENOM if the new F is less than FOPT. */

		if (f < fopt) {
			ksav = knew;
			densav = denom;
			delsq = delta * delta;
			scaden = zero;
			biglsq = zero;
			knew = 0;
			i__1 = npt;
			for (k = 1; k <= i__1; ++k) {
				hdiag = zero;
				i__2 = nptm;
				for (jj = 1; jj <= i__2; ++jj) {
					hdiag += zmat[k + jj * zmat_dim1] * zmat[k + jj * zmat_dim1];
				}
				den = beta * hdiag + vlag[k] * vlag[k];
				distsq = zero;
				for (j = 1; j <= n; ++j) {
					distsq += (xpt[k + j * xpt_dim1] - xnew[j]) * (xpt[k + j * xpt_dim1] - xnew[j]);
				}
				{
					double d3 = distsq / delsq;
					temp = d3*d3;
					if(one > temp){ temp = one; }
				}
				if (temp * den > scaden) {
					scaden = temp * den;
					knew = k;
					denom = den;
				}
				{
					double d2 = temp * (vlag[k] * vlag[k]);
					if(d2 > biglsq){ biglsq = d2; }
				}
			}
			if (scaden <= half * biglsq) {
			knew = ksav;
			denom = densav;
			}
		}
    }

/*     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be */
/*     moved. Also update the second derivative terms of the model. */

    update_(n, npt, &bmat[bmat_offset], &zmat[zmat_offset], ndim, &vlag[1], 
	    beta, denom, knew-1, &w[1]);
    ih = 0;
    pqold = pq[knew];
    pq[knew] = zero;
    for (i__ = 1; i__ <= n; ++i__) {
		temp = pqold * xpt[knew + i__ * xpt_dim1];
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
			++ih;
			hq[ih] += temp * xpt[knew + j * xpt_dim1];
		}
    }
    for (jj = 1; jj <= nptm; ++jj) {
		temp = diff * zmat[knew + jj * zmat_dim1];
		for (k = 1; k <= npt; ++k) {
			pq[k] += temp * zmat[k + jj * zmat_dim1];
		}
    }

/*     Include the new interpolation point, and make the changes to GOPT at */
/*     the old XOPT that are caused by the updating of the quadratic model. */

    fval[knew] = f;
    for (i__ = 1; i__ <= n; ++i__) {
		xpt[knew + i__ * xpt_dim1] = xnew[i__];
		w[i__] = bmat[knew + i__ * bmat_dim1];
    }
    for (k = 1; k <= npt; ++k) {
		suma = zero;
		for (jj = 1; jj <= nptm; ++jj) {
			suma += zmat[knew + jj * zmat_dim1] * zmat[k + jj * zmat_dim1];
		}
		sumb = zero;
		for (j = 1; j <= n; ++j) {
			sumb += xpt[k + j * xpt_dim1] * xopt[j];
		}
		temp = suma * sumb;
		for (i__ = 1; i__ <= n; ++i__) {
			w[i__] += temp * xpt[k + i__ * xpt_dim1];
		}
    }
    for (i__ = 1; i__ <= n; ++i__) {
		gopt[i__] += diff * w[i__];
    }

/*     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT. */

    if (f < fopt) {
		kopt = knew;
		xoptsq = zero;
		ih = 0;
		for (j = 1; j <= n; ++j) {
			xopt[j] = xnew[j];
			xoptsq += xopt[j] * xopt[j];
			for (i__ = 1; i__ <= j; ++i__) {
				++ih;
				if (i__ < j) {
					gopt[j] += hq[ih] * d__[i__];
				}
				gopt[i__] += hq[ih] * d__[j];
			}
		}
		for (k = 1; k <= npt; ++k) {
			temp = zero;
			for (j = 1; j <= n; ++j) {
				temp += xpt[k + j * xpt_dim1] * d__[j];
			}
			temp = pq[k] * temp;
			for (i__ = 1; i__ <= n; ++i__) {
				gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
			}
		}
    }

/*     Calculate the parameters of the least Frobenius norm interpolant to */
/*     the current data, the gradient of this interpolant at XOPT being put */
/*     into VLAG(NPT+I), I=1,2,...,N. */

    if (ntrits > 0) {
		for (k = 1; k <= npt; ++k) {
			vlag[k] = fval[k] - fval[kopt];
			w[k] = zero;
		}
		for (j = 1; j <= nptm; ++j) {
			sum = zero;
			for (k = 1; k <= npt; ++k) {
				sum += zmat[k + j * zmat_dim1] * vlag[k];
			}
			for (k = 1; k <= npt; ++k) {
				w[k] += sum * zmat[k + j * zmat_dim1];
			}
		}
		for (k = 1; k <= npt; ++k) {
			sum = zero;
			for (j = 1; j <= n; ++j) {
				sum += xpt[k + j * xpt_dim1] * xopt[j];
			}
			w[k + npt] = w[k];
			w[k] = sum * w[k];
		}
		gqsq = zero;
		gisq = zero;
		for (i__ = 1; i__ <= n; ++i__) {
			sum = zero;
			for (k = 1; k <= npt; ++k) {
				sum = sum + bmat[k + i__ * bmat_dim1] * vlag[k] + xpt[k + i__ 
				* xpt_dim1] * w[k];
			}
			if (xopt[i__] == sl[i__]) {
				double d1 = gopt[i__];
				if(zero < d1){ d1 = zero; }
				gqsq += d1 * d1;
				d1 = sum;
				if(zero < sum){ sum = zero; }
				gisq += d1 * d1;
			} else if (xopt[i__] == su[i__]) {
				double d1 = gopt[i__];
				if(zero > d1){ d1 = zero; }
				gqsq += d1 * d1;
				d1 = sum;
				if(zero > sum){ sum = zero; }
				gisq += d1 * d1;
			} else {
				gqsq += gopt[i__] * gopt[i__];
				gisq += sum * sum;
			}
			vlag[npt + i__] = sum;
		}

	/*     Test whether to replace the new quadratic model by the least Frobenius */
	/*     norm interpolant, making the replacement if the test is satisfied. */

		++itest;
		if (gqsq < ten * gisq) {
			itest = 0;
		}
		if (itest >= 3) {
			i__1 = imax(npt,nh);
			for (i__ = 1; i__ <= i__1; ++i__) {
			if (i__ <= n) {
				gopt[i__] = vlag[npt + i__];
			}
			if (i__ <= npt) {
				pq[i__] = w[npt + i__];
			}
			if (i__ <= nh) {
				hq[i__] = zero;
			}
			itest = 0;
			}
		}
    }

/*     If a trust region step has provided a sufficient decrease in F, then */
/*     branch for another trust region calculation. The case NTRITS=0 occurs */
/*     when the new interpolation point was reached by an alternative step. */

    if (ntrits == 0) {
		goto L60;
    }
    if (f <= fopt + tenth * vquad) {
		goto L60;
    }

/*     Alternatively, find out if the interpolation points are close enough */
/*       to the best point so far. */
	{
		double d3 = two*delta;
		double d4 = ten*rho;
		double d1 = d3*d3;
		double d2 = d4*d4;
		distsq = (d1 > d2) ? d1 : d2;
	}
L650:
    knew = 0;
    for (k = 1; k <= npt; ++k) {
		sum = zero;
		for (j = 1; j <= n; ++j) {
			sum += (xpt[k + j * xpt_dim1] - xopt[j]) * (xpt[k + j * xpt_dim1] - xopt[j]);
		}
		if (sum > distsq) {
			knew = k;
			distsq = sum;
		}
    }

/*     If KNEW is positive, then ALTMOV finds alternative new positions for */
/*     the KNEW-th interpolation point within distance ADELT of XOPT. It is */
/*     reached via label 90. Otherwise, there is a branch to label 60 for */
/*     another trust region iteration, unless the calculations with the */
/*     current RHO are complete. */

    if (knew > 0) {
		dist = sqrt(distsq);
		if (ntrits == -1) {
			delta = tenth * delta;
			if(half * dist < delta){ delta = half * dist; }
			if (delta <= rho * 1.5) {
				delta = rho;
			}
		}
		ntrits = 0;
		adelt = tenth * dist;
		if(delta < adelt){ adelt = delta; }
		if(rho > adelt){ adelt = rho; }
		dsq = adelt * adelt;
		goto L90;
    }
    if (ntrits == -1) {
		goto L680;
    }
    if (ratio > zero) {
		goto L60;
	}
    if ((delta > dnorm ? delta : dnorm) > rho) {
		goto L60;
    }

/*     The calculations with the current value of RHO are complete. Pick the */
/*       next values of RHO and DELTA. */

L680:
    if (rho > rhoend) {
		delta = half * rho;
		ratio = rho / rhoend;
		if (ratio <= 16.) {
			rho = rhoend;
		} else if (ratio <= 250.) {
			rho = sqrt(ratio) * rhoend;
		} else {
			rho = tenth * rho;
		}
		delta = (delta > rho ? delta : rho);
		if (iprint >= 2) {
			printf("New RHO = %g, Number of function values = %d", rho, nf);
			printf("Least value of F = %g\n    The corresponding X is:\n", fval[kopt]);
			for (i__ = 1; i__ <= n; ++i__) {
				printf("\t%g\n", xbase[i__] + xopt[i__]);
			}
		}
		ntrits = 0;
		nfsav = nf;
		goto L60;
    }

/*     Return from the calculation, after another Newton-Raphson step, if */
/*       it is too short to have been tried before. */

    if (ntrits == -1) {
		goto L360;
    }
L720:
    if (fval[kopt] <= fsave) {
		for (i__ = 1; i__ <= n; ++i__) {
			x[i__] = xbase[i__] + xopt[i__];
			if(xl[i__] > x[i__]){ x[i__] = xl[i__]; }
			if(xu[i__] < x[i__]){ x[i__] = xu[i__]; }
			if (xopt[i__] == sl[i__]) {
				x[i__] = xl[i__];
			}
			if (xopt[i__] == su[i__]) {
				x[i__] = xu[i__];
			}
		}
		f = fval[kopt];
    }
    if (iprint >= 1) {
		printf("At the return from BOBYQA\nNumber of function values = %d\n", nf);
		printf("Least value of F = %g\n    The corresponding X is:\n", f);
		for (i__ = 1; i__ <= n; ++i__) {
			printf("\t%g\n", x[i__]);
		}
    }
}


void BOBYQA(
	int n, int npt,
	double *x, double *xl, double *xu,
	double (*func)(void *data, int n, const double *x), void *data,
	double rhobeg, double rhoend, int iprint, int maxfun, double *w
){
    /* Local variables */
    int j, id, np, iw, igo, ihq, ixb, ixa, ifv, isl, jsl, ipq, ivl,
	     ixn, ixo, ixp, isu, jsu, ndim;
    int ibmat, izmat;
    extern void bobyqb_(
		int n, int npt,
		double *x, double *xl, double *xu,
		double (*func)(void *data, int n, const double *x), void *data,
		double rhobeg, double rhoend, int iprint, int maxfun, double *xbase, 
		double *xpt, double *fval, double *xopt, double *gopt,
		double *hq, double *pq, double *bmat, double *zmat, 
		int ndim, double *sl, double *su, double *xnew, 
		double *xalt, double *d__, double *vlag, double *w
	);

/*     This subroutine seeks the least value of a function of many variables, */
/*     by applying a trust region method that forms quadratic models by */
/*     interpolation. There is usually some freedom in the interpolation */
/*     conditions, which is taken up by minimizing the Frobenius norm of */
/*     the change to the second derivative of the model, beginning with the */
/*     zero matrix. The values of the variables are constrained by upper and */
/*     lower bounds. The arguments of the subroutine are as follows. */

/*     N must be set to the number of variables and must be at least two. */
/*     NPT is the number of interpolation conditions. Its value must be in */
/*       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not */
/*       recommended. */
/*     Initial values of the variables must be set in X(1),X(2),...,X(N). They */
/*       will be changed to the values that give the least calculated F. */
/*     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper */
/*       bounds, respectively, on X(I). The construction of quadratic models */
/*       requires XL(I) to be strictly less than XU(I) for each I. Further, */
/*       the contribution to a model from changes to the I-th variable is */
/*       damaged severely by rounding errors if XU(I)-XL(I) is too small. */
/*     RHOBEG and RHOEND must be set to the initial and final values of a trust */
/*       region radius, so both must be positive with RHOEND no greater than */
/*       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest */
/*       expected change to a variable, while RHOEND should indicate the */
/*       accuracy that is required in the final values of the variables. An */
/*       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N, */
/*       is less than 2*RHOBEG. */
/*     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the */
/*       amount of printing. Specifically, there is no output if IPRINT=0 and */
/*       there is output only at the return if IPRINT=1. Otherwise, each new */
/*       value of RHO is printed, with the best vector of variables so far and */
/*       the corresponding value of the objective function. Further, each new */
/*       value of F with its variables are output if IPRINT=3. */
/*     MAXFUN must be set to an upper bound on the number of calls of CALFUN. */
/*     The array W will be used for working space. Its length must be at least */
/*       (NPT+5)*(NPT+N)+3*N*(N+5)/2. */

/*     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set */
/*     F to the value of the objective function for the current values of the */
/*     variables X(1),X(2),...,X(N), which are generated automatically in a */
/*     way that satisfies the bounds given in XL and XU. */

/*     Return if the value of NPT is unacceptable. */

    np = n + 1;
    if (npt < n + 2 || npt > (n + 2) * np / 2) {
		printf("Return from BOBYQA because NPT is not in the required interval\n");
		return;
    }

	/* Partition the working space array, so that different parts of it can */
	/* be treated separately during the calculation of BOBYQB. The partition */
	/* requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the */
	/* space that is taken by the last array in the argument list of BOBYQB. */

    ndim = npt + n;
    ixb = 0;
    ixp = ixb + n;
    ifv = ixp + n * npt;
    ixo = ifv + npt;
    igo = ixo + n;
    ihq = igo + n;
    ipq = ihq + n * np / 2;
    ibmat = ipq + npt;
    izmat = ibmat + ndim * n;
    isl = izmat + npt * (npt - np);
    isu = isl + n;
    ixn = isu + n;
    ixa = ixn + n;
    id = ixa + n;
    ivl = id + n;
    iw = ivl + ndim;

	/* Return if there is insufficient space between the bounds. Modify the */
	/* initial X if necessary in order to avoid conflicts between the bounds */
	/* and the construction of the first quadratic model. The lower and upper */
	/* bounds on moves from the updated X are set now, in the ISL and ISU */
	/* partitions of W, in order to provide useful and exact information about */
	/* components of X that become within distance RHOBEG from their bounds. */

    for (j = 0; j < n; ++j) {
		double temp = xu[j] - xl[j];
		if (temp < 2*rhobeg) {
			printf("Return from BOBYQA because one of the differences XU(I)-XL(I) is less than 2*RHOBEG.\n");
			return;
		}
		jsl = isl + j;
		jsu = jsl + n;
		w[jsl] = xl[j] - x[j];
		w[jsu] = xu[j] - x[j];
		if (w[jsl] >= -rhobeg) {
			if (w[jsl] >= 0) {
				x[j] = xl[j];
				w[jsl] = 0;
				w[jsu] = temp;
			} else {
				x[j] = xl[j] + rhobeg;
				w[jsl] = -rhobeg;
				w[jsu] = xu[j] - x[j];
				if(rhobeg > w[jsu]){ w[jsu] = rhobeg; }
			}
		} else if (w[jsu] <= rhobeg) {
			if (w[jsu] <= 0) {
				x[j] = xu[j];
				w[jsl] = -temp;
				w[jsu] = 0;
			} else {
				x[j] = xu[j] - rhobeg;
				w[jsl] = xl[j] - x[j];
				if(-rhobeg < w[jsl]){ w[jsl] = -rhobeg; }
				w[jsu] = rhobeg;
			}
		}
    }

    bobyqb_(
		n, npt, x, xl, xu, func, data, rhobeg, rhoend, iprint, maxfun,
		&w[ixb], &w[ixp], &w[ifv], &w[ixo], &w[igo], &w[ihq], &w[ipq],
		&w[ibmat], &w[izmat], ndim, &w[isl], &w[isu], &w[ixn],
		&w[ixa], &w[id], &w[ivl], &w[iw]
	);
}
