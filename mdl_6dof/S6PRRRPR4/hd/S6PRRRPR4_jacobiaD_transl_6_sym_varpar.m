% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:20
% EndTime: 2019-02-26 20:12:21
% DurationCPUTime: 0.42s
% Computational Cost: add. (675->79), mult. (1208->127), div. (0->0), fcn. (1207->14), ass. (0->65)
t350 = sin(qJ(3));
t353 = cos(qJ(3));
t346 = qJ(4) + pkin(12);
t332 = pkin(5) * cos(t346) + cos(qJ(4)) * pkin(4);
t342 = qJ(6) + t346;
t338 = sin(t342);
t339 = cos(t342);
t369 = t339 * r_i_i_C(1) - t338 * r_i_i_C(2);
t365 = pkin(3) + t332 + t369;
t391 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(9);
t331 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t346);
t328 = t331 * qJD(4);
t345 = qJD(4) + qJD(6);
t368 = t338 * r_i_i_C(1) + t339 * r_i_i_C(2);
t393 = t368 * t345 + t328;
t355 = t393 * t353 + (t365 * t350 - t391 * t353) * qJD(3) - t350 * qJD(5);
t347 = sin(pkin(11));
t351 = sin(qJ(2));
t354 = cos(qJ(2));
t389 = cos(pkin(11));
t390 = cos(pkin(6));
t366 = t390 * t389;
t323 = t347 * t354 + t351 * t366;
t348 = sin(pkin(6));
t377 = t348 * t389;
t313 = t323 * t353 - t350 * t377;
t378 = t347 * t390;
t325 = -t351 * t378 + t389 * t354;
t387 = t348 * t350;
t386 = t348 * t353;
t385 = t348 * t354;
t319 = t323 * qJD(2);
t373 = t313 * t345 - t319;
t364 = t354 * t366;
t381 = qJD(2) * t351;
t318 = -qJD(2) * t364 + t347 * t381;
t360 = -t323 * t350 - t353 * t377;
t309 = t360 * qJD(3) - t318 * t353;
t322 = t347 * t351 - t364;
t375 = -t322 * t345 - t309;
t384 = (t375 * t338 - t373 * t339) * r_i_i_C(1) + (t373 * t338 + t375 * t339) * r_i_i_C(2);
t315 = t325 * t353 + t347 * t387;
t321 = t325 * qJD(2);
t372 = t315 * t345 - t321;
t324 = t389 * t351 + t354 * t378;
t320 = t324 * qJD(2);
t363 = -t325 * t350 + t347 * t386;
t311 = t363 * qJD(3) - t320 * t353;
t374 = -t324 * t345 - t311;
t383 = (t374 * t338 - t372 * t339) * r_i_i_C(1) + (t372 * t338 + t374 * t339) * r_i_i_C(2);
t327 = t390 * t350 + t351 * t386;
t362 = -t327 * t345 + t348 * t381;
t359 = -t351 * t387 + t390 * t353;
t379 = qJD(2) * t385;
t317 = t359 * qJD(3) + t353 * t379;
t367 = t345 * t385 - t317;
t382 = (t367 * t338 + t362 * t339) * r_i_i_C(1) + (-t362 * t338 + t367 * t339) * r_i_i_C(2);
t361 = pkin(8) + t331 + t368;
t329 = t332 * qJD(4);
t357 = t369 * t345 + t329;
t356 = -t391 * t350 - t365 * t353 - pkin(2);
t316 = t327 * qJD(3) + t350 * t379;
t310 = t315 * qJD(3) - t320 * t350;
t308 = t313 * qJD(3) - t318 * t350;
t1 = [0, -t361 * t320 + t356 * t321 + t355 * t324 + t357 * t325, t315 * qJD(5) - t365 * t310 + t391 * t311 - t363 * t393, -t311 * t331 - t315 * t329 + t321 * t332 - t324 * t328 + t383, t310, t383; 0, -t361 * t318 + t356 * t319 + t355 * t322 + t357 * t323, t313 * qJD(5) - t365 * t308 + t391 * t309 - t360 * t393, -t309 * t331 - t313 * t329 + t319 * t332 - t322 * t328 + t384, t308, t384; 0 ((t356 * qJD(2) + t357) * t351 + (t361 * qJD(2) - t355) * t354) * t348, t327 * qJD(5) - t365 * t316 + t391 * t317 - t359 * t393, -t317 * t331 - t327 * t329 + (t354 * t328 + t332 * t381) * t348 + t382, t316, t382;];
JaD_transl  = t1;
