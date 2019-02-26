% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:47
% EndTime: 2019-02-26 21:55:47
% DurationCPUTime: 0.50s
% Computational Cost: add. (535->94), mult. (1256->158), div. (0->0), fcn. (1310->12), ass. (0->71)
t350 = qJ(4) + qJ(5);
t347 = sin(t350);
t354 = cos(pkin(6));
t351 = sin(pkin(12));
t353 = cos(pkin(12));
t359 = cos(qJ(2));
t382 = qJD(2) * t359;
t356 = sin(qJ(2));
t383 = qJD(2) * t356;
t399 = t351 * t383 - t353 * t382;
t327 = t399 * t354;
t368 = t359 * t351 + t356 * t353;
t333 = t368 * t354;
t338 = t356 * t351 - t359 * t353;
t360 = cos(qJ(1));
t357 = sin(qJ(1));
t384 = qJD(1) * t357;
t335 = -t351 * t382 - t353 * t383;
t390 = t357 * t335;
t317 = t390 - t333 * t384 + (-qJD(1) * t338 - t327) * t360;
t349 = qJD(4) + qJD(5);
t352 = sin(pkin(6));
t393 = t352 * t360;
t400 = -t349 * t393 + t317;
t401 = t400 * t347;
t398 = pkin(2) * t354;
t397 = r_i_i_C(3) + pkin(10) + pkin(9);
t396 = t347 * t349;
t348 = cos(t350);
t395 = t348 * t349;
t394 = t352 * t357;
t392 = t356 * t357;
t391 = t356 * t360;
t389 = t357 * t359;
t388 = t359 * t360;
t369 = t357 * t333 + t360 * t338;
t378 = qJD(1) * t393;
t365 = t349 * t369 + t378;
t321 = t360 * t333 - t357 * t338;
t314 = -t321 * qJD(1) + t357 * t327 + t360 * t335;
t372 = t349 * t394 + t314;
t310 = -t372 * t347 + t365 * t348;
t311 = t365 * t347 + t372 * t348;
t387 = t310 * r_i_i_C(1) - t311 * r_i_i_C(2);
t379 = t352 * t384;
t366 = -t321 * t349 + t379;
t375 = t400 * t348;
t386 = (t366 * t348 - t401) * r_i_i_C(1) + (-t366 * t347 - t375) * r_i_i_C(2);
t331 = t368 * t352;
t325 = t399 * t352;
t374 = -t349 * t354 + t325;
t385 = (-t331 * t395 + t374 * t347) * r_i_i_C(1) + (t331 * t396 + t374 * t348) * r_i_i_C(2);
t381 = pkin(2) * t383;
t373 = -r_i_i_C(1) * t347 - r_i_i_C(2) * t348;
t332 = t338 * t354;
t371 = -t360 * t332 - t357 * t368;
t370 = t357 * t332 - t360 * t368;
t358 = cos(qJ(4));
t345 = t358 * pkin(4) + pkin(3);
t367 = t348 * r_i_i_C(1) - t347 * r_i_i_C(2) + t345;
t364 = qJD(2) * t368;
t363 = t338 * qJD(2);
t355 = sin(qJ(4));
t362 = -qJD(4) * t355 * pkin(4) + t373 * t349;
t346 = t359 * pkin(2) + pkin(1);
t336 = -t352 * qJD(3) + t382 * t398;
t334 = t356 * t398 + (-pkin(8) - qJ(3)) * t352;
t328 = t354 * t364;
t316 = t370 * qJD(1) - t360 * t328 + t357 * t363;
t313 = t371 * qJD(1) - t357 * t328 - t360 * t363;
t1 = [(t321 * t396 - t375) * r_i_i_C(1) + (t321 * t395 + t401) * r_i_i_C(2) - t317 * t345 + t357 * t381 - t360 * t336 + t397 * t316 + (-t360 * t346 + (t373 * t352 + t334) * t357) * qJD(1) + (-t355 * t379 + (t321 * t355 + t358 * t393) * qJD(4)) * pkin(4), t397 * t314 + t362 * t370 - t367 * t313 + ((t354 * t392 - t388) * qJD(2) + (-t354 * t388 + t392) * qJD(1)) * pkin(2), t378 (t358 * t378 - t314 * t355 + (-t355 * t394 + t358 * t369) * qJD(4)) * pkin(4) + t387, t387, 0; -t360 * t381 + t311 * r_i_i_C(1) + t310 * r_i_i_C(2) + t314 * t345 - t357 * t336 + t397 * t313 + (-t334 * t360 - t346 * t357) * qJD(1) + (t355 * t378 + (t355 * t369 + t358 * t394) * qJD(4)) * pkin(4), -t397 * (t369 * qJD(1) + t360 * t327 - t390) + t362 * t371 + t367 * t316 + ((-t354 * t391 - t389) * qJD(2) + (-t354 * t389 - t391) * qJD(1)) * pkin(2), t379 (t358 * t379 - t317 * t355 + (-t321 * t358 + t355 * t393) * qJD(4)) * pkin(4) + t386, t386, 0; 0, -t397 * t325 + (-t338 * t362 - t364 * t367 - t381) * t352, 0 (t325 * t355 + (-t331 * t358 - t354 * t355) * qJD(4)) * pkin(4) + t385, t385, 0;];
JaD_transl  = t1;
