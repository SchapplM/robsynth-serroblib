% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR9_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR9_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:29
% EndTime: 2019-02-26 22:20:29
% DurationCPUTime: 0.46s
% Computational Cost: add. (398->123), mult. (1234->207), div. (0->0), fcn. (1244->12), ass. (0->66)
t349 = sin(pkin(13));
t357 = cos(qJ(3));
t390 = cos(pkin(13));
t370 = qJD(3) * t390;
t354 = sin(qJ(3));
t379 = qJD(3) * t354;
t393 = t349 * t379 - t357 * t370;
t392 = pkin(3) * t354;
t391 = pkin(10) + qJ(4);
t350 = sin(pkin(7));
t351 = sin(pkin(6));
t389 = t350 * t351;
t358 = cos(qJ(2));
t388 = t354 * t358;
t355 = sin(qJ(2));
t356 = sin(qJ(1));
t387 = t355 * t356;
t386 = t355 * t357;
t359 = cos(qJ(1));
t385 = t355 * t359;
t384 = t356 * t358;
t383 = t358 * t359;
t382 = qJD(1) * t356;
t381 = qJD(1) * t359;
t380 = qJD(2) * t355;
t378 = qJD(3) * t357;
t377 = pkin(3) * t379;
t376 = pkin(3) * t378;
t353 = cos(pkin(6));
t375 = t353 * t387;
t374 = t351 * t382;
t373 = t351 * t381;
t319 = t393 * t350;
t352 = cos(pkin(7));
t371 = -r_i_i_C(1) * t319 + qJD(4) * t352 + t350 * t376;
t338 = -t357 * t349 - t354 * t390;
t324 = t338 * t350;
t369 = -r_i_i_C(1) * t324 + t350 * t392 + t391 * t352 + pkin(9);
t331 = -t353 * t383 + t387;
t367 = t353 * t384 + t385;
t332 = t353 * t385 + t384;
t366 = t375 - t383;
t365 = -t354 * t349 + t357 * t390;
t318 = -qJD(1) * t375 - t356 * t380 + (qJD(2) * t353 + qJD(1)) * t383;
t321 = t393 * t352;
t330 = -t349 * t378 - t354 * t370;
t364 = -t318 * t365 - t331 * t321 - t332 * t330;
t361 = qJD(3) * t338;
t322 = t352 * t361;
t329 = t365 * qJD(3);
t363 = -t318 * t338 + t331 * t322 + t332 * t329;
t325 = t365 * t352;
t326 = t338 * t352;
t328 = -t391 * t350 + t352 * t392;
t362 = -r_i_i_C(1) * t326 + r_i_i_C(2) * t325 - r_i_i_C(3) * t350 + t328;
t315 = t331 * qJD(1) + t366 * qJD(2);
t316 = t332 * qJD(1) + t367 * qJD(2);
t360 = -t315 * t326 - t316 * t365 + t321 * t367 - t330 * t366;
t348 = pkin(3) * t357 + pkin(2);
t336 = -qJD(4) * t350 + t352 * t376;
t323 = t365 * t350;
t320 = t350 * t361;
t317 = t367 * qJD(1) + t332 * qJD(2);
t314 = -t315 * t350 + t352 * t373;
t313 = t315 * t325 - t316 * t338 - t367 * t322 + t366 * t329 + (t320 * t356 + t323 * t381) * t351;
t1 = [t364 * r_i_i_C(1) + t363 * r_i_i_C(2) - t318 * t348 + t332 * t377 + t331 * t336 - pkin(1) * t381 + t362 * t317 + ((r_i_i_C(2) * t320 + t371) * t359 + (-r_i_i_C(2) * t323 - r_i_i_C(3) * t352 - t369) * t382) * t351 (t315 * t365 - t321 * t366 - t330 * t367) * r_i_i_C(1) + (t315 * t338 + t322 * t366 + t329 * t367) * r_i_i_C(2) + t315 * t348 + t367 * t377 + t366 * t336 + t362 * t316, t313 * r_i_i_C(1) + ((t319 * t356 + t324 * t381) * t351 - t360) * r_i_i_C(2) + (t316 * t354 + (t315 * t352 + t350 * t373) * t357 + (t366 * t357 + (t352 * t367 - t356 * t389) * t354) * qJD(3)) * pkin(3), t314, 0, 0; t360 * r_i_i_C(1) + t313 * r_i_i_C(2) + t314 * r_i_i_C(3) - t316 * t348 + t366 * t377 + t315 * t328 - t367 * t336 - pkin(1) * t382 + (t371 * t356 + t369 * t381) * t351 (-t317 * t365 + t321 * t332 - t330 * t331) * r_i_i_C(1) + (-t317 * t338 - t322 * t332 + t329 * t331) * r_i_i_C(2) - t317 * t348 + t331 * t377 - t332 * t336 - t362 * t318 (-t317 * t325 - t363) * r_i_i_C(1) + (-t317 * t326 + t364) * r_i_i_C(2) + ((-t320 * t359 + t323 * t382) * r_i_i_C(1) + (-t319 * t359 + t324 * t382) * r_i_i_C(2)) * t351 + (-t318 * t354 + (-t317 * t352 + t350 * t374) * t357 + (-t332 * t357 + (t331 * t352 + t359 * t389) * t354) * qJD(3)) * pkin(3), t317 * t350 + t352 * t374, 0, 0; 0 ((t321 * t355 + t330 * t358) * r_i_i_C(1) + (-t322 * t355 - t329 * t358) * r_i_i_C(2) - t358 * t377 - t355 * t336 + ((-r_i_i_C(1) * t365 - r_i_i_C(2) * t338 - t348) * t355 - t362 * t358) * qJD(2)) * t351 (r_i_i_C(1) * t320 + r_i_i_C(2) * t319 - t350 * t377) * t353 + ((t322 * t358 - t329 * t355) * r_i_i_C(1) + (t321 * t358 - t330 * t355) * r_i_i_C(2) + ((-t325 * t355 + t338 * t358) * r_i_i_C(1) + (-t326 * t355 - t358 * t365) * r_i_i_C(2)) * qJD(2) + ((-t352 * t388 - t386) * qJD(3) + (-t352 * t386 - t388) * qJD(2)) * pkin(3)) * t351, t380 * t389, 0, 0;];
JaD_transl  = t1;
