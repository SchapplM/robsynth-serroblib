% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:43
% EndTime: 2019-02-26 20:04:44
% DurationCPUTime: 0.48s
% Computational Cost: add. (695->102), mult. (1179->168), div. (0->0), fcn. (1178->14), ass. (0->67)
t350 = qJ(3) + pkin(12);
t345 = sin(t350);
t346 = cos(t350);
t356 = sin(qJ(3));
t358 = cos(qJ(5));
t351 = qJ(5) + qJ(6);
t347 = sin(t351);
t348 = cos(t351);
t375 = r_i_i_C(1) * t348 - r_i_i_C(2) * t347;
t371 = pkin(5) * t358 + pkin(4) + t375;
t400 = r_i_i_C(3) + pkin(10) + pkin(9);
t407 = (t356 * pkin(3) + t371 * t345 - t400 * t346) * qJD(3);
t357 = sin(qJ(2));
t360 = cos(qJ(2));
t352 = sin(pkin(11));
t399 = cos(pkin(6));
t384 = t352 * t399;
t398 = cos(pkin(11));
t339 = -t357 * t384 + t398 * t360;
t374 = t347 * r_i_i_C(1) + t348 * r_i_i_C(2);
t349 = qJD(5) + qJD(6);
t355 = sin(qJ(5));
t402 = t355 * pkin(5);
t406 = qJD(5) * t402 + t374 * t349;
t397 = t347 * t349;
t396 = t348 * t349;
t353 = sin(pkin(6));
t395 = t352 * t353;
t394 = t353 * t357;
t393 = t353 * t360;
t372 = t399 * t398;
t337 = t352 * t360 + t357 * t372;
t333 = t337 * qJD(2);
t383 = t353 * t398;
t365 = -t337 * t346 + t345 * t383;
t379 = -t349 * t365 - t333;
t370 = t360 * t372;
t389 = qJD(2) * t357;
t332 = -qJD(2) * t370 + t352 * t389;
t367 = -t337 * t345 - t346 * t383;
t321 = t367 * qJD(3) - t332 * t346;
t336 = t352 * t357 - t370;
t381 = -t336 * t349 - t321;
t392 = (t381 * t347 - t379 * t348) * r_i_i_C(1) + (t379 * t347 + t381 * t348) * r_i_i_C(2);
t329 = t339 * t346 + t345 * t395;
t335 = t339 * qJD(2);
t378 = t329 * t349 - t335;
t338 = t398 * t357 + t360 * t384;
t334 = t338 * qJD(2);
t369 = -t339 * t345 + t346 * t395;
t323 = t369 * qJD(3) - t334 * t346;
t380 = -t338 * t349 - t323;
t391 = (t380 * t347 - t378 * t348) * r_i_i_C(1) + (t378 * t347 + t380 * t348) * r_i_i_C(2);
t331 = t399 * t345 + t346 * t394;
t385 = t353 * t389;
t368 = -t331 * t349 + t385;
t366 = -t345 * t394 + t399 * t346;
t386 = qJD(2) * t393;
t325 = t366 * qJD(3) + t346 * t386;
t373 = t349 * t393 - t325;
t390 = (t373 * t347 + t368 * t348) * r_i_i_C(1) + (-t368 * t347 + t373 * t348) * r_i_i_C(2);
t388 = qJD(5) * t358;
t359 = cos(qJ(3));
t363 = -pkin(3) * t359 - t400 * t345 - t371 * t346 - pkin(2);
t362 = t406 * t346 + t407;
t354 = -qJ(4) - pkin(8);
t1 = [0 (-t334 * t347 + t339 * t396) * r_i_i_C(1) + (-t334 * t348 - t339 * t397) * r_i_i_C(2) + t334 * t354 + t339 * qJD(4) + (-t334 * t355 + t339 * t388) * pkin(5) + t363 * t335 + t362 * t338, t400 * t323 - t406 * t369 + t371 * (-t329 * qJD(3) + t334 * t345) + (t334 * t356 + (-t339 * t359 - t356 * t395) * qJD(3)) * pkin(3), t335 (-t323 * t355 + t335 * t358 + (-t329 * t358 - t338 * t355) * qJD(5)) * pkin(5) + t391, t391; 0 (-t332 * t347 + t337 * t396) * r_i_i_C(1) + (-t332 * t348 - t337 * t397) * r_i_i_C(2) + t332 * t354 + t337 * qJD(4) + (-t332 * t355 + t337 * t388) * pkin(5) + t363 * t333 + t362 * t336, t400 * t321 - t406 * t367 + t371 * (t365 * qJD(3) + t332 * t345) + (t332 * t356 + (-t337 * t359 + t356 * t383) * qJD(3)) * pkin(3), t333 (-t321 * t355 + t333 * t358 + (-t336 * t355 + t358 * t365) * qJD(5)) * pkin(5) + t392, t392; 0 ((pkin(5) * t388 + t363 * qJD(2) + t375 * t349 + qJD(4)) * t357 + (-qJD(2) * t354 + (-qJD(5) * t346 + qJD(2)) * t402 - t407 + t374 * (-t346 * t349 + qJD(2))) * t360) * t353, t400 * t325 - t406 * t366 + t371 * (-t331 * qJD(3) - t345 * t386) + (-t356 * t386 + (-t399 * t356 - t359 * t394) * qJD(3)) * pkin(3), t385 (t358 * t385 - t325 * t355 + (-t331 * t358 + t355 * t393) * qJD(5)) * pkin(5) + t390, t390;];
JaD_transl  = t1;
