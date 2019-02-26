% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR1
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

function JaD_transl = S6PRRPRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:07
% EndTime: 2019-02-26 20:04:08
% DurationCPUTime: 0.46s
% Computational Cost: add. (813->94), mult. (1084->159), div. (0->0), fcn. (1067->14), ass. (0->64)
t416 = pkin(10) + r_i_i_C(3);
t373 = sin(qJ(6));
t376 = cos(qJ(6));
t389 = t376 * r_i_i_C(1) - t373 * r_i_i_C(2);
t422 = pkin(5) + t389;
t401 = qJD(6) * t376;
t402 = qJD(6) * t373;
t421 = -r_i_i_C(1) * t402 - t401 * r_i_i_C(2);
t370 = qJ(3) + pkin(12);
t359 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t370);
t356 = t359 * qJD(3);
t367 = qJ(5) + t370;
t363 = sin(t367);
t364 = cos(t367);
t369 = qJD(3) + qJD(5);
t420 = -(t363 * t422 - t416 * t364) * t369 + t356;
t375 = sin(qJ(2));
t378 = cos(qJ(2));
t371 = sin(pkin(11));
t412 = cos(pkin(6));
t396 = t371 * t412;
t411 = cos(pkin(11));
t353 = -t375 * t396 + t411 * t378;
t418 = t373 * r_i_i_C(1) + t376 * r_i_i_C(2);
t410 = t363 * t369;
t409 = t364 * t369;
t372 = sin(pkin(6));
t408 = t371 * t372;
t407 = t372 * t375;
t406 = t372 * t378;
t405 = qJD(2) * t375;
t404 = qJD(2) * t378;
t403 = qJD(6) * t364;
t399 = t363 * t407;
t398 = t364 * t407;
t397 = t372 * t405;
t395 = t372 * t411;
t391 = t364 * t395;
t390 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t370);
t352 = t411 * t375 + t378 * t396;
t348 = t352 * qJD(2);
t388 = t369 * t408 - t348;
t387 = t412 * t411;
t385 = t378 * t387;
t384 = t412 * t369 + t372 * t404;
t351 = t371 * t378 + t375 * t387;
t383 = -t416 * t363 - t364 * t422 - pkin(2) + t390;
t332 = -t353 * t410 + t388 * t364;
t382 = t421 * (-t353 * t363 + t364 * t408) + t416 * t332 + t422 * (-t353 * t409 - t388 * t363);
t346 = -qJD(2) * t385 + t371 * t405;
t330 = -t346 * t364 - t351 * t410 - t369 * t391;
t381 = t421 * (-t351 * t363 - t391) + t416 * t330 + t422 * (-t351 * t409 + (t369 * t395 + t346) * t363);
t337 = t384 * t364 - t369 * t399;
t380 = t421 * (t412 * t364 - t399) + t416 * t337 + t422 * (-t384 * t363 - t369 * t398);
t379 = t418 * t403 - t420;
t368 = -pkin(9) - qJ(4) - pkin(8);
t357 = t390 * qJD(3);
t350 = t371 * t375 - t385;
t349 = t353 * qJD(2);
t347 = t351 * qJD(2);
t345 = t412 * t363 + t398;
t341 = t353 * t364 + t363 * t408;
t339 = t351 * t364 - t363 * t395;
t1 = [0 (-t348 * t373 + t353 * t401) * r_i_i_C(1) + (-t348 * t376 - t353 * t402) * r_i_i_C(2) + t348 * t368 + t353 * qJD(4) + t383 * t349 + t379 * t352, -t348 * t359 + t353 * t357 + t356 * t408 + t382, t349, t382 (-t332 * t373 + t349 * t376) * r_i_i_C(1) + (-t332 * t376 - t349 * t373) * r_i_i_C(2) + ((-t341 * t376 - t352 * t373) * r_i_i_C(1) + (t341 * t373 - t352 * t376) * r_i_i_C(2)) * qJD(6); 0 (-t346 * t373 + t351 * t401) * r_i_i_C(1) + (-t346 * t376 - t351 * t402) * r_i_i_C(2) + t346 * t368 + t351 * qJD(4) + t383 * t347 + t379 * t350, -t346 * t359 + t351 * t357 - t356 * t395 + t381, t347, t381 (-t330 * t373 + t347 * t376) * r_i_i_C(1) + (-t330 * t376 - t347 * t373) * r_i_i_C(2) + ((-t339 * t376 - t350 * t373) * r_i_i_C(1) + (t339 * t373 - t350 * t376) * r_i_i_C(2)) * qJD(6); 0 ((t383 * qJD(2) + t389 * qJD(6) + qJD(4)) * t375 + (-qJD(2) * t368 + t418 * (qJD(2) - t403) + t420) * t378) * t372, t412 * t356 + (t357 * t375 + t359 * t404) * t372 + t380, t397, t380 (-t337 * t373 + t376 * t397) * r_i_i_C(1) + (-t337 * t376 - t373 * t397) * r_i_i_C(2) + ((-t345 * t376 + t373 * t406) * r_i_i_C(1) + (t345 * t373 + t376 * t406) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
