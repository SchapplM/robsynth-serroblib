% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR1
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
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:40
% EndTime: 2019-02-26 20:10:41
% DurationCPUTime: 0.50s
% Computational Cost: add. (850->102), mult. (1129->172), div. (0->0), fcn. (1105->14), ass. (0->71)
t419 = pkin(10) + r_i_i_C(3);
t373 = sin(qJ(6));
t376 = cos(qJ(6));
t390 = t376 * r_i_i_C(1) - t373 * r_i_i_C(2);
t425 = pkin(5) + t390;
t401 = qJD(6) * t376;
t402 = qJD(6) * t373;
t424 = -r_i_i_C(1) * t402 - t401 * r_i_i_C(2);
t369 = qJD(3) + qJD(4);
t374 = sin(qJ(3));
t414 = pkin(3) * qJD(3);
t370 = qJ(3) + qJ(4);
t366 = sin(t370);
t418 = pkin(4) * t366;
t355 = -t369 * t418 - t374 * t414;
t365 = pkin(12) + t370;
t363 = sin(t365);
t364 = cos(t365);
t423 = -(t363 * t425 - t419 * t364) * t369 + t355;
t375 = sin(qJ(2));
t378 = cos(qJ(2));
t371 = sin(pkin(11));
t413 = cos(pkin(6));
t396 = t371 * t413;
t412 = cos(pkin(11));
t353 = -t375 * t396 + t412 * t378;
t421 = t373 * r_i_i_C(1) + t376 * r_i_i_C(2);
t411 = t363 * t369;
t410 = t364 * t369;
t367 = cos(t370);
t409 = t367 * t369;
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
t395 = t372 * t412;
t391 = t364 * t395;
t352 = t412 * t375 + t378 * t396;
t348 = t352 * qJD(2);
t389 = t369 * t408 - t348;
t388 = t413 * t412;
t386 = t378 * t388;
t346 = -qJD(2) * t386 + t371 * t405;
t385 = t369 * t395 + t346;
t384 = t413 * t369 + t372 * t404;
t351 = t371 * t378 + t375 * t388;
t377 = cos(qJ(3));
t383 = -t377 * pkin(3) - pkin(4) * t367 - t419 * t363 - t364 * t425 - pkin(2);
t330 = -t346 * t364 - t351 * t411 - t369 * t391;
t382 = t424 * (-t351 * t363 - t391) + t419 * t330 + t425 * (-t351 * t410 + t385 * t363);
t332 = -t353 * t411 + t389 * t364;
t381 = t424 * (-t353 * t363 + t364 * t408) + t419 * t332 + t425 * (-t353 * t410 - t389 * t363);
t337 = t384 * t364 - t369 * t399;
t380 = t424 * (t413 * t364 - t399) + t419 * t337 + t425 * (-t384 * t363 - t369 * t398);
t379 = t421 * t403 - t423;
t368 = -qJ(5) - pkin(9) - pkin(8);
t359 = -t374 * pkin(3) - t418;
t356 = -pkin(4) * t409 - t377 * t414;
t350 = t371 * t375 - t386;
t349 = t353 * qJD(2);
t347 = t351 * qJD(2);
t345 = t413 * t363 + t398;
t341 = t353 * t364 + t363 * t408;
t339 = t351 * t364 - t363 * t395;
t1 = [0 (-t348 * t373 + t353 * t401) * r_i_i_C(1) + (-t348 * t376 - t353 * t402) * r_i_i_C(2) + t348 * t368 + t353 * qJD(5) + t383 * t349 + t379 * t352, -t348 * t359 + t353 * t356 + t355 * t408 + t381 (-t353 * t409 - t389 * t366) * pkin(4) + t381, t349 (-t332 * t373 + t349 * t376) * r_i_i_C(1) + (-t332 * t376 - t349 * t373) * r_i_i_C(2) + ((-t341 * t376 - t352 * t373) * r_i_i_C(1) + (t341 * t373 - t352 * t376) * r_i_i_C(2)) * qJD(6); 0 (-t346 * t373 + t351 * t401) * r_i_i_C(1) + (-t346 * t376 - t351 * t402) * r_i_i_C(2) + t346 * t368 + t351 * qJD(5) + t383 * t347 + t379 * t350, -t346 * t359 + t351 * t356 - t355 * t395 + t382 (-t351 * t409 + t385 * t366) * pkin(4) + t382, t347 (-t330 * t373 + t347 * t376) * r_i_i_C(1) + (-t330 * t376 - t347 * t373) * r_i_i_C(2) + ((-t339 * t376 - t350 * t373) * r_i_i_C(1) + (t339 * t373 - t350 * t376) * r_i_i_C(2)) * qJD(6); 0 ((t383 * qJD(2) + t390 * qJD(6) + qJD(5)) * t375 + (-qJD(2) * t368 + t421 * (qJD(2) - t403) + t423) * t378) * t372, t413 * t355 + (t356 * t375 + t359 * t404) * t372 + t380 (-t384 * t366 - t407 * t409) * pkin(4) + t380, t397 (-t337 * t373 + t376 * t397) * r_i_i_C(1) + (-t337 * t376 - t373 * t397) * r_i_i_C(2) + ((-t345 * t376 + t373 * t406) * r_i_i_C(1) + (t345 * t373 + t376 * t406) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
