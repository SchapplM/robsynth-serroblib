% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR7_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:32
% EndTime: 2019-02-26 19:57:32
% DurationCPUTime: 0.73s
% Computational Cost: add. (271->109), mult. (972->215), div. (0->0), fcn. (1026->14), ass. (0->74)
t373 = cos(pkin(8));
t412 = pkin(10) + r_i_i_C(3);
t413 = t412 * t373 + qJ(3);
t366 = sin(pkin(14));
t374 = cos(pkin(7));
t411 = t366 * t374;
t368 = sin(pkin(8));
t376 = sin(qJ(4));
t410 = t368 * t376;
t378 = cos(qJ(4));
t409 = t368 * t378;
t369 = sin(pkin(7));
t370 = sin(pkin(6));
t408 = t369 * t370;
t375 = cos(pkin(6));
t407 = t369 * t375;
t406 = t370 * t374;
t371 = cos(pkin(14));
t405 = t371 * t374;
t377 = sin(qJ(2));
t404 = t371 * t377;
t403 = t373 * t376;
t402 = t373 * t378;
t401 = t374 * t377;
t379 = cos(qJ(2));
t400 = t374 * t379;
t399 = t375 * t377;
t398 = t375 * t379;
t396 = qJD(2) * t377;
t395 = qJD(2) * t379;
t394 = qJD(4) * t376;
t393 = qJD(4) * t378;
t392 = t369 * t410;
t391 = t369 * t409;
t372 = cos(pkin(13));
t390 = t372 * t398;
t389 = t396 * t408;
t387 = t368 * t389;
t386 = t378 * r_i_i_C(1) - t376 * r_i_i_C(2) + pkin(3);
t367 = sin(pkin(13));
t360 = -t367 * t377 + t390;
t385 = t360 * t374 - t372 * t408;
t382 = t367 * t398 + t372 * t377;
t384 = t367 * t408 - t374 * t382;
t383 = -t366 * t377 + t371 * t400;
t361 = t367 * t379 + t372 * t399;
t381 = t367 * t399 - t372 * t379;
t353 = (-t366 * t379 - t371 * t401) * t370;
t354 = (-t366 * t401 + t371 * t379) * t370;
t380 = (t376 * r_i_i_C(1) + t378 * r_i_i_C(2)) * t373 - t412 * t368;
t359 = t375 * t374 - t379 * t408;
t358 = t381 * qJD(2);
t357 = t382 * qJD(2);
t356 = t361 * qJD(2);
t355 = -qJD(2) * t390 + t367 * t396;
t352 = qJD(2) * t354;
t351 = qJD(2) * t353;
t348 = t367 * t406 + t369 * t382;
t347 = -t360 * t369 - t372 * t406;
t346 = t370 * t404 + (t370 * t400 + t407) * t366;
t345 = t383 * t370 + t371 * t407;
t344 = -t371 * t382 + t381 * t411;
t343 = t366 * t382 + t381 * t405;
t342 = t360 * t371 - t361 * t411;
t341 = -t360 * t366 - t361 * t405;
t338 = -t357 * t371 + t358 * t411;
t337 = t357 * t366 + t358 * t405;
t334 = -t355 * t371 - t356 * t411;
t333 = t355 * t366 - t356 * t405;
t332 = t384 * t366 - t371 * t381;
t331 = t366 * t381 + t384 * t371;
t330 = t361 * t371 + t385 * t366;
t329 = -t361 * t366 + t385 * t371;
t1 = [0, t358 * pkin(2) + t386 * (t357 * t411 + t358 * t371) + t380 * (t357 * t405 - t358 * t366) + ((t343 * t402 - t344 * t376) * r_i_i_C(1) + (-t343 * t403 - t344 * t378) * r_i_i_C(2)) * qJD(4) + (-t381 * qJD(3) - t413 * t357 + ((-t357 * t376 - t381 * t393) * r_i_i_C(1) + (-t357 * t378 + t381 * t394) * r_i_i_C(2)) * t368) * t369, -t358 * t369 (t337 * t402 - t338 * t376 - t358 * t391) * r_i_i_C(1) + (-t337 * t403 - t338 * t378 + t358 * t392) * r_i_i_C(2) + ((-t331 * t403 - t332 * t378 - t348 * t410) * r_i_i_C(1) + (-t331 * t402 + t332 * t376 - t348 * t409) * r_i_i_C(2)) * qJD(4), 0, 0; 0, -t356 * pkin(2) + t386 * (t355 * t411 - t356 * t371) + t380 * (t355 * t405 + t356 * t366) + ((t341 * t402 - t342 * t376) * r_i_i_C(1) + (-t341 * t403 - t342 * t378) * r_i_i_C(2)) * qJD(4) + (t361 * qJD(3) - t413 * t355 + ((-t355 * t376 + t361 * t393) * r_i_i_C(1) + (-t355 * t378 - t361 * t394) * r_i_i_C(2)) * t368) * t369, t356 * t369 (t333 * t402 - t334 * t376 + t356 * t391) * r_i_i_C(1) + (-t333 * t403 - t334 * t378 - t356 * t392) * r_i_i_C(2) + ((-t329 * t403 - t330 * t378 - t347 * t410) * r_i_i_C(1) + (-t329 * t402 + t330 * t376 - t347 * t409) * r_i_i_C(2)) * qJD(4), 0, 0; 0 ((t353 * t402 - t354 * t376) * r_i_i_C(1) + (-t353 * t403 - t354 * t378) * r_i_i_C(2)) * qJD(4) + (-pkin(2) * t396 + (t377 * qJD(3) + t413 * t395 + ((t376 * t395 + t377 * t393) * r_i_i_C(1) + (-t377 * t394 + t378 * t395) * r_i_i_C(2)) * t368) * t369 + (t386 * (-t366 * t400 - t404) - t380 * t383) * qJD(2)) * t370, t389 (t351 * t402 - t352 * t376 + t378 * t387) * r_i_i_C(1) + (-t351 * t403 - t352 * t378 - t376 * t387) * r_i_i_C(2) + ((-t345 * t403 - t346 * t378 - t359 * t410) * r_i_i_C(1) + (-t345 * t402 + t346 * t376 - t359 * t409) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
