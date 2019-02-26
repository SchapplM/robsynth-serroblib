% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:38
% EndTime: 2019-02-26 20:18:39
% DurationCPUTime: 0.51s
% Computational Cost: add. (1182->103), mult. (1356->173), div. (0->0), fcn. (1329->14), ass. (0->72)
t422 = pkin(11) + r_i_i_C(3);
t376 = sin(qJ(6));
t379 = cos(qJ(6));
t392 = t379 * r_i_i_C(1) - t376 * r_i_i_C(2);
t428 = pkin(5) + t392;
t404 = qJD(6) * t379;
t405 = qJD(6) * t376;
t427 = -r_i_i_C(1) * t405 - t404 * r_i_i_C(2);
t371 = qJD(3) + qJD(4);
t377 = sin(qJ(3));
t417 = pkin(3) * qJD(3);
t373 = qJ(3) + qJ(4);
t368 = sin(t373);
t421 = pkin(4) * t368;
t358 = -t371 * t421 - t377 * t417;
t370 = qJ(5) + t373;
t365 = sin(t370);
t366 = cos(t370);
t367 = qJD(5) + t371;
t426 = -(t365 * t428 - t422 * t366) * t367 + t358;
t378 = sin(qJ(2));
t381 = cos(qJ(2));
t374 = sin(pkin(12));
t416 = cos(pkin(6));
t398 = t374 * t416;
t415 = cos(pkin(12));
t357 = -t378 * t398 + t415 * t381;
t424 = t376 * r_i_i_C(1) + t379 * r_i_i_C(2);
t414 = t365 * t367;
t413 = t366 * t367;
t369 = cos(t373);
t412 = t369 * t371;
t375 = sin(pkin(6));
t411 = t374 * t375;
t410 = t375 * t378;
t409 = t375 * t381;
t408 = qJD(2) * t378;
t407 = qJD(2) * t381;
t406 = qJD(6) * t366;
t402 = t365 * t410;
t401 = t366 * t410;
t400 = t375 * t408;
t399 = t375 * t407;
t397 = t375 * t415;
t393 = t366 * t397;
t356 = t415 * t378 + t381 * t398;
t350 = t356 * qJD(2);
t391 = t367 * t411 - t350;
t390 = t416 * t415;
t388 = t381 * t390;
t387 = t416 * t367 + t399;
t355 = t374 * t381 + t378 * t390;
t380 = cos(qJ(3));
t386 = -t380 * pkin(3) - pkin(4) * t369 - t422 * t365 - t366 * t428 - pkin(2);
t348 = -qJD(2) * t388 + t374 * t408;
t332 = -t348 * t366 - t355 * t414 - t367 * t393;
t385 = t427 * (-t355 * t365 - t393) + t422 * t332 + t428 * (-t355 * t413 + (t367 * t397 + t348) * t365);
t334 = -t357 * t414 + t391 * t366;
t384 = t427 * (-t357 * t365 + t366 * t411) + t422 * t334 + t428 * (-t357 * t413 - t391 * t365);
t339 = t387 * t366 - t367 * t402;
t383 = t427 * (t416 * t366 - t402) + t422 * t339 + t428 * (-t387 * t365 - t367 * t401);
t382 = t424 * t406 - t426;
t372 = -pkin(10) - pkin(9) - pkin(8);
t361 = -t377 * pkin(3) - t421;
t359 = -pkin(4) * t412 - t380 * t417;
t354 = t374 * t378 - t388;
t351 = t357 * qJD(2);
t349 = t355 * qJD(2);
t347 = t416 * t365 + t401;
t343 = t357 * t366 + t365 * t411;
t341 = t355 * t366 - t365 * t397;
t1 = [0 (-t350 * t376 + t357 * t404) * r_i_i_C(1) + (-t350 * t379 - t357 * t405) * r_i_i_C(2) + t350 * t372 + t386 * t351 + t382 * t356, -t350 * t361 + t357 * t359 + t358 * t411 + t384 (-t357 * t412 + (-t371 * t411 + t350) * t368) * pkin(4) + t384, t384 (-t334 * t376 + t351 * t379) * r_i_i_C(1) + (-t334 * t379 - t351 * t376) * r_i_i_C(2) + ((-t343 * t379 - t356 * t376) * r_i_i_C(1) + (t343 * t376 - t356 * t379) * r_i_i_C(2)) * qJD(6); 0 (-t348 * t376 + t355 * t404) * r_i_i_C(1) + (-t348 * t379 - t355 * t405) * r_i_i_C(2) + t348 * t372 + t386 * t349 + t382 * t354, -t348 * t361 + t355 * t359 - t358 * t397 + t385 (-t355 * t412 + (t371 * t397 + t348) * t368) * pkin(4) + t385, t385 (-t332 * t376 + t349 * t379) * r_i_i_C(1) + (-t332 * t379 - t349 * t376) * r_i_i_C(2) + ((-t341 * t379 - t354 * t376) * r_i_i_C(1) + (t341 * t376 - t354 * t379) * r_i_i_C(2)) * qJD(6); 0 ((t386 * qJD(2) + t392 * qJD(6)) * t378 + (-qJD(2) * t372 + t424 * (qJD(2) - t406) + t426) * t381) * t375, t416 * t358 + (t359 * t378 + t361 * t407) * t375 + t383 (-t410 * t412 + (-t416 * t371 - t399) * t368) * pkin(4) + t383, t383 (-t339 * t376 + t379 * t400) * r_i_i_C(1) + (-t339 * t379 - t376 * t400) * r_i_i_C(2) + ((-t347 * t379 + t376 * t409) * r_i_i_C(1) + (t347 * t376 + t379 * t409) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
