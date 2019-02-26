% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:09
% EndTime: 2019-02-26 20:08:09
% DurationCPUTime: 0.52s
% Computational Cost: add. (471->110), mult. (1558->203), div. (0->0), fcn. (1646->12), ass. (0->72)
t403 = sin(qJ(3));
t406 = cos(qJ(3));
t396 = sin(pkin(12));
t399 = cos(pkin(12));
t407 = cos(qJ(2));
t401 = cos(pkin(6));
t404 = sin(qJ(2));
t436 = t401 * t404;
t414 = t396 * t436 - t399 * t407;
t400 = cos(pkin(7));
t435 = t401 * t407;
t415 = t396 * t435 + t399 * t404;
t397 = sin(pkin(7));
t398 = sin(pkin(6));
t444 = t397 * t398;
t416 = t396 * t444 - t400 * t415;
t448 = t403 * t414 + t406 * t416;
t388 = t396 * t407 + t399 * t436;
t426 = t399 * t435;
t387 = -t396 * t404 + t426;
t440 = t398 * t399;
t417 = -t387 * t400 + t397 * t440;
t367 = t388 * t403 + t406 * t417;
t446 = t388 * t406;
t402 = sin(qJ(5));
t443 = t397 * t402;
t442 = t397 * t403;
t405 = cos(qJ(5));
t441 = t397 * t405;
t439 = t398 * t400;
t438 = t400 * t403;
t437 = t400 * t406;
t434 = t403 * t404;
t433 = t403 * t407;
t432 = t404 * t406;
t431 = t406 * t407;
t430 = qJD(2) * t404;
t429 = qJD(5) * t402;
t428 = qJD(5) * t405;
t427 = r_i_i_C(3) + pkin(10) + pkin(3);
t425 = t400 * t431;
t424 = t401 * t397 * t406;
t423 = qJD(3) * t442;
t422 = t430 * t444;
t421 = t405 * r_i_i_C(1) - t402 * r_i_i_C(2);
t420 = -t402 * r_i_i_C(1) - t405 * r_i_i_C(2);
t419 = qJ(4) - t420;
t418 = pkin(4) + pkin(9) + t421;
t373 = t387 * t403 + t388 * t437;
t374 = -t403 * t415 - t414 * t437;
t413 = t425 - t434;
t412 = t400 * t432 + t433;
t411 = t400 * t433 + t432;
t410 = qJD(5) * t420;
t409 = qJD(5) * t421 + qJD(4);
t408 = t403 * t416 - t406 * t414;
t386 = t401 * t400 - t407 * t444;
t385 = t414 * qJD(2);
t384 = t415 * qJD(2);
t383 = t388 * qJD(2);
t382 = -qJD(2) * t426 + t396 * t430;
t381 = t412 * t398;
t378 = t396 * t439 + t397 * t415;
t377 = -t387 * t397 - t399 * t439;
t375 = -t398 * t413 - t424;
t371 = (-qJD(2) * t425 - qJD(3) * t431 + (qJD(3) * t400 + qJD(2)) * t434) * t398;
t365 = t401 * t423 + (qJD(2) * t412 + qJD(3) * t411) * t398;
t363 = -t384 * t437 + t385 * t403 + (-t406 * t415 + t414 * t438) * qJD(3);
t361 = -t382 * t437 - t383 * t403 + (t387 * t406 - t388 * t438) * qJD(3);
t359 = qJD(3) * t408 - t384 * t403 - t385 * t437;
t357 = -t382 * t403 + t383 * t437 - t423 * t440 + (t387 * t438 + t446) * qJD(3);
t1 = [0 (t363 * t402 + t374 * t428) * r_i_i_C(1) + (t363 * t405 - t374 * t429) * r_i_i_C(2) + t363 * qJ(4) + t374 * qJD(4) + t385 * pkin(2) + t427 * (-qJD(3) * t374 + t384 * t438 + t385 * t406) + (-t384 * t418 - t410 * t414) * t397, t409 * t408 + t419 * (t448 * qJD(3) - t384 * t406 + t385 * t438) - t427 * t359, t359 (t359 * t405 + t385 * t443) * r_i_i_C(1) + (-t359 * t402 + t385 * t441) * r_i_i_C(2) + ((-t378 * t405 + t402 * t448) * r_i_i_C(1) + (t378 * t402 + t405 * t448) * r_i_i_C(2)) * qJD(5), 0; 0 (t361 * t402 + t373 * t428) * r_i_i_C(1) + (t361 * t405 - t373 * t429) * r_i_i_C(2) + t361 * qJ(4) + t373 * qJD(4) - t383 * pkin(2) + t427 * (-qJD(3) * t373 + t382 * t438 - t383 * t406) + (-t382 * t418 + t388 * t410) * t397, t409 * (-t403 * t417 + t446) + t419 * (-t367 * qJD(3) - t382 * t406 - t383 * t438) - t427 * t357, t357 (t357 * t405 - t383 * t443) * r_i_i_C(1) + (-t357 * t402 - t383 * t441) * r_i_i_C(2) + ((-t367 * t402 - t377 * t405) * r_i_i_C(1) + (-t367 * t405 + t377 * t402) * r_i_i_C(2)) * qJD(5), 0; 0 (-t371 * t402 + t381 * t428) * r_i_i_C(1) + (-t371 * t405 - t381 * t429) * r_i_i_C(2) - t371 * qJ(4) + t381 * qJD(4) + (-t427 * (qJD(2) * t411 + qJD(3) * t412) - pkin(2) * t430 + (qJD(2) * t407 * t418 + t404 * t410) * t397) * t398, t409 * (t398 * t411 + t401 * t442) + t419 * (qJD(3) * t424 + (t413 * qJD(3) + (-t400 * t434 + t431) * qJD(2)) * t398) - t427 * t365, t365 (t365 * t405 - t402 * t422) * r_i_i_C(1) + (-t365 * t402 - t405 * t422) * r_i_i_C(2) + ((-t375 * t402 - t386 * t405) * r_i_i_C(1) + (-t375 * t405 + t386 * t402) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
