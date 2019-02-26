% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR2
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
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:17
% EndTime: 2019-02-26 20:19:17
% DurationCPUTime: 0.54s
% Computational Cost: add. (966->101), mult. (1433->167), div. (0->0), fcn. (1431->14), ass. (0->74)
t400 = qJ(3) + qJ(4);
t396 = cos(t400);
t398 = qJD(3) + qJD(4);
t404 = sin(qJ(3));
t397 = qJD(5) + qJD(6);
t399 = qJ(5) + qJ(6);
t393 = sin(t399);
t395 = cos(t399);
t425 = r_i_i_C(1) * t393 + r_i_i_C(2) * t395;
t403 = sin(qJ(5));
t456 = t403 * pkin(5);
t416 = -qJD(5) * t456 - t425 * t397;
t394 = sin(t400);
t451 = t394 * t398;
t455 = r_i_i_C(3) + pkin(11) + pkin(10);
t406 = cos(qJ(5));
t426 = t395 * r_i_i_C(1) - t393 * r_i_i_C(2);
t463 = t406 * pkin(5) + pkin(4) + t426;
t411 = qJD(3) * t404 * pkin(3) - (t455 * t398 + t416) * t396 + t451 * t463;
t405 = sin(qJ(2));
t408 = cos(qJ(2));
t401 = sin(pkin(12));
t454 = cos(pkin(6));
t434 = t401 * t454;
t453 = cos(pkin(12));
t387 = -t405 * t434 + t453 * t408;
t452 = t393 * t397;
t450 = t395 * t397;
t449 = t396 * t398;
t402 = sin(pkin(6));
t448 = t401 * t402;
t447 = t402 * t405;
t446 = t402 * t408;
t422 = t454 * t453;
t385 = t401 * t408 + t405 * t422;
t433 = t402 * t453;
t375 = t385 * t396 - t394 * t433;
t381 = t385 * qJD(2);
t429 = t375 * t397 - t381;
t420 = t408 * t422;
t442 = qJD(2) * t405;
t380 = -qJD(2) * t420 + t401 * t442;
t419 = t398 * t433 + t380;
t365 = -t385 * t451 - t419 * t396;
t384 = t401 * t405 - t420;
t431 = -t384 * t397 - t365;
t445 = (t431 * t393 - t429 * t395) * r_i_i_C(1) + (t429 * t393 + t431 * t395) * r_i_i_C(2);
t377 = t387 * t396 + t394 * t448;
t383 = t387 * qJD(2);
t428 = t377 * t397 - t383;
t386 = t453 * t405 + t408 * t434;
t382 = t386 * qJD(2);
t423 = t398 * t448 - t382;
t367 = -t387 * t451 + t423 * t396;
t430 = -t386 * t397 - t367;
t444 = (t430 * t393 - t428 * t395) * r_i_i_C(1) + (t428 * t393 + t430 * t395) * r_i_i_C(2);
t437 = t396 * t447;
t379 = t454 * t394 + t437;
t436 = t402 * t442;
t418 = -t379 * t397 + t436;
t435 = qJD(2) * t446;
t417 = t454 * t398 + t435;
t438 = t394 * t447;
t373 = t417 * t396 - t398 * t438;
t424 = t397 * t446 - t373;
t443 = (t424 * t393 + t418 * t395) * r_i_i_C(1) + (-t418 * t393 + t424 * t395) * r_i_i_C(2);
t441 = qJD(5) * t406;
t407 = cos(qJ(3));
t415 = -t407 * pkin(3) - t455 * t394 - t396 * t463 - pkin(2);
t414 = t416 * (-t385 * t394 - t396 * t433) + t455 * t365 + t463 * (-t385 * t449 + t419 * t394);
t413 = t416 * (-t387 * t394 + t396 * t448) + t455 * t367 + t463 * (-t387 * t449 - t423 * t394);
t412 = t416 * (t454 * t396 - t438) + t455 * t373 + t463 * (-t417 * t394 - t398 * t437);
t410 = -pkin(9) - pkin(8);
t1 = [0 (-t382 * t393 + t387 * t450) * r_i_i_C(1) + (-t382 * t395 - t387 * t452) * r_i_i_C(2) + t382 * t410 + (-t382 * t403 + t387 * t441) * pkin(5) + t415 * t383 + t411 * t386 (t382 * t404 + (-t387 * t407 - t404 * t448) * qJD(3)) * pkin(3) + t413, t413 (-t367 * t403 + t383 * t406 + (-t377 * t406 - t386 * t403) * qJD(5)) * pkin(5) + t444, t444; 0 (-t380 * t393 + t385 * t450) * r_i_i_C(1) + (-t380 * t395 - t385 * t452) * r_i_i_C(2) + t380 * t410 + (-t380 * t403 + t385 * t441) * pkin(5) + t415 * t381 + t411 * t384 (t380 * t404 + (-t385 * t407 + t404 * t433) * qJD(3)) * pkin(3) + t414, t414 (-t365 * t403 + t381 * t406 + (-t375 * t406 - t384 * t403) * qJD(5)) * pkin(5) + t445, t445; 0 ((pkin(5) * t441 + t415 * qJD(2) + t426 * t397) * t405 + ((-t410 + t425 + t456) * qJD(2) - t411) * t408) * t402 (-t404 * t435 + (-t454 * t404 - t407 * t447) * qJD(3)) * pkin(3) + t412, t412 (t406 * t436 - t373 * t403 + (-t379 * t406 + t403 * t446) * qJD(5)) * pkin(5) + t443, t443;];
JaD_transl  = t1;
