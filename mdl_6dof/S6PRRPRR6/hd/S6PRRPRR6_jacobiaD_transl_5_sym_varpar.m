% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:07
% EndTime: 2019-02-26 20:07:08
% DurationCPUTime: 0.55s
% Computational Cost: add. (500->113), mult. (1426->206), div. (0->0), fcn. (1512->14), ass. (0->80)
t461 = r_i_i_C(3) + pkin(10) + qJ(4);
t410 = sin(pkin(12));
t413 = cos(pkin(12));
t420 = cos(qJ(2));
t415 = cos(pkin(6));
t418 = sin(qJ(2));
t451 = t415 * t418;
t399 = t410 * t420 + t413 * t451;
t419 = cos(qJ(3));
t460 = t399 * t419;
t408 = pkin(13) + qJ(5);
t406 = sin(t408);
t411 = sin(pkin(7));
t459 = t406 * t411;
t407 = cos(t408);
t458 = t407 * t411;
t412 = sin(pkin(6));
t457 = t411 * t412;
t417 = sin(qJ(3));
t456 = t411 * t417;
t455 = t412 * t413;
t414 = cos(pkin(7));
t454 = t412 * t414;
t453 = t414 * t417;
t452 = t414 * t419;
t450 = t415 * t420;
t449 = t417 * t418;
t448 = t417 * t420;
t447 = t418 * t419;
t446 = t419 * t420;
t445 = qJD(2) * t418;
t444 = qJD(5) * t406;
t443 = qJD(5) * t407;
t442 = t413 * t450;
t441 = t415 * t411 * t419;
t440 = qJD(3) * t456;
t439 = t445 * t457;
t438 = r_i_i_C(1) * t407 - r_i_i_C(2) * t406;
t437 = -r_i_i_C(1) * t406 - r_i_i_C(2) * t407;
t405 = cos(pkin(13)) * pkin(4) + pkin(3);
t436 = -t405 - t438;
t398 = -t410 * t418 + t442;
t435 = -t398 * t417 - t399 * t452;
t384 = t398 * t419 - t399 * t453;
t434 = t398 * t414 - t411 * t455;
t430 = t410 * t451 - t413 * t420;
t431 = t410 * t450 + t413 * t418;
t433 = t417 * t431 + t430 * t452;
t385 = -t419 * t431 + t430 * t453;
t432 = t410 * t457 - t414 * t431;
t429 = t414 * t446 - t449;
t428 = t414 * t447 + t448;
t427 = t414 * t448 + t447;
t426 = t414 * t449 - t446;
t425 = qJD(5) * t438;
t424 = qJD(5) * t437;
t423 = sin(pkin(13)) * pkin(4) + pkin(9) - t437;
t422 = -t399 * t417 + t419 * t434;
t421 = t417 * t430 + t419 * t432;
t381 = t417 * t432 - t419 * t430;
t397 = t414 * t415 - t420 * t457;
t396 = t430 * qJD(2);
t395 = t431 * qJD(2);
t394 = t399 * qJD(2);
t393 = -qJD(2) * t442 + t410 * t445;
t392 = t426 * t412;
t389 = t410 * t454 + t411 * t431;
t388 = -t398 * t411 - t413 * t454;
t387 = t412 * t427 + t415 * t456;
t383 = (-qJD(2) * t427 - qJD(3) * t428) * t412;
t379 = t417 * t434 + t460;
t377 = qJD(3) * t441 + (-qJD(2) * t426 + qJD(3) * t429) * t412;
t376 = t415 * t440 + (qJD(2) * t428 + qJD(3) * t427) * t412;
t375 = qJD(3) * t433 + t395 * t453 + t396 * t419;
t373 = qJD(3) * t435 + t393 * t453 - t394 * t419;
t371 = qJD(3) * t421 - t395 * t419 + t396 * t453;
t370 = qJD(3) * t381 - t395 * t417 - t396 * t452;
t369 = qJD(3) * t422 - t393 * t419 - t394 * t453;
t368 = -t393 * t417 + t394 * t452 - t440 * t455 + (t398 * t453 + t460) * qJD(3);
t1 = [0 (t375 * t407 - t385 * t444) * r_i_i_C(1) + (-t375 * t406 - t385 * t443) * r_i_i_C(2) + t375 * t405 - t433 * qJD(4) + t396 * pkin(2) + t461 * (qJD(3) * t385 - t395 * t452 + t396 * t417) + (-t395 * t423 - t425 * t430) * t411, qJD(4) * t381 + t370 * t436 + t371 * t461 + t421 * t424, t370 (-t371 * t406 - t396 * t458) * r_i_i_C(1) + (-t371 * t407 + t396 * t459) * r_i_i_C(2) + ((-t381 * t407 - t389 * t406) * r_i_i_C(1) + (t381 * t406 - t389 * t407) * r_i_i_C(2)) * qJD(5), 0; 0 (t373 * t407 - t384 * t444) * r_i_i_C(1) + (-t373 * t406 - t384 * t443) * r_i_i_C(2) + t373 * t405 - t435 * qJD(4) - t394 * pkin(2) + t461 * (qJD(3) * t384 - t393 * t452 - t394 * t417) + (-t393 * t423 + t399 * t425) * t411, qJD(4) * t379 + t368 * t436 + t369 * t461 + t422 * t424, t368 (-t369 * t406 + t394 * t458) * r_i_i_C(1) + (-t369 * t407 - t394 * t459) * r_i_i_C(2) + ((-t379 * t407 - t388 * t406) * r_i_i_C(1) + (t379 * t406 - t388 * t407) * r_i_i_C(2)) * qJD(5), 0; 0 (t383 * t407 + t392 * t444) * r_i_i_C(1) + (-t383 * t406 + t392 * t443) * r_i_i_C(2) + t383 * t405 + (-t461 * (-qJD(2) * t429 + qJD(3) * t426) + t428 * qJD(4) - pkin(2) * t445 + (qJD(2) * t420 * t423 + t418 * t425) * t411) * t412, qJD(4) * t387 + t461 * t377 + (t412 * t429 + t441) * t424 + t436 * t376, t376 (-t377 * t406 + t407 * t439) * r_i_i_C(1) + (-t377 * t407 - t406 * t439) * r_i_i_C(2) + ((-t387 * t407 - t397 * t406) * r_i_i_C(1) + (t387 * t406 - t397 * t407) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
