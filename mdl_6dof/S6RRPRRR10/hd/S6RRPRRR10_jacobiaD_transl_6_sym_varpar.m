% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR10
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
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:19
% EndTime: 2019-02-26 21:59:20
% DurationCPUTime: 0.74s
% Computational Cost: add. (1063->129), mult. (1896->202), div. (0->0), fcn. (1907->14), ass. (0->80)
t431 = sin(qJ(1));
t484 = cos(pkin(6));
t445 = qJD(2) * t484 + qJD(1);
t430 = sin(qJ(2));
t460 = t431 * t484;
t450 = t430 * t460;
t470 = qJD(2) * t430;
t433 = cos(qJ(2));
t434 = cos(qJ(1));
t475 = t434 * t433;
t396 = -qJD(1) * t450 - t431 * t470 + t445 * t475;
t459 = t434 * t484;
t407 = t430 * t459 + t431 * t433;
t423 = pkin(12) + qJ(4);
t419 = sin(t423);
t420 = cos(t423);
t427 = sin(pkin(6));
t471 = qJD(1) * t427;
t464 = t431 * t471;
t477 = t427 * t434;
t466 = t420 * t477;
t390 = (-qJD(4) * t407 + t464) * t419 - qJD(4) * t466 + t396 * t420;
t449 = t433 * t459;
t476 = t431 * t430;
t406 = -t449 + t476;
t424 = qJD(5) + qJD(6);
t456 = t406 * t424 + t390;
t432 = cos(qJ(5));
t418 = t432 * pkin(5) + pkin(4);
t425 = qJ(5) + qJ(6);
t421 = sin(t425);
t422 = cos(t425);
t448 = t422 * r_i_i_C(1) - t421 * r_i_i_C(2);
t444 = t418 + t448;
t485 = r_i_i_C(3) + pkin(11) + pkin(10);
t493 = (t444 * t419 - t485 * t420) * qJD(4);
t408 = t434 * t430 + t433 * t460;
t395 = t408 * qJD(1) + t407 * qJD(2);
t401 = -t407 * t420 + t419 * t477;
t492 = -t401 * t424 - t395;
t447 = t421 * r_i_i_C(1) + t422 * r_i_i_C(2);
t429 = sin(qJ(5));
t486 = t429 * pkin(5);
t491 = qJD(5) * t486 + t447 * t424;
t417 = cos(pkin(12)) * pkin(3) + pkin(2);
t489 = t485 * t419 + t444 * t420 + t417;
t482 = t421 * t424;
t481 = t422 * t424;
t480 = t427 * t430;
t479 = t427 * t431;
t478 = t427 * t433;
t469 = qJD(2) * t433;
t393 = -qJD(1) * t449 - t434 * t469 + t445 * t476;
t409 = -t450 + t475;
t403 = t409 * t420 + t419 * t479;
t455 = -t403 * t424 - t393;
t394 = t407 * qJD(1) + t408 * qJD(2);
t443 = -t409 * t419 + t420 * t479;
t463 = t434 * t471;
t388 = t443 * qJD(4) - t394 * t420 + t419 * t463;
t458 = t408 * t424 + t388;
t383 = -t458 * t421 + t455 * t422;
t384 = t455 * t421 + t458 * t422;
t474 = t383 * r_i_i_C(1) - t384 * r_i_i_C(2);
t473 = (-t456 * t421 - t422 * t492) * r_i_i_C(1) + (t421 * t492 - t456 * t422) * r_i_i_C(2);
t405 = t484 * t419 + t420 * t480;
t461 = t427 * t470;
t442 = -t405 * t424 + t461;
t440 = -t419 * t480 + t484 * t420;
t462 = t427 * t469;
t398 = t440 * qJD(4) + t420 * t462;
t446 = t424 * t478 - t398;
t472 = (t446 * t421 + t442 * t422) * r_i_i_C(1) + (-t442 * t421 + t446 * t422) * r_i_i_C(2);
t468 = qJD(5) * t432;
t465 = pkin(3) * sin(pkin(12)) + pkin(8);
t437 = t401 * qJD(4) - t396 * t419 + t420 * t464;
t436 = t491 * t420 + t493;
t428 = -pkin(9) - qJ(3);
t387 = t403 * qJD(4) - t394 * t419 - t420 * t463;
t1 = [-t406 * qJD(3) - t390 * t418 + t395 * t428 - t396 * t417 + (-t456 * r_i_i_C(1) + r_i_i_C(2) * t492) * t422 + (r_i_i_C(1) * t492 + t456 * r_i_i_C(2)) * t421 + t485 * t437 + (-t434 * pkin(1) - t465 * t479) * qJD(1) + (-t395 * t429 + (-t401 * t429 - t406 * t432) * qJD(5)) * pkin(5) (-t394 * t421 + t409 * t481) * r_i_i_C(1) + (-t394 * t422 - t409 * t482) * r_i_i_C(2) + t394 * t428 + t409 * qJD(3) + (-t394 * t429 + t409 * t468) * pkin(5) + t489 * t393 + t436 * t408, -t393, -t444 * t387 + t485 * t388 - t443 * t491 (-t388 * t429 - t393 * t432 + (-t403 * t432 - t408 * t429) * qJD(5)) * pkin(5) + t474, t474; t384 * r_i_i_C(1) + t383 * r_i_i_C(2) + t408 * qJD(3) + t388 * t418 + t393 * t428 - t394 * t417 + t485 * t387 + (-pkin(1) * t431 + t465 * t477) * qJD(1) + (-t393 * t429 + (-t403 * t429 + t408 * t432) * qJD(5)) * pkin(5) (t396 * t421 + t407 * t481) * r_i_i_C(1) + (t396 * t422 - t407 * t482) * r_i_i_C(2) - t396 * t428 + t407 * qJD(3) + (t396 * t429 + t407 * t468) * pkin(5) - t489 * t395 + t436 * t406, t395, t485 * t390 - t491 * (-t407 * t419 - t466) + t444 * t437 (-t390 * t429 + t395 * t432 + (t401 * t432 - t406 * t429) * qJD(5)) * pkin(5) + t473, t473; 0 ((pkin(5) * t468 - qJD(2) * t489 + t448 * t424 + qJD(3)) * t430 + (-qJD(2) * t428 + (-qJD(5) * t420 + qJD(2)) * t486 - t493 + t447 * (-t420 * t424 + qJD(2))) * t433) * t427, t461, t485 * t398 - t491 * t440 + t444 * (-t405 * qJD(4) - t419 * t462) (t432 * t461 - t398 * t429 + (-t405 * t432 + t429 * t478) * qJD(5)) * pkin(5) + t472, t472;];
JaD_transl  = t1;
