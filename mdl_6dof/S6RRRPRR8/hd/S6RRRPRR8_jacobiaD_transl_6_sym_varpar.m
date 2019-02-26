% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:59
% EndTime: 2019-02-26 22:20:00
% DurationCPUTime: 0.90s
% Computational Cost: add. (1093->146), mult. (1995->230), div. (0->0), fcn. (1994->14), ass. (0->82)
t432 = sin(qJ(1));
t486 = cos(pkin(6));
t447 = qJD(2) * t486 + qJD(1);
t431 = sin(qJ(2));
t462 = t432 * t486;
t452 = t431 * t462;
t471 = qJD(2) * t431;
t435 = cos(qJ(2));
t436 = cos(qJ(1));
t476 = t436 * t435;
t397 = -qJD(1) * t452 - t432 * t471 + t447 * t476;
t461 = t436 * t486;
t408 = t431 * t461 + t432 * t435;
t425 = qJ(3) + pkin(12);
t420 = sin(t425);
t421 = cos(t425);
t427 = sin(pkin(6));
t472 = qJD(1) * t427;
t466 = t432 * t472;
t478 = t427 * t436;
t467 = t421 * t478;
t391 = t420 * (-qJD(3) * t408 + t466) - qJD(3) * t467 + t397 * t421;
t451 = t435 * t461;
t477 = t432 * t431;
t407 = -t451 + t477;
t424 = qJD(5) + qJD(6);
t458 = t407 * t424 + t391;
t430 = sin(qJ(3));
t433 = cos(qJ(5));
t418 = t433 * pkin(5) + pkin(4);
t426 = qJ(5) + qJ(6);
t422 = sin(t426);
t423 = cos(t426);
t450 = r_i_i_C(1) * t423 - r_i_i_C(2) * t422;
t446 = t418 + t450;
t487 = r_i_i_C(3) + pkin(11) + pkin(10);
t496 = (t430 * pkin(3) + t420 * t446 - t421 * t487) * qJD(3);
t409 = t436 * t431 + t435 * t462;
t396 = qJD(1) * t409 + qJD(2) * t408;
t402 = -t408 * t421 + t420 * t478;
t495 = -t402 * t424 - t396;
t449 = r_i_i_C(1) * t422 + r_i_i_C(2) * t423;
t429 = sin(qJ(5));
t489 = t429 * pkin(5);
t494 = qJD(5) * t489 + t424 * t449;
t434 = cos(qJ(3));
t419 = t434 * pkin(3) + pkin(2);
t493 = t420 * t487 + t421 * t446 + t419;
t484 = t422 * t424;
t483 = t423 * t424;
t482 = t427 * t431;
t481 = t427 * t432;
t480 = t427 * t434;
t479 = t427 * t435;
t470 = qJD(2) * t435;
t394 = -qJD(1) * t451 - t436 * t470 + t447 * t477;
t410 = -t452 + t476;
t404 = t410 * t421 + t420 * t481;
t457 = -t404 * t424 - t394;
t395 = qJD(1) * t408 + qJD(2) * t409;
t445 = -t410 * t420 + t421 * t481;
t465 = t436 * t472;
t389 = qJD(3) * t445 - t395 * t421 + t420 * t465;
t460 = t409 * t424 + t389;
t384 = -t422 * t460 + t423 * t457;
t385 = t422 * t457 + t423 * t460;
t475 = t384 * r_i_i_C(1) - t385 * r_i_i_C(2);
t474 = (-t458 * t422 - t423 * t495) * r_i_i_C(1) + (t422 * t495 - t458 * t423) * r_i_i_C(2);
t406 = t420 * t486 + t421 * t482;
t463 = t427 * t471;
t444 = -t406 * t424 + t463;
t442 = -t420 * t482 + t421 * t486;
t464 = t427 * t470;
t399 = qJD(3) * t442 + t421 * t464;
t448 = t424 * t479 - t399;
t473 = (t422 * t448 + t423 * t444) * r_i_i_C(1) + (-t422 * t444 + t423 * t448) * r_i_i_C(2);
t469 = qJD(5) * t433;
t439 = qJD(3) * t402 - t397 * t420 + t421 * t466;
t438 = t494 * t421 + t496;
t428 = -qJ(4) - pkin(9);
t388 = qJD(3) * t404 - t395 * t420 - t421 * t465;
t1 = [-t407 * qJD(4) - t391 * t418 + t396 * t428 - t397 * t419 + (-r_i_i_C(1) * t458 + r_i_i_C(2) * t495) * t423 + (r_i_i_C(1) * t495 + r_i_i_C(2) * t458) * t422 + t487 * t439 + (-t436 * pkin(1) - pkin(8) * t481) * qJD(1) + (-t396 * t429 + (-t402 * t429 - t407 * t433) * qJD(5)) * pkin(5) + (-t430 * t466 + (t408 * t430 + t434 * t478) * qJD(3)) * pkin(3) (-t395 * t422 + t410 * t483) * r_i_i_C(1) + (-t395 * t423 - t410 * t484) * r_i_i_C(2) + t395 * t428 + t410 * qJD(4) + (-t395 * t429 + t410 * t469) * pkin(5) + t493 * t394 + t438 * t409, t487 * t389 - t494 * t445 - t446 * t388 + (t434 * t465 + t395 * t430 + (-t410 * t434 - t430 * t481) * qJD(3)) * pkin(3), -t394 (-t389 * t429 - t394 * t433 + (-t404 * t433 - t409 * t429) * qJD(5)) * pkin(5) + t475, t475; t385 * r_i_i_C(1) + t384 * r_i_i_C(2) + t409 * qJD(4) + t389 * t418 + t394 * t428 - t395 * t419 + t487 * t388 + (-t432 * pkin(1) + pkin(8) * t478) * qJD(1) + (-t394 * t429 + (-t404 * t429 + t409 * t433) * qJD(5)) * pkin(5) + (t430 * t465 + (-t410 * t430 + t432 * t480) * qJD(3)) * pkin(3) (t397 * t422 + t408 * t483) * r_i_i_C(1) + (t397 * t423 - t408 * t484) * r_i_i_C(2) - t397 * t428 + t408 * qJD(4) + (t397 * t429 + t408 * t469) * pkin(5) - t493 * t396 + t438 * t407, t487 * t391 - t494 * (-t408 * t420 - t467) + t446 * t439 + (t434 * t466 - t397 * t430 + (-t408 * t434 + t430 * t478) * qJD(3)) * pkin(3), t396 (-t391 * t429 + t396 * t433 + (t402 * t433 - t407 * t429) * qJD(5)) * pkin(5) + t474, t474; 0 ((pkin(5) * t469 - qJD(2) * t493 + t424 * t450 + qJD(4)) * t431 + (-qJD(2) * t428 + (-qJD(5) * t421 + qJD(2)) * t489 - t496 + t449 * (-t421 * t424 + qJD(2))) * t435) * t427, t487 * t399 - t494 * t442 + t446 * (-qJD(3) * t406 - t420 * t464) + (-t430 * t464 + (-t430 * t486 - t431 * t480) * qJD(3)) * pkin(3), t463 (t433 * t463 - t399 * t429 + (-t406 * t433 + t429 * t479) * qJD(5)) * pkin(5) + t473, t473;];
JaD_transl  = t1;
