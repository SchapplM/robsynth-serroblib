% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR13_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR13_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:11
% EndTime: 2019-02-26 22:01:12
% DurationCPUTime: 0.48s
% Computational Cost: add. (821->112), mult. (1936->173), div. (0->0), fcn. (1942->12), ass. (0->75)
t423 = qJD(5) + qJD(6);
t431 = cos(qJ(5));
t424 = qJ(5) + qJ(6);
t421 = sin(t424);
t422 = cos(t424);
t453 = t422 * r_i_i_C(1) - t421 * r_i_i_C(2);
t479 = pkin(5) * qJD(5);
t438 = -t423 * t453 - t431 * t479;
t426 = cos(pkin(6));
t429 = sin(qJ(2));
t434 = cos(qJ(1));
t472 = t434 * t429;
t430 = sin(qJ(1));
t433 = cos(qJ(2));
t473 = t430 * t433;
t412 = t426 * t472 + t473;
t413 = t426 * t473 + t472;
t400 = qJD(1) * t413 + qJD(2) * t412;
t471 = t434 * t433;
t474 = t430 * t429;
t411 = -t426 * t471 + t474;
t428 = sin(qJ(4));
t432 = cos(qJ(4));
t425 = sin(pkin(6));
t467 = qJD(1) * t425;
t462 = t430 * t467;
t475 = t425 * t434;
t463 = t432 * t475;
t483 = t428 * (qJD(4) * t411 + t462) - qJD(4) * t463 - t400 * t432;
t448 = t411 * t432 + t428 * t475;
t393 = qJD(4) * t448 + t400 * t428 + t432 * t462;
t420 = t431 * pkin(5) + pkin(4);
t450 = t420 + t453;
t480 = r_i_i_C(3) + pkin(11) + pkin(10);
t437 = t428 * t450 - t432 * t480 + qJ(3);
t482 = -pkin(2) - pkin(9);
t481 = pkin(3) + pkin(8);
t478 = t425 * t429;
t477 = t425 * t430;
t476 = t425 * t433;
t464 = t426 * t474;
t466 = qJD(2) * t429;
t401 = -qJD(1) * t464 - t430 * t466 + (qJD(2) * t426 + qJD(1)) * t471;
t407 = -t411 * t428 + t463;
t454 = -t407 * t423 - t401;
t458 = -t412 * t423 - t393;
t470 = (t421 * t458 - t454 * t422) * r_i_i_C(1) + (t421 * t454 + t422 * t458) * r_i_i_C(2);
t399 = qJD(1) * t412 + qJD(2) * t413;
t405 = t413 * t428 + t432 * t477;
t456 = -t405 * t423 - t399;
t444 = t464 - t471;
t398 = qJD(1) * t411 + qJD(2) * t444;
t447 = t413 * t432 - t428 * t477;
t461 = t434 * t467;
t396 = qJD(4) * t447 - t398 * t428 + t432 * t461;
t457 = -t423 * t444 + t396;
t387 = -t421 * t457 + t422 * t456;
t388 = t421 * t456 + t422 * t457;
t469 = t387 * r_i_i_C(1) - t388 * r_i_i_C(2);
t445 = -t426 * t432 + t428 * t476;
t460 = qJD(2) * t476;
t443 = t423 * t445 + t460;
t446 = t426 * t428 + t432 * t476;
t459 = t425 * t466;
t402 = qJD(4) * t446 - t428 * t459;
t451 = -t423 * t478 + t402;
t468 = (t421 * t451 + t422 * t443) * r_i_i_C(1) + (-t421 * t443 + t422 * t451) * r_i_i_C(2);
t452 = -t421 * r_i_i_C(1) - t422 * r_i_i_C(2);
t449 = t452 + t482;
t427 = sin(qJ(5));
t440 = t427 * pkin(5) - t449;
t439 = t423 * t452 - t427 * t479;
t436 = qJD(3) + t439 * t428 + (t428 * t480 + t432 * t450) * qJD(4);
t395 = qJD(4) * t405 + t398 * t432 + t428 * t461;
t1 = [-t400 * qJ(3) - t411 * qJD(3) - t450 * t393 + ((-t407 * t421 - t412 * t422) * r_i_i_C(1) + (-t407 * t422 + t412 * t421) * r_i_i_C(2)) * t423 + t449 * t401 - t480 * t483 + (-t434 * pkin(1) - t477 * t481) * qJD(1) + (-t401 * t427 + (-t407 * t427 - t412 * t431) * qJD(5)) * pkin(5), t440 * t398 - t437 * t399 + t438 * t413 - t436 * t444, -t398, -t395 * t450 + t396 * t480 + t439 * t447 (-t396 * t427 - t399 * t431 + (-t405 * t431 + t427 * t444) * qJD(5)) * pkin(5) + t469, t469; t388 * r_i_i_C(1) + t387 * r_i_i_C(2) - t398 * qJ(3) + t413 * qJD(3) + t396 * t420 + t482 * t399 + t480 * t395 + (-pkin(1) * t430 + t475 * t481) * qJD(1) + (-t399 * t427 + (-t405 * t427 - t431 * t444) * qJD(5)) * pkin(5), -t400 * t440 + t401 * t437 + t411 * t438 + t412 * t436, t400, t480 * t393 + t439 * t448 - t450 * t483 (-t393 * t427 + t401 * t431 + (t407 * t431 - t412 * t427) * qJD(5)) * pkin(5) + t470, t470; 0 ((qJD(2) * t437 - t438) * t433 + (-qJD(2) * t440 + t436) * t429) * t425, t459, -t480 * t402 - t439 * t446 + t450 * (qJD(4) * t445 + t432 * t459) (t431 * t460 + t402 * t427 + (-t427 * t478 + t431 * t445) * qJD(5)) * pkin(5) + t468, t468;];
JaD_transl  = t1;
