% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR9
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
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:48
% EndTime: 2019-02-26 21:58:49
% DurationCPUTime: 0.70s
% Computational Cost: add. (1257->133), mult. (1819->207), div. (0->0), fcn. (1811->14), ass. (0->84)
t497 = pkin(11) + r_i_i_C(3);
t442 = sin(qJ(6));
t445 = cos(qJ(6));
t457 = t445 * r_i_i_C(1) - t442 * r_i_i_C(2);
t455 = pkin(5) + t457;
t472 = qJD(6) * t445;
t473 = qJD(6) * t442;
t505 = -r_i_i_C(1) * t473 - t472 * r_i_i_C(2);
t441 = cos(pkin(6));
t443 = sin(qJ(2));
t447 = cos(qJ(1));
t479 = t447 * t443;
t444 = sin(qJ(1));
t446 = cos(qJ(2));
t480 = t444 * t446;
t421 = t441 * t479 + t480;
t438 = pkin(12) + qJ(4);
t436 = qJ(5) + t438;
t432 = sin(t436);
t433 = cos(t436);
t440 = sin(pkin(6));
t482 = t440 * t447;
t411 = -t421 * t433 + t432 * t482;
t478 = t447 * t446;
t465 = t441 * t478;
t481 = t444 * t443;
t420 = -t465 + t481;
t504 = -t411 * t442 - t420 * t445;
t503 = t411 * t445 - t420 * t442;
t434 = sin(t438);
t439 = qJD(4) + qJD(5);
t492 = pkin(4) * qJD(4);
t502 = -t434 * t492 - (t455 * t432 - t497 * t433) * t439;
t500 = t442 * r_i_i_C(1) + t445 * r_i_i_C(2);
t435 = cos(t438);
t425 = pkin(4) * t435 + cos(pkin(12)) * pkin(3) + pkin(2);
t498 = t497 * t432 + t455 * t433 + t425;
t493 = pkin(8) + pkin(4) * t434 + sin(pkin(12)) * pkin(3);
t422 = t441 * t480 + t479;
t407 = t422 * qJD(1) + t421 * qJD(2);
t491 = t407 * t442;
t490 = t407 * t445;
t467 = t441 * t481;
t423 = -t467 + t478;
t487 = t423 * t433;
t486 = t432 * t439;
t485 = t440 * t443;
t484 = t440 * t444;
t483 = t440 * t446;
t477 = qJD(1) * t440;
t476 = qJD(2) * t443;
t475 = qJD(2) * t446;
t474 = qJD(6) * t433;
t469 = t432 * t485;
t468 = t433 * t485;
t466 = t433 * t482;
t464 = t444 * t477;
t463 = t447 * t477;
t462 = t440 * t475;
t461 = t440 * t476;
t459 = qJD(2) * t441 + qJD(1);
t408 = -qJD(1) * t467 - t444 * t476 + t459 * t478;
t460 = -t408 * t433 + t439 * t466;
t406 = t421 * qJD(1) + t422 * qJD(2);
t456 = t439 * t484 - t406;
t454 = -t421 * t439 + t464;
t453 = t439 * t441 + t462;
t396 = t454 * t433 + (t439 * t482 - t408) * t432;
t397 = -t421 * t486 + t432 * t464 - t460;
t451 = t505 * (-t421 * t432 - t466) + t497 * t397 + t455 * t396;
t394 = t456 * t432 - t433 * t463 + t439 * t487;
t395 = -t423 * t486 + t432 * t463 + t456 * t433;
t450 = t505 * (-t423 * t432 + t433 * t484) + t497 * t395 - t455 * t394;
t404 = t453 * t433 - t439 * t469;
t449 = t505 * (t441 * t433 - t469) + t497 * t404 + t455 * (-t453 * t432 - t439 * t468);
t448 = t500 * t474 - t502;
t437 = -pkin(10) - pkin(9) - qJ(3);
t417 = t441 * t432 + t468;
t413 = t432 * t484 + t487;
t405 = -qJD(1) * t465 - t447 * t475 + t459 * t481;
t399 = -t454 * t432 + t460;
t387 = t395 * t445 - t405 * t442 + (-t413 * t442 + t422 * t445) * qJD(6);
t386 = -t395 * t442 - t405 * t445 + (-t413 * t445 - t422 * t442) * qJD(6);
t1 = [(t399 * t445 - t491) * r_i_i_C(1) + (-t399 * t442 - t490) * r_i_i_C(2) + t399 * pkin(5) - t408 * t425 + t407 * t437 - t420 * qJD(3) + t497 * t396 + (t504 * r_i_i_C(1) - t503 * r_i_i_C(2)) * qJD(6) + (t421 * t434 + t435 * t482) * t492 + (-t447 * pkin(1) - t493 * t484) * qJD(1) (-t406 * t442 + t423 * t472) * r_i_i_C(1) + (-t406 * t445 - t423 * t473) * r_i_i_C(2) + t406 * t437 + t423 * qJD(3) + t498 * t405 + t448 * t422, -t405 (t435 * t463 + t406 * t434 + (-t423 * t435 - t434 * t484) * qJD(4)) * pkin(4) + t450, t450, t386 * r_i_i_C(1) - t387 * r_i_i_C(2); t395 * pkin(5) + t387 * r_i_i_C(1) + t386 * r_i_i_C(2) + t422 * qJD(3) + t405 * t437 - t406 * t425 + t497 * t394 + (-t423 * t434 + t435 * t484) * t492 + (-pkin(1) * t444 + t493 * t482) * qJD(1) (t408 * t442 + t421 * t472) * r_i_i_C(1) + (t408 * t445 - t421 * t473) * r_i_i_C(2) - t408 * t437 + t421 * qJD(3) - t498 * t407 + t448 * t420, t407 (t435 * t464 - t408 * t434 + (-t421 * t435 + t434 * t482) * qJD(4)) * pkin(4) + t451, t451 (-t397 * t442 + t490) * r_i_i_C(1) + (-t397 * t445 - t491) * r_i_i_C(2) + (t503 * r_i_i_C(1) + t504 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t498 + t457 * qJD(6) + qJD(3)) * t443 + (-qJD(2) * t437 + t500 * (qJD(2) - t474) + t502) * t446) * t440, t461 (-t434 * t462 + (-t434 * t441 - t435 * t485) * qJD(4)) * pkin(4) + t449, t449 (-t404 * t442 + t445 * t461) * r_i_i_C(1) + (-t404 * t445 - t442 * t461) * r_i_i_C(2) + ((-t417 * t445 + t442 * t483) * r_i_i_C(1) + (t417 * t442 + t445 * t483) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
