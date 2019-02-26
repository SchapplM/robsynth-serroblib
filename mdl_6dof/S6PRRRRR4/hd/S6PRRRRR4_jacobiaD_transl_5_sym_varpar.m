% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:28
% EndTime: 2019-02-26 20:20:28
% DurationCPUTime: 0.66s
% Computational Cost: add. (695->129), mult. (1858->229), div. (0->0), fcn. (1972->14), ass. (0->83)
t442 = sin(qJ(3));
t445 = cos(qJ(3));
t435 = sin(pkin(13));
t438 = cos(pkin(13));
t446 = cos(qJ(2));
t440 = cos(pkin(6));
t443 = sin(qJ(2));
t481 = t440 * t443;
t456 = t435 * t481 - t438 * t446;
t439 = cos(pkin(7));
t480 = t440 * t446;
t457 = t435 * t480 + t438 * t443;
t436 = sin(pkin(7));
t437 = sin(pkin(6));
t487 = t436 * t437;
t458 = t435 * t487 - t439 * t457;
t409 = t442 * t458 - t445 * t456;
t425 = t435 * t446 + t438 * t481;
t468 = t438 * t480;
t424 = -t435 * t443 + t468;
t459 = -t424 * t439 + t438 * t487;
t494 = -t425 * t445 + t442 * t459;
t493 = r_i_i_C(3) + pkin(11) + pkin(10);
t434 = qJ(4) + qJ(5);
t431 = sin(t434);
t433 = qJD(4) + qJD(5);
t490 = t431 * t433;
t432 = cos(t434);
t489 = t432 * t433;
t488 = t433 * t443;
t486 = t436 * t440;
t444 = cos(qJ(4));
t485 = t436 * t444;
t484 = t437 * t439;
t483 = t439 * t442;
t482 = t439 * t445;
t479 = t442 * t443;
t478 = t442 * t446;
t477 = t443 * t445;
t476 = t445 * t446;
t420 = t425 * qJD(2);
t462 = -t420 * t436 - t433 * t494;
t472 = qJD(2) * t443;
t419 = -qJD(2) * t468 + t435 * t472;
t449 = -t425 * t442 - t445 * t459;
t397 = qJD(3) * t449 - t419 * t445 - t420 * t483;
t416 = -t424 * t436 - t438 * t484;
t466 = -t416 * t433 - t397;
t475 = (t431 * t466 - t432 * t462) * r_i_i_C(1) + (t431 * t462 + t432 * t466) * r_i_i_C(2);
t422 = t456 * qJD(2);
t461 = t409 * t433 + t422 * t436;
t421 = t457 * qJD(2);
t448 = t442 * t456 + t445 * t458;
t399 = qJD(3) * t448 - t421 * t445 + t422 * t483;
t417 = t435 * t484 + t436 * t457;
t465 = -t417 * t433 - t399;
t474 = (t431 * t465 - t432 * t461) * r_i_i_C(1) + (t431 * t461 + t432 * t465) * r_i_i_C(2);
t453 = t439 * t478 + t477;
t415 = t437 * t453 + t442 * t486;
t463 = t472 * t487;
t451 = -t415 * t433 + t463;
t452 = t439 * t479 - t476;
t455 = t439 * t476 - t479;
t467 = qJD(3) * t486;
t405 = t445 * t467 + (-qJD(2) * t452 + qJD(3) * t455) * t437;
t423 = t440 * t439 - t446 * t487;
t464 = -t423 * t433 - t405;
t473 = (t431 * t464 + t432 * t451) * r_i_i_C(1) + (-t431 * t451 + t432 * t464) * r_i_i_C(2);
t471 = qJD(2) * t446;
t470 = qJD(4) * t444;
t441 = sin(qJ(4));
t469 = qJD(4) * t441 * pkin(4);
t430 = t444 * pkin(4) + pkin(3);
t460 = t432 * r_i_i_C(1) - t431 * r_i_i_C(2) + t430;
t412 = t424 * t445 - t425 * t483;
t413 = -t445 * t457 + t456 * t483;
t454 = -t439 * t477 - t478;
t450 = -t469 + (-t431 * r_i_i_C(1) - t432 * r_i_i_C(2)) * t433;
t418 = t452 * t437;
t411 = (-qJD(2) * t453 + qJD(3) * t454) * t437;
t403 = t421 * t483 + t422 * t445 + (t442 * t457 + t456 * t482) * qJD(3);
t401 = t419 * t483 - t420 * t445 + (-t424 * t442 - t425 * t482) * qJD(3);
t1 = [0 (t403 * t432 - t413 * t490) * r_i_i_C(1) + (-t403 * t431 - t413 * t489) * r_i_i_C(2) + t403 * t430 - t413 * t469 + t422 * pkin(2) + t493 * (qJD(3) * t413 - t421 * t482 + t422 * t442) + ((-t421 * t431 - t456 * t489) * r_i_i_C(1) + (-t421 * t432 + t456 * t490) * r_i_i_C(2) - t421 * pkin(9) + (-t421 * t441 - t456 * t470) * pkin(4)) * t436, t493 * t399 + t450 * t448 + t460 * (-t409 * qJD(3) + t421 * t442 + t422 * t482) (-t422 * t485 - t399 * t441 + (-t409 * t444 - t417 * t441) * qJD(4)) * pkin(4) + t474, t474, 0; 0 (t401 * t432 - t412 * t490) * r_i_i_C(1) + (-t401 * t431 - t412 * t489) * r_i_i_C(2) + t401 * t430 - t412 * t469 - t420 * pkin(2) + t493 * (qJD(3) * t412 - t419 * t482 - t420 * t442) + ((-t419 * t431 + t425 * t489) * r_i_i_C(1) + (-t419 * t432 - t425 * t490) * r_i_i_C(2) - t419 * pkin(9) + (-t419 * t441 + t425 * t470) * pkin(4)) * t436, t493 * t397 + t450 * t449 + t460 * (t494 * qJD(3) + t419 * t442 - t420 * t482) (t420 * t485 - t397 * t441 + (-t416 * t441 + t444 * t494) * qJD(4)) * pkin(4) + t475, t475, 0; 0 (t411 * t432 + t418 * t490) * r_i_i_C(1) + (-t411 * t431 + t418 * t489) * r_i_i_C(2) + t411 * t430 + t418 * t469 + (-t493 * (-qJD(2) * t455 + qJD(3) * t452) - pkin(2) * t472 + ((t431 * t471 + t432 * t488) * r_i_i_C(1) + (-t431 * t488 + t432 * t471) * r_i_i_C(2) + pkin(9) * t471 + (t441 * t471 + t443 * t470) * pkin(4)) * t436) * t437, t493 * t405 + t450 * (t437 * t455 + t445 * t486) + t460 * (-t442 * t467 + (qJD(2) * t454 - qJD(3) * t453) * t437) (t444 * t463 - t405 * t441 + (-t415 * t444 - t423 * t441) * qJD(4)) * pkin(4) + t473, t473, 0;];
JaD_transl  = t1;
