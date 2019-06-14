% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 18:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:28:51
% EndTime: 2019-05-06 18:29:01
% DurationCPUTime: 7.76s
% Computational Cost: add. (105440->339), mult. (238956->429), div. (0->0), fcn. (191970->12), ass. (0->143)
t576 = Ifges(6,1) + Ifges(7,1);
t569 = Ifges(6,4) + Ifges(7,4);
t568 = Ifges(6,5) + Ifges(7,5);
t575 = Ifges(6,2) + Ifges(7,2);
t567 = Ifges(6,6) + Ifges(7,6);
t574 = Ifges(6,3) + Ifges(7,3);
t533 = cos(pkin(6));
t527 = qJD(1) * t533 + qJD(2);
t530 = sin(pkin(11));
t532 = cos(pkin(11));
t536 = sin(qJ(2));
t531 = sin(pkin(6));
t557 = qJD(1) * t531;
t552 = t536 * t557;
t509 = t527 * t532 - t530 * t552;
t510 = t527 * t530 + t532 * t552;
t535 = sin(qJ(4));
t539 = cos(qJ(4));
t491 = t509 * t539 - t510 * t535;
t540 = cos(qJ(2));
t556 = qJD(1) * t540;
t520 = (qJD(2) * t556 + qJDD(1) * t536) * t531;
t526 = qJDD(1) * t533 + qJDD(2);
t497 = -t520 * t530 + t526 * t532;
t498 = t520 * t532 + t526 * t530;
t460 = qJD(4) * t491 + t497 * t535 + t498 * t539;
t492 = t509 * t535 + t510 * t539;
t551 = t531 * t556;
t523 = qJD(4) - t551;
t534 = sin(qJ(5));
t538 = cos(qJ(5));
t479 = -t492 * t534 + t523 * t538;
t555 = qJDD(1) * t531;
t521 = -qJD(2) * t552 + t540 * t555;
t513 = qJDD(4) - t521;
t441 = qJD(5) * t479 + t460 * t538 + t513 * t534;
t480 = t492 * t538 + t523 * t534;
t454 = -mrSges(7,1) * t479 + mrSges(7,2) * t480;
t518 = (-pkin(2) * t540 - qJ(3) * t536) * t557;
t525 = t527 ^ 2;
t542 = qJD(1) ^ 2;
t537 = sin(qJ(1));
t541 = cos(qJ(1));
t546 = -g(1) * t541 - g(2) * t537;
t517 = -pkin(1) * t542 + pkin(8) * t555 + t546;
t550 = g(1) * t537 - g(2) * t541;
t572 = pkin(8) * t531;
t516 = qJDD(1) * pkin(1) + t542 * t572 + t550;
t565 = t516 * t533;
t558 = t517 * t540 + t536 * t565;
t476 = -t525 * pkin(2) + t526 * qJ(3) + (-g(3) * t536 + t518 * t556) * t531 + t558;
t571 = t533 * g(3);
t477 = -t521 * pkin(2) - t571 - t520 * qJ(3) + (-t516 + (pkin(2) * t536 - qJ(3) * t540) * t527 * qJD(1)) * t531;
t435 = -0.2e1 * qJD(3) * t510 - t530 * t476 + t477 * t532;
t431 = (-t509 * t551 - t498) * pkin(9) + (t509 * t510 - t521) * pkin(3) + t435;
t436 = 0.2e1 * qJD(3) * t509 + t476 * t532 + t477 * t530;
t499 = -pkin(3) * t551 - pkin(9) * t510;
t508 = t509 ^ 2;
t434 = -pkin(3) * t508 + pkin(9) * t497 + t499 * t551 + t436;
t426 = t431 * t535 + t434 * t539;
t471 = -pkin(4) * t491 - pkin(10) * t492;
t522 = t523 ^ 2;
t424 = -pkin(4) * t522 + pkin(10) * t513 + t471 * t491 + t426;
t563 = t531 * t540;
t488 = -g(3) * t563 - t517 * t536 + t540 * t565;
t475 = -t526 * pkin(2) - t525 * qJ(3) + t518 * t552 + qJDD(3) - t488;
t442 = -t497 * pkin(3) - t508 * pkin(9) + t499 * t510 + t475;
t459 = -qJD(4) * t492 + t497 * t539 - t498 * t535;
t429 = (-t491 * t523 - t460) * pkin(10) + (t492 * t523 - t459) * pkin(4) + t442;
t419 = -t534 * t424 + t429 * t538;
t458 = qJDD(5) - t459;
t490 = qJD(5) - t491;
t416 = -0.2e1 * qJD(6) * t480 + (t479 * t490 - t441) * qJ(6) + (t479 * t480 + t458) * pkin(5) + t419;
t461 = -mrSges(7,2) * t490 + mrSges(7,3) * t479;
t554 = m(7) * t416 + mrSges(7,1) * t458 + t461 * t490;
t413 = -t441 * mrSges(7,3) - t480 * t454 + t554;
t420 = t424 * t538 + t429 * t534;
t440 = -qJD(5) * t480 - t460 * t534 + t513 * t538;
t463 = pkin(5) * t490 - qJ(6) * t480;
t478 = t479 ^ 2;
t418 = -pkin(5) * t478 + qJ(6) * t440 + 0.2e1 * qJD(6) * t479 - t463 * t490 + t420;
t560 = t479 * t569 + t480 * t576 + t490 * t568;
t561 = -t479 * t575 - t480 * t569 - t490 * t567;
t573 = mrSges(6,1) * t419 + mrSges(7,1) * t416 - mrSges(6,2) * t420 - mrSges(7,2) * t418 + pkin(5) * t413 + t440 * t567 + t441 * t568 + t458 * t574 - t479 * t560 - t480 * t561;
t570 = -mrSges(6,2) - mrSges(7,2);
t564 = t531 * t536;
t455 = -mrSges(6,1) * t479 + mrSges(6,2) * t480;
t462 = -mrSges(6,2) * t490 + mrSges(6,3) * t479;
t406 = m(6) * t419 + t458 * mrSges(6,1) + t490 * t462 + (-t454 - t455) * t480 + (-mrSges(6,3) - mrSges(7,3)) * t441 + t554;
t553 = m(7) * t418 + mrSges(7,3) * t440 + t454 * t479;
t464 = mrSges(7,1) * t490 - mrSges(7,3) * t480;
t559 = -mrSges(6,1) * t490 + mrSges(6,3) * t480 - t464;
t411 = m(6) * t420 + t440 * mrSges(6,3) + t479 * t455 + t458 * t570 + t490 * t559 + t553;
t404 = -t406 * t534 + t411 * t538;
t470 = -mrSges(5,1) * t491 + mrSges(5,2) * t492;
t482 = mrSges(5,1) * t523 - mrSges(5,3) * t492;
t400 = m(5) * t426 - mrSges(5,2) * t513 + mrSges(5,3) * t459 + t470 * t491 - t482 * t523 + t404;
t425 = t431 * t539 - t434 * t535;
t423 = -pkin(4) * t513 - pkin(10) * t522 + t471 * t492 - t425;
t421 = -pkin(5) * t440 - qJ(6) * t478 + t463 * t480 + qJDD(6) + t423;
t547 = -m(7) * t421 + mrSges(7,1) * t440 + t461 * t479;
t412 = -m(6) * t423 + mrSges(6,1) * t440 + t441 * t570 + t462 * t479 + t480 * t559 + t547;
t481 = -mrSges(5,2) * t523 + mrSges(5,3) * t491;
t408 = m(5) * t425 + t513 * mrSges(5,1) - t460 * mrSges(5,3) - t492 * t470 + t523 * t481 + t412;
t395 = t400 * t535 + t408 * t539;
t493 = -mrSges(4,1) * t509 + mrSges(4,2) * t510;
t495 = mrSges(4,2) * t551 + mrSges(4,3) * t509;
t393 = m(4) * t435 - mrSges(4,1) * t521 - mrSges(4,3) * t498 - t493 * t510 - t495 * t551 + t395;
t496 = -mrSges(4,1) * t551 - mrSges(4,3) * t510;
t548 = t400 * t539 - t408 * t535;
t394 = m(4) * t436 + mrSges(4,2) * t521 + mrSges(4,3) * t497 + t493 * t509 + t496 * t551 + t548;
t388 = t393 * t532 + t394 * t530;
t403 = t406 * t538 + t411 * t534;
t562 = -t479 * t567 - t480 * t568 - t490 * t574;
t549 = -t393 * t530 + t394 * t532;
t545 = m(5) * t442 - mrSges(5,1) * t459 + mrSges(5,2) * t460 - t481 * t491 + t482 * t492 + t403;
t401 = m(4) * t475 - mrSges(4,1) * t497 + mrSges(4,2) * t498 - t495 * t509 + t496 * t510 + t545;
t414 = t441 * mrSges(7,2) + t480 * t464 - t547;
t396 = -mrSges(6,1) * t423 + mrSges(6,3) * t420 - mrSges(7,1) * t421 + mrSges(7,3) * t418 - pkin(5) * t414 + qJ(6) * t553 + (-qJ(6) * t464 + t560) * t490 + t562 * t480 + (-mrSges(7,2) * qJ(6) + t567) * t458 + t569 * t441 + t575 * t440;
t402 = mrSges(6,2) * t423 + mrSges(7,2) * t421 - mrSges(6,3) * t419 - mrSges(7,3) * t416 - qJ(6) * t413 + t440 * t569 + t441 * t576 + t458 * t568 - t479 * t562 + t490 * t561;
t467 = Ifges(5,4) * t492 + Ifges(5,2) * t491 + Ifges(5,6) * t523;
t468 = Ifges(5,1) * t492 + Ifges(5,4) * t491 + Ifges(5,5) * t523;
t543 = mrSges(5,1) * t425 - mrSges(5,2) * t426 + Ifges(5,5) * t460 + Ifges(5,6) * t459 + Ifges(5,3) * t513 + pkin(4) * t412 + pkin(10) * t404 + t538 * t396 + t534 * t402 + t492 * t467 - t491 * t468;
t519 = (-mrSges(3,1) * t540 + mrSges(3,2) * t536) * t557;
t515 = -mrSges(3,2) * t527 + mrSges(3,3) * t551;
t514 = mrSges(3,1) * t527 - mrSges(3,3) * t552;
t503 = -t531 * t516 - t571;
t502 = Ifges(3,5) * t527 + (Ifges(3,1) * t536 + Ifges(3,4) * t540) * t557;
t501 = Ifges(3,6) * t527 + (t536 * Ifges(3,4) + Ifges(3,2) * t540) * t557;
t500 = Ifges(3,3) * t527 + (Ifges(3,5) * t536 + Ifges(3,6) * t540) * t557;
t489 = -g(3) * t564 + t558;
t485 = Ifges(4,1) * t510 + Ifges(4,4) * t509 - Ifges(4,5) * t551;
t484 = Ifges(4,4) * t510 + Ifges(4,2) * t509 - Ifges(4,6) * t551;
t483 = Ifges(4,5) * t510 + Ifges(4,6) * t509 - Ifges(4,3) * t551;
t466 = Ifges(5,5) * t492 + Ifges(5,6) * t491 + Ifges(5,3) * t523;
t397 = m(3) * t488 + mrSges(3,1) * t526 - mrSges(3,3) * t520 + t515 * t527 - t519 * t552 - t401;
t389 = -mrSges(5,1) * t442 + mrSges(5,3) * t426 + Ifges(5,4) * t460 + Ifges(5,2) * t459 + Ifges(5,6) * t513 - pkin(4) * t403 - t492 * t466 + t523 * t468 - t573;
t387 = m(3) * t489 - mrSges(3,2) * t526 + mrSges(3,3) * t521 - t514 * t527 + t519 * t551 + t549;
t386 = mrSges(5,2) * t442 - mrSges(5,3) * t425 + Ifges(5,1) * t460 + Ifges(5,4) * t459 + Ifges(5,5) * t513 - pkin(10) * t403 - t396 * t534 + t402 * t538 + t466 * t491 - t467 * t523;
t385 = mrSges(4,2) * t475 - mrSges(4,3) * t435 + Ifges(4,1) * t498 + Ifges(4,4) * t497 - Ifges(4,5) * t521 - pkin(9) * t395 + t386 * t539 - t389 * t535 + t483 * t509 + t484 * t551;
t384 = -mrSges(4,1) * t475 + mrSges(4,3) * t436 + Ifges(4,4) * t498 + Ifges(4,2) * t497 - Ifges(4,6) * t521 - pkin(3) * t545 + pkin(9) * t548 + t535 * t386 + t539 * t389 - t510 * t483 - t485 * t551;
t383 = Ifges(3,5) * t520 + Ifges(3,6) * t521 + Ifges(3,3) * t526 + mrSges(3,1) * t488 - mrSges(3,2) * t489 + t530 * t385 + t532 * t384 - pkin(2) * t401 + qJ(3) * t549 + (t501 * t536 - t502 * t540) * t557;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t550 - mrSges(2,2) * t546 + (mrSges(3,2) * t503 - mrSges(3,3) * t488 + Ifges(3,1) * t520 + Ifges(3,4) * t521 + Ifges(3,5) * t526 - qJ(3) * t388 - t384 * t530 + t385 * t532 + t500 * t551 - t501 * t527) * t564 + (-pkin(2) * t388 - t500 * t552 - t543 - pkin(3) * t395 + t527 * t502 + Ifges(3,6) * t526 + Ifges(3,4) * t520 + t509 * t485 - t510 * t484 - Ifges(4,5) * t498 - mrSges(3,1) * t503 - Ifges(4,6) * t497 + mrSges(3,3) * t489 + (Ifges(3,2) + Ifges(4,3)) * t521 - mrSges(4,1) * t435 + mrSges(4,2) * t436) * t563 + t533 * t383 + pkin(1) * ((t387 * t536 + t397 * t540) * t533 + (-m(3) * t503 + t521 * mrSges(3,1) - t520 * mrSges(3,2) + (-t514 * t536 + t515 * t540) * t557 - t388) * t531) + (t387 * t540 - t397 * t536) * t572; t383; t401; t543; t573; t414;];
tauJ  = t1;
