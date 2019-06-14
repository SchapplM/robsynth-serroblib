% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:33:05
% EndTime: 2019-05-06 16:33:14
% DurationCPUTime: 7.47s
% Computational Cost: add. (76457->345), mult. (173332->433), div. (0->0), fcn. (127825->12), ass. (0->149)
t614 = -2 * qJD(3);
t613 = Ifges(3,1) + Ifges(4,2);
t605 = Ifges(3,4) + Ifges(4,6);
t604 = Ifges(3,5) - Ifges(4,4);
t612 = Ifges(3,2) + Ifges(4,3);
t603 = Ifges(3,6) - Ifges(4,5);
t611 = Ifges(3,3) + Ifges(4,1);
t563 = cos(pkin(6));
t555 = qJD(1) * t563 + qJD(2);
t566 = sin(qJ(2));
t561 = sin(pkin(6));
t594 = qJD(1) * t561;
t588 = t566 * t594;
t610 = (pkin(2) * t555 + t614) * t588;
t572 = qJD(1) ^ 2;
t567 = sin(qJ(1));
t571 = cos(qJ(1));
t582 = -g(1) * t571 - g(2) * t567;
t591 = qJDD(1) * t561;
t538 = -pkin(1) * t572 + pkin(8) * t591 + t582;
t570 = cos(qJ(2));
t600 = t561 * t566;
t586 = t567 * g(1) - g(2) * t571;
t608 = pkin(8) * t561;
t537 = qJDD(1) * pkin(1) + t572 * t608 + t586;
t602 = t537 * t563;
t499 = -g(3) * t600 + t570 * t538 + t566 * t602;
t539 = (-pkin(2) * t570 - qJ(3) * t566) * t594;
t553 = t555 ^ 2;
t554 = qJDD(1) * t563 + qJDD(2);
t593 = qJD(1) * t570;
t587 = t561 * t593;
t470 = pkin(2) * t553 - t554 * qJ(3) - t539 * t587 + t555 * t614 - t499;
t609 = -pkin(2) - pkin(9);
t607 = g(3) * t563;
t606 = mrSges(3,1) - mrSges(4,2);
t601 = t561 ^ 2 * t572;
t599 = t561 * t570;
t542 = pkin(3) * t588 - pkin(9) * t555;
t543 = (qJD(2) * t593 + qJDD(1) * t566) * t561;
t544 = -qJD(2) * t588 + t570 * t591;
t589 = t570 ^ 2 * t601;
t464 = -pkin(3) * t589 - t607 - qJ(3) * t543 + t609 * t544 + (-t537 + (-qJ(3) * t555 * t570 - t542 * t566) * qJD(1)) * t561 + t610;
t595 = g(3) * t599 + t566 * t538;
t580 = -qJ(3) * t553 + t539 * t588 + qJDD(3) + t595;
t466 = pkin(3) * t543 + t609 * t554 + (-pkin(3) * t555 * t594 - pkin(9) * t566 * t601 - t602) * t570 + t580;
t565 = sin(qJ(4));
t569 = cos(qJ(4));
t451 = t569 * t464 + t565 * t466;
t524 = t555 * t565 + t569 * t587;
t525 = t555 * t569 - t565 * t587;
t500 = pkin(4) * t524 - qJ(5) * t525;
t532 = qJDD(4) + t543;
t548 = qJD(4) + t588;
t546 = t548 ^ 2;
t445 = -pkin(4) * t546 + qJ(5) * t532 - t500 * t524 + t451;
t463 = pkin(3) * t544 - pkin(9) * t589 + t555 * t542 - t470;
t496 = qJD(4) * t525 + t569 * t544 + t554 * t565;
t497 = -qJD(4) * t524 - t544 * t565 + t554 * t569;
t448 = (t524 * t548 - t497) * qJ(5) + (t525 * t548 + t496) * pkin(4) + t463;
t560 = sin(pkin(11));
t562 = cos(pkin(11));
t506 = t525 * t562 + t548 * t560;
t440 = -0.2e1 * qJD(5) * t506 - t445 * t560 + t562 * t448;
t481 = t497 * t562 + t532 * t560;
t505 = -t525 * t560 + t548 * t562;
t438 = (t505 * t524 - t481) * pkin(10) + (t505 * t506 + t496) * pkin(5) + t440;
t441 = 0.2e1 * qJD(5) * t505 + t562 * t445 + t560 * t448;
t480 = -t497 * t560 + t532 * t562;
t487 = pkin(5) * t524 - pkin(10) * t506;
t504 = t505 ^ 2;
t439 = -pkin(5) * t504 + pkin(10) * t480 - t487 * t524 + t441;
t564 = sin(qJ(6));
t568 = cos(qJ(6));
t436 = t438 * t568 - t439 * t564;
t477 = t505 * t568 - t506 * t564;
t454 = qJD(6) * t477 + t480 * t564 + t481 * t568;
t478 = t505 * t564 + t506 * t568;
t459 = -mrSges(7,1) * t477 + mrSges(7,2) * t478;
t523 = qJD(6) + t524;
t467 = -mrSges(7,2) * t523 + mrSges(7,3) * t477;
t494 = qJDD(6) + t496;
t432 = m(7) * t436 + mrSges(7,1) * t494 - mrSges(7,3) * t454 - t459 * t478 + t467 * t523;
t437 = t438 * t564 + t439 * t568;
t453 = -qJD(6) * t478 + t480 * t568 - t481 * t564;
t468 = mrSges(7,1) * t523 - mrSges(7,3) * t478;
t433 = m(7) * t437 - mrSges(7,2) * t494 + mrSges(7,3) * t453 + t459 * t477 - t468 * t523;
t425 = t568 * t432 + t564 * t433;
t482 = -mrSges(6,1) * t505 + mrSges(6,2) * t506;
t485 = -mrSges(6,2) * t524 + mrSges(6,3) * t505;
t423 = m(6) * t440 + mrSges(6,1) * t496 - mrSges(6,3) * t481 - t482 * t506 + t485 * t524 + t425;
t486 = mrSges(6,1) * t524 - mrSges(6,3) * t506;
t583 = -t432 * t564 + t568 * t433;
t424 = m(6) * t441 - mrSges(6,2) * t496 + mrSges(6,3) * t480 + t482 * t505 - t486 * t524 + t583;
t419 = t562 * t423 + t560 * t424;
t598 = (t566 * t604 + t570 * t603) * t594 + t611 * t555;
t597 = (t566 * t613 + t570 * t605) * t594 + t604 * t555;
t596 = (-t605 * t566 - t612 * t570) * t594 - t603 * t555;
t590 = t570 * t602;
t501 = mrSges(5,1) * t524 + mrSges(5,2) * t525;
t508 = mrSges(5,1) * t548 - mrSges(5,3) * t525;
t584 = -t423 * t560 + t562 * t424;
t417 = m(5) * t451 - mrSges(5,2) * t532 - mrSges(5,3) * t496 - t501 * t524 - t508 * t548 + t584;
t450 = -t565 * t464 + t466 * t569;
t444 = -pkin(4) * t532 - qJ(5) * t546 + t525 * t500 + qJDD(5) - t450;
t442 = -pkin(5) * t480 - pkin(10) * t504 + t487 * t506 + t444;
t578 = m(7) * t442 - t453 * mrSges(7,1) + mrSges(7,2) * t454 - t477 * t467 + t468 * t478;
t435 = m(6) * t444 - t480 * mrSges(6,1) + mrSges(6,2) * t481 - t505 * t485 + t486 * t506 + t578;
t507 = -mrSges(5,2) * t548 - mrSges(5,3) * t524;
t428 = m(5) * t450 + mrSges(5,1) * t532 - mrSges(5,3) * t497 - t501 * t525 + t507 * t548 - t435;
t585 = t569 * t417 - t565 * t428;
t515 = -t537 * t561 - t607;
t410 = t417 * t565 + t428 * t569;
t471 = -pkin(2) * t544 + (-t555 * t587 - t543) * qJ(3) + t515 + t610;
t535 = -mrSges(4,1) * t587 - mrSges(4,3) * t555;
t581 = -m(4) * t471 + t543 * mrSges(4,3) - t535 * t587 - t585;
t476 = -pkin(2) * t554 + t580 - t590;
t579 = -m(4) * t476 - t543 * mrSges(4,1) - t410;
t576 = m(5) * t463 + t496 * mrSges(5,1) + t497 * mrSges(5,2) + t524 * t507 + t525 * t508 + t419;
t455 = Ifges(7,5) * t478 + Ifges(7,6) * t477 + Ifges(7,3) * t523;
t457 = Ifges(7,1) * t478 + Ifges(7,4) * t477 + Ifges(7,5) * t523;
t426 = -mrSges(7,1) * t442 + mrSges(7,3) * t437 + Ifges(7,4) * t454 + Ifges(7,2) * t453 + Ifges(7,6) * t494 - t455 * t478 + t457 * t523;
t456 = Ifges(7,4) * t478 + Ifges(7,2) * t477 + Ifges(7,6) * t523;
t427 = mrSges(7,2) * t442 - mrSges(7,3) * t436 + Ifges(7,1) * t454 + Ifges(7,4) * t453 + Ifges(7,5) * t494 + t455 * t477 - t456 * t523;
t472 = Ifges(6,5) * t506 + Ifges(6,6) * t505 + Ifges(6,3) * t524;
t474 = Ifges(6,1) * t506 + Ifges(6,4) * t505 + Ifges(6,5) * t524;
t412 = -mrSges(6,1) * t444 + mrSges(6,3) * t441 + Ifges(6,4) * t481 + Ifges(6,2) * t480 + Ifges(6,6) * t496 - pkin(5) * t578 + pkin(10) * t583 + t568 * t426 + t564 * t427 - t506 * t472 + t524 * t474;
t473 = Ifges(6,4) * t506 + Ifges(6,2) * t505 + Ifges(6,6) * t524;
t414 = mrSges(6,2) * t444 - mrSges(6,3) * t440 + Ifges(6,1) * t481 + Ifges(6,4) * t480 + Ifges(6,5) * t496 - pkin(10) * t425 - t426 * t564 + t427 * t568 + t472 * t505 - t473 * t524;
t489 = Ifges(5,4) * t525 - Ifges(5,2) * t524 + Ifges(5,6) * t548;
t490 = Ifges(5,1) * t525 - Ifges(5,4) * t524 + Ifges(5,5) * t548;
t575 = mrSges(5,1) * t450 - mrSges(5,2) * t451 + Ifges(5,5) * t497 - Ifges(5,6) * t496 + Ifges(5,3) * t532 - pkin(4) * t435 + qJ(5) * t584 + t562 * t412 + t560 * t414 + t525 * t489 + t524 * t490;
t536 = mrSges(4,1) * t588 + mrSges(4,2) * t555;
t540 = (mrSges(4,2) * t570 - mrSges(4,3) * t566) * t594;
t574 = -m(4) * t470 + t554 * mrSges(4,3) + t555 * t536 + t540 * t587 + t576;
t573 = mrSges(7,1) * t436 - mrSges(7,2) * t437 + Ifges(7,5) * t454 + Ifges(7,6) * t453 + Ifges(7,3) * t494 + t478 * t456 - t477 * t457;
t541 = (-mrSges(3,1) * t570 + mrSges(3,2) * t566) * t594;
t534 = -mrSges(3,2) * t555 + mrSges(3,3) * t587;
t533 = mrSges(3,1) * t555 - mrSges(3,3) * t588;
t498 = t590 - t595;
t488 = Ifges(5,5) * t525 - Ifges(5,6) * t524 + Ifges(5,3) * t548;
t415 = t574 + m(3) * t499 - mrSges(3,2) * t554 - t533 * t555 + (mrSges(3,3) + mrSges(4,1)) * t544 + t541 * t587;
t409 = mrSges(4,2) * t554 + t535 * t555 + t540 * t588 - t579;
t408 = t544 * mrSges(4,2) - t536 * t588 - t581;
t407 = m(3) * t498 - mrSges(3,3) * t543 + (t534 - t535) * t555 + t606 * t554 + (-t540 - t541) * t588 + t579;
t406 = -t573 + t548 * t490 + Ifges(5,6) * t532 - t525 * t488 + t505 * t474 - t506 * t473 + Ifges(5,4) * t497 - Ifges(6,6) * t480 - Ifges(6,5) * t481 - mrSges(5,1) * t463 + mrSges(5,3) * t451 + mrSges(6,2) * t441 - mrSges(6,1) * t440 - pkin(5) * t425 - pkin(4) * t419 + (-Ifges(5,2) - Ifges(6,3)) * t496;
t405 = mrSges(5,2) * t463 - mrSges(5,3) * t450 + Ifges(5,1) * t497 - Ifges(5,4) * t496 + Ifges(5,5) * t532 - qJ(5) * t419 - t412 * t560 + t414 * t562 - t488 * t524 - t489 * t548;
t404 = mrSges(3,1) * t498 - mrSges(3,2) * t499 + mrSges(4,2) * t476 - mrSges(4,3) * t470 + t569 * t405 - t565 * t406 - pkin(9) * t410 - pkin(2) * t409 + qJ(3) * t574 + t611 * t554 + (qJ(3) * mrSges(4,1) + t603) * t544 + t604 * t543 + (-t566 * t596 - t570 * t597) * t594;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t586 - mrSges(2,2) * t582 + (mrSges(4,1) * t476 + mrSges(3,2) * t515 - mrSges(3,3) * t498 - mrSges(4,3) * t471 + pkin(3) * t410 - qJ(3) * t408 + t613 * t543 + t605 * t544 + t604 * t554 + t596 * t555 + t598 * t587 + t575) * t600 + (-mrSges(3,1) * t515 - mrSges(4,1) * t470 + mrSges(4,2) * t471 + mrSges(3,3) * t499 - pkin(2) * t408 + pkin(3) * t576 - pkin(9) * t585 - t565 * t405 - t569 * t406 + t605 * t543 + t612 * t544 + t603 * t554 + t597 * t555 - t598 * t588) * t599 + t563 * t404 + pkin(1) * ((t407 * t570 + t415 * t566) * t563 + (-m(3) * t515 - t543 * mrSges(3,2) + t606 * t544 + (t534 * t570 + (-t533 + t536) * t566) * t594 + t581) * t561) + (-t407 * t566 + t415 * t570) * t608; t404; t409; t575; t435; t573;];
tauJ  = t1;
