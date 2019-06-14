% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR12
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
% Datum: 2019-05-06 16:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:14:59
% EndTime: 2019-05-06 16:15:08
% DurationCPUTime: 7.43s
% Computational Cost: add. (74002->346), mult. (169883->433), div. (0->0), fcn. (126121->12), ass. (0->151)
t619 = -2 * qJD(3);
t618 = Ifges(3,1) + Ifges(4,2);
t609 = Ifges(3,4) + Ifges(4,6);
t608 = Ifges(3,5) - Ifges(4,4);
t617 = Ifges(3,2) + Ifges(4,3);
t607 = Ifges(3,6) - Ifges(4,5);
t616 = Ifges(3,3) + Ifges(4,1);
t566 = cos(pkin(6));
t558 = qJD(1) * t566 + qJD(2);
t569 = sin(qJ(2));
t564 = sin(pkin(6));
t598 = qJD(1) * t564;
t592 = t569 * t598;
t615 = (pkin(2) * t558 + t619) * t592;
t575 = qJD(1) ^ 2;
t570 = sin(qJ(1));
t574 = cos(qJ(1));
t586 = -g(1) * t574 - t570 * g(2);
t595 = qJDD(1) * t564;
t541 = -pkin(1) * t575 + pkin(8) * t595 + t586;
t573 = cos(qJ(2));
t604 = t564 * t569;
t590 = t570 * g(1) - g(2) * t574;
t611 = t564 * pkin(8);
t540 = qJDD(1) * pkin(1) + t575 * t611 + t590;
t606 = t540 * t566;
t503 = -g(3) * t604 + t573 * t541 + t569 * t606;
t542 = (-t573 * pkin(2) - t569 * qJ(3)) * t598;
t556 = t558 ^ 2;
t557 = qJDD(1) * t566 + qJDD(2);
t597 = qJD(1) * t573;
t591 = t564 * t597;
t482 = pkin(2) * t556 - t557 * qJ(3) - t542 * t591 + t558 * t619 - t503;
t614 = 2 * qJD(5);
t613 = -pkin(2) - pkin(9);
t612 = g(3) * t566;
t610 = mrSges(3,1) - mrSges(4,2);
t605 = t564 ^ 2 * t575;
t603 = t564 * t573;
t545 = pkin(3) * t592 - pkin(9) * t558;
t546 = (qJD(2) * t597 + qJDD(1) * t569) * t564;
t547 = -qJD(2) * t592 + t573 * t595;
t593 = t573 ^ 2 * t605;
t464 = -pkin(3) * t593 - t612 - qJ(3) * t546 + t613 * t547 + (-t540 + (-qJ(3) * t558 * t573 - t545 * t569) * qJD(1)) * t564 + t615;
t599 = g(3) * t603 + t569 * t541;
t583 = -qJ(3) * t556 + t542 * t592 + qJDD(3) + t599;
t467 = pkin(3) * t546 + t613 * t557 + (-pkin(3) * t558 * t598 - pkin(9) * t569 * t605 - t606) * t573 + t583;
t568 = sin(qJ(4));
t572 = cos(qJ(4));
t449 = -t464 * t568 + t572 * t467;
t527 = -t558 * t568 - t572 * t591;
t501 = qJD(4) * t527 - t547 * t568 + t557 * t572;
t528 = t558 * t572 - t568 * t591;
t535 = qJDD(4) + t546;
t551 = qJD(4) + t592;
t445 = (t527 * t551 - t501) * qJ(5) + (t527 * t528 + t535) * pkin(4) + t449;
t450 = t572 * t464 + t568 * t467;
t500 = -qJD(4) * t528 - t547 * t572 - t557 * t568;
t510 = pkin(4) * t551 - qJ(5) * t528;
t526 = t527 ^ 2;
t447 = -pkin(4) * t526 + qJ(5) * t500 - t510 * t551 + t450;
t563 = sin(pkin(11));
t565 = cos(pkin(11));
t506 = t527 * t565 - t528 * t563;
t442 = t563 * t445 + t565 * t447 + t506 * t614;
t475 = t500 * t565 - t501 * t563;
t507 = t527 * t563 + t528 * t565;
t484 = -mrSges(6,1) * t506 + mrSges(6,2) * t507;
t491 = mrSges(6,1) * t551 - mrSges(6,3) * t507;
t485 = -pkin(5) * t506 - pkin(10) * t507;
t549 = t551 ^ 2;
t439 = -pkin(5) * t549 + pkin(10) * t535 + t485 * t506 + t442;
t463 = pkin(3) * t547 - pkin(9) * t593 + t558 * t545 - t482;
t452 = -pkin(4) * t500 - qJ(5) * t526 + t528 * t510 + qJDD(5) + t463;
t476 = t500 * t563 + t501 * t565;
t443 = (-t506 * t551 - t476) * pkin(10) + (t507 * t551 - t475) * pkin(5) + t452;
t567 = sin(qJ(6));
t571 = cos(qJ(6));
t436 = -t439 * t567 + t443 * t571;
t488 = -t507 * t567 + t551 * t571;
t455 = qJD(6) * t488 + t476 * t571 + t535 * t567;
t489 = t507 * t571 + t551 * t567;
t468 = -mrSges(7,1) * t488 + mrSges(7,2) * t489;
t505 = qJD(6) - t506;
t469 = -mrSges(7,2) * t505 + mrSges(7,3) * t488;
t474 = qJDD(6) - t475;
t433 = m(7) * t436 + mrSges(7,1) * t474 - mrSges(7,3) * t455 - t468 * t489 + t469 * t505;
t437 = t439 * t571 + t443 * t567;
t454 = -qJD(6) * t489 - t476 * t567 + t535 * t571;
t470 = mrSges(7,1) * t505 - mrSges(7,3) * t489;
t434 = m(7) * t437 - mrSges(7,2) * t474 + mrSges(7,3) * t454 + t468 * t488 - t470 * t505;
t587 = -t433 * t567 + t571 * t434;
t420 = m(6) * t442 - mrSges(6,2) * t535 + mrSges(6,3) * t475 + t484 * t506 - t491 * t551 + t587;
t585 = -t445 * t565 + t447 * t563;
t441 = -0.2e1 * qJD(5) * t507 - t585;
t490 = -mrSges(6,2) * t551 + mrSges(6,3) * t506;
t438 = -pkin(5) * t535 - pkin(10) * t549 + (t614 + t485) * t507 + t585;
t581 = -m(7) * t438 + t454 * mrSges(7,1) - mrSges(7,2) * t455 + t488 * t469 - t470 * t489;
t429 = m(6) * t441 + mrSges(6,1) * t535 - mrSges(6,3) * t476 - t484 * t507 + t490 * t551 + t581;
t416 = t563 * t420 + t565 * t429;
t423 = t571 * t433 + t567 * t434;
t602 = (t569 * t608 + t573 * t607) * t598 + t616 * t558;
t601 = (t569 * t609 + t573 * t617) * t598 + t607 * t558;
t600 = (t569 * t618 + t573 * t609) * t598 + t608 * t558;
t594 = t573 * t606;
t508 = -mrSges(5,1) * t527 + mrSges(5,2) * t528;
t509 = -mrSges(5,2) * t551 + mrSges(5,3) * t527;
t413 = m(5) * t449 + mrSges(5,1) * t535 - mrSges(5,3) * t501 - t508 * t528 + t509 * t551 + t416;
t511 = mrSges(5,1) * t551 - mrSges(5,3) * t528;
t588 = t565 * t420 - t429 * t563;
t414 = m(5) * t450 - mrSges(5,2) * t535 + mrSges(5,3) * t500 + t508 * t527 - t511 * t551 + t588;
t589 = -t568 * t413 + t572 * t414;
t518 = -t540 * t564 - t612;
t409 = t413 * t572 + t414 * t568;
t483 = -pkin(2) * t547 + (-t558 * t591 - t546) * qJ(3) + t518 + t615;
t538 = -mrSges(4,1) * t591 - mrSges(4,3) * t558;
t584 = -m(4) * t483 + t546 * mrSges(4,3) - t538 * t591 - t589;
t486 = -pkin(2) * t557 + t583 - t594;
t582 = -m(4) * t486 - t546 * mrSges(4,1) - t409;
t421 = m(6) * t452 - t475 * mrSges(6,1) + t476 * mrSges(6,2) - t506 * t490 + t507 * t491 + t423;
t457 = Ifges(7,4) * t489 + Ifges(7,2) * t488 + Ifges(7,6) * t505;
t458 = Ifges(7,1) * t489 + Ifges(7,4) * t488 + Ifges(7,5) * t505;
t579 = mrSges(7,1) * t436 - mrSges(7,2) * t437 + Ifges(7,5) * t455 + Ifges(7,6) * t454 + Ifges(7,3) * t474 + t457 * t489 - t458 * t488;
t578 = -m(5) * t463 + t500 * mrSges(5,1) - t501 * mrSges(5,2) + t527 * t509 - t528 * t511 - t421;
t539 = mrSges(4,1) * t592 + mrSges(4,2) * t558;
t543 = (t573 * mrSges(4,2) - t569 * mrSges(4,3)) * t598;
t577 = -m(4) * t482 + t557 * mrSges(4,3) + t558 * t539 + t543 * t591 - t578;
t456 = Ifges(7,5) * t489 + Ifges(7,6) * t488 + Ifges(7,3) * t505;
t426 = -mrSges(7,1) * t438 + mrSges(7,3) * t437 + Ifges(7,4) * t455 + Ifges(7,2) * t454 + Ifges(7,6) * t474 - t456 * t489 + t458 * t505;
t427 = mrSges(7,2) * t438 - mrSges(7,3) * t436 + Ifges(7,1) * t455 + Ifges(7,4) * t454 + Ifges(7,5) * t474 + t456 * t488 - t457 * t505;
t478 = Ifges(6,4) * t507 + Ifges(6,2) * t506 + Ifges(6,6) * t551;
t479 = Ifges(6,1) * t507 + Ifges(6,4) * t506 + Ifges(6,5) * t551;
t493 = Ifges(5,4) * t528 + Ifges(5,2) * t527 + Ifges(5,6) * t551;
t494 = Ifges(5,1) * t528 + Ifges(5,4) * t527 + Ifges(5,5) * t551;
t576 = mrSges(5,1) * t449 + mrSges(6,1) * t441 - mrSges(5,2) * t450 - mrSges(6,2) * t442 + pkin(4) * t416 + pkin(5) * t581 + pkin(10) * t587 + t571 * t426 + t567 * t427 + t507 * t478 - t506 * t479 - t527 * t494 + Ifges(6,6) * t475 + Ifges(6,5) * t476 + t528 * t493 + Ifges(5,6) * t500 + Ifges(5,5) * t501 + (Ifges(6,3) + Ifges(5,3)) * t535;
t544 = (-t573 * mrSges(3,1) + t569 * mrSges(3,2)) * t598;
t537 = -mrSges(3,2) * t558 + mrSges(3,3) * t591;
t536 = mrSges(3,1) * t558 - mrSges(3,3) * t592;
t502 = t594 - t599;
t492 = Ifges(5,5) * t528 + Ifges(5,6) * t527 + Ifges(5,3) * t551;
t477 = Ifges(6,5) * t507 + Ifges(6,6) * t506 + Ifges(6,3) * t551;
t417 = (mrSges(3,3) + mrSges(4,1)) * t547 + t544 * t591 + t577 - t558 * t536 - t557 * mrSges(3,2) + m(3) * t503;
t411 = -mrSges(6,1) * t452 + mrSges(6,3) * t442 + Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t535 - pkin(5) * t423 - t477 * t507 + t479 * t551 - t579;
t410 = mrSges(6,2) * t452 - mrSges(6,3) * t441 + Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t535 - pkin(10) * t423 - t426 * t567 + t427 * t571 + t477 * t506 - t478 * t551;
t408 = mrSges(4,2) * t557 + t538 * t558 + t543 * t592 - t582;
t407 = t547 * mrSges(4,2) - t539 * t592 - t584;
t406 = m(3) * t502 - mrSges(3,3) * t546 + (t537 - t538) * t558 + t610 * t557 + (-t543 - t544) * t592 + t582;
t405 = mrSges(5,2) * t463 - mrSges(5,3) * t449 + Ifges(5,1) * t501 + Ifges(5,4) * t500 + Ifges(5,5) * t535 - qJ(5) * t416 + t410 * t565 - t411 * t563 + t492 * t527 - t493 * t551;
t404 = -mrSges(5,1) * t463 + mrSges(5,3) * t450 + Ifges(5,4) * t501 + Ifges(5,2) * t500 + Ifges(5,6) * t535 - pkin(4) * t421 + qJ(5) * t588 + t563 * t410 + t565 * t411 - t528 * t492 + t551 * t494;
t403 = mrSges(3,1) * t502 - mrSges(3,2) * t503 + mrSges(4,2) * t486 - mrSges(4,3) * t482 + t572 * t405 - t568 * t404 - pkin(9) * t409 - pkin(2) * t408 + qJ(3) * t577 + t616 * t557 + (qJ(3) * mrSges(4,1) + t607) * t547 + t608 * t546 + (t601 * t569 - t600 * t573) * t598;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t590 - mrSges(2,2) * t586 + (mrSges(4,1) * t486 + mrSges(3,2) * t518 - mrSges(3,3) * t502 - mrSges(4,3) * t483 + pkin(3) * t409 - qJ(3) * t407 + t618 * t546 + t609 * t547 + t608 * t557 - t601 * t558 + t602 * t591 + t576) * t604 + (-mrSges(3,1) * t518 - mrSges(4,1) * t482 + mrSges(4,2) * t483 + mrSges(3,3) * t503 - pkin(2) * t407 - pkin(3) * t578 - pkin(9) * t589 - t572 * t404 - t568 * t405 + t609 * t546 + t617 * t547 + t607 * t557 + t600 * t558 - t602 * t592) * t603 + t566 * t403 + pkin(1) * ((t573 * t406 + t569 * t417) * t566 + (-m(3) * t518 - t546 * mrSges(3,2) + t610 * t547 + (t537 * t573 + (-t536 + t539) * t569) * t598 + t584) * t564) + (-t406 * t569 + t417 * t573) * t611; t403; t408; t576; t421; t579;];
tauJ  = t1;
