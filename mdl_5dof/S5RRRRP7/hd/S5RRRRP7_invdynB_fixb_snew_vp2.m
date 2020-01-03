% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:28
% EndTime: 2019-12-31 21:56:34
% DurationCPUTime: 4.22s
% Computational Cost: add. (41735->290), mult. (84270->357), div. (0->0), fcn. (56406->8), ass. (0->113)
t598 = Ifges(5,1) + Ifges(6,1);
t591 = Ifges(5,4) - Ifges(6,5);
t597 = -Ifges(5,5) - Ifges(6,4);
t596 = Ifges(5,2) + Ifges(6,3);
t589 = Ifges(5,6) - Ifges(6,6);
t595 = -Ifges(5,3) - Ifges(6,2);
t562 = sin(qJ(3));
t563 = sin(qJ(2));
t565 = cos(qJ(3));
t566 = cos(qJ(2));
t539 = (t562 * t563 - t565 * t566) * qJD(1);
t594 = cos(qJ(4));
t568 = qJD(1) ^ 2;
t593 = pkin(2) * t568;
t592 = -mrSges(5,3) - mrSges(6,2);
t564 = sin(qJ(1));
t567 = cos(qJ(1));
t553 = -t567 * g(1) - t564 * g(2);
t542 = -t568 * pkin(1) + qJDD(1) * pkin(6) + t553;
t588 = t563 * t542;
t580 = qJD(1) * qJD(2);
t547 = t563 * qJDD(1) + t566 * t580;
t507 = qJDD(2) * pkin(2) - t547 * pkin(7) - t588 + (pkin(7) * t580 + t563 * t593 - g(3)) * t566;
t530 = -t563 * g(3) + t566 * t542;
t548 = t566 * qJDD(1) - t563 * t580;
t582 = qJD(1) * t563;
t551 = qJD(2) * pkin(2) - pkin(7) * t582;
t560 = t566 ^ 2;
t508 = t548 * pkin(7) - qJD(2) * t551 - t560 * t593 + t530;
t486 = t562 * t507 + t565 * t508;
t540 = (t562 * t566 + t563 * t565) * qJD(1);
t514 = -t540 * qJD(3) - t562 * t547 + t565 * t548;
t525 = t539 * mrSges(4,1) + t540 * mrSges(4,2);
t559 = qJD(2) + qJD(3);
t532 = t559 * mrSges(4,1) - t540 * mrSges(4,3);
t558 = qJDD(2) + qJDD(3);
t515 = -t539 * qJD(3) + t565 * t547 + t562 * t548;
t552 = t564 * g(1) - t567 * g(2);
t572 = -qJDD(1) * pkin(1) - t552;
t516 = -t548 * pkin(2) + t551 * t582 + (-pkin(7) * t560 - pkin(6)) * t568 + t572;
t481 = (t539 * t559 - t515) * pkin(8) + (t540 * t559 - t514) * pkin(3) + t516;
t526 = t539 * pkin(3) - t540 * pkin(8);
t557 = t559 ^ 2;
t484 = -t557 * pkin(3) + t558 * pkin(8) - t539 * t526 + t486;
t561 = sin(qJ(4));
t479 = t561 * t481 + t594 * t484;
t528 = t594 * t540 + t561 * t559;
t489 = t528 * qJD(4) + t561 * t515 - t594 * t558;
t513 = qJDD(4) - t514;
t535 = qJD(4) + t539;
t519 = t535 * mrSges(5,1) - t528 * mrSges(5,3);
t527 = t561 * t540 - t594 * t559;
t504 = t527 * pkin(4) - t528 * qJ(5);
t534 = t535 ^ 2;
t475 = -t534 * pkin(4) + t513 * qJ(5) + 0.2e1 * qJD(5) * t535 - t527 * t504 + t479;
t520 = -t535 * mrSges(6,1) + t528 * mrSges(6,2);
t579 = m(6) * t475 + t513 * mrSges(6,3) + t535 * t520;
t505 = t527 * mrSges(6,1) - t528 * mrSges(6,3);
t583 = -t527 * mrSges(5,1) - t528 * mrSges(5,2) - t505;
t470 = m(5) * t479 - t513 * mrSges(5,2) + t592 * t489 - t535 * t519 + t583 * t527 + t579;
t478 = t594 * t481 - t561 * t484;
t490 = -t527 * qJD(4) + t594 * t515 + t561 * t558;
t518 = -t535 * mrSges(5,2) - t527 * mrSges(5,3);
t476 = -t513 * pkin(4) - t534 * qJ(5) + t528 * t504 + qJDD(5) - t478;
t517 = -t527 * mrSges(6,2) + t535 * mrSges(6,3);
t574 = -m(6) * t476 + t513 * mrSges(6,1) + t535 * t517;
t472 = m(5) * t478 + t513 * mrSges(5,1) + t592 * t490 + t535 * t518 + t583 * t528 + t574;
t575 = t594 * t470 - t561 * t472;
t462 = m(4) * t486 - t558 * mrSges(4,2) + t514 * mrSges(4,3) - t539 * t525 - t559 * t532 + t575;
t485 = t565 * t507 - t562 * t508;
t531 = -t559 * mrSges(4,2) - t539 * mrSges(4,3);
t483 = -t558 * pkin(3) - t557 * pkin(8) + t540 * t526 - t485;
t477 = -0.2e1 * qJD(5) * t528 + (t527 * t535 - t490) * qJ(5) + (t528 * t535 + t489) * pkin(4) + t483;
t473 = m(6) * t477 + t489 * mrSges(6,1) - t490 * mrSges(6,3) + t527 * t517 - t528 * t520;
t570 = -m(5) * t483 - t489 * mrSges(5,1) - t490 * mrSges(5,2) - t527 * t518 - t528 * t519 - t473;
t467 = m(4) * t485 + t558 * mrSges(4,1) - t515 * mrSges(4,3) - t540 * t525 + t559 * t531 + t570;
t457 = t562 * t462 + t565 * t467;
t529 = -t566 * g(3) - t588;
t546 = (-mrSges(3,1) * t566 + mrSges(3,2) * t563) * qJD(1);
t581 = qJD(1) * t566;
t550 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t581;
t455 = m(3) * t529 + qJDD(2) * mrSges(3,1) - t547 * mrSges(3,3) + qJD(2) * t550 - t546 * t582 + t457;
t549 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t582;
t576 = t565 * t462 - t562 * t467;
t456 = m(3) * t530 - qJDD(2) * mrSges(3,2) + t548 * mrSges(3,3) - qJD(2) * t549 + t546 * t581 + t576;
t577 = -t563 * t455 + t566 * t456;
t449 = m(2) * t553 - t568 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t577;
t541 = -t568 * pkin(6) + t572;
t465 = t561 * t470 + t594 * t472;
t571 = m(4) * t516 - t514 * mrSges(4,1) + t515 * mrSges(4,2) + t539 * t531 + t540 * t532 + t465;
t569 = -m(3) * t541 + t548 * mrSges(3,1) - t547 * mrSges(3,2) - t549 * t582 + t550 * t581 - t571;
t459 = m(2) * t552 + qJDD(1) * mrSges(2,1) - t568 * mrSges(2,2) + t569;
t587 = t564 * t449 + t567 * t459;
t450 = t566 * t455 + t563 * t456;
t586 = t596 * t527 - t591 * t528 - t589 * t535;
t585 = t589 * t527 + t597 * t528 + t595 * t535;
t584 = -t591 * t527 + t598 * t528 - t597 * t535;
t578 = t567 * t449 - t564 * t459;
t538 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t563 + Ifges(3,4) * t566) * qJD(1);
t537 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t563 + Ifges(3,2) * t566) * qJD(1);
t536 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t563 + Ifges(3,6) * t566) * qJD(1);
t523 = Ifges(4,1) * t540 - Ifges(4,4) * t539 + Ifges(4,5) * t559;
t522 = Ifges(4,4) * t540 - Ifges(4,2) * t539 + Ifges(4,6) * t559;
t521 = Ifges(4,5) * t540 - Ifges(4,6) * t539 + Ifges(4,3) * t559;
t464 = mrSges(5,2) * t483 + mrSges(6,2) * t476 - mrSges(5,3) * t478 - mrSges(6,3) * t477 - qJ(5) * t473 - t591 * t489 + t598 * t490 - t513 * t597 + t585 * t527 + t586 * t535;
t463 = -mrSges(5,1) * t483 - mrSges(6,1) * t477 + mrSges(6,2) * t475 + mrSges(5,3) * t479 - pkin(4) * t473 - t596 * t489 + t591 * t490 + t589 * t513 + t585 * t528 + t584 * t535;
t451 = Ifges(4,4) * t515 + Ifges(4,2) * t514 + Ifges(4,6) * t558 - t540 * t521 + t559 * t523 - mrSges(4,1) * t516 + mrSges(4,3) * t486 - mrSges(5,1) * t478 + mrSges(5,2) * t479 + mrSges(6,1) * t476 - mrSges(6,3) * t475 - pkin(4) * t574 - qJ(5) * t579 - pkin(3) * t465 + (pkin(4) * t505 + t586) * t528 + (qJ(5) * t505 - t584) * t527 + t595 * t513 + (pkin(4) * mrSges(6,2) + t597) * t490 + (qJ(5) * mrSges(6,2) + t589) * t489;
t446 = mrSges(4,2) * t516 - mrSges(4,3) * t485 + Ifges(4,1) * t515 + Ifges(4,4) * t514 + Ifges(4,5) * t558 - pkin(8) * t465 - t561 * t463 + t594 * t464 - t539 * t521 - t559 * t522;
t445 = mrSges(3,2) * t541 - mrSges(3,3) * t529 + Ifges(3,1) * t547 + Ifges(3,4) * t548 + Ifges(3,5) * qJDD(2) - pkin(7) * t457 - qJD(2) * t537 + t565 * t446 - t562 * t451 + t536 * t581;
t444 = -mrSges(3,1) * t541 + mrSges(3,3) * t530 + Ifges(3,4) * t547 + Ifges(3,2) * t548 + Ifges(3,6) * qJDD(2) - pkin(2) * t571 + pkin(7) * t576 + qJD(2) * t538 + t562 * t446 + t565 * t451 - t536 * t582;
t443 = -pkin(1) * t450 + mrSges(2,3) * t553 - pkin(2) * t457 - Ifges(3,5) * t547 - Ifges(3,6) * t548 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t529 + mrSges(3,2) * t530 - pkin(8) * t575 - Ifges(4,5) * t515 - Ifges(4,6) * t514 - Ifges(4,3) * t558 - mrSges(4,1) * t485 + mrSges(4,2) * t486 - t561 * t464 - t594 * t463 - pkin(3) * t570 + mrSges(2,1) * g(3) + t568 * Ifges(2,5) - t540 * t522 - t539 * t523 + Ifges(2,6) * qJDD(1) + (-t563 * t537 + t566 * t538) * qJD(1);
t442 = -mrSges(2,2) * g(3) - mrSges(2,3) * t552 + Ifges(2,5) * qJDD(1) - t568 * Ifges(2,6) - pkin(6) * t450 - t563 * t444 + t566 * t445;
t1 = [-m(1) * g(1) + t578; -m(1) * g(2) + t587; (-m(1) - m(2)) * g(3) + t450; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t587 + t567 * t442 - t564 * t443; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t578 + t564 * t442 + t567 * t443; -mrSges(1,1) * g(2) + mrSges(2,1) * t552 + mrSges(1,2) * g(1) - mrSges(2,2) * t553 + Ifges(2,3) * qJDD(1) + pkin(1) * t569 + pkin(6) * t577 + t566 * t444 + t563 * t445;];
tauB = t1;
