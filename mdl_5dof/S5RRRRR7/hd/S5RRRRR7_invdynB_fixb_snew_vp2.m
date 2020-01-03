% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:54
% EndTime: 2019-12-31 22:22:01
% DurationCPUTime: 6.76s
% Computational Cost: add. (101233->314), mult. (219704->399), div. (0->0), fcn. (158228->10), ass. (0->123)
t582 = qJD(1) ^ 2;
t598 = pkin(2) * t582;
t576 = sin(qJ(1));
t581 = cos(qJ(1));
t564 = -t581 * g(1) - t576 * g(2);
t553 = -t582 * pkin(1) + qJDD(1) * pkin(6) + t564;
t575 = sin(qJ(2));
t597 = t575 * t553;
t580 = cos(qJ(2));
t593 = qJD(1) * qJD(2);
t558 = t575 * qJDD(1) + t580 * t593;
t521 = qJDD(2) * pkin(2) - t558 * pkin(7) - t597 + (pkin(7) * t593 + t575 * t598 - g(3)) * t580;
t541 = -t575 * g(3) + t580 * t553;
t559 = t580 * qJDD(1) - t575 * t593;
t595 = qJD(1) * t575;
t562 = qJD(2) * pkin(2) - pkin(7) * t595;
t571 = t580 ^ 2;
t522 = t559 * pkin(7) - qJD(2) * t562 - t571 * t598 + t541;
t574 = sin(qJ(3));
t579 = cos(qJ(3));
t506 = t579 * t521 - t574 * t522;
t550 = (-t574 * t575 + t579 * t580) * qJD(1);
t527 = t550 * qJD(3) + t579 * t558 + t574 * t559;
t551 = (t574 * t580 + t575 * t579) * qJD(1);
t569 = qJDD(2) + qJDD(3);
t570 = qJD(2) + qJD(3);
t489 = (t550 * t570 - t527) * pkin(8) + (t550 * t551 + t569) * pkin(3) + t506;
t507 = t574 * t521 + t579 * t522;
t526 = -t551 * qJD(3) - t574 * t558 + t579 * t559;
t544 = t570 * pkin(3) - t551 * pkin(8);
t546 = t550 ^ 2;
t492 = -t546 * pkin(3) + t526 * pkin(8) - t570 * t544 + t507;
t573 = sin(qJ(4));
t578 = cos(qJ(4));
t487 = t573 * t489 + t578 * t492;
t538 = t573 * t550 + t578 * t551;
t504 = -t538 * qJD(4) + t578 * t526 - t573 * t527;
t537 = t578 * t550 - t573 * t551;
t516 = -t537 * mrSges(5,1) + t538 * mrSges(5,2);
t567 = qJD(4) + t570;
t530 = t567 * mrSges(5,1) - t538 * mrSges(5,3);
t566 = qJDD(4) + t569;
t517 = -t537 * pkin(4) - t538 * pkin(9);
t565 = t567 ^ 2;
t484 = -t565 * pkin(4) + t566 * pkin(9) + t537 * t517 + t487;
t563 = t576 * g(1) - t581 * g(2);
t587 = -qJDD(1) * pkin(1) - t563;
t528 = -t559 * pkin(2) + t562 * t595 + (-pkin(7) * t571 - pkin(6)) * t582 + t587;
t496 = -t526 * pkin(3) - t546 * pkin(8) + t551 * t544 + t528;
t505 = t537 * qJD(4) + t573 * t526 + t578 * t527;
t485 = (-t537 * t567 - t505) * pkin(9) + (t538 * t567 - t504) * pkin(4) + t496;
t572 = sin(qJ(5));
t577 = cos(qJ(5));
t481 = -t572 * t484 + t577 * t485;
t523 = -t572 * t538 + t577 * t567;
t494 = t523 * qJD(5) + t577 * t505 + t572 * t566;
t502 = qJDD(5) - t504;
t524 = t577 * t538 + t572 * t567;
t508 = -t523 * mrSges(6,1) + t524 * mrSges(6,2);
t534 = qJD(5) - t537;
t509 = -t534 * mrSges(6,2) + t523 * mrSges(6,3);
t479 = m(6) * t481 + t502 * mrSges(6,1) - t494 * mrSges(6,3) - t524 * t508 + t534 * t509;
t482 = t577 * t484 + t572 * t485;
t493 = -t524 * qJD(5) - t572 * t505 + t577 * t566;
t510 = t534 * mrSges(6,1) - t524 * mrSges(6,3);
t480 = m(6) * t482 - t502 * mrSges(6,2) + t493 * mrSges(6,3) + t523 * t508 - t534 * t510;
t588 = -t572 * t479 + t577 * t480;
t470 = m(5) * t487 - t566 * mrSges(5,2) + t504 * mrSges(5,3) + t537 * t516 - t567 * t530 + t588;
t486 = t578 * t489 - t573 * t492;
t529 = -t567 * mrSges(5,2) + t537 * mrSges(5,3);
t483 = -t566 * pkin(4) - t565 * pkin(9) + t538 * t517 - t486;
t585 = -m(6) * t483 + t493 * mrSges(6,1) - t494 * mrSges(6,2) + t523 * t509 - t524 * t510;
t475 = m(5) * t486 + t566 * mrSges(5,1) - t505 * mrSges(5,3) - t538 * t516 + t567 * t529 + t585;
t465 = t573 * t470 + t578 * t475;
t539 = -t550 * mrSges(4,1) + t551 * mrSges(4,2);
t542 = -t570 * mrSges(4,2) + t550 * mrSges(4,3);
t463 = m(4) * t506 + t569 * mrSges(4,1) - t527 * mrSges(4,3) - t551 * t539 + t570 * t542 + t465;
t543 = t570 * mrSges(4,1) - t551 * mrSges(4,3);
t589 = t578 * t470 - t573 * t475;
t464 = m(4) * t507 - t569 * mrSges(4,2) + t526 * mrSges(4,3) + t550 * t539 - t570 * t543 + t589;
t457 = t579 * t463 + t574 * t464;
t540 = -t580 * g(3) - t597;
t557 = (-mrSges(3,1) * t580 + mrSges(3,2) * t575) * qJD(1);
t594 = qJD(1) * t580;
t561 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t594;
t455 = m(3) * t540 + qJDD(2) * mrSges(3,1) - t558 * mrSges(3,3) + qJD(2) * t561 - t557 * t595 + t457;
t560 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t595;
t590 = -t574 * t463 + t579 * t464;
t456 = m(3) * t541 - qJDD(2) * mrSges(3,2) + t559 * mrSges(3,3) - qJD(2) * t560 + t557 * t594 + t590;
t591 = -t575 * t455 + t580 * t456;
t449 = m(2) * t564 - t582 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t591;
t552 = -t582 * pkin(6) + t587;
t471 = t577 * t479 + t572 * t480;
t586 = m(5) * t496 - t504 * mrSges(5,1) + t505 * mrSges(5,2) - t537 * t529 + t538 * t530 + t471;
t584 = m(4) * t528 - t526 * mrSges(4,1) + t527 * mrSges(4,2) - t550 * t542 + t551 * t543 + t586;
t583 = -m(3) * t552 + t559 * mrSges(3,1) - t558 * mrSges(3,2) - t560 * t595 + t561 * t594 - t584;
t467 = m(2) * t563 + qJDD(1) * mrSges(2,1) - t582 * mrSges(2,2) + t583;
t596 = t576 * t449 + t581 * t467;
t450 = t580 * t455 + t575 * t456;
t592 = t581 * t449 - t576 * t467;
t549 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t575 + Ifges(3,4) * t580) * qJD(1);
t548 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t575 + Ifges(3,2) * t580) * qJD(1);
t547 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t575 + Ifges(3,6) * t580) * qJD(1);
t533 = Ifges(4,1) * t551 + Ifges(4,4) * t550 + Ifges(4,5) * t570;
t532 = Ifges(4,4) * t551 + Ifges(4,2) * t550 + Ifges(4,6) * t570;
t531 = Ifges(4,5) * t551 + Ifges(4,6) * t550 + Ifges(4,3) * t570;
t513 = Ifges(5,1) * t538 + Ifges(5,4) * t537 + Ifges(5,5) * t567;
t512 = Ifges(5,4) * t538 + Ifges(5,2) * t537 + Ifges(5,6) * t567;
t511 = Ifges(5,5) * t538 + Ifges(5,6) * t537 + Ifges(5,3) * t567;
t499 = Ifges(6,1) * t524 + Ifges(6,4) * t523 + Ifges(6,5) * t534;
t498 = Ifges(6,4) * t524 + Ifges(6,2) * t523 + Ifges(6,6) * t534;
t497 = Ifges(6,5) * t524 + Ifges(6,6) * t523 + Ifges(6,3) * t534;
t473 = mrSges(6,2) * t483 - mrSges(6,3) * t481 + Ifges(6,1) * t494 + Ifges(6,4) * t493 + Ifges(6,5) * t502 + t523 * t497 - t534 * t498;
t472 = -mrSges(6,1) * t483 + mrSges(6,3) * t482 + Ifges(6,4) * t494 + Ifges(6,2) * t493 + Ifges(6,6) * t502 - t524 * t497 + t534 * t499;
t459 = -mrSges(5,1) * t496 - mrSges(6,1) * t481 + mrSges(6,2) * t482 + mrSges(5,3) * t487 + Ifges(5,4) * t505 - Ifges(6,5) * t494 + Ifges(5,2) * t504 + Ifges(5,6) * t566 - Ifges(6,6) * t493 - Ifges(6,3) * t502 - pkin(4) * t471 - t524 * t498 + t523 * t499 - t538 * t511 + t567 * t513;
t458 = mrSges(5,2) * t496 - mrSges(5,3) * t486 + Ifges(5,1) * t505 + Ifges(5,4) * t504 + Ifges(5,5) * t566 - pkin(9) * t471 - t572 * t472 + t577 * t473 + t537 * t511 - t567 * t512;
t451 = mrSges(4,2) * t528 - mrSges(4,3) * t506 + Ifges(4,1) * t527 + Ifges(4,4) * t526 + Ifges(4,5) * t569 - pkin(8) * t465 + t578 * t458 - t573 * t459 + t550 * t531 - t570 * t532;
t446 = -mrSges(4,1) * t528 + mrSges(4,3) * t507 + Ifges(4,4) * t527 + Ifges(4,2) * t526 + Ifges(4,6) * t569 - pkin(3) * t586 + pkin(8) * t589 + t573 * t458 + t578 * t459 - t551 * t531 + t570 * t533;
t445 = Ifges(2,6) * qJDD(1) - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) + (-t575 * t548 + t580 * t549) * qJD(1) - t577 * t472 + t582 * Ifges(2,5) - t572 * t473 - Ifges(5,3) * t566 - Ifges(4,3) * t569 - t551 * t532 - Ifges(3,5) * t558 - Ifges(3,6) * t559 + mrSges(2,3) * t564 - mrSges(3,1) * t540 + mrSges(3,2) * t541 + t550 * t533 + t537 * t513 - t538 * t512 - Ifges(4,6) * t526 - Ifges(4,5) * t527 - Ifges(5,6) * t504 - Ifges(5,5) * t505 - mrSges(4,1) * t506 + mrSges(4,2) * t507 + mrSges(5,2) * t487 - mrSges(5,1) * t486 - pkin(3) * t465 - pkin(9) * t588 - pkin(4) * t585 - pkin(2) * t457 - pkin(1) * t450;
t444 = mrSges(3,2) * t552 - mrSges(3,3) * t540 + Ifges(3,1) * t558 + Ifges(3,4) * t559 + Ifges(3,5) * qJDD(2) - pkin(7) * t457 - qJD(2) * t548 - t574 * t446 + t579 * t451 + t547 * t594;
t443 = -mrSges(3,1) * t552 + mrSges(3,3) * t541 + Ifges(3,4) * t558 + Ifges(3,2) * t559 + Ifges(3,6) * qJDD(2) - pkin(2) * t584 + pkin(7) * t590 + qJD(2) * t549 + t579 * t446 + t574 * t451 - t547 * t595;
t442 = -mrSges(2,2) * g(3) - mrSges(2,3) * t563 + Ifges(2,5) * qJDD(1) - t582 * Ifges(2,6) - pkin(6) * t450 - t575 * t443 + t580 * t444;
t1 = [-m(1) * g(1) + t592; -m(1) * g(2) + t596; (-m(1) - m(2)) * g(3) + t450; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t596 + t581 * t442 - t576 * t445; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t592 + t576 * t442 + t581 * t445; -mrSges(1,1) * g(2) + mrSges(2,1) * t563 + mrSges(1,2) * g(1) - mrSges(2,2) * t564 + Ifges(2,3) * qJDD(1) + pkin(1) * t583 + pkin(6) * t591 + t580 * t443 + t575 * t444;];
tauB = t1;
