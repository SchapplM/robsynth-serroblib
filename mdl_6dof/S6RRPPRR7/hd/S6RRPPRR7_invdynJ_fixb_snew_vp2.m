% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:03:17
% EndTime: 2019-05-06 11:03:22
% DurationCPUTime: 3.46s
% Computational Cost: add. (23421->330), mult. (53927->402), div. (0->0), fcn. (36652->10), ass. (0->147)
t654 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t631 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t630 = Ifges(5,5) + Ifges(3,6) - Ifges(4,6);
t629 = Ifges(5,6) + Ifges(3,5) + Ifges(4,4);
t653 = Ifges(4,3) + Ifges(5,1) + Ifges(3,2);
t652 = -Ifges(5,3) - Ifges(3,3) - Ifges(4,2);
t599 = qJD(1) ^ 2;
t594 = sin(qJ(1));
t598 = cos(qJ(1));
t620 = t594 * g(1) - g(2) * t598;
t589 = sin(pkin(6));
t646 = pkin(8) * t589;
t557 = qJDD(1) * pkin(1) + t599 * t646 + t620;
t590 = cos(pkin(6));
t524 = -t590 * g(3) - t589 * t557;
t593 = sin(qJ(2));
t597 = cos(qJ(2));
t635 = qJD(1) * t597;
t564 = (qJD(2) * t635 + qJDD(1) * t593) * t589;
t636 = qJD(1) * t589;
t622 = t593 * t636;
t632 = qJDD(1) * t589;
t565 = -qJD(2) * t622 + t597 * t632;
t582 = qJD(1) * t590 + qJD(2);
t621 = t589 * t635;
t617 = t582 * t621;
t615 = -t565 * pkin(2) + t524 + (-t564 - t617) * qJ(3);
t647 = pkin(2) * t582;
t487 = (-(2 * qJD(3)) + t647) * t622 + t615;
t651 = m(4) * t487 - t565 * mrSges(4,1);
t650 = t582 ^ 2;
t649 = pkin(3) + pkin(9);
t581 = qJDD(1) * t590 + qJDD(2);
t648 = pkin(2) * t581;
t645 = t593 * g(3);
t644 = mrSges(3,3) + mrSges(4,2);
t643 = t557 * t590;
t642 = t589 ^ 2 * t599;
t641 = t589 * t593;
t640 = t589 * t597;
t550 = -pkin(3) * t582 - qJ(4) * t622;
t627 = t597 ^ 2 * t642;
t605 = -qJ(4) * t627 + qJDD(4) - t615 + ((2 * qJD(3)) + t550) * t622;
t473 = pkin(4) * t564 + t649 * t565 + (pkin(4) * t597 + (-pkin(2) - pkin(9)) * t593) * t582 * t636 + t605;
t563 = (t593 * pkin(4) + t597 * pkin(9)) * t636;
t616 = -g(1) * t598 - g(2) * t594;
t558 = -pkin(1) * t599 + pkin(8) * t632 + t616;
t507 = -g(3) * t640 - t593 * t558 + t597 * t643;
t559 = (-t597 * pkin(2) - t593 * qJ(3)) * t636;
t609 = -qJ(3) * t650 + t559 * t622 + qJDD(3) - t507;
t633 = qJD(1) * qJD(4);
t618 = -0.2e1 * t589 * t633;
t602 = t593 * t618 + t609 + (-t564 + t617) * qJ(4);
t626 = t597 * t642;
t477 = -pkin(4) * t650 + (-pkin(3) * t626 - t563 * t636) * t593 + (-pkin(2) - t649) * t581 + t602;
t592 = sin(qJ(5));
t596 = cos(qJ(5));
t470 = t592 * t473 + t596 * t477;
t539 = -t582 * t592 - t596 * t621;
t505 = -qJD(5) * t539 + t565 * t592 - t581 * t596;
t538 = -t582 * t596 + t592 * t621;
t509 = -mrSges(6,1) * t538 + mrSges(6,2) * t539;
t570 = qJD(5) + t622;
t514 = mrSges(6,1) * t570 - mrSges(6,3) * t539;
t548 = qJDD(5) + t564;
t510 = -pkin(5) * t538 - pkin(10) * t539;
t568 = t570 ^ 2;
t467 = -pkin(5) * t568 + pkin(10) * t548 + t510 * t538 + t470;
t634 = qJD(3) * t582;
t569 = 0.2e1 * t634;
t638 = t597 * t558 + t593 * t643;
t614 = pkin(2) * t650 - t581 * qJ(3) - t559 * t621 - t638;
t606 = pkin(3) * t627 + t565 * qJ(4) - t582 * t550 + t614;
t475 = t581 * pkin(4) - t650 * pkin(9) + t569 + t597 * t618 + (-t563 * t635 - t645) * t589 - t606;
t506 = qJD(5) * t538 - t565 * t596 - t581 * t592;
t471 = (t539 * t570 - t505) * pkin(5) + (-t538 * t570 - t506) * pkin(10) + t475;
t591 = sin(qJ(6));
t595 = cos(qJ(6));
t464 = -t467 * t591 + t471 * t595;
t511 = -t539 * t591 + t570 * t595;
t483 = qJD(6) * t511 + t506 * t595 + t548 * t591;
t512 = t539 * t595 + t570 * t591;
t493 = -mrSges(7,1) * t511 + mrSges(7,2) * t512;
t536 = qJD(6) - t538;
t495 = -mrSges(7,2) * t536 + mrSges(7,3) * t511;
t503 = qJDD(6) - t505;
t461 = m(7) * t464 + mrSges(7,1) * t503 - mrSges(7,3) * t483 - t493 * t512 + t495 * t536;
t465 = t467 * t595 + t471 * t591;
t482 = -qJD(6) * t512 - t506 * t591 + t548 * t595;
t496 = mrSges(7,1) * t536 - mrSges(7,3) * t512;
t462 = m(7) * t465 - mrSges(7,2) * t503 + mrSges(7,3) * t482 + t493 * t511 - t496 * t536;
t619 = -t461 * t591 + t595 * t462;
t449 = m(6) * t470 - mrSges(6,2) * t548 + mrSges(6,3) * t505 + t509 * t538 - t514 * t570 + t619;
t469 = t473 * t596 - t477 * t592;
t513 = -mrSges(6,2) * t570 + mrSges(6,3) * t538;
t466 = -pkin(5) * t548 - pkin(10) * t568 + t510 * t539 - t469;
t612 = -m(7) * t466 + t482 * mrSges(7,1) - mrSges(7,2) * t483 + t511 * t495 - t496 * t512;
t457 = m(6) * t469 + mrSges(6,1) * t548 - mrSges(6,3) * t506 - t509 * t539 + t513 * t570 + t612;
t639 = t596 * t449 - t592 * t457;
t451 = t595 * t461 + t591 * t462;
t551 = mrSges(5,2) * t582 - mrSges(5,3) * t622;
t637 = -mrSges(3,1) * t582 + mrSges(3,3) * t622 + t551;
t628 = g(3) * t641;
t625 = (-t593 * t629 - t597 * t630) * t636 + t652 * t582;
t624 = (-t593 * t631 - t597 * t653) * t636 - t630 * t582;
t623 = (-t593 * t654 - t631 * t597) * t636 - t629 * t582;
t446 = t592 * t449 + t596 * t457;
t479 = -t648 + (-t593 * t626 - t581) * pkin(3) + t602;
t554 = -mrSges(5,1) * t582 + mrSges(5,3) * t621;
t562 = (t593 * mrSges(5,1) - t597 * mrSges(5,2)) * t636;
t613 = m(5) * t479 + t581 * mrSges(5,2) + t582 * t554 - t562 * t622 + t639;
t480 = pkin(3) * t565 - t622 * t647 + t605;
t611 = m(5) * t480 + t564 * mrSges(5,1) - t554 * t621 + t446;
t610 = -m(6) * t475 + t505 * mrSges(6,1) - t506 * mrSges(6,2) + t538 * t513 - t539 * t514 - t451;
t492 = t609 - t648;
t556 = mrSges(4,2) * t621 + mrSges(4,3) * t582;
t608 = m(4) * t492 - t581 * mrSges(4,1) - t582 * t556 + t613;
t607 = t565 * mrSges(5,2) - t611;
t478 = -0.2e1 * t634 + (0.2e1 * t597 * t633 + t645) * t589 + t606;
t604 = -m(5) * t478 - t565 * mrSges(5,3) - t610;
t488 = Ifges(7,5) * t512 + Ifges(7,6) * t511 + Ifges(7,3) * t536;
t490 = Ifges(7,1) * t512 + Ifges(7,4) * t511 + Ifges(7,5) * t536;
t454 = -mrSges(7,1) * t466 + mrSges(7,3) * t465 + Ifges(7,4) * t483 + Ifges(7,2) * t482 + Ifges(7,6) * t503 - t488 * t512 + t490 * t536;
t489 = Ifges(7,4) * t512 + Ifges(7,2) * t511 + Ifges(7,6) * t536;
t455 = mrSges(7,2) * t466 - mrSges(7,3) * t464 + Ifges(7,1) * t483 + Ifges(7,4) * t482 + Ifges(7,5) * t503 + t488 * t511 - t489 * t536;
t498 = Ifges(6,4) * t539 + Ifges(6,2) * t538 + Ifges(6,6) * t570;
t499 = Ifges(6,1) * t539 + Ifges(6,4) * t538 + Ifges(6,5) * t570;
t603 = mrSges(6,1) * t469 - mrSges(6,2) * t470 + Ifges(6,5) * t506 + Ifges(6,6) * t505 + Ifges(6,3) * t548 + pkin(5) * t612 + pkin(10) * t619 + t595 * t454 + t591 * t455 + t539 * t498 - t538 * t499;
t601 = mrSges(7,1) * t464 - mrSges(7,2) * t465 + Ifges(7,5) * t483 + Ifges(7,6) * t482 + Ifges(7,3) * t503 + t489 * t512 - t490 * t511;
t486 = t569 - t614 - t628;
t553 = -mrSges(4,1) * t582 + mrSges(4,2) * t622;
t560 = (-t597 * mrSges(4,1) - t593 * mrSges(4,3)) * t636;
t600 = m(4) * t486 + t581 * mrSges(4,3) + t582 * t553 + t560 * t621 + t604;
t561 = (-mrSges(3,1) * t597 + mrSges(3,2) * t593) * t636;
t555 = -mrSges(3,2) * t582 + mrSges(3,3) * t621;
t508 = -t628 + t638;
t497 = Ifges(6,5) * t539 + Ifges(6,6) * t538 + Ifges(6,3) * t570;
t447 = t637 * t582 + (-mrSges(3,2) + mrSges(5,1)) * t581 + t644 * t565 + (t561 - t562) * t621 + m(3) * t508 + t600;
t445 = -mrSges(5,3) * t564 + t613;
t444 = t551 * t622 - t607;
t443 = t560 * t622 + (mrSges(4,2) - mrSges(5,3)) * t564 + t608;
t442 = -t564 * mrSges(4,3) + (-t556 * t597 + (-t551 - t553) * t593) * t636 + t607 + t651;
t441 = m(3) * t507 + mrSges(3,1) * t581 + t555 * t582 + (-t560 - t561) * t622 + (mrSges(5,3) - t644) * t564 - t608;
t440 = -mrSges(6,1) * t475 + mrSges(6,3) * t470 + Ifges(6,4) * t506 + Ifges(6,2) * t505 + Ifges(6,6) * t548 - pkin(5) * t451 - t497 * t539 + t499 * t570 - t601;
t439 = mrSges(6,2) * t475 - mrSges(6,3) * t469 + Ifges(6,1) * t506 + Ifges(6,4) * t505 + Ifges(6,5) * t548 - pkin(10) * t451 - t454 * t591 + t455 * t595 + t497 * t538 - t498 * t570;
t438 = -t596 * t440 + qJ(3) * (t582 * t551 + t600) - t592 * t439 - pkin(4) * t610 + mrSges(3,1) * t507 - mrSges(3,2) * t508 + mrSges(4,3) * t486 - mrSges(4,1) * t492 - mrSges(5,1) * t478 + mrSges(5,2) * t479 - pkin(9) * t639 - pkin(3) * t445 - pkin(2) * t443 + (qJ(3) * mrSges(5,1) - t652) * t581 + (qJ(3) * mrSges(4,2) + t630) * t565 + t629 * t564 + (-t624 * t593 + (-qJ(3) * t562 + t623) * t597) * t636;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t620 - mrSges(2,2) * t616 + (mrSges(5,1) * t480 + mrSges(3,2) * t524 + mrSges(4,2) * t492 - mrSges(3,3) * t507 - mrSges(4,3) * t487 - mrSges(5,3) * t479 + pkin(4) * t446 - qJ(3) * t442 - qJ(4) * t445 + t564 * t654 + t631 * t565 + t629 * t581 + t624 * t582 - t625 * t621 + t603) * t641 + (-t596 * t439 - qJ(4) * t604 + t592 * t440 - mrSges(3,1) * t524 + mrSges(3,3) * t508 + mrSges(4,2) * t486 - mrSges(4,1) * t487 - mrSges(5,2) * t480 + mrSges(5,3) * t478 + pkin(9) * t446 + pkin(3) * t444 - pkin(2) * t442 + (-qJ(4) * t551 - t623) * t582 + (-qJ(4) * mrSges(5,1) + t630) * t581 + t653 * t565 + t631 * t564 + (qJ(4) * t562 * t597 + t625 * t593) * t636) * t640 + t590 * t438 + pkin(1) * ((t441 * t597 + t447 * t593) * t590 + (-m(3) * t524 + (mrSges(3,1) - mrSges(5,2)) * t565 + (-mrSges(3,2) + mrSges(4,3)) * t564 + ((t555 + t556) * t597 + (t553 + t637) * t593) * t636 + t611 - t651) * t589) + (-t441 * t593 + t447 * t597) * t646; t438; t443; t444; t603; t601;];
tauJ  = t1;
