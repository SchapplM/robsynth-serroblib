% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:57:45
% EndTime: 2019-05-06 02:58:01
% DurationCPUTime: 15.48s
% Computational Cost: add. (271352->343), mult. (534921->428), div. (0->0), fcn. (368397->12), ass. (0->136)
t645 = sin(qJ(1));
t650 = cos(qJ(1));
t630 = t645 * g(1) - g(2) * t650;
t621 = qJDD(1) * pkin(1) + t630;
t631 = -g(1) * t650 - g(2) * t645;
t652 = qJD(1) ^ 2;
t623 = -pkin(1) * t652 + t631;
t639 = sin(pkin(11));
t640 = cos(pkin(11));
t601 = t639 * t621 + t640 * t623;
t588 = -pkin(2) * t652 + qJDD(1) * pkin(7) + t601;
t638 = -g(3) + qJDD(2);
t644 = sin(qJ(3));
t649 = cos(qJ(3));
t581 = t649 * t588 + t644 * t638;
t622 = (-mrSges(4,1) * t649 + mrSges(4,2) * t644) * qJD(1);
t665 = qJD(1) * qJD(3);
t635 = t644 * t665;
t626 = qJDD(1) * t649 - t635;
t667 = qJD(1) * t644;
t628 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t667;
t600 = t640 * t621 - t639 * t623;
t587 = -qJDD(1) * pkin(2) - t652 * pkin(7) - t600;
t663 = t649 * t665;
t625 = qJDD(1) * t644 + t663;
t572 = (-t625 - t663) * pkin(8) + (-t626 + t635) * pkin(3) + t587;
t624 = (-pkin(3) * t649 - pkin(8) * t644) * qJD(1);
t651 = qJD(3) ^ 2;
t666 = qJD(1) * t649;
t578 = -pkin(3) * t651 + qJDD(3) * pkin(8) + t624 * t666 + t581;
t643 = sin(qJ(4));
t648 = cos(qJ(4));
t559 = t648 * t572 - t643 * t578;
t619 = qJD(3) * t648 - t643 * t667;
t596 = qJD(4) * t619 + qJDD(3) * t643 + t625 * t648;
t618 = qJDD(4) - t626;
t620 = qJD(3) * t643 + t648 * t667;
t633 = qJD(4) - t666;
t550 = (t619 * t633 - t596) * pkin(9) + (t619 * t620 + t618) * pkin(4) + t559;
t560 = t643 * t572 + t648 * t578;
t595 = -qJD(4) * t620 + qJDD(3) * t648 - t625 * t643;
t605 = pkin(4) * t633 - pkin(9) * t620;
t617 = t619 ^ 2;
t552 = -pkin(4) * t617 + pkin(9) * t595 - t605 * t633 + t560;
t642 = sin(qJ(5));
t647 = cos(qJ(5));
t540 = t647 * t550 - t642 * t552;
t598 = t619 * t647 - t620 * t642;
t564 = qJD(5) * t598 + t595 * t642 + t596 * t647;
t599 = t619 * t642 + t620 * t647;
t614 = qJDD(5) + t618;
t632 = qJD(5) + t633;
t538 = (t598 * t632 - t564) * pkin(10) + (t598 * t599 + t614) * pkin(5) + t540;
t541 = t642 * t550 + t647 * t552;
t563 = -qJD(5) * t599 + t595 * t647 - t596 * t642;
t584 = pkin(5) * t632 - pkin(10) * t599;
t597 = t598 ^ 2;
t539 = -pkin(5) * t597 + pkin(10) * t563 - t584 * t632 + t541;
t641 = sin(qJ(6));
t646 = cos(qJ(6));
t536 = t538 * t646 - t539 * t641;
t575 = t598 * t646 - t599 * t641;
t547 = qJD(6) * t575 + t563 * t641 + t564 * t646;
t576 = t598 * t641 + t599 * t646;
t561 = -mrSges(7,1) * t575 + mrSges(7,2) * t576;
t627 = qJD(6) + t632;
t565 = -mrSges(7,2) * t627 + mrSges(7,3) * t575;
t607 = qJDD(6) + t614;
t534 = m(7) * t536 + mrSges(7,1) * t607 - mrSges(7,3) * t547 - t561 * t576 + t565 * t627;
t537 = t538 * t641 + t539 * t646;
t546 = -qJD(6) * t576 + t563 * t646 - t564 * t641;
t566 = mrSges(7,1) * t627 - mrSges(7,3) * t576;
t535 = m(7) * t537 - mrSges(7,2) * t607 + mrSges(7,3) * t546 + t561 * t575 - t566 * t627;
t526 = t646 * t534 + t641 * t535;
t579 = -mrSges(6,1) * t598 + mrSges(6,2) * t599;
t582 = -mrSges(6,2) * t632 + mrSges(6,3) * t598;
t524 = m(6) * t540 + mrSges(6,1) * t614 - mrSges(6,3) * t564 - t579 * t599 + t582 * t632 + t526;
t583 = mrSges(6,1) * t632 - mrSges(6,3) * t599;
t657 = -t534 * t641 + t646 * t535;
t525 = m(6) * t541 - mrSges(6,2) * t614 + mrSges(6,3) * t563 + t579 * t598 - t583 * t632 + t657;
t520 = t647 * t524 + t642 * t525;
t602 = -mrSges(5,1) * t619 + mrSges(5,2) * t620;
t603 = -mrSges(5,2) * t633 + mrSges(5,3) * t619;
t518 = m(5) * t559 + mrSges(5,1) * t618 - mrSges(5,3) * t596 - t602 * t620 + t603 * t633 + t520;
t604 = mrSges(5,1) * t633 - mrSges(5,3) * t620;
t658 = -t524 * t642 + t647 * t525;
t519 = m(5) * t560 - mrSges(5,2) * t618 + mrSges(5,3) * t595 + t602 * t619 - t604 * t633 + t658;
t659 = -t518 * t643 + t648 * t519;
t513 = m(4) * t581 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t626 - qJD(3) * t628 + t622 * t666 + t659;
t580 = -t644 * t588 + t638 * t649;
t629 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t666;
t577 = -qJDD(3) * pkin(3) - pkin(8) * t651 + t624 * t667 - t580;
t558 = -pkin(4) * t595 - pkin(9) * t617 + t620 * t605 + t577;
t543 = -pkin(5) * t563 - pkin(10) * t597 + t584 * t599 + t558;
t656 = m(7) * t543 - t546 * mrSges(7,1) + t547 * mrSges(7,2) - t575 * t565 + t576 * t566;
t655 = m(6) * t558 - t563 * mrSges(6,1) + t564 * mrSges(6,2) - t598 * t582 + t599 * t583 + t656;
t653 = -m(5) * t577 + t595 * mrSges(5,1) - t596 * mrSges(5,2) + t619 * t603 - t620 * t604 - t655;
t530 = m(4) * t580 + qJDD(3) * mrSges(4,1) - t625 * mrSges(4,3) + qJD(3) * t629 - t622 * t667 + t653;
t660 = t649 * t513 - t530 * t644;
t507 = m(3) * t601 - mrSges(3,1) * t652 - qJDD(1) * mrSges(3,2) + t660;
t514 = t518 * t648 + t519 * t643;
t654 = -m(4) * t587 + t626 * mrSges(4,1) - mrSges(4,2) * t625 - t628 * t667 + t629 * t666 - t514;
t510 = m(3) * t600 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t652 + t654;
t501 = t639 * t507 + t640 * t510;
t499 = m(2) * t630 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t652 + t501;
t661 = t640 * t507 - t510 * t639;
t500 = m(2) * t631 - mrSges(2,1) * t652 - qJDD(1) * mrSges(2,2) + t661;
t668 = t650 * t499 + t645 * t500;
t508 = t644 * t513 + t649 * t530;
t664 = m(3) * t638 + t508;
t662 = -t499 * t645 + t650 * t500;
t613 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t644 + Ifges(4,4) * t649) * qJD(1);
t612 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t644 + Ifges(4,2) * t649) * qJD(1);
t611 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t644 + Ifges(4,6) * t649) * qJD(1);
t591 = Ifges(5,1) * t620 + Ifges(5,4) * t619 + Ifges(5,5) * t633;
t590 = Ifges(5,4) * t620 + Ifges(5,2) * t619 + Ifges(5,6) * t633;
t589 = Ifges(5,5) * t620 + Ifges(5,6) * t619 + Ifges(5,3) * t633;
t571 = Ifges(6,1) * t599 + Ifges(6,4) * t598 + Ifges(6,5) * t632;
t570 = Ifges(6,4) * t599 + Ifges(6,2) * t598 + Ifges(6,6) * t632;
t569 = Ifges(6,5) * t599 + Ifges(6,6) * t598 + Ifges(6,3) * t632;
t555 = Ifges(7,1) * t576 + Ifges(7,4) * t575 + Ifges(7,5) * t627;
t554 = Ifges(7,4) * t576 + Ifges(7,2) * t575 + Ifges(7,6) * t627;
t553 = Ifges(7,5) * t576 + Ifges(7,6) * t575 + Ifges(7,3) * t627;
t528 = mrSges(7,2) * t543 - mrSges(7,3) * t536 + Ifges(7,1) * t547 + Ifges(7,4) * t546 + Ifges(7,5) * t607 + t553 * t575 - t554 * t627;
t527 = -mrSges(7,1) * t543 + mrSges(7,3) * t537 + Ifges(7,4) * t547 + Ifges(7,2) * t546 + Ifges(7,6) * t607 - t553 * t576 + t555 * t627;
t516 = mrSges(6,2) * t558 - mrSges(6,3) * t540 + Ifges(6,1) * t564 + Ifges(6,4) * t563 + Ifges(6,5) * t614 - pkin(10) * t526 - t527 * t641 + t528 * t646 + t569 * t598 - t570 * t632;
t515 = -mrSges(6,1) * t558 + mrSges(6,3) * t541 + Ifges(6,4) * t564 + Ifges(6,2) * t563 + Ifges(6,6) * t614 - pkin(5) * t656 + pkin(10) * t657 + t646 * t527 + t641 * t528 - t599 * t569 + t632 * t571;
t504 = mrSges(5,2) * t577 - mrSges(5,3) * t559 + Ifges(5,1) * t596 + Ifges(5,4) * t595 + Ifges(5,5) * t618 - pkin(9) * t520 - t515 * t642 + t516 * t647 + t589 * t619 - t590 * t633;
t503 = -mrSges(5,1) * t577 + mrSges(5,3) * t560 + Ifges(5,4) * t596 + Ifges(5,2) * t595 + Ifges(5,6) * t618 - pkin(4) * t655 + pkin(9) * t658 + t647 * t515 + t642 * t516 - t620 * t589 + t633 * t591;
t502 = -mrSges(6,1) * t540 - mrSges(7,1) * t536 + mrSges(7,2) * t537 + Ifges(4,6) * qJDD(3) - pkin(5) * t526 - mrSges(5,1) * t559 + mrSges(5,2) * t560 - mrSges(4,1) * t587 - Ifges(5,6) * t595 - Ifges(5,5) * t596 - pkin(4) * t520 - pkin(3) * t514 + t575 * t555 - t576 * t554 + mrSges(4,3) * t581 + t619 * t591 - t620 * t590 + Ifges(4,4) * t625 + Ifges(4,2) * t626 - Ifges(7,3) * t607 + qJD(3) * t613 - Ifges(6,3) * t614 - Ifges(5,3) * t618 - Ifges(6,6) * t563 - Ifges(6,5) * t564 + mrSges(6,2) * t541 - Ifges(7,6) * t546 - t611 * t667 + t598 * t571 - t599 * t570 - Ifges(7,5) * t547;
t495 = mrSges(4,2) * t587 - mrSges(4,3) * t580 + Ifges(4,1) * t625 + Ifges(4,4) * t626 + Ifges(4,5) * qJDD(3) - pkin(8) * t514 - qJD(3) * t612 - t503 * t643 + t504 * t648 + t611 * t666;
t494 = Ifges(3,6) * qJDD(1) + t652 * Ifges(3,5) - mrSges(3,1) * t638 + mrSges(3,3) * t601 - Ifges(4,5) * t625 - Ifges(4,6) * t626 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t580 + mrSges(4,2) * t581 - t643 * t504 - t648 * t503 - pkin(3) * t653 - pkin(8) * t659 - pkin(2) * t508 + (-t612 * t644 + t613 * t649) * qJD(1);
t493 = mrSges(3,2) * t638 - mrSges(3,3) * t600 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t652 - pkin(7) * t508 + t495 * t649 - t502 * t644;
t492 = -mrSges(2,2) * g(3) - mrSges(2,3) * t630 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t652 - qJ(2) * t501 + t493 * t640 - t494 * t639;
t491 = mrSges(2,1) * g(3) + mrSges(2,3) * t631 + t652 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t664 + qJ(2) * t661 + t639 * t493 + t640 * t494;
t1 = [-m(1) * g(1) + t662; -m(1) * g(2) + t668; (-m(1) - m(2)) * g(3) + t664; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t668 - t645 * t491 + t650 * t492; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t662 + t650 * t491 + t645 * t492; pkin(1) * t501 + mrSges(2,1) * t630 - mrSges(2,2) * t631 + pkin(7) * t660 + t644 * t495 + t649 * t502 + pkin(2) * t654 + mrSges(3,1) * t600 - mrSges(3,2) * t601 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
