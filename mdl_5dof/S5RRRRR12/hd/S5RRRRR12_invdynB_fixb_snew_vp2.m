% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR12_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:48:36
% EndTime: 2019-12-31 22:49:10
% DurationCPUTime: 33.12s
% Computational Cost: add. (524744->343), mult. (1303988->465), div. (0->0), fcn. (1077567->14), ass. (0->153)
t639 = cos(pkin(5));
t634 = t639 * qJD(1) + qJD(2);
t636 = sin(pkin(6));
t638 = cos(pkin(6));
t637 = sin(pkin(5));
t648 = cos(qJ(2));
t667 = qJD(1) * t648;
t663 = t637 * t667;
t619 = (t634 * t636 + t638 * t663) * pkin(9);
t643 = sin(qJ(2));
t669 = qJD(1) * t637;
t682 = pkin(9) * t636;
t623 = (-pkin(2) * t648 - t643 * t682) * t669;
t666 = qJD(1) * qJD(2);
t629 = (qJDD(1) * t643 + t648 * t666) * t637;
t633 = t639 * qJDD(1) + qJDD(2);
t644 = sin(qJ(1));
t649 = cos(qJ(1));
t631 = t644 * g(1) - t649 * g(2);
t650 = qJD(1) ^ 2;
t683 = pkin(8) * t637;
t626 = qJDD(1) * pkin(1) + t650 * t683 + t631;
t632 = -t649 * g(1) - t644 * g(2);
t627 = -t650 * pkin(1) + qJDD(1) * t683 + t632;
t672 = t639 * t648;
t659 = t626 * t672 - t643 * t627;
t668 = qJD(1) * t643;
t681 = pkin(9) * t638;
t580 = -t629 * t681 + t633 * pkin(2) + t634 * t619 + (-g(3) * t648 - t623 * t668) * t637 + t659;
t664 = t637 * t668;
t622 = t634 * pkin(2) - t664 * t681;
t630 = (qJDD(1) * t648 - t643 * t666) * t637;
t657 = t630 * t638 + t633 * t636;
t673 = t639 * t643;
t670 = t626 * t673 + t648 * t627;
t581 = -t634 * t622 + (-g(3) * t643 + t623 * t667) * t637 + t657 * pkin(9) + t670;
t680 = t639 * g(3);
t586 = -t629 * t682 - t630 * pkin(2) - t680 + (-t626 + (-t619 * t648 + t622 * t643) * qJD(1)) * t637;
t642 = sin(qJ(3));
t647 = cos(qJ(3));
t559 = -t642 * t581 + (t580 * t638 + t586 * t636) * t647;
t674 = t638 * t648;
t679 = t636 * t642;
t610 = t634 * t679 + (t642 * t674 + t643 * t647) * t669;
t596 = -t610 * qJD(3) - t642 * t629 + t657 * t647;
t678 = t636 * t647;
t609 = (-t642 * t643 + t647 * t674) * t669 + t634 * t678;
t677 = t637 * t643;
t676 = t637 * t648;
t675 = t638 * t642;
t560 = t580 * t675 + t647 * t581 + t586 * t679;
t598 = -t609 * mrSges(4,1) + t610 * mrSges(4,2);
t620 = t638 * t634 - t636 * t663 + qJD(3);
t604 = t620 * mrSges(4,1) - t610 * mrSges(4,3);
t611 = -t636 * t630 + t638 * t633 + qJDD(3);
t599 = -t609 * pkin(3) - t610 * pkin(10);
t618 = t620 ^ 2;
t553 = -t618 * pkin(3) + t611 * pkin(10) + t609 * t599 + t560;
t565 = -t636 * t580 + t638 * t586;
t597 = t609 * qJD(3) + t647 * t629 + t657 * t642;
t555 = (-t609 * t620 - t597) * pkin(10) + (t610 * t620 - t596) * pkin(3) + t565;
t641 = sin(qJ(4));
t646 = cos(qJ(4));
t550 = t646 * t553 + t641 * t555;
t602 = t646 * t610 + t641 * t620;
t569 = -t602 * qJD(4) - t641 * t597 + t646 * t611;
t601 = -t641 * t610 + t646 * t620;
t582 = -t601 * mrSges(5,1) + t602 * mrSges(5,2);
t608 = qJD(4) - t609;
t591 = t608 * mrSges(5,1) - t602 * mrSges(5,3);
t595 = qJDD(4) - t596;
t583 = -t601 * pkin(4) - t602 * pkin(11);
t607 = t608 ^ 2;
t547 = -t607 * pkin(4) + t595 * pkin(11) + t601 * t583 + t550;
t552 = -t611 * pkin(3) - t618 * pkin(10) + t610 * t599 - t559;
t570 = t601 * qJD(4) + t646 * t597 + t641 * t611;
t548 = (-t601 * t608 - t570) * pkin(11) + (t602 * t608 - t569) * pkin(4) + t552;
t640 = sin(qJ(5));
t645 = cos(qJ(5));
t544 = -t640 * t547 + t645 * t548;
t588 = -t640 * t602 + t645 * t608;
t558 = t588 * qJD(5) + t645 * t570 + t640 * t595;
t589 = t645 * t602 + t640 * t608;
t566 = -t588 * mrSges(6,1) + t589 * mrSges(6,2);
t568 = qJDD(5) - t569;
t600 = qJD(5) - t601;
t571 = -t600 * mrSges(6,2) + t588 * mrSges(6,3);
t542 = m(6) * t544 + t568 * mrSges(6,1) - t558 * mrSges(6,3) - t589 * t566 + t600 * t571;
t545 = t645 * t547 + t640 * t548;
t557 = -t589 * qJD(5) - t640 * t570 + t645 * t595;
t572 = t600 * mrSges(6,1) - t589 * mrSges(6,3);
t543 = m(6) * t545 - t568 * mrSges(6,2) + t557 * mrSges(6,3) + t588 * t566 - t600 * t572;
t660 = -t640 * t542 + t645 * t543;
t535 = m(5) * t550 - t595 * mrSges(5,2) + t569 * mrSges(5,3) + t601 * t582 - t608 * t591 + t660;
t549 = -t641 * t553 + t646 * t555;
t590 = -t608 * mrSges(5,2) + t601 * mrSges(5,3);
t546 = -t595 * pkin(4) - t607 * pkin(11) + t602 * t583 - t549;
t652 = -m(6) * t546 + t557 * mrSges(6,1) - t558 * mrSges(6,2) + t588 * t571 - t589 * t572;
t540 = m(5) * t549 + t595 * mrSges(5,1) - t570 * mrSges(5,3) - t602 * t582 + t608 * t590 + t652;
t661 = t646 * t535 - t641 * t540;
t526 = m(4) * t560 - t611 * mrSges(4,2) + t596 * mrSges(4,3) + t609 * t598 - t620 * t604 + t661;
t529 = t641 * t535 + t646 * t540;
t603 = -t620 * mrSges(4,2) + t609 * mrSges(4,3);
t528 = m(4) * t565 - t596 * mrSges(4,1) + t597 * mrSges(4,2) - t609 * t603 + t610 * t604 + t529;
t536 = t645 * t542 + t640 * t543;
t651 = -m(5) * t552 + t569 * mrSges(5,1) - t570 * mrSges(5,2) + t601 * t590 - t602 * t591 - t536;
t532 = m(4) * t559 + t611 * mrSges(4,1) - t597 * mrSges(4,3) - t610 * t598 + t620 * t603 + t651;
t515 = t638 * t647 * t532 + t526 * t675 - t636 * t528;
t605 = -g(3) * t676 + t659;
t625 = -t634 * mrSges(3,2) + mrSges(3,3) * t663;
t628 = (-mrSges(3,1) * t648 + mrSges(3,2) * t643) * t669;
t511 = m(3) * t605 + t633 * mrSges(3,1) - t629 * mrSges(3,3) + t634 * t625 - t628 * t664 + t515;
t514 = t526 * t679 + t638 * t528 + t532 * t678;
t615 = -t637 * t626 - t680;
t624 = t634 * mrSges(3,1) - mrSges(3,3) * t664;
t513 = m(3) * t615 - t630 * mrSges(3,1) + t629 * mrSges(3,2) + (t624 * t643 - t625 * t648) * t669 + t514;
t520 = t647 * t526 - t642 * t532;
t606 = -g(3) * t677 + t670;
t519 = m(3) * t606 - t633 * mrSges(3,2) + t630 * mrSges(3,3) - t634 * t624 + t628 * t663 + t520;
t501 = t511 * t672 - t637 * t513 + t519 * t673;
t499 = m(2) * t631 + qJDD(1) * mrSges(2,1) - t650 * mrSges(2,2) + t501;
t505 = -t643 * t511 + t648 * t519;
t504 = m(2) * t632 - t650 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t505;
t671 = t649 * t499 + t644 * t504;
t500 = t511 * t676 + t639 * t513 + t519 * t677;
t662 = -t644 * t499 + t649 * t504;
t561 = Ifges(6,5) * t589 + Ifges(6,6) * t588 + Ifges(6,3) * t600;
t563 = Ifges(6,1) * t589 + Ifges(6,4) * t588 + Ifges(6,5) * t600;
t537 = -mrSges(6,1) * t546 + mrSges(6,3) * t545 + Ifges(6,4) * t558 + Ifges(6,2) * t557 + Ifges(6,6) * t568 - t589 * t561 + t600 * t563;
t562 = Ifges(6,4) * t589 + Ifges(6,2) * t588 + Ifges(6,6) * t600;
t538 = mrSges(6,2) * t546 - mrSges(6,3) * t544 + Ifges(6,1) * t558 + Ifges(6,4) * t557 + Ifges(6,5) * t568 + t588 * t561 - t600 * t562;
t573 = Ifges(5,5) * t602 + Ifges(5,6) * t601 + Ifges(5,3) * t608;
t574 = Ifges(5,4) * t602 + Ifges(5,2) * t601 + Ifges(5,6) * t608;
t521 = mrSges(5,2) * t552 - mrSges(5,3) * t549 + Ifges(5,1) * t570 + Ifges(5,4) * t569 + Ifges(5,5) * t595 - pkin(11) * t536 - t640 * t537 + t645 * t538 + t601 * t573 - t608 * t574;
t575 = Ifges(5,1) * t602 + Ifges(5,4) * t601 + Ifges(5,5) * t608;
t522 = -mrSges(5,1) * t552 - mrSges(6,1) * t544 + mrSges(6,2) * t545 + mrSges(5,3) * t550 + Ifges(5,4) * t570 - Ifges(6,5) * t558 + Ifges(5,2) * t569 + Ifges(5,6) * t595 - Ifges(6,6) * t557 - Ifges(6,3) * t568 - pkin(4) * t536 - t589 * t562 + t588 * t563 - t602 * t573 + t608 * t575;
t593 = Ifges(4,4) * t610 + Ifges(4,2) * t609 + Ifges(4,6) * t620;
t594 = Ifges(4,1) * t610 + Ifges(4,4) * t609 + Ifges(4,5) * t620;
t506 = mrSges(4,1) * t559 - mrSges(4,2) * t560 + Ifges(4,5) * t597 + Ifges(4,6) * t596 + Ifges(4,3) * t611 + pkin(3) * t651 + pkin(10) * t661 + t641 * t521 + t646 * t522 + t610 * t593 - t609 * t594;
t612 = Ifges(3,3) * t634 + (Ifges(3,5) * t643 + Ifges(3,6) * t648) * t669;
t614 = Ifges(3,5) * t634 + (Ifges(3,1) * t643 + Ifges(3,4) * t648) * t669;
t592 = Ifges(4,5) * t610 + Ifges(4,6) * t609 + Ifges(4,3) * t620;
t507 = mrSges(4,2) * t565 - mrSges(4,3) * t559 + Ifges(4,1) * t597 + Ifges(4,4) * t596 + Ifges(4,5) * t611 - pkin(10) * t529 + t646 * t521 - t641 * t522 + t609 * t592 - t620 * t593;
t508 = Ifges(4,4) * t597 + Ifges(4,2) * t596 + Ifges(4,6) * t611 - t610 * t592 + t620 * t594 - mrSges(4,1) * t565 + mrSges(4,3) * t560 - Ifges(5,5) * t570 - Ifges(5,6) * t569 - Ifges(5,3) * t595 - t602 * t574 + t601 * t575 - mrSges(5,1) * t549 + mrSges(5,2) * t550 - t640 * t538 - t645 * t537 - pkin(4) * t652 - pkin(11) * t660 - pkin(3) * t529;
t653 = pkin(9) * t520 + t507 * t642 + t508 * t647;
t496 = -mrSges(3,1) * t615 + mrSges(3,3) * t606 + Ifges(3,4) * t629 + Ifges(3,2) * t630 + Ifges(3,6) * t633 - pkin(2) * t514 - t636 * t506 - t612 * t664 + t634 * t614 + t653 * t638;
t613 = Ifges(3,6) * t634 + (Ifges(3,4) * t643 + Ifges(3,2) * t648) * t669;
t497 = t612 * t663 + mrSges(3,2) * t615 - mrSges(3,3) * t605 + Ifges(3,1) * t629 + Ifges(3,4) * t630 + Ifges(3,5) * t633 + t647 * t507 - t642 * t508 - t634 * t613 + (-t514 * t636 - t515 * t638) * pkin(9);
t654 = pkin(8) * t505 + t496 * t648 + t497 * t643;
t495 = mrSges(3,1) * t605 - mrSges(3,2) * t606 + Ifges(3,5) * t629 + Ifges(3,6) * t630 + Ifges(3,3) * t633 + pkin(2) * t515 + t638 * t506 + (t613 * t643 - t614 * t648) * t669 + t653 * t636;
t494 = -mrSges(2,2) * g(3) - mrSges(2,3) * t631 + Ifges(2,5) * qJDD(1) - t650 * Ifges(2,6) - t643 * t496 + t648 * t497 + (-t500 * t637 - t501 * t639) * pkin(8);
t493 = mrSges(2,1) * g(3) + mrSges(2,3) * t632 + t650 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t500 - t637 * t495 + t654 * t639;
t1 = [-m(1) * g(1) + t662; -m(1) * g(2) + t671; (-m(1) - m(2)) * g(3) + t500; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t671 - t644 * t493 + t649 * t494; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t662 + t649 * t493 + t644 * t494; -mrSges(1,1) * g(2) + mrSges(2,1) * t631 + mrSges(1,2) * g(1) - mrSges(2,2) * t632 + Ifges(2,3) * qJDD(1) + pkin(1) * t501 + t639 * t495 + t654 * t637;];
tauB = t1;
