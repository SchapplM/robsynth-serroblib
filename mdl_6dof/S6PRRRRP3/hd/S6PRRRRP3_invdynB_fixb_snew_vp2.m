% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:43:20
% EndTime: 2019-05-05 09:43:34
% DurationCPUTime: 10.96s
% Computational Cost: add. (181697->319), mult. (352961->395), div. (0->0), fcn. (251386->12), ass. (0->133)
t707 = Ifges(6,1) + Ifges(7,1);
t704 = Ifges(6,4) + Ifges(7,4);
t703 = Ifges(6,5) + Ifges(7,5);
t706 = Ifges(6,2) + Ifges(7,2);
t702 = -Ifges(6,6) - Ifges(7,6);
t705 = -Ifges(6,3) - Ifges(7,3);
t665 = sin(pkin(11));
t667 = cos(pkin(11));
t656 = g(1) * t665 - g(2) * t667;
t657 = -g(1) * t667 - g(2) * t665;
t664 = -g(3) + qJDD(1);
t666 = sin(pkin(6));
t668 = cos(pkin(6));
t672 = sin(qJ(2));
t676 = cos(qJ(2));
t618 = -t672 * t657 + (t656 * t668 + t664 * t666) * t676;
t678 = qJD(2) ^ 2;
t699 = t668 * t672;
t700 = t666 * t672;
t619 = t656 * t699 + t676 * t657 + t664 * t700;
t612 = -pkin(2) * t678 + qJDD(2) * pkin(8) + t619;
t636 = -t656 * t666 + t664 * t668;
t671 = sin(qJ(3));
t675 = cos(qJ(3));
t605 = t675 * t612 + t671 * t636;
t653 = (-pkin(3) * t675 - pkin(9) * t671) * qJD(2);
t677 = qJD(3) ^ 2;
t693 = qJD(2) * t675;
t591 = -pkin(3) * t677 + qJDD(3) * pkin(9) + t653 * t693 + t605;
t611 = -qJDD(2) * pkin(2) - t678 * pkin(8) - t618;
t692 = qJD(2) * qJD(3);
t689 = t675 * t692;
t654 = qJDD(2) * t671 + t689;
t663 = t671 * t692;
t655 = qJDD(2) * t675 - t663;
t594 = (-t654 - t689) * pkin(9) + (-t655 + t663) * pkin(3) + t611;
t670 = sin(qJ(4));
t674 = cos(qJ(4));
t579 = -t670 * t591 + t674 * t594;
t694 = qJD(2) * t671;
t650 = qJD(3) * t674 - t670 * t694;
t627 = qJD(4) * t650 + qJDD(3) * t670 + t654 * t674;
t647 = qJDD(4) - t655;
t651 = qJD(3) * t670 + t674 * t694;
t662 = qJD(4) - t693;
t576 = (t650 * t662 - t627) * pkin(10) + (t650 * t651 + t647) * pkin(4) + t579;
t580 = t674 * t591 + t670 * t594;
t626 = -qJD(4) * t651 + qJDD(3) * t674 - t654 * t670;
t635 = pkin(4) * t662 - pkin(10) * t651;
t646 = t650 ^ 2;
t578 = -pkin(4) * t646 + pkin(10) * t626 - t635 * t662 + t580;
t669 = sin(qJ(5));
t673 = cos(qJ(5));
t570 = t673 * t576 - t669 * t578;
t629 = t650 * t673 - t651 * t669;
t587 = qJD(5) * t629 + t626 * t669 + t627 * t673;
t630 = t650 * t669 + t651 * t673;
t606 = -mrSges(7,1) * t629 + mrSges(7,2) * t630;
t607 = -mrSges(6,1) * t629 + mrSges(6,2) * t630;
t661 = qJD(5) + t662;
t614 = -mrSges(6,2) * t661 + mrSges(6,3) * t629;
t643 = qJDD(5) + t647;
t567 = -0.2e1 * qJD(6) * t630 + (t629 * t661 - t587) * qJ(6) + (t629 * t630 + t643) * pkin(5) + t570;
t613 = -mrSges(7,2) * t661 + mrSges(7,3) * t629;
t691 = m(7) * t567 + t643 * mrSges(7,1) + t661 * t613;
t559 = m(6) * t570 + t643 * mrSges(6,1) + t661 * t614 + (-t606 - t607) * t630 + (-mrSges(6,3) - mrSges(7,3)) * t587 + t691;
t571 = t669 * t576 + t673 * t578;
t586 = -qJD(5) * t630 + t626 * t673 - t627 * t669;
t616 = mrSges(7,1) * t661 - mrSges(7,3) * t630;
t617 = mrSges(6,1) * t661 - mrSges(6,3) * t630;
t615 = pkin(5) * t661 - qJ(6) * t630;
t628 = t629 ^ 2;
t569 = -pkin(5) * t628 + qJ(6) * t586 + 0.2e1 * qJD(6) * t629 - t615 * t661 + t571;
t690 = m(7) * t569 + t586 * mrSges(7,3) + t629 * t606;
t562 = m(6) * t571 + t586 * mrSges(6,3) + t629 * t607 + (-t616 - t617) * t661 + (-mrSges(6,2) - mrSges(7,2)) * t643 + t690;
t557 = t673 * t559 + t669 * t562;
t631 = -mrSges(5,1) * t650 + mrSges(5,2) * t651;
t633 = -mrSges(5,2) * t662 + mrSges(5,3) * t650;
t554 = m(5) * t579 + mrSges(5,1) * t647 - mrSges(5,3) * t627 - t631 * t651 + t633 * t662 + t557;
t634 = mrSges(5,1) * t662 - mrSges(5,3) * t651;
t685 = -t559 * t669 + t673 * t562;
t555 = m(5) * t580 - mrSges(5,2) * t647 + mrSges(5,3) * t626 + t631 * t650 - t634 * t662 + t685;
t551 = t554 * t674 + t555 * t670;
t658 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t694;
t659 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t693;
t680 = -m(4) * t611 + t655 * mrSges(4,1) - mrSges(4,2) * t654 - t658 * t694 + t659 * t693 - t551;
t547 = m(3) * t618 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t678 + t680;
t701 = t547 * t676;
t652 = (-mrSges(4,1) * t675 + mrSges(4,2) * t671) * qJD(2);
t686 = -t554 * t670 + t674 * t555;
t550 = m(4) * t605 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t655 - qJD(3) * t658 + t652 * t693 + t686;
t604 = -t671 * t612 + t636 * t675;
t590 = -qJDD(3) * pkin(3) - pkin(9) * t677 + t653 * t694 - t604;
t581 = -pkin(4) * t626 - pkin(10) * t646 + t651 * t635 + t590;
t573 = -pkin(5) * t586 - qJ(6) * t628 + t615 * t630 + qJDD(6) + t581;
t684 = m(7) * t573 - t586 * mrSges(7,1) + t587 * mrSges(7,2) - t629 * t613 + t630 * t616;
t681 = m(6) * t581 - t586 * mrSges(6,1) + t587 * mrSges(6,2) - t629 * t614 + t630 * t617 + t684;
t679 = -m(5) * t590 + t626 * mrSges(5,1) - t627 * mrSges(5,2) + t650 * t633 - t651 * t634 - t681;
t564 = m(4) * t604 + qJDD(3) * mrSges(4,1) - t654 * mrSges(4,3) + qJD(3) * t659 - t652 * t694 + t679;
t687 = t675 * t550 - t564 * t671;
t541 = m(3) * t619 - mrSges(3,1) * t678 - qJDD(2) * mrSges(3,2) + t687;
t544 = t671 * t550 + t675 * t564;
t543 = m(3) * t636 + t544;
t530 = t541 * t699 - t543 * t666 + t668 * t701;
t528 = m(2) * t656 + t530;
t534 = t676 * t541 - t547 * t672;
t533 = m(2) * t657 + t534;
t698 = t667 * t528 + t665 * t533;
t697 = t702 * t629 - t703 * t630 + t705 * t661;
t696 = -t706 * t629 - t704 * t630 + t702 * t661;
t695 = t704 * t629 + t707 * t630 + t703 * t661;
t529 = t541 * t700 + t668 * t543 + t666 * t701;
t688 = -t528 * t665 + t667 * t533;
t552 = -mrSges(6,1) * t581 + mrSges(6,3) * t571 - mrSges(7,1) * t573 + mrSges(7,3) * t569 - pkin(5) * t684 + qJ(6) * t690 + (-qJ(6) * t616 + t695) * t661 + (-mrSges(7,2) * qJ(6) - t702) * t643 + t697 * t630 + t704 * t587 + t706 * t586;
t565 = -t587 * mrSges(7,3) - t630 * t606 + t691;
t556 = mrSges(6,2) * t581 + mrSges(7,2) * t573 - mrSges(6,3) * t570 - mrSges(7,3) * t567 - qJ(6) * t565 + t704 * t586 + t707 * t587 - t697 * t629 + t703 * t643 + t696 * t661;
t620 = Ifges(5,5) * t651 + Ifges(5,6) * t650 + Ifges(5,3) * t662;
t622 = Ifges(5,1) * t651 + Ifges(5,4) * t650 + Ifges(5,5) * t662;
t536 = -mrSges(5,1) * t590 + mrSges(5,3) * t580 + Ifges(5,4) * t627 + Ifges(5,2) * t626 + Ifges(5,6) * t647 - pkin(4) * t681 + pkin(10) * t685 + t673 * t552 + t669 * t556 - t651 * t620 + t662 * t622;
t621 = Ifges(5,4) * t651 + Ifges(5,2) * t650 + Ifges(5,6) * t662;
t537 = mrSges(5,2) * t590 - mrSges(5,3) * t579 + Ifges(5,1) * t627 + Ifges(5,4) * t626 + Ifges(5,5) * t647 - pkin(10) * t557 - t552 * t669 + t556 * t673 + t620 * t650 - t621 * t662;
t640 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t671 + Ifges(4,6) * t675) * qJD(2);
t641 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t671 + Ifges(4,2) * t675) * qJD(2);
t526 = mrSges(4,2) * t611 - mrSges(4,3) * t604 + Ifges(4,1) * t654 + Ifges(4,4) * t655 + Ifges(4,5) * qJDD(3) - pkin(9) * t551 - qJD(3) * t641 - t536 * t670 + t537 * t674 + t640 * t693;
t642 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t671 + Ifges(4,4) * t675) * qJD(2);
t535 = t650 * t622 - t651 * t621 + Ifges(4,4) * t654 + Ifges(4,2) * t655 - Ifges(5,3) * t647 + qJD(3) * t642 - Ifges(5,6) * t626 - Ifges(5,5) * t627 + mrSges(4,3) * t605 - mrSges(4,1) * t611 + mrSges(5,2) * t580 - mrSges(5,1) * t579 + mrSges(7,2) * t569 - mrSges(6,1) * t570 + mrSges(6,2) * t571 - mrSges(7,1) * t567 - pkin(5) * t565 - pkin(4) * t557 - pkin(3) * t551 + t705 * t643 + Ifges(4,6) * qJDD(3) - t640 * t694 + t695 * t629 + t696 * t630 + t702 * t586 - t703 * t587;
t524 = mrSges(3,2) * t636 - mrSges(3,3) * t618 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t678 - pkin(8) * t544 + t526 * t675 - t535 * t671;
t525 = Ifges(3,6) * qJDD(2) + t678 * Ifges(3,5) - mrSges(3,1) * t636 + mrSges(3,3) * t619 - Ifges(4,5) * t654 - Ifges(4,6) * t655 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t604 + mrSges(4,2) * t605 - t670 * t537 - t674 * t536 - pkin(3) * t679 - pkin(9) * t686 - pkin(2) * t544 + (-t641 * t671 + t642 * t675) * qJD(2);
t682 = pkin(7) * t534 + t524 * t672 + t525 * t676;
t523 = mrSges(3,1) * t618 - mrSges(3,2) * t619 + Ifges(3,3) * qJDD(2) + pkin(2) * t680 + pkin(8) * t687 + t671 * t526 + t675 * t535;
t522 = mrSges(2,2) * t664 - mrSges(2,3) * t656 + t676 * t524 - t672 * t525 + (-t529 * t666 - t530 * t668) * pkin(7);
t521 = -mrSges(2,1) * t664 + mrSges(2,3) * t657 - pkin(1) * t529 - t666 * t523 + t682 * t668;
t1 = [-m(1) * g(1) + t688; -m(1) * g(2) + t698; -m(1) * g(3) + m(2) * t664 + t529; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t698 - t665 * t521 + t667 * t522; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t688 + t667 * t521 + t665 * t522; -mrSges(1,1) * g(2) + mrSges(2,1) * t656 + mrSges(1,2) * g(1) - mrSges(2,2) * t657 + pkin(1) * t530 + t668 * t523 + t682 * t666;];
tauB  = t1;
