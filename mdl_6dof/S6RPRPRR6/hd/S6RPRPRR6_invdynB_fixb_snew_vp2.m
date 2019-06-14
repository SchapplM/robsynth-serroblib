% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:02:30
% EndTime: 2019-05-05 19:02:53
% DurationCPUTime: 23.12s
% Computational Cost: add. (361256->365), mult. (892356->460), div. (0->0), fcn. (686634->12), ass. (0->151)
t730 = qJD(1) ^ 2;
t760 = cos(qJ(3));
t721 = cos(pkin(10));
t759 = pkin(2) * t721;
t719 = sin(pkin(10));
t758 = mrSges(3,2) * t719;
t717 = t721 ^ 2;
t757 = t717 * t730;
t725 = sin(qJ(1));
t728 = cos(qJ(1));
t707 = -t728 * g(1) - t725 * g(2);
t703 = -t730 * pkin(1) + qJDD(1) * qJ(2) + t707;
t752 = qJD(1) * qJD(2);
t748 = -t721 * g(3) - 0.2e1 * t719 * t752;
t668 = (-pkin(7) * qJDD(1) + t730 * t759 - t703) * t719 + t748;
t689 = -t719 * g(3) + (t703 + 0.2e1 * t752) * t721;
t750 = qJDD(1) * t721;
t675 = -pkin(2) * t757 + pkin(7) * t750 + t689;
t724 = sin(qJ(3));
t651 = t724 * t668 + t760 * t675;
t749 = t721 * t760;
t753 = t719 * qJD(1);
t701 = -qJD(1) * t749 + t724 * t753;
t735 = t760 * t719 + t721 * t724;
t702 = t735 * qJD(1);
t682 = t701 * mrSges(4,1) + t702 * mrSges(4,2);
t751 = qJDD(1) * t719;
t754 = t702 * qJD(3);
t686 = -qJDD(1) * t749 + t724 * t751 + t754;
t696 = qJD(3) * mrSges(4,1) - t702 * mrSges(4,3);
t681 = t701 * pkin(3) - t702 * qJ(4);
t729 = qJD(3) ^ 2;
t637 = -t729 * pkin(3) + qJDD(3) * qJ(4) - t701 * t681 + t651;
t716 = t719 ^ 2;
t706 = t725 * g(1) - t728 * g(2);
t741 = qJDD(2) - t706;
t685 = (-pkin(1) - t759) * qJDD(1) + (-qJ(2) + (-t716 - t717) * pkin(7)) * t730 + t741;
t755 = t701 * qJD(3);
t687 = t735 * qJDD(1) - t755;
t642 = (-t687 + t755) * qJ(4) + (t686 + t754) * pkin(3) + t685;
t718 = sin(pkin(11));
t720 = cos(pkin(11));
t694 = t718 * qJD(3) + t720 * t702;
t623 = -0.2e1 * qJD(4) * t694 - t718 * t637 + t720 * t642;
t674 = t718 * qJDD(3) + t720 * t687;
t693 = t720 * qJD(3) - t718 * t702;
t616 = (t693 * t701 - t674) * pkin(8) + (t693 * t694 + t686) * pkin(4) + t623;
t624 = 0.2e1 * qJD(4) * t693 + t720 * t637 + t718 * t642;
t672 = t701 * pkin(4) - t694 * pkin(8);
t673 = t720 * qJDD(3) - t718 * t687;
t692 = t693 ^ 2;
t618 = -t692 * pkin(4) + t673 * pkin(8) - t701 * t672 + t624;
t723 = sin(qJ(5));
t727 = cos(qJ(5));
t610 = t727 * t616 - t723 * t618;
t661 = t727 * t693 - t723 * t694;
t636 = t661 * qJD(5) + t723 * t673 + t727 * t674;
t662 = t723 * t693 + t727 * t694;
t684 = qJDD(5) + t686;
t699 = qJD(5) + t701;
t608 = (t661 * t699 - t636) * pkin(9) + (t661 * t662 + t684) * pkin(5) + t610;
t611 = t723 * t616 + t727 * t618;
t635 = -t662 * qJD(5) + t727 * t673 - t723 * t674;
t654 = t699 * pkin(5) - t662 * pkin(9);
t660 = t661 ^ 2;
t609 = -t660 * pkin(5) + t635 * pkin(9) - t699 * t654 + t611;
t722 = sin(qJ(6));
t726 = cos(qJ(6));
t606 = t726 * t608 - t722 * t609;
t647 = t726 * t661 - t722 * t662;
t622 = t647 * qJD(6) + t722 * t635 + t726 * t636;
t648 = t722 * t661 + t726 * t662;
t631 = -t647 * mrSges(7,1) + t648 * mrSges(7,2);
t698 = qJD(6) + t699;
t640 = -t698 * mrSges(7,2) + t647 * mrSges(7,3);
t680 = qJDD(6) + t684;
t602 = m(7) * t606 + t680 * mrSges(7,1) - t622 * mrSges(7,3) - t648 * t631 + t698 * t640;
t607 = t722 * t608 + t726 * t609;
t621 = -t648 * qJD(6) + t726 * t635 - t722 * t636;
t641 = t698 * mrSges(7,1) - t648 * mrSges(7,3);
t603 = m(7) * t607 - t680 * mrSges(7,2) + t621 * mrSges(7,3) + t647 * t631 - t698 * t641;
t596 = t726 * t602 + t722 * t603;
t649 = -t661 * mrSges(6,1) + t662 * mrSges(6,2);
t652 = -t699 * mrSges(6,2) + t661 * mrSges(6,3);
t594 = m(6) * t610 + t684 * mrSges(6,1) - t636 * mrSges(6,3) - t662 * t649 + t699 * t652 + t596;
t653 = t699 * mrSges(6,1) - t662 * mrSges(6,3);
t742 = -t722 * t602 + t726 * t603;
t595 = m(6) * t611 - t684 * mrSges(6,2) + t635 * mrSges(6,3) + t661 * t649 - t699 * t653 + t742;
t590 = t727 * t594 + t723 * t595;
t663 = -t693 * mrSges(5,1) + t694 * mrSges(5,2);
t670 = -t701 * mrSges(5,2) + t693 * mrSges(5,3);
t588 = m(5) * t623 + t686 * mrSges(5,1) - t674 * mrSges(5,3) - t694 * t663 + t701 * t670 + t590;
t671 = t701 * mrSges(5,1) - t694 * mrSges(5,3);
t743 = -t723 * t594 + t727 * t595;
t589 = m(5) * t624 - t686 * mrSges(5,2) + t673 * mrSges(5,3) + t693 * t663 - t701 * t671 + t743;
t744 = -t718 * t588 + t720 * t589;
t581 = m(4) * t651 - qJDD(3) * mrSges(4,2) - t686 * mrSges(4,3) - qJD(3) * t696 - t701 * t682 + t744;
t650 = t760 * t668 - t724 * t675;
t695 = -qJD(3) * mrSges(4,2) - t701 * mrSges(4,3);
t634 = -qJDD(3) * pkin(3) - t729 * qJ(4) + t702 * t681 + qJDD(4) - t650;
t625 = -t673 * pkin(4) - t692 * pkin(8) + t694 * t672 + t634;
t613 = -t635 * pkin(5) - t660 * pkin(9) + t662 * t654 + t625;
t737 = m(7) * t613 - t621 * mrSges(7,1) + t622 * mrSges(7,2) - t647 * t640 + t648 * t641;
t733 = m(6) * t625 - t635 * mrSges(6,1) + t636 * mrSges(6,2) - t661 * t652 + t662 * t653 + t737;
t731 = -m(5) * t634 + t673 * mrSges(5,1) - t674 * mrSges(5,2) + t693 * t670 - t694 * t671 - t733;
t605 = m(4) * t650 + qJDD(3) * mrSges(4,1) - t687 * mrSges(4,3) + qJD(3) * t695 - t702 * t682 + t731;
t576 = t724 * t581 + t760 * t605;
t688 = -t719 * t703 + t748;
t736 = mrSges(3,3) * qJDD(1) + t730 * (-mrSges(3,1) * t721 + t758);
t574 = m(3) * t688 - t736 * t719 + t576;
t745 = t760 * t581 - t724 * t605;
t575 = m(3) * t689 + t736 * t721 + t745;
t746 = -t719 * t574 + t721 * t575;
t566 = m(2) * t707 - t730 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t746;
t700 = -qJDD(1) * pkin(1) - t730 * qJ(2) + t741;
t582 = t720 * t588 + t718 * t589;
t734 = m(4) * t685 + t686 * mrSges(4,1) + t687 * mrSges(4,2) + t701 * t695 + t702 * t696 + t582;
t732 = -m(3) * t700 + mrSges(3,1) * t750 - t734 + (t716 * t730 + t757) * mrSges(3,3);
t578 = -t730 * mrSges(2,2) + t732 + m(2) * t706 + (mrSges(2,1) - t758) * qJDD(1);
t756 = t725 * t566 + t728 * t578;
t567 = t721 * t574 + t719 * t575;
t747 = t728 * t566 - t725 * t578;
t740 = Ifges(3,1) * t719 + Ifges(3,4) * t721;
t739 = Ifges(3,4) * t719 + Ifges(3,2) * t721;
t738 = Ifges(3,5) * t719 + Ifges(3,6) * t721;
t705 = t738 * qJD(1);
t678 = Ifges(4,1) * t702 - Ifges(4,4) * t701 + Ifges(4,5) * qJD(3);
t677 = Ifges(4,4) * t702 - Ifges(4,2) * t701 + Ifges(4,6) * qJD(3);
t676 = Ifges(4,5) * t702 - Ifges(4,6) * t701 + Ifges(4,3) * qJD(3);
t657 = Ifges(5,1) * t694 + Ifges(5,4) * t693 + Ifges(5,5) * t701;
t656 = Ifges(5,4) * t694 + Ifges(5,2) * t693 + Ifges(5,6) * t701;
t655 = Ifges(5,5) * t694 + Ifges(5,6) * t693 + Ifges(5,3) * t701;
t645 = Ifges(6,1) * t662 + Ifges(6,4) * t661 + Ifges(6,5) * t699;
t644 = Ifges(6,4) * t662 + Ifges(6,2) * t661 + Ifges(6,6) * t699;
t643 = Ifges(6,5) * t662 + Ifges(6,6) * t661 + Ifges(6,3) * t699;
t628 = Ifges(7,1) * t648 + Ifges(7,4) * t647 + Ifges(7,5) * t698;
t627 = Ifges(7,4) * t648 + Ifges(7,2) * t647 + Ifges(7,6) * t698;
t626 = Ifges(7,5) * t648 + Ifges(7,6) * t647 + Ifges(7,3) * t698;
t598 = mrSges(7,2) * t613 - mrSges(7,3) * t606 + Ifges(7,1) * t622 + Ifges(7,4) * t621 + Ifges(7,5) * t680 + t647 * t626 - t698 * t627;
t597 = -mrSges(7,1) * t613 + mrSges(7,3) * t607 + Ifges(7,4) * t622 + Ifges(7,2) * t621 + Ifges(7,6) * t680 - t648 * t626 + t698 * t628;
t584 = mrSges(6,2) * t625 - mrSges(6,3) * t610 + Ifges(6,1) * t636 + Ifges(6,4) * t635 + Ifges(6,5) * t684 - pkin(9) * t596 - t722 * t597 + t726 * t598 + t661 * t643 - t699 * t644;
t583 = -mrSges(6,1) * t625 + mrSges(6,3) * t611 + Ifges(6,4) * t636 + Ifges(6,2) * t635 + Ifges(6,6) * t684 - pkin(5) * t737 + pkin(9) * t742 + t726 * t597 + t722 * t598 - t662 * t643 + t699 * t645;
t570 = mrSges(5,2) * t634 - mrSges(5,3) * t623 + Ifges(5,1) * t674 + Ifges(5,4) * t673 + Ifges(5,5) * t686 - pkin(8) * t590 - t723 * t583 + t727 * t584 + t693 * t655 - t701 * t656;
t569 = -mrSges(5,1) * t634 + mrSges(5,3) * t624 + Ifges(5,4) * t674 + Ifges(5,2) * t673 + Ifges(5,6) * t686 - pkin(4) * t733 + pkin(8) * t743 + t727 * t583 + t723 * t584 - t694 * t655 + t701 * t657;
t568 = mrSges(6,2) * t611 - mrSges(6,1) * t610 + Ifges(4,4) * t687 + t693 * t657 - t694 * t656 + mrSges(7,2) * t607 - mrSges(7,1) * t606 + t661 * t645 - t662 * t644 - pkin(5) * t596 - pkin(4) * t590 - Ifges(6,6) * t635 - Ifges(6,5) * t636 + qJD(3) * t678 - Ifges(7,3) * t680 - Ifges(6,3) * t684 - mrSges(4,1) * t685 - pkin(3) * t582 - Ifges(5,6) * t673 - Ifges(5,5) * t674 + mrSges(4,3) * t651 - t702 * t676 + mrSges(5,2) * t624 - Ifges(7,6) * t621 - Ifges(7,5) * t622 - mrSges(5,1) * t623 + Ifges(4,6) * qJDD(3) + t647 * t628 - t648 * t627 + (-Ifges(5,3) - Ifges(4,2)) * t686;
t563 = mrSges(4,2) * t685 - mrSges(4,3) * t650 + Ifges(4,1) * t687 - Ifges(4,4) * t686 + Ifges(4,5) * qJDD(3) - qJ(4) * t582 - qJD(3) * t677 - t718 * t569 + t720 * t570 - t701 * t676;
t562 = t721 * qJD(1) * t705 + mrSges(3,2) * t700 - mrSges(3,3) * t688 - pkin(7) * t576 + t740 * qJDD(1) + t760 * t563 - t724 * t568;
t561 = mrSges(2,1) * g(3) - pkin(1) * t567 + mrSges(2,3) * t707 - pkin(2) * t576 - mrSges(3,1) * t688 + mrSges(3,2) * t689 - pkin(3) * t731 - qJ(4) * t744 - t718 * t570 - t720 * t569 - mrSges(4,1) * t650 + mrSges(4,2) * t651 - Ifges(4,5) * t687 + Ifges(4,6) * t686 - Ifges(4,3) * qJDD(3) - t702 * t677 - t701 * t678 + (Ifges(2,6) - t738) * qJDD(1) + (-t719 * t739 + t721 * t740 + Ifges(2,5)) * t730;
t560 = -mrSges(3,1) * t700 + mrSges(3,3) * t689 - pkin(2) * t734 + pkin(7) * t745 + t739 * qJDD(1) + t724 * t563 + t760 * t568 - t705 * t753;
t559 = -mrSges(2,2) * g(3) - mrSges(2,3) * t706 + Ifges(2,5) * qJDD(1) - t730 * Ifges(2,6) - qJ(2) * t567 - t719 * t560 + t721 * t562;
t1 = [-m(1) * g(1) + t747; -m(1) * g(2) + t756; (-m(1) - m(2)) * g(3) + t567; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t756 + t728 * t559 - t725 * t561; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t747 + t725 * t559 + t728 * t561; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t706 - mrSges(2,2) * t707 + t719 * t562 + t721 * t560 + pkin(1) * (-mrSges(3,2) * t751 + t732) + qJ(2) * t746;];
tauB  = t1;
