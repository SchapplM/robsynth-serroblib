% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR4
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
% Datum: 2019-05-06 03:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:13:50
% EndTime: 2019-05-06 03:14:17
% DurationCPUTime: 27.82s
% Computational Cost: add. (438862->365), mult. (1096892->456), div. (0->0), fcn. (888076->12), ass. (0->150)
t732 = qJD(1) ^ 2;
t721 = cos(pkin(11));
t759 = pkin(2) * t721;
t720 = sin(pkin(11));
t758 = mrSges(3,2) * t720;
t718 = t721 ^ 2;
t757 = t718 * t732;
t726 = sin(qJ(1));
t731 = cos(qJ(1));
t706 = -t731 * g(1) - t726 * g(2);
t702 = -t732 * pkin(1) + qJDD(1) * qJ(2) + t706;
t753 = qJD(1) * qJD(2);
t751 = -t721 * g(3) - 0.2e1 * t720 * t753;
t677 = (-pkin(7) * qJDD(1) + t732 * t759 - t702) * t720 + t751;
t693 = -t720 * g(3) + (t702 + 0.2e1 * t753) * t721;
t752 = qJDD(1) * t721;
t678 = -pkin(2) * t757 + pkin(7) * t752 + t693;
t725 = sin(qJ(3));
t730 = cos(qJ(3));
t657 = t730 * t677 - t725 * t678;
t740 = t720 * t730 + t721 * t725;
t739 = -t720 * t725 + t721 * t730;
t700 = t739 * qJD(1);
t754 = t700 * qJD(3);
t691 = qJDD(1) * t740 + t754;
t701 = t740 * qJD(1);
t639 = (-t691 + t754) * pkin(8) + (t700 * t701 + qJDD(3)) * pkin(3) + t657;
t658 = t725 * t677 + t730 * t678;
t690 = -t701 * qJD(3) + qJDD(1) * t739;
t696 = qJD(3) * pkin(3) - t701 * pkin(8);
t699 = t700 ^ 2;
t648 = -t699 * pkin(3) + t690 * pkin(8) - qJD(3) * t696 + t658;
t724 = sin(qJ(4));
t729 = cos(qJ(4));
t623 = t729 * t639 - t724 * t648;
t683 = t729 * t700 - t724 * t701;
t654 = t683 * qJD(4) + t724 * t690 + t729 * t691;
t684 = t724 * t700 + t729 * t701;
t716 = qJDD(3) + qJDD(4);
t719 = qJD(3) + qJD(4);
t615 = (t683 * t719 - t654) * pkin(9) + (t683 * t684 + t716) * pkin(4) + t623;
t624 = t724 * t639 + t729 * t648;
t653 = -t684 * qJD(4) + t729 * t690 - t724 * t691;
t675 = t719 * pkin(4) - t684 * pkin(9);
t679 = t683 ^ 2;
t617 = -t679 * pkin(4) + t653 * pkin(9) - t719 * t675 + t624;
t723 = sin(qJ(5));
t728 = cos(qJ(5));
t613 = t723 * t615 + t728 * t617;
t669 = t723 * t683 + t728 * t684;
t628 = -t669 * qJD(5) + t728 * t653 - t723 * t654;
t668 = t728 * t683 - t723 * t684;
t645 = -t668 * mrSges(6,1) + t669 * mrSges(6,2);
t714 = qJD(5) + t719;
t660 = t714 * mrSges(6,1) - t669 * mrSges(6,3);
t713 = qJDD(5) + t716;
t647 = -t668 * pkin(5) - t669 * pkin(10);
t712 = t714 ^ 2;
t610 = -t712 * pkin(5) + t713 * pkin(10) + t668 * t647 + t613;
t717 = t720 ^ 2;
t705 = t726 * g(1) - t731 * g(2);
t744 = qJDD(2) - t705;
t689 = (-pkin(1) - t759) * qJDD(1) + (-qJ(2) + (-t717 - t718) * pkin(7)) * t732 + t744;
t650 = -t690 * pkin(3) - t699 * pkin(8) + t701 * t696 + t689;
t622 = -t653 * pkin(4) - t679 * pkin(9) + t684 * t675 + t650;
t629 = t668 * qJD(5) + t723 * t653 + t728 * t654;
t611 = t622 + (-t668 * t714 - t629) * pkin(10) + (t669 * t714 - t628) * pkin(5);
t722 = sin(qJ(6));
t727 = cos(qJ(6));
t607 = -t722 * t610 + t727 * t611;
t655 = -t722 * t669 + t727 * t714;
t620 = t655 * qJD(6) + t727 * t629 + t722 * t713;
t627 = qJDD(6) - t628;
t656 = t727 * t669 + t722 * t714;
t634 = -t655 * mrSges(7,1) + t656 * mrSges(7,2);
t664 = qJD(6) - t668;
t637 = -t664 * mrSges(7,2) + t655 * mrSges(7,3);
t605 = m(7) * t607 + t627 * mrSges(7,1) - t620 * mrSges(7,3) - t656 * t634 + t664 * t637;
t608 = t727 * t610 + t722 * t611;
t619 = -t656 * qJD(6) - t722 * t629 + t727 * t713;
t638 = t664 * mrSges(7,1) - t656 * mrSges(7,3);
t606 = m(7) * t608 - t627 * mrSges(7,2) + t619 * mrSges(7,3) + t655 * t634 - t664 * t638;
t745 = -t722 * t605 + t727 * t606;
t596 = m(6) * t613 - t713 * mrSges(6,2) + t628 * mrSges(6,3) + t668 * t645 - t714 * t660 + t745;
t612 = t728 * t615 - t723 * t617;
t659 = -t714 * mrSges(6,2) + t668 * mrSges(6,3);
t609 = -t713 * pkin(5) - t712 * pkin(10) + t669 * t647 - t612;
t736 = -m(7) * t609 + t619 * mrSges(7,1) - t620 * mrSges(7,2) + t655 * t637 - t656 * t638;
t601 = m(6) * t612 + t713 * mrSges(6,1) - t629 * mrSges(6,3) - t669 * t645 + t714 * t659 + t736;
t591 = t723 * t596 + t728 * t601;
t670 = -t683 * mrSges(5,1) + t684 * mrSges(5,2);
t673 = -t719 * mrSges(5,2) + t683 * mrSges(5,3);
t589 = m(5) * t623 + t716 * mrSges(5,1) - t654 * mrSges(5,3) - t684 * t670 + t719 * t673 + t591;
t674 = t719 * mrSges(5,1) - t684 * mrSges(5,3);
t746 = t728 * t596 - t723 * t601;
t590 = m(5) * t624 - t716 * mrSges(5,2) + t653 * mrSges(5,3) + t683 * t670 - t719 * t674 + t746;
t583 = t729 * t589 + t724 * t590;
t687 = -t700 * mrSges(4,1) + t701 * mrSges(4,2);
t694 = -qJD(3) * mrSges(4,2) + t700 * mrSges(4,3);
t581 = m(4) * t657 + qJDD(3) * mrSges(4,1) - t691 * mrSges(4,3) + qJD(3) * t694 - t701 * t687 + t583;
t695 = qJD(3) * mrSges(4,1) - t701 * mrSges(4,3);
t747 = -t724 * t589 + t729 * t590;
t582 = m(4) * t658 - qJDD(3) * mrSges(4,2) + t690 * mrSges(4,3) - qJD(3) * t695 + t700 * t687 + t747;
t576 = t730 * t581 + t725 * t582;
t692 = -t720 * t702 + t751;
t738 = mrSges(3,3) * qJDD(1) + t732 * (-mrSges(3,1) * t721 + t758);
t574 = m(3) * t692 - t720 * t738 + t576;
t748 = -t725 * t581 + t730 * t582;
t575 = m(3) * t693 + t721 * t738 + t748;
t749 = -t720 * t574 + t721 * t575;
t568 = m(2) * t706 - t732 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t749;
t698 = -qJDD(1) * pkin(1) - t732 * qJ(2) + t744;
t597 = t727 * t605 + t722 * t606;
t737 = m(6) * t622 - t628 * mrSges(6,1) + t629 * mrSges(6,2) - t668 * t659 + t669 * t660 + t597;
t735 = m(5) * t650 - t653 * mrSges(5,1) + t654 * mrSges(5,2) - t683 * t673 + t684 * t674 + t737;
t734 = m(4) * t689 - t690 * mrSges(4,1) + t691 * mrSges(4,2) - t700 * t694 + t701 * t695 + t735;
t733 = -m(3) * t698 + mrSges(3,1) * t752 - t734 + (t717 * t732 + t757) * mrSges(3,3);
t593 = t733 + (mrSges(2,1) - t758) * qJDD(1) + m(2) * t705 - t732 * mrSges(2,2);
t756 = t726 * t568 + t731 * t593;
t569 = t721 * t574 + t720 * t575;
t741 = Ifges(3,5) * t720 + Ifges(3,6) * t721;
t755 = t732 * t741;
t750 = t731 * t568 - t726 * t593;
t743 = Ifges(3,1) * t720 + Ifges(3,4) * t721;
t742 = Ifges(3,4) * t720 + Ifges(3,2) * t721;
t682 = Ifges(4,1) * t701 + Ifges(4,4) * t700 + Ifges(4,5) * qJD(3);
t681 = Ifges(4,4) * t701 + Ifges(4,2) * t700 + Ifges(4,6) * qJD(3);
t680 = Ifges(4,5) * t701 + Ifges(4,6) * t700 + Ifges(4,3) * qJD(3);
t663 = Ifges(5,1) * t684 + Ifges(5,4) * t683 + Ifges(5,5) * t719;
t662 = Ifges(5,4) * t684 + Ifges(5,2) * t683 + Ifges(5,6) * t719;
t661 = Ifges(5,5) * t684 + Ifges(5,6) * t683 + Ifges(5,3) * t719;
t642 = Ifges(6,1) * t669 + Ifges(6,4) * t668 + Ifges(6,5) * t714;
t641 = Ifges(6,4) * t669 + Ifges(6,2) * t668 + Ifges(6,6) * t714;
t640 = Ifges(6,5) * t669 + Ifges(6,6) * t668 + Ifges(6,3) * t714;
t632 = Ifges(7,1) * t656 + Ifges(7,4) * t655 + Ifges(7,5) * t664;
t631 = Ifges(7,4) * t656 + Ifges(7,2) * t655 + Ifges(7,6) * t664;
t630 = Ifges(7,5) * t656 + Ifges(7,6) * t655 + Ifges(7,3) * t664;
t599 = mrSges(7,2) * t609 - mrSges(7,3) * t607 + Ifges(7,1) * t620 + Ifges(7,4) * t619 + Ifges(7,5) * t627 + t655 * t630 - t664 * t631;
t598 = -mrSges(7,1) * t609 + mrSges(7,3) * t608 + Ifges(7,4) * t620 + Ifges(7,2) * t619 + Ifges(7,6) * t627 - t656 * t630 + t664 * t632;
t585 = -mrSges(6,1) * t622 - mrSges(7,1) * t607 + mrSges(7,2) * t608 + mrSges(6,3) * t613 + Ifges(6,4) * t629 - Ifges(7,5) * t620 + Ifges(6,2) * t628 + Ifges(6,6) * t713 - Ifges(7,6) * t619 - Ifges(7,3) * t627 - pkin(5) * t597 - t656 * t631 + t655 * t632 - t669 * t640 + t714 * t642;
t584 = mrSges(6,2) * t622 - mrSges(6,3) * t612 + Ifges(6,1) * t629 + Ifges(6,4) * t628 + Ifges(6,5) * t713 - pkin(10) * t597 - t722 * t598 + t727 * t599 + t668 * t640 - t714 * t641;
t577 = mrSges(5,2) * t650 - mrSges(5,3) * t623 + Ifges(5,1) * t654 + Ifges(5,4) * t653 + Ifges(5,5) * t716 - pkin(9) * t591 + t728 * t584 - t723 * t585 + t683 * t661 - t719 * t662;
t570 = -mrSges(5,1) * t650 + mrSges(5,3) * t624 + Ifges(5,4) * t654 + Ifges(5,2) * t653 + Ifges(5,6) * t716 - pkin(4) * t737 + pkin(9) * t746 + t723 * t584 + t728 * t585 - t684 * t661 + t719 * t663;
t565 = mrSges(4,2) * t689 - mrSges(4,3) * t657 + Ifges(4,1) * t691 + Ifges(4,4) * t690 + Ifges(4,5) * qJDD(3) - pkin(8) * t583 - qJD(3) * t681 - t724 * t570 + t729 * t577 + t700 * t680;
t564 = -mrSges(4,1) * t689 + mrSges(4,3) * t658 + Ifges(4,4) * t691 + Ifges(4,2) * t690 + Ifges(4,6) * qJDD(3) - pkin(3) * t735 + pkin(8) * t747 + qJD(3) * t682 + t729 * t570 + t724 * t577 - t701 * t680;
t563 = -pkin(10) * t745 + (Ifges(2,6) - t741) * qJDD(1) + (-t720 * t742 + t721 * t743 + Ifges(2,5)) * t732 - pkin(5) * t736 - pkin(2) * t576 - Ifges(4,3) * qJDD(3) + mrSges(2,1) * g(3) + t700 * t682 + t683 * t663 - t684 * t662 - Ifges(6,3) * t713 - Ifges(5,3) * t716 - mrSges(5,1) * t623 - t701 * t681 + mrSges(2,3) * t706 - mrSges(4,1) * t657 + mrSges(4,2) * t658 - Ifges(5,6) * t653 - Ifges(5,5) * t654 - Ifges(4,6) * t690 - Ifges(4,5) * t691 - mrSges(3,1) * t692 + mrSges(3,2) * t693 - t669 * t641 + t668 * t642 - mrSges(6,1) * t612 - pkin(4) * t591 + mrSges(5,2) * t624 - Ifges(6,6) * t628 - Ifges(6,5) * t629 - pkin(1) * t569 - pkin(3) * t583 - t727 * t598 - t722 * t599 + mrSges(6,2) * t613;
t562 = mrSges(3,2) * t698 - mrSges(3,3) * t692 - pkin(7) * t576 + qJDD(1) * t743 - t725 * t564 + t730 * t565 + t721 * t755;
t561 = -mrSges(3,1) * t698 + mrSges(3,3) * t693 - pkin(2) * t734 + pkin(7) * t748 + qJDD(1) * t742 + t730 * t564 + t725 * t565 - t720 * t755;
t560 = -mrSges(2,2) * g(3) - mrSges(2,3) * t705 + Ifges(2,5) * qJDD(1) - t732 * Ifges(2,6) - qJ(2) * t569 - t720 * t561 + t721 * t562;
t1 = [-m(1) * g(1) + t750; -m(1) * g(2) + t756; (-m(1) - m(2)) * g(3) + t569; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t756 + t731 * t560 - t726 * t563; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t750 + t726 * t560 + t731 * t563; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t705 - mrSges(2,2) * t706 + t720 * t562 + t721 * t561 + pkin(1) * (-qJDD(1) * t758 + t733) + qJ(2) * t749;];
tauB  = t1;
