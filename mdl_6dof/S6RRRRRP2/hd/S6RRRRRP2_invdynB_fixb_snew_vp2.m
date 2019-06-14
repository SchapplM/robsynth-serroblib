% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:30:41
% EndTime: 2019-05-08 04:31:01
% DurationCPUTime: 14.07s
% Computational Cost: add. (212853->364), mult. (460992->447), div. (0->0), fcn. (342628->10), ass. (0->141)
t755 = Ifges(6,1) + Ifges(7,1);
t748 = Ifges(6,4) - Ifges(7,5);
t754 = -Ifges(6,5) - Ifges(7,4);
t753 = Ifges(6,2) + Ifges(7,3);
t746 = Ifges(6,6) - Ifges(7,6);
t752 = -Ifges(6,3) - Ifges(7,2);
t751 = cos(qJ(5));
t724 = qJD(1) ^ 2;
t750 = pkin(2) * t724;
t749 = -mrSges(6,3) - mrSges(7,2);
t719 = sin(qJ(1));
t723 = cos(qJ(1));
t705 = -t723 * g(1) - t719 * g(2);
t694 = -t724 * pkin(1) + qJDD(1) * pkin(7) + t705;
t718 = sin(qJ(2));
t745 = t718 * t694;
t722 = cos(qJ(2));
t737 = qJD(1) * qJD(2);
t699 = t718 * qJDD(1) + t722 * t737;
t661 = qJDD(2) * pkin(2) - t699 * pkin(8) - t745 + (pkin(8) * t737 + t718 * t750 - g(3)) * t722;
t682 = -t718 * g(3) + t722 * t694;
t700 = t722 * qJDD(1) - t718 * t737;
t739 = qJD(1) * t718;
t703 = qJD(2) * pkin(2) - pkin(8) * t739;
t714 = t722 ^ 2;
t662 = t700 * pkin(8) - qJD(2) * t703 - t714 * t750 + t682;
t717 = sin(qJ(3));
t721 = cos(qJ(3));
t642 = t721 * t661 - t717 * t662;
t691 = (-t717 * t718 + t721 * t722) * qJD(1);
t667 = t691 * qJD(3) + t721 * t699 + t717 * t700;
t692 = (t717 * t722 + t718 * t721) * qJD(1);
t712 = qJDD(2) + qJDD(3);
t713 = qJD(2) + qJD(3);
t614 = (t691 * t713 - t667) * pkin(9) + (t691 * t692 + t712) * pkin(3) + t642;
t643 = t717 * t661 + t721 * t662;
t666 = -t692 * qJD(3) - t717 * t699 + t721 * t700;
t685 = t713 * pkin(3) - t692 * pkin(9);
t687 = t691 ^ 2;
t619 = -t687 * pkin(3) + t666 * pkin(9) - t713 * t685 + t643;
t716 = sin(qJ(4));
t720 = cos(qJ(4));
t612 = t716 * t614 + t720 * t619;
t679 = t716 * t691 + t720 * t692;
t637 = -t679 * qJD(4) + t720 * t666 - t716 * t667;
t678 = t720 * t691 - t716 * t692;
t656 = -t678 * mrSges(5,1) + t679 * mrSges(5,2);
t710 = qJD(4) + t713;
t670 = t710 * mrSges(5,1) - t679 * mrSges(5,3);
t709 = qJDD(4) + t712;
t657 = -t678 * pkin(4) - t679 * pkin(10);
t708 = t710 ^ 2;
t608 = -t708 * pkin(4) + t709 * pkin(10) + t678 * t657 + t612;
t704 = t719 * g(1) - t723 * g(2);
t729 = -qJDD(1) * pkin(1) - t704;
t668 = -t700 * pkin(2) + t703 * t739 + (-pkin(8) * t714 - pkin(7)) * t724 + t729;
t623 = -t666 * pkin(3) - t687 * pkin(9) + t692 * t685 + t668;
t638 = t678 * qJD(4) + t716 * t666 + t720 * t667;
t610 = (-t678 * t710 - t638) * pkin(10) + (t679 * t710 - t637) * pkin(4) + t623;
t715 = sin(qJ(5));
t605 = t751 * t608 + t715 * t610;
t664 = t751 * t679 + t715 * t710;
t620 = t664 * qJD(5) + t715 * t638 - t751 * t709;
t634 = qJDD(5) - t637;
t675 = qJD(5) - t678;
t649 = t675 * mrSges(6,1) - t664 * mrSges(6,3);
t663 = t715 * t679 - t751 * t710;
t644 = t663 * pkin(5) - t664 * qJ(6);
t674 = t675 ^ 2;
t601 = -t674 * pkin(5) + t634 * qJ(6) + 0.2e1 * qJD(6) * t675 - t663 * t644 + t605;
t650 = -t675 * mrSges(7,1) + t664 * mrSges(7,2);
t736 = m(7) * t601 + t634 * mrSges(7,3) + t675 * t650;
t645 = t663 * mrSges(7,1) - t664 * mrSges(7,3);
t740 = -t663 * mrSges(6,1) - t664 * mrSges(6,2) - t645;
t596 = m(6) * t605 - t634 * mrSges(6,2) + t749 * t620 - t675 * t649 + t740 * t663 + t736;
t604 = -t715 * t608 + t751 * t610;
t621 = -t663 * qJD(5) + t751 * t638 + t715 * t709;
t648 = -t675 * mrSges(6,2) - t663 * mrSges(6,3);
t602 = -t634 * pkin(5) - t674 * qJ(6) + t664 * t644 + qJDD(6) - t604;
t647 = -t663 * mrSges(7,2) + t675 * mrSges(7,3);
t730 = -m(7) * t602 + t634 * mrSges(7,1) + t675 * t647;
t598 = m(6) * t604 + t634 * mrSges(6,1) + t749 * t621 + t675 * t648 + t740 * t664 + t730;
t731 = t751 * t596 - t715 * t598;
t588 = m(5) * t612 - t709 * mrSges(5,2) + t637 * mrSges(5,3) + t678 * t656 - t710 * t670 + t731;
t611 = t720 * t614 - t716 * t619;
t669 = -t710 * mrSges(5,2) + t678 * mrSges(5,3);
t607 = -t709 * pkin(4) - t708 * pkin(10) + t679 * t657 - t611;
t603 = -0.2e1 * qJD(6) * t664 + (t663 * t675 - t621) * qJ(6) + (t664 * t675 + t620) * pkin(5) + t607;
t599 = m(7) * t603 + t620 * mrSges(7,1) - t621 * mrSges(7,3) + t663 * t647 - t664 * t650;
t726 = -m(6) * t607 - t620 * mrSges(6,1) - t621 * mrSges(6,2) - t663 * t648 - t664 * t649 - t599;
t593 = m(5) * t611 + t709 * mrSges(5,1) - t638 * mrSges(5,3) - t679 * t656 + t710 * t669 + t726;
t583 = t716 * t588 + t720 * t593;
t680 = -t691 * mrSges(4,1) + t692 * mrSges(4,2);
t683 = -t713 * mrSges(4,2) + t691 * mrSges(4,3);
t581 = m(4) * t642 + t712 * mrSges(4,1) - t667 * mrSges(4,3) - t692 * t680 + t713 * t683 + t583;
t684 = t713 * mrSges(4,1) - t692 * mrSges(4,3);
t732 = t720 * t588 - t716 * t593;
t582 = m(4) * t643 - t712 * mrSges(4,2) + t666 * mrSges(4,3) + t691 * t680 - t713 * t684 + t732;
t575 = t721 * t581 + t717 * t582;
t681 = -t722 * g(3) - t745;
t698 = (-mrSges(3,1) * t722 + mrSges(3,2) * t718) * qJD(1);
t738 = qJD(1) * t722;
t702 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t738;
t573 = m(3) * t681 + qJDD(2) * mrSges(3,1) - t699 * mrSges(3,3) + qJD(2) * t702 - t698 * t739 + t575;
t701 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t739;
t733 = -t717 * t581 + t721 * t582;
t574 = m(3) * t682 - qJDD(2) * mrSges(3,2) + t700 * mrSges(3,3) - qJD(2) * t701 + t698 * t738 + t733;
t734 = -t718 * t573 + t722 * t574;
t568 = m(2) * t705 - t724 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t734;
t693 = -t724 * pkin(7) + t729;
t591 = t715 * t596 + t751 * t598;
t728 = m(5) * t623 - t637 * mrSges(5,1) + t638 * mrSges(5,2) - t678 * t669 + t679 * t670 + t591;
t727 = m(4) * t668 - t666 * mrSges(4,1) + t667 * mrSges(4,2) - t691 * t683 + t692 * t684 + t728;
t725 = -m(3) * t693 + t700 * mrSges(3,1) - t699 * mrSges(3,2) - t701 * t739 + t702 * t738 - t727;
t585 = m(2) * t704 + qJDD(1) * mrSges(2,1) - t724 * mrSges(2,2) + t725;
t744 = t719 * t568 + t723 * t585;
t569 = t722 * t573 + t718 * t574;
t743 = t663 * t753 - t664 * t748 - t675 * t746;
t742 = t663 * t746 + t664 * t754 + t675 * t752;
t741 = -t748 * t663 + t664 * t755 - t754 * t675;
t735 = t723 * t568 - t719 * t585;
t690 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t718 + Ifges(3,4) * t722) * qJD(1);
t689 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t718 + Ifges(3,2) * t722) * qJD(1);
t688 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t718 + Ifges(3,6) * t722) * qJD(1);
t673 = Ifges(4,1) * t692 + Ifges(4,4) * t691 + Ifges(4,5) * t713;
t672 = Ifges(4,4) * t692 + Ifges(4,2) * t691 + Ifges(4,6) * t713;
t671 = Ifges(4,5) * t692 + Ifges(4,6) * t691 + Ifges(4,3) * t713;
t653 = Ifges(5,1) * t679 + Ifges(5,4) * t678 + Ifges(5,5) * t710;
t652 = Ifges(5,4) * t679 + Ifges(5,2) * t678 + Ifges(5,6) * t710;
t651 = Ifges(5,5) * t679 + Ifges(5,6) * t678 + Ifges(5,3) * t710;
t590 = mrSges(6,2) * t607 + mrSges(7,2) * t602 - mrSges(6,3) * t604 - mrSges(7,3) * t603 - qJ(6) * t599 - t748 * t620 + t621 * t755 - t634 * t754 + t742 * t663 + t743 * t675;
t589 = -mrSges(6,1) * t607 - mrSges(7,1) * t603 + mrSges(7,2) * t601 + mrSges(6,3) * t605 - pkin(5) * t599 - t620 * t753 + t748 * t621 + t746 * t634 + t742 * t664 + t741 * t675;
t577 = Ifges(5,4) * t638 + Ifges(5,2) * t637 + Ifges(5,6) * t709 - t679 * t651 + t710 * t653 - mrSges(5,1) * t623 + mrSges(5,3) * t612 - mrSges(6,1) * t604 + mrSges(6,2) * t605 + mrSges(7,1) * t602 - mrSges(7,3) * t601 - pkin(5) * t730 - qJ(6) * t736 - pkin(4) * t591 + (pkin(5) * t645 + t743) * t664 + (qJ(6) * t645 - t741) * t663 + t752 * t634 + (pkin(5) * mrSges(7,2) + t754) * t621 + (qJ(6) * mrSges(7,2) + t746) * t620;
t576 = mrSges(5,2) * t623 - mrSges(5,3) * t611 + Ifges(5,1) * t638 + Ifges(5,4) * t637 + Ifges(5,5) * t709 - pkin(10) * t591 - t715 * t589 + t751 * t590 + t678 * t651 - t710 * t652;
t565 = mrSges(4,2) * t668 - mrSges(4,3) * t642 + Ifges(4,1) * t667 + Ifges(4,4) * t666 + Ifges(4,5) * t712 - pkin(9) * t583 + t720 * t576 - t716 * t577 + t691 * t671 - t713 * t672;
t564 = -mrSges(4,1) * t668 + mrSges(4,3) * t643 + Ifges(4,4) * t667 + Ifges(4,2) * t666 + Ifges(4,6) * t712 - pkin(3) * t728 + pkin(9) * t732 + t716 * t576 + t720 * t577 - t692 * t671 + t713 * t673;
t563 = -pkin(4) * t726 - pkin(10) * t731 - t751 * t589 - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) - pkin(1) * t569 + Ifges(2,6) * qJDD(1) + (-t718 * t689 + t722 * t690) * qJD(1) - pkin(3) * t583 + t724 * Ifges(2,5) - t715 * t590 + mrSges(2,3) * t705 - Ifges(5,3) * t709 - Ifges(4,3) * t712 - t692 * t672 - Ifges(3,5) * t699 - Ifges(3,6) * t700 + mrSges(3,2) * t682 + t691 * t673 + t678 * t653 - t679 * t652 - mrSges(3,1) * t681 - Ifges(4,6) * t666 - Ifges(4,5) * t667 + mrSges(4,2) * t643 - Ifges(5,6) * t637 - Ifges(5,5) * t638 - mrSges(4,1) * t642 + mrSges(5,2) * t612 - mrSges(5,1) * t611 - pkin(2) * t575;
t562 = mrSges(3,2) * t693 - mrSges(3,3) * t681 + Ifges(3,1) * t699 + Ifges(3,4) * t700 + Ifges(3,5) * qJDD(2) - pkin(8) * t575 - qJD(2) * t689 - t717 * t564 + t721 * t565 + t688 * t738;
t561 = -mrSges(3,1) * t693 + mrSges(3,3) * t682 + Ifges(3,4) * t699 + Ifges(3,2) * t700 + Ifges(3,6) * qJDD(2) - pkin(2) * t727 + pkin(8) * t733 + qJD(2) * t690 + t721 * t564 + t717 * t565 - t688 * t739;
t560 = -mrSges(2,2) * g(3) - mrSges(2,3) * t704 + Ifges(2,5) * qJDD(1) - t724 * Ifges(2,6) - pkin(7) * t569 - t718 * t561 + t722 * t562;
t1 = [-m(1) * g(1) + t735; -m(1) * g(2) + t744; (-m(1) - m(2)) * g(3) + t569; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t744 + t723 * t560 - t719 * t563; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t735 + t719 * t560 + t723 * t563; -mrSges(1,1) * g(2) + mrSges(2,1) * t704 + mrSges(1,2) * g(1) - mrSges(2,2) * t705 + Ifges(2,3) * qJDD(1) + pkin(1) * t725 + pkin(7) * t734 + t722 * t561 + t718 * t562;];
tauB  = t1;
