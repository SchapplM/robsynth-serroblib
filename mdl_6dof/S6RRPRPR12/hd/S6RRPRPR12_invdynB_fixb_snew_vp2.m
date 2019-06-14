% Calculate vector of inverse dynamics base forces with Newton-Euler for
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 16:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:15:08
% EndTime: 2019-05-06 16:15:26
% DurationCPUTime: 17.48s
% Computational Cost: add. (269273->383), mult. (618545->479), div. (0->0), fcn. (458677->12), ass. (0->161)
t840 = -2 * qJD(3);
t839 = Ifges(3,1) + Ifges(4,2);
t831 = Ifges(3,4) + Ifges(4,6);
t830 = Ifges(3,5) - Ifges(4,4);
t838 = Ifges(3,2) + Ifges(4,3);
t829 = Ifges(3,6) - Ifges(4,5);
t837 = Ifges(3,3) + Ifges(4,1);
t787 = cos(pkin(6));
t779 = qJD(1) * t787 + qJD(2);
t790 = sin(qJ(2));
t785 = sin(pkin(6));
t818 = qJD(1) * t785;
t812 = t790 * t818;
t836 = (pkin(2) * t779 + t840) * t812;
t791 = sin(qJ(1));
t795 = cos(qJ(1));
t774 = t791 * g(1) - g(2) * t795;
t796 = qJD(1) ^ 2;
t759 = pkin(8) * t785 * t796 + qJDD(1) * pkin(1) + t774;
t775 = -g(1) * t795 - g(2) * t791;
t815 = qJDD(1) * t785;
t760 = -pkin(1) * t796 + pkin(8) * t815 + t775;
t794 = cos(qJ(2));
t825 = t787 * t790;
t827 = t785 * t790;
t724 = -g(3) * t827 + t759 * t825 + t794 * t760;
t761 = (-pkin(2) * t794 - qJ(3) * t790) * t818;
t777 = t779 ^ 2;
t778 = qJDD(1) * t787 + qJDD(2);
t817 = qJD(1) * t794;
t811 = t785 * t817;
t706 = t777 * pkin(2) - t778 * qJ(3) - t761 * t811 + t779 * t840 - t724;
t835 = 2 * qJD(5);
t834 = -pkin(2) - pkin(9);
t833 = t787 * g(3);
t832 = mrSges(3,1) - mrSges(4,2);
t828 = t785 ^ 2 * t796;
t826 = t785 * t794;
t824 = t787 * t794;
t739 = -t785 * t759 - t833;
t755 = mrSges(3,1) * t779 - mrSges(3,3) * t812;
t756 = -mrSges(3,2) * t779 + mrSges(3,3) * t811;
t758 = mrSges(4,1) * t812 + mrSges(4,2) * t779;
t765 = (qJD(2) * t817 + qJDD(1) * t790) * t785;
t766 = -qJD(2) * t812 + t794 * t815;
t707 = -t766 * pkin(2) + (-t779 * t811 - t765) * qJ(3) + t739 + t836;
t757 = -mrSges(4,1) * t811 - mrSges(4,3) * t779;
t764 = pkin(3) * t812 - pkin(9) * t779;
t814 = t794 ^ 2 * t828;
t690 = -pkin(3) * t814 - t833 - t765 * qJ(3) + t834 * t766 + (-t759 + (-qJ(3) * t779 * t794 - t764 * t790) * qJD(1)) * t785 + t836;
t819 = g(3) * t826 + t790 * t760;
t804 = -t777 * qJ(3) + t761 * t812 + qJDD(3) + t819;
t693 = t765 * pkin(3) + t834 * t778 + (-pkin(3) * t779 * t818 - pkin(9) * t790 * t828 - t759 * t787) * t794 + t804;
t789 = sin(qJ(4));
t793 = cos(qJ(4));
t676 = -t789 * t690 + t793 * t693;
t748 = -t779 * t789 - t793 * t811;
t722 = qJD(4) * t748 - t766 * t789 + t778 * t793;
t749 = t779 * t793 - t789 * t811;
t754 = qJDD(4) + t765;
t770 = qJD(4) + t812;
t673 = (t748 * t770 - t722) * qJ(5) + (t748 * t749 + t754) * pkin(4) + t676;
t677 = t793 * t690 + t789 * t693;
t721 = -qJD(4) * t749 - t766 * t793 - t778 * t789;
t731 = pkin(4) * t770 - qJ(5) * t749;
t747 = t748 ^ 2;
t675 = -pkin(4) * t747 + qJ(5) * t721 - t731 * t770 + t677;
t784 = sin(pkin(11));
t786 = cos(pkin(11));
t727 = t748 * t786 - t749 * t784;
t670 = t784 * t673 + t786 * t675 + t727 * t835;
t699 = t721 * t786 - t722 * t784;
t728 = t748 * t784 + t749 * t786;
t708 = -mrSges(6,1) * t727 + mrSges(6,2) * t728;
t714 = mrSges(6,1) * t770 - mrSges(6,3) * t728;
t709 = -pkin(5) * t727 - pkin(10) * t728;
t768 = t770 ^ 2;
t668 = -pkin(5) * t768 + pkin(10) * t754 + t709 * t727 + t670;
t689 = t766 * pkin(3) - pkin(9) * t814 + t779 * t764 - t706;
t679 = -t721 * pkin(4) - t747 * qJ(5) + t749 * t731 + qJDD(5) + t689;
t700 = t721 * t784 + t722 * t786;
t671 = (-t727 * t770 - t700) * pkin(10) + (t728 * t770 - t699) * pkin(5) + t679;
t788 = sin(qJ(6));
t792 = cos(qJ(6));
t665 = -t668 * t788 + t671 * t792;
t711 = -t728 * t788 + t770 * t792;
t682 = qJD(6) * t711 + t700 * t792 + t754 * t788;
t712 = t728 * t792 + t770 * t788;
t694 = -mrSges(7,1) * t711 + mrSges(7,2) * t712;
t726 = qJD(6) - t727;
t695 = -mrSges(7,2) * t726 + mrSges(7,3) * t711;
t698 = qJDD(6) - t699;
t663 = m(7) * t665 + mrSges(7,1) * t698 - mrSges(7,3) * t682 - t694 * t712 + t695 * t726;
t666 = t668 * t792 + t671 * t788;
t681 = -qJD(6) * t712 - t700 * t788 + t754 * t792;
t696 = mrSges(7,1) * t726 - mrSges(7,3) * t712;
t664 = m(7) * t666 - mrSges(7,2) * t698 + mrSges(7,3) * t681 + t694 * t711 - t696 * t726;
t807 = -t663 * t788 + t792 * t664;
t654 = m(6) * t670 - mrSges(6,2) * t754 + mrSges(6,3) * t699 + t708 * t727 - t714 * t770 + t807;
t806 = -t786 * t673 + t784 * t675;
t669 = -0.2e1 * qJD(5) * t728 - t806;
t713 = -mrSges(6,2) * t770 + mrSges(6,3) * t727;
t667 = -t754 * pkin(5) - t768 * pkin(10) + (t835 + t709) * t728 + t806;
t801 = -m(7) * t667 + t681 * mrSges(7,1) - mrSges(7,2) * t682 + t711 * t695 - t696 * t712;
t659 = m(6) * t669 + mrSges(6,1) * t754 - mrSges(6,3) * t700 - t708 * t728 + t713 * t770 + t801;
t647 = t784 * t654 + t786 * t659;
t729 = -mrSges(5,1) * t748 + mrSges(5,2) * t749;
t730 = -mrSges(5,2) * t770 + mrSges(5,3) * t748;
t645 = m(5) * t676 + mrSges(5,1) * t754 - mrSges(5,3) * t722 - t729 * t749 + t730 * t770 + t647;
t732 = mrSges(5,1) * t770 - mrSges(5,3) * t749;
t808 = t786 * t654 - t659 * t784;
t646 = m(5) * t677 - mrSges(5,2) * t754 + mrSges(5,3) * t721 + t729 * t748 - t732 * t770 + t808;
t809 = -t789 * t645 + t793 * t646;
t805 = m(4) * t707 - t765 * mrSges(4,3) + t757 * t811 + t809;
t638 = m(3) * t739 + t765 * mrSges(3,2) - t832 * t766 + (-t756 * t794 + (t755 - t758) * t790) * t818 + t805;
t813 = t759 * t824;
t723 = t813 - t819;
t762 = (mrSges(4,2) * t794 - mrSges(4,3) * t790) * t818;
t763 = (-mrSges(3,1) * t794 + mrSges(3,2) * t790) * t818;
t641 = t793 * t645 + t789 * t646;
t710 = -t778 * pkin(2) + t804 - t813;
t802 = -m(4) * t710 - t765 * mrSges(4,1) - t641;
t639 = m(3) * t723 - t765 * mrSges(3,3) + (t756 - t757) * t779 + t832 * t778 + (-t762 - t763) * t812 + t802;
t655 = t792 * t663 + t788 * t664;
t799 = m(6) * t679 - t699 * mrSges(6,1) + t700 * mrSges(6,2) - t727 * t713 + t728 * t714 + t655;
t798 = -m(5) * t689 + t721 * mrSges(5,1) - t722 * mrSges(5,2) + t748 * t730 - t749 * t732 - t799;
t797 = -m(4) * t706 + t778 * mrSges(4,3) + t779 * t758 + t762 * t811 - t798;
t651 = t797 + t763 * t811 + (mrSges(3,3) + mrSges(4,1)) * t766 - t778 * mrSges(3,2) - t779 * t755 + m(3) * t724;
t628 = -t638 * t785 + t639 * t824 + t651 * t825;
t626 = m(2) * t774 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t796 + t628;
t634 = -t639 * t790 + t794 * t651;
t633 = m(2) * t775 - mrSges(2,1) * t796 - qJDD(1) * mrSges(2,2) + t634;
t823 = t795 * t626 + t791 * t633;
t822 = (t830 * t790 + t829 * t794) * t818 + t837 * t779;
t821 = (-t831 * t790 - t838 * t794) * t818 - t829 * t779;
t820 = (t839 * t790 + t831 * t794) * t818 + t830 * t779;
t627 = t787 * t638 + t639 * t826 + t651 * t827;
t810 = -t626 * t791 + t795 * t633;
t683 = Ifges(7,5) * t712 + Ifges(7,6) * t711 + Ifges(7,3) * t726;
t685 = Ifges(7,1) * t712 + Ifges(7,4) * t711 + Ifges(7,5) * t726;
t656 = -mrSges(7,1) * t667 + mrSges(7,3) * t666 + Ifges(7,4) * t682 + Ifges(7,2) * t681 + Ifges(7,6) * t698 - t683 * t712 + t685 * t726;
t684 = Ifges(7,4) * t712 + Ifges(7,2) * t711 + Ifges(7,6) * t726;
t657 = mrSges(7,2) * t667 - mrSges(7,3) * t665 + Ifges(7,1) * t682 + Ifges(7,4) * t681 + Ifges(7,5) * t698 + t683 * t711 - t684 * t726;
t701 = Ifges(6,5) * t728 + Ifges(6,6) * t727 + Ifges(6,3) * t770;
t702 = Ifges(6,4) * t728 + Ifges(6,2) * t727 + Ifges(6,6) * t770;
t642 = mrSges(6,2) * t679 - mrSges(6,3) * t669 + Ifges(6,1) * t700 + Ifges(6,4) * t699 + Ifges(6,5) * t754 - pkin(10) * t655 - t656 * t788 + t657 * t792 + t701 * t727 - t702 * t770;
t703 = Ifges(6,1) * t728 + Ifges(6,4) * t727 + Ifges(6,5) * t770;
t643 = -mrSges(6,1) * t679 - mrSges(7,1) * t665 + mrSges(7,2) * t666 + mrSges(6,3) * t670 + Ifges(6,4) * t700 - Ifges(7,5) * t682 + Ifges(6,2) * t699 + Ifges(6,6) * t754 - Ifges(7,6) * t681 - Ifges(7,3) * t698 - pkin(5) * t655 - t684 * t712 + t685 * t711 - t701 * t728 + t703 * t770;
t715 = Ifges(5,5) * t749 + Ifges(5,6) * t748 + Ifges(5,3) * t770;
t717 = Ifges(5,1) * t749 + Ifges(5,4) * t748 + Ifges(5,5) * t770;
t629 = -mrSges(5,1) * t689 + mrSges(5,3) * t677 + Ifges(5,4) * t722 + Ifges(5,2) * t721 + Ifges(5,6) * t754 - pkin(4) * t799 + qJ(5) * t808 + t784 * t642 + t786 * t643 - t749 * t715 + t770 * t717;
t716 = Ifges(5,4) * t749 + Ifges(5,2) * t748 + Ifges(5,6) * t770;
t630 = mrSges(5,2) * t689 - mrSges(5,3) * t676 + Ifges(5,1) * t722 + Ifges(5,4) * t721 + Ifges(5,5) * t754 - qJ(5) * t647 + t642 * t786 - t643 * t784 + t715 * t748 - t716 * t770;
t640 = t766 * mrSges(4,2) - t758 * t812 + t805;
t623 = -mrSges(3,1) * t739 - mrSges(4,1) * t706 + mrSges(4,2) * t707 + mrSges(3,3) * t724 - pkin(2) * t640 - pkin(3) * t798 - pkin(9) * t809 - t793 * t629 - t789 * t630 + t831 * t765 + t838 * t766 + t829 * t778 + t820 * t779 - t822 * t812;
t624 = t822 * t811 + t821 * t779 + pkin(3) * t641 - qJ(3) * t640 + t839 * t765 + pkin(5) * t801 + pkin(10) * t807 + pkin(4) * t647 + (Ifges(6,3) + Ifges(5,3)) * t754 + t830 * t778 + t831 * t766 + t788 * t657 + t792 * t656 - t748 * t717 + t749 * t716 + mrSges(3,2) * t739 - t727 * t703 + t728 * t702 + Ifges(5,6) * t721 + Ifges(5,5) * t722 - mrSges(3,3) * t723 - mrSges(4,3) * t707 + mrSges(4,1) * t710 + Ifges(6,6) * t699 + Ifges(6,5) * t700 + mrSges(6,1) * t669 - mrSges(6,2) * t670 + mrSges(5,1) * t676 - mrSges(5,2) * t677;
t803 = pkin(8) * t634 + t623 * t794 + t624 * t790;
t622 = mrSges(3,1) * t723 - mrSges(3,2) * t724 + mrSges(4,2) * t710 - mrSges(4,3) * t706 + t793 * t630 - t789 * t629 - pkin(9) * t641 + pkin(2) * (-t779 * t757 + t802) + qJ(3) * t797 + (-mrSges(4,2) * pkin(2) + t837) * t778 + (mrSges(4,1) * qJ(3) + t829) * t766 + t830 * t765 + (-t820 * t794 + (-pkin(2) * t762 - t821) * t790) * t818;
t621 = -mrSges(2,2) * g(3) - mrSges(2,3) * t774 + Ifges(2,5) * qJDD(1) - t796 * Ifges(2,6) - t790 * t623 + t794 * t624 + (-t627 * t785 - t628 * t787) * pkin(8);
t620 = mrSges(2,1) * g(3) + mrSges(2,3) * t775 + t796 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t627 - t785 * t622 + t787 * t803;
t1 = [-m(1) * g(1) + t810; -m(1) * g(2) + t823; (-m(1) - m(2)) * g(3) + t627; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t823 - t791 * t620 + t795 * t621; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t810 + t795 * t620 + t791 * t621; -mrSges(1,1) * g(2) + mrSges(2,1) * t774 + mrSges(1,2) * g(1) - mrSges(2,2) * t775 + Ifges(2,3) * qJDD(1) + pkin(1) * t628 + t787 * t622 + t785 * t803;];
tauB  = t1;
