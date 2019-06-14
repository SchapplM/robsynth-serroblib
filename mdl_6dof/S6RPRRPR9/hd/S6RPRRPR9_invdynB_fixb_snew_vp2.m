% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 23:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:22:52
% EndTime: 2019-05-05 23:24:38
% DurationCPUTime: 108.69s
% Computational Cost: add. (1668535->403), mult. (5237182->551), div. (0->0), fcn. (4506346->16), ass. (0->184)
t806 = sin(pkin(12));
t808 = sin(pkin(6));
t810 = cos(pkin(12));
t812 = cos(pkin(6));
t815 = sin(qJ(3));
t811 = cos(pkin(7));
t819 = cos(qJ(3));
t854 = t811 * t819;
t807 = sin(pkin(7));
t859 = t807 * t819;
t823 = (-t806 * t815 + t810 * t854) * t808 + t812 * t859;
t780 = t823 * qJD(1);
t855 = t811 * t815;
t860 = t807 * t815;
t825 = t812 * t860 + (t806 * t819 + t810 * t855) * t808;
t781 = t825 * qJD(1);
t769 = -t781 * qJD(3) + t823 * qJDD(1);
t857 = t808 * t811;
t792 = (t807 * t812 + t810 * t857) * qJD(1) * pkin(9);
t816 = sin(qJ(1));
t820 = cos(qJ(1));
t803 = -t820 * g(1) - t816 * g(2);
t821 = qJD(1) ^ 2;
t863 = qJ(2) * t808;
t796 = -t821 * pkin(1) + qJDD(1) * t863 + t803;
t866 = pkin(9) * t806;
t837 = -pkin(2) * t810 - t807 * t866;
t852 = qJD(1) * t808;
t864 = pkin(9) * qJDD(1);
t832 = qJD(1) * t837 * t852 + t811 * t864;
t802 = t816 * g(1) - t820 * g(2);
t795 = qJDD(1) * pkin(1) + t821 * t863 + t802;
t848 = qJD(2) * t852;
t856 = t810 * t812;
t858 = t808 * t810;
t838 = -g(3) * t858 + t795 * t856 - 0.2e1 * t806 * t848;
t749 = (pkin(2) * qJDD(1) + qJD(1) * t792) * t812 + (-t832 * t808 - t796) * t806 + t838;
t797 = (pkin(2) * t812 - t857 * t866) * qJD(1);
t861 = t806 * t812;
t849 = t795 * t861 + (t796 + 0.2e1 * t848) * t810;
t750 = (-qJD(1) * t797 + t807 * t864) * t812 + (-g(3) * t806 + t832 * t810) * t808 + t849;
t847 = -t812 * g(3) + qJDD(2);
t759 = (-t795 + t837 * qJDD(1) + (-t792 * t810 + t797 * t806) * qJD(1)) * t808 + t847;
t718 = -t815 * t750 + (t749 * t811 + t759 * t807) * t819;
t867 = 2 * qJD(5);
t865 = Ifges(3,3) * t812;
t862 = t806 * t808;
t719 = t749 * t855 + t819 * t750 + t759 * t860;
t767 = -t780 * mrSges(4,1) + t781 * mrSges(4,2);
t833 = -t807 * t858 + t811 * t812;
t793 = t833 * qJD(1) + qJD(3);
t777 = t793 * mrSges(4,1) - t781 * mrSges(4,3);
t790 = t833 * qJDD(1) + qJDD(3);
t768 = -t780 * pkin(3) - t781 * pkin(10);
t789 = t793 ^ 2;
t710 = -t789 * pkin(3) + t790 * pkin(10) + t780 * t768 + t719;
t733 = -t807 * t749 + t811 * t759;
t770 = t780 * qJD(3) + t825 * qJDD(1);
t713 = (-t780 * t793 - t770) * pkin(10) + (t781 * t793 - t769) * pkin(3) + t733;
t814 = sin(qJ(4));
t818 = cos(qJ(4));
t702 = -t814 * t710 + t818 * t713;
t774 = -t814 * t781 + t818 * t793;
t745 = t774 * qJD(4) + t818 * t770 + t814 * t790;
t766 = qJDD(4) - t769;
t775 = t818 * t781 + t814 * t793;
t779 = qJD(4) - t780;
t699 = (t774 * t779 - t745) * qJ(5) + (t774 * t775 + t766) * pkin(4) + t702;
t703 = t818 * t710 + t814 * t713;
t744 = -t775 * qJD(4) - t814 * t770 + t818 * t790;
t761 = t779 * pkin(4) - t775 * qJ(5);
t773 = t774 ^ 2;
t701 = -t773 * pkin(4) + t744 * qJ(5) - t779 * t761 + t703;
t805 = sin(pkin(13));
t809 = cos(pkin(13));
t753 = t809 * t774 - t805 * t775;
t696 = t805 * t699 + t809 * t701 + t753 * t867;
t723 = t809 * t744 - t805 * t745;
t754 = t805 * t774 + t809 * t775;
t731 = -t753 * mrSges(6,1) + t754 * mrSges(6,2);
t737 = t779 * mrSges(6,1) - t754 * mrSges(6,3);
t732 = -t753 * pkin(5) - t754 * pkin(11);
t778 = t779 ^ 2;
t694 = -t778 * pkin(5) + t766 * pkin(11) + t753 * t732 + t696;
t709 = -t790 * pkin(3) - t789 * pkin(10) + t781 * t768 - t718;
t704 = -t744 * pkin(4) - t773 * qJ(5) + t775 * t761 + qJDD(5) + t709;
t724 = t805 * t744 + t809 * t745;
t697 = (-t753 * t779 - t724) * pkin(11) + (t754 * t779 - t723) * pkin(5) + t704;
t813 = sin(qJ(6));
t817 = cos(qJ(6));
t691 = -t813 * t694 + t817 * t697;
t734 = -t813 * t754 + t817 * t779;
t707 = t734 * qJD(6) + t817 * t724 + t813 * t766;
t735 = t817 * t754 + t813 * t779;
t720 = -t734 * mrSges(7,1) + t735 * mrSges(7,2);
t722 = qJDD(6) - t723;
t752 = qJD(6) - t753;
t725 = -t752 * mrSges(7,2) + t734 * mrSges(7,3);
t689 = m(7) * t691 + t722 * mrSges(7,1) - t707 * mrSges(7,3) - t735 * t720 + t752 * t725;
t692 = t817 * t694 + t813 * t697;
t706 = -t735 * qJD(6) - t813 * t724 + t817 * t766;
t726 = t752 * mrSges(7,1) - t735 * mrSges(7,3);
t690 = m(7) * t692 - t722 * mrSges(7,2) + t706 * mrSges(7,3) + t734 * t720 - t752 * t726;
t843 = -t813 * t689 + t817 * t690;
t679 = m(6) * t696 - t766 * mrSges(6,2) + t723 * mrSges(6,3) + t753 * t731 - t779 * t737 + t843;
t840 = -t809 * t699 + t805 * t701;
t695 = -0.2e1 * qJD(5) * t754 - t840;
t736 = -t779 * mrSges(6,2) + t753 * mrSges(6,3);
t693 = -t766 * pkin(5) - t778 * pkin(11) + (t867 + t732) * t754 + t840;
t829 = -m(7) * t693 + t706 * mrSges(7,1) - t707 * mrSges(7,2) + t734 * t725 - t735 * t726;
t685 = m(6) * t695 + t766 * mrSges(6,1) - t724 * mrSges(6,3) - t754 * t731 + t779 * t736 + t829;
t674 = t805 * t679 + t809 * t685;
t755 = -t774 * mrSges(5,1) + t775 * mrSges(5,2);
t760 = -t779 * mrSges(5,2) + t774 * mrSges(5,3);
t672 = m(5) * t702 + t766 * mrSges(5,1) - t745 * mrSges(5,3) - t775 * t755 + t779 * t760 + t674;
t762 = t779 * mrSges(5,1) - t775 * mrSges(5,3);
t844 = t809 * t679 - t805 * t685;
t673 = m(5) * t703 - t766 * mrSges(5,2) + t744 * mrSges(5,3) + t774 * t755 - t779 * t762 + t844;
t845 = -t814 * t672 + t818 * t673;
t663 = m(4) * t719 - t790 * mrSges(4,2) + t769 * mrSges(4,3) + t780 * t767 - t793 * t777 + t845;
t666 = t818 * t672 + t814 * t673;
t776 = -t793 * mrSges(4,2) + t780 * mrSges(4,3);
t665 = m(4) * t733 - t769 * mrSges(4,1) + t770 * mrSges(4,2) - t780 * t776 + t781 * t777 + t666;
t681 = t817 * t689 + t813 * t690;
t826 = m(6) * t704 - t723 * mrSges(6,1) + t724 * mrSges(6,2) - t753 * t736 + t754 * t737 + t681;
t822 = -m(5) * t709 + t744 * mrSges(5,1) - t745 * mrSges(5,2) + t774 * t760 - t775 * t762 - t826;
t680 = m(4) * t718 + t790 * mrSges(4,1) - t770 * mrSges(4,3) - t781 * t767 + t793 * t776 + t822;
t652 = t663 * t855 - t807 * t665 + t680 * t854;
t771 = -t806 * t796 + t838;
t842 = -mrSges(3,1) * t810 + mrSges(3,2) * t806;
t794 = t842 * t852;
t835 = -mrSges(3,2) * t812 + mrSges(3,3) * t858;
t799 = t835 * qJD(1);
t836 = mrSges(3,1) * t812 - mrSges(3,3) * t862;
t648 = m(3) * t771 + t836 * qJDD(1) + (-t794 * t862 + t799 * t812) * qJD(1) + t652;
t651 = t663 * t860 + t811 * t665 + t680 * t859;
t782 = -t808 * t795 + t847;
t798 = t836 * qJD(1);
t650 = m(3) * t782 + (t842 * qJDD(1) + (t798 * t806 - t799 * t810) * qJD(1)) * t808 + t651;
t659 = t819 * t663 - t815 * t680;
t772 = -g(3) * t862 + t849;
t658 = m(3) * t772 + t835 * qJDD(1) + (t794 * t858 - t798 * t812) * qJD(1) + t659;
t638 = t648 * t856 - t808 * t650 + t658 * t861;
t636 = m(2) * t802 + qJDD(1) * mrSges(2,1) - t821 * mrSges(2,2) + t638;
t644 = -t806 * t648 + t810 * t658;
t643 = m(2) * t803 - t821 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t644;
t853 = t820 * t636 + t816 * t643;
t637 = t648 * t858 + t812 * t650 + t658 * t862;
t846 = -t816 * t636 + t820 * t643;
t841 = Ifges(3,5) * t806 + Ifges(3,6) * t810;
t714 = Ifges(7,5) * t735 + Ifges(7,6) * t734 + Ifges(7,3) * t752;
t716 = Ifges(7,1) * t735 + Ifges(7,4) * t734 + Ifges(7,5) * t752;
t682 = -mrSges(7,1) * t693 + mrSges(7,3) * t692 + Ifges(7,4) * t707 + Ifges(7,2) * t706 + Ifges(7,6) * t722 - t735 * t714 + t752 * t716;
t715 = Ifges(7,4) * t735 + Ifges(7,2) * t734 + Ifges(7,6) * t752;
t683 = mrSges(7,2) * t693 - mrSges(7,3) * t691 + Ifges(7,1) * t707 + Ifges(7,4) * t706 + Ifges(7,5) * t722 + t734 * t714 - t752 * t715;
t727 = Ifges(6,5) * t754 + Ifges(6,6) * t753 + Ifges(6,3) * t779;
t728 = Ifges(6,4) * t754 + Ifges(6,2) * t753 + Ifges(6,6) * t779;
t667 = mrSges(6,2) * t704 - mrSges(6,3) * t695 + Ifges(6,1) * t724 + Ifges(6,4) * t723 + Ifges(6,5) * t766 - pkin(11) * t681 - t813 * t682 + t817 * t683 + t753 * t727 - t779 * t728;
t729 = Ifges(6,1) * t754 + Ifges(6,4) * t753 + Ifges(6,5) * t779;
t668 = -mrSges(6,1) * t704 - mrSges(7,1) * t691 + mrSges(7,2) * t692 + mrSges(6,3) * t696 + Ifges(6,4) * t724 - Ifges(7,5) * t707 + Ifges(6,2) * t723 + Ifges(6,6) * t766 - Ifges(7,6) * t706 - Ifges(7,3) * t722 - pkin(5) * t681 - t735 * t715 + t734 * t716 - t754 * t727 + t779 * t729;
t738 = Ifges(5,5) * t775 + Ifges(5,6) * t774 + Ifges(5,3) * t779;
t740 = Ifges(5,1) * t775 + Ifges(5,4) * t774 + Ifges(5,5) * t779;
t653 = -mrSges(5,1) * t709 + mrSges(5,3) * t703 + Ifges(5,4) * t745 + Ifges(5,2) * t744 + Ifges(5,6) * t766 - pkin(4) * t826 + qJ(5) * t844 + t805 * t667 + t809 * t668 - t775 * t738 + t779 * t740;
t739 = Ifges(5,4) * t775 + Ifges(5,2) * t774 + Ifges(5,6) * t779;
t654 = mrSges(5,2) * t709 - mrSges(5,3) * t702 + Ifges(5,1) * t745 + Ifges(5,4) * t744 + Ifges(5,5) * t766 - qJ(5) * t674 + t809 * t667 - t805 * t668 + t774 * t738 - t779 * t739;
t763 = Ifges(4,5) * t781 + Ifges(4,6) * t780 + Ifges(4,3) * t793;
t764 = Ifges(4,4) * t781 + Ifges(4,2) * t780 + Ifges(4,6) * t793;
t640 = mrSges(4,2) * t733 - mrSges(4,3) * t718 + Ifges(4,1) * t770 + Ifges(4,4) * t769 + Ifges(4,5) * t790 - pkin(10) * t666 - t814 * t653 + t818 * t654 + t780 * t763 - t793 * t764;
t765 = Ifges(4,1) * t781 + Ifges(4,4) * t780 + Ifges(4,5) * t793;
t645 = -mrSges(5,1) * t702 + mrSges(5,2) * t703 + mrSges(6,2) * t696 - mrSges(6,1) * t695 - pkin(3) * t666 - pkin(11) * t843 - pkin(5) * t829 + (-Ifges(5,3) - Ifges(6,3)) * t766 - t817 * t682 - t813 * t683 + Ifges(4,6) * t790 + t793 * t765 - t781 * t763 + t774 * t740 - t775 * t739 + Ifges(4,2) * t769 + Ifges(4,4) * t770 - t754 * t728 + t753 * t729 - Ifges(5,5) * t745 - Ifges(5,6) * t744 - mrSges(4,1) * t733 - Ifges(6,6) * t723 - Ifges(6,5) * t724 + mrSges(4,3) * t719 - pkin(4) * t674;
t831 = pkin(9) * t659 + t640 * t815 + t645 * t819;
t639 = mrSges(4,1) * t718 - mrSges(4,2) * t719 + Ifges(4,5) * t770 + Ifges(4,6) * t769 + Ifges(4,3) * t790 + pkin(3) * t822 + pkin(10) * t845 + t818 * t653 + t814 * t654 + t781 * t764 - t780 * t765;
t785 = (t841 * t808 + t865) * qJD(1);
t828 = Ifges(3,5) * t812 + (Ifges(3,1) * t806 + Ifges(3,4) * t810) * t808;
t787 = t828 * qJD(1);
t827 = Ifges(3,6) * t812 + (Ifges(3,4) * t806 + Ifges(3,2) * t810) * t808;
t633 = -mrSges(3,1) * t782 + mrSges(3,3) * t772 - pkin(2) * t651 - t807 * t639 + (-t785 * t862 + t787 * t812) * qJD(1) + t831 * t811 + t827 * qJDD(1);
t786 = t827 * qJD(1);
t634 = mrSges(3,2) * t782 - mrSges(3,3) * t771 + t819 * t640 - t815 * t645 + (t785 * t858 - t786 * t812) * qJD(1) + (-t651 * t807 - t652 * t811) * pkin(9) + t828 * qJDD(1);
t830 = qJ(2) * t644 + t633 * t810 + t634 * t806;
t632 = qJDD(1) * t865 + mrSges(3,1) * t771 - mrSges(3,2) * t772 + pkin(2) * t652 + t811 * t639 + t831 * t807 + (t841 * qJDD(1) + (t786 * t806 - t787 * t810) * qJD(1)) * t808;
t631 = -mrSges(2,2) * g(3) - mrSges(2,3) * t802 + Ifges(2,5) * qJDD(1) - t821 * Ifges(2,6) - t806 * t633 + t810 * t634 + (-t637 * t808 - t638 * t812) * qJ(2);
t630 = mrSges(2,1) * g(3) + mrSges(2,3) * t803 + t821 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t637 - t808 * t632 + t830 * t812;
t1 = [-m(1) * g(1) + t846; -m(1) * g(2) + t853; (-m(1) - m(2)) * g(3) + t637; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t853 - t816 * t630 + t820 * t631; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t846 + t820 * t630 + t816 * t631; -mrSges(1,1) * g(2) + mrSges(2,1) * t802 + mrSges(1,2) * g(1) - mrSges(2,2) * t803 + Ifges(2,3) * qJDD(1) + pkin(1) * t638 + t812 * t632 + t830 * t808;];
tauB  = t1;
