% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:17:24
% EndTime: 2019-05-06 02:18:21
% DurationCPUTime: 54.99s
% Computational Cost: add. (791821->379), mult. (2456931->508), div. (0->0), fcn. (2088823->14), ass. (0->174)
t790 = sin(pkin(12));
t792 = sin(pkin(6));
t793 = cos(pkin(12));
t795 = cos(pkin(6));
t798 = sin(qJ(3));
t794 = cos(pkin(7));
t801 = cos(qJ(3));
t839 = t794 * t801;
t791 = sin(pkin(7));
t844 = t791 * t801;
t806 = (-t790 * t798 + t793 * t839) * t792 + t795 * t844;
t765 = t806 * qJD(1);
t840 = t794 * t798;
t845 = t791 * t798;
t808 = t795 * t845 + (t790 * t801 + t793 * t840) * t792;
t766 = t808 * qJD(1);
t752 = -t766 * qJD(3) + t806 * qJDD(1);
t860 = Ifges(6,1) + Ifges(7,1);
t853 = Ifges(6,4) - Ifges(7,5);
t859 = -Ifges(6,5) - Ifges(7,4);
t858 = Ifges(6,2) + Ifges(7,3);
t851 = Ifges(6,6) - Ifges(7,6);
t857 = -Ifges(6,3) - Ifges(7,2);
t842 = t792 * t794;
t777 = (t791 * t795 + t793 * t842) * qJD(1) * pkin(9);
t799 = sin(qJ(1));
t802 = cos(qJ(1));
t788 = -t802 * g(1) - t799 * g(2);
t803 = qJD(1) ^ 2;
t848 = qJ(2) * t792;
t781 = -t803 * pkin(1) + qJDD(1) * t848 + t788;
t855 = pkin(9) * t790;
t818 = -pkin(2) * t793 - t791 * t855;
t833 = qJD(1) * t792;
t849 = pkin(9) * qJDD(1);
t813 = qJD(1) * t818 * t833 + t794 * t849;
t787 = t799 * g(1) - t802 * g(2);
t780 = qJDD(1) * pkin(1) + t803 * t848 + t787;
t828 = qJD(2) * t833;
t841 = t793 * t795;
t843 = t792 * t793;
t819 = -g(3) * t843 + t780 * t841 - 0.2e1 * t790 * t828;
t733 = (pkin(2) * qJDD(1) + qJD(1) * t777) * t795 + (-t813 * t792 - t781) * t790 + t819;
t782 = (pkin(2) * t795 - t842 * t855) * qJD(1);
t846 = t790 * t795;
t829 = t780 * t846 + (t781 + 0.2e1 * t828) * t793;
t734 = (-qJD(1) * t782 + t791 * t849) * t795 + (-g(3) * t790 + t813 * t793) * t792 + t829;
t827 = -t795 * g(3) + qJDD(2);
t742 = (-t780 + t818 * qJDD(1) + (-t777 * t793 + t782 * t790) * qJD(1)) * t792 + t827;
t695 = -t798 * t734 + (t733 * t794 + t742 * t791) * t801;
t856 = cos(qJ(5));
t854 = -mrSges(6,3) - mrSges(7,2);
t850 = Ifges(3,3) * t795;
t847 = t790 * t792;
t696 = t733 * t840 + t801 * t734 + t742 * t845;
t750 = -t765 * mrSges(4,1) + t766 * mrSges(4,2);
t814 = -t791 * t843 + t794 * t795;
t778 = t814 * qJD(1) + qJD(3);
t761 = t778 * mrSges(4,1) - t766 * mrSges(4,3);
t775 = t814 * qJDD(1) + qJDD(3);
t751 = -t765 * pkin(3) - t766 * pkin(10);
t774 = t778 ^ 2;
t692 = -t774 * pkin(3) + t775 * pkin(10) + t765 * t751 + t696;
t709 = -t791 * t733 + t794 * t742;
t753 = t765 * qJD(3) + t808 * qJDD(1);
t694 = (-t765 * t778 - t753) * pkin(10) + (t766 * t778 - t752) * pkin(3) + t709;
t797 = sin(qJ(4));
t800 = cos(qJ(4));
t688 = t800 * t692 + t797 * t694;
t759 = t800 * t766 + t797 * t778;
t728 = -t759 * qJD(4) - t797 * t753 + t800 * t775;
t758 = -t797 * t766 + t800 * t778;
t735 = -t758 * mrSges(5,1) + t759 * mrSges(5,2);
t764 = qJD(4) - t765;
t744 = t764 * mrSges(5,1) - t759 * mrSges(5,3);
t749 = qJDD(4) - t752;
t736 = -t758 * pkin(4) - t759 * pkin(11);
t763 = t764 ^ 2;
t684 = -t763 * pkin(4) + t749 * pkin(11) + t758 * t736 + t688;
t691 = -t775 * pkin(3) - t774 * pkin(10) + t766 * t751 - t695;
t729 = t758 * qJD(4) + t800 * t753 + t797 * t775;
t686 = (-t758 * t764 - t729) * pkin(11) + (t759 * t764 - t728) * pkin(4) + t691;
t796 = sin(qJ(5));
t681 = t856 * t684 + t796 * t686;
t741 = t856 * t759 + t796 * t764;
t699 = t741 * qJD(5) + t796 * t729 - t856 * t749;
t755 = qJD(5) - t758;
t718 = t755 * mrSges(6,1) - t741 * mrSges(6,3);
t726 = qJDD(5) - t728;
t740 = t796 * t759 - t856 * t764;
t712 = t740 * pkin(5) - t741 * qJ(6);
t754 = t755 ^ 2;
t677 = -t754 * pkin(5) + t726 * qJ(6) + 0.2e1 * qJD(6) * t755 - t740 * t712 + t681;
t719 = -t755 * mrSges(7,1) + t741 * mrSges(7,2);
t830 = m(7) * t677 + t726 * mrSges(7,3) + t755 * t719;
t713 = t740 * mrSges(7,1) - t741 * mrSges(7,3);
t834 = -t740 * mrSges(6,1) - t741 * mrSges(6,2) - t713;
t673 = m(6) * t681 - t726 * mrSges(6,2) + t854 * t699 - t755 * t718 + t834 * t740 + t830;
t680 = -t796 * t684 + t856 * t686;
t700 = -t740 * qJD(5) + t856 * t729 + t796 * t749;
t717 = -t755 * mrSges(6,2) - t740 * mrSges(6,3);
t678 = -t726 * pkin(5) - t754 * qJ(6) + t741 * t712 + qJDD(6) - t680;
t716 = -t740 * mrSges(7,2) + t755 * mrSges(7,3);
t823 = -m(7) * t678 + t726 * mrSges(7,1) + t755 * t716;
t674 = m(6) * t680 + t726 * mrSges(6,1) + t854 * t700 + t755 * t717 + t834 * t741 + t823;
t824 = t856 * t673 - t796 * t674;
t666 = m(5) * t688 - t749 * mrSges(5,2) + t728 * mrSges(5,3) + t758 * t735 - t764 * t744 + t824;
t687 = -t797 * t692 + t800 * t694;
t743 = -t764 * mrSges(5,2) + t758 * mrSges(5,3);
t683 = -t749 * pkin(4) - t763 * pkin(11) + t759 * t736 - t687;
t679 = -0.2e1 * qJD(6) * t741 + (t740 * t755 - t700) * qJ(6) + (t741 * t755 + t699) * pkin(5) + t683;
t675 = m(7) * t679 + t699 * mrSges(7,1) - t700 * mrSges(7,3) + t740 * t716 - t741 * t719;
t804 = -m(6) * t683 - t699 * mrSges(6,1) - t700 * mrSges(6,2) - t740 * t717 - t741 * t718 - t675;
t671 = m(5) * t687 + t749 * mrSges(5,1) - t729 * mrSges(5,3) - t759 * t735 + t764 * t743 + t804;
t825 = t800 * t666 - t797 * t671;
t657 = m(4) * t696 - t775 * mrSges(4,2) + t752 * mrSges(4,3) + t765 * t750 - t778 * t761 + t825;
t660 = t797 * t666 + t800 * t671;
t760 = -t778 * mrSges(4,2) + t765 * mrSges(4,3);
t659 = m(4) * t709 - t752 * mrSges(4,1) + t753 * mrSges(4,2) - t765 * t760 + t766 * t761 + t660;
t669 = t796 * t673 + t856 * t674;
t805 = -m(5) * t691 + t728 * mrSges(5,1) - t729 * mrSges(5,2) + t758 * t743 - t759 * t744 - t669;
t663 = m(4) * t695 + t775 * mrSges(4,1) - t753 * mrSges(4,3) - t766 * t750 + t778 * t760 + t805;
t646 = t657 * t840 - t791 * t659 + t663 * t839;
t756 = -t790 * t781 + t819;
t822 = -mrSges(3,1) * t793 + mrSges(3,2) * t790;
t779 = t822 * t833;
t816 = -mrSges(3,2) * t795 + mrSges(3,3) * t843;
t784 = t816 * qJD(1);
t817 = mrSges(3,1) * t795 - mrSges(3,3) * t847;
t642 = m(3) * t756 + t817 * qJDD(1) + (-t779 * t847 + t784 * t795) * qJD(1) + t646;
t645 = t657 * t845 + t794 * t659 + t663 * t844;
t767 = -t792 * t780 + t827;
t783 = t817 * qJD(1);
t644 = m(3) * t767 + (t822 * qJDD(1) + (t783 * t790 - t784 * t793) * qJD(1)) * t792 + t645;
t652 = t801 * t657 - t798 * t663;
t757 = -g(3) * t847 + t829;
t651 = m(3) * t757 + t816 * qJDD(1) + (t779 * t843 - t783 * t795) * qJD(1) + t652;
t632 = t642 * t841 - t792 * t644 + t651 * t846;
t630 = m(2) * t787 + qJDD(1) * mrSges(2,1) - t803 * mrSges(2,2) + t632;
t638 = -t790 * t642 + t793 * t651;
t637 = m(2) * t788 - t803 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t638;
t838 = t802 * t630 + t799 * t637;
t837 = t858 * t740 - t853 * t741 - t851 * t755;
t836 = t851 * t740 + t859 * t741 + t857 * t755;
t835 = -t853 * t740 + t860 * t741 - t859 * t755;
t631 = t642 * t843 + t795 * t644 + t651 * t847;
t826 = -t799 * t630 + t802 * t637;
t821 = Ifges(3,5) * t790 + Ifges(3,6) * t793;
t667 = -mrSges(6,1) * t683 - mrSges(7,1) * t679 + mrSges(7,2) * t677 + mrSges(6,3) * t681 - pkin(5) * t675 - t858 * t699 + t853 * t700 + t851 * t726 + t836 * t741 + t835 * t755;
t668 = mrSges(6,2) * t683 + mrSges(7,2) * t678 - mrSges(6,3) * t680 - mrSges(7,3) * t679 - qJ(6) * t675 - t853 * t699 + t860 * t700 - t726 * t859 + t836 * t740 + t837 * t755;
t722 = Ifges(5,5) * t759 + Ifges(5,6) * t758 + Ifges(5,3) * t764;
t723 = Ifges(5,4) * t759 + Ifges(5,2) * t758 + Ifges(5,6) * t764;
t647 = mrSges(5,2) * t691 - mrSges(5,3) * t687 + Ifges(5,1) * t729 + Ifges(5,4) * t728 + Ifges(5,5) * t749 - pkin(11) * t669 - t796 * t667 + t856 * t668 + t758 * t722 - t764 * t723;
t724 = Ifges(5,1) * t759 + Ifges(5,4) * t758 + Ifges(5,5) * t764;
t653 = Ifges(5,4) * t729 + Ifges(5,2) * t728 + Ifges(5,6) * t749 - t759 * t722 + t764 * t724 - mrSges(5,1) * t691 + mrSges(5,3) * t688 - mrSges(6,1) * t680 + mrSges(6,2) * t681 + mrSges(7,1) * t678 - mrSges(7,3) * t677 - pkin(5) * t823 - qJ(6) * t830 - pkin(4) * t669 + (pkin(5) * t713 + t837) * t741 + (qJ(6) * t713 - t835) * t740 + t857 * t726 + (pkin(5) * mrSges(7,2) + t859) * t700 + (qJ(6) * mrSges(7,2) + t851) * t699;
t745 = Ifges(4,5) * t766 + Ifges(4,6) * t765 + Ifges(4,3) * t778;
t746 = Ifges(4,4) * t766 + Ifges(4,2) * t765 + Ifges(4,6) * t778;
t634 = mrSges(4,2) * t709 - mrSges(4,3) * t695 + Ifges(4,1) * t753 + Ifges(4,4) * t752 + Ifges(4,5) * t775 - pkin(10) * t660 + t800 * t647 - t797 * t653 + t765 * t745 - t778 * t746;
t747 = Ifges(4,1) * t766 + Ifges(4,4) * t765 + Ifges(4,5) * t778;
t639 = Ifges(4,4) * t753 + Ifges(4,2) * t752 + Ifges(4,6) * t775 - t766 * t745 + t778 * t747 - mrSges(4,1) * t709 + mrSges(4,3) * t696 - Ifges(5,5) * t729 - Ifges(5,6) * t728 - Ifges(5,3) * t749 - t759 * t723 + t758 * t724 - mrSges(5,1) * t687 + mrSges(5,2) * t688 - t796 * t668 - t856 * t667 - pkin(4) * t804 - pkin(11) * t824 - pkin(3) * t660;
t812 = pkin(9) * t652 + t634 * t798 + t639 * t801;
t633 = mrSges(4,1) * t695 - mrSges(4,2) * t696 + Ifges(4,5) * t753 + Ifges(4,6) * t752 + Ifges(4,3) * t775 + pkin(3) * t805 + pkin(10) * t825 + t797 * t647 + t800 * t653 + t766 * t746 - t765 * t747;
t770 = (t821 * t792 + t850) * qJD(1);
t810 = Ifges(3,5) * t795 + (Ifges(3,1) * t790 + Ifges(3,4) * t793) * t792;
t772 = t810 * qJD(1);
t809 = Ifges(3,6) * t795 + (Ifges(3,4) * t790 + Ifges(3,2) * t793) * t792;
t627 = -mrSges(3,1) * t767 + mrSges(3,3) * t757 - pkin(2) * t645 - t791 * t633 + (-t770 * t847 + t772 * t795) * qJD(1) + t812 * t794 + t809 * qJDD(1);
t771 = t809 * qJD(1);
t628 = mrSges(3,2) * t767 - mrSges(3,3) * t756 + t801 * t634 - t798 * t639 + (t770 * t843 - t771 * t795) * qJD(1) + (-t645 * t791 - t646 * t794) * pkin(9) + t810 * qJDD(1);
t811 = qJ(2) * t638 + t627 * t793 + t628 * t790;
t626 = qJDD(1) * t850 + mrSges(3,1) * t756 - mrSges(3,2) * t757 + pkin(2) * t646 + t794 * t633 + t812 * t791 + (t821 * qJDD(1) + (t771 * t790 - t772 * t793) * qJD(1)) * t792;
t625 = -mrSges(2,2) * g(3) - mrSges(2,3) * t787 + Ifges(2,5) * qJDD(1) - t803 * Ifges(2,6) - t790 * t627 + t793 * t628 + (-t631 * t792 - t632 * t795) * qJ(2);
t624 = mrSges(2,1) * g(3) + mrSges(2,3) * t788 + t803 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t631 - t792 * t626 + t811 * t795;
t1 = [-m(1) * g(1) + t826; -m(1) * g(2) + t838; (-m(1) - m(2)) * g(3) + t631; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t838 - t799 * t624 + t802 * t625; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t826 + t802 * t624 + t799 * t625; -mrSges(1,1) * g(2) + mrSges(2,1) * t787 + mrSges(1,2) * g(1) - mrSges(2,2) * t788 + Ifges(2,3) * qJDD(1) + pkin(1) * t632 + t795 * t626 + t811 * t792;];
tauB  = t1;
