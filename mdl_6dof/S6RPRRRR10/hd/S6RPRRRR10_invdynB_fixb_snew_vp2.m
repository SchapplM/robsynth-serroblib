% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 04:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:45:22
% EndTime: 2019-05-06 04:47:09
% DurationCPUTime: 110.60s
% Computational Cost: add. (1743631->404), mult. (5407276->549), div. (0->0), fcn. (4674466->16), ass. (0->187)
t824 = sin(pkin(13));
t826 = sin(pkin(6));
t827 = cos(pkin(13));
t829 = cos(pkin(6));
t838 = cos(qJ(3));
t828 = cos(pkin(7));
t833 = sin(qJ(3));
t872 = t828 * t833;
t825 = sin(pkin(7));
t876 = t825 * t833;
t842 = t829 * t876 + (t824 * t838 + t827 * t872) * t826;
t798 = t842 * qJD(1);
t874 = t826 * t828;
t853 = t825 * t829 + t827 * t874;
t849 = t853 * t838;
t878 = t824 * t826;
t868 = t833 * t878;
t784 = -t798 * qJD(3) + (t849 - t868) * qJDD(1);
t847 = t853 * qJD(1);
t809 = pkin(9) * t847;
t834 = sin(qJ(1));
t839 = cos(qJ(1));
t822 = -t839 * g(1) - t834 * g(2);
t840 = qJD(1) ^ 2;
t880 = qJ(2) * t826;
t813 = -t840 * pkin(1) + qJDD(1) * t880 + t822;
t883 = pkin(9) * t824;
t856 = -pkin(2) * t827 - t825 * t883;
t869 = qJD(1) * t826;
t881 = pkin(9) * qJDD(1);
t851 = qJD(1) * t856 * t869 + t828 * t881;
t821 = t834 * g(1) - t839 * g(2);
t812 = qJDD(1) * pkin(1) + t840 * t880 + t821;
t866 = qJD(2) * t869;
t873 = t827 * t829;
t875 = t826 * t827;
t857 = -g(3) * t875 + t812 * t873 - 0.2e1 * t824 * t866;
t764 = (pkin(2) * qJDD(1) + qJD(1) * t809) * t829 + (-t826 * t851 - t813) * t824 + t857;
t814 = (pkin(2) * t829 - t874 * t883) * qJD(1);
t877 = t824 * t829;
t867 = t812 * t877 + (t813 + 0.2e1 * t866) * t827;
t765 = (-qJD(1) * t814 + t825 * t881) * t829 + (-g(3) * t824 + t827 * t851) * t826 + t867;
t865 = -t829 * g(3) + qJDD(2);
t773 = (-t812 + t856 * qJDD(1) + (-t809 * t827 + t814 * t824) * qJD(1)) * t826 + t865;
t737 = -t833 * t765 + (t764 * t828 + t773 * t825) * t838;
t882 = Ifges(3,3) * t829;
t818 = qJD(1) * t868;
t797 = t838 * t847 - t818;
t782 = -mrSges(4,1) * t797 + mrSges(4,2) * t798;
t785 = t797 * qJD(3) + qJDD(1) * t842;
t852 = -t825 * t875 + t828 * t829;
t810 = qJD(1) * t852 + qJD(3);
t791 = -mrSges(4,2) * t810 + mrSges(4,3) * t797;
t807 = qJDD(1) * t852 + qJDD(3);
t783 = -pkin(3) * t797 - pkin(10) * t798;
t806 = t810 ^ 2;
t724 = -t807 * pkin(3) - t806 * pkin(10) + t798 * t783 - t737;
t832 = sin(qJ(4));
t837 = cos(qJ(4));
t790 = t798 * t837 + t810 * t832;
t759 = -qJD(4) * t790 - t785 * t832 + t807 * t837;
t789 = -t798 * t832 + t810 * t837;
t760 = qJD(4) * t789 + t785 * t837 + t807 * t832;
t796 = -qJD(1) * t849 + qJD(4) + t818;
t774 = -mrSges(5,2) * t796 + mrSges(5,3) * t789;
t775 = mrSges(5,1) * t796 - mrSges(5,3) * t790;
t738 = t764 * t872 + t838 * t765 + t773 * t876;
t725 = -pkin(3) * t806 + pkin(10) * t807 + t783 * t797 + t738;
t748 = -t825 * t764 + t828 * t773;
t728 = (-t797 * t810 - t785) * pkin(10) + (t798 * t810 - t784) * pkin(3) + t748;
t717 = -t832 * t725 + t837 * t728;
t781 = qJDD(4) - t784;
t714 = (t789 * t796 - t760) * pkin(11) + (t789 * t790 + t781) * pkin(4) + t717;
t718 = t837 * t725 + t832 * t728;
t776 = pkin(4) * t796 - pkin(11) * t790;
t788 = t789 ^ 2;
t716 = -pkin(4) * t788 + pkin(11) * t759 - t776 * t796 + t718;
t831 = sin(qJ(5));
t836 = cos(qJ(5));
t711 = t831 * t714 + t836 * t716;
t767 = t789 * t836 - t790 * t831;
t768 = t789 * t831 + t790 * t836;
t747 = -pkin(5) * t767 - pkin(12) * t768;
t780 = qJDD(5) + t781;
t794 = qJD(5) + t796;
t793 = t794 ^ 2;
t709 = -pkin(5) * t793 + pkin(12) * t780 + t747 * t767 + t711;
t719 = -t759 * pkin(4) - t788 * pkin(11) + t790 * t776 + t724;
t735 = -qJD(5) * t768 + t759 * t836 - t760 * t831;
t736 = qJD(5) * t767 + t759 * t831 + t760 * t836;
t712 = (-t767 * t794 - t736) * pkin(12) + (t768 * t794 - t735) * pkin(5) + t719;
t830 = sin(qJ(6));
t835 = cos(qJ(6));
t706 = -t709 * t830 + t712 * t835;
t749 = -t768 * t830 + t794 * t835;
t722 = qJD(6) * t749 + t736 * t835 + t780 * t830;
t734 = qJDD(6) - t735;
t750 = t768 * t835 + t794 * t830;
t739 = -mrSges(7,1) * t749 + mrSges(7,2) * t750;
t766 = qJD(6) - t767;
t740 = -mrSges(7,2) * t766 + mrSges(7,3) * t749;
t704 = m(7) * t706 + mrSges(7,1) * t734 - mrSges(7,3) * t722 - t739 * t750 + t740 * t766;
t707 = t709 * t835 + t712 * t830;
t721 = -qJD(6) * t750 - t736 * t830 + t780 * t835;
t741 = mrSges(7,1) * t766 - mrSges(7,3) * t750;
t705 = m(7) * t707 - mrSges(7,2) * t734 + mrSges(7,3) * t721 + t739 * t749 - t741 * t766;
t696 = t835 * t704 + t830 * t705;
t751 = -mrSges(6,2) * t794 + mrSges(6,3) * t767;
t752 = mrSges(6,1) * t794 - mrSges(6,3) * t768;
t843 = m(6) * t719 - t735 * mrSges(6,1) + mrSges(6,2) * t736 - t767 * t751 + t752 * t768 + t696;
t841 = -m(5) * t724 + t759 * mrSges(5,1) - mrSges(5,2) * t760 + t789 * t774 - t775 * t790 - t843;
t692 = m(4) * t737 + mrSges(4,1) * t807 - mrSges(4,3) * t785 - t782 * t798 + t791 * t810 + t841;
t879 = t692 * t838;
t792 = mrSges(4,1) * t810 - mrSges(4,3) * t798;
t746 = -mrSges(6,1) * t767 + mrSges(6,2) * t768;
t861 = -t704 * t830 + t835 * t705;
t695 = m(6) * t711 - mrSges(6,2) * t780 + mrSges(6,3) * t735 + t746 * t767 - t752 * t794 + t861;
t710 = t714 * t836 - t716 * t831;
t708 = -pkin(5) * t780 - pkin(12) * t793 + t747 * t768 - t710;
t846 = -m(7) * t708 + t721 * mrSges(7,1) - mrSges(7,2) * t722 + t749 * t740 - t741 * t750;
t700 = m(6) * t710 + mrSges(6,1) * t780 - mrSges(6,3) * t736 - t746 * t768 + t751 * t794 + t846;
t689 = t831 * t695 + t836 * t700;
t769 = -mrSges(5,1) * t789 + mrSges(5,2) * t790;
t687 = m(5) * t717 + mrSges(5,1) * t781 - mrSges(5,3) * t760 - t769 * t790 + t774 * t796 + t689;
t862 = t836 * t695 - t700 * t831;
t688 = m(5) * t718 - mrSges(5,2) * t781 + mrSges(5,3) * t759 + t769 * t789 - t775 * t796 + t862;
t863 = -t687 * t832 + t837 * t688;
t678 = m(4) * t738 - mrSges(4,2) * t807 + mrSges(4,3) * t784 + t782 * t797 - t792 * t810 + t863;
t681 = t837 * t687 + t832 * t688;
t680 = m(4) * t748 - mrSges(4,1) * t784 + mrSges(4,2) * t785 - t791 * t797 + t792 * t798 + t681;
t667 = t678 * t872 - t825 * t680 + t828 * t879;
t786 = -t824 * t813 + t857;
t860 = -mrSges(3,1) * t827 + mrSges(3,2) * t824;
t811 = t860 * t869;
t854 = -mrSges(3,2) * t829 + mrSges(3,3) * t875;
t816 = t854 * qJD(1);
t855 = mrSges(3,1) * t829 - mrSges(3,3) * t878;
t663 = m(3) * t786 + t855 * qJDD(1) + (-t811 * t878 + t816 * t829) * qJD(1) + t667;
t666 = t678 * t876 + t828 * t680 + t825 * t879;
t799 = -t826 * t812 + t865;
t815 = t855 * qJD(1);
t665 = m(3) * t799 + (t860 * qJDD(1) + (t815 * t824 - t816 * t827) * qJD(1)) * t826 + t666;
t674 = t838 * t678 - t833 * t692;
t787 = -g(3) * t878 + t867;
t673 = m(3) * t787 + t854 * qJDD(1) + (t811 * t875 - t815 * t829) * qJD(1) + t674;
t653 = t663 * t873 - t665 * t826 + t673 * t877;
t651 = m(2) * t821 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t840 + t653;
t659 = -t663 * t824 + t827 * t673;
t658 = m(2) * t822 - mrSges(2,1) * t840 - qJDD(1) * mrSges(2,2) + t659;
t871 = t839 * t651 + t834 * t658;
t652 = t663 * t875 + t829 * t665 + t673 * t878;
t864 = -t834 * t651 + t839 * t658;
t859 = Ifges(3,5) * t824 + Ifges(3,6) * t827;
t729 = Ifges(7,5) * t750 + Ifges(7,6) * t749 + Ifges(7,3) * t766;
t731 = Ifges(7,1) * t750 + Ifges(7,4) * t749 + Ifges(7,5) * t766;
t697 = -mrSges(7,1) * t708 + mrSges(7,3) * t707 + Ifges(7,4) * t722 + Ifges(7,2) * t721 + Ifges(7,6) * t734 - t729 * t750 + t731 * t766;
t730 = Ifges(7,4) * t750 + Ifges(7,2) * t749 + Ifges(7,6) * t766;
t698 = mrSges(7,2) * t708 - mrSges(7,3) * t706 + Ifges(7,1) * t722 + Ifges(7,4) * t721 + Ifges(7,5) * t734 + t729 * t749 - t730 * t766;
t742 = Ifges(6,5) * t768 + Ifges(6,6) * t767 + Ifges(6,3) * t794;
t743 = Ifges(6,4) * t768 + Ifges(6,2) * t767 + Ifges(6,6) * t794;
t682 = mrSges(6,2) * t719 - mrSges(6,3) * t710 + Ifges(6,1) * t736 + Ifges(6,4) * t735 + Ifges(6,5) * t780 - pkin(12) * t696 - t697 * t830 + t698 * t835 + t742 * t767 - t743 * t794;
t744 = Ifges(6,1) * t768 + Ifges(6,4) * t767 + Ifges(6,5) * t794;
t683 = -mrSges(6,1) * t719 - mrSges(7,1) * t706 + mrSges(7,2) * t707 + mrSges(6,3) * t711 + Ifges(6,4) * t736 - Ifges(7,5) * t722 + Ifges(6,2) * t735 + Ifges(6,6) * t780 - Ifges(7,6) * t721 - Ifges(7,3) * t734 - pkin(5) * t696 - t730 * t750 + t731 * t749 - t742 * t768 + t744 * t794;
t753 = Ifges(5,5) * t790 + Ifges(5,6) * t789 + Ifges(5,3) * t796;
t755 = Ifges(5,1) * t790 + Ifges(5,4) * t789 + Ifges(5,5) * t796;
t668 = -mrSges(5,1) * t724 + mrSges(5,3) * t718 + Ifges(5,4) * t760 + Ifges(5,2) * t759 + Ifges(5,6) * t781 - pkin(4) * t843 + pkin(11) * t862 + t831 * t682 + t836 * t683 - t790 * t753 + t796 * t755;
t754 = Ifges(5,4) * t790 + Ifges(5,2) * t789 + Ifges(5,6) * t796;
t669 = mrSges(5,2) * t724 - mrSges(5,3) * t717 + Ifges(5,1) * t760 + Ifges(5,4) * t759 + Ifges(5,5) * t781 - pkin(11) * t689 + t682 * t836 - t683 * t831 + t753 * t789 - t754 * t796;
t777 = Ifges(4,5) * t798 + Ifges(4,6) * t797 + Ifges(4,3) * t810;
t778 = Ifges(4,4) * t798 + Ifges(4,2) * t797 + Ifges(4,6) * t810;
t655 = mrSges(4,2) * t748 - mrSges(4,3) * t737 + Ifges(4,1) * t785 + Ifges(4,4) * t784 + Ifges(4,5) * t807 - pkin(10) * t681 - t668 * t832 + t669 * t837 + t777 * t797 - t778 * t810;
t779 = Ifges(4,1) * t798 + Ifges(4,4) * t797 + Ifges(4,5) * t810;
t660 = -pkin(5) * t846 - pkin(12) * t861 - mrSges(6,1) * t710 - t768 * t743 - Ifges(6,6) * t735 - Ifges(6,5) * t736 + mrSges(4,3) * t738 - t830 * t698 - t835 * t697 + t767 * t744 - t798 * t777 + Ifges(4,6) * t807 + t810 * t779 + t789 * t755 - t790 * t754 - pkin(3) * t681 - Ifges(5,6) * t759 - Ifges(5,5) * t760 - mrSges(4,1) * t748 - Ifges(6,3) * t780 - Ifges(5,3) * t781 + Ifges(4,2) * t784 + Ifges(4,4) * t785 - pkin(4) * t689 - mrSges(5,1) * t717 + mrSges(5,2) * t718 + mrSges(6,2) * t711;
t850 = pkin(9) * t674 + t655 * t833 + t660 * t838;
t654 = mrSges(4,1) * t737 - mrSges(4,2) * t738 + Ifges(4,5) * t785 + Ifges(4,6) * t784 + Ifges(4,3) * t807 + pkin(3) * t841 + pkin(10) * t863 + t837 * t668 + t832 * t669 + t798 * t778 - t797 * t779;
t802 = (t826 * t859 + t882) * qJD(1);
t845 = Ifges(3,5) * t829 + (Ifges(3,1) * t824 + Ifges(3,4) * t827) * t826;
t804 = t845 * qJD(1);
t844 = Ifges(3,6) * t829 + (Ifges(3,4) * t824 + Ifges(3,2) * t827) * t826;
t648 = -mrSges(3,1) * t799 + mrSges(3,3) * t787 - pkin(2) * t666 - t825 * t654 + (-t802 * t878 + t804 * t829) * qJD(1) + t850 * t828 + t844 * qJDD(1);
t803 = t844 * qJD(1);
t649 = mrSges(3,2) * t799 - mrSges(3,3) * t786 + t838 * t655 - t833 * t660 + (t802 * t875 - t803 * t829) * qJD(1) + (-t666 * t825 - t667 * t828) * pkin(9) + t845 * qJDD(1);
t848 = qJ(2) * t659 + t648 * t827 + t649 * t824;
t647 = qJDD(1) * t882 + mrSges(3,1) * t786 - mrSges(3,2) * t787 + pkin(2) * t667 + t828 * t654 + t850 * t825 + (t859 * qJDD(1) + (t803 * t824 - t804 * t827) * qJD(1)) * t826;
t646 = -mrSges(2,2) * g(3) - mrSges(2,3) * t821 + Ifges(2,5) * qJDD(1) - t840 * Ifges(2,6) - t824 * t648 + t827 * t649 + (-t652 * t826 - t653 * t829) * qJ(2);
t645 = mrSges(2,1) * g(3) + mrSges(2,3) * t822 + t840 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t652 - t826 * t647 + t829 * t848;
t1 = [-m(1) * g(1) + t864; -m(1) * g(2) + t871; (-m(1) - m(2)) * g(3) + t652; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t871 - t834 * t645 + t839 * t646; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t864 + t839 * t645 + t834 * t646; -mrSges(1,1) * g(2) + mrSges(2,1) * t821 + mrSges(1,2) * g(1) - mrSges(2,2) * t822 + Ifges(2,3) * qJDD(1) + pkin(1) * t653 + t829 * t647 + t826 * t848;];
tauB  = t1;
