% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_invdynB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 18:22:37
% EndTime: 2019-05-08 18:28:55
% DurationCPUTime: 375.13s
% Computational Cost: add. (5997720->432), mult. (15377734->585), div. (0->0), fcn. (13327090->18), ass. (0->196)
t828 = cos(pkin(6));
t821 = qJD(1) * t828 + qJD(2);
t824 = sin(pkin(7));
t827 = cos(pkin(7));
t825 = sin(pkin(6));
t839 = cos(qJ(2));
t859 = qJD(1) * t839;
t856 = t825 * t859;
t806 = (t821 * t824 + t827 * t856) * pkin(11);
t833 = sin(qJ(2));
t861 = qJD(1) * t825;
t881 = pkin(11) * t824;
t810 = (-pkin(2) * t839 - t833 * t881) * t861;
t858 = qJD(1) * qJD(2);
t816 = (qJDD(1) * t833 + t839 * t858) * t825;
t820 = qJDD(1) * t828 + qJDD(2);
t834 = sin(qJ(1));
t840 = cos(qJ(1));
t818 = t834 * g(1) - g(2) * t840;
t841 = qJD(1) ^ 2;
t882 = pkin(10) * t825;
t813 = qJDD(1) * pkin(1) + t841 * t882 + t818;
t819 = -g(1) * t840 - g(2) * t834;
t814 = -pkin(1) * t841 + qJDD(1) * t882 + t819;
t865 = t828 * t839;
t852 = t813 * t865 - t833 * t814;
t860 = qJD(1) * t833;
t880 = pkin(11) * t827;
t771 = -t816 * t880 + t820 * pkin(2) + t821 * t806 + (-g(3) * t839 - t810 * t860) * t825 + t852;
t857 = t825 * t860;
t809 = pkin(2) * t821 - t857 * t880;
t817 = (qJDD(1) * t839 - t833 * t858) * t825;
t848 = t817 * t827 + t820 * t824;
t866 = t828 * t833;
t862 = t813 * t866 + t839 * t814;
t772 = -t821 * t809 + (-g(3) * t833 + t810 * t859) * t825 + t848 * pkin(11) + t862;
t877 = t828 * g(3);
t777 = -t816 * t881 - t817 * pkin(2) - t877 + (-t813 + (-t806 * t839 + t809 * t833) * qJD(1)) * t825;
t832 = sin(qJ(3));
t838 = cos(qJ(3));
t868 = t827 * t838;
t873 = t824 * t838;
t745 = t771 * t868 - t772 * t832 + t777 * t873;
t867 = t827 * t839;
t797 = t821 * t873 + (-t832 * t833 + t838 * t867) * t861;
t786 = t797 * qJD(3) + t838 * t816 + t832 * t848;
t874 = t824 * t832;
t798 = t821 * t874 + (t832 * t867 + t833 * t838) * t861;
t823 = sin(pkin(8));
t879 = pkin(12) * t823;
t787 = -pkin(3) * t797 - t798 * t879;
t807 = t821 * t827 - t824 * t856 + qJD(3);
t826 = cos(pkin(8));
t849 = t797 * t826 + t807 * t823;
t790 = t849 * pkin(12);
t799 = -t817 * t824 + t820 * t827 + qJDD(3);
t878 = pkin(12) * t826;
t729 = pkin(3) * t799 - t786 * t878 - t787 * t798 + t790 * t807 + t745;
t869 = t827 * t832;
t746 = t771 * t869 + t838 * t772 + t777 * t874;
t792 = pkin(3) * t807 - t798 * t878;
t785 = -t798 * qJD(3) - t832 * t816 + t838 * t848;
t850 = t785 * t826 + t799 * t823;
t730 = pkin(12) * t850 + t797 * t787 - t807 * t792 + t746;
t760 = -t771 * t824 + t827 * t777;
t737 = -pkin(3) * t785 - t786 * t879 - t790 * t797 + t792 * t798 + t760;
t831 = sin(qJ(4));
t837 = cos(qJ(4));
t716 = -t831 * t730 + (t729 * t826 + t737 * t823) * t837;
t781 = t837 * t798 + t831 * t849;
t750 = -t781 * qJD(4) - t831 * t786 + t837 * t850;
t780 = -t831 * t798 + t837 * t849;
t751 = t780 * qJD(4) + t837 * t786 + t831 * t850;
t761 = -mrSges(5,1) * t780 + mrSges(5,2) * t781;
t791 = -t797 * t823 + t807 * t826 + qJD(4);
t766 = -mrSges(5,2) * t791 + mrSges(5,3) * t780;
t773 = -t785 * t823 + t799 * t826 + qJDD(4);
t870 = t826 * t831;
t875 = t823 * t831;
t717 = t729 * t870 + t837 * t730 + t737 * t875;
t762 = -pkin(4) * t780 - pkin(13) * t781;
t789 = t791 ^ 2;
t713 = -pkin(4) * t789 + pkin(13) * t773 + t762 * t780 + t717;
t718 = -t823 * t729 + t826 * t737;
t715 = (-t780 * t791 - t751) * pkin(13) + (t781 * t791 - t750) * pkin(4) + t718;
t830 = sin(qJ(5));
t836 = cos(qJ(5));
t709 = t836 * t713 + t830 * t715;
t764 = -t781 * t830 + t791 * t836;
t765 = t781 * t836 + t791 * t830;
t748 = -pkin(5) * t764 - pkin(14) * t765;
t749 = qJDD(5) - t750;
t779 = qJD(5) - t780;
t778 = t779 ^ 2;
t707 = -pkin(5) * t778 + pkin(14) * t749 + t748 * t764 + t709;
t712 = -t773 * pkin(4) - t789 * pkin(13) + t781 * t762 - t716;
t733 = -qJD(5) * t765 - t751 * t830 + t773 * t836;
t734 = qJD(5) * t764 + t751 * t836 + t773 * t830;
t710 = (-t764 * t779 - t734) * pkin(14) + (t765 * t779 - t733) * pkin(5) + t712;
t829 = sin(qJ(6));
t835 = cos(qJ(6));
t704 = -t707 * t829 + t710 * t835;
t753 = -t765 * t829 + t779 * t835;
t721 = qJD(6) * t753 + t734 * t835 + t749 * t829;
t732 = qJDD(6) - t733;
t754 = t765 * t835 + t779 * t829;
t738 = -mrSges(7,1) * t753 + mrSges(7,2) * t754;
t763 = qJD(6) - t764;
t739 = -mrSges(7,2) * t763 + mrSges(7,3) * t753;
t702 = m(7) * t704 + mrSges(7,1) * t732 - mrSges(7,3) * t721 - t738 * t754 + t739 * t763;
t705 = t707 * t835 + t710 * t829;
t720 = -qJD(6) * t754 - t734 * t829 + t749 * t835;
t740 = mrSges(7,1) * t763 - mrSges(7,3) * t754;
t703 = m(7) * t705 - mrSges(7,2) * t732 + mrSges(7,3) * t720 + t738 * t753 - t740 * t763;
t696 = t702 * t835 + t703 * t829;
t755 = -mrSges(6,2) * t779 + mrSges(6,3) * t764;
t756 = mrSges(6,1) * t779 - mrSges(6,3) * t765;
t842 = -m(6) * t712 + t733 * mrSges(6,1) - mrSges(6,2) * t734 + t764 * t755 - t756 * t765 - t696;
t692 = m(5) * t716 + mrSges(5,1) * t773 - mrSges(5,3) * t751 - t761 * t781 + t766 * t791 + t842;
t876 = t692 * t837;
t872 = t825 * t833;
t871 = t825 * t839;
t767 = mrSges(5,1) * t791 - mrSges(5,3) * t781;
t747 = -mrSges(6,1) * t764 + mrSges(6,2) * t765;
t853 = -t702 * t829 + t835 * t703;
t695 = m(6) * t709 - mrSges(6,2) * t749 + mrSges(6,3) * t733 + t747 * t764 - t756 * t779 + t853;
t708 = -t713 * t830 + t715 * t836;
t706 = -pkin(5) * t749 - pkin(14) * t778 + t748 * t765 - t708;
t843 = -m(7) * t706 + t720 * mrSges(7,1) - mrSges(7,2) * t721 + t753 * t739 - t740 * t754;
t700 = m(6) * t708 + mrSges(6,1) * t749 - mrSges(6,3) * t734 - t747 * t765 + t755 * t779 + t843;
t854 = t836 * t695 - t700 * t830;
t686 = m(5) * t717 - mrSges(5,2) * t773 + mrSges(5,3) * t750 + t761 * t780 - t767 * t791 + t854;
t689 = t830 * t695 + t836 * t700;
t688 = m(5) * t718 - mrSges(5,1) * t750 + mrSges(5,2) * t751 - t766 * t780 + t767 * t781 + t689;
t675 = t686 * t870 - t688 * t823 + t826 * t876;
t788 = -mrSges(4,1) * t797 + mrSges(4,2) * t798;
t793 = -mrSges(4,2) * t807 + mrSges(4,3) * t797;
t671 = m(4) * t745 + mrSges(4,1) * t799 - mrSges(4,3) * t786 - t788 * t798 + t793 * t807 + t675;
t674 = t686 * t875 + t826 * t688 + t823 * t876;
t794 = mrSges(4,1) * t807 - mrSges(4,3) * t798;
t673 = m(4) * t760 - mrSges(4,1) * t785 + mrSges(4,2) * t786 - t793 * t797 + t794 * t798 + t674;
t680 = t837 * t686 - t692 * t831;
t679 = m(4) * t746 - mrSges(4,2) * t799 + mrSges(4,3) * t785 + t788 * t797 - t794 * t807 + t680;
t660 = t671 * t868 - t673 * t824 + t679 * t869;
t795 = -g(3) * t871 + t852;
t812 = -mrSges(3,2) * t821 + mrSges(3,3) * t856;
t815 = (-mrSges(3,1) * t839 + mrSges(3,2) * t833) * t861;
t656 = m(3) * t795 + mrSges(3,1) * t820 - mrSges(3,3) * t816 + t812 * t821 - t815 * t857 + t660;
t659 = t671 * t873 + t827 * t673 + t679 * t874;
t803 = -t825 * t813 - t877;
t811 = mrSges(3,1) * t821 - mrSges(3,3) * t857;
t658 = m(3) * t803 - t817 * mrSges(3,1) + t816 * mrSges(3,2) + (t811 * t833 - t812 * t839) * t861 + t659;
t665 = -t671 * t832 + t838 * t679;
t796 = -g(3) * t872 + t862;
t664 = m(3) * t796 - mrSges(3,2) * t820 + mrSges(3,3) * t817 - t811 * t821 + t815 * t856 + t665;
t646 = t656 * t865 - t658 * t825 + t664 * t866;
t644 = m(2) * t818 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t841 + t646;
t653 = -t656 * t833 + t839 * t664;
t652 = m(2) * t819 - mrSges(2,1) * t841 - qJDD(1) * mrSges(2,2) + t653;
t863 = t840 * t644 + t834 * t652;
t645 = t656 * t871 + t828 * t658 + t664 * t872;
t855 = -t644 * t834 + t840 * t652;
t722 = Ifges(7,5) * t754 + Ifges(7,6) * t753 + Ifges(7,3) * t763;
t724 = Ifges(7,1) * t754 + Ifges(7,4) * t753 + Ifges(7,5) * t763;
t697 = -mrSges(7,1) * t706 + mrSges(7,3) * t705 + Ifges(7,4) * t721 + Ifges(7,2) * t720 + Ifges(7,6) * t732 - t722 * t754 + t724 * t763;
t723 = Ifges(7,4) * t754 + Ifges(7,2) * t753 + Ifges(7,6) * t763;
t698 = mrSges(7,2) * t706 - mrSges(7,3) * t704 + Ifges(7,1) * t721 + Ifges(7,4) * t720 + Ifges(7,5) * t732 + t722 * t753 - t723 * t763;
t741 = Ifges(6,5) * t765 + Ifges(6,6) * t764 + Ifges(6,3) * t779;
t742 = Ifges(6,4) * t765 + Ifges(6,2) * t764 + Ifges(6,6) * t779;
t681 = mrSges(6,2) * t712 - mrSges(6,3) * t708 + Ifges(6,1) * t734 + Ifges(6,4) * t733 + Ifges(6,5) * t749 - pkin(14) * t696 - t697 * t829 + t698 * t835 + t741 * t764 - t742 * t779;
t743 = Ifges(6,1) * t765 + Ifges(6,4) * t764 + Ifges(6,5) * t779;
t682 = -mrSges(6,1) * t712 - mrSges(7,1) * t704 + mrSges(7,2) * t705 + mrSges(6,3) * t709 + Ifges(6,4) * t734 - Ifges(7,5) * t721 + Ifges(6,2) * t733 + Ifges(6,6) * t749 - Ifges(7,6) * t720 - Ifges(7,3) * t732 - pkin(5) * t696 - t723 * t754 + t724 * t753 - t741 * t765 + t743 * t779;
t758 = Ifges(5,4) * t781 + Ifges(5,2) * t780 + Ifges(5,6) * t791;
t759 = Ifges(5,1) * t781 + Ifges(5,4) * t780 + Ifges(5,5) * t791;
t666 = mrSges(5,1) * t716 - mrSges(5,2) * t717 + Ifges(5,5) * t751 + Ifges(5,6) * t750 + Ifges(5,3) * t773 + pkin(4) * t842 + pkin(13) * t854 + t830 * t681 + t836 * t682 + t781 * t758 - t780 * t759;
t783 = Ifges(4,4) * t798 + Ifges(4,2) * t797 + Ifges(4,6) * t807;
t784 = Ifges(4,1) * t798 + Ifges(4,4) * t797 + Ifges(4,5) * t807;
t757 = Ifges(5,5) * t781 + Ifges(5,6) * t780 + Ifges(5,3) * t791;
t667 = mrSges(5,2) * t718 - mrSges(5,3) * t716 + Ifges(5,1) * t751 + Ifges(5,4) * t750 + Ifges(5,5) * t773 - pkin(13) * t689 + t681 * t836 - t682 * t830 + t757 * t780 - t758 * t791;
t668 = Ifges(5,4) * t751 + Ifges(5,2) * t750 + Ifges(5,6) * t773 - t781 * t757 + t791 * t759 - mrSges(5,1) * t718 + mrSges(5,3) * t717 - Ifges(6,5) * t734 - Ifges(6,6) * t733 - Ifges(6,3) * t749 - t765 * t742 + t764 * t743 - mrSges(6,1) * t708 + mrSges(6,2) * t709 - t829 * t698 - t835 * t697 - pkin(5) * t843 - pkin(14) * t853 - pkin(4) * t689;
t844 = pkin(12) * t680 + t667 * t831 + t668 * t837;
t647 = mrSges(4,1) * t745 - mrSges(4,2) * t746 + Ifges(4,5) * t786 + Ifges(4,6) * t785 + Ifges(4,3) * t799 + pkin(3) * t675 + t826 * t666 + t798 * t783 - t797 * t784 + t823 * t844;
t800 = Ifges(3,3) * t821 + (Ifges(3,5) * t833 + Ifges(3,6) * t839) * t861;
t802 = Ifges(3,5) * t821 + (Ifges(3,1) * t833 + Ifges(3,4) * t839) * t861;
t782 = Ifges(4,5) * t798 + Ifges(4,6) * t797 + Ifges(4,3) * t807;
t648 = -mrSges(4,1) * t760 + mrSges(4,3) * t746 + Ifges(4,4) * t786 + Ifges(4,2) * t785 + Ifges(4,6) * t799 - pkin(3) * t674 - t823 * t666 - t798 * t782 + t807 * t784 + t826 * t844;
t649 = mrSges(4,2) * t760 - mrSges(4,3) * t745 + Ifges(4,1) * t786 + Ifges(4,4) * t785 + Ifges(4,5) * t799 + t837 * t667 - t831 * t668 + t797 * t782 - t807 * t783 + (-t674 * t823 - t675 * t826) * pkin(12);
t845 = pkin(11) * t665 + t648 * t838 + t649 * t832;
t641 = -mrSges(3,1) * t803 + mrSges(3,3) * t796 + Ifges(3,4) * t816 + Ifges(3,2) * t817 + Ifges(3,6) * t820 - pkin(2) * t659 - t824 * t647 - t800 * t857 + t821 * t802 + t827 * t845;
t801 = Ifges(3,6) * t821 + (Ifges(3,4) * t833 + Ifges(3,2) * t839) * t861;
t642 = t800 * t856 + mrSges(3,2) * t803 - mrSges(3,3) * t795 + Ifges(3,1) * t816 + Ifges(3,4) * t817 + Ifges(3,5) * t820 - t832 * t648 + t838 * t649 - t821 * t801 + (-t659 * t824 - t660 * t827) * pkin(11);
t846 = pkin(10) * t653 + t641 * t839 + t642 * t833;
t640 = mrSges(3,1) * t795 - mrSges(3,2) * t796 + Ifges(3,5) * t816 + Ifges(3,6) * t817 + Ifges(3,3) * t820 + pkin(2) * t660 + t827 * t647 + (t801 * t833 - t802 * t839) * t861 + t845 * t824;
t639 = -mrSges(2,2) * g(3) - mrSges(2,3) * t818 + Ifges(2,5) * qJDD(1) - t841 * Ifges(2,6) - t833 * t641 + t839 * t642 + (-t645 * t825 - t646 * t828) * pkin(10);
t638 = mrSges(2,1) * g(3) + mrSges(2,3) * t819 + t841 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t645 - t825 * t640 + t828 * t846;
t1 = [-m(1) * g(1) + t855; -m(1) * g(2) + t863; (-m(1) - m(2)) * g(3) + t645; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(9) * t863 - t834 * t638 + t840 * t639; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(9) * t855 + t840 * t638 + t834 * t639; -mrSges(1,1) * g(2) + mrSges(2,1) * t818 + mrSges(1,2) * g(1) - mrSges(2,2) * t819 + Ifges(2,3) * qJDD(1) + pkin(1) * t646 + t828 * t640 + t825 * t846;];
tauB  = t1;
