% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-09 12:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_invdynB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-09 11:30:14
% EndTime: 2019-05-09 11:35:59
% DurationCPUTime: 336.66s
% Computational Cost: add. (5264401->431), mult. (14277968->587), div. (0->0), fcn. (12309530->18), ass. (0->196)
t839 = cos(pkin(6));
t830 = qJD(1) * t839 + qJD(2);
t834 = sin(pkin(7));
t838 = cos(pkin(7));
t835 = sin(pkin(6));
t848 = cos(qJ(2));
t868 = qJD(1) * t848;
t865 = t835 * t868;
t813 = (t830 * t834 + t838 * t865) * qJ(3);
t843 = sin(qJ(2));
t870 = qJD(1) * t835;
t887 = qJ(3) * t834;
t819 = (-pkin(2) * t848 - t843 * t887) * t870;
t867 = qJD(1) * qJD(2);
t825 = (qJDD(1) * t843 + t848 * t867) * t835;
t829 = qJDD(1) * t839 + qJDD(2);
t844 = sin(qJ(1));
t849 = cos(qJ(1));
t827 = t844 * g(1) - g(2) * t849;
t850 = qJD(1) ^ 2;
t891 = pkin(10) * t835;
t822 = qJDD(1) * pkin(1) + t850 * t891 + t827;
t828 = -g(1) * t849 - g(2) * t844;
t823 = -pkin(1) * t850 + qJDD(1) * t891 + t828;
t874 = t839 * t848;
t861 = t822 * t874 - t843 * t823;
t869 = qJD(1) * t843;
t886 = qJ(3) * t838;
t778 = -t825 * t886 + t829 * pkin(2) + t830 * t813 + (-g(3) * t848 - t819 * t869) * t835 + t861;
t866 = t835 * t869;
t818 = pkin(2) * t830 - t866 * t886;
t826 = (qJDD(1) * t848 - t843 * t867) * t835;
t857 = t826 * t838 + t829 * t834;
t875 = t839 * t843;
t871 = t822 * t875 + t848 * t823;
t779 = -t830 * t818 + (-g(3) * t843 + t819 * t868) * t835 + t857 * qJ(3) + t871;
t888 = t839 * g(3);
t783 = -t825 * t887 - t826 * pkin(2) - t888 + (-t822 + (-t813 * t848 + t818 * t843) * qJD(1)) * t835;
t832 = sin(pkin(14));
t836 = cos(pkin(14));
t876 = t838 * t848;
t884 = t832 * t834;
t807 = t830 * t884 + (t832 * t876 + t836 * t843) * t870;
t878 = t836 * t838;
t881 = t834 * t836;
t752 = -0.2e1 * qJD(3) * t807 + t778 * t878 - t779 * t832 + t783 * t881;
t806 = t830 * t881 + (-t832 * t843 + t836 * t876) * t870;
t833 = sin(pkin(8));
t890 = pkin(11) * t833;
t792 = -pkin(3) * t806 - t807 * t890;
t816 = t830 * t838 - t834 * t865;
t837 = cos(pkin(8));
t858 = t806 * t837 + t816 * t833;
t795 = t858 * pkin(11);
t801 = t836 * t825 + t832 * t857;
t808 = -t826 * t834 + t829 * t838;
t889 = pkin(11) * t837;
t736 = pkin(3) * t808 - t792 * t807 + t795 * t816 - t801 * t889 + t752;
t883 = t832 * t838;
t753 = 0.2e1 * qJD(3) * t806 + t778 * t883 + t836 * t779 + t783 * t884;
t797 = pkin(3) * t816 - t807 * t889;
t800 = -t832 * t825 + t836 * t857;
t859 = t800 * t837 + t808 * t833;
t737 = pkin(11) * t859 + t806 * t792 - t816 * t797 + t753;
t767 = -t778 * t834 + t838 * t783 + qJDD(3);
t741 = -pkin(3) * t800 - t795 * t806 + t797 * t807 - t801 * t890 + t767;
t842 = sin(qJ(4));
t847 = cos(qJ(4));
t723 = -t842 * t737 + (t736 * t837 + t741 * t833) * t847;
t787 = t847 * t807 + t842 * t858;
t765 = -t787 * qJD(4) - t842 * t801 + t847 * t859;
t786 = -t842 * t807 + t847 * t858;
t766 = t786 * qJD(4) + t847 * t801 + t842 * t859;
t768 = -mrSges(5,1) * t786 + mrSges(5,2) * t787;
t796 = -t806 * t833 + t816 * t837 + qJD(4);
t773 = -mrSges(5,2) * t796 + mrSges(5,3) * t786;
t791 = -t800 * t833 + t808 * t837 + qJDD(4);
t877 = t837 * t842;
t882 = t833 * t842;
t724 = t736 * t877 + t847 * t737 + t741 * t882;
t769 = -pkin(4) * t786 - pkin(12) * t787;
t794 = t796 ^ 2;
t720 = -pkin(4) * t794 + pkin(12) * t791 + t769 * t786 + t724;
t725 = -t833 * t736 + t837 * t741;
t722 = (-t786 * t796 - t766) * pkin(12) + (t787 * t796 - t765) * pkin(4) + t725;
t841 = sin(qJ(5));
t846 = cos(qJ(5));
t716 = t846 * t720 + t841 * t722;
t771 = -t787 * t841 + t796 * t846;
t772 = t787 * t846 + t796 * t841;
t755 = -pkin(5) * t771 - pkin(13) * t772;
t764 = qJDD(5) - t765;
t785 = qJD(5) - t786;
t784 = t785 ^ 2;
t714 = -pkin(5) * t784 + pkin(13) * t764 + t755 * t771 + t716;
t719 = -t791 * pkin(4) - t794 * pkin(12) + t787 * t769 - t723;
t744 = -qJD(5) * t772 - t766 * t841 + t791 * t846;
t745 = qJD(5) * t771 + t766 * t846 + t791 * t841;
t717 = (-t771 * t785 - t745) * pkin(13) + (t772 * t785 - t744) * pkin(5) + t719;
t840 = sin(qJ(6));
t845 = cos(qJ(6));
t711 = -t714 * t840 + t717 * t845;
t757 = -t772 * t840 + t785 * t845;
t728 = qJD(6) * t757 + t745 * t845 + t764 * t840;
t758 = t772 * t845 + t785 * t840;
t738 = -mrSges(7,1) * t757 + mrSges(7,2) * t758;
t743 = qJDD(6) - t744;
t770 = qJD(6) - t771;
t746 = -mrSges(7,2) * t770 + mrSges(7,3) * t757;
t709 = m(7) * t711 + mrSges(7,1) * t743 - mrSges(7,3) * t728 - t738 * t758 + t746 * t770;
t712 = t714 * t845 + t717 * t840;
t727 = -qJD(6) * t758 - t745 * t840 + t764 * t845;
t747 = mrSges(7,1) * t770 - mrSges(7,3) * t758;
t710 = m(7) * t712 - mrSges(7,2) * t743 + mrSges(7,3) * t727 + t738 * t757 - t747 * t770;
t703 = t709 * t845 + t710 * t840;
t759 = -mrSges(6,2) * t785 + mrSges(6,3) * t771;
t760 = mrSges(6,1) * t785 - mrSges(6,3) * t772;
t851 = -m(6) * t719 + t744 * mrSges(6,1) - mrSges(6,2) * t745 + t771 * t759 - t760 * t772 - t703;
t699 = m(5) * t723 + mrSges(5,1) * t791 - mrSges(5,3) * t766 - t768 * t787 + t773 * t796 + t851;
t885 = t699 * t847;
t880 = t835 * t843;
t879 = t835 * t848;
t774 = mrSges(5,1) * t796 - mrSges(5,3) * t787;
t754 = -mrSges(6,1) * t771 + mrSges(6,2) * t772;
t862 = -t709 * t840 + t845 * t710;
t702 = m(6) * t716 - mrSges(6,2) * t764 + mrSges(6,3) * t744 + t754 * t771 - t760 * t785 + t862;
t715 = -t720 * t841 + t722 * t846;
t713 = -pkin(5) * t764 - pkin(13) * t784 + t755 * t772 - t715;
t852 = -m(7) * t713 + t727 * mrSges(7,1) - mrSges(7,2) * t728 + t757 * t746 - t747 * t758;
t707 = m(6) * t715 + mrSges(6,1) * t764 - mrSges(6,3) * t745 - t754 * t772 + t759 * t785 + t852;
t863 = t846 * t702 - t707 * t841;
t693 = m(5) * t724 - mrSges(5,2) * t791 + mrSges(5,3) * t765 + t768 * t786 - t774 * t796 + t863;
t696 = t841 * t702 + t846 * t707;
t695 = m(5) * t725 - mrSges(5,1) * t765 + mrSges(5,2) * t766 - t773 * t786 + t774 * t787 + t696;
t682 = t693 * t877 - t695 * t833 + t837 * t885;
t793 = -mrSges(4,1) * t806 + mrSges(4,2) * t807;
t798 = -mrSges(4,2) * t816 + mrSges(4,3) * t806;
t678 = m(4) * t752 + mrSges(4,1) * t808 - mrSges(4,3) * t801 - t793 * t807 + t798 * t816 + t682;
t681 = t693 * t882 + t837 * t695 + t833 * t885;
t799 = mrSges(4,1) * t816 - mrSges(4,3) * t807;
t680 = m(4) * t767 - mrSges(4,1) * t800 + mrSges(4,2) * t801 - t798 * t806 + t799 * t807 + t681;
t687 = t847 * t693 - t699 * t842;
t686 = m(4) * t753 - mrSges(4,2) * t808 + mrSges(4,3) * t800 + t793 * t806 - t799 * t816 + t687;
t667 = t678 * t878 - t680 * t834 + t686 * t883;
t802 = -g(3) * t879 + t861;
t821 = -mrSges(3,2) * t830 + mrSges(3,3) * t865;
t824 = (-mrSges(3,1) * t848 + mrSges(3,2) * t843) * t870;
t663 = m(3) * t802 + mrSges(3,1) * t829 - mrSges(3,3) * t825 + t821 * t830 - t824 * t866 + t667;
t666 = t678 * t881 + t838 * t680 + t686 * t884;
t812 = -t835 * t822 - t888;
t820 = mrSges(3,1) * t830 - mrSges(3,3) * t866;
t665 = m(3) * t812 - t826 * mrSges(3,1) + t825 * mrSges(3,2) + (t820 * t843 - t821 * t848) * t870 + t666;
t672 = -t678 * t832 + t836 * t686;
t803 = -g(3) * t880 + t871;
t671 = m(3) * t803 - mrSges(3,2) * t829 + mrSges(3,3) * t826 - t820 * t830 + t824 * t865 + t672;
t653 = t663 * t874 - t665 * t835 + t671 * t875;
t651 = m(2) * t827 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t850 + t653;
t660 = -t663 * t843 + t848 * t671;
t659 = m(2) * t828 - mrSges(2,1) * t850 - qJDD(1) * mrSges(2,2) + t660;
t872 = t849 * t651 + t844 * t659;
t652 = t663 * t879 + t839 * t665 + t671 * t880;
t864 = -t651 * t844 + t849 * t659;
t729 = Ifges(7,5) * t758 + Ifges(7,6) * t757 + Ifges(7,3) * t770;
t731 = Ifges(7,1) * t758 + Ifges(7,4) * t757 + Ifges(7,5) * t770;
t704 = -mrSges(7,1) * t713 + mrSges(7,3) * t712 + Ifges(7,4) * t728 + Ifges(7,2) * t727 + Ifges(7,6) * t743 - t729 * t758 + t731 * t770;
t730 = Ifges(7,4) * t758 + Ifges(7,2) * t757 + Ifges(7,6) * t770;
t705 = mrSges(7,2) * t713 - mrSges(7,3) * t711 + Ifges(7,1) * t728 + Ifges(7,4) * t727 + Ifges(7,5) * t743 + t729 * t757 - t730 * t770;
t748 = Ifges(6,5) * t772 + Ifges(6,6) * t771 + Ifges(6,3) * t785;
t749 = Ifges(6,4) * t772 + Ifges(6,2) * t771 + Ifges(6,6) * t785;
t688 = mrSges(6,2) * t719 - mrSges(6,3) * t715 + Ifges(6,1) * t745 + Ifges(6,4) * t744 + Ifges(6,5) * t764 - pkin(13) * t703 - t704 * t840 + t705 * t845 + t748 * t771 - t749 * t785;
t750 = Ifges(6,1) * t772 + Ifges(6,4) * t771 + Ifges(6,5) * t785;
t689 = -mrSges(6,1) * t719 - mrSges(7,1) * t711 + mrSges(7,2) * t712 + mrSges(6,3) * t716 + Ifges(6,4) * t745 - Ifges(7,5) * t728 + Ifges(6,2) * t744 + Ifges(6,6) * t764 - Ifges(7,6) * t727 - Ifges(7,3) * t743 - pkin(5) * t703 - t730 * t758 + t731 * t757 - t748 * t772 + t750 * t785;
t762 = Ifges(5,4) * t787 + Ifges(5,2) * t786 + Ifges(5,6) * t796;
t763 = Ifges(5,1) * t787 + Ifges(5,4) * t786 + Ifges(5,5) * t796;
t673 = mrSges(5,1) * t723 - mrSges(5,2) * t724 + Ifges(5,5) * t766 + Ifges(5,6) * t765 + Ifges(5,3) * t791 + pkin(4) * t851 + pkin(12) * t863 + t841 * t688 + t846 * t689 + t787 * t762 - t786 * t763;
t789 = Ifges(4,4) * t807 + Ifges(4,2) * t806 + Ifges(4,6) * t816;
t790 = Ifges(4,1) * t807 + Ifges(4,4) * t806 + Ifges(4,5) * t816;
t761 = Ifges(5,5) * t787 + Ifges(5,6) * t786 + Ifges(5,3) * t796;
t674 = mrSges(5,2) * t725 - mrSges(5,3) * t723 + Ifges(5,1) * t766 + Ifges(5,4) * t765 + Ifges(5,5) * t791 - pkin(12) * t696 + t688 * t846 - t689 * t841 + t761 * t786 - t762 * t796;
t675 = Ifges(5,4) * t766 + Ifges(5,2) * t765 + Ifges(5,6) * t791 - t787 * t761 + t796 * t763 - mrSges(5,1) * t725 + mrSges(5,3) * t724 - Ifges(6,5) * t745 - Ifges(6,6) * t744 - Ifges(6,3) * t764 - t772 * t749 + t771 * t750 - mrSges(6,1) * t715 + mrSges(6,2) * t716 - t840 * t705 - t845 * t704 - pkin(5) * t852 - pkin(13) * t862 - pkin(4) * t696;
t854 = pkin(11) * t687 + t674 * t842 + t675 * t847;
t654 = mrSges(4,1) * t752 - mrSges(4,2) * t753 + Ifges(4,5) * t801 + Ifges(4,6) * t800 + Ifges(4,3) * t808 + pkin(3) * t682 + t837 * t673 + t807 * t789 - t806 * t790 + t833 * t854;
t809 = Ifges(3,3) * t830 + (Ifges(3,5) * t843 + Ifges(3,6) * t848) * t870;
t811 = Ifges(3,5) * t830 + (Ifges(3,1) * t843 + Ifges(3,4) * t848) * t870;
t788 = Ifges(4,5) * t807 + Ifges(4,6) * t806 + Ifges(4,3) * t816;
t655 = -mrSges(4,1) * t767 + mrSges(4,3) * t753 + Ifges(4,4) * t801 + Ifges(4,2) * t800 + Ifges(4,6) * t808 - pkin(3) * t681 - t833 * t673 - t807 * t788 + t816 * t790 + t837 * t854;
t656 = mrSges(4,2) * t767 - mrSges(4,3) * t752 + Ifges(4,1) * t801 + Ifges(4,4) * t800 + Ifges(4,5) * t808 + t847 * t674 - t842 * t675 + t806 * t788 - t816 * t789 + (-t681 * t833 - t682 * t837) * pkin(11);
t853 = qJ(3) * t672 + t655 * t836 + t656 * t832;
t648 = -mrSges(3,1) * t812 + mrSges(3,3) * t803 + Ifges(3,4) * t825 + Ifges(3,2) * t826 + Ifges(3,6) * t829 - pkin(2) * t666 - t834 * t654 - t809 * t866 + t830 * t811 + t838 * t853;
t810 = Ifges(3,6) * t830 + (Ifges(3,4) * t843 + Ifges(3,2) * t848) * t870;
t649 = t809 * t865 + mrSges(3,2) * t812 - mrSges(3,3) * t802 + Ifges(3,1) * t825 + Ifges(3,4) * t826 + Ifges(3,5) * t829 - t832 * t655 + t836 * t656 - t830 * t810 + (-t666 * t834 - t667 * t838) * qJ(3);
t855 = pkin(10) * t660 + t648 * t848 + t649 * t843;
t647 = mrSges(3,1) * t802 - mrSges(3,2) * t803 + Ifges(3,5) * t825 + Ifges(3,6) * t826 + Ifges(3,3) * t829 + pkin(2) * t667 + t838 * t654 + (t810 * t843 - t811 * t848) * t870 + t853 * t834;
t646 = -mrSges(2,2) * g(3) - mrSges(2,3) * t827 + Ifges(2,5) * qJDD(1) - t850 * Ifges(2,6) - t843 * t648 + t848 * t649 + (-t652 * t835 - t653 * t839) * pkin(10);
t645 = mrSges(2,1) * g(3) + mrSges(2,3) * t828 + t850 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t652 - t835 * t647 + t839 * t855;
t1 = [-m(1) * g(1) + t864; -m(1) * g(2) + t872; (-m(1) - m(2)) * g(3) + t652; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(9) * t872 - t844 * t645 + t849 * t646; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(9) * t864 + t849 * t645 + t844 * t646; -mrSges(1,1) * g(2) + mrSges(2,1) * t827 + mrSges(1,2) * g(1) - mrSges(2,2) * t828 + Ifges(2,3) * qJDD(1) + pkin(1) * t653 + t839 * t647 + t835 * t855;];
tauB  = t1;
