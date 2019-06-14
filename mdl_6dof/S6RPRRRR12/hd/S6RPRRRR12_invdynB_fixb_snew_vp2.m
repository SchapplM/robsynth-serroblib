% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 07:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_invdynB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 06:32:42
% EndTime: 2019-05-06 06:37:41
% DurationCPUTime: 299.72s
% Computational Cost: add. (4386647->418), mult. (14307186->580), div. (0->0), fcn. (12579308->18), ass. (0->198)
t833 = sin(pkin(7));
t835 = cos(pkin(14));
t838 = cos(pkin(6));
t834 = sin(pkin(6));
t837 = cos(pkin(7));
t884 = t834 * t837;
t818 = (t833 * t838 + t835 * t884) * qJD(1) * pkin(10);
t843 = sin(qJ(1));
t848 = cos(qJ(1));
t829 = -g(1) * t848 - g(2) * t843;
t849 = qJD(1) ^ 2;
t892 = qJ(2) * t834;
t822 = -pkin(1) * t849 + qJDD(1) * t892 + t829;
t831 = sin(pkin(14));
t897 = pkin(10) * t831;
t863 = -pkin(2) * t835 - t833 * t897;
t877 = qJD(1) * t834;
t893 = pkin(10) * qJDD(1);
t859 = qJD(1) * t863 * t877 + t837 * t893;
t828 = t843 * g(1) - g(2) * t848;
t821 = qJDD(1) * pkin(1) + t849 * t892 + t828;
t875 = qJD(2) * t877;
t883 = t835 * t838;
t885 = t834 * t835;
t864 = -g(3) * t885 + t821 * t883 - 0.2e1 * t831 * t875;
t781 = (pkin(2) * qJDD(1) + qJD(1) * t818) * t838 + (-t859 * t834 - t822) * t831 + t864;
t823 = (pkin(2) * t838 - t884 * t897) * qJD(1);
t889 = t831 * t838;
t876 = t821 * t889 + (t822 + 0.2e1 * t875) * t835;
t782 = (-qJD(1) * t823 + t833 * t893) * t838 + (-g(3) * t831 + t859 * t835) * t834 + t876;
t874 = -g(3) * t838 + qJDD(2);
t790 = (-t821 + t863 * qJDD(1) + (-t818 * t835 + t823 * t831) * qJD(1)) * t834 + t874;
t842 = sin(qJ(3));
t847 = cos(qJ(3));
t880 = t837 * t847;
t886 = t833 * t847;
t755 = t781 * t880 - t782 * t842 + t790 * t886;
t851 = t838 * t886 + (-t831 * t842 + t835 * t880) * t834;
t807 = t851 * qJD(1);
t881 = t837 * t842;
t887 = t833 * t842;
t852 = t838 * t887 + (t831 * t847 + t835 * t881) * t834;
t808 = t852 * qJD(1);
t832 = sin(pkin(8));
t896 = pkin(11) * t832;
t795 = -pkin(3) * t807 - t808 * t896;
t798 = qJD(3) * t807 + t852 * qJDD(1);
t860 = -t833 * t885 + t837 * t838;
t819 = t860 * qJD(1) + qJD(3);
t836 = cos(pkin(8));
t866 = t807 * t836 + t819 * t832;
t800 = t866 * pkin(11);
t816 = t860 * qJDD(1) + qJDD(3);
t895 = pkin(11) * t836;
t739 = pkin(3) * t816 - t795 * t808 - t798 * t895 + t800 * t819 + t755;
t756 = t781 * t881 + t847 * t782 + t790 * t887;
t804 = pkin(3) * t819 - t808 * t895;
t797 = -qJD(3) * t808 + t851 * qJDD(1);
t867 = t797 * t836 + t816 * t832;
t740 = t867 * pkin(11) + t795 * t807 - t804 * t819 + t756;
t770 = -t781 * t833 + t837 * t790;
t743 = -pkin(3) * t797 - t798 * t896 - t800 * t807 + t804 * t808 + t770;
t841 = sin(qJ(4));
t846 = cos(qJ(4));
t726 = -t841 * t740 + (t739 * t836 + t743 * t832) * t846;
t789 = t846 * t808 + t866 * t841;
t763 = -qJD(4) * t789 - t798 * t841 + t867 * t846;
t788 = -t841 * t808 + t866 * t846;
t894 = Ifges(3,3) * t838;
t764 = qJD(4) * t788 + t798 * t846 + t867 * t841;
t771 = -mrSges(5,1) * t788 + mrSges(5,2) * t789;
t801 = -t807 * t832 + t819 * t836 + qJD(4);
t776 = -mrSges(5,2) * t801 + mrSges(5,3) * t788;
t791 = -t797 * t832 + t816 * t836 + qJDD(4);
t882 = t836 * t841;
t888 = t832 * t841;
t727 = t739 * t882 + t846 * t740 + t743 * t888;
t772 = -pkin(4) * t788 - pkin(12) * t789;
t799 = t801 ^ 2;
t723 = -pkin(4) * t799 + pkin(12) * t791 + t772 * t788 + t727;
t728 = -t739 * t832 + t836 * t743;
t725 = (-t788 * t801 - t764) * pkin(12) + (t789 * t801 - t763) * pkin(4) + t728;
t840 = sin(qJ(5));
t845 = cos(qJ(5));
t719 = t845 * t723 + t840 * t725;
t774 = -t789 * t840 + t801 * t845;
t775 = t789 * t845 + t801 * t840;
t758 = -pkin(5) * t774 - pkin(13) * t775;
t762 = qJDD(5) - t763;
t786 = qJD(5) - t788;
t785 = t786 ^ 2;
t717 = -pkin(5) * t785 + pkin(13) * t762 + t758 * t774 + t719;
t722 = -pkin(4) * t791 - pkin(12) * t799 + t789 * t772 - t726;
t747 = -qJD(5) * t775 - t764 * t840 + t791 * t845;
t748 = qJD(5) * t774 + t764 * t845 + t791 * t840;
t720 = (-t774 * t786 - t748) * pkin(13) + (t775 * t786 - t747) * pkin(5) + t722;
t839 = sin(qJ(6));
t844 = cos(qJ(6));
t714 = -t717 * t839 + t720 * t844;
t760 = -t775 * t839 + t786 * t844;
t731 = qJD(6) * t760 + t748 * t844 + t762 * t839;
t761 = t775 * t844 + t786 * t839;
t744 = -mrSges(7,1) * t760 + mrSges(7,2) * t761;
t746 = qJDD(6) - t747;
t773 = qJD(6) - t774;
t749 = -mrSges(7,2) * t773 + mrSges(7,3) * t760;
t712 = m(7) * t714 + mrSges(7,1) * t746 - mrSges(7,3) * t731 - t744 * t761 + t749 * t773;
t715 = t717 * t844 + t720 * t839;
t730 = -qJD(6) * t761 - t748 * t839 + t762 * t844;
t750 = mrSges(7,1) * t773 - mrSges(7,3) * t761;
t713 = m(7) * t715 - mrSges(7,2) * t746 + mrSges(7,3) * t730 + t744 * t760 - t750 * t773;
t706 = t712 * t844 + t713 * t839;
t765 = -mrSges(6,2) * t786 + mrSges(6,3) * t774;
t766 = mrSges(6,1) * t786 - mrSges(6,3) * t775;
t850 = -m(6) * t722 + t747 * mrSges(6,1) - mrSges(6,2) * t748 + t774 * t765 - t766 * t775 - t706;
t702 = m(5) * t726 + mrSges(5,1) * t791 - mrSges(5,3) * t764 - t771 * t789 + t776 * t801 + t850;
t891 = t702 * t846;
t890 = t831 * t834;
t777 = mrSges(5,1) * t801 - mrSges(5,3) * t789;
t757 = -mrSges(6,1) * t774 + mrSges(6,2) * t775;
t871 = -t712 * t839 + t844 * t713;
t705 = m(6) * t719 - mrSges(6,2) * t762 + mrSges(6,3) * t747 + t757 * t774 - t766 * t786 + t871;
t718 = -t723 * t840 + t725 * t845;
t716 = -pkin(5) * t762 - pkin(13) * t785 + t758 * t775 - t718;
t855 = -m(7) * t716 + t730 * mrSges(7,1) - mrSges(7,2) * t731 + t760 * t749 - t750 * t761;
t710 = m(6) * t718 + mrSges(6,1) * t762 - mrSges(6,3) * t748 - t757 * t775 + t765 * t786 + t855;
t872 = t845 * t705 - t710 * t840;
t696 = m(5) * t727 - mrSges(5,2) * t791 + mrSges(5,3) * t763 + t771 * t788 - t777 * t801 + t872;
t699 = t840 * t705 + t845 * t710;
t698 = m(5) * t728 - mrSges(5,1) * t763 + mrSges(5,2) * t764 - t776 * t788 + t777 * t789 + t699;
t685 = t696 * t882 - t698 * t832 + t836 * t891;
t796 = -mrSges(4,1) * t807 + mrSges(4,2) * t808;
t805 = -mrSges(4,2) * t819 + mrSges(4,3) * t807;
t681 = m(4) * t755 + mrSges(4,1) * t816 - mrSges(4,3) * t798 - t796 * t808 + t805 * t819 + t685;
t684 = t696 * t888 + t836 * t698 + t832 * t891;
t806 = mrSges(4,1) * t819 - mrSges(4,3) * t808;
t683 = m(4) * t770 - mrSges(4,1) * t797 + mrSges(4,2) * t798 - t805 * t807 + t806 * t808 + t684;
t690 = t846 * t696 - t702 * t841;
t689 = m(4) * t756 - mrSges(4,2) * t816 + mrSges(4,3) * t797 + t796 * t807 - t806 * t819 + t690;
t670 = t681 * t880 - t683 * t833 + t689 * t881;
t802 = -t822 * t831 + t864;
t870 = -mrSges(3,1) * t835 + mrSges(3,2) * t831;
t820 = t870 * t877;
t861 = -mrSges(3,2) * t838 + mrSges(3,3) * t885;
t825 = t861 * qJD(1);
t862 = mrSges(3,1) * t838 - mrSges(3,3) * t890;
t666 = m(3) * t802 + t862 * qJDD(1) + (-t820 * t890 + t825 * t838) * qJD(1) + t670;
t669 = t681 * t886 + t837 * t683 + t689 * t887;
t809 = -t821 * t834 + t874;
t824 = t862 * qJD(1);
t668 = m(3) * t809 + (t870 * qJDD(1) + (t824 * t831 - t825 * t835) * qJD(1)) * t834 + t669;
t675 = -t681 * t842 + t847 * t689;
t803 = -g(3) * t890 + t876;
t674 = m(3) * t803 + t861 * qJDD(1) + (t820 * t885 - t824 * t838) * qJD(1) + t675;
t656 = t666 * t883 - t668 * t834 + t674 * t889;
t654 = m(2) * t828 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t849 + t656;
t663 = -t666 * t831 + t835 * t674;
t662 = m(2) * t829 - mrSges(2,1) * t849 - qJDD(1) * mrSges(2,2) + t663;
t878 = t848 * t654 + t843 * t662;
t655 = t666 * t885 + t838 * t668 + t674 * t890;
t873 = -t843 * t654 + t848 * t662;
t869 = Ifges(3,5) * t831 + Ifges(3,6) * t835;
t732 = Ifges(7,5) * t761 + Ifges(7,6) * t760 + Ifges(7,3) * t773;
t734 = Ifges(7,1) * t761 + Ifges(7,4) * t760 + Ifges(7,5) * t773;
t707 = -mrSges(7,1) * t716 + mrSges(7,3) * t715 + Ifges(7,4) * t731 + Ifges(7,2) * t730 + Ifges(7,6) * t746 - t732 * t761 + t734 * t773;
t733 = Ifges(7,4) * t761 + Ifges(7,2) * t760 + Ifges(7,6) * t773;
t708 = mrSges(7,2) * t716 - mrSges(7,3) * t714 + Ifges(7,1) * t731 + Ifges(7,4) * t730 + Ifges(7,5) * t746 + t732 * t760 - t733 * t773;
t751 = Ifges(6,5) * t775 + Ifges(6,6) * t774 + Ifges(6,3) * t786;
t752 = Ifges(6,4) * t775 + Ifges(6,2) * t774 + Ifges(6,6) * t786;
t691 = mrSges(6,2) * t722 - mrSges(6,3) * t718 + Ifges(6,1) * t748 + Ifges(6,4) * t747 + Ifges(6,5) * t762 - pkin(13) * t706 - t707 * t839 + t708 * t844 + t751 * t774 - t752 * t786;
t753 = Ifges(6,1) * t775 + Ifges(6,4) * t774 + Ifges(6,5) * t786;
t692 = -mrSges(6,1) * t722 - mrSges(7,1) * t714 + mrSges(7,2) * t715 + mrSges(6,3) * t719 + Ifges(6,4) * t748 - Ifges(7,5) * t731 + Ifges(6,2) * t747 + Ifges(6,6) * t762 - Ifges(7,6) * t730 - Ifges(7,3) * t746 - pkin(5) * t706 - t733 * t761 + t734 * t760 - t751 * t775 + t753 * t786;
t768 = Ifges(5,4) * t789 + Ifges(5,2) * t788 + Ifges(5,6) * t801;
t769 = Ifges(5,1) * t789 + Ifges(5,4) * t788 + Ifges(5,5) * t801;
t676 = mrSges(5,1) * t726 - mrSges(5,2) * t727 + Ifges(5,5) * t764 + Ifges(5,6) * t763 + Ifges(5,3) * t791 + pkin(4) * t850 + pkin(12) * t872 + t840 * t691 + t845 * t692 + t789 * t768 - t788 * t769;
t792 = Ifges(4,5) * t808 + Ifges(4,6) * t807 + Ifges(4,3) * t819;
t794 = Ifges(4,1) * t808 + Ifges(4,4) * t807 + Ifges(4,5) * t819;
t767 = Ifges(5,5) * t789 + Ifges(5,6) * t788 + Ifges(5,3) * t801;
t677 = mrSges(5,2) * t728 - mrSges(5,3) * t726 + Ifges(5,1) * t764 + Ifges(5,4) * t763 + Ifges(5,5) * t791 - pkin(12) * t699 + t691 * t845 - t692 * t840 + t767 * t788 - t768 * t801;
t678 = Ifges(5,4) * t764 + Ifges(5,2) * t763 + Ifges(5,6) * t791 - t789 * t767 + t801 * t769 - mrSges(5,1) * t728 + mrSges(5,3) * t727 - Ifges(6,5) * t748 - Ifges(6,6) * t747 - Ifges(6,3) * t762 - t775 * t752 + t774 * t753 - mrSges(6,1) * t718 + mrSges(6,2) * t719 - t839 * t708 - t844 * t707 - pkin(5) * t855 - pkin(13) * t871 - pkin(4) * t699;
t857 = pkin(11) * t690 + t677 * t841 + t678 * t846;
t658 = -mrSges(4,1) * t770 + mrSges(4,3) * t756 + Ifges(4,4) * t798 + Ifges(4,2) * t797 + Ifges(4,6) * t816 - pkin(3) * t684 - t676 * t832 - t792 * t808 + t794 * t819 + t857 * t836;
t793 = Ifges(4,4) * t808 + Ifges(4,2) * t807 + Ifges(4,6) * t819;
t659 = mrSges(4,2) * t770 - mrSges(4,3) * t755 + Ifges(4,1) * t798 + Ifges(4,4) * t797 + Ifges(4,5) * t816 + t677 * t846 - t678 * t841 + t792 * t807 - t793 * t819 + (-t684 * t832 - t685 * t836) * pkin(11);
t858 = pkin(10) * t675 + t658 * t847 + t659 * t842;
t657 = mrSges(4,1) * t755 - mrSges(4,2) * t756 + Ifges(4,5) * t798 + Ifges(4,6) * t797 + Ifges(4,3) * t816 + pkin(3) * t685 + t676 * t836 + t793 * t808 - t794 * t807 + t857 * t832;
t812 = (t869 * t834 + t894) * qJD(1);
t854 = Ifges(3,5) * t838 + (Ifges(3,1) * t831 + Ifges(3,4) * t835) * t834;
t814 = t854 * qJD(1);
t853 = Ifges(3,6) * t838 + (Ifges(3,4) * t831 + Ifges(3,2) * t835) * t834;
t651 = -mrSges(3,1) * t809 + mrSges(3,3) * t803 - pkin(2) * t669 - t657 * t833 + (-t812 * t890 + t814 * t838) * qJD(1) + t858 * t837 + t853 * qJDD(1);
t813 = t853 * qJD(1);
t652 = mrSges(3,2) * t809 - mrSges(3,3) * t802 - t658 * t842 + t659 * t847 + (t812 * t885 - t813 * t838) * qJD(1) + (-t669 * t833 - t670 * t837) * pkin(10) + t854 * qJDD(1);
t856 = qJ(2) * t663 + t651 * t835 + t652 * t831;
t650 = qJDD(1) * t894 + mrSges(3,1) * t802 - mrSges(3,2) * t803 + pkin(2) * t670 + t657 * t837 + t858 * t833 + (t869 * qJDD(1) + (t813 * t831 - t814 * t835) * qJD(1)) * t834;
t649 = -mrSges(2,2) * g(3) - mrSges(2,3) * t828 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t849 - t651 * t831 + t652 * t835 + (-t655 * t834 - t656 * t838) * qJ(2);
t648 = mrSges(2,1) * g(3) + mrSges(2,3) * t829 + Ifges(2,5) * t849 + Ifges(2,6) * qJDD(1) - pkin(1) * t655 - t650 * t834 + t856 * t838;
t1 = [-m(1) * g(1) + t873; -m(1) * g(2) + t878; (-m(1) - m(2)) * g(3) + t655; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(9) * t878 - t843 * t648 + t848 * t649; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(9) * t873 + t848 * t648 + t843 * t649; -mrSges(1,1) * g(2) + mrSges(2,1) * t828 + mrSges(1,2) * g(1) - mrSges(2,2) * t829 + Ifges(2,3) * qJDD(1) + pkin(1) * t656 + t650 * t838 + t856 * t834;];
tauB  = t1;
