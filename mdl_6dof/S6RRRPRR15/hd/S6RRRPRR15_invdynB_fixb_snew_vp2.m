% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 17:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR15_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR15_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:08:14
% EndTime: 2019-05-07 17:09:02
% DurationCPUTime: 45.92s
% Computational Cost: add. (748233->403), mult. (1873715->521), div. (0->0), fcn. (1555760->14), ass. (0->176)
t885 = Ifges(4,5) - Ifges(5,4);
t894 = -Ifges(4,2) - Ifges(5,3);
t893 = Ifges(5,2) + Ifges(4,1);
t884 = Ifges(4,6) - Ifges(5,5);
t883 = -Ifges(5,6) - Ifges(4,4);
t892 = Ifges(4,3) + Ifges(5,1);
t838 = cos(qJ(2));
t829 = sin(pkin(7));
t882 = cos(pkin(6));
t851 = t882 * qJD(1) + qJD(2);
t847 = t851 * t829;
t830 = sin(pkin(6));
t870 = qJD(1) * t830;
t881 = cos(pkin(7));
t852 = t881 * t870;
t891 = t838 * t852 + t847;
t807 = t891 * pkin(10);
t834 = sin(qJ(2));
t887 = pkin(10) * t834;
t813 = (-pkin(2) * t838 - t829 * t887) * t870;
t866 = qJD(1) * qJD(2);
t819 = (qJDD(1) * t834 + t838 * t866) * t830;
t827 = t882 * qJDD(1) + qJDD(2);
t835 = sin(qJ(1));
t839 = cos(qJ(1));
t825 = t835 * g(1) - g(2) * t839;
t840 = qJD(1) ^ 2;
t888 = pkin(9) * t830;
t816 = qJDD(1) * pkin(1) + t840 * t888 + t825;
t826 = -g(1) * t839 - g(2) * t835;
t817 = -pkin(1) * t840 + qJDD(1) * t888 + t826;
t858 = t838 * t882;
t854 = t816 * t858 - t834 * t817;
t863 = t881 * pkin(10);
t869 = qJD(1) * t834;
t748 = -t819 * t863 + t827 * pkin(2) + t851 * t807 + (-g(3) * t838 - t813 * t869) * t830 + t854;
t812 = t851 * pkin(2) - t852 * t887;
t820 = (qJDD(1) * t838 - t834 * t866) * t830;
t849 = t881 * t820 + t827 * t829;
t868 = qJD(1) * t838;
t859 = t834 * t882;
t871 = t816 * t859 + t838 * t817;
t749 = -t851 * t812 + (-g(3) * t834 + t813 * t868) * t830 + t849 * pkin(10) + t871;
t864 = t882 * g(3);
t754 = -t819 * t829 * pkin(10) - t864 - t820 * pkin(2) + (-t816 + (-t807 * t838 + t812 * t834) * qJD(1)) * t830;
t833 = sin(qJ(3));
t860 = t833 * t881;
t879 = t829 * t833;
t890 = cos(qJ(3));
t728 = t748 * t860 + t890 * t749 + t754 * t879;
t862 = t830 * t869;
t794 = t833 * t862 - t891 * t890;
t795 = t833 * t847 + (t890 * t834 + t838 * t860) * t870;
t776 = pkin(3) * t794 - qJ(4) * t795;
t798 = -t829 * t820 + t881 * t827 + qJDD(3);
t861 = t830 * t868;
t809 = t829 * t861 - t881 * t851 - qJD(3);
t806 = t809 ^ 2;
t722 = t806 * pkin(3) - t798 * qJ(4) + 0.2e1 * qJD(4) * t809 + t794 * t776 - t728;
t733 = -t829 * t748 + t881 * t754;
t853 = t881 * t890;
t865 = t829 * t890;
t772 = qJD(3) * t795 + t819 * t833 - t820 * t853 - t827 * t865;
t773 = -t794 * qJD(3) + t890 * t819 + t849 * t833;
t785 = -mrSges(4,1) * t809 - mrSges(4,3) * t795;
t880 = t794 * t809;
t842 = (-t773 - t880) * qJ(4) + t733 + (-t809 * pkin(3) - 0.2e1 * qJD(4)) * t795;
t721 = t772 * pkin(3) + t842;
t784 = mrSges(5,1) * t795 - mrSges(5,2) * t809;
t727 = t748 * t853 - t833 * t749 + t754 * t865;
t723 = -t798 * pkin(3) - t806 * qJ(4) + t795 * t776 + qJDD(4) - t727;
t716 = (t794 * t795 - t798) * pkin(11) + (t773 - t880) * pkin(4) + t723;
t786 = pkin(4) * t795 + pkin(11) * t809;
t793 = t794 ^ 2;
t717 = -t793 * pkin(4) - t795 * t786 + (pkin(3) + pkin(11)) * t772 + t842;
t832 = sin(qJ(5));
t837 = cos(qJ(5));
t712 = t832 * t716 + t837 * t717;
t781 = t794 * t832 - t809 * t837;
t737 = -qJD(5) * t781 + t772 * t837 - t798 * t832;
t780 = t794 * t837 + t809 * t832;
t750 = -mrSges(6,1) * t780 + mrSges(6,2) * t781;
t792 = qJD(5) + t795;
t761 = mrSges(6,1) * t792 - mrSges(6,3) * t781;
t771 = qJDD(5) + t773;
t751 = -pkin(5) * t780 - pkin(12) * t781;
t791 = t792 ^ 2;
t710 = -pkin(5) * t791 + pkin(12) * t771 + t751 * t780 + t712;
t719 = -t772 * pkin(4) - t793 * pkin(11) - t809 * t786 - t722;
t738 = qJD(5) * t780 + t772 * t832 + t798 * t837;
t713 = (-t780 * t792 - t738) * pkin(12) + (t781 * t792 - t737) * pkin(5) + t719;
t831 = sin(qJ(6));
t836 = cos(qJ(6));
t707 = -t710 * t831 + t713 * t836;
t758 = -t781 * t831 + t792 * t836;
t726 = qJD(6) * t758 + t738 * t836 + t771 * t831;
t759 = t781 * t836 + t792 * t831;
t734 = -mrSges(7,1) * t758 + mrSges(7,2) * t759;
t736 = qJDD(6) - t737;
t779 = qJD(6) - t780;
t739 = -mrSges(7,2) * t779 + mrSges(7,3) * t758;
t705 = m(7) * t707 + mrSges(7,1) * t736 - mrSges(7,3) * t726 - t734 * t759 + t739 * t779;
t708 = t710 * t836 + t713 * t831;
t725 = -qJD(6) * t759 - t738 * t831 + t771 * t836;
t740 = mrSges(7,1) * t779 - mrSges(7,3) * t759;
t706 = m(7) * t708 - mrSges(7,2) * t736 + mrSges(7,3) * t725 + t734 * t758 - t740 * t779;
t855 = -t705 * t831 + t836 * t706;
t697 = m(6) * t712 - mrSges(6,2) * t771 + mrSges(6,3) * t737 + t750 * t780 - t761 * t792 + t855;
t711 = t716 * t837 - t717 * t832;
t760 = -mrSges(6,2) * t792 + mrSges(6,3) * t780;
t709 = -pkin(5) * t771 - pkin(12) * t791 + t751 * t781 - t711;
t845 = -m(7) * t709 + t725 * mrSges(7,1) - mrSges(7,2) * t726 + t758 * t739 - t740 * t759;
t701 = m(6) * t711 + mrSges(6,1) * t771 - mrSges(6,3) * t738 - t750 * t781 + t760 * t792 + t845;
t856 = t837 * t697 - t832 * t701;
t848 = m(5) * t721 - t773 * mrSges(5,3) - t795 * t784 + t856;
t783 = mrSges(5,1) * t794 + mrSges(5,3) * t809;
t872 = mrSges(4,2) * t809 - mrSges(4,3) * t794 - t783;
t886 = mrSges(4,1) - mrSges(5,2);
t686 = m(4) * t733 + t773 * mrSges(4,2) + t886 * t772 + t795 * t785 + t872 * t794 + t848;
t777 = mrSges(4,1) * t794 + mrSges(4,2) * t795;
t691 = t832 * t697 + t837 * t701;
t778 = -mrSges(5,2) * t794 - mrSges(5,3) * t795;
t844 = -m(5) * t723 - t773 * mrSges(5,1) - t795 * t778 - t691;
t689 = m(4) * t727 - t773 * mrSges(4,3) - t795 * t777 + t886 * t798 - t872 * t809 + t844;
t698 = t836 * t705 + t831 * t706;
t843 = -m(6) * t719 + t737 * mrSges(6,1) - t738 * mrSges(6,2) + t780 * t760 - t781 * t761 - t698;
t841 = -m(5) * t722 + t798 * mrSges(5,3) - t809 * t784 - t843;
t695 = t841 + (-t777 - t778) * t794 + (-mrSges(4,3) - mrSges(5,1)) * t772 + t809 * t785 - t798 * mrSges(4,2) + m(4) * t728;
t677 = -t829 * t686 + t689 * t853 + t695 * t860;
t877 = t830 * t838;
t788 = -g(3) * t877 + t854;
t815 = -t851 * mrSges(3,2) + mrSges(3,3) * t861;
t818 = (-mrSges(3,1) * t838 + mrSges(3,2) * t834) * t870;
t673 = m(3) * t788 + t827 * mrSges(3,1) - t819 * mrSges(3,3) + t851 * t815 - t818 * t862 + t677;
t682 = -t833 * t689 + t890 * t695;
t878 = t830 * t834;
t789 = -g(3) * t878 + t871;
t814 = t851 * mrSges(3,1) - mrSges(3,3) * t862;
t681 = m(3) * t789 - t827 * mrSges(3,2) + t820 * mrSges(3,3) - t851 * t814 + t818 * t861 + t682;
t669 = -t673 * t834 + t838 * t681;
t889 = pkin(9) * t669;
t676 = t881 * t686 + t689 * t865 + t695 * t879;
t802 = -t830 * t816 - t864;
t675 = m(3) * t802 - t820 * mrSges(3,1) + t819 * mrSges(3,2) + (t814 * t834 - t815 * t838) * t870 + t676;
t663 = t673 * t858 - t675 * t830 + t681 * t859;
t661 = m(2) * t825 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t840 + t663;
t668 = m(2) * t826 - mrSges(2,1) * t840 - qJDD(1) * mrSges(2,2) + t669;
t876 = t839 * t661 + t835 * t668;
t875 = t884 * t794 - t885 * t795 + t892 * t809;
t874 = t894 * t794 - t883 * t795 - t884 * t809;
t873 = -t883 * t794 - t893 * t795 + t885 * t809;
t662 = t673 * t877 + t882 * t675 + t681 * t878;
t857 = -t661 * t835 + t839 * t668;
t801 = Ifges(3,5) * qJD(2) + (Ifges(3,5) * t882 + (Ifges(3,1) * t834 + Ifges(3,4) * t838) * t830) * qJD(1);
t800 = Ifges(3,6) * qJD(2) + (Ifges(3,6) * t882 + (Ifges(3,4) * t834 + Ifges(3,2) * t838) * t830) * qJD(1);
t799 = Ifges(3,3) * qJD(2) + (Ifges(3,3) * t882 + (Ifges(3,5) * t834 + Ifges(3,6) * t838) * t830) * qJD(1);
t743 = Ifges(6,1) * t781 + Ifges(6,4) * t780 + Ifges(6,5) * t792;
t742 = Ifges(6,4) * t781 + Ifges(6,2) * t780 + Ifges(6,6) * t792;
t741 = Ifges(6,5) * t781 + Ifges(6,6) * t780 + Ifges(6,3) * t792;
t731 = Ifges(7,1) * t759 + Ifges(7,4) * t758 + Ifges(7,5) * t779;
t730 = Ifges(7,4) * t759 + Ifges(7,2) * t758 + Ifges(7,6) * t779;
t729 = Ifges(7,5) * t759 + Ifges(7,6) * t758 + Ifges(7,3) * t779;
t700 = mrSges(7,2) * t709 - mrSges(7,3) * t707 + Ifges(7,1) * t726 + Ifges(7,4) * t725 + Ifges(7,5) * t736 + t729 * t758 - t730 * t779;
t699 = -mrSges(7,1) * t709 + mrSges(7,3) * t708 + Ifges(7,4) * t726 + Ifges(7,2) * t725 + Ifges(7,6) * t736 - t729 * t759 + t731 * t779;
t690 = -t772 * mrSges(5,2) - t794 * t783 + t848;
t684 = -mrSges(6,1) * t719 - mrSges(7,1) * t707 + mrSges(7,2) * t708 + mrSges(6,3) * t712 + Ifges(6,4) * t738 - Ifges(7,5) * t726 + Ifges(6,2) * t737 + Ifges(6,6) * t771 - Ifges(7,6) * t725 - Ifges(7,3) * t736 - pkin(5) * t698 - t730 * t759 + t731 * t758 - t741 * t781 + t743 * t792;
t683 = mrSges(6,2) * t719 - mrSges(6,3) * t711 + Ifges(6,1) * t738 + Ifges(6,4) * t737 + Ifges(6,5) * t771 - pkin(12) * t698 - t699 * t831 + t700 * t836 + t741 * t780 - t742 * t792;
t670 = -qJ(4) * t690 + t836 * t699 + t831 * t700 + pkin(12) * t855 - t780 * t743 + t781 * t742 + Ifges(6,3) * t771 + pkin(5) * t845 + Ifges(6,5) * t738 + Ifges(6,6) * t737 + mrSges(4,2) * t733 - mrSges(4,3) * t727 - mrSges(5,3) * t721 + mrSges(5,1) * t723 - mrSges(6,2) * t712 + mrSges(6,1) * t711 + pkin(4) * t691 + t874 * t809 + t885 * t798 + t875 * t794 + t893 * t773 + t883 * t772;
t665 = -mrSges(4,1) * t733 - mrSges(5,1) * t722 + mrSges(5,2) * t721 + mrSges(4,3) * t728 - pkin(3) * t690 - pkin(4) * t843 - pkin(11) * t856 - t832 * t683 - t837 * t684 + t894 * t772 - t883 * t773 + t875 * t795 + t884 * t798 + t873 * t809;
t664 = mrSges(4,1) * t727 - mrSges(4,2) * t728 + mrSges(5,2) * t723 - mrSges(5,3) * t722 + t837 * t683 - t832 * t684 - pkin(11) * t691 + pkin(3) * (t809 * t783 + t844) + qJ(4) * t841 + (-mrSges(5,2) * pkin(3) + t892) * t798 + t874 * t795 + (-qJ(4) * t778 - t873) * t794 + t885 * t773 + (-mrSges(5,1) * qJ(4) - t884) * t772;
t659 = Ifges(3,1) * t819 + Ifges(3,4) * t820 + Ifges(3,5) * t827 + t799 * t861 - t851 * t800 + mrSges(3,2) * t802 - mrSges(3,3) * t788 + t890 * t670 - t833 * t665 + (-t829 * t676 - t881 * t677) * pkin(10);
t658 = -mrSges(3,1) * t802 + mrSges(3,3) * t789 + Ifges(3,4) * t819 + Ifges(3,2) * t820 + Ifges(3,6) * t827 - pkin(2) * t676 - t829 * t664 + t665 * t853 + t670 * t860 + t682 * t863 - t799 * t862 + t851 * t801;
t657 = mrSges(3,1) * t788 - mrSges(3,2) * t789 + t881 * t664 + Ifges(3,5) * t819 + Ifges(3,6) * t820 + Ifges(3,3) * t827 + pkin(2) * t677 + (t800 * t834 - t801 * t838) * t870 + (pkin(10) * t682 + t890 * t665 + t670 * t833) * t829;
t656 = -mrSges(2,2) * g(3) - mrSges(2,3) * t825 + Ifges(2,5) * qJDD(1) - t840 * Ifges(2,6) - t834 * t658 + t838 * t659 + (-t662 * t830 - t882 * t663) * pkin(9);
t655 = mrSges(2,1) * g(3) + mrSges(2,3) * t826 + t840 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t662 - t830 * t657 + t658 * t858 + t659 * t859 + t882 * t889;
t1 = [-m(1) * g(1) + t857; -m(1) * g(2) + t876; (-m(1) - m(2)) * g(3) + t662; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t876 - t835 * t655 + t839 * t656; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t857 + t839 * t655 + t835 * t656; -mrSges(1,1) * g(2) + mrSges(2,1) * t825 + mrSges(1,2) * g(1) - mrSges(2,2) * t826 + t882 * t657 + Ifges(2,3) * qJDD(1) + pkin(1) * t663 + (t658 * t838 + t659 * t834 + t889) * t830;];
tauB  = t1;
