% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-06 00:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:43:00
% EndTime: 2019-05-06 00:43:45
% DurationCPUTime: 43.20s
% Computational Cost: add. (617194->386), mult. (1925370->508), div. (0->0), fcn. (1631441->14), ass. (0->179)
t897 = Ifges(5,1) + Ifges(6,2);
t889 = Ifges(5,4) + Ifges(6,6);
t888 = Ifges(5,5) - Ifges(6,4);
t896 = Ifges(5,2) + Ifges(6,3);
t887 = Ifges(5,6) - Ifges(6,5);
t895 = -Ifges(5,3) - Ifges(6,1);
t823 = sin(pkin(12));
t826 = cos(pkin(12));
t831 = sin(qJ(3));
t827 = cos(pkin(7));
t834 = cos(qJ(3));
t873 = t827 * t834;
t894 = -t823 * t831 + t826 * t873;
t825 = sin(pkin(6));
t828 = cos(pkin(6));
t874 = t827 * t831;
t824 = sin(pkin(7));
t879 = t824 * t831;
t840 = t828 * t879 + (t823 * t834 + t826 * t874) * t825;
t795 = t840 * qJD(1);
t878 = t824 * t834;
t865 = t828 * t878;
t781 = -t795 * qJD(3) + qJDD(1) * (t825 * t894 + t865);
t876 = t825 * t827;
t808 = (t824 * t828 + t826 * t876) * qJD(1) * pkin(9);
t832 = sin(qJ(1));
t835 = cos(qJ(1));
t821 = -t835 * g(1) - t832 * g(2);
t836 = qJD(1) ^ 2;
t884 = qJ(2) * t825;
t812 = -t836 * pkin(1) + qJDD(1) * t884 + t821;
t890 = pkin(9) * t823;
t854 = -pkin(2) * t826 - t824 * t890;
t867 = qJD(1) * t825;
t885 = pkin(9) * qJDD(1);
t849 = qJD(1) * t854 * t867 + t827 * t885;
t820 = t832 * g(1) - t835 * g(2);
t811 = qJDD(1) * pkin(1) + t836 * t884 + t820;
t862 = qJD(2) * t867;
t875 = t826 * t828;
t877 = t825 * t826;
t855 = -g(3) * t877 + t811 * t875 - 0.2e1 * t823 * t862;
t757 = (pkin(2) * qJDD(1) + qJD(1) * t808) * t828 + (-t825 * t849 - t812) * t823 + t855;
t813 = (pkin(2) * t828 - t876 * t890) * qJD(1);
t881 = t823 * t828;
t863 = t811 * t881 + (t812 + 0.2e1 * t862) * t826;
t758 = (-qJD(1) * t813 + t824 * t885) * t828 + (-g(3) * t823 + t826 * t849) * t825 + t863;
t861 = -t828 * g(3) + qJDD(2);
t767 = (-t811 + t854 * qJDD(1) + (-t808 * t826 + t813 * t823) * qJD(1)) * t825 + t861;
t724 = -t831 * t758 + (t757 * t827 + t767 * t824) * t834;
t892 = -2 * qJD(5);
t891 = cos(qJ(4));
t886 = Ifges(3,3) * t828;
t850 = -t824 * t877 + t827 * t828;
t809 = qJD(1) * t850 + qJD(3);
t830 = sin(qJ(4));
t788 = t830 * t795 - t809 * t891;
t794 = qJD(1) * t865 + t867 * t894;
t793 = qJD(4) - t794;
t883 = t788 * t793;
t882 = t823 * t825;
t725 = t757 * t874 + t834 * t758 + t767 * t879;
t779 = -t794 * mrSges(4,1) + t795 * mrSges(4,2);
t791 = t809 * mrSges(4,1) - t795 * mrSges(4,3);
t806 = qJDD(1) * t850 + qJDD(3);
t780 = -t794 * pkin(3) - t795 * pkin(10);
t805 = t809 ^ 2;
t721 = -t805 * pkin(3) + t806 * pkin(10) + t794 * t780 + t725;
t733 = -t824 * t757 + t827 * t767;
t782 = t794 * qJD(3) + qJDD(1) * t840;
t723 = (-t794 * t809 - t782) * pkin(10) + (t795 * t809 - t781) * pkin(3) + t733;
t716 = -t830 * t721 + t723 * t891;
t752 = -t788 * qJD(4) + t782 * t891 + t830 * t806;
t789 = t795 * t891 + t830 * t809;
t760 = t788 * mrSges(5,1) + t789 * mrSges(5,2);
t768 = t788 * mrSges(6,1) - t793 * mrSges(6,3);
t770 = -t793 * mrSges(5,2) - t788 * mrSges(5,3);
t778 = qJDD(4) - t781;
t759 = t788 * pkin(4) - t789 * qJ(5);
t792 = t793 ^ 2;
t714 = -t778 * pkin(4) - t792 * qJ(5) + t789 * t759 + qJDD(5) - t716;
t709 = (t788 * t789 - t778) * pkin(11) + (t752 + t883) * pkin(5) + t714;
t751 = t789 * qJD(4) + t830 * t782 - t806 * t891;
t772 = t789 * pkin(5) - t793 * pkin(11);
t787 = t788 ^ 2;
t720 = -t806 * pkin(3) - t805 * pkin(10) + t795 * t780 - t724;
t837 = (-t752 + t883) * qJ(5) + t720 + (pkin(4) * t793 + t892) * t789;
t712 = -t787 * pkin(5) - t789 * t772 + (pkin(4) + pkin(11)) * t751 + t837;
t829 = sin(qJ(6));
t833 = cos(qJ(6));
t707 = t833 * t709 - t829 * t712;
t765 = t833 * t788 - t829 * t793;
t728 = t765 * qJD(6) + t829 * t751 + t833 * t778;
t766 = t829 * t788 + t833 * t793;
t734 = -t765 * mrSges(7,1) + t766 * mrSges(7,2);
t784 = qJD(6) + t789;
t737 = -t784 * mrSges(7,2) + t765 * mrSges(7,3);
t748 = qJDD(6) + t752;
t705 = m(7) * t707 + t748 * mrSges(7,1) - t728 * mrSges(7,3) - t766 * t734 + t784 * t737;
t708 = t829 * t709 + t833 * t712;
t727 = -t766 * qJD(6) + t833 * t751 - t829 * t778;
t738 = t784 * mrSges(7,1) - t766 * mrSges(7,3);
t706 = m(7) * t708 - t748 * mrSges(7,2) + t727 * mrSges(7,3) + t765 * t734 - t784 * t738;
t697 = t833 * t705 + t829 * t706;
t761 = -t788 * mrSges(6,2) - t789 * mrSges(6,3);
t842 = -m(6) * t714 - t752 * mrSges(6,1) - t789 * t761 - t697;
t695 = m(5) * t716 - t752 * mrSges(5,3) - t789 * t760 + (-t768 + t770) * t793 + (mrSges(5,1) - mrSges(6,2)) * t778 + t842;
t717 = t891 * t721 + t830 * t723;
t771 = t793 * mrSges(5,1) - t789 * mrSges(5,3);
t841 = -t792 * pkin(4) + t778 * qJ(5) - t788 * t759 + t717;
t713 = t793 * t892 - t841;
t769 = t789 * mrSges(6,1) + t793 * mrSges(6,2);
t711 = -t751 * pkin(5) - t787 * pkin(11) + ((2 * qJD(5)) + t772) * t793 + t841;
t845 = -m(7) * t711 + t727 * mrSges(7,1) - t728 * mrSges(7,2) + t765 * t737 - t766 * t738;
t839 = -m(6) * t713 + t778 * mrSges(6,3) + t793 * t769 - t845;
t702 = m(5) * t717 - t778 * mrSges(5,2) - t793 * t771 + (-t760 - t761) * t788 + (-mrSges(5,3) - mrSges(6,1)) * t751 + t839;
t859 = -t830 * t695 + t891 * t702;
t687 = m(4) * t725 - t806 * mrSges(4,2) + t781 * mrSges(4,3) + t794 * t779 - t809 * t791 + t859;
t690 = t891 * t695 + t830 * t702;
t790 = -t809 * mrSges(4,2) + t794 * mrSges(4,3);
t689 = m(4) * t733 - t781 * mrSges(4,1) + t782 * mrSges(4,2) - t794 * t790 + t795 * t791 + t690;
t715 = t751 * pkin(4) + t837;
t871 = -t829 * t705 + t833 * t706;
t848 = -m(6) * t715 + t751 * mrSges(6,2) + t788 * t768 - t871;
t838 = -m(5) * t720 - t751 * mrSges(5,1) - t788 * t770 + (t769 - t771) * t789 + (-mrSges(5,2) + mrSges(6,3)) * t752 + t848;
t693 = m(4) * t724 + t806 * mrSges(4,1) - t782 * mrSges(4,3) - t795 * t779 + t809 * t790 + t838;
t676 = t687 * t874 - t824 * t689 + t693 * t873;
t785 = -t823 * t812 + t855;
t858 = -mrSges(3,1) * t826 + mrSges(3,2) * t823;
t810 = t858 * t867;
t852 = -mrSges(3,2) * t828 + mrSges(3,3) * t877;
t815 = t852 * qJD(1);
t853 = mrSges(3,1) * t828 - mrSges(3,3) * t882;
t672 = m(3) * t785 + t853 * qJDD(1) + (-t810 * t882 + t815 * t828) * qJD(1) + t676;
t675 = t687 * t879 + t827 * t689 + t693 * t878;
t796 = -t825 * t811 + t861;
t814 = t853 * qJD(1);
t674 = m(3) * t796 + (t858 * qJDD(1) + (t814 * t823 - t815 * t826) * qJD(1)) * t825 + t675;
t683 = t834 * t687 - t831 * t693;
t786 = -g(3) * t882 + t863;
t682 = m(3) * t786 + t852 * qJDD(1) + (t810 * t877 - t814 * t828) * qJD(1) + t683;
t662 = t672 * t875 - t825 * t674 + t682 * t881;
t660 = m(2) * t820 + qJDD(1) * mrSges(2,1) - t836 * mrSges(2,2) + t662;
t668 = -t823 * t672 + t826 * t682;
t667 = m(2) * t821 - t836 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t668;
t872 = t835 * t660 + t832 * t667;
t870 = t788 * t887 - t789 * t888 + t793 * t895;
t869 = t788 * t896 - t789 * t889 - t793 * t887;
t868 = -t889 * t788 + t789 * t897 + t888 * t793;
t661 = t672 * t877 + t828 * t674 + t682 * t882;
t860 = -t832 * t660 + t835 * t667;
t857 = Ifges(3,5) * t823 + Ifges(3,6) * t826;
t696 = -t752 * mrSges(6,3) - t789 * t769 - t848;
t729 = Ifges(7,5) * t766 + Ifges(7,6) * t765 + Ifges(7,3) * t784;
t731 = Ifges(7,1) * t766 + Ifges(7,4) * t765 + Ifges(7,5) * t784;
t698 = -mrSges(7,1) * t711 + mrSges(7,3) * t708 + Ifges(7,4) * t728 + Ifges(7,2) * t727 + Ifges(7,6) * t748 - t766 * t729 + t784 * t731;
t730 = Ifges(7,4) * t766 + Ifges(7,2) * t765 + Ifges(7,6) * t784;
t699 = mrSges(7,2) * t711 - mrSges(7,3) * t707 + Ifges(7,1) * t728 + Ifges(7,4) * t727 + Ifges(7,5) * t748 + t765 * t729 - t784 * t730;
t677 = -mrSges(5,1) * t720 - mrSges(6,1) * t713 + mrSges(6,2) * t715 + mrSges(5,3) * t717 - pkin(4) * t696 - pkin(5) * t845 - pkin(11) * t871 - t833 * t698 - t829 * t699 - t751 * t896 + t889 * t752 + t887 * t778 + t870 * t789 + t868 * t793;
t678 = mrSges(6,1) * t714 + mrSges(7,1) * t707 + mrSges(5,2) * t720 - mrSges(7,2) * t708 - mrSges(5,3) * t716 - mrSges(6,3) * t715 + Ifges(7,5) * t728 + Ifges(7,6) * t727 + Ifges(7,3) * t748 + pkin(5) * t697 - qJ(5) * t696 + t766 * t730 - t765 * t731 + t869 * t793 + t870 * t788 + t888 * t778 + t897 * t752 - t889 * t751;
t774 = Ifges(4,5) * t795 + Ifges(4,6) * t794 + Ifges(4,3) * t809;
t775 = Ifges(4,4) * t795 + Ifges(4,2) * t794 + Ifges(4,6) * t809;
t664 = mrSges(4,2) * t733 - mrSges(4,3) * t724 + Ifges(4,1) * t782 + Ifges(4,4) * t781 + Ifges(4,5) * t806 - pkin(10) * t690 - t830 * t677 + t678 * t891 + t794 * t774 - t809 * t775;
t776 = Ifges(4,1) * t795 + Ifges(4,4) * t794 + Ifges(4,5) * t809;
t669 = -pkin(3) * t690 - qJ(5) * t839 - pkin(4) * (-t793 * t768 + t842) - t833 * t699 + t829 * t698 + t809 * t776 + Ifges(4,6) * t806 - t795 * t774 + Ifges(4,2) * t781 + Ifges(4,4) * t782 - mrSges(4,1) * t733 + mrSges(4,3) * t725 + mrSges(5,2) * t717 - mrSges(5,1) * t716 + mrSges(6,3) * t713 - mrSges(6,2) * t714 + pkin(11) * t697 + t869 * t789 + (qJ(5) * t761 - t868) * t788 + (pkin(4) * mrSges(6,2) + t895) * t778 - t888 * t752 + (qJ(5) * mrSges(6,1) + t887) * t751;
t847 = pkin(9) * t683 + t664 * t831 + t669 * t834;
t663 = mrSges(4,1) * t724 - mrSges(4,2) * t725 + Ifges(4,5) * t782 + Ifges(4,6) * t781 + Ifges(4,3) * t806 + pkin(3) * t838 + pkin(10) * t859 + t677 * t891 + t830 * t678 + t795 * t775 - t794 * t776;
t799 = (t825 * t857 + t886) * qJD(1);
t844 = Ifges(3,5) * t828 + (Ifges(3,1) * t823 + Ifges(3,4) * t826) * t825;
t801 = t844 * qJD(1);
t843 = Ifges(3,6) * t828 + (Ifges(3,4) * t823 + Ifges(3,2) * t826) * t825;
t657 = -mrSges(3,1) * t796 + mrSges(3,3) * t786 - pkin(2) * t675 - t824 * t663 + (-t799 * t882 + t801 * t828) * qJD(1) + t847 * t827 + t843 * qJDD(1);
t800 = t843 * qJD(1);
t658 = mrSges(3,2) * t796 - mrSges(3,3) * t785 + t834 * t664 - t831 * t669 + (t799 * t877 - t800 * t828) * qJD(1) + (-t675 * t824 - t676 * t827) * pkin(9) + t844 * qJDD(1);
t846 = qJ(2) * t668 + t657 * t826 + t658 * t823;
t656 = qJDD(1) * t886 + mrSges(3,1) * t785 - mrSges(3,2) * t786 + pkin(2) * t676 + t827 * t663 + t847 * t824 + (t857 * qJDD(1) + (t800 * t823 - t801 * t826) * qJD(1)) * t825;
t655 = -mrSges(2,2) * g(3) - mrSges(2,3) * t820 + Ifges(2,5) * qJDD(1) - t836 * Ifges(2,6) - t823 * t657 + t826 * t658 + (-t661 * t825 - t662 * t828) * qJ(2);
t654 = mrSges(2,1) * g(3) + mrSges(2,3) * t821 + t836 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t661 - t825 * t656 + t828 * t846;
t1 = [-m(1) * g(1) + t860; -m(1) * g(2) + t872; (-m(1) - m(2)) * g(3) + t661; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t872 - t832 * t654 + t835 * t655; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t860 + t835 * t654 + t832 * t655; -mrSges(1,1) * g(2) + mrSges(2,1) * t820 + mrSges(1,2) * g(1) - mrSges(2,2) * t821 + Ifges(2,3) * qJDD(1) + pkin(1) * t662 + t828 * t656 + t825 * t846;];
tauB  = t1;
