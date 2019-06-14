% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 21:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR13_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:54:55
% EndTime: 2019-05-05 20:55:35
% DurationCPUTime: 38.83s
% Computational Cost: add. (547958->386), mult. (1736231->510), div. (0->0), fcn. (1459104->14), ass. (0->178)
t925 = Ifges(4,1) + Ifges(5,2);
t917 = Ifges(4,5) - Ifges(5,4);
t924 = -Ifges(4,2) - Ifges(5,3);
t916 = Ifges(4,6) - Ifges(5,5);
t915 = -Ifges(5,6) - Ifges(4,4);
t923 = Ifges(4,3) + Ifges(5,1);
t853 = cos(pkin(12));
t851 = sin(pkin(7));
t914 = cos(pkin(6));
t889 = t914 * t851;
t852 = sin(pkin(6));
t913 = cos(pkin(7));
t892 = t852 * t913;
t831 = (t853 * t892 + t889) * qJD(1) * pkin(9);
t857 = sin(qJ(1));
t860 = cos(qJ(1));
t847 = -g(1) * t860 - g(2) * t857;
t861 = qJD(1) ^ 2;
t911 = qJ(2) * t852;
t836 = -pkin(1) * t861 + qJDD(1) * t911 + t847;
t850 = sin(pkin(12));
t919 = pkin(9) * t851;
t876 = -pkin(2) * t853 - t850 * t919;
t895 = t913 * pkin(9);
t901 = qJD(1) * t852;
t870 = qJD(1) * t876 * t901 + qJDD(1) * t895;
t846 = t857 * g(1) - g(2) * t860;
t835 = qJDD(1) * pkin(1) + t861 * t911 + t846;
t891 = t853 * t914;
t894 = qJD(2) * t901;
t907 = t852 * t853;
t878 = -g(3) * t907 + t835 * t891 - 0.2e1 * t850 * t894;
t884 = qJDD(1) * t914;
t888 = qJD(1) * t914;
t770 = pkin(2) * t884 + t831 * t888 + (-t870 * t852 - t836) * t850 + t878;
t837 = (-pkin(9) * t850 * t892 + t914 * pkin(2)) * qJD(1);
t893 = t850 * t914;
t897 = t835 * t893 + (t836 + 0.2e1 * t894) * t853;
t771 = t884 * t919 - t837 * t888 + (-g(3) * t850 + t870 * t853) * t852 + t897;
t879 = -t914 * g(3) + qJDD(2);
t781 = (-t835 + t876 * qJDD(1) + (-t831 * t853 + t837 * t850) * qJD(1)) * t852 + t879;
t856 = sin(qJ(3));
t890 = t856 * t913;
t908 = t851 * t856;
t920 = cos(qJ(3));
t747 = t770 * t890 + t920 * t771 + t781 * t908;
t883 = t913 * t920;
t909 = t850 * t852;
t921 = t856 * t909 - t883 * t907 - t920 * t889;
t816 = t921 * qJD(1);
t862 = t856 * t889 + (t920 * t850 + t853 * t890) * t852;
t817 = t862 * qJD(1);
t796 = pkin(3) * t816 - qJ(4) * t817;
t922 = -t851 * t907 + t914 * t913;
t832 = -t922 * qJD(1) - qJD(3);
t828 = t832 ^ 2;
t829 = t922 * qJDD(1) + qJDD(3);
t742 = pkin(3) * t828 - t829 * qJ(4) + 0.2e1 * qJD(4) * t832 + t816 * t796 - t747;
t918 = mrSges(4,1) - mrSges(5,2);
t896 = t851 * t920;
t746 = t770 * t883 - t856 * t771 + t781 * t896;
t797 = mrSges(4,1) * t816 + mrSges(4,2) * t817;
t800 = -t816 * qJD(3) + t862 * qJDD(1);
t744 = -t829 * pkin(3) - t828 * qJ(4) + t817 * t796 + qJDD(4) - t746;
t910 = t816 * t832;
t737 = (t816 * t817 - t829) * pkin(10) + (t800 - t910) * pkin(4) + t744;
t799 = qJD(3) * t817 + t921 * qJDD(1);
t810 = pkin(4) * t817 + pkin(10) * t832;
t815 = t816 ^ 2;
t755 = -t770 * t851 + t913 * t781;
t864 = (-t800 - t910) * qJ(4) + t755 + (-t832 * pkin(3) - 0.2e1 * qJD(4)) * t817;
t741 = -pkin(4) * t815 - t810 * t817 + (pkin(3) + pkin(10)) * t799 + t864;
t855 = sin(qJ(5));
t859 = cos(qJ(5));
t734 = t855 * t737 + t859 * t741;
t805 = t816 * t855 - t832 * t859;
t765 = -qJD(5) * t805 + t799 * t859 - t829 * t855;
t804 = t816 * t859 + t832 * t855;
t772 = -mrSges(6,1) * t804 + mrSges(6,2) * t805;
t814 = qJD(5) + t817;
t783 = mrSges(6,1) * t814 - mrSges(6,3) * t805;
t795 = qJDD(5) + t800;
t773 = -pkin(5) * t804 - pkin(11) * t805;
t813 = t814 ^ 2;
t732 = -pkin(5) * t813 + pkin(11) * t795 + t773 * t804 + t734;
t739 = -pkin(4) * t799 - pkin(10) * t815 - t832 * t810 - t742;
t766 = qJD(5) * t804 + t799 * t855 + t829 * t859;
t735 = (-t804 * t814 - t766) * pkin(11) + (t805 * t814 - t765) * pkin(5) + t739;
t854 = sin(qJ(6));
t858 = cos(qJ(6));
t729 = -t732 * t854 + t735 * t858;
t779 = -t805 * t854 + t814 * t858;
t750 = qJD(6) * t779 + t766 * t858 + t795 * t854;
t780 = t805 * t858 + t814 * t854;
t756 = -mrSges(7,1) * t779 + mrSges(7,2) * t780;
t801 = qJD(6) - t804;
t757 = -mrSges(7,2) * t801 + mrSges(7,3) * t779;
t763 = qJDD(6) - t765;
t727 = m(7) * t729 + mrSges(7,1) * t763 - mrSges(7,3) * t750 - t756 * t780 + t757 * t801;
t730 = t732 * t858 + t735 * t854;
t749 = -qJD(6) * t780 - t766 * t854 + t795 * t858;
t758 = mrSges(7,1) * t801 - mrSges(7,3) * t780;
t728 = m(7) * t730 - mrSges(7,2) * t763 + mrSges(7,3) * t749 + t756 * t779 - t758 * t801;
t885 = -t727 * t854 + t858 * t728;
t719 = m(6) * t734 - mrSges(6,2) * t795 + mrSges(6,3) * t765 + t772 * t804 - t783 * t814 + t885;
t733 = t737 * t859 - t741 * t855;
t782 = -mrSges(6,2) * t814 + mrSges(6,3) * t804;
t731 = -pkin(5) * t795 - pkin(11) * t813 + t773 * t805 - t733;
t869 = -m(7) * t731 + t749 * mrSges(7,1) - mrSges(7,2) * t750 + t779 * t757 - t758 * t780;
t723 = m(6) * t733 + mrSges(6,1) * t795 - mrSges(6,3) * t766 - t772 * t805 + t782 * t814 + t869;
t713 = t719 * t855 + t723 * t859;
t798 = -mrSges(5,2) * t816 - mrSges(5,3) * t817;
t868 = -m(5) * t744 - t800 * mrSges(5,1) - t817 * t798 - t713;
t807 = mrSges(5,1) * t816 + mrSges(5,3) * t832;
t902 = mrSges(4,2) * t832 - mrSges(4,3) * t816 - t807;
t709 = m(4) * t746 - mrSges(4,3) * t800 - t797 * t817 + t918 * t829 - t902 * t832 + t868;
t809 = -mrSges(4,1) * t832 - mrSges(4,3) * t817;
t745 = pkin(3) * t799 + t864;
t808 = mrSges(5,1) * t817 - mrSges(5,2) * t832;
t886 = t859 * t719 - t723 * t855;
t874 = m(5) * t745 - t800 * mrSges(5,3) - t817 * t808 + t886;
t711 = m(4) * t755 + mrSges(4,2) * t800 + t918 * t799 + t809 * t817 + t902 * t816 + t874;
t720 = t858 * t727 + t854 * t728;
t867 = -m(6) * t739 + mrSges(6,1) * t765 - t766 * mrSges(6,2) + t782 * t804 - t805 * t783 - t720;
t863 = -m(5) * t742 + t829 * mrSges(5,3) - t832 * t808 - t867;
t717 = t863 + (-t797 - t798) * t816 + (-mrSges(4,3) - mrSges(5,1)) * t799 + t809 * t832 - mrSges(4,2) * t829 + m(4) * t747;
t699 = t709 * t883 - t851 * t711 + t717 * t890;
t802 = -t850 * t836 + t878;
t881 = -mrSges(3,1) * t853 + mrSges(3,2) * t850;
t834 = t881 * t901;
t872 = -mrSges(3,2) * t914 + mrSges(3,3) * t907;
t839 = t872 * qJD(1);
t873 = mrSges(3,1) * t914 - mrSges(3,3) * t909;
t695 = m(3) * t802 + t873 * qJDD(1) + (-t834 * t909 + t914 * t839) * qJD(1) + t699;
t704 = -t856 * t709 + t920 * t717;
t803 = -g(3) * t909 + t897;
t838 = t873 * qJD(1);
t703 = m(3) * t803 + t872 * qJDD(1) + (t834 * t907 - t914 * t838) * qJD(1) + t704;
t691 = -t695 * t850 + t853 * t703;
t912 = qJ(2) * t691;
t698 = t709 * t896 + t913 * t711 + t717 * t908;
t818 = -t852 * t835 + t879;
t697 = m(3) * t818 + (t881 * qJDD(1) + (t838 * t850 - t839 * t853) * qJD(1)) * t852 + t698;
t685 = t695 * t891 - t697 * t852 + t703 * t893;
t683 = m(2) * t846 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t861 + t685;
t690 = m(2) * t847 - mrSges(2,1) * t861 - qJDD(1) * mrSges(2,2) + t691;
t906 = t860 * t683 + t857 * t690;
t905 = t916 * t816 - t917 * t817 + t923 * t832;
t904 = t924 * t816 - t915 * t817 - t916 * t832;
t903 = -t915 * t816 - t925 * t817 + t917 * t832;
t684 = t695 * t907 + t914 * t697 + t703 * t909;
t887 = -t683 * t857 + t860 * t690;
t880 = Ifges(3,5) * t850 + Ifges(3,6) * t853;
t866 = t914 * Ifges(3,5) + (Ifges(3,1) * t850 + Ifges(3,4) * t853) * t852;
t865 = t914 * Ifges(3,6) + (Ifges(3,4) * t850 + Ifges(3,2) * t853) * t852;
t823 = t866 * qJD(1);
t822 = t865 * qJD(1);
t821 = (t914 * Ifges(3,3) + t880 * t852) * qJD(1);
t761 = Ifges(6,1) * t805 + Ifges(6,4) * t804 + Ifges(6,5) * t814;
t760 = Ifges(6,4) * t805 + Ifges(6,2) * t804 + Ifges(6,6) * t814;
t759 = Ifges(6,5) * t805 + Ifges(6,6) * t804 + Ifges(6,3) * t814;
t753 = Ifges(7,1) * t780 + Ifges(7,4) * t779 + Ifges(7,5) * t801;
t752 = Ifges(7,4) * t780 + Ifges(7,2) * t779 + Ifges(7,6) * t801;
t751 = Ifges(7,5) * t780 + Ifges(7,6) * t779 + Ifges(7,3) * t801;
t722 = mrSges(7,2) * t731 - mrSges(7,3) * t729 + Ifges(7,1) * t750 + Ifges(7,4) * t749 + Ifges(7,5) * t763 + t751 * t779 - t752 * t801;
t721 = -mrSges(7,1) * t731 + mrSges(7,3) * t730 + Ifges(7,4) * t750 + Ifges(7,2) * t749 + Ifges(7,6) * t763 - t751 * t780 + t753 * t801;
t712 = -mrSges(5,2) * t799 - t807 * t816 + t874;
t706 = -mrSges(6,1) * t739 - mrSges(7,1) * t729 + mrSges(7,2) * t730 + mrSges(6,3) * t734 + Ifges(6,4) * t766 - Ifges(7,5) * t750 + Ifges(6,2) * t765 + Ifges(6,6) * t795 - Ifges(7,6) * t749 - Ifges(7,3) * t763 - pkin(5) * t720 - t752 * t780 + t753 * t779 - t759 * t805 + t761 * t814;
t705 = mrSges(6,2) * t739 - mrSges(6,3) * t733 + Ifges(6,1) * t766 + Ifges(6,4) * t765 + Ifges(6,5) * t795 - pkin(11) * t720 - t721 * t854 + t722 * t858 + t759 * t804 - t760 * t814;
t692 = pkin(4) * t713 - qJ(4) * t712 + t858 * t721 + pkin(11) * t885 + t854 * t722 + t805 * t760 - t804 * t761 + Ifges(6,3) * t795 + pkin(5) * t869 + Ifges(6,6) * t765 + Ifges(6,5) * t766 + mrSges(4,2) * t755 - mrSges(4,3) * t746 + mrSges(5,1) * t744 - mrSges(5,3) * t745 - mrSges(6,2) * t734 + mrSges(6,1) * t733 + t904 * t832 + t917 * t829 + t905 * t816 + t925 * t800 + t915 * t799;
t687 = -mrSges(4,1) * t755 - mrSges(5,1) * t742 + mrSges(5,2) * t745 + mrSges(4,3) * t747 - pkin(3) * t712 - pkin(4) * t867 - pkin(10) * t886 - t855 * t705 - t859 * t706 + t924 * t799 - t915 * t800 + t905 * t817 + t916 * t829 + t903 * t832;
t686 = mrSges(4,1) * t746 - mrSges(4,2) * t747 + mrSges(5,2) * t744 - mrSges(5,3) * t742 + t859 * t705 - t855 * t706 - pkin(10) * t713 + pkin(3) * (t807 * t832 + t868) + qJ(4) * t863 + (-pkin(3) * mrSges(5,2) + t923) * t829 + t904 * t817 + (-qJ(4) * t798 - t903) * t816 + t917 * t800 + (-qJ(4) * mrSges(5,1) - t916) * t799;
t681 = mrSges(3,2) * t818 - mrSges(3,3) * t802 + t920 * t692 - t856 * t687 + (t821 * t907 - t914 * t822) * qJD(1) + (-t698 * t851 - t913 * t699) * pkin(9) + t866 * qJDD(1);
t680 = t687 * t883 + t704 * t895 + t692 * t890 - mrSges(3,1) * t818 + mrSges(3,3) * t803 - pkin(2) * t698 - t851 * t686 + (-t821 * t909 + t914 * t823) * qJD(1) + t865 * qJDD(1);
t679 = Ifges(3,3) * t884 + mrSges(3,1) * t802 - mrSges(3,2) * t803 + t913 * t686 + pkin(2) * t699 + (pkin(9) * t704 + t920 * t687 + t692 * t856) * t851 + (t880 * qJDD(1) + (t822 * t850 - t823 * t853) * qJD(1)) * t852;
t678 = -mrSges(2,2) * g(3) - mrSges(2,3) * t846 + Ifges(2,5) * qJDD(1) - t861 * Ifges(2,6) - t850 * t680 + t853 * t681 + (-t684 * t852 - t914 * t685) * qJ(2);
t677 = mrSges(2,1) * g(3) + mrSges(2,3) * t847 + t861 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t684 - t852 * t679 + t680 * t891 + t681 * t893 + t914 * t912;
t1 = [-m(1) * g(1) + t887; -m(1) * g(2) + t906; (-m(1) - m(2)) * g(3) + t684; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t906 - t857 * t677 + t860 * t678; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t887 + t860 * t677 + t857 * t678; -mrSges(1,1) * g(2) + mrSges(2,1) * t846 + mrSges(1,2) * g(1) - mrSges(2,2) * t847 + t914 * t679 + Ifges(2,3) * qJDD(1) + pkin(1) * t685 + (t680 * t853 + t681 * t850 + t912) * t852;];
tauB  = t1;
