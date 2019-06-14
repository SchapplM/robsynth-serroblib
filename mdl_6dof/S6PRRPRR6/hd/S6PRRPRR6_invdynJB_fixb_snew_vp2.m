% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:42:03
% EndTime: 2019-05-05 05:42:50
% DurationCPUTime: 45.68s
% Computational Cost: add. (797474->352), mult. (1730667->468), div. (0->0), fcn. (1402480->16), ass. (0->157)
t899 = sin(pkin(12));
t903 = cos(pkin(12));
t888 = g(1) * t899 - g(2) * t903;
t889 = -g(1) * t903 - g(2) * t899;
t897 = -g(3) + qJDD(1);
t909 = sin(qJ(2));
t905 = cos(pkin(6));
t913 = cos(qJ(2));
t931 = t905 * t913;
t901 = sin(pkin(6));
t934 = t901 * t913;
t856 = t888 * t931 - t889 * t909 + t897 * t934;
t900 = sin(pkin(7));
t914 = qJD(2) ^ 2;
t854 = pkin(9) * t900 * t914 + qJDD(2) * pkin(2) + t856;
t932 = t905 * t909;
t935 = t901 * t909;
t857 = t888 * t932 + t913 * t889 + t897 * t935;
t927 = qJDD(2) * t900;
t855 = -pkin(2) * t914 + pkin(9) * t927 + t857;
t875 = -t888 * t901 + t897 * t905;
t904 = cos(pkin(7));
t908 = sin(qJ(3));
t912 = cos(qJ(3));
t823 = -t908 * t855 + (t854 * t904 + t875 * t900) * t912;
t929 = qJD(2) * t900;
t879 = (-pkin(3) * t912 - qJ(4) * t908) * t929;
t895 = qJD(2) * t904 + qJD(3);
t893 = t895 ^ 2;
t894 = qJDD(2) * t904 + qJDD(3);
t926 = t908 * t929;
t817 = -t894 * pkin(3) - t893 * qJ(4) + t879 * t926 + qJDD(4) - t823;
t898 = sin(pkin(13));
t902 = cos(pkin(13));
t873 = t895 * t902 - t898 * t926;
t928 = qJD(2) * t912;
t925 = t900 * t928;
t859 = mrSges(5,2) * t925 + mrSges(5,3) * t873;
t874 = t895 * t898 + t902 * t926;
t860 = -mrSges(5,1) * t925 - mrSges(5,3) * t874;
t881 = (qJD(3) * t928 + qJDD(2) * t908) * t900;
t861 = -t881 * t898 + t894 * t902;
t862 = t881 * t902 + t894 * t898;
t933 = t904 * t908;
t936 = t900 * t908;
t824 = t854 * t933 + t912 * t855 + t875 * t936;
t818 = -pkin(3) * t893 + qJ(4) * t894 + t879 * t925 + t824;
t871 = t904 * t875;
t882 = -qJD(3) * t926 + t912 * t927;
t821 = -t882 * pkin(3) - t881 * qJ(4) + t871 + (-t854 + (pkin(3) * t908 - qJ(4) * t912) * t895 * qJD(2)) * t900;
t806 = -0.2e1 * qJD(4) * t874 - t898 * t818 + t902 * t821;
t803 = (-t873 * t925 - t862) * pkin(10) + (t873 * t874 - t882) * pkin(4) + t806;
t807 = 0.2e1 * qJD(4) * t873 + t902 * t818 + t898 * t821;
t863 = -pkin(4) * t925 - pkin(10) * t874;
t872 = t873 ^ 2;
t805 = -pkin(4) * t872 + pkin(10) * t861 + t863 * t925 + t807;
t907 = sin(qJ(5));
t911 = cos(qJ(5));
t800 = t907 * t803 + t911 * t805;
t849 = t873 * t911 - t874 * t907;
t850 = t873 * t907 + t874 * t911;
t836 = -pkin(5) * t849 - pkin(11) * t850;
t876 = qJDD(5) - t882;
t887 = qJD(5) - t925;
t886 = t887 ^ 2;
t798 = -pkin(5) * t886 + pkin(11) * t876 + t836 * t849 + t800;
t808 = -t861 * pkin(4) - t872 * pkin(10) + t874 * t863 + t817;
t827 = -qJD(5) * t850 + t861 * t911 - t862 * t907;
t828 = qJD(5) * t849 + t861 * t907 + t862 * t911;
t801 = (-t849 * t887 - t828) * pkin(11) + (t850 * t887 - t827) * pkin(5) + t808;
t906 = sin(qJ(6));
t910 = cos(qJ(6));
t795 = -t798 * t906 + t801 * t910;
t838 = -t850 * t906 + t887 * t910;
t811 = qJD(6) * t838 + t828 * t910 + t876 * t906;
t839 = t850 * t910 + t887 * t906;
t822 = -mrSges(7,1) * t838 + mrSges(7,2) * t839;
t826 = qJDD(6) - t827;
t848 = qJD(6) - t849;
t829 = -mrSges(7,2) * t848 + mrSges(7,3) * t838;
t792 = m(7) * t795 + mrSges(7,1) * t826 - mrSges(7,3) * t811 - t822 * t839 + t829 * t848;
t796 = t798 * t910 + t801 * t906;
t810 = -qJD(6) * t839 - t828 * t906 + t876 * t910;
t830 = mrSges(7,1) * t848 - mrSges(7,3) * t839;
t793 = m(7) * t796 - mrSges(7,2) * t826 + mrSges(7,3) * t810 + t822 * t838 - t830 * t848;
t783 = t910 * t792 + t906 * t793;
t840 = -mrSges(6,2) * t887 + mrSges(6,3) * t849;
t841 = mrSges(6,1) * t887 - mrSges(6,3) * t850;
t917 = m(6) * t808 - t827 * mrSges(6,1) + mrSges(6,2) * t828 - t849 * t840 + t841 * t850 + t783;
t782 = m(5) * t817 - t861 * mrSges(5,1) + mrSges(5,2) * t862 - t873 * t859 + t860 * t874 + t917;
t878 = -mrSges(4,2) * t895 + mrSges(4,3) * t925;
t880 = (-mrSges(4,1) * t912 + mrSges(4,2) * t908) * t929;
t778 = m(4) * t823 + mrSges(4,1) * t894 - mrSges(4,3) * t881 + t878 * t895 - t880 * t926 - t782;
t937 = t778 * t912;
t877 = mrSges(4,1) * t895 - mrSges(4,3) * t926;
t784 = -t792 * t906 + t910 * t793;
t835 = -mrSges(6,1) * t849 + mrSges(6,2) * t850;
t781 = m(6) * t800 - mrSges(6,2) * t876 + mrSges(6,3) * t827 + t835 * t849 - t841 * t887 + t784;
t799 = t803 * t911 - t805 * t907;
t797 = -pkin(5) * t876 - pkin(11) * t886 + t836 * t850 - t799;
t794 = -m(7) * t797 + t810 * mrSges(7,1) - mrSges(7,2) * t811 + t838 * t829 - t830 * t839;
t788 = m(6) * t799 + mrSges(6,1) * t876 - mrSges(6,3) * t828 - t835 * t850 + t840 * t887 + t794;
t775 = t907 * t781 + t911 * t788;
t853 = -mrSges(5,1) * t873 + mrSges(5,2) * t874;
t773 = m(5) * t806 - mrSges(5,1) * t882 - mrSges(5,3) * t862 - t853 * t874 - t859 * t925 + t775;
t922 = t911 * t781 - t788 * t907;
t774 = m(5) * t807 + mrSges(5,2) * t882 + mrSges(5,3) * t861 + t853 * t873 + t860 * t925 + t922;
t923 = -t773 * t898 + t902 * t774;
t764 = m(4) * t824 - mrSges(4,2) * t894 + mrSges(4,3) * t882 - t877 * t895 + t880 * t925 + t923;
t767 = t902 * t773 + t898 * t774;
t837 = -t900 * t854 + t871;
t766 = m(4) * t837 - t882 * mrSges(4,1) + t881 * mrSges(4,2) + (t877 * t908 - t878 * t912) * t929 + t767;
t753 = t764 * t933 - t766 * t900 + t904 * t937;
t749 = m(3) * t856 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t914 + t753;
t752 = t764 * t936 + t904 * t766 + t900 * t937;
t751 = m(3) * t875 + t752;
t760 = t912 * t764 - t778 * t908;
t759 = m(3) * t857 - mrSges(3,1) * t914 - qJDD(2) * mrSges(3,2) + t760;
t739 = t749 * t931 - t751 * t901 + t759 * t932;
t737 = m(2) * t888 + t739;
t745 = -t749 * t909 + t913 * t759;
t744 = m(2) * t889 + t745;
t930 = t903 * t737 + t899 * t744;
t738 = t749 * t934 + t905 * t751 + t759 * t935;
t924 = -t737 * t899 + t903 * t744;
t921 = m(2) * t897 + t738;
t812 = Ifges(7,5) * t839 + Ifges(7,6) * t838 + Ifges(7,3) * t848;
t814 = Ifges(7,1) * t839 + Ifges(7,4) * t838 + Ifges(7,5) * t848;
t785 = -mrSges(7,1) * t797 + mrSges(7,3) * t796 + Ifges(7,4) * t811 + Ifges(7,2) * t810 + Ifges(7,6) * t826 - t812 * t839 + t814 * t848;
t813 = Ifges(7,4) * t839 + Ifges(7,2) * t838 + Ifges(7,6) * t848;
t786 = mrSges(7,2) * t797 - mrSges(7,3) * t795 + Ifges(7,1) * t811 + Ifges(7,4) * t810 + Ifges(7,5) * t826 + t812 * t838 - t813 * t848;
t831 = Ifges(6,5) * t850 + Ifges(6,6) * t849 + Ifges(6,3) * t887;
t832 = Ifges(6,4) * t850 + Ifges(6,2) * t849 + Ifges(6,6) * t887;
t768 = mrSges(6,2) * t808 - mrSges(6,3) * t799 + Ifges(6,1) * t828 + Ifges(6,4) * t827 + Ifges(6,5) * t876 - pkin(11) * t783 - t785 * t906 + t786 * t910 + t831 * t849 - t832 * t887;
t833 = Ifges(6,1) * t850 + Ifges(6,4) * t849 + Ifges(6,5) * t887;
t916 = mrSges(7,1) * t795 - mrSges(7,2) * t796 + Ifges(7,5) * t811 + Ifges(7,6) * t810 + Ifges(7,3) * t826 + t813 * t839 - t814 * t838;
t769 = -mrSges(6,1) * t808 + mrSges(6,3) * t800 + Ifges(6,4) * t828 + Ifges(6,2) * t827 + Ifges(6,6) * t876 - pkin(5) * t783 - t831 * t850 + t833 * t887 - t916;
t842 = Ifges(5,5) * t874 + Ifges(5,6) * t873 - Ifges(5,3) * t925;
t844 = Ifges(5,1) * t874 + Ifges(5,4) * t873 - Ifges(5,5) * t925;
t754 = -mrSges(5,1) * t817 + mrSges(5,3) * t807 + Ifges(5,4) * t862 + Ifges(5,2) * t861 - Ifges(5,6) * t882 - pkin(4) * t917 + pkin(10) * t922 + t907 * t768 + t911 * t769 - t874 * t842 - t844 * t925;
t843 = Ifges(5,4) * t874 + Ifges(5,2) * t873 - Ifges(5,6) * t925;
t755 = mrSges(5,2) * t817 - mrSges(5,3) * t806 + Ifges(5,1) * t862 + Ifges(5,4) * t861 - Ifges(5,5) * t882 - pkin(10) * t775 + t768 * t911 - t769 * t907 + t842 * t873 + t843 * t925;
t866 = Ifges(4,6) * t895 + (Ifges(4,4) * t908 + Ifges(4,2) * t912) * t929;
t867 = Ifges(4,5) * t895 + (Ifges(4,1) * t908 + Ifges(4,4) * t912) * t929;
t740 = Ifges(4,5) * t881 + Ifges(4,6) * t882 + Ifges(4,3) * t894 + mrSges(4,1) * t823 - mrSges(4,2) * t824 + t898 * t755 + t902 * t754 - pkin(3) * t782 + qJ(4) * t923 + (t866 * t908 - t867 * t912) * t929;
t865 = Ifges(4,3) * t895 + (Ifges(4,5) * t908 + Ifges(4,6) * t912) * t929;
t741 = mrSges(4,2) * t837 - mrSges(4,3) * t823 + Ifges(4,1) * t881 + Ifges(4,4) * t882 + Ifges(4,5) * t894 - qJ(4) * t767 - t754 * t898 + t755 * t902 + t865 * t925 - t866 * t895;
t915 = mrSges(6,1) * t799 - mrSges(6,2) * t800 + Ifges(6,5) * t828 + Ifges(6,6) * t827 + Ifges(6,3) * t876 + pkin(5) * t794 + pkin(11) * t784 + t910 * t785 + t906 * t786 + t850 * t832 - t849 * t833;
t746 = -mrSges(4,1) * t837 - mrSges(5,1) * t806 + mrSges(5,2) * t807 - t865 * t926 + Ifges(4,6) * t894 + t895 * t867 + Ifges(4,4) * t881 + mrSges(4,3) * t824 + (Ifges(4,2) + Ifges(5,3)) * t882 + t873 * t844 - t874 * t843 - Ifges(5,6) * t861 - Ifges(5,5) * t862 - pkin(4) * t775 - pkin(3) * t767 - t915;
t918 = pkin(9) * t760 + t741 * t908 + t746 * t912;
t734 = -mrSges(3,1) * t875 + mrSges(3,3) * t857 + t914 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t752 - t900 * t740 + t904 * t918;
t735 = mrSges(3,2) * t875 - mrSges(3,3) * t856 + Ifges(3,5) * qJDD(2) - t914 * Ifges(3,6) + t912 * t741 - t908 * t746 + (-t752 * t900 - t753 * t904) * pkin(9);
t919 = pkin(8) * t745 + t734 * t913 + t735 * t909;
t733 = mrSges(3,1) * t856 - mrSges(3,2) * t857 + Ifges(3,3) * qJDD(2) + pkin(2) * t753 + t904 * t740 + t900 * t918;
t732 = mrSges(2,2) * t897 - mrSges(2,3) * t888 - t909 * t734 + t913 * t735 + (-t738 * t901 - t739 * t905) * pkin(8);
t731 = -mrSges(2,1) * t897 + mrSges(2,3) * t889 - pkin(1) * t738 - t901 * t733 + t905 * t919;
t1 = [-m(1) * g(1) + t924; -m(1) * g(2) + t930; -m(1) * g(3) + t921; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t930 - t899 * t731 + t903 * t732; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t924 + t903 * t731 + t899 * t732; -mrSges(1,1) * g(2) + mrSges(2,1) * t888 + mrSges(1,2) * g(1) - mrSges(2,2) * t889 + pkin(1) * t739 + t905 * t733 + t901 * t919; t921; t733; t740; t782; t915; t916;];
tauJB  = t1;
