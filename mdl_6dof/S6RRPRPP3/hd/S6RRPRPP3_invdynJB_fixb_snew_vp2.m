% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-05-06 12:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:29:40
% EndTime: 2019-05-06 12:29:52
% DurationCPUTime: 7.09s
% Computational Cost: add. (91625->341), mult. (196596->398), div. (0->0), fcn. (133238->8), ass. (0->134)
t906 = sin(qJ(1));
t908 = cos(qJ(1));
t894 = -g(1) * t908 - g(2) * t906;
t910 = qJD(1) ^ 2;
t876 = -pkin(1) * t910 + qJDD(1) * pkin(7) + t894;
t905 = sin(qJ(2));
t907 = cos(qJ(2));
t859 = -t907 * g(3) - t905 * t876;
t886 = (-pkin(2) * t907 - qJ(3) * t905) * qJD(1);
t909 = qJD(2) ^ 2;
t941 = qJD(1) * t905;
t838 = -qJDD(2) * pkin(2) - t909 * qJ(3) + t886 * t941 + qJDD(3) - t859;
t939 = qJD(1) * qJD(2);
t929 = t907 * t939;
t888 = qJDD(1) * t905 + t929;
t902 = sin(pkin(9));
t903 = cos(pkin(9));
t863 = qJDD(2) * t903 - t888 * t902;
t883 = qJD(2) * t902 + t903 * t941;
t940 = qJD(1) * t907;
t865 = -pkin(3) * t940 - pkin(8) * t883;
t882 = qJD(2) * t903 - t902 * t941;
t881 = t882 ^ 2;
t804 = -t863 * pkin(3) - t881 * pkin(8) + t883 * t865 + t838;
t904 = sin(qJ(4));
t949 = cos(qJ(4));
t856 = t904 * t882 + t883 * t949;
t864 = qJDD(2) * t902 + t888 * t903;
t812 = t856 * qJD(4) - t863 * t949 + t904 * t864;
t855 = -t882 * t949 + t904 * t883;
t813 = -t855 * qJD(4) + t904 * t863 + t864 * t949;
t896 = -qJD(4) + t940;
t840 = mrSges(5,2) * t896 - mrSges(5,3) * t855;
t841 = -mrSges(5,1) * t896 - mrSges(5,3) * t856;
t846 = mrSges(6,1) * t856 - mrSges(6,2) * t896;
t945 = t855 * t896;
t951 = -2 * qJD(5);
t912 = (-t813 - t945) * qJ(5) + t804 + (-t896 * pkin(4) + t951) * t856;
t794 = t812 * pkin(4) + t912;
t844 = mrSges(6,1) * t855 + mrSges(6,3) * t896;
t842 = pkin(5) * t856 + qJ(6) * t896;
t854 = t855 ^ 2;
t950 = 2 * qJD(6);
t791 = t912 - t854 * pkin(5) + t855 * t950 - t856 * t842 + (pkin(4) + qJ(6)) * t812;
t843 = mrSges(7,1) * t856 + mrSges(7,3) * t896;
t845 = -mrSges(7,1) * t855 - mrSges(7,2) * t896;
t923 = m(7) * t791 - t813 * mrSges(7,2) + t812 * mrSges(7,3) - t856 * t843 + t855 * t845;
t918 = m(6) * t794 - t812 * mrSges(6,2) - t855 * t844 + t923;
t959 = m(5) * t804 + t812 * mrSges(5,1) + t855 * t840 + t918 + (t841 - t846) * t856 + (mrSges(5,2) - mrSges(6,3)) * t813;
t862 = -mrSges(4,1) * t940 - mrSges(4,3) * t883;
t921 = mrSges(4,2) * t940 + mrSges(4,3) * t882;
t781 = m(4) * t838 - t863 * mrSges(4,1) + t864 * mrSges(4,2) + t883 * t862 - t882 * t921 + t959;
t958 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t938 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t937 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t957 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t936 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t956 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t893 = t906 * g(1) - t908 * g(2);
t875 = -qJDD(1) * pkin(1) - t910 * pkin(7) - t893;
t898 = t905 * t939;
t889 = qJDD(1) * t907 - t898;
t833 = (-t888 - t929) * qJ(3) + (-t889 + t898) * pkin(2) + t875;
t860 = -g(3) * t905 + t907 * t876;
t839 = -pkin(2) * t909 + qJDD(2) * qJ(3) + t886 * t940 + t860;
t802 = -0.2e1 * qJD(3) * t883 + t903 * t833 - t902 * t839;
t798 = (-t882 * t940 - t864) * pkin(8) + (t882 * t883 - t889) * pkin(3) + t802;
t803 = 0.2e1 * qJD(3) * t882 + t902 * t833 + t903 * t839;
t801 = -pkin(3) * t881 + pkin(8) * t863 + t865 * t940 + t803;
t795 = t798 * t949 - t904 * t801;
t828 = pkin(4) * t855 - qJ(5) * t856;
t885 = qJDD(4) - t889;
t895 = t896 ^ 2;
t793 = -t885 * pkin(4) - t895 * qJ(5) + t856 * t828 + qJDD(5) - t795;
t830 = -mrSges(6,2) * t855 - mrSges(6,3) * t856;
t954 = -m(6) * t793 - t813 * mrSges(6,1) - t856 * t830;
t784 = -t813 * mrSges(6,3) - t856 * t846 + t918;
t827 = -mrSges(7,2) * t856 + mrSges(7,3) * t855;
t796 = t904 * t798 + t949 * t801;
t917 = -t895 * pkin(4) + t885 * qJ(5) - t855 * t828 + t796;
t789 = -t812 * pkin(5) - t854 * qJ(6) + qJDD(6) + (t951 - t842) * t896 + t917;
t934 = m(7) * t789 + t885 * mrSges(7,2) - t896 * t843;
t786 = -t812 * mrSges(7,1) - t855 * t827 + t934;
t792 = 0.2e1 * qJD(5) * t896 - t917;
t931 = -t938 * t855 + t958 * t856 - t937 * t896;
t933 = t936 * t855 - t937 * t856 + t956 * t896;
t768 = -mrSges(5,1) * t804 - mrSges(6,1) * t792 + mrSges(7,1) * t789 + mrSges(6,2) * t794 + mrSges(5,3) * t796 - mrSges(7,3) * t791 - pkin(4) * t784 + pkin(5) * t786 - qJ(6) * t923 + t957 * t812 + t938 * t813 + t933 * t856 + t936 * t885 - t931 * t896;
t787 = t896 * t950 + (t855 * t856 - t885) * qJ(6) + (t813 - t945) * pkin(5) + t793;
t924 = -m(7) * t787 + t885 * mrSges(7,3) - t896 * t845;
t785 = t813 * mrSges(7,1) + t856 * t827 - t924;
t932 = t957 * t855 + t938 * t856 - t936 * t896;
t769 = mrSges(6,1) * t793 + mrSges(7,1) * t787 + mrSges(5,2) * t804 - mrSges(7,2) * t791 - mrSges(5,3) * t795 - mrSges(6,3) * t794 + pkin(5) * t785 - qJ(5) * t784 - t938 * t812 + t958 * t813 + t933 * t855 + t937 * t885 + t932 * t896;
t849 = Ifges(4,5) * t883 + Ifges(4,6) * t882 - Ifges(4,3) * t940;
t851 = Ifges(4,1) * t883 + Ifges(4,4) * t882 - Ifges(4,5) * t940;
t829 = mrSges(5,1) * t855 + mrSges(5,2) * t856;
t947 = -mrSges(7,1) - mrSges(5,3);
t775 = m(5) * t795 + (-t840 + t844) * t896 + (mrSges(5,1) - mrSges(6,2)) * t885 + (-t827 - t829) * t856 + t947 * t813 + t924 + t954;
t920 = -m(6) * t792 + t885 * mrSges(6,3) - t896 * t846 + t934;
t943 = -t827 - t830;
t778 = m(5) * t796 - t885 * mrSges(5,2) + t896 * t841 + (-t829 + t943) * t855 + (-mrSges(6,1) + t947) * t812 + t920;
t925 = -t775 * t904 + t949 * t778;
t753 = -mrSges(4,1) * t838 + mrSges(4,3) * t803 + Ifges(4,4) * t864 + Ifges(4,2) * t863 - Ifges(4,6) * t889 - pkin(3) * t959 + pkin(8) * t925 + t949 * t768 + t904 * t769 - t883 * t849 - t851 * t940;
t773 = t949 * t775 + t904 * t778;
t850 = Ifges(4,4) * t883 + Ifges(4,2) * t882 - Ifges(4,6) * t940;
t754 = mrSges(4,2) * t838 - mrSges(4,3) * t802 + Ifges(4,1) * t864 + Ifges(4,4) * t863 - Ifges(4,5) * t889 - pkin(8) * t773 - t904 * t768 + t769 * t949 + t882 * t849 + t850 * t940;
t857 = -mrSges(4,1) * t882 + mrSges(4,2) * t883;
t771 = m(4) * t802 - t889 * mrSges(4,1) - t864 * mrSges(4,3) - t883 * t857 - t921 * t940 + t773;
t772 = m(4) * t803 + mrSges(4,2) * t889 + mrSges(4,3) * t863 + t857 * t882 + t862 * t940 + t925;
t767 = -t771 * t902 + t903 * t772;
t873 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t905 + Ifges(3,2) * t907) * qJD(1);
t874 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t905 + Ifges(3,4) * t907) * qJD(1);
t953 = mrSges(3,1) * t859 - mrSges(3,2) * t860 + Ifges(3,5) * t888 + Ifges(3,6) * t889 + Ifges(3,3) * qJDD(2) - pkin(2) * t781 + qJ(3) * t767 + t903 * t753 + t902 * t754 + (t873 * t905 - t874 * t907) * qJD(1);
t782 = t885 * mrSges(6,2) - t896 * t844 + t785 - t954;
t952 = t936 * t812 - t937 * t813 - t931 * t855 - t932 * t856 - t956 * t885 - mrSges(5,1) * t795 + mrSges(5,2) * t796 - mrSges(6,2) * t793 - mrSges(7,2) * t789 + mrSges(6,3) * t792 + mrSges(7,3) * t787 + pkin(4) * t782 - qJ(5) * (t943 * t855 + (-mrSges(6,1) - mrSges(7,1)) * t812 + t920) + qJ(6) * t785;
t887 = (-mrSges(3,1) * t907 + mrSges(3,2) * t905) * qJD(1);
t891 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t941;
t765 = m(3) * t860 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t889 - qJD(2) * t891 + t887 * t940 + t767;
t892 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t940;
t780 = m(3) * t859 + qJDD(2) * mrSges(3,1) - t888 * mrSges(3,3) + qJD(2) * t892 - t887 * t941 - t781;
t926 = t907 * t765 - t780 * t905;
t757 = m(2) * t894 - mrSges(2,1) * t910 - qJDD(1) * mrSges(2,2) + t926;
t766 = t771 * t903 + t772 * t902;
t915 = -m(3) * t875 + t889 * mrSges(3,1) - mrSges(3,2) * t888 - t891 * t941 + t892 * t940 - t766;
t761 = m(2) * t893 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t910 + t915;
t944 = t906 * t757 + t908 * t761;
t759 = t905 * t765 + t907 * t780;
t927 = t908 * t757 - t761 * t906;
t872 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t905 + Ifges(3,6) * t907) * qJD(1);
t750 = mrSges(3,2) * t875 - mrSges(3,3) * t859 + Ifges(3,1) * t888 + Ifges(3,4) * t889 + Ifges(3,5) * qJDD(2) - qJ(3) * t766 - qJD(2) * t873 - t753 * t902 + t754 * t903 + t872 * t940;
t752 = t882 * t851 - t883 * t850 - mrSges(4,1) * t802 + mrSges(4,2) * t803 + t952 - t872 * t941 - pkin(3) * t773 - Ifges(4,6) * t863 - Ifges(4,5) * t864 + qJD(2) * t874 - mrSges(3,1) * t875 + mrSges(3,3) * t860 + Ifges(3,4) * t888 - pkin(2) * t766 + Ifges(3,6) * qJDD(2) + (Ifges(3,2) + Ifges(4,3)) * t889;
t919 = mrSges(2,1) * t893 - mrSges(2,2) * t894 + Ifges(2,3) * qJDD(1) + pkin(1) * t915 + pkin(7) * t926 + t905 * t750 + t907 * t752;
t748 = mrSges(2,1) * g(3) + mrSges(2,3) * t894 + t910 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t759 - t953;
t747 = -mrSges(2,2) * g(3) - mrSges(2,3) * t893 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t910 - pkin(7) * t759 + t750 * t907 - t752 * t905;
t1 = [-m(1) * g(1) + t927; -m(1) * g(2) + t944; (-m(1) - m(2)) * g(3) + t759; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t944 + t908 * t747 - t906 * t748; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t927 + t906 * t747 + t908 * t748; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t919; t919; t953; t781; -t952; t782; t786;];
tauJB  = t1;
