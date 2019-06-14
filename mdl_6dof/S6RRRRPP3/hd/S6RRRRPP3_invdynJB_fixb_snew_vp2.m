% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:08:31
% EndTime: 2019-05-07 18:08:42
% DurationCPUTime: 7.08s
% Computational Cost: add. (101112->344), mult. (202553->403), div. (0->0), fcn. (139475->8), ass. (0->137)
t964 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t946 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t945 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t963 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t944 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t962 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t913 = sin(qJ(3));
t914 = sin(qJ(2));
t916 = cos(qJ(3));
t917 = cos(qJ(2));
t887 = (t913 * t917 + t914 * t916) * qJD(1);
t947 = qJD(1) * qJD(2);
t894 = qJDD(1) * t914 + t917 * t947;
t895 = qJDD(1) * t917 - t914 * t947;
t855 = -qJD(3) * t887 - t894 * t913 + t895 * t916;
t948 = qJD(1) * t917;
t949 = qJD(1) * t914;
t886 = -t913 * t949 + t916 * t948;
t856 = qJD(3) * t886 + t894 * t916 + t895 * t913;
t899 = qJD(2) * pkin(2) - pkin(8) * t949;
t911 = t917 ^ 2;
t919 = qJD(1) ^ 2;
t915 = sin(qJ(1));
t918 = cos(qJ(1));
t900 = t915 * g(1) - t918 * g(2);
t931 = -qJDD(1) * pkin(1) - t900;
t857 = -t895 * pkin(2) + t899 * t949 + (-pkin(8) * t911 - pkin(7)) * t919 + t931;
t909 = qJD(2) + qJD(3);
t804 = (-t886 * t909 - t856) * pkin(9) + (t887 * t909 - t855) * pkin(3) + t857;
t901 = -g(1) * t918 - g(2) * t915;
t889 = -pkin(1) * t919 + qJDD(1) * pkin(7) + t901;
t952 = t914 * t889;
t955 = pkin(2) * t919;
t845 = qJDD(2) * pkin(2) - t894 * pkin(8) - t952 + (pkin(8) * t947 + t914 * t955 - g(3)) * t917;
t877 = -g(3) * t914 + t917 * t889;
t846 = pkin(8) * t895 - qJD(2) * t899 - t911 * t955 + t877;
t811 = t913 * t845 + t916 * t846;
t871 = -pkin(3) * t886 - pkin(9) * t887;
t907 = t909 ^ 2;
t908 = qJDD(2) + qJDD(3);
t808 = -pkin(3) * t907 + pkin(9) * t908 + t871 * t886 + t811;
t912 = sin(qJ(4));
t957 = cos(qJ(4));
t801 = t804 * t957 - t912 * t808;
t874 = t912 * t887 - t909 * t957;
t875 = t887 * t957 + t912 * t909;
t841 = pkin(4) * t874 - qJ(5) * t875;
t854 = qJDD(4) - t855;
t882 = qJD(4) - t886;
t881 = t882 ^ 2;
t799 = -t854 * pkin(4) - t881 * qJ(5) + t875 * t841 + qJDD(5) - t801;
t818 = -t874 * qJD(4) + t856 * t957 + t912 * t908;
t843 = -mrSges(6,2) * t874 - mrSges(6,3) * t875;
t961 = -m(6) * t799 - t818 * mrSges(6,1) - t875 * t843;
t870 = -mrSges(4,1) * t886 + mrSges(4,2) * t887;
t879 = mrSges(4,1) * t909 - mrSges(4,3) * t887;
t840 = -mrSges(7,2) * t875 + mrSges(7,3) * t874;
t842 = mrSges(5,1) * t874 + mrSges(5,2) * t875;
t860 = mrSges(6,1) * t874 - mrSges(6,3) * t882;
t863 = -mrSges(5,2) * t882 - mrSges(5,3) * t874;
t953 = t874 * t882;
t793 = -0.2e1 * qJD(6) * t882 + (t874 * t875 - t854) * qJ(6) + (t818 + t953) * pkin(5) + t799;
t861 = -mrSges(7,1) * t874 + mrSges(7,2) * t882;
t934 = -m(7) * t793 + t854 * mrSges(7,3) + t882 * t861;
t954 = -mrSges(7,1) - mrSges(5,3);
t782 = m(5) * t801 + (-t860 + t863) * t882 + (-t840 - t842) * t875 + (mrSges(5,1) - mrSges(6,2)) * t854 + t954 * t818 + t934 + t961;
t802 = t912 * t804 + t957 * t808;
t817 = t875 * qJD(4) + t912 * t856 - t908 * t957;
t864 = mrSges(5,1) * t882 - mrSges(5,3) * t875;
t927 = -t881 * pkin(4) + t854 * qJ(5) - t874 * t841 + t802;
t958 = -2 * qJD(5);
t798 = t882 * t958 - t927;
t862 = mrSges(6,1) * t875 + mrSges(6,2) * t882;
t858 = pkin(5) * t875 - qJ(6) * t882;
t873 = t874 ^ 2;
t795 = -t817 * pkin(5) - t873 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t858) * t882 + t927;
t859 = mrSges(7,1) * t875 - mrSges(7,3) * t882;
t942 = m(7) * t795 + t854 * mrSges(7,2) + t882 * t859;
t930 = -m(6) * t798 + t854 * mrSges(6,3) + t882 * t862 + t942;
t950 = -t840 - t843;
t785 = m(5) * t802 - t854 * mrSges(5,2) - t882 * t864 + (-t842 + t950) * t874 + (-mrSges(6,1) + t954) * t817 + t930;
t935 = -t782 * t912 + t957 * t785;
t775 = m(4) * t811 - mrSges(4,2) * t908 + mrSges(4,3) * t855 + t870 * t886 - t879 * t909 + t935;
t810 = t916 * t845 - t913 * t846;
t878 = -mrSges(4,2) * t909 + mrSges(4,3) * t886;
t807 = -t908 * pkin(3) - t907 * pkin(9) + t887 * t871 - t810;
t924 = (-t818 + t953) * qJ(5) + t807 + (pkin(4) * t882 + t958) * t875;
t800 = t817 * pkin(4) + t924;
t797 = -t873 * pkin(5) + 0.2e1 * qJD(6) * t874 - t875 * t858 + (pkin(4) + qJ(6)) * t817 + t924;
t933 = m(7) * t797 - t818 * mrSges(7,2) + t817 * mrSges(7,3) - t875 * t859 + t874 * t861;
t928 = -m(6) * t800 + t817 * mrSges(6,2) + t874 * t860 - t933;
t923 = -m(5) * t807 - t817 * mrSges(5,1) - t874 * t863 + (t862 - t864) * t875 + (-mrSges(5,2) + mrSges(6,3)) * t818 + t928;
t780 = m(4) * t810 + t908 * mrSges(4,1) - t856 * mrSges(4,3) - t887 * t870 + t909 * t878 + t923;
t767 = t913 * t775 + t916 * t780;
t876 = -t917 * g(3) - t952;
t884 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t914 + Ifges(3,2) * t917) * qJD(1);
t885 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t914 + Ifges(3,4) * t917) * qJD(1);
t790 = -t818 * mrSges(6,3) - t875 * t862 - t928;
t792 = -t817 * mrSges(7,1) - t874 * t840 + t942;
t939 = -t946 * t874 + t875 * t964 + t945 * t882;
t941 = t874 * t944 - t875 * t945 - t882 * t962;
t761 = -mrSges(5,1) * t807 - mrSges(6,1) * t798 + mrSges(7,1) * t795 + mrSges(6,2) * t800 + mrSges(5,3) * t802 - mrSges(7,3) * t797 - pkin(4) * t790 + pkin(5) * t792 - qJ(6) * t933 + t817 * t963 + t946 * t818 + t944 * t854 + t941 * t875 + t939 * t882;
t791 = t818 * mrSges(7,1) + t875 * t840 - t934;
t940 = t874 * t963 + t875 * t946 + t882 * t944;
t769 = mrSges(6,1) * t799 + mrSges(7,1) * t793 + mrSges(5,2) * t807 - mrSges(7,2) * t797 - mrSges(5,3) * t801 - mrSges(6,3) * t800 + pkin(5) * t791 - qJ(5) * t790 - t946 * t817 + t818 * t964 + t945 * t854 + t941 * t874 - t940 * t882;
t867 = Ifges(4,4) * t887 + Ifges(4,2) * t886 + Ifges(4,6) * t909;
t868 = Ifges(4,1) * t887 + Ifges(4,4) * t886 + Ifges(4,5) * t909;
t925 = -mrSges(4,1) * t810 + mrSges(4,2) * t811 - Ifges(4,5) * t856 - Ifges(4,6) * t855 - Ifges(4,3) * t908 - pkin(3) * t923 - pkin(9) * t935 - t957 * t761 - t912 * t769 - t887 * t867 + t886 * t868;
t960 = mrSges(3,1) * t876 - mrSges(3,2) * t877 + Ifges(3,5) * t894 + Ifges(3,6) * t895 + Ifges(3,3) * qJDD(2) + pkin(2) * t767 + (t884 * t914 - t885 * t917) * qJD(1) - t925;
t788 = t854 * mrSges(6,2) + t882 * t860 + t791 - t961;
t959 = -t817 * t944 + t818 * t945 + t962 * t854 + t874 * t939 + t875 * t940 + mrSges(5,1) * t801 - mrSges(5,2) * t802 + mrSges(6,2) * t799 + mrSges(7,2) * t795 - mrSges(6,3) * t798 - mrSges(7,3) * t793 - pkin(4) * t788 + qJ(5) * (t950 * t874 + (-mrSges(6,1) - mrSges(7,1)) * t817 + t930) - qJ(6) * t791;
t893 = (-mrSges(3,1) * t917 + mrSges(3,2) * t914) * qJD(1);
t898 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t948;
t765 = m(3) * t876 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t894 + qJD(2) * t898 - t893 * t949 + t767;
t897 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t949;
t936 = t916 * t775 - t780 * t913;
t766 = m(3) * t877 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t895 - qJD(2) * t897 + t893 * t948 + t936;
t937 = -t765 * t914 + t917 * t766;
t756 = m(2) * t901 - mrSges(2,1) * t919 - qJDD(1) * mrSges(2,2) + t937;
t888 = -t919 * pkin(7) + t931;
t777 = t957 * t782 + t912 * t785;
t926 = m(4) * t857 - t855 * mrSges(4,1) + mrSges(4,2) * t856 - t886 * t878 + t879 * t887 + t777;
t922 = -m(3) * t888 + t895 * mrSges(3,1) - mrSges(3,2) * t894 - t897 * t949 + t898 * t948 - t926;
t771 = m(2) * t900 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t919 + t922;
t951 = t915 * t756 + t918 * t771;
t758 = t917 * t765 + t914 * t766;
t938 = t918 * t756 - t771 * t915;
t866 = Ifges(4,5) * t887 + Ifges(4,6) * t886 + Ifges(4,3) * t909;
t753 = mrSges(4,2) * t857 - mrSges(4,3) * t810 + Ifges(4,1) * t856 + Ifges(4,4) * t855 + Ifges(4,5) * t908 - pkin(9) * t777 - t912 * t761 + t769 * t957 + t886 * t866 - t909 * t867;
t759 = -mrSges(4,1) * t857 + mrSges(4,3) * t811 + Ifges(4,4) * t856 + Ifges(4,2) * t855 + Ifges(4,6) * t908 - pkin(3) * t777 - t887 * t866 + t909 * t868 - t959;
t883 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t914 + Ifges(3,6) * t917) * qJD(1);
t749 = -mrSges(3,1) * t888 + mrSges(3,3) * t877 + Ifges(3,4) * t894 + Ifges(3,2) * t895 + Ifges(3,6) * qJDD(2) - pkin(2) * t926 + pkin(8) * t936 + qJD(2) * t885 + t913 * t753 + t916 * t759 - t883 * t949;
t752 = mrSges(3,2) * t888 - mrSges(3,3) * t876 + Ifges(3,1) * t894 + Ifges(3,4) * t895 + Ifges(3,5) * qJDD(2) - pkin(8) * t767 - qJD(2) * t884 + t753 * t916 - t759 * t913 + t883 * t948;
t929 = mrSges(2,1) * t900 - mrSges(2,2) * t901 + Ifges(2,3) * qJDD(1) + pkin(1) * t922 + pkin(7) * t937 + t917 * t749 + t914 * t752;
t750 = mrSges(2,1) * g(3) + mrSges(2,3) * t901 + t919 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t758 - t960;
t747 = -mrSges(2,2) * g(3) - mrSges(2,3) * t900 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t919 - pkin(7) * t758 - t749 * t914 + t752 * t917;
t1 = [-m(1) * g(1) + t938; -m(1) * g(2) + t951; (-m(1) - m(2)) * g(3) + t758; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t951 + t918 * t747 - t915 * t750; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t938 + t915 * t747 + t918 * t750; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t929; t929; t960; -t925; t959; t788; t792;];
tauJB  = t1;
