% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:09:38
% EndTime: 2019-05-07 20:10:33
% DurationCPUTime: 40.39s
% Computational Cost: add. (686362->389), mult. (1386197->489), div. (0->0), fcn. (1025063->12), ass. (0->154)
t957 = sin(qJ(3));
t958 = sin(qJ(2));
t962 = cos(qJ(3));
t963 = cos(qJ(2));
t929 = (t957 * t963 + t958 * t962) * qJD(1);
t983 = qJD(1) * qJD(2);
t937 = qJDD(1) * t958 + t963 * t983;
t938 = qJDD(1) * t963 - t958 * t983;
t902 = -t929 * qJD(3) - t957 * t937 + t938 * t962;
t984 = qJD(1) * t963;
t985 = qJD(1) * t958;
t928 = -t957 * t985 + t962 * t984;
t903 = qJD(3) * t928 + t937 * t962 + t938 * t957;
t942 = qJD(2) * pkin(2) - pkin(8) * t985;
t952 = t963 ^ 2;
t965 = qJD(1) ^ 2;
t959 = sin(qJ(1));
t964 = cos(qJ(1));
t943 = t959 * g(1) - t964 * g(2);
t974 = -qJDD(1) * pkin(1) - t943;
t904 = -t938 * pkin(2) + t942 * t985 + (-pkin(8) * t952 - pkin(7)) * t965 + t974;
t950 = qJD(2) + qJD(3);
t853 = (-t928 * t950 - t903) * pkin(9) + (t929 * t950 - t902) * pkin(3) + t904;
t944 = -g(1) * t964 - g(2) * t959;
t931 = -pkin(1) * t965 + qJDD(1) * pkin(7) + t944;
t987 = t958 * t931;
t989 = pkin(2) * t965;
t893 = qJDD(2) * pkin(2) - t937 * pkin(8) - t987 + (pkin(8) * t983 + t958 * t989 - g(3)) * t963;
t918 = -g(3) * t958 + t963 * t931;
t894 = pkin(8) * t938 - qJD(2) * t942 - t952 * t989 + t918;
t871 = t957 * t893 + t962 * t894;
t913 = -pkin(3) * t928 - pkin(9) * t929;
t948 = t950 ^ 2;
t949 = qJDD(2) + qJDD(3);
t859 = -pkin(3) * t948 + pkin(9) * t949 + t913 * t928 + t871;
t956 = sin(qJ(4));
t961 = cos(qJ(4));
t842 = t961 * t853 - t956 * t859;
t915 = -t929 * t956 + t950 * t961;
t874 = qJD(4) * t915 + t903 * t961 + t949 * t956;
t901 = qJDD(4) - t902;
t916 = t929 * t961 + t950 * t956;
t924 = qJD(4) - t928;
t832 = (t915 * t924 - t874) * qJ(5) + (t915 * t916 + t901) * pkin(4) + t842;
t843 = t956 * t853 + t961 * t859;
t873 = -qJD(4) * t916 - t903 * t956 + t949 * t961;
t906 = pkin(4) * t924 - qJ(5) * t916;
t914 = t915 ^ 2;
t834 = -pkin(4) * t914 + qJ(5) * t873 - t906 * t924 + t843;
t953 = sin(pkin(11));
t954 = cos(pkin(11));
t887 = t915 * t953 + t916 * t954;
t826 = -0.2e1 * qJD(5) * t887 + t954 * t832 - t953 * t834;
t856 = t873 * t953 + t874 * t954;
t886 = t915 * t954 - t916 * t953;
t824 = (t886 * t924 - t856) * pkin(10) + (t886 * t887 + t901) * pkin(5) + t826;
t827 = 0.2e1 * qJD(5) * t886 + t953 * t832 + t954 * t834;
t855 = t873 * t954 - t874 * t953;
t877 = pkin(5) * t924 - pkin(10) * t887;
t885 = t886 ^ 2;
t825 = -pkin(5) * t885 + pkin(10) * t855 - t877 * t924 + t827;
t955 = sin(qJ(6));
t960 = cos(qJ(6));
t822 = t824 * t960 - t825 * t955;
t866 = t886 * t960 - t887 * t955;
t840 = qJD(6) * t866 + t855 * t955 + t856 * t960;
t867 = t886 * t955 + t887 * t960;
t850 = -mrSges(7,1) * t866 + mrSges(7,2) * t867;
t922 = qJD(6) + t924;
t860 = -mrSges(7,2) * t922 + mrSges(7,3) * t866;
t896 = qJDD(6) + t901;
t814 = m(7) * t822 + mrSges(7,1) * t896 - mrSges(7,3) * t840 - t850 * t867 + t860 * t922;
t823 = t824 * t955 + t825 * t960;
t839 = -qJD(6) * t867 + t855 * t960 - t856 * t955;
t861 = mrSges(7,1) * t922 - mrSges(7,3) * t867;
t815 = m(7) * t823 - mrSges(7,2) * t896 + mrSges(7,3) * t839 + t850 * t866 - t861 * t922;
t808 = t960 * t814 + t955 * t815;
t868 = -mrSges(6,1) * t886 + mrSges(6,2) * t887;
t875 = -mrSges(6,2) * t924 + mrSges(6,3) * t886;
t806 = m(6) * t826 + mrSges(6,1) * t901 - mrSges(6,3) * t856 - t868 * t887 + t875 * t924 + t808;
t876 = mrSges(6,1) * t924 - mrSges(6,3) * t887;
t977 = -t814 * t955 + t960 * t815;
t807 = m(6) * t827 - mrSges(6,2) * t901 + mrSges(6,3) * t855 + t868 * t886 - t876 * t924 + t977;
t802 = t954 * t806 + t953 * t807;
t863 = Ifges(6,4) * t887 + Ifges(6,2) * t886 + Ifges(6,6) * t924;
t864 = Ifges(6,1) * t887 + Ifges(6,4) * t886 + Ifges(6,5) * t924;
t879 = Ifges(5,4) * t916 + Ifges(5,2) * t915 + Ifges(5,6) * t924;
t880 = Ifges(5,1) * t916 + Ifges(5,4) * t915 + Ifges(5,5) * t924;
t846 = Ifges(7,4) * t867 + Ifges(7,2) * t866 + Ifges(7,6) * t922;
t847 = Ifges(7,1) * t867 + Ifges(7,4) * t866 + Ifges(7,5) * t922;
t972 = -mrSges(7,1) * t822 + mrSges(7,2) * t823 - Ifges(7,5) * t840 - Ifges(7,6) * t839 - Ifges(7,3) * t896 - t867 * t846 + t866 * t847;
t991 = mrSges(5,1) * t842 + mrSges(6,1) * t826 - mrSges(5,2) * t843 - mrSges(6,2) * t827 + Ifges(5,5) * t874 + Ifges(6,5) * t856 + Ifges(5,6) * t873 + Ifges(6,6) * t855 + pkin(4) * t802 + pkin(5) * t808 + t887 * t863 - t886 * t864 + t916 * t879 - t915 * t880 - (-Ifges(6,3) - Ifges(5,3)) * t901 - t972;
t912 = -mrSges(4,1) * t928 + mrSges(4,2) * t929;
t920 = mrSges(4,1) * t950 - mrSges(4,3) * t929;
t891 = -mrSges(5,1) * t915 + mrSges(5,2) * t916;
t905 = -mrSges(5,2) * t924 + mrSges(5,3) * t915;
t800 = m(5) * t842 + mrSges(5,1) * t901 - mrSges(5,3) * t874 - t891 * t916 + t905 * t924 + t802;
t907 = mrSges(5,1) * t924 - mrSges(5,3) * t916;
t978 = -t806 * t953 + t954 * t807;
t801 = m(5) * t843 - mrSges(5,2) * t901 + mrSges(5,3) * t873 + t891 * t915 - t907 * t924 + t978;
t979 = -t800 * t956 + t961 * t801;
t791 = m(4) * t871 - mrSges(4,2) * t949 + mrSges(4,3) * t902 + t912 * t928 - t920 * t950 + t979;
t870 = t893 * t962 - t957 * t894;
t919 = -mrSges(4,2) * t950 + mrSges(4,3) * t928;
t858 = -pkin(3) * t949 - pkin(9) * t948 + t929 * t913 - t870;
t844 = -pkin(4) * t873 - qJ(5) * t914 + t916 * t906 + qJDD(5) + t858;
t829 = -pkin(5) * t855 - pkin(10) * t885 + t877 * t887 + t844;
t976 = m(7) * t829 - t839 * mrSges(7,1) + t840 * mrSges(7,2) - t866 * t860 + t867 * t861;
t820 = m(6) * t844 - t855 * mrSges(6,1) + mrSges(6,2) * t856 - t886 * t875 + t876 * t887 + t976;
t968 = -m(5) * t858 + t873 * mrSges(5,1) - mrSges(5,2) * t874 + t915 * t905 - t907 * t916 - t820;
t817 = m(4) * t870 + mrSges(4,1) * t949 - mrSges(4,3) * t903 - t912 * t929 + t919 * t950 + t968;
t785 = t957 * t791 + t962 * t817;
t917 = -t963 * g(3) - t987;
t926 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t958 + Ifges(3,2) * t963) * qJD(1);
t927 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t958 + Ifges(3,4) * t963) * qJD(1);
t845 = Ifges(7,5) * t867 + Ifges(7,6) * t866 + Ifges(7,3) * t922;
t809 = -mrSges(7,1) * t829 + mrSges(7,3) * t823 + Ifges(7,4) * t840 + Ifges(7,2) * t839 + Ifges(7,6) * t896 - t845 * t867 + t847 * t922;
t810 = mrSges(7,2) * t829 - mrSges(7,3) * t822 + Ifges(7,1) * t840 + Ifges(7,4) * t839 + Ifges(7,5) * t896 + t845 * t866 - t846 * t922;
t862 = Ifges(6,5) * t887 + Ifges(6,6) * t886 + Ifges(6,3) * t924;
t795 = -mrSges(6,1) * t844 + mrSges(6,3) * t827 + Ifges(6,4) * t856 + Ifges(6,2) * t855 + Ifges(6,6) * t901 - pkin(5) * t976 + pkin(10) * t977 + t960 * t809 + t955 * t810 - t887 * t862 + t924 * t864;
t796 = mrSges(6,2) * t844 - mrSges(6,3) * t826 + Ifges(6,1) * t856 + Ifges(6,4) * t855 + Ifges(6,5) * t901 - pkin(10) * t808 - t809 * t955 + t810 * t960 + t862 * t886 - t863 * t924;
t878 = Ifges(5,5) * t916 + Ifges(5,6) * t915 + Ifges(5,3) * t924;
t777 = -mrSges(5,1) * t858 + mrSges(5,3) * t843 + Ifges(5,4) * t874 + Ifges(5,2) * t873 + Ifges(5,6) * t901 - pkin(4) * t820 + qJ(5) * t978 + t954 * t795 + t953 * t796 - t916 * t878 + t924 * t880;
t779 = mrSges(5,2) * t858 - mrSges(5,3) * t842 + Ifges(5,1) * t874 + Ifges(5,4) * t873 + Ifges(5,5) * t901 - qJ(5) * t802 - t795 * t953 + t796 * t954 + t878 * t915 - t879 * t924;
t909 = Ifges(4,4) * t929 + Ifges(4,2) * t928 + Ifges(4,6) * t950;
t910 = Ifges(4,1) * t929 + Ifges(4,4) * t928 + Ifges(4,5) * t950;
t970 = -mrSges(4,1) * t870 + mrSges(4,2) * t871 - Ifges(4,5) * t903 - Ifges(4,6) * t902 - Ifges(4,3) * t949 - pkin(3) * t968 - pkin(9) * t979 - t961 * t777 - t956 * t779 - t929 * t909 + t928 * t910;
t990 = mrSges(3,1) * t917 - mrSges(3,2) * t918 + Ifges(3,5) * t937 + Ifges(3,6) * t938 + Ifges(3,3) * qJDD(2) + pkin(2) * t785 + (t926 * t958 - t927 * t963) * qJD(1) - t970;
t936 = (-mrSges(3,1) * t963 + mrSges(3,2) * t958) * qJD(1);
t941 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t984;
t783 = m(3) * t917 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t937 + qJD(2) * t941 - t936 * t985 + t785;
t940 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t985;
t980 = t962 * t791 - t817 * t957;
t784 = m(3) * t918 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t938 - qJD(2) * t940 + t936 * t984 + t980;
t981 = -t783 * t958 + t963 * t784;
t772 = m(2) * t944 - mrSges(2,1) * t965 - qJDD(1) * mrSges(2,2) + t981;
t930 = -t965 * pkin(7) + t974;
t793 = t961 * t800 + t956 * t801;
t971 = m(4) * t904 - t902 * mrSges(4,1) + mrSges(4,2) * t903 - t928 * t919 + t920 * t929 + t793;
t969 = -m(3) * t930 + t938 * mrSges(3,1) - mrSges(3,2) * t937 - t940 * t985 + t941 * t984 - t971;
t787 = m(2) * t943 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t965 + t969;
t986 = t959 * t772 + t964 * t787;
t774 = t963 * t783 + t958 * t784;
t982 = t964 * t772 - t787 * t959;
t908 = Ifges(4,5) * t929 + Ifges(4,6) * t928 + Ifges(4,3) * t950;
t769 = mrSges(4,2) * t904 - mrSges(4,3) * t870 + Ifges(4,1) * t903 + Ifges(4,4) * t902 + Ifges(4,5) * t949 - pkin(9) * t793 - t777 * t956 + t779 * t961 + t908 * t928 - t909 * t950;
t775 = -mrSges(4,1) * t904 + mrSges(4,3) * t871 + Ifges(4,4) * t903 + Ifges(4,2) * t902 + Ifges(4,6) * t949 - pkin(3) * t793 - t929 * t908 + t950 * t910 - t991;
t925 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t958 + Ifges(3,6) * t963) * qJD(1);
t765 = -mrSges(3,1) * t930 + mrSges(3,3) * t918 + Ifges(3,4) * t937 + Ifges(3,2) * t938 + Ifges(3,6) * qJDD(2) - pkin(2) * t971 + pkin(8) * t980 + qJD(2) * t927 + t957 * t769 + t962 * t775 - t925 * t985;
t768 = mrSges(3,2) * t930 - mrSges(3,3) * t917 + Ifges(3,1) * t937 + Ifges(3,4) * t938 + Ifges(3,5) * qJDD(2) - pkin(8) * t785 - qJD(2) * t926 + t769 * t962 - t775 * t957 + t925 * t984;
t973 = mrSges(2,1) * t943 - mrSges(2,2) * t944 + Ifges(2,3) * qJDD(1) + pkin(1) * t969 + pkin(7) * t981 + t963 * t765 + t958 * t768;
t766 = mrSges(2,1) * g(3) + mrSges(2,3) * t944 + t965 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t774 - t990;
t763 = -mrSges(2,2) * g(3) - mrSges(2,3) * t943 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t965 - pkin(7) * t774 - t765 * t958 + t768 * t963;
t1 = [-m(1) * g(1) + t982; -m(1) * g(2) + t986; (-m(1) - m(2)) * g(3) + t774; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t986 + t964 * t763 - t959 * t766; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t982 + t959 * t763 + t964 * t766; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t973; t973; t990; -t970; t991; t820; -t972;];
tauJB  = t1;
