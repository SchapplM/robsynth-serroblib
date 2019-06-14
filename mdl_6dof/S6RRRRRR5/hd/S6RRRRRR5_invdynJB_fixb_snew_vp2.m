% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 10:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 10:09:42
% EndTime: 2019-05-08 10:11:38
% DurationCPUTime: 81.80s
% Computational Cost: add. (1399498->403), mult. (3032122->517), div. (0->0), fcn. (2503429->14), ass. (0->169)
t1003 = sin(pkin(6));
t1044 = pkin(8) * t1003;
t1004 = cos(pkin(6));
t1043 = t1004 * g(3);
t1010 = sin(qJ(1));
t1016 = cos(qJ(1));
t1017 = qJD(1) ^ 2;
t1015 = cos(qJ(2));
t1035 = t1004 * t1015;
t1009 = sin(qJ(2));
t1036 = t1004 * t1009;
t1008 = sin(qJ(3));
t1014 = cos(qJ(3));
t1007 = sin(qJ(4));
t1013 = cos(qJ(4));
t1006 = sin(qJ(5));
t1012 = cos(qJ(5));
t1005 = sin(qJ(6));
t1011 = cos(qJ(6));
t1039 = qJD(1) * t1015;
t994 = t1010 * g(1) - g(2) * t1016;
t983 = qJDD(1) * pkin(1) + t1017 * t1044 + t994;
t1034 = qJDD(1) * t1003;
t995 = -g(1) * t1016 - g(2) * t1010;
t984 = -pkin(1) * t1017 + pkin(8) * t1034 + t995;
t1041 = t1015 * t984 + t983 * t1036;
t1040 = qJD(1) * t1003;
t986 = (-pkin(2) * t1015 - pkin(9) * t1009) * t1040;
t999 = qJD(1) * t1004 + qJD(2);
t997 = t999 ^ 2;
t998 = qJDD(1) * t1004 + qJDD(2);
t943 = -t997 * pkin(2) + t998 * pkin(9) + (-g(3) * t1009 + t1039 * t986) * t1003 + t1041;
t987 = (qJD(2) * t1039 + qJDD(1) * t1009) * t1003;
t1038 = t1003 * t1009;
t1033 = qJD(1) * t1038;
t988 = -qJD(2) * t1033 + t1015 * t1034;
t944 = -t988 * pkin(2) - t987 * pkin(9) - t1043 + (-t983 + (pkin(2) * t1009 - pkin(9) * t1015) * t999 * qJD(1)) * t1003;
t918 = -t1008 * t943 + t1014 * t944;
t975 = -t1008 * t1033 + t1014 * t999;
t955 = qJD(3) * t975 + t1008 * t998 + t1014 * t987;
t976 = t1008 * t999 + t1014 * t1033;
t980 = qJDD(3) - t988;
t1037 = t1003 * t1015;
t1032 = qJD(1) * t1037;
t993 = qJD(3) - t1032;
t903 = (t975 * t993 - t955) * pkin(10) + (t975 * t976 + t980) * pkin(3) + t918;
t919 = t1008 * t944 + t1014 * t943;
t954 = -qJD(3) * t976 - t1008 * t987 + t1014 * t998;
t964 = pkin(3) * t993 - pkin(10) * t976;
t974 = t975 ^ 2;
t906 = -pkin(3) * t974 + pkin(10) * t954 - t964 * t993 + t919;
t883 = -t1007 * t906 + t1013 * t903;
t959 = -t1007 * t976 + t1013 * t975;
t924 = qJD(4) * t959 + t1007 * t954 + t1013 * t955;
t960 = t1007 * t975 + t1013 * t976;
t979 = qJDD(4) + t980;
t991 = qJD(4) + t993;
t879 = (t959 * t991 - t924) * pkin(11) + (t959 * t960 + t979) * pkin(4) + t883;
t884 = t1007 * t903 + t1013 * t906;
t923 = -qJD(4) * t960 - t1007 * t955 + t1013 * t954;
t947 = pkin(4) * t991 - pkin(11) * t960;
t958 = t959 ^ 2;
t881 = -pkin(4) * t958 + pkin(11) * t923 - t947 * t991 + t884;
t876 = t1006 * t879 + t1012 * t881;
t936 = -t1006 * t960 + t1012 * t959;
t937 = t1006 * t959 + t1012 * t960;
t917 = -pkin(5) * t936 - pkin(12) * t937;
t973 = qJDD(5) + t979;
t990 = qJD(5) + t991;
t989 = t990 ^ 2;
t873 = -pkin(5) * t989 + pkin(12) * t973 + t917 * t936 + t876;
t956 = -g(3) * t1037 - t1009 * t984 + t1035 * t983;
t942 = -t998 * pkin(2) - t997 * pkin(9) + t986 * t1033 - t956;
t910 = -t954 * pkin(3) - t974 * pkin(10) + t976 * t964 + t942;
t889 = -t923 * pkin(4) - t958 * pkin(11) + t960 * t947 + t910;
t895 = -qJD(5) * t937 - t1006 * t924 + t1012 * t923;
t896 = qJD(5) * t936 + t1006 * t923 + t1012 * t924;
t877 = t889 + (-t936 * t990 - t896) * pkin(12) + (t937 * t990 - t895) * pkin(5);
t870 = -t1005 * t873 + t1011 * t877;
t926 = -t1005 * t937 + t1011 * t990;
t887 = qJD(6) * t926 + t1005 * t973 + t1011 * t896;
t894 = qJDD(6) - t895;
t927 = t1005 * t990 + t1011 * t937;
t907 = -mrSges(7,1) * t926 + mrSges(7,2) * t927;
t935 = qJD(6) - t936;
t908 = -mrSges(7,2) * t935 + mrSges(7,3) * t926;
t866 = m(7) * t870 + mrSges(7,1) * t894 - mrSges(7,3) * t887 - t907 * t927 + t908 * t935;
t871 = t1005 * t877 + t1011 * t873;
t886 = -qJD(6) * t927 - t1005 * t896 + t1011 * t973;
t909 = mrSges(7,1) * t935 - mrSges(7,3) * t927;
t867 = m(7) * t871 - mrSges(7,2) * t894 + mrSges(7,3) * t886 + t907 * t926 - t909 * t935;
t1031 = -t1005 * t866 + t1011 * t867;
t916 = -mrSges(6,1) * t936 + mrSges(6,2) * t937;
t929 = mrSges(6,1) * t990 - mrSges(6,3) * t937;
t853 = m(6) * t876 - mrSges(6,2) * t973 + mrSges(6,3) * t895 + t916 * t936 - t929 * t990 + t1031;
t875 = -t1006 * t881 + t1012 * t879;
t872 = -pkin(5) * t973 - pkin(12) * t989 + t917 * t937 - t875;
t1025 = -m(7) * t872 + t886 * mrSges(7,1) - mrSges(7,2) * t887 + t926 * t908 - t909 * t927;
t928 = -mrSges(6,2) * t990 + mrSges(6,3) * t936;
t862 = m(6) * t875 + mrSges(6,1) * t973 - mrSges(6,3) * t896 - t916 * t937 + t928 * t990 + t1025;
t847 = t1006 * t853 + t1012 * t862;
t938 = -mrSges(5,1) * t959 + mrSges(5,2) * t960;
t945 = -mrSges(5,2) * t991 + mrSges(5,3) * t959;
t844 = m(5) * t883 + mrSges(5,1) * t979 - mrSges(5,3) * t924 - t938 * t960 + t945 * t991 + t847;
t1030 = -t1006 * t862 + t1012 * t853;
t946 = mrSges(5,1) * t991 - mrSges(5,3) * t960;
t845 = m(5) * t884 - mrSges(5,2) * t979 + mrSges(5,3) * t923 + t938 * t959 - t946 * t991 + t1030;
t838 = t1007 * t845 + t1013 * t844;
t961 = -mrSges(4,1) * t975 + mrSges(4,2) * t976;
t962 = -mrSges(4,2) * t993 + mrSges(4,3) * t975;
t836 = m(4) * t918 + mrSges(4,1) * t980 - mrSges(4,3) * t955 - t961 * t976 + t962 * t993 + t838;
t1029 = -t1007 * t844 + t1013 * t845;
t963 = mrSges(4,1) * t993 - mrSges(4,3) * t976;
t837 = m(4) * t919 - mrSges(4,2) * t980 + mrSges(4,3) * t954 + t961 * t975 - t963 * t993 + t1029;
t1028 = -t1008 * t836 + t1014 * t837;
t957 = -g(3) * t1038 + t1041;
t981 = mrSges(3,1) * t999 - mrSges(3,3) * t1033;
t985 = (-mrSges(3,1) * t1015 + mrSges(3,2) * t1009) * t1040;
t828 = m(3) * t957 - mrSges(3,2) * t998 + mrSges(3,3) * t988 + t1032 * t985 - t981 * t999 + t1028;
t831 = t1008 * t837 + t1014 * t836;
t968 = -t1003 * t983 - t1043;
t982 = -mrSges(3,2) * t999 + mrSges(3,3) * t1032;
t830 = m(3) * t968 - t988 * mrSges(3,1) + t987 * mrSges(3,2) + (t1009 * t981 - t1015 * t982) * t1040 + t831;
t855 = t1005 * t867 + t1011 * t866;
t1026 = m(6) * t889 - t895 * mrSges(6,1) + t896 * mrSges(6,2) - t936 * t928 + t937 * t929 + t855;
t1022 = m(5) * t910 - t923 * mrSges(5,1) + t924 * mrSges(5,2) - t959 * t945 + t960 * t946 + t1026;
t1019 = -m(4) * t942 + t954 * mrSges(4,1) - t955 * mrSges(4,2) + t975 * t962 - t976 * t963 - t1022;
t850 = m(3) * t956 + t998 * mrSges(3,1) - t987 * mrSges(3,3) - t1033 * t985 + t999 * t982 + t1019;
t818 = -t1003 * t830 + t850 * t1035 + t828 * t1036;
t815 = m(2) * t994 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1017 + t818;
t823 = -t1009 * t850 + t1015 * t828;
t821 = m(2) * t995 - mrSges(2,1) * t1017 - qJDD(1) * mrSges(2,2) + t823;
t1042 = t1010 * t821 + t1016 * t815;
t817 = t1004 * t830 + t850 * t1037 + t828 * t1038;
t1027 = -t1010 * t815 + t1016 * t821;
t897 = Ifges(7,5) * t927 + Ifges(7,6) * t926 + Ifges(7,3) * t935;
t899 = Ifges(7,1) * t927 + Ifges(7,4) * t926 + Ifges(7,5) * t935;
t859 = -mrSges(7,1) * t872 + mrSges(7,3) * t871 + Ifges(7,4) * t887 + Ifges(7,2) * t886 + Ifges(7,6) * t894 - t897 * t927 + t899 * t935;
t898 = Ifges(7,4) * t927 + Ifges(7,2) * t926 + Ifges(7,6) * t935;
t860 = mrSges(7,2) * t872 - mrSges(7,3) * t870 + Ifges(7,1) * t887 + Ifges(7,4) * t886 + Ifges(7,5) * t894 + t897 * t926 - t898 * t935;
t911 = Ifges(6,5) * t937 + Ifges(6,6) * t936 + Ifges(6,3) * t990;
t912 = Ifges(6,4) * t937 + Ifges(6,2) * t936 + Ifges(6,6) * t990;
t839 = mrSges(6,2) * t889 - mrSges(6,3) * t875 + Ifges(6,1) * t896 + Ifges(6,4) * t895 + Ifges(6,5) * t973 - pkin(12) * t855 - t1005 * t859 + t1011 * t860 + t911 * t936 - t912 * t990;
t1021 = mrSges(7,1) * t870 - mrSges(7,2) * t871 + Ifges(7,5) * t887 + Ifges(7,6) * t886 + Ifges(7,3) * t894 + t898 * t927 - t899 * t926;
t913 = Ifges(6,1) * t937 + Ifges(6,4) * t936 + Ifges(6,5) * t990;
t840 = -mrSges(6,1) * t889 + mrSges(6,3) * t876 + Ifges(6,4) * t896 + Ifges(6,2) * t895 + Ifges(6,6) * t973 - pkin(5) * t855 - t911 * t937 + t913 * t990 - t1021;
t930 = Ifges(5,5) * t960 + Ifges(5,6) * t959 + Ifges(5,3) * t991;
t932 = Ifges(5,1) * t960 + Ifges(5,4) * t959 + Ifges(5,5) * t991;
t824 = -mrSges(5,1) * t910 + mrSges(5,3) * t884 + Ifges(5,4) * t924 + Ifges(5,2) * t923 + Ifges(5,6) * t979 - pkin(4) * t1026 + pkin(11) * t1030 + t1006 * t839 + t1012 * t840 - t960 * t930 + t991 * t932;
t931 = Ifges(5,4) * t960 + Ifges(5,2) * t959 + Ifges(5,6) * t991;
t832 = mrSges(5,2) * t910 - mrSges(5,3) * t883 + Ifges(5,1) * t924 + Ifges(5,4) * t923 + Ifges(5,5) * t979 - pkin(11) * t847 - t1006 * t840 + t1012 * t839 + t930 * t959 - t931 * t991;
t948 = Ifges(4,5) * t976 + Ifges(4,6) * t975 + Ifges(4,3) * t993;
t950 = Ifges(4,1) * t976 + Ifges(4,4) * t975 + Ifges(4,5) * t993;
t810 = -mrSges(4,1) * t942 + mrSges(4,3) * t919 + Ifges(4,4) * t955 + Ifges(4,2) * t954 + Ifges(4,6) * t980 - pkin(3) * t1022 + pkin(10) * t1029 + t1007 * t832 + t1013 * t824 - t976 * t948 + t993 * t950;
t949 = Ifges(4,4) * t976 + Ifges(4,2) * t975 + Ifges(4,6) * t993;
t811 = mrSges(4,2) * t942 - mrSges(4,3) * t918 + Ifges(4,1) * t955 + Ifges(4,4) * t954 + Ifges(4,5) * t980 - pkin(10) * t838 - t1007 * t824 + t1013 * t832 + t948 * t975 - t949 * t993;
t966 = Ifges(3,6) * t999 + (Ifges(3,4) * t1009 + Ifges(3,2) * t1015) * t1040;
t967 = Ifges(3,5) * t999 + (Ifges(3,1) * t1009 + Ifges(3,4) * t1015) * t1040;
t807 = Ifges(3,5) * t987 + Ifges(3,6) * t988 + Ifges(3,3) * t998 + mrSges(3,1) * t956 - mrSges(3,2) * t957 + t1008 * t811 + t1014 * t810 + pkin(2) * t1019 + pkin(9) * t1028 + (t1009 * t966 - t1015 * t967) * t1040;
t965 = Ifges(3,3) * t999 + (Ifges(3,5) * t1009 + Ifges(3,6) * t1015) * t1040;
t809 = mrSges(3,2) * t968 - mrSges(3,3) * t956 + Ifges(3,1) * t987 + Ifges(3,4) * t988 + Ifges(3,5) * t998 - pkin(9) * t831 - t1008 * t810 + t1014 * t811 + t1032 * t965 - t966 * t999;
t1023 = -mrSges(6,1) * t875 + mrSges(6,2) * t876 - Ifges(6,5) * t896 - Ifges(6,6) * t895 - Ifges(6,3) * t973 - pkin(5) * t1025 - pkin(12) * t1031 - t1005 * t860 - t1011 * t859 - t937 * t912 + t936 * t913;
t1020 = -mrSges(5,1) * t883 + mrSges(5,2) * t884 - Ifges(5,5) * t924 - Ifges(5,6) * t923 - Ifges(5,3) * t979 - pkin(4) * t847 - t960 * t931 + t959 * t932 + t1023;
t1018 = mrSges(4,1) * t918 - mrSges(4,2) * t919 + Ifges(4,5) * t955 + Ifges(4,6) * t954 + Ifges(4,3) * t980 + pkin(3) * t838 + t976 * t949 - t975 * t950 - t1020;
t813 = -mrSges(3,1) * t968 + mrSges(3,3) * t957 + Ifges(3,4) * t987 + Ifges(3,2) * t988 + Ifges(3,6) * t998 - pkin(2) * t831 - t1033 * t965 + t999 * t967 - t1018;
t1024 = mrSges(2,1) * t994 - mrSges(2,2) * t995 + Ifges(2,3) * qJDD(1) + pkin(1) * t818 + t1004 * t807 + t813 * t1037 + t809 * t1038 + t823 * t1044;
t805 = -mrSges(2,2) * g(3) - mrSges(2,3) * t994 + Ifges(2,5) * qJDD(1) - t1017 * Ifges(2,6) - t1009 * t813 + t1015 * t809 + (-t1003 * t817 - t1004 * t818) * pkin(8);
t804 = mrSges(2,1) * g(3) + mrSges(2,3) * t995 + t1017 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t817 - t1003 * t807 + (pkin(8) * t823 + t1009 * t809 + t1015 * t813) * t1004;
t1 = [-m(1) * g(1) + t1027; -m(1) * g(2) + t1042; (-m(1) - m(2)) * g(3) + t817; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1042 - t1010 * t804 + t1016 * t805; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1027 + t1010 * t805 + t1016 * t804; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1024; t1024; t807; t1018; -t1020; -t1023; t1021;];
tauJB  = t1;
