% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 00:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 00:20:11
% EndTime: 2019-05-08 00:23:05
% DurationCPUTime: 169.57s
% Computational Cost: add. (2888568->419), mult. (7157883->557), div. (0->0), fcn. (6078626->16), ass. (0->182)
t1036 = sin(pkin(7));
t1039 = cos(pkin(7));
t1043 = sin(qJ(3));
t1048 = cos(qJ(3));
t1040 = cos(pkin(6));
t1032 = qJD(1) * t1040 + qJD(2);
t1037 = sin(pkin(6));
t1049 = cos(qJ(2));
t1074 = t1037 * t1049;
t1067 = qJD(1) * t1074;
t1016 = (t1032 * t1036 + t1039 * t1067) * pkin(10);
t1044 = sin(qJ(2));
t1080 = qJD(1) * t1037;
t1084 = pkin(10) * t1036;
t1020 = (-pkin(2) * t1049 - t1044 * t1084) * t1080;
t1078 = qJD(1) * t1049;
t1026 = (qJD(2) * t1078 + qJDD(1) * t1044) * t1037;
t1031 = qJDD(1) * t1040 + qJDD(2);
t1045 = sin(qJ(1));
t1050 = cos(qJ(1));
t1029 = t1045 * g(1) - g(2) * t1050;
t1051 = qJD(1) ^ 2;
t1085 = pkin(9) * t1037;
t1023 = qJDD(1) * pkin(1) + t1051 * t1085 + t1029;
t1030 = -g(1) * t1050 - g(2) * t1045;
t1024 = -pkin(1) * t1051 + qJDD(1) * t1085 + t1030;
t1070 = t1040 * t1049;
t1062 = t1023 * t1070 - t1044 * t1024;
t1079 = qJD(1) * t1044;
t1083 = pkin(10) * t1039;
t977 = -t1026 * t1083 + t1031 * pkin(2) + t1032 * t1016 + (-g(3) * t1049 - t1020 * t1079) * t1037 + t1062;
t1075 = t1037 * t1044;
t1068 = qJD(1) * t1075;
t1019 = pkin(2) * t1032 - t1068 * t1083;
t1027 = (-qJD(2) * t1079 + qJDD(1) * t1049) * t1037;
t1058 = t1027 * t1039 + t1031 * t1036;
t1071 = t1040 * t1044;
t1069 = t1023 * t1071 + t1049 * t1024;
t978 = -t1032 * t1019 + (-g(3) * t1044 + t1020 * t1078) * t1037 + t1058 * pkin(10) + t1069;
t1082 = t1040 * g(3);
t984 = -t1026 * t1084 - t1027 * pkin(2) - t1082 + (-t1023 + (-t1016 * t1049 + t1019 * t1044) * qJD(1)) * t1037;
t944 = -t1043 * t978 + (t1036 * t984 + t1039 * t977) * t1048;
t1041 = sin(qJ(6));
t1046 = cos(qJ(6));
t1035 = sin(pkin(13));
t1038 = cos(pkin(13));
t1072 = t1039 * t1049;
t1076 = t1036 * t1048;
t1006 = (-t1043 * t1044 + t1048 * t1072) * t1080 + t1032 * t1076;
t1005 = qJD(4) - t1006;
t1004 = t1005 ^ 2;
t1087 = 2 * qJD(5);
t1042 = sin(qJ(4));
t1047 = cos(qJ(4));
t1008 = -t1027 * t1036 + t1031 * t1039 + qJDD(3);
t1017 = t1032 * t1039 - t1036 * t1067 + qJD(3);
t1015 = t1017 ^ 2;
t1073 = t1039 * t1043;
t1077 = t1036 * t1043;
t945 = t1048 * t978 + t977 * t1073 + t984 * t1077;
t1007 = t1032 * t1077 + (t1043 * t1072 + t1044 * t1048) * t1080;
t996 = -pkin(3) * t1006 - pkin(11) * t1007;
t936 = -pkin(3) * t1015 + pkin(11) * t1008 + t1006 * t996 + t945;
t959 = -t1036 * t977 + t1039 * t984;
t993 = -t1007 * qJD(3) - t1043 * t1026 + t1048 * t1058;
t994 = t1006 * qJD(3) + t1048 * t1026 + t1043 * t1058;
t939 = (-t1006 * t1017 - t994) * pkin(11) + (t1007 * t1017 - t993) * pkin(3) + t959;
t928 = -t1042 * t936 + t1047 * t939;
t998 = -t1007 * t1042 + t1017 * t1047;
t962 = qJD(4) * t998 + t1008 * t1042 + t1047 * t994;
t992 = qJDD(4) - t993;
t999 = t1007 * t1047 + t1017 * t1042;
t925 = (t1005 * t998 - t962) * qJ(5) + (t998 * t999 + t992) * pkin(4) + t928;
t929 = t1042 * t939 + t1047 * t936;
t961 = -qJD(4) * t999 + t1008 * t1047 - t1042 * t994;
t987 = pkin(4) * t1005 - qJ(5) * t999;
t997 = t998 ^ 2;
t927 = -pkin(4) * t997 + qJ(5) * t961 - t1005 * t987 + t929;
t979 = -t1035 * t999 + t1038 * t998;
t922 = t1035 * t925 + t1038 * t927 + t979 * t1087;
t980 = t1035 * t998 + t1038 * t999;
t958 = -pkin(5) * t979 - pkin(12) * t980;
t920 = -pkin(5) * t1004 + pkin(12) * t992 + t958 * t979 + t922;
t935 = -t1008 * pkin(3) - t1015 * pkin(11) + t1007 * t996 - t944;
t930 = -t961 * pkin(4) - t997 * qJ(5) + t999 * t987 + qJDD(5) + t935;
t948 = -t1035 * t962 + t1038 * t961;
t949 = t1035 * t961 + t1038 * t962;
t923 = (-t1005 * t979 - t949) * pkin(12) + (t1005 * t980 - t948) * pkin(5) + t930;
t917 = -t1041 * t920 + t1046 * t923;
t963 = t1005 * t1046 - t1041 * t980;
t933 = qJD(6) * t963 + t1041 * t992 + t1046 * t949;
t947 = qJDD(6) - t948;
t964 = t1005 * t1041 + t1046 * t980;
t950 = -mrSges(7,1) * t963 + mrSges(7,2) * t964;
t976 = qJD(6) - t979;
t951 = -mrSges(7,2) * t976 + mrSges(7,3) * t963;
t914 = m(7) * t917 + mrSges(7,1) * t947 - mrSges(7,3) * t933 - t950 * t964 + t951 * t976;
t918 = t1041 * t923 + t1046 * t920;
t932 = -qJD(6) * t964 - t1041 * t949 + t1046 * t992;
t952 = mrSges(7,1) * t976 - mrSges(7,3) * t964;
t915 = m(7) * t918 - mrSges(7,2) * t947 + mrSges(7,3) * t932 + t950 * t963 - t952 * t976;
t906 = -t1041 * t914 + t1046 * t915;
t957 = -mrSges(6,1) * t979 + mrSges(6,2) * t980;
t966 = mrSges(6,1) * t1005 - mrSges(6,3) * t980;
t900 = m(6) * t922 - mrSges(6,2) * t992 + mrSges(6,3) * t948 - t1005 * t966 + t957 * t979 + t906;
t1061 = t1035 * t927 - t1038 * t925;
t919 = -t992 * pkin(5) - t1004 * pkin(12) + (t1087 + t958) * t980 + t1061;
t916 = -m(7) * t919 + t932 * mrSges(7,1) - mrSges(7,2) * t933 + t963 * t951 - t952 * t964;
t921 = -0.2e1 * qJD(5) * t980 - t1061;
t965 = -mrSges(6,2) * t1005 + mrSges(6,3) * t979;
t910 = m(6) * t921 + mrSges(6,1) * t992 - mrSges(6,3) * t949 + t1005 * t965 - t957 * t980 + t916;
t897 = t1035 * t900 + t1038 * t910;
t940 = Ifges(7,5) * t964 + Ifges(7,6) * t963 + Ifges(7,3) * t976;
t942 = Ifges(7,1) * t964 + Ifges(7,4) * t963 + Ifges(7,5) * t976;
t907 = -mrSges(7,1) * t919 + mrSges(7,3) * t918 + Ifges(7,4) * t933 + Ifges(7,2) * t932 + Ifges(7,6) * t947 - t940 * t964 + t942 * t976;
t941 = Ifges(7,4) * t964 + Ifges(7,2) * t963 + Ifges(7,6) * t976;
t908 = mrSges(7,2) * t919 - mrSges(7,3) * t917 + Ifges(7,1) * t933 + Ifges(7,4) * t932 + Ifges(7,5) * t947 + t940 * t963 - t941 * t976;
t954 = Ifges(6,4) * t980 + Ifges(6,2) * t979 + Ifges(6,6) * t1005;
t955 = Ifges(6,1) * t980 + Ifges(6,4) * t979 + Ifges(6,5) * t1005;
t968 = Ifges(5,4) * t999 + Ifges(5,2) * t998 + Ifges(5,6) * t1005;
t969 = Ifges(5,1) * t999 + Ifges(5,4) * t998 + Ifges(5,5) * t1005;
t1088 = Ifges(5,5) * t962 + Ifges(5,6) * t961 + t999 * t968 - t998 * t969 + mrSges(5,1) * t928 - mrSges(5,2) * t929 + Ifges(6,5) * t949 + Ifges(6,6) * t948 + t980 * t954 - t979 * t955 + mrSges(6,1) * t921 - mrSges(6,2) * t922 + t1041 * t908 + t1046 * t907 + pkin(5) * t916 + pkin(12) * t906 + pkin(4) * t897 + (Ifges(5,3) + Ifges(6,3)) * t992;
t1002 = -g(3) * t1074 + t1062;
t1022 = -mrSges(3,2) * t1032 + mrSges(3,3) * t1067;
t1025 = (-mrSges(3,1) * t1049 + mrSges(3,2) * t1044) * t1080;
t1001 = mrSges(4,1) * t1017 - mrSges(4,3) * t1007;
t981 = -mrSges(5,1) * t998 + mrSges(5,2) * t999;
t986 = -mrSges(5,2) * t1005 + mrSges(5,3) * t998;
t895 = m(5) * t928 + mrSges(5,1) * t992 - mrSges(5,3) * t962 + t1005 * t986 - t981 * t999 + t897;
t1065 = -t1035 * t910 + t1038 * t900;
t988 = mrSges(5,1) * t1005 - mrSges(5,3) * t999;
t896 = m(5) * t929 - mrSges(5,2) * t992 + mrSges(5,3) * t961 - t1005 * t988 + t981 * t998 + t1065;
t1064 = -t1042 * t895 + t1047 * t896;
t995 = -mrSges(4,1) * t1006 + mrSges(4,2) * t1007;
t886 = m(4) * t945 - mrSges(4,2) * t1008 + mrSges(4,3) * t993 - t1001 * t1017 + t1006 * t995 + t1064;
t1000 = -mrSges(4,2) * t1017 + mrSges(4,3) * t1006;
t889 = t1042 * t896 + t1047 * t895;
t888 = m(4) * t959 - mrSges(4,1) * t993 + mrSges(4,2) * t994 - t1000 * t1006 + t1001 * t1007 + t889;
t905 = t1041 * t915 + t1046 * t914;
t904 = m(6) * t930 - t948 * mrSges(6,1) + mrSges(6,2) * t949 - t979 * t965 + t966 * t980 + t905;
t1053 = -m(5) * t935 + t961 * mrSges(5,1) - mrSges(5,2) * t962 + t998 * t986 - t988 * t999 - t904;
t903 = m(4) * t944 + mrSges(4,1) * t1008 - mrSges(4,3) * t994 + t1000 * t1017 - t1007 * t995 + t1053;
t875 = t1039 * t1048 * t903 - t1036 * t888 + t886 * t1073;
t871 = m(3) * t1002 + mrSges(3,1) * t1031 - mrSges(3,3) * t1026 + t1022 * t1032 - t1025 * t1068 + t875;
t1012 = -t1037 * t1023 - t1082;
t1021 = mrSges(3,1) * t1032 - mrSges(3,3) * t1068;
t874 = t1039 * t888 + t903 * t1076 + t886 * t1077;
t873 = m(3) * t1012 - t1027 * mrSges(3,1) + t1026 * mrSges(3,2) + (t1021 * t1044 - t1022 * t1049) * t1080 + t874;
t1003 = -g(3) * t1075 + t1069;
t882 = -t1043 * t903 + t1048 * t886;
t881 = m(3) * t1003 - mrSges(3,2) * t1031 + mrSges(3,3) * t1027 - t1021 * t1032 + t1025 * t1067 + t882;
t860 = -t1037 * t873 + t871 * t1070 + t881 * t1071;
t857 = m(2) * t1029 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1051 + t860;
t867 = -t1044 * t871 + t1049 * t881;
t865 = m(2) * t1030 - mrSges(2,1) * t1051 - qJDD(1) * mrSges(2,2) + t867;
t1081 = t1045 * t865 + t1050 * t857;
t859 = t1040 * t873 + t871 * t1074 + t881 * t1075;
t1063 = -t1045 * t857 + t1050 * t865;
t953 = Ifges(6,5) * t980 + Ifges(6,6) * t979 + Ifges(6,3) * t1005;
t890 = mrSges(6,2) * t930 - mrSges(6,3) * t921 + Ifges(6,1) * t949 + Ifges(6,4) * t948 + Ifges(6,5) * t992 - pkin(12) * t905 - t1005 * t954 - t1041 * t907 + t1046 * t908 + t953 * t979;
t1054 = mrSges(7,1) * t917 - mrSges(7,2) * t918 + Ifges(7,5) * t933 + Ifges(7,6) * t932 + Ifges(7,3) * t947 + t941 * t964 - t942 * t963;
t891 = -mrSges(6,1) * t930 + mrSges(6,3) * t922 + Ifges(6,4) * t949 + Ifges(6,2) * t948 + Ifges(6,6) * t992 - pkin(5) * t905 + t1005 * t955 - t953 * t980 - t1054;
t967 = Ifges(5,5) * t999 + Ifges(5,6) * t998 + Ifges(5,3) * t1005;
t876 = -mrSges(5,1) * t935 + mrSges(5,3) * t929 + Ifges(5,4) * t962 + Ifges(5,2) * t961 + Ifges(5,6) * t992 - pkin(4) * t904 + qJ(5) * t1065 + t1005 * t969 + t1035 * t890 + t1038 * t891 - t999 * t967;
t877 = mrSges(5,2) * t935 - mrSges(5,3) * t928 + Ifges(5,1) * t962 + Ifges(5,4) * t961 + Ifges(5,5) * t992 - qJ(5) * t897 - t1005 * t968 - t1035 * t891 + t1038 * t890 + t967 * t998;
t989 = Ifges(4,5) * t1007 + Ifges(4,6) * t1006 + Ifges(4,3) * t1017;
t990 = Ifges(4,4) * t1007 + Ifges(4,2) * t1006 + Ifges(4,6) * t1017;
t862 = mrSges(4,2) * t959 - mrSges(4,3) * t944 + Ifges(4,1) * t994 + Ifges(4,4) * t993 + Ifges(4,5) * t1008 - pkin(11) * t889 + t1006 * t989 - t1017 * t990 - t1042 * t876 + t1047 * t877;
t991 = Ifges(4,1) * t1007 + Ifges(4,4) * t1006 + Ifges(4,5) * t1017;
t868 = -mrSges(4,1) * t959 + mrSges(4,3) * t945 + Ifges(4,4) * t994 + Ifges(4,2) * t993 + Ifges(4,6) * t1008 - pkin(3) * t889 - t1007 * t989 + t1017 * t991 - t1088;
t1056 = pkin(10) * t882 + t1043 * t862 + t1048 * t868;
t1010 = Ifges(3,6) * t1032 + (Ifges(3,4) * t1044 + Ifges(3,2) * t1049) * t1080;
t1011 = Ifges(3,5) * t1032 + (Ifges(3,1) * t1044 + Ifges(3,4) * t1049) * t1080;
t861 = mrSges(4,1) * t944 - mrSges(4,2) * t945 + Ifges(4,5) * t994 + Ifges(4,6) * t993 + Ifges(4,3) * t1008 + pkin(3) * t1053 + pkin(11) * t1064 - t1006 * t991 + t1007 * t990 + t1042 * t877 + t1047 * t876;
t851 = mrSges(3,1) * t1002 - mrSges(3,2) * t1003 + Ifges(3,5) * t1026 + Ifges(3,6) * t1027 + Ifges(3,3) * t1031 + pkin(2) * t875 + t1039 * t861 + (t1010 * t1044 - t1011 * t1049) * t1080 + t1056 * t1036;
t1009 = Ifges(3,3) * t1032 + (Ifges(3,5) * t1044 + Ifges(3,6) * t1049) * t1080;
t853 = -mrSges(3,1) * t1012 + mrSges(3,3) * t1003 + Ifges(3,4) * t1026 + Ifges(3,2) * t1027 + Ifges(3,6) * t1031 - pkin(2) * t874 - t1009 * t1068 + t1032 * t1011 - t1036 * t861 + t1039 * t1056;
t855 = t1009 * t1067 + mrSges(3,2) * t1012 - mrSges(3,3) * t1002 + Ifges(3,1) * t1026 + Ifges(3,4) * t1027 + Ifges(3,5) * t1031 - t1032 * t1010 - t1043 * t868 + t1048 * t862 + (-t1036 * t874 - t1039 * t875) * pkin(10);
t1055 = mrSges(2,1) * t1029 - mrSges(2,2) * t1030 + Ifges(2,3) * qJDD(1) + pkin(1) * t860 + t1040 * t851 + t853 * t1074 + t855 * t1075 + t867 * t1085;
t849 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1029 + Ifges(2,5) * qJDD(1) - t1051 * Ifges(2,6) - t1044 * t853 + t1049 * t855 + (-t1037 * t859 - t1040 * t860) * pkin(9);
t848 = mrSges(2,1) * g(3) + mrSges(2,3) * t1030 + t1051 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t859 - t1037 * t851 + (pkin(9) * t867 + t1044 * t855 + t1049 * t853) * t1040;
t1 = [-m(1) * g(1) + t1063; -m(1) * g(2) + t1081; (-m(1) - m(2)) * g(3) + t859; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1081 - t1045 * t848 + t1050 * t849; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1063 + t1045 * t849 + t1050 * t848; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1055; t1055; t851; t861; t1088; t904; t1054;];
tauJB  = t1;
