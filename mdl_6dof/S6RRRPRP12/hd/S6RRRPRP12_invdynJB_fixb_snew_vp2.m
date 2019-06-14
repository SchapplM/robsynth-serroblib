% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:24:55
% EndTime: 2019-05-07 09:25:21
% DurationCPUTime: 12.55s
% Computational Cost: add. (182714->358), mult. (391068->433), div. (0->0), fcn. (296817->10), ass. (0->150)
t1067 = Ifges(6,1) + Ifges(7,1);
t1053 = Ifges(6,4) - Ifges(7,5);
t1064 = Ifges(7,4) + Ifges(6,5);
t1066 = Ifges(6,2) + Ifges(7,3);
t1062 = Ifges(6,6) - Ifges(7,6);
t1065 = Ifges(4,1) + Ifges(5,2);
t1052 = -Ifges(4,5) + Ifges(5,4);
t1063 = -Ifges(4,2) - Ifges(5,3);
t1051 = Ifges(4,6) - Ifges(5,5);
t1050 = -Ifges(5,6) - Ifges(4,4);
t1061 = Ifges(4,3) + Ifges(5,1);
t1060 = Ifges(6,3) + Ifges(7,2);
t1007 = sin(qJ(5));
t1055 = cos(qJ(5));
t1006 = cos(pkin(6));
t1001 = qJD(1) * t1006 + qJD(2);
t1008 = sin(qJ(3));
t1005 = sin(pkin(6));
t1009 = sin(qJ(2));
t1034 = t1005 * t1009;
t1029 = qJD(1) * t1034;
t1056 = cos(qJ(3));
t975 = -t1001 * t1056 + t1008 * t1029;
t1011 = cos(qJ(2));
t1033 = t1005 * t1011;
t1028 = qJD(1) * t1033;
t995 = -qJD(3) + t1028;
t956 = -t1007 * t995 - t1055 * t975;
t957 = t1007 * t975 - t1055 * t995;
t976 = t1008 * t1001 + t1029 * t1056;
t973 = qJD(5) + t976;
t1059 = -t1053 * t957 - t1062 * t973 + t1066 * t956;
t1058 = -t1053 * t956 + t1064 * t973 + t1067 * t957;
t989 = (qJD(1) * qJD(2) * t1009 - qJDD(1) * t1011) * t1005;
t1000 = qJDD(1) * t1006 + qJDD(2);
t1035 = qJD(1) * t1011;
t1032 = t1006 * t1009;
t1013 = qJD(1) ^ 2;
t1048 = pkin(8) * t1005;
t1010 = sin(qJ(1));
t1012 = cos(qJ(1));
t997 = t1010 * g(1) - g(2) * t1012;
t984 = qJDD(1) * pkin(1) + t1013 * t1048 + t997;
t998 = -g(1) * t1012 - g(2) * t1010;
t985 = -pkin(1) * t1013 + qJDD(1) * t1048 + t998;
t1038 = t1011 * t985 + t984 * t1032;
t1036 = qJD(1) * t1005;
t987 = (-pkin(2) * t1011 - pkin(9) * t1009) * t1036;
t999 = t1001 ^ 2;
t916 = -t999 * pkin(2) + t1000 * pkin(9) + (-g(3) * t1009 + t1035 * t987) * t1005 + t1038;
t1047 = t1006 * g(3);
t988 = (qJD(2) * t1035 + qJDD(1) * t1009) * t1005;
t917 = t989 * pkin(2) - t988 * pkin(9) - t1047 + (-t984 + (pkin(2) * t1009 - pkin(9) * t1011) * t1001 * qJD(1)) * t1005;
t893 = t1008 * t917 + t1056 * t916;
t951 = pkin(3) * t975 - qJ(4) * t976;
t981 = qJDD(3) + t989;
t994 = t995 ^ 2;
t889 = t994 * pkin(3) - t981 * qJ(4) + 0.2e1 * qJD(4) * t995 + t975 * t951 - t893;
t947 = t976 * qJD(3) - t1000 * t1056 + t1008 * t988;
t962 = pkin(4) * t976 + pkin(10) * t995;
t974 = t975 ^ 2;
t887 = -t947 * pkin(4) - t974 * pkin(10) - t995 * t962 - t889;
t902 = t957 * qJD(5) + t1007 * t981 - t1055 * t947;
t903 = -t956 * qJD(5) + t1007 * t947 + t1055 * t981;
t882 = -0.2e1 * qJD(6) * t957 + (t956 * t973 - t903) * qJ(6) + (t957 * t973 + t902) * pkin(5) + t887;
t926 = -mrSges(7,2) * t956 + mrSges(7,3) * t973;
t929 = -mrSges(7,1) * t973 + mrSges(7,2) * t957;
t873 = m(7) * t882 + t902 * mrSges(7,1) - t903 * mrSges(7,3) + t956 * t926 - t957 * t929;
t927 = -mrSges(6,2) * t973 - mrSges(6,3) * t956;
t928 = mrSges(6,1) * t973 - mrSges(6,3) * t957;
t1019 = m(6) * t887 + t902 * mrSges(6,1) + t903 * mrSges(6,2) + t956 * t927 + t957 * t928 + t873;
t961 = mrSges(5,1) * t976 - mrSges(5,2) * t995;
t1016 = -m(5) * t889 + t981 * mrSges(5,3) - t995 * t961 + t1019;
t1039 = t1050 * t975 + t1052 * t995 + t1065 * t976;
t1040 = -t1050 * t976 - t1051 * t995 + t1063 * t975;
t1046 = t975 * t995;
t892 = -t1008 * t916 + t1056 * t917;
t890 = -t981 * pkin(3) - t994 * qJ(4) + t976 * t951 + qJDD(4) - t892;
t948 = -t975 * qJD(3) + t1008 * t1000 + t1056 * t988;
t884 = (t975 * t976 - t981) * pkin(10) + (t948 - t1046) * pkin(4) + t890;
t1031 = t1006 * t1011;
t949 = -g(3) * t1033 - t1009 * t985 + t1031 * t984;
t915 = -t1000 * pkin(2) - t999 * pkin(9) + t987 * t1029 - t949;
t1017 = (-t948 - t1046) * qJ(4) + t915 + (-pkin(3) * t995 - 0.2e1 * qJD(4)) * t976;
t888 = -t974 * pkin(4) - t976 * t962 + (pkin(3) + pkin(10)) * t947 + t1017;
t880 = t1007 * t884 + t1055 * t888;
t920 = pkin(5) * t956 - qJ(6) * t957;
t944 = qJDD(5) + t948;
t972 = t973 ^ 2;
t876 = -pkin(5) * t972 + qJ(6) * t944 + 0.2e1 * qJD(6) * t973 - t920 * t956 + t880;
t1030 = m(7) * t876 + t944 * mrSges(7,3) + t973 * t929;
t921 = mrSges(7,1) * t956 - mrSges(7,3) * t957;
t1042 = -mrSges(6,1) * t956 - mrSges(6,2) * t957 - t921;
t1054 = -mrSges(6,3) - mrSges(7,2);
t866 = m(6) * t880 - t944 * mrSges(6,2) + t1042 * t956 + t1054 * t902 - t973 * t928 + t1030;
t879 = -t1007 * t888 + t1055 * t884;
t877 = -t944 * pkin(5) - t972 * qJ(6) + t957 * t920 + qJDD(6) - t879;
t1025 = -m(7) * t877 + t944 * mrSges(7,1) + t973 * t926;
t868 = m(6) * t879 + t944 * mrSges(6,1) + t1042 * t957 + t1054 * t903 + t973 * t927 + t1025;
t861 = t1007 * t866 + t1055 * t868;
t953 = -mrSges(5,2) * t975 - mrSges(5,3) * t976;
t1020 = m(5) * t890 + t948 * mrSges(5,1) + t976 * t953 + t861;
t960 = mrSges(5,1) * t975 + mrSges(5,3) * t995;
t858 = t981 * mrSges(5,2) - t995 * t960 + t1020;
t1043 = -t1060 * t973 + t1062 * t956 - t1064 * t957;
t859 = -mrSges(6,1) * t887 - mrSges(7,1) * t882 + mrSges(7,2) * t876 + mrSges(6,3) * t880 - pkin(5) * t873 + t1043 * t957 + t1053 * t903 + t1058 * t973 + t1062 * t944 - t1066 * t902;
t860 = mrSges(6,2) * t887 + mrSges(7,2) * t877 - mrSges(6,3) * t879 - mrSges(7,3) * t882 - qJ(6) * t873 + t1043 * t956 - t1053 * t902 + t1059 * t973 + t1064 * t944 + t1067 * t903;
t1057 = t1039 * t975 + t1040 * t976 + t1061 * t981 - t1051 * t947 - t1052 * t948 + mrSges(4,1) * t892 - mrSges(4,2) * t893 + mrSges(5,2) * t890 - mrSges(5,3) * t889 - pkin(3) * t858 - pkin(10) * t861 + qJ(4) * (-t947 * mrSges(5,1) - t975 * t953 + t1016) - t1007 * t859 + t1055 * t860;
t952 = mrSges(4,1) * t975 + mrSges(4,2) * t976;
t958 = mrSges(4,2) * t995 - mrSges(4,3) * t975;
t856 = m(4) * t892 - t948 * mrSges(4,3) - t976 * t952 + (-t958 + t960) * t995 + (mrSges(4,1) - mrSges(5,2)) * t981 - t1020;
t959 = -mrSges(4,1) * t995 - mrSges(4,3) * t976;
t864 = t1016 - t981 * mrSges(4,2) + t995 * t959 + (-t952 - t953) * t975 + (-mrSges(4,3) - mrSges(5,1)) * t947 + m(4) * t893;
t1027 = -t1008 * t856 + t1056 * t864;
t950 = -g(3) * t1034 + t1038;
t982 = mrSges(3,1) * t1001 - mrSges(3,3) * t1029;
t986 = (-mrSges(3,1) * t1011 + mrSges(3,2) * t1009) * t1036;
t848 = m(3) * t950 - mrSges(3,2) * t1000 - mrSges(3,3) * t989 - t1001 * t982 + t1028 * t986 + t1027;
t851 = t1008 * t864 + t1056 * t856;
t967 = -t1005 * t984 - t1047;
t983 = -mrSges(3,2) * t1001 + mrSges(3,3) * t1028;
t850 = m(3) * t967 + t989 * mrSges(3,1) + t988 * mrSges(3,2) + (t1009 * t982 - t1011 * t983) * t1036 + t851;
t1044 = -t1007 * t868 + t1055 * t866;
t891 = t947 * pkin(3) + t1017;
t1024 = -m(5) * t891 + t947 * mrSges(5,2) + t975 * t960 - t1044;
t1018 = -m(4) * t915 - t947 * mrSges(4,1) - t975 * t958 + (-t959 + t961) * t976 + (-mrSges(4,2) + mrSges(5,3)) * t948 + t1024;
t854 = m(3) * t949 + t1000 * mrSges(3,1) - t988 * mrSges(3,3) + t1001 * t983 - t1029 * t986 + t1018;
t838 = -t1005 * t850 + t854 * t1031 + t848 * t1032;
t835 = m(2) * t997 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1013 + t838;
t844 = -t1009 * t854 + t1011 * t848;
t842 = m(2) * t998 - mrSges(2,1) * t1013 - qJDD(1) * mrSges(2,2) + t844;
t1045 = t1010 * t842 + t1012 * t835;
t1041 = t1051 * t975 + t1052 * t976 + t1061 * t995;
t837 = t1006 * t850 + t854 * t1033 + t848 * t1034;
t1026 = -t1010 * t835 + t1012 * t842;
t857 = -t948 * mrSges(5,3) - t976 * t961 - t1024;
t833 = -mrSges(4,1) * t915 - mrSges(5,1) * t889 + mrSges(5,2) * t891 + mrSges(4,3) * t893 - pkin(3) * t857 + pkin(4) * t1019 - pkin(10) * t1044 - t1007 * t860 - t1039 * t995 + t1041 * t976 - t1050 * t948 + t1051 * t981 - t1055 * t859 + t1063 * t947;
t872 = t903 * mrSges(7,2) + t957 * t921 - t1025;
t1015 = mrSges(6,1) * t879 - mrSges(7,1) * t877 - mrSges(6,2) * t880 + mrSges(7,3) * t876 - pkin(5) * t872 + qJ(6) * t1030 - t1059 * t957 + (-qJ(6) * t921 + t1058) * t956 + t1060 * t944 + t1064 * t903 + (-mrSges(7,2) * qJ(6) - t1062) * t902;
t839 = mrSges(5,1) * t890 + mrSges(4,2) * t915 - mrSges(4,3) * t892 - mrSges(5,3) * t891 + pkin(4) * t861 - qJ(4) * t857 + t1040 * t995 + t1041 * t975 + t1050 * t947 - t1052 * t981 + t1065 * t948 + t1015;
t965 = Ifges(3,6) * t1001 + (Ifges(3,4) * t1009 + Ifges(3,2) * t1011) * t1036;
t966 = Ifges(3,5) * t1001 + (Ifges(3,1) * t1009 + Ifges(3,4) * t1011) * t1036;
t828 = Ifges(3,5) * t988 - Ifges(3,6) * t989 + Ifges(3,3) * t1000 + mrSges(3,1) * t949 - mrSges(3,2) * t950 + t1008 * t839 + t1056 * t833 + pkin(2) * t1018 + pkin(9) * t1027 + (t1009 * t965 - t1011 * t966) * t1036;
t964 = Ifges(3,3) * t1001 + (Ifges(3,5) * t1009 + Ifges(3,6) * t1011) * t1036;
t830 = mrSges(3,2) * t967 - mrSges(3,3) * t949 + Ifges(3,1) * t988 - Ifges(3,4) * t989 + Ifges(3,5) * t1000 - pkin(9) * t851 - t1001 * t965 - t1008 * t833 + t964 * t1028 + t1056 * t839;
t832 = -mrSges(3,1) * t967 + mrSges(3,3) * t950 + Ifges(3,4) * t988 - Ifges(3,2) * t989 + Ifges(3,6) * t1000 - pkin(2) * t851 + t1001 * t966 - t1029 * t964 - t1057;
t1021 = mrSges(2,1) * t997 - mrSges(2,2) * t998 + Ifges(2,3) * qJDD(1) + pkin(1) * t838 + t1006 * t828 + t832 * t1033 + t830 * t1034 + t844 * t1048;
t826 = -mrSges(2,2) * g(3) - mrSges(2,3) * t997 + Ifges(2,5) * qJDD(1) - t1013 * Ifges(2,6) - t1009 * t832 + t1011 * t830 + (-t1005 * t837 - t1006 * t838) * pkin(8);
t825 = mrSges(2,1) * g(3) + mrSges(2,3) * t998 + t1013 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t837 - t1005 * t828 + (pkin(8) * t844 + t1009 * t830 + t1011 * t832) * t1006;
t1 = [-m(1) * g(1) + t1026; -m(1) * g(2) + t1045; (-m(1) - m(2)) * g(3) + t837; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1045 - t1010 * t825 + t1012 * t826; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1026 + t1010 * t826 + t1012 * t825; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1021; t1021; t828; t1057; t858; t1015; t872;];
tauJB  = t1;
