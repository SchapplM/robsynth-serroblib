% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-05-06 09:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:09:37
% EndTime: 2019-05-06 09:09:50
% DurationCPUTime: 7.61s
% Computational Cost: add. (71945->348), mult. (165847->404), div. (0->0), fcn. (109055->8), ass. (0->140)
t1031 = Ifges(6,4) + Ifges(7,4);
t1046 = Ifges(6,2) + Ifges(7,2);
t1039 = Ifges(6,6) + Ifges(7,6);
t1045 = -2 * qJD(4);
t1044 = Ifges(4,1) + Ifges(5,2);
t1043 = -Ifges(5,1) - Ifges(4,3);
t1042 = Ifges(6,1) + Ifges(7,1);
t1030 = Ifges(4,5) - Ifges(5,4);
t1041 = Ifges(6,5) + Ifges(7,5);
t1040 = -Ifges(4,2) - Ifges(5,3);
t1029 = Ifges(4,6) - Ifges(5,5);
t1028 = -Ifges(5,6) - Ifges(4,4);
t1038 = Ifges(6,3) + Ifges(7,3);
t987 = cos(qJ(2));
t1015 = qJD(1) * t987;
t984 = sin(qJ(2));
t1016 = qJD(1) * t984;
t1027 = cos(pkin(9));
t982 = sin(pkin(9));
t955 = -t1027 * t1015 + t982 * t1016;
t983 = sin(qJ(5));
t986 = cos(qJ(5));
t937 = -qJD(2) * t983 + t955 * t986;
t938 = qJD(2) * t986 + t955 * t983;
t956 = (t1027 * t984 + t982 * t987) * qJD(1);
t953 = qJD(5) + t956;
t1037 = t1031 * t938 + t1039 * t953 + t1046 * t937;
t1011 = qJD(1) * qJD(2);
t985 = sin(qJ(1));
t988 = cos(qJ(1));
t974 = -g(1) * t988 - g(2) * t985;
t990 = qJD(1) ^ 2;
t962 = -pkin(1) * t990 + qJDD(1) * pkin(7) + t974;
t1026 = t962 * t984;
t1033 = pkin(2) * t990;
t967 = qJDD(1) * t984 + t987 * t1011;
t901 = qJDD(2) * pkin(2) - qJ(3) * t967 - t1026 + (qJ(3) * t1011 + t984 * t1033 - g(3)) * t987;
t941 = -g(3) * t984 + t987 * t962;
t968 = qJDD(1) * t987 - t984 * t1011;
t970 = qJD(2) * pkin(2) - qJ(3) * t1016;
t981 = t987 ^ 2;
t902 = qJ(3) * t968 - qJD(2) * t970 - t981 * t1033 + t941;
t1021 = t1027 * t902 + t982 * t901;
t923 = pkin(3) * t955 - qJ(4) * t956;
t989 = qJD(2) ^ 2;
t1036 = pkin(3) * t989 - qJDD(2) * qJ(4) + qJD(2) * t1045 + t955 * t923 - t1021;
t875 = -0.2e1 * qJD(3) * t956 + t1027 * t901 - t982 * t902;
t934 = -t1027 * t968 + t967 * t982;
t947 = pkin(4) * t956 - qJD(2) * pkin(8);
t1013 = qJD(3) * t955;
t950 = -0.2e1 * t1013;
t954 = t955 ^ 2;
t868 = -pkin(4) * t934 - pkin(8) * t954 + qJD(2) * t947 - t1036 + t950;
t895 = -qJD(5) * t938 - qJDD(2) * t983 + t934 * t986;
t911 = pkin(5) * t953 - qJ(6) * t938;
t936 = t937 ^ 2;
t863 = -t895 * pkin(5) - qJ(6) * t936 + t911 * t938 + qJDD(6) + t868;
t896 = qJD(5) * t937 + qJDD(2) * t986 + t934 * t983;
t912 = mrSges(7,1) * t953 - mrSges(7,3) * t938;
t1006 = m(7) * t863 + t896 * mrSges(7,2) + t938 * t912;
t909 = -mrSges(7,2) * t953 + mrSges(7,3) * t937;
t910 = -mrSges(6,2) * t953 + mrSges(6,3) * t937;
t913 = mrSges(6,1) * t953 - mrSges(6,3) * t938;
t1035 = -m(6) * t868 - t896 * mrSges(6,2) + (t909 + t910) * t937 + (mrSges(6,1) + mrSges(7,1)) * t895 - t938 * t913 - t1006;
t1017 = t1030 * qJD(2) + t1028 * t955 + t1044 * t956;
t1018 = t1029 * qJD(2) - t1028 * t956 + t1040 * t955;
t924 = mrSges(4,1) * t955 + mrSges(4,2) * t956;
t935 = t1027 * t967 + t982 * t968;
t943 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t955;
t945 = mrSges(5,1) * t955 - qJD(2) * mrSges(5,3);
t1014 = qJD(2) * t955;
t872 = -qJDD(2) * pkin(3) - t989 * qJ(4) + t956 * t923 + qJDD(4) - t875;
t866 = (t955 * t956 - qJDD(2)) * pkin(8) + (t935 + t1014) * pkin(4) + t872;
t973 = g(1) * t985 - t988 * g(2);
t1000 = -qJDD(1) * pkin(1) - t973;
t908 = -pkin(2) * t968 + qJDD(3) + t970 * t1016 + (-qJ(3) * t981 - pkin(7)) * t990 + t1000;
t993 = (-t935 + t1014) * qJ(4) + t908 + (qJD(2) * pkin(3) + t1045) * t956;
t870 = -pkin(4) * t954 - t947 * t956 + (pkin(3) + pkin(8)) * t934 + t993;
t860 = t986 * t866 - t870 * t983;
t933 = qJDD(5) + t935;
t855 = -0.2e1 * qJD(6) * t938 + (t937 * t953 - t896) * qJ(6) + (t937 * t938 + t933) * pkin(5) + t860;
t1008 = m(7) * t855 + t933 * mrSges(7,1) + t953 * t909;
t903 = -mrSges(7,1) * t937 + mrSges(7,2) * t938;
t904 = -mrSges(6,1) * t937 + mrSges(6,2) * t938;
t847 = m(6) * t860 + t933 * mrSges(6,1) + t910 * t953 + (-t903 - t904) * t938 + (-mrSges(6,3) - mrSges(7,3)) * t896 + t1008;
t861 = t983 * t866 + t986 * t870;
t857 = -pkin(5) * t936 + t895 * qJ(6) + 0.2e1 * qJD(6) * t937 - t911 * t953 + t861;
t1007 = m(7) * t857 + t895 * mrSges(7,3) + t937 * t903;
t849 = m(6) * t861 + t895 * mrSges(6,3) + t904 * t937 + (-t912 - t913) * t953 + (-mrSges(6,2) - mrSges(7,2)) * t933 + t1007;
t842 = t847 * t986 + t849 * t983;
t925 = -mrSges(5,2) * t955 - mrSges(5,3) * t956;
t996 = -m(5) * t872 - t935 * mrSges(5,1) - t956 * t925 - t842;
t837 = m(4) * t875 - mrSges(4,3) * t935 - t924 * t956 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t943 - t945) * qJD(2) + t996;
t876 = t950 + t1021;
t944 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t956;
t871 = 0.2e1 * t1013 + t1036;
t946 = mrSges(5,1) * t956 + qJD(2) * mrSges(5,2);
t995 = -m(5) * t871 + qJDD(2) * mrSges(5,3) + qJD(2) * t946 - t1035;
t845 = -qJDD(2) * mrSges(4,2) + t995 + (-t924 - t925) * t955 + (-mrSges(4,3) - mrSges(5,1)) * t934 - qJD(2) * t944 + m(4) * t876;
t831 = t1027 * t837 + t982 * t845;
t1022 = -t1031 * t937 - t1041 * t953 - t1042 * t938;
t1023 = -t1038 * t953 - t1039 * t937 - t1041 * t938;
t858 = -t895 * mrSges(7,1) - t909 * t937 + t1006;
t832 = -mrSges(6,1) * t868 + mrSges(6,3) * t861 - mrSges(7,1) * t863 + mrSges(7,3) * t857 - pkin(5) * t858 + qJ(6) * t1007 + (-qJ(6) * t912 - t1022) * t953 + t1023 * t938 + (-mrSges(7,2) * qJ(6) + t1039) * t933 + t1031 * t896 + t1046 * t895;
t840 = qJDD(2) * mrSges(5,2) + qJD(2) * t945 - t996;
t852 = -t896 * mrSges(7,3) - t903 * t938 + t1008;
t841 = mrSges(6,2) * t868 + mrSges(7,2) * t863 - mrSges(6,3) * t860 - mrSges(7,3) * t855 - qJ(6) * t852 - t1023 * t937 + t1031 * t895 - t1037 * t953 + t1041 * t933 + t1042 * t896;
t940 = -g(3) * t987 - t1026;
t958 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t984 + Ifges(3,2) * t987) * qJD(1);
t959 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t984 + Ifges(3,4) * t987) * qJD(1);
t1034 = (t984 * t958 - t987 * t959) * qJD(1) + (Ifges(3,3) - t1043) * qJDD(2) + t1017 * t955 + t1018 * t956 - t1029 * t934 + t1030 * t935 + mrSges(3,1) * t940 + mrSges(4,1) * t875 - mrSges(3,2) * t941 - mrSges(4,2) * t876 + mrSges(5,2) * t872 - mrSges(5,3) * t871 + Ifges(3,5) * t967 + Ifges(3,6) * t968 + pkin(2) * t831 - pkin(3) * t840 - pkin(8) * t842 + qJ(4) * (-mrSges(5,1) * t934 - t925 * t955 + t995) - t983 * t832 + t986 * t841;
t966 = (-mrSges(3,1) * t987 + mrSges(3,2) * t984) * qJD(1);
t972 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1015;
t829 = m(3) * t940 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t967 + qJD(2) * t972 - t966 * t1016 + t831;
t1003 = t1027 * t845 - t837 * t982;
t971 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1016;
t830 = m(3) * t941 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t968 - qJD(2) * t971 + t966 * t1015 + t1003;
t1004 = -t829 * t984 + t987 * t830;
t823 = m(2) * t974 - mrSges(2,1) * t990 - qJDD(1) * mrSges(2,2) + t1004;
t1024 = -t983 * t847 + t986 * t849;
t874 = pkin(3) * t934 + t993;
t839 = m(5) * t874 - t934 * mrSges(5,2) - t935 * mrSges(5,3) - t955 * t945 - t956 * t946 + t1024;
t838 = m(4) * t908 + t934 * mrSges(4,1) + mrSges(4,2) * t935 + t955 * t943 + t944 * t956 + t839;
t961 = -pkin(7) * t990 + t1000;
t992 = -m(3) * t961 + t968 * mrSges(3,1) - mrSges(3,2) * t967 + t972 * t1015 - t971 * t1016 - t838;
t834 = m(2) * t973 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t990 + t992;
t1025 = t985 * t823 + t988 * t834;
t825 = t987 * t829 + t984 * t830;
t1019 = t1043 * qJD(2) + t1029 * t955 - t1030 * t956;
t1005 = t988 * t823 - t834 * t985;
t819 = -mrSges(4,1) * t908 - mrSges(5,1) * t871 + mrSges(5,2) * t874 + mrSges(4,3) * t876 - pkin(3) * t839 - pkin(4) * t1035 - pkin(8) * t1024 + t1017 * qJD(2) + t1029 * qJDD(2) + t1019 * t956 - t1028 * t935 + t1040 * t934 - t986 * t832 - t983 * t841;
t994 = mrSges(6,1) * t860 + mrSges(7,1) * t855 - mrSges(6,2) * t861 - mrSges(7,2) * t857 + pkin(5) * t852 + t1022 * t937 + t1037 * t938 + t1038 * t933 + t1039 * t895 + t1041 * t896;
t820 = mrSges(5,1) * t872 + mrSges(4,2) * t908 - mrSges(4,3) * t875 - mrSges(5,3) * t874 + pkin(4) * t842 - qJ(4) * t839 - t1018 * qJD(2) + t1030 * qJDD(2) + t1019 * t955 + t1028 * t934 + t1044 * t935 + t994;
t957 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t984 + Ifges(3,6) * t987) * qJD(1);
t815 = -mrSges(3,1) * t961 + mrSges(3,3) * t941 + Ifges(3,4) * t967 + Ifges(3,2) * t968 + Ifges(3,6) * qJDD(2) - pkin(2) * t838 + qJ(3) * t1003 + qJD(2) * t959 - t957 * t1016 + t1027 * t819 + t982 * t820;
t817 = mrSges(3,2) * t961 - mrSges(3,3) * t940 + Ifges(3,1) * t967 + Ifges(3,4) * t968 + Ifges(3,5) * qJDD(2) - qJ(3) * t831 - qJD(2) * t958 + t957 * t1015 + t1027 * t820 - t982 * t819;
t997 = mrSges(2,1) * t973 - mrSges(2,2) * t974 + Ifges(2,3) * qJDD(1) + pkin(1) * t992 + pkin(7) * t1004 + t987 * t815 + t984 * t817;
t818 = mrSges(2,1) * g(3) + mrSges(2,3) * t974 + t990 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t825 - t1034;
t813 = -mrSges(2,2) * g(3) - mrSges(2,3) * t973 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t990 - pkin(7) * t825 - t815 * t984 + t817 * t987;
t1 = [-m(1) * g(1) + t1005; -m(1) * g(2) + t1025; (-m(1) - m(2)) * g(3) + t825; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1025 + t988 * t813 - t985 * t818; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1005 + t985 * t813 + t988 * t818; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t997; t997; t1034; t838; t840; t994; t858;];
tauJB  = t1;
