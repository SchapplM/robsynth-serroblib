% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:53:38
% EndTime: 2019-05-06 18:53:55
% DurationCPUTime: 8.67s
% Computational Cost: add. (90335->349), mult. (185543->404), div. (0->0), fcn. (114028->8), ass. (0->140)
t1032 = Ifges(6,4) + Ifges(7,4);
t1051 = Ifges(6,2) + Ifges(7,2);
t1043 = Ifges(6,6) + Ifges(7,6);
t991 = cos(qJ(2));
t1020 = qJD(1) * t991;
t1048 = -2 * qJD(3);
t988 = sin(qJ(1));
t992 = cos(qJ(1));
t971 = -t992 * g(1) - t988 * g(2);
t994 = qJD(1) ^ 2;
t945 = -t994 * pkin(1) + qJDD(1) * pkin(7) + t971;
t987 = sin(qJ(2));
t928 = -t987 * g(3) + t991 * t945;
t957 = (-pkin(2) * t991 - qJ(3) * t987) * qJD(1);
t993 = qJD(2) ^ 2;
t904 = t993 * pkin(2) - qJDD(2) * qJ(3) + (qJD(2) * t1048) - t957 * t1020 - t928;
t1019 = qJD(1) * qJD(2);
t1013 = t987 * t1019;
t961 = t991 * qJDD(1) - t1013;
t977 = t987 * qJD(1);
t969 = pkin(3) * t977 - (qJD(2) * pkin(8));
t984 = t991 ^ 2;
t889 = -t984 * t994 * pkin(8) + t961 * pkin(3) + qJD(2) * t969 - t904;
t986 = sin(qJ(4));
t990 = cos(qJ(4));
t956 = t990 * qJD(2) - t986 * t1020;
t919 = -t956 * qJD(4) - t986 * qJDD(2) - t990 * t961;
t974 = t977 + qJD(4);
t929 = t974 * pkin(4) - t956 * pkin(9);
t955 = -t986 * qJD(2) - t990 * t1020;
t953 = t955 ^ 2;
t871 = -t919 * pkin(4) - t953 * pkin(9) + t956 * t929 + t889;
t920 = t955 * qJD(4) + t990 * qJDD(2) - t986 * t961;
t985 = sin(qJ(5));
t989 = cos(qJ(5));
t923 = t985 * t955 + t989 * t956;
t881 = -t923 * qJD(5) + t989 * t919 - t985 * t920;
t972 = qJD(5) + t974;
t908 = t972 * pkin(5) - t923 * qJ(6);
t922 = t989 * t955 - t985 * t956;
t921 = t922 ^ 2;
t861 = -t881 * pkin(5) - t921 * qJ(6) + t923 * t908 + qJDD(6) + t871;
t882 = t922 * qJD(5) + t985 * t919 + t989 * t920;
t909 = t972 * mrSges(7,1) - t923 * mrSges(7,3);
t1015 = m(7) * t861 + t882 * mrSges(7,2) + t923 * t909;
t906 = -t972 * mrSges(7,2) + t922 * mrSges(7,3);
t907 = -t972 * mrSges(6,2) + t922 * mrSges(6,3);
t910 = t972 * mrSges(6,1) - t923 * mrSges(6,3);
t1050 = m(6) * t871 + t882 * mrSges(6,2) + t923 * t910 + t1015 - (t906 + t907) * t922 - (mrSges(6,1) + mrSges(7,1)) * t881;
t925 = -t974 * mrSges(5,2) + t955 * mrSges(5,3);
t926 = t974 * mrSges(5,1) - t956 * mrSges(5,3);
t1049 = -m(5) * t889 + t919 * mrSges(5,1) - t920 * mrSges(5,2) + t955 * t925 - t956 * t926 - t1050;
t1047 = Ifges(3,1) + Ifges(4,2);
t1046 = Ifges(6,1) + Ifges(7,1);
t1033 = Ifges(3,4) + Ifges(4,6);
t1031 = (Ifges(3,5) - Ifges(4,4));
t1045 = Ifges(6,5) + Ifges(7,5);
t1044 = Ifges(3,2) + Ifges(4,3);
t1030 = (Ifges(3,6) - Ifges(4,5));
t1042 = Ifges(3,3) + Ifges(4,1);
t1041 = Ifges(6,3) + Ifges(7,3);
t1040 = t1032 * t923 + t1043 * t972 + t1051 * t922;
t1021 = (t1031 * qJD(2)) + (t1033 * t991 + t1047 * t987) * qJD(1);
t1022 = (t1030 * qJD(2)) + (t1033 * t987 + t1044 * t991) * qJD(1);
t970 = t988 * g(1) - t992 * g(2);
t1007 = -qJDD(1) * pkin(1) - t970;
t1012 = t991 * t1019;
t960 = t987 * qJDD(1) + t1012;
t1000 = pkin(2) * t1013 + t977 * t1048 + (-t960 - t1012) * qJ(3) + t1007;
t885 = -t969 * t977 + (-pkin(3) * t984 - pkin(7)) * t994 + (-pkin(2) - pkin(8)) * t961 + t1000;
t927 = -t991 * g(3) - t987 * t945;
t905 = -qJDD(2) * pkin(2) - t993 * qJ(3) + t957 * t977 + qJDD(3) - t927;
t890 = (-t987 * t991 * t994 - qJDD(2)) * pkin(8) + (t960 - t1012) * pkin(3) + t905;
t868 = -t986 * t885 + t990 * t890;
t954 = qJDD(4) + t960;
t864 = (t955 * t974 - t920) * pkin(9) + (t955 * t956 + t954) * pkin(4) + t868;
t869 = t990 * t885 + t986 * t890;
t866 = -t953 * pkin(4) + t919 * pkin(9) - t974 * t929 + t869;
t858 = t989 * t864 - t985 * t866;
t947 = qJDD(5) + t954;
t853 = -0.2e1 * qJD(6) * t923 + (t922 * t972 - t882) * qJ(6) + (t922 * t923 + t947) * pkin(5) + t858;
t1017 = m(7) * t853 + t947 * mrSges(7,1) + t972 * t906;
t899 = -t922 * mrSges(7,1) + t923 * mrSges(7,2);
t900 = -t922 * mrSges(6,1) + t923 * mrSges(6,2);
t841 = m(6) * t858 + t947 * mrSges(6,1) + t972 * t907 + (-t899 - t900) * t923 + (-mrSges(6,3) - mrSges(7,3)) * t882 + t1017;
t859 = t985 * t864 + t989 * t866;
t855 = -t921 * pkin(5) + t881 * qJ(6) + 0.2e1 * qJD(6) * t922 - t972 * t908 + t859;
t1016 = m(7) * t855 + t881 * mrSges(7,3) + t922 * t899;
t844 = m(6) * t859 + t881 * mrSges(6,3) + t922 * t900 + (-t909 - t910) * t972 + (-mrSges(6,2) - mrSges(7,2)) * t947 + t1016;
t1008 = -t985 * t841 + t989 * t844;
t1025 = -t1032 * t922 - t1045 * t972 - t1046 * t923;
t1026 = -t1041 * t972 - t1043 * t922 - t1045 * t923;
t856 = -t881 * mrSges(7,1) - t922 * t906 + t1015;
t832 = -mrSges(6,1) * t871 + mrSges(6,3) * t859 - mrSges(7,1) * t861 + mrSges(7,3) * t855 - pkin(5) * t856 + qJ(6) * t1016 + (-qJ(6) * t909 - t1025) * t972 + (-qJ(6) * mrSges(7,2) + t1043) * t947 + t1026 * t923 + t1032 * t882 + t1051 * t881;
t850 = -t882 * mrSges(7,3) - t923 * t899 + t1017;
t837 = mrSges(6,2) * t871 + mrSges(7,2) * t861 - mrSges(6,3) * t858 - mrSges(7,3) * t853 - qJ(6) * t850 - t1026 * t922 + t1032 * t881 - t1040 * t972 + t1045 * t947 + t1046 * t882;
t911 = Ifges(5,5) * t956 + Ifges(5,6) * t955 + Ifges(5,3) * t974;
t913 = Ifges(5,1) * t956 + Ifges(5,4) * t955 + Ifges(5,5) * t974;
t817 = -mrSges(5,1) * t889 + mrSges(5,3) * t869 + Ifges(5,4) * t920 + Ifges(5,2) * t919 + Ifges(5,6) * t954 - pkin(4) * t1050 + pkin(9) * t1008 + t989 * t832 + t985 * t837 - t956 * t911 + t974 * t913;
t839 = t989 * t841 + t985 * t844;
t912 = Ifges(5,4) * t956 + Ifges(5,2) * t955 + Ifges(5,6) * t974;
t818 = mrSges(5,2) * t889 - mrSges(5,3) * t868 + Ifges(5,1) * t920 + Ifges(5,4) * t919 + Ifges(5,5) * t954 - pkin(9) * t839 - t985 * t832 + t989 * t837 + t955 * t911 - t974 * t912;
t924 = -t955 * mrSges(5,1) + t956 * mrSges(5,2);
t835 = m(5) * t868 + t954 * mrSges(5,1) - t920 * mrSges(5,3) - t956 * t924 + t974 * t925 + t839;
t836 = m(5) * t869 - t954 * mrSges(5,2) + t919 * mrSges(5,3) + t955 * t924 - t974 * t926 + t1008;
t831 = t990 * t835 + t986 * t836;
t1003 = -m(4) * t905 - t960 * mrSges(4,1) - t831;
t958 = (mrSges(4,2) * t991 - mrSges(4,3) * t987) * qJD(1);
t967 = -mrSges(4,1) * t1020 - qJD(2) * mrSges(4,3);
t830 = qJDD(2) * mrSges(4,2) + qJD(2) * t967 + t958 * t977 - t1003;
t968 = mrSges(4,1) * t977 + qJD(2) * mrSges(4,2);
t996 = -m(4) * t904 + qJDD(2) * mrSges(4,3) + qJD(2) * t968 + t958 * t1020 - t1049;
t1039 = -(t1021 * t991 - t1022 * t987) * qJD(1) + t1042 * qJDD(2) + t1030 * t961 + t1031 * t960 + mrSges(3,1) * t927 - mrSges(3,2) * t928 + mrSges(4,2) * t905 - mrSges(4,3) * t904 - pkin(2) * t830 - pkin(8) * t831 + qJ(3) * (t961 * mrSges(4,1) + t996) - t986 * t817 + t990 * t818;
t1035 = t994 * pkin(7);
t959 = (-mrSges(3,1) * t991 + mrSges(3,2) * t987) * qJD(1);
t966 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1020;
t828 = m(3) * t927 - t960 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t966 - t967) * qJD(2) + (-t958 - t959) * t977 + t1003;
t965 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t977;
t847 = t996 - qJDD(2) * mrSges(3,2) + (mrSges(4,1) + mrSges(3,3)) * t961 - qJD(2) * t965 + t959 * t1020 + m(3) * t928;
t1009 = -t987 * t828 + t991 * t847;
t821 = m(2) * t971 - t994 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t1009;
t1027 = -t986 * t835 + t990 * t836;
t901 = -t961 * pkin(2) + t1000 - t1035;
t1005 = -m(4) * t901 - t961 * mrSges(4,2) + t968 * t977 - t1027;
t944 = t1007 - t1035;
t998 = -m(3) * t944 + t966 * t1020 + t961 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t960 + (-t965 * t987 - t967 * t991) * qJD(1) + t1005;
t825 = m(2) * t970 + qJDD(1) * mrSges(2,1) - t994 * mrSges(2,2) + t998;
t1028 = t988 * t821 + t992 * t825;
t823 = t991 * t828 + t987 * t847;
t1023 = t1042 * qJD(2) + (t1030 * t991 + t1031 * t987) * qJD(1);
t1010 = t992 * t821 - t988 * t825;
t829 = -t960 * mrSges(4,3) + t967 * t1020 - t1005;
t814 = -mrSges(3,1) * t944 - mrSges(4,1) * t904 + mrSges(4,2) * t901 + mrSges(3,3) * t928 - pkin(2) * t829 - pkin(3) * t1049 - pkin(8) * t1027 + t1021 * qJD(2) + t1030 * qJDD(2) - t1023 * t977 + t1033 * t960 + t1044 * t961 - t990 * t817 - t986 * t818;
t999 = mrSges(6,1) * t858 + mrSges(7,1) * t853 - mrSges(6,2) * t859 - mrSges(7,2) * t855 + pkin(5) * t850 + t1025 * t922 + t1040 * t923 + t1041 * t947 + t1043 * t881 + t1045 * t882;
t995 = mrSges(5,1) * t868 - mrSges(5,2) * t869 + Ifges(5,5) * t920 + Ifges(5,6) * t919 + Ifges(5,3) * t954 + pkin(4) * t839 + t956 * t912 - t955 * t913 + t999;
t816 = mrSges(4,1) * t905 + mrSges(3,2) * t944 - mrSges(3,3) * t927 - mrSges(4,3) * t901 + pkin(3) * t831 - qJ(3) * t829 - t1022 * qJD(2) + t1031 * qJDD(2) + t1023 * t1020 + t1033 * t961 + t1047 * t960 + t995;
t1002 = mrSges(2,1) * t970 - mrSges(2,2) * t971 + Ifges(2,3) * qJDD(1) + pkin(1) * t998 + pkin(7) * t1009 + t991 * t814 + t987 * t816;
t812 = mrSges(2,1) * g(3) + mrSges(2,3) * t971 + t994 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t823 - t1039;
t811 = -mrSges(2,2) * g(3) - mrSges(2,3) * t970 + Ifges(2,5) * qJDD(1) - t994 * Ifges(2,6) - pkin(7) * t823 - t987 * t814 + t991 * t816;
t1 = [-m(1) * g(1) + t1010; -m(1) * g(2) + t1028; (-m(1) - m(2)) * g(3) + t823; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1028 + t992 * t811 - t988 * t812; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1010 + t988 * t811 + t992 * t812; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1002; t1002; t1039; t830; t995; t999; t856;];
tauJB  = t1;
