% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:04:23
% EndTime: 2019-05-06 23:05:41
% DurationCPUTime: 74.99s
% Computational Cost: add. (1223263->401), mult. (2784719->518), div. (0->0), fcn. (2289384->14), ass. (0->166)
t992 = sin(pkin(6));
t1030 = pkin(8) * t992;
t994 = cos(pkin(6));
t1029 = t994 * g(3);
t998 = sin(qJ(2));
t1028 = t992 * t998;
t1027 = t994 * t998;
t1004 = cos(qJ(1));
t1005 = qJD(1) ^ 2;
t1003 = cos(qJ(2));
t1022 = t1003 * t994;
t1021 = qJD(1) * t1003;
t1018 = t992 * t1021;
t1002 = cos(qJ(4));
t1001 = cos(qJ(5));
t1000 = cos(qJ(6));
t999 = sin(qJ(1));
t982 = t999 * g(1) - g(2) * t1004;
t972 = qJDD(1) * pkin(1) + t1005 * t1030 + t982;
t1020 = qJDD(1) * t992;
t983 = -g(1) * t1004 - g(2) * t999;
t973 = -pkin(1) * t1005 + pkin(8) * t1020 + t983;
t1025 = t1003 * t973 + t972 * t1027;
t1024 = qJD(1) * t992;
t974 = (-pkin(2) * t1003 - qJ(3) * t998) * t1024;
t987 = qJD(1) * t994 + qJD(2);
t985 = t987 ^ 2;
t986 = qJDD(1) * t994 + qJDD(2);
t932 = -t985 * pkin(2) + t986 * qJ(3) + (-g(3) * t998 + t1021 * t974) * t992 + t1025;
t976 = (qJD(2) * t1021 + qJDD(1) * t998) * t992;
t1019 = t998 * t1024;
t977 = -qJD(2) * t1019 + t1003 * t1020;
t933 = -t977 * pkin(2) - t1029 - t976 * qJ(3) + (-t972 + (pkin(2) * t998 - qJ(3) * t1003) * t987 * qJD(1)) * t992;
t991 = sin(pkin(12));
t993 = cos(pkin(12));
t965 = t1019 * t993 + t987 * t991;
t907 = -0.2e1 * qJD(3) * t965 - t991 * t932 + t993 * t933;
t952 = t976 * t993 + t986 * t991;
t964 = -t1019 * t991 + t987 * t993;
t895 = (-t1018 * t964 - t952) * pkin(9) + (t964 * t965 - t977) * pkin(3) + t907;
t908 = 0.2e1 * qJD(3) * t964 + t993 * t932 + t991 * t933;
t951 = -t976 * t991 + t986 * t993;
t953 = -pkin(3) * t1018 - pkin(9) * t965;
t962 = t964 ^ 2;
t898 = -pkin(3) * t962 + pkin(9) * t951 + t1018 * t953 + t908;
t997 = sin(qJ(4));
t875 = t1002 * t895 - t898 * t997;
t945 = t1002 * t964 - t965 * t997;
t916 = qJD(4) * t945 + t1002 * t952 + t951 * t997;
t946 = t1002 * t965 + t964 * t997;
t969 = qJDD(4) - t977;
t981 = qJD(4) - t1018;
t872 = (t945 * t981 - t916) * pkin(10) + (t945 * t946 + t969) * pkin(4) + t875;
t876 = t1002 * t898 + t997 * t895;
t915 = -qJD(4) * t946 + t1002 * t951 - t952 * t997;
t936 = pkin(4) * t981 - pkin(10) * t946;
t944 = t945 ^ 2;
t874 = -pkin(4) * t944 + pkin(10) * t915 - t936 * t981 + t876;
t996 = sin(qJ(5));
t869 = t1001 * t874 + t996 * t872;
t925 = t1001 * t945 - t946 * t996;
t926 = t1001 * t946 + t945 * t996;
t910 = -pkin(5) * t925 - pkin(11) * t926;
t968 = qJDD(5) + t969;
t979 = qJD(5) + t981;
t978 = t979 ^ 2;
t866 = -pkin(5) * t978 + pkin(11) * t968 + t910 * t925 + t869;
t1023 = t1003 * t992;
t942 = -g(3) * t1023 + t1022 * t972 - t998 * t973;
t931 = -t986 * pkin(2) - t985 * qJ(3) + t974 * t1019 + qJDD(3) - t942;
t911 = -t951 * pkin(3) - t962 * pkin(9) + t965 * t953 + t931;
t881 = -t915 * pkin(4) - t944 * pkin(10) + t946 * t936 + t911;
t891 = -qJD(5) * t926 + t1001 * t915 - t916 * t996;
t892 = qJD(5) * t925 + t1001 * t916 + t915 * t996;
t870 = t881 + (-t925 * t979 - t892) * pkin(11) + (t926 * t979 - t891) * pkin(5);
t995 = sin(qJ(6));
t863 = t1000 * t870 - t866 * t995;
t912 = t1000 * t979 - t926 * t995;
t879 = qJD(6) * t912 + t1000 * t892 + t968 * t995;
t886 = qJDD(6) - t891;
t913 = t1000 * t926 + t979 * t995;
t899 = -mrSges(7,1) * t912 + mrSges(7,2) * t913;
t924 = qJD(6) - t925;
t900 = -mrSges(7,2) * t924 + mrSges(7,3) * t912;
t859 = m(7) * t863 + mrSges(7,1) * t886 - mrSges(7,3) * t879 - t899 * t913 + t900 * t924;
t864 = t1000 * t866 + t870 * t995;
t878 = -qJD(6) * t913 + t1000 * t968 - t892 * t995;
t901 = mrSges(7,1) * t924 - mrSges(7,3) * t913;
t860 = m(7) * t864 - mrSges(7,2) * t886 + mrSges(7,3) * t878 + t899 * t912 - t901 * t924;
t1013 = t1000 * t860 - t859 * t995;
t909 = -mrSges(6,1) * t925 + mrSges(6,2) * t926;
t918 = mrSges(6,1) * t979 - mrSges(6,3) * t926;
t845 = m(6) * t869 - mrSges(6,2) * t968 + mrSges(6,3) * t891 + t909 * t925 - t918 * t979 + t1013;
t868 = t1001 * t872 - t874 * t996;
t865 = -pkin(5) * t968 - pkin(11) * t978 + t910 * t926 - t868;
t1011 = -m(7) * t865 + t878 * mrSges(7,1) - mrSges(7,2) * t879 + t912 * t900 - t901 * t913;
t917 = -mrSges(6,2) * t979 + mrSges(6,3) * t925;
t855 = m(6) * t868 + mrSges(6,1) * t968 - mrSges(6,3) * t892 - t909 * t926 + t917 * t979 + t1011;
t839 = t1001 * t855 + t996 * t845;
t927 = -mrSges(5,1) * t945 + mrSges(5,2) * t946;
t934 = -mrSges(5,2) * t981 + mrSges(5,3) * t945;
t837 = m(5) * t875 + mrSges(5,1) * t969 - mrSges(5,3) * t916 - t927 * t946 + t934 * t981 + t839;
t1014 = t1001 * t845 - t855 * t996;
t935 = mrSges(5,1) * t981 - mrSges(5,3) * t946;
t838 = m(5) * t876 - mrSges(5,2) * t969 + mrSges(5,3) * t915 + t927 * t945 - t935 * t981 + t1014;
t831 = t1002 * t837 + t997 * t838;
t947 = -mrSges(4,1) * t964 + mrSges(4,2) * t965;
t949 = mrSges(4,2) * t1018 + mrSges(4,3) * t964;
t829 = m(4) * t907 - mrSges(4,1) * t977 - mrSges(4,3) * t952 - t1018 * t949 - t947 * t965 + t831;
t1015 = t1002 * t838 - t837 * t997;
t950 = -mrSges(4,1) * t1018 - mrSges(4,3) * t965;
t830 = m(4) * t908 + mrSges(4,2) * t977 + mrSges(4,3) * t951 + t1018 * t950 + t947 * t964 + t1015;
t1016 = -t829 * t991 + t993 * t830;
t943 = -g(3) * t1028 + t1025;
t970 = mrSges(3,1) * t987 - mrSges(3,3) * t1019;
t975 = (-mrSges(3,1) * t1003 + mrSges(3,2) * t998) * t1024;
t821 = m(3) * t943 - mrSges(3,2) * t986 + mrSges(3,3) * t977 + t1018 * t975 - t970 * t987 + t1016;
t824 = t993 * t829 + t991 * t830;
t957 = -t992 * t972 - t1029;
t971 = -mrSges(3,2) * t987 + mrSges(3,3) * t1018;
t823 = m(3) * t957 - t977 * mrSges(3,1) + t976 * mrSges(3,2) + (-t1003 * t971 + t970 * t998) * t1024 + t824;
t848 = t1000 * t859 + t995 * t860;
t1012 = m(6) * t881 - t891 * mrSges(6,1) + t892 * mrSges(6,2) - t925 * t917 + t926 * t918 + t848;
t1008 = m(5) * t911 - t915 * mrSges(5,1) + t916 * mrSges(5,2) - t945 * t934 + t946 * t935 + t1012;
t846 = m(4) * t931 - t951 * mrSges(4,1) + t952 * mrSges(4,2) - t964 * t949 + t965 * t950 + t1008;
t842 = m(3) * t942 + t986 * mrSges(3,1) - t976 * mrSges(3,3) - t1019 * t975 + t987 * t971 - t846;
t811 = t842 * t1022 + t821 * t1027 - t823 * t992;
t808 = m(2) * t982 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1005 + t811;
t816 = t1003 * t821 - t842 * t998;
t814 = m(2) * t983 - mrSges(2,1) * t1005 - qJDD(1) * mrSges(2,2) + t816;
t1026 = t1004 * t808 + t999 * t814;
t810 = t842 * t1023 + t821 * t1028 + t994 * t823;
t1017 = t1004 * t814 - t808 * t999;
t887 = Ifges(7,5) * t913 + Ifges(7,6) * t912 + Ifges(7,3) * t924;
t889 = Ifges(7,1) * t913 + Ifges(7,4) * t912 + Ifges(7,5) * t924;
t852 = -mrSges(7,1) * t865 + mrSges(7,3) * t864 + Ifges(7,4) * t879 + Ifges(7,2) * t878 + Ifges(7,6) * t886 - t887 * t913 + t889 * t924;
t888 = Ifges(7,4) * t913 + Ifges(7,2) * t912 + Ifges(7,6) * t924;
t853 = mrSges(7,2) * t865 - mrSges(7,3) * t863 + Ifges(7,1) * t879 + Ifges(7,4) * t878 + Ifges(7,5) * t886 + t887 * t912 - t888 * t924;
t902 = Ifges(6,5) * t926 + Ifges(6,6) * t925 + Ifges(6,3) * t979;
t903 = Ifges(6,4) * t926 + Ifges(6,2) * t925 + Ifges(6,6) * t979;
t832 = mrSges(6,2) * t881 - mrSges(6,3) * t868 + Ifges(6,1) * t892 + Ifges(6,4) * t891 + Ifges(6,5) * t968 - pkin(11) * t848 + t1000 * t853 - t852 * t995 + t902 * t925 - t903 * t979;
t1007 = mrSges(7,1) * t863 - mrSges(7,2) * t864 + Ifges(7,5) * t879 + Ifges(7,6) * t878 + Ifges(7,3) * t886 + t888 * t913 - t889 * t912;
t904 = Ifges(6,1) * t926 + Ifges(6,4) * t925 + Ifges(6,5) * t979;
t833 = -mrSges(6,1) * t881 + mrSges(6,3) * t869 + Ifges(6,4) * t892 + Ifges(6,2) * t891 + Ifges(6,6) * t968 - pkin(5) * t848 - t902 * t926 + t904 * t979 - t1007;
t919 = Ifges(5,5) * t946 + Ifges(5,6) * t945 + Ifges(5,3) * t981;
t921 = Ifges(5,1) * t946 + Ifges(5,4) * t945 + Ifges(5,5) * t981;
t817 = -mrSges(5,1) * t911 + mrSges(5,3) * t876 + Ifges(5,4) * t916 + Ifges(5,2) * t915 + Ifges(5,6) * t969 - pkin(4) * t1012 + pkin(10) * t1014 + t1001 * t833 + t996 * t832 - t946 * t919 + t981 * t921;
t920 = Ifges(5,4) * t946 + Ifges(5,2) * t945 + Ifges(5,6) * t981;
t825 = mrSges(5,2) * t911 - mrSges(5,3) * t875 + Ifges(5,1) * t916 + Ifges(5,4) * t915 + Ifges(5,5) * t969 - pkin(10) * t839 + t1001 * t832 - t833 * t996 + t919 * t945 - t920 * t981;
t937 = Ifges(4,5) * t965 + Ifges(4,6) * t964 - Ifges(4,3) * t1018;
t939 = Ifges(4,1) * t965 + Ifges(4,4) * t964 - Ifges(4,5) * t1018;
t803 = -mrSges(4,1) * t931 + mrSges(4,3) * t908 + Ifges(4,4) * t952 + Ifges(4,2) * t951 - Ifges(4,6) * t977 - pkin(3) * t1008 + pkin(9) * t1015 + t1002 * t817 - t1018 * t939 + t997 * t825 - t965 * t937;
t938 = Ifges(4,4) * t965 + Ifges(4,2) * t964 - Ifges(4,6) * t1018;
t804 = mrSges(4,2) * t931 - mrSges(4,3) * t907 + Ifges(4,1) * t952 + Ifges(4,4) * t951 - Ifges(4,5) * t977 - pkin(9) * t831 + t1002 * t825 + t1018 * t938 - t817 * t997 + t937 * t964;
t955 = Ifges(3,6) * t987 + (Ifges(3,4) * t998 + Ifges(3,2) * t1003) * t1024;
t956 = Ifges(3,5) * t987 + (Ifges(3,1) * t998 + Ifges(3,4) * t1003) * t1024;
t800 = Ifges(3,5) * t976 + Ifges(3,6) * t977 + Ifges(3,3) * t986 + mrSges(3,1) * t942 - mrSges(3,2) * t943 + t991 * t804 + t993 * t803 - pkin(2) * t846 + qJ(3) * t1016 + (-t1003 * t956 + t955 * t998) * t1024;
t954 = Ifges(3,3) * t987 + (Ifges(3,5) * t998 + Ifges(3,6) * t1003) * t1024;
t802 = mrSges(3,2) * t957 - mrSges(3,3) * t942 + Ifges(3,1) * t976 + Ifges(3,4) * t977 + Ifges(3,5) * t986 - qJ(3) * t824 + t1018 * t954 - t803 * t991 + t804 * t993 - t955 * t987;
t1009 = -mrSges(6,1) * t868 + mrSges(6,2) * t869 - Ifges(6,5) * t892 - Ifges(6,6) * t891 - Ifges(6,3) * t968 - pkin(5) * t1011 - pkin(11) * t1013 - t1000 * t852 - t995 * t853 - t926 * t903 + t925 * t904;
t1006 = mrSges(5,1) * t875 - mrSges(5,2) * t876 + Ifges(5,5) * t916 + Ifges(5,6) * t915 + Ifges(5,3) * t969 + pkin(4) * t839 + t946 * t920 - t945 * t921 - t1009;
t806 = -t954 * t1019 - pkin(3) * t831 + (Ifges(3,2) + Ifges(4,3)) * t977 - pkin(2) * t824 - mrSges(4,1) * t907 + mrSges(4,2) * t908 + mrSges(3,3) * t943 - Ifges(4,6) * t951 - Ifges(4,5) * t952 - mrSges(3,1) * t957 + t964 * t939 - t965 * t938 + Ifges(3,4) * t976 + Ifges(3,6) * t986 + t987 * t956 - t1006;
t1010 = mrSges(2,1) * t982 - mrSges(2,2) * t983 + Ifges(2,3) * qJDD(1) + pkin(1) * t811 + t806 * t1023 + t802 * t1028 + t816 * t1030 + t994 * t800;
t798 = -mrSges(2,2) * g(3) - mrSges(2,3) * t982 + Ifges(2,5) * qJDD(1) - t1005 * Ifges(2,6) + t1003 * t802 - t998 * t806 + (-t810 * t992 - t811 * t994) * pkin(8);
t797 = mrSges(2,1) * g(3) + mrSges(2,3) * t983 + t1005 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t810 - t992 * t800 + (pkin(8) * t816 + t1003 * t806 + t802 * t998) * t994;
t1 = [-m(1) * g(1) + t1017; -m(1) * g(2) + t1026; (-m(1) - m(2)) * g(3) + t810; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1026 + t1004 * t798 - t999 * t797; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1017 + t1004 * t797 + t999 * t798; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1010; t1010; t800; t846; t1006; -t1009; t1007;];
tauJB  = t1;
