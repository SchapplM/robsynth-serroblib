% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:40:52
% EndTime: 2019-05-08 05:41:36
% DurationCPUTime: 33.17s
% Computational Cost: add. (566359->378), mult. (1203790->475), div. (0->0), fcn. (970422->12), ass. (0->159)
t1033 = Ifges(6,1) + Ifges(7,1);
t1025 = Ifges(6,4) - Ifges(7,5);
t1024 = -Ifges(6,5) - Ifges(7,4);
t1032 = Ifges(6,2) + Ifges(7,3);
t1023 = Ifges(6,6) - Ifges(7,6);
t1031 = -Ifges(6,3) - Ifges(7,2);
t1029 = cos(qJ(5));
t991 = cos(qJ(2));
t1010 = qJD(1) * t991;
t983 = cos(pkin(6));
t987 = sin(qJ(2));
t1019 = t983 * t987;
t982 = sin(pkin(6));
t1028 = pkin(8) * t982;
t988 = sin(qJ(1));
t992 = cos(qJ(1));
t973 = t988 * g(1) - t992 * g(2);
t993 = qJD(1) ^ 2;
t962 = qJDD(1) * pkin(1) + t993 * t1028 + t973;
t1009 = qJDD(1) * t982;
t974 = -t992 * g(1) - t988 * g(2);
t963 = -t993 * pkin(1) + pkin(8) * t1009 + t974;
t1012 = t962 * t1019 + t991 * t963;
t1011 = qJD(1) * t982;
t965 = (-pkin(2) * t991 - pkin(9) * t987) * t1011;
t978 = t983 * qJD(1) + qJD(2);
t976 = t978 ^ 2;
t977 = t983 * qJDD(1) + qJDD(2);
t921 = -t976 * pkin(2) + t977 * pkin(9) + (-g(3) * t987 + t965 * t1010) * t982 + t1012;
t1027 = t983 * g(3);
t966 = (qJD(2) * t1010 + qJDD(1) * t987) * t982;
t1007 = t987 * t1011;
t967 = -qJD(2) * t1007 + t991 * t1009;
t922 = -t967 * pkin(2) - t966 * pkin(9) - t1027 + (-t962 + (pkin(2) * t987 - pkin(9) * t991) * t978 * qJD(1)) * t982;
t986 = sin(qJ(3));
t990 = cos(qJ(3));
t882 = -t986 * t921 + t990 * t922;
t954 = -t986 * t1007 + t990 * t978;
t934 = t954 * qJD(3) + t990 * t966 + t986 * t977;
t955 = t990 * t1007 + t986 * t978;
t959 = qJDD(3) - t967;
t1006 = t982 * t1010;
t972 = qJD(3) - t1006;
t873 = (t954 * t972 - t934) * pkin(10) + (t954 * t955 + t959) * pkin(3) + t882;
t883 = t990 * t921 + t986 * t922;
t933 = -t955 * qJD(3) - t986 * t966 + t990 * t977;
t944 = t972 * pkin(3) - t955 * pkin(10);
t953 = t954 ^ 2;
t876 = -t953 * pkin(3) + t933 * pkin(10) - t972 * t944 + t883;
t985 = sin(qJ(4));
t989 = cos(qJ(4));
t871 = t985 * t873 + t989 * t876;
t939 = t989 * t954 - t985 * t955;
t940 = t985 * t954 + t989 * t955;
t916 = -pkin(4) * t939 - pkin(11) * t940;
t958 = qJDD(4) + t959;
t970 = qJD(4) + t972;
t969 = t970 ^ 2;
t866 = -t969 * pkin(4) + t958 * pkin(11) + t939 * t916 + t871;
t1018 = t983 * t991;
t1020 = t982 * t991;
t935 = -g(3) * t1020 + t962 * t1018 - t987 * t963;
t920 = -t977 * pkin(2) - t976 * pkin(9) + t965 * t1007 - t935;
t881 = -t933 * pkin(3) - t953 * pkin(10) + t955 * t944 + t920;
t900 = -t940 * qJD(4) + t989 * t933 - t985 * t934;
t901 = t939 * qJD(4) + t985 * t933 + t989 * t934;
t868 = (-t939 * t970 - t901) * pkin(11) + (t940 * t970 - t900) * pkin(4) + t881;
t984 = sin(qJ(5));
t863 = t1029 * t866 + t984 * t868;
t899 = qJDD(5) - t900;
t923 = -t1029 * t970 + t984 * t940;
t924 = t1029 * t940 + t984 * t970;
t904 = pkin(5) * t923 - qJ(6) * t924;
t938 = qJD(5) - t939;
t937 = t938 ^ 2;
t859 = -pkin(5) * t937 + qJ(6) * t899 + 0.2e1 * qJD(6) * t938 - t904 * t923 + t863;
t910 = -mrSges(7,1) * t938 + mrSges(7,2) * t924;
t1008 = m(7) * t859 + t899 * mrSges(7,3) + t938 * t910;
t1014 = t1024 * t938 + t1025 * t923 - t1033 * t924;
t1015 = -t1023 * t938 - t1025 * t924 + t1032 * t923;
t862 = t1029 * t868 - t984 * t866;
t860 = -t899 * pkin(5) - t937 * qJ(6) + t924 * t904 + qJDD(6) - t862;
t907 = -mrSges(7,2) * t923 + mrSges(7,3) * t938;
t1001 = -m(7) * t860 + t899 * mrSges(7,1) + t938 * t907;
t880 = -t923 * qJD(5) + t1029 * t901 + t984 * t958;
t905 = mrSges(7,1) * t923 - mrSges(7,3) * t924;
t856 = mrSges(7,2) * t880 + t905 * t924 - t1001;
t879 = t924 * qJD(5) - t1029 * t958 + t984 * t901;
t1030 = -t1014 * t923 - t1015 * t924 - t1031 * t899 - t1023 * t879 - t1024 * t880 + mrSges(6,1) * t862 - mrSges(7,1) * t860 - mrSges(6,2) * t863 + mrSges(7,3) * t859 - pkin(5) * t856 + qJ(6) * (-mrSges(7,2) * t879 - t905 * t923 + t1008);
t1026 = -mrSges(6,3) - mrSges(7,2);
t1021 = t982 * t987;
t1013 = -mrSges(6,1) * t923 - mrSges(6,2) * t924 - t905;
t909 = mrSges(6,1) * t938 - mrSges(6,3) * t924;
t850 = m(6) * t863 - mrSges(6,2) * t899 + t1013 * t923 + t1026 * t879 - t909 * t938 + t1008;
t908 = -mrSges(6,2) * t938 - mrSges(6,3) * t923;
t852 = m(6) * t862 + mrSges(6,1) * t899 + t1013 * t924 + t1026 * t880 + t908 * t938 + t1001;
t1002 = t1029 * t850 - t984 * t852;
t915 = -mrSges(5,1) * t939 + mrSges(5,2) * t940;
t926 = t970 * mrSges(5,1) - t940 * mrSges(5,3);
t838 = m(5) * t871 - t958 * mrSges(5,2) + t900 * mrSges(5,3) + t939 * t915 - t970 * t926 + t1002;
t870 = t989 * t873 - t985 * t876;
t925 = -t970 * mrSges(5,2) + t939 * mrSges(5,3);
t865 = -t958 * pkin(4) - t969 * pkin(11) + t940 * t916 - t870;
t861 = -0.2e1 * qJD(6) * t924 + (t923 * t938 - t880) * qJ(6) + (t924 * t938 + t879) * pkin(5) + t865;
t857 = m(7) * t861 + t879 * mrSges(7,1) - t880 * mrSges(7,3) + t923 * t907 - t924 * t910;
t996 = -m(6) * t865 - t879 * mrSges(6,1) - t880 * mrSges(6,2) - t923 * t908 - t924 * t909 - t857;
t847 = m(5) * t870 + t958 * mrSges(5,1) - t901 * mrSges(5,3) - t940 * t915 + t970 * t925 + t996;
t832 = t985 * t838 + t989 * t847;
t941 = -mrSges(4,1) * t954 + mrSges(4,2) * t955;
t942 = -t972 * mrSges(4,2) + t954 * mrSges(4,3);
t830 = m(4) * t882 + t959 * mrSges(4,1) - t934 * mrSges(4,3) - t955 * t941 + t972 * t942 + t832;
t1003 = t989 * t838 - t985 * t847;
t943 = t972 * mrSges(4,1) - t955 * mrSges(4,3);
t831 = m(4) * t883 - t959 * mrSges(4,2) + t933 * mrSges(4,3) + t954 * t941 - t972 * t943 + t1003;
t1004 = -t986 * t830 + t990 * t831;
t936 = -g(3) * t1021 + t1012;
t960 = t978 * mrSges(3,1) - mrSges(3,3) * t1007;
t964 = (-mrSges(3,1) * t991 + mrSges(3,2) * t987) * t1011;
t821 = m(3) * t936 - t977 * mrSges(3,2) + t967 * mrSges(3,3) + t964 * t1006 - t978 * t960 + t1004;
t824 = t990 * t830 + t986 * t831;
t948 = -t982 * t962 - t1027;
t961 = -t978 * mrSges(3,2) + mrSges(3,3) * t1006;
t823 = m(3) * t948 - t967 * mrSges(3,1) + t966 * mrSges(3,2) + (t960 * t987 - t961 * t991) * t1011 + t824;
t844 = t1029 * t852 + t984 * t850;
t999 = m(5) * t881 - t900 * mrSges(5,1) + t901 * mrSges(5,2) - t939 * t925 + t940 * t926 + t844;
t995 = -m(4) * t920 + t933 * mrSges(4,1) - t934 * mrSges(4,2) + t954 * t942 - t955 * t943 - t999;
t835 = m(3) * t935 + t977 * mrSges(3,1) - t966 * mrSges(3,3) - t964 * t1007 + t978 * t961 + t995;
t811 = t835 * t1018 + t821 * t1019 - t982 * t823;
t808 = m(2) * t973 + qJDD(1) * mrSges(2,1) - t993 * mrSges(2,2) + t811;
t817 = t991 * t821 - t987 * t835;
t815 = m(2) * t974 - t993 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t817;
t1017 = t992 * t808 + t988 * t815;
t1016 = t1023 * t923 + t1024 * t924 + t1031 * t938;
t810 = t835 * t1020 + t821 * t1021 + t983 * t823;
t1005 = -t988 * t808 + t992 * t815;
t840 = -mrSges(6,1) * t865 - mrSges(7,1) * t861 + mrSges(7,2) * t859 + mrSges(6,3) * t863 - pkin(5) * t857 - t1014 * t938 + t1016 * t924 + t1023 * t899 + t1025 * t880 - t1032 * t879;
t842 = mrSges(6,2) * t865 + mrSges(7,2) * t860 - mrSges(6,3) * t862 - mrSges(7,3) * t861 - qJ(6) * t857 + t1015 * t938 + t1016 * t923 - t1024 * t899 - t1025 * t879 + t1033 * t880;
t911 = Ifges(5,5) * t940 + Ifges(5,6) * t939 + Ifges(5,3) * t970;
t912 = Ifges(5,4) * t940 + Ifges(5,2) * t939 + Ifges(5,6) * t970;
t825 = mrSges(5,2) * t881 - mrSges(5,3) * t870 + Ifges(5,1) * t901 + Ifges(5,4) * t900 + Ifges(5,5) * t958 - pkin(11) * t844 + t1029 * t842 - t984 * t840 + t939 * t911 - t970 * t912;
t913 = Ifges(5,1) * t940 + Ifges(5,4) * t939 + Ifges(5,5) * t970;
t826 = -mrSges(5,1) * t881 + mrSges(5,3) * t871 + Ifges(5,4) * t901 + Ifges(5,2) * t900 + Ifges(5,6) * t958 - pkin(4) * t844 - t940 * t911 + t970 * t913 - t1030;
t927 = Ifges(4,5) * t955 + Ifges(4,6) * t954 + Ifges(4,3) * t972;
t929 = Ifges(4,1) * t955 + Ifges(4,4) * t954 + Ifges(4,5) * t972;
t806 = -mrSges(4,1) * t920 + mrSges(4,3) * t883 + Ifges(4,4) * t934 + Ifges(4,2) * t933 + Ifges(4,6) * t959 - pkin(3) * t999 + pkin(10) * t1003 + t985 * t825 + t989 * t826 - t955 * t927 + t972 * t929;
t928 = Ifges(4,4) * t955 + Ifges(4,2) * t954 + Ifges(4,6) * t972;
t812 = mrSges(4,2) * t920 - mrSges(4,3) * t882 + Ifges(4,1) * t934 + Ifges(4,4) * t933 + Ifges(4,5) * t959 - pkin(10) * t832 + t989 * t825 - t985 * t826 + t954 * t927 - t972 * t928;
t946 = Ifges(3,6) * t978 + (Ifges(3,4) * t987 + Ifges(3,2) * t991) * t1011;
t947 = Ifges(3,5) * t978 + (Ifges(3,1) * t987 + Ifges(3,4) * t991) * t1011;
t801 = Ifges(3,5) * t966 + Ifges(3,6) * t967 + Ifges(3,3) * t977 + mrSges(3,1) * t935 - mrSges(3,2) * t936 + t986 * t812 + t990 * t806 + pkin(2) * t995 + pkin(9) * t1004 + (t946 * t987 - t947 * t991) * t1011;
t945 = Ifges(3,3) * t978 + (Ifges(3,5) * t987 + Ifges(3,6) * t991) * t1011;
t803 = mrSges(3,2) * t948 - mrSges(3,3) * t935 + Ifges(3,1) * t966 + Ifges(3,4) * t967 + Ifges(3,5) * t977 - pkin(9) * t824 + t945 * t1006 - t986 * t806 + t990 * t812 - t978 * t946;
t998 = -mrSges(5,1) * t870 + mrSges(5,2) * t871 - Ifges(5,5) * t901 - Ifges(5,6) * t900 - Ifges(5,3) * t958 - pkin(4) * t996 - pkin(11) * t1002 - t1029 * t840 - t984 * t842 - t940 * t912 + t939 * t913;
t994 = mrSges(4,1) * t882 - mrSges(4,2) * t883 + Ifges(4,5) * t934 + Ifges(4,6) * t933 + Ifges(4,3) * t959 + pkin(3) * t832 + t955 * t928 - t954 * t929 - t998;
t805 = -mrSges(3,1) * t948 + mrSges(3,3) * t936 + Ifges(3,4) * t966 + Ifges(3,2) * t967 + Ifges(3,6) * t977 - pkin(2) * t824 - t945 * t1007 + t978 * t947 - t994;
t1000 = mrSges(2,1) * t973 - mrSges(2,2) * t974 + Ifges(2,3) * qJDD(1) + pkin(1) * t811 + t805 * t1020 + t803 * t1021 + t817 * t1028 + t983 * t801;
t799 = -mrSges(2,2) * g(3) - mrSges(2,3) * t973 + Ifges(2,5) * qJDD(1) - t993 * Ifges(2,6) + t991 * t803 - t987 * t805 + (-t810 * t982 - t811 * t983) * pkin(8);
t798 = mrSges(2,1) * g(3) + mrSges(2,3) * t974 + t993 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t810 - t982 * t801 + (pkin(8) * t817 + t803 * t987 + t805 * t991) * t983;
t1 = [-m(1) * g(1) + t1005; -m(1) * g(2) + t1017; (-m(1) - m(2)) * g(3) + t810; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1017 - t988 * t798 + t992 * t799; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1005 + t992 * t798 + t988 * t799; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1000; t1000; t801; t994; -t998; t1030; t856;];
tauJB  = t1;
