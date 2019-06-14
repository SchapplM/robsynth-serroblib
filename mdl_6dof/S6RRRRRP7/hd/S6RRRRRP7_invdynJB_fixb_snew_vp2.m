% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP7
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
% Datum: 2019-05-08 05:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:23:32
% EndTime: 2019-05-08 05:24:16
% DurationCPUTime: 34.08s
% Computational Cost: add. (577611->380), mult. (1227329->475), div. (0->0), fcn. (990837->12), ass. (0->160)
t1036 = Ifges(6,1) + Ifges(7,1);
t1029 = Ifges(6,4) + Ifges(7,4);
t1028 = Ifges(6,5) + Ifges(7,5);
t1035 = Ifges(6,2) + Ifges(7,2);
t1027 = Ifges(6,6) + Ifges(7,6);
t1034 = Ifges(6,3) + Ifges(7,3);
t984 = sin(pkin(6));
t1015 = qJD(1) * t984;
t989 = sin(qJ(2));
t1010 = t989 * t1015;
t985 = cos(pkin(6));
t980 = t985 * qJD(1) + qJD(2);
t988 = sin(qJ(3));
t993 = cos(qJ(3));
t957 = -t1010 * t988 + t993 * t980;
t958 = t1010 * t993 + t988 * t980;
t987 = sin(qJ(4));
t992 = cos(qJ(4));
t944 = t987 * t957 + t992 * t958;
t994 = cos(qJ(2));
t1014 = qJD(1) * t994;
t1009 = t984 * t1014;
t974 = qJD(3) - t1009;
t972 = qJD(4) + t974;
t986 = sin(qJ(5));
t991 = cos(qJ(5));
t928 = -t986 * t944 + t991 * t972;
t929 = t991 * t944 + t986 * t972;
t943 = t992 * t957 - t987 * t958;
t942 = qJD(5) - t943;
t1018 = t1028 * t942 + t1029 * t928 + t1036 * t929;
t1019 = -t1027 * t942 - t1029 * t929 - t1035 * t928;
t1023 = t985 * t989;
t1032 = pkin(8) * t984;
t990 = sin(qJ(1));
t995 = cos(qJ(1));
t975 = t990 * g(1) - t995 * g(2);
t996 = qJD(1) ^ 2;
t965 = qJDD(1) * pkin(1) + t1032 * t996 + t975;
t1013 = qJDD(1) * t984;
t976 = -t995 * g(1) - t990 * g(2);
t966 = -t996 * pkin(1) + pkin(8) * t1013 + t976;
t1016 = t965 * t1023 + t994 * t966;
t968 = (-pkin(2) * t994 - pkin(9) * t989) * t1015;
t978 = t980 ^ 2;
t979 = t985 * qJDD(1) + qJDD(2);
t925 = -t978 * pkin(2) + t979 * pkin(9) + (-g(3) * t989 + t1014 * t968) * t984 + t1016;
t1031 = t985 * g(3);
t969 = (qJD(2) * t1014 + qJDD(1) * t989) * t984;
t970 = -qJD(2) * t1010 + t1013 * t994;
t926 = -t970 * pkin(2) - t969 * pkin(9) - t1031 + (-t965 + (pkin(2) * t989 - pkin(9) * t994) * t980 * qJD(1)) * t984;
t888 = -t988 * t925 + t993 * t926;
t939 = t957 * qJD(3) + t993 * t969 + t988 * t979;
t962 = qJDD(3) - t970;
t877 = (t957 * t974 - t939) * pkin(10) + (t957 * t958 + t962) * pkin(3) + t888;
t889 = t993 * t925 + t988 * t926;
t938 = -t958 * qJD(3) - t988 * t969 + t993 * t979;
t948 = t974 * pkin(3) - t958 * pkin(10);
t956 = t957 ^ 2;
t880 = -t956 * pkin(3) + t938 * pkin(10) - t974 * t948 + t889;
t875 = t987 * t877 + t992 * t880;
t920 = -t943 * pkin(4) - t944 * pkin(11);
t961 = qJDD(4) + t962;
t971 = t972 ^ 2;
t869 = -t971 * pkin(4) + t961 * pkin(11) + t943 * t920 + t875;
t1022 = t985 * t994;
t1024 = t984 * t994;
t940 = -g(3) * t1024 + t1022 * t965 - t989 * t966;
t924 = -t979 * pkin(2) - t978 * pkin(9) + t968 * t1010 - t940;
t886 = -t938 * pkin(3) - t956 * pkin(10) + t958 * t948 + t924;
t904 = -t944 * qJD(4) + t992 * t938 - t987 * t939;
t905 = t943 * qJD(4) + t987 * t938 + t992 * t939;
t872 = (-t943 * t972 - t905) * pkin(11) + (t944 * t972 - t904) * pkin(4) + t886;
t864 = -t986 * t869 + t991 * t872;
t885 = t928 * qJD(5) + t991 * t905 + t986 * t961;
t903 = qJDD(5) - t904;
t860 = -0.2e1 * qJD(6) * t929 + (t928 * t942 - t885) * qJ(6) + (t928 * t929 + t903) * pkin(5) + t864;
t910 = -t942 * mrSges(7,2) + t928 * mrSges(7,3);
t1012 = m(7) * t860 + t903 * mrSges(7,1) + t942 * t910;
t908 = -t928 * mrSges(7,1) + t929 * mrSges(7,2);
t858 = -t885 * mrSges(7,3) - t929 * t908 + t1012;
t865 = t991 * t869 + t986 * t872;
t884 = -t929 * qJD(5) - t986 * t905 + t991 * t961;
t912 = t942 * pkin(5) - t929 * qJ(6);
t927 = t928 ^ 2;
t863 = -t927 * pkin(5) + t884 * qJ(6) + 0.2e1 * qJD(6) * t928 - t942 * t912 + t865;
t1033 = mrSges(6,1) * t864 + mrSges(7,1) * t860 - mrSges(6,2) * t865 - mrSges(7,2) * t863 + pkin(5) * t858 - t1018 * t928 - t1019 * t929 + t1027 * t884 + t1028 * t885 + t1034 * t903;
t1030 = -mrSges(6,2) - mrSges(7,2);
t1025 = t984 * t989;
t909 = -t928 * mrSges(6,1) + t929 * mrSges(6,2);
t911 = -t942 * mrSges(6,2) + t928 * mrSges(6,3);
t850 = m(6) * t864 + t903 * mrSges(6,1) + t942 * t911 + (-t908 - t909) * t929 + (-mrSges(6,3) - mrSges(7,3)) * t885 + t1012;
t1011 = m(7) * t863 + t884 * mrSges(7,3) + t928 * t908;
t913 = t942 * mrSges(7,1) - t929 * mrSges(7,3);
t1017 = -t942 * mrSges(6,1) + t929 * mrSges(6,3) - t913;
t853 = m(6) * t865 + t884 * mrSges(6,3) + t1017 * t942 + t1030 * t903 + t928 * t909 + t1011;
t1005 = -t986 * t850 + t991 * t853;
t919 = -t943 * mrSges(5,1) + t944 * mrSges(5,2);
t931 = t972 * mrSges(5,1) - t944 * mrSges(5,3);
t843 = m(5) * t875 - t961 * mrSges(5,2) + t904 * mrSges(5,3) + t943 * t919 - t972 * t931 + t1005;
t874 = t992 * t877 - t987 * t880;
t930 = -t972 * mrSges(5,2) + t943 * mrSges(5,3);
t868 = -t961 * pkin(4) - t971 * pkin(11) + t944 * t920 - t874;
t866 = -t884 * pkin(5) - t927 * qJ(6) + t929 * t912 + qJDD(6) + t868;
t1004 = -m(7) * t866 + t884 * mrSges(7,1) + t928 * t910;
t999 = -m(6) * t868 + t884 * mrSges(6,1) + t1017 * t929 + t1030 * t885 + t928 * t911 + t1004;
t855 = m(5) * t874 + t961 * mrSges(5,1) - t905 * mrSges(5,3) - t944 * t919 + t972 * t930 + t999;
t835 = t987 * t843 + t992 * t855;
t945 = -t957 * mrSges(4,1) + t958 * mrSges(4,2);
t946 = -t974 * mrSges(4,2) + t957 * mrSges(4,3);
t833 = m(4) * t888 + t962 * mrSges(4,1) - t939 * mrSges(4,3) - t958 * t945 + t974 * t946 + t835;
t1006 = t992 * t843 - t987 * t855;
t947 = t974 * mrSges(4,1) - t958 * mrSges(4,3);
t834 = m(4) * t889 - t962 * mrSges(4,2) + t938 * mrSges(4,3) + t957 * t945 - t974 * t947 + t1006;
t1007 = -t988 * t833 + t993 * t834;
t941 = -g(3) * t1025 + t1016;
t963 = t980 * mrSges(3,1) - mrSges(3,3) * t1010;
t967 = (-mrSges(3,1) * t994 + mrSges(3,2) * t989) * t1015;
t825 = m(3) * t941 - t979 * mrSges(3,2) + t970 * mrSges(3,3) + t1009 * t967 - t980 * t963 + t1007;
t828 = t993 * t833 + t988 * t834;
t952 = -t984 * t965 - t1031;
t964 = -t980 * mrSges(3,2) + mrSges(3,3) * t1009;
t827 = m(3) * t952 - t970 * mrSges(3,1) + t969 * mrSges(3,2) + (t963 * t989 - t964 * t994) * t1015 + t828;
t847 = t991 * t850 + t986 * t853;
t1002 = m(5) * t886 - t904 * mrSges(5,1) + t905 * mrSges(5,2) - t943 * t930 + t944 * t931 + t847;
t998 = -m(4) * t924 + t938 * mrSges(4,1) - t939 * mrSges(4,2) + t957 * t946 - t958 * t947 - t1002;
t840 = m(3) * t940 + t979 * mrSges(3,1) - t969 * mrSges(3,3) - t1010 * t967 + t980 * t964 + t998;
t814 = t840 * t1022 + t825 * t1023 - t984 * t827;
t811 = m(2) * t975 + qJDD(1) * mrSges(2,1) - t996 * mrSges(2,2) + t814;
t820 = t994 * t825 - t989 * t840;
t818 = m(2) * t976 - t996 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t820;
t1021 = t995 * t811 + t990 * t818;
t1020 = -t1027 * t928 - t1028 * t929 - t1034 * t942;
t813 = t840 * t1024 + t825 * t1025 + t985 * t827;
t1008 = -t990 * t811 + t995 * t818;
t861 = t885 * mrSges(7,2) + t929 * t913 - t1004;
t837 = -mrSges(6,1) * t868 + mrSges(6,3) * t865 - mrSges(7,1) * t866 + mrSges(7,3) * t863 - pkin(5) * t861 + qJ(6) * t1011 + (-qJ(6) * t913 + t1018) * t942 + t1020 * t929 + (-qJ(6) * mrSges(7,2) + t1027) * t903 + t1029 * t885 + t1035 * t884;
t845 = mrSges(6,2) * t868 + mrSges(7,2) * t866 - mrSges(6,3) * t864 - mrSges(7,3) * t860 - qJ(6) * t858 + t1019 * t942 - t1020 * t928 + t1028 * t903 + t1029 * t884 + t1036 * t885;
t915 = Ifges(5,5) * t944 + Ifges(5,6) * t943 + Ifges(5,3) * t972;
t916 = Ifges(5,4) * t944 + Ifges(5,2) * t943 + Ifges(5,6) * t972;
t821 = mrSges(5,2) * t886 - mrSges(5,3) * t874 + Ifges(5,1) * t905 + Ifges(5,4) * t904 + Ifges(5,5) * t961 - pkin(11) * t847 - t986 * t837 + t991 * t845 + t943 * t915 - t972 * t916;
t917 = Ifges(5,1) * t944 + Ifges(5,4) * t943 + Ifges(5,5) * t972;
t829 = -mrSges(5,1) * t886 + mrSges(5,3) * t875 + Ifges(5,4) * t905 + Ifges(5,2) * t904 + Ifges(5,6) * t961 - pkin(4) * t847 - t944 * t915 + t972 * t917 - t1033;
t932 = Ifges(4,5) * t958 + Ifges(4,6) * t957 + Ifges(4,3) * t974;
t934 = Ifges(4,1) * t958 + Ifges(4,4) * t957 + Ifges(4,5) * t974;
t809 = -mrSges(4,1) * t924 + mrSges(4,3) * t889 + Ifges(4,4) * t939 + Ifges(4,2) * t938 + Ifges(4,6) * t962 - pkin(3) * t1002 + pkin(10) * t1006 + t987 * t821 + t992 * t829 - t958 * t932 + t974 * t934;
t933 = Ifges(4,4) * t958 + Ifges(4,2) * t957 + Ifges(4,6) * t974;
t815 = mrSges(4,2) * t924 - mrSges(4,3) * t888 + Ifges(4,1) * t939 + Ifges(4,4) * t938 + Ifges(4,5) * t962 - pkin(10) * t835 + t992 * t821 - t987 * t829 + t957 * t932 - t974 * t933;
t950 = Ifges(3,6) * t980 + (Ifges(3,4) * t989 + Ifges(3,2) * t994) * t1015;
t951 = Ifges(3,5) * t980 + (Ifges(3,1) * t989 + Ifges(3,4) * t994) * t1015;
t804 = Ifges(3,5) * t969 + Ifges(3,6) * t970 + Ifges(3,3) * t979 + mrSges(3,1) * t940 - mrSges(3,2) * t941 + t988 * t815 + t993 * t809 + pkin(2) * t998 + pkin(9) * t1007 + (t950 * t989 - t951 * t994) * t1015;
t949 = Ifges(3,3) * t980 + (Ifges(3,5) * t989 + Ifges(3,6) * t994) * t1015;
t806 = mrSges(3,2) * t952 - mrSges(3,3) * t940 + Ifges(3,1) * t969 + Ifges(3,4) * t970 + Ifges(3,5) * t979 - pkin(9) * t828 + t1009 * t949 - t988 * t809 + t993 * t815 - t980 * t950;
t1000 = -mrSges(5,1) * t874 + mrSges(5,2) * t875 - Ifges(5,5) * t905 - Ifges(5,6) * t904 - Ifges(5,3) * t961 - pkin(4) * t999 - pkin(11) * t1005 - t991 * t837 - t986 * t845 - t944 * t916 + t943 * t917;
t997 = mrSges(4,1) * t888 - mrSges(4,2) * t889 + Ifges(4,5) * t939 + Ifges(4,6) * t938 + Ifges(4,3) * t962 + pkin(3) * t835 + t958 * t933 - t957 * t934 - t1000;
t808 = -mrSges(3,1) * t952 + mrSges(3,3) * t941 + Ifges(3,4) * t969 + Ifges(3,2) * t970 + Ifges(3,6) * t979 - pkin(2) * t828 - t1010 * t949 + t980 * t951 - t997;
t1003 = mrSges(2,1) * t975 - mrSges(2,2) * t976 + Ifges(2,3) * qJDD(1) + pkin(1) * t814 + t808 * t1024 + t806 * t1025 + t1032 * t820 + t985 * t804;
t802 = -mrSges(2,2) * g(3) - mrSges(2,3) * t975 + Ifges(2,5) * qJDD(1) - t996 * Ifges(2,6) + t994 * t806 - t989 * t808 + (-t813 * t984 - t814 * t985) * pkin(8);
t801 = mrSges(2,1) * g(3) + mrSges(2,3) * t976 + t996 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t813 - t984 * t804 + (pkin(8) * t820 + t806 * t989 + t808 * t994) * t985;
t1 = [-m(1) * g(1) + t1008; -m(1) * g(2) + t1021; (-m(1) - m(2)) * g(3) + t813; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1021 - t990 * t801 + t995 * t802; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1008 + t995 * t801 + t990 * t802; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1003; t1003; t804; t997; -t1000; t1033; t861;];
tauJB  = t1;
