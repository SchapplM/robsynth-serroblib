% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-05-07 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 03:58:34
% EndTime: 2019-05-07 03:58:56
% DurationCPUTime: 16.08s
% Computational Cost: add. (256826->364), mult. (542286->437), div. (0->0), fcn. (395838->10), ass. (0->153)
t1041 = -2 * qJD(4);
t1040 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t1039 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t1017 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t1038 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t1016 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t1037 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t1026 = cos(pkin(10));
t975 = sin(pkin(6));
t1005 = t975 * t1026;
t1027 = cos(pkin(6));
t1007 = qJ(4) * t1027;
t1031 = cos(qJ(2));
t1009 = t1031 * qJD(1);
t1000 = qJD(2) * t1009;
t977 = sin(qJ(2));
t1020 = qJD(1) * t977;
t1008 = qJD(2) * t1020;
t978 = sin(qJ(1));
t980 = cos(qJ(1));
t969 = t978 * g(1) - t980 * g(2);
t982 = qJD(1) ^ 2;
t953 = -qJDD(1) * pkin(1) - t982 * pkin(8) - t969;
t963 = t977 * qJDD(1) + t1000;
t964 = qJDD(1) * t1031 - t1008;
t916 = (-t1000 - t963) * pkin(9) + (-t964 + t1008) * pkin(2) + t953;
t970 = -g(1) * t980 - g(2) * t978;
t954 = -pkin(1) * t982 + qJDD(1) * pkin(8) + t970;
t944 = -t977 * g(3) + t1031 * t954;
t962 = (-pkin(2) * t1031 - pkin(9) * t977) * qJD(1);
t981 = qJD(2) ^ 2;
t925 = -t981 * pkin(2) + qJDD(2) * pkin(9) + t1009 * t962 + t944;
t976 = sin(qJ(3));
t979 = cos(qJ(3));
t891 = t979 * t916 - t976 * t925;
t959 = qJD(2) * t979 - t1020 * t976;
t933 = qJD(3) * t959 + qJDD(2) * t976 + t963 * t979;
t1025 = qJ(4) * t975;
t960 = qJD(2) * t976 + t1020 * t979;
t934 = -pkin(3) * t959 - t1025 * t960;
t998 = t1009 - qJD(3);
t995 = t998 * t975;
t988 = t1027 * t959 - t995;
t935 = t988 * qJ(4);
t958 = qJDD(3) - t964;
t870 = t958 * pkin(3) - t1007 * t933 - t960 * t934 - t935 * t998 + t891;
t892 = t976 * t916 + t979 * t925;
t940 = -pkin(3) * t998 - t1007 * t960;
t932 = -qJD(3) * t960 + qJDD(2) * t979 - t963 * t976;
t994 = t1027 * t932 + t958 * t975;
t871 = qJ(4) * t994 + t959 * t934 + t940 * t998 + t892;
t943 = -t1031 * g(3) - t977 * t954;
t924 = -qJDD(2) * pkin(2) - pkin(9) * t981 + t962 * t1020 - t943;
t874 = -pkin(3) * t932 - t1025 * t933 - t935 * t959 + t940 * t960 + t924;
t974 = sin(pkin(10));
t922 = t1026 * t960 + t974 * t988;
t999 = t1027 * t1026;
t864 = t1005 * t874 + t1041 * t922 + t870 * t999 - t974 * t871;
t921 = t1026 * t995 - t959 * t999 + t974 * t960;
t894 = pkin(4) * t921 - qJ(5) * t922;
t915 = t1027 * t958 - t975 * t932;
t938 = t1027 * t998 + t975 * t959;
t937 = t938 ^ 2;
t861 = -t915 * pkin(4) - t937 * qJ(5) + t922 * t894 + qJDD(5) - t864;
t896 = -mrSges(6,2) * t921 - mrSges(6,3) * t922;
t901 = t1026 * t933 + t974 * t994;
t1036 = -m(6) * t861 - t901 * mrSges(6,1) - t922 * t896;
t1010 = t1017 * t938 + t1039 * t921 - t1040 * t922;
t1011 = -t1016 * t938 + t1038 * t921 + t1039 * t922;
t893 = -mrSges(7,2) * t922 + mrSges(7,3) * t921;
t1022 = -t893 - t896;
t1024 = t921 * t938;
t1032 = 2 * qJD(6);
t855 = t938 * t1032 + (t921 * t922 - t915) * qJ(6) + (t901 - t1024) * pkin(5) + t861;
t907 = -mrSges(7,1) * t921 - mrSges(7,2) * t938;
t1001 = -m(7) * t855 + t915 * mrSges(7,3) - t938 * t907;
t853 = t901 * mrSges(7,1) + t922 * t893 - t1001;
t906 = mrSges(6,1) * t921 + mrSges(6,3) * t938;
t851 = t915 * mrSges(6,2) - t938 * t906 - t1036 + t853;
t1033 = -2 * qJD(5);
t900 = -t1005 * t958 - t932 * t999 + t974 * t933;
t904 = pkin(5) * t922 + qJ(6) * t938;
t918 = t921 * t1041;
t920 = t921 ^ 2;
t1006 = t974 * t1027;
t1013 = t975 * t974 * t874 + t870 * t1006 + t1026 * t871;
t991 = t937 * pkin(4) - t915 * qJ(5) - t1013;
t857 = -t900 * pkin(5) - t920 * qJ(6) - t921 * t894 + qJDD(6) + t918 + (t1033 - t904) * t938 - t991;
t860 = 0.2e1 * qJD(5) * t938 + ((2 * qJD(4)) + t894) * t921 + t991;
t865 = t918 + t1013;
t905 = mrSges(7,1) * t922 + mrSges(7,3) * t938;
t1015 = m(7) * t857 + t915 * mrSges(7,2) - t938 * t905;
t908 = mrSges(6,1) * t922 - mrSges(6,2) * t938;
t992 = -m(6) * t860 + t915 * mrSges(6,3) - t938 * t908 + t1015;
t833 = mrSges(5,1) * t864 - mrSges(5,2) * t865 + mrSges(6,2) * t861 - mrSges(6,3) * t860 + mrSges(7,2) * t857 - mrSges(7,3) * t855 - qJ(6) * t853 - pkin(4) * t851 + qJ(5) * t992 + t1011 * t922 + t1037 * t915 + t1017 * t901 + (qJ(5) * t1022 - t1010) * t921 + (qJ(5) * (-mrSges(6,1) - mrSges(7,1)) - t1016) * t900;
t1012 = t1016 * t921 - t1017 * t922 + t1037 * t938;
t866 = t1027 * t874 - t975 * t870 + qJDD(4);
t986 = (-t901 - t1024) * qJ(5) + t866 + (-pkin(4) * t938 + t1033) * t922;
t859 = -t920 * pkin(5) + t921 * t1032 - t922 * t904 + (pkin(4) + qJ(6)) * t900 + t986;
t1014 = m(7) * t859 + t900 * mrSges(7,3) + t921 * t907;
t863 = t900 * pkin(4) + t986;
t997 = m(6) * t863 - t901 * mrSges(6,3) - t922 * t908 + t1014;
t852 = -t900 * mrSges(6,2) - t901 * mrSges(7,2) - t922 * t905 - t921 * t906 + t997;
t854 = -t900 * mrSges(7,1) - t921 * t893 + t1015;
t834 = -mrSges(5,1) * t866 + mrSges(5,3) * t865 - mrSges(6,1) * t860 + mrSges(6,2) * t863 + mrSges(7,1) * t857 - mrSges(7,3) * t859 + pkin(5) * t854 - qJ(6) * t1014 - pkin(4) * t852 + t1010 * t938 + (qJ(6) * t905 + t1012) * t922 + t1016 * t915 + (mrSges(7,2) * qJ(6) + t1039) * t901 + t1038 * t900;
t1021 = mrSges(5,2) * t938 - mrSges(5,3) * t921 - t906;
t1028 = -mrSges(7,1) - mrSges(5,3);
t1029 = mrSges(5,1) - mrSges(6,2);
t895 = mrSges(5,1) * t921 + mrSges(5,2) * t922;
t846 = m(5) * t864 - t1021 * t938 + (-t893 - t895) * t922 + t1029 * t915 + t1028 * t901 + t1001 + t1036;
t903 = -mrSges(5,1) * t938 - mrSges(5,3) * t922;
t849 = m(5) * t865 - t915 * mrSges(5,2) + t938 * t903 + (-t895 + t1022) * t921 + (-mrSges(6,1) + t1028) * t900 + t992;
t850 = m(5) * t866 + (t903 - t905) * t922 + t1021 * t921 + (mrSges(5,2) - mrSges(7,2)) * t901 + t1029 * t900 + t997;
t840 = t849 * t1006 + t846 * t999 - t975 * t850;
t841 = mrSges(6,1) * t861 + mrSges(7,1) * t855 + mrSges(5,2) * t866 - mrSges(7,2) * t859 - mrSges(5,3) * t864 - mrSges(6,3) * t863 + pkin(5) * t853 - qJ(5) * t852 + t1011 * t938 + t1012 * t921 + t1017 * t915 - t1039 * t900 + t1040 * t901;
t844 = t1026 * t849 - t974 * t846;
t928 = Ifges(4,4) * t960 + Ifges(4,2) * t959 - Ifges(4,6) * t998;
t929 = Ifges(4,1) * t960 + Ifges(4,4) * t959 - Ifges(4,5) * t998;
t1035 = mrSges(4,1) * t891 - mrSges(4,2) * t892 + Ifges(4,5) * t933 + Ifges(4,6) * t932 + Ifges(4,3) * t958 + pkin(3) * t840 + t1027 * t833 + t960 * t928 - t959 * t929 + t975 * (qJ(4) * t844 + t1026 * t834 + t841 * t974);
t839 = t975 * (t1026 * t846 + t849 * t974) + t1027 * t850;
t927 = Ifges(4,5) * t960 + Ifges(4,6) * t959 - Ifges(4,3) * t998;
t818 = -mrSges(4,1) * t924 + mrSges(4,3) * t892 + Ifges(4,4) * t933 + Ifges(4,2) * t932 + Ifges(4,6) * t958 - pkin(3) * t839 + t1006 * t841 + t1007 * t844 - t975 * t833 + t834 * t999 - t960 * t927 - t929 * t998;
t819 = Ifges(4,1) * t933 + Ifges(4,4) * t932 + Ifges(4,5) * t958 + t959 * t927 + t998 * t928 + mrSges(4,2) * t924 - mrSges(4,3) * t891 + t1026 * t841 - t974 * t834 + (-t1027 * t840 - t975 * t839) * qJ(4);
t936 = -mrSges(4,1) * t959 + mrSges(4,2) * t960;
t941 = mrSges(4,2) * t998 + t959 * mrSges(4,3);
t837 = m(4) * t891 + t958 * mrSges(4,1) - t933 * mrSges(4,3) - t960 * t936 - t941 * t998 + t840;
t942 = -mrSges(4,1) * t998 - t960 * mrSges(4,3);
t843 = m(4) * t892 - t958 * mrSges(4,2) + t932 * mrSges(4,3) + t959 * t936 + t942 * t998 + t844;
t832 = -t976 * t837 + t979 * t843;
t838 = -m(4) * t924 + t932 * mrSges(4,1) - t933 * mrSges(4,2) + t959 * t941 - t960 * t942 - t839;
t951 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t977 + Ifges(3,2) * t1031) * qJD(1);
t952 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t977 + Ifges(3,4) * t1031) * qJD(1);
t1034 = mrSges(3,1) * t943 - mrSges(3,2) * t944 + Ifges(3,5) * t963 + Ifges(3,6) * t964 + Ifges(3,3) * qJDD(2) + pkin(2) * t838 + pkin(9) * t832 - qJD(1) * (t1031 * t952 - t951 * t977) + t979 * t818 + t976 * t819;
t961 = (-mrSges(3,1) * t1031 + mrSges(3,2) * t977) * qJD(1);
t967 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1020;
t830 = m(3) * t944 - qJDD(2) * mrSges(3,2) + t964 * mrSges(3,3) - qJD(2) * t967 + t1009 * t961 + t832;
t968 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1009;
t836 = m(3) * t943 + qJDD(2) * mrSges(3,1) - t963 * mrSges(3,3) + qJD(2) * t968 - t1020 * t961 + t838;
t1002 = t1031 * t830 - t836 * t977;
t822 = m(2) * t970 - mrSges(2,1) * t982 - qJDD(1) * mrSges(2,2) + t1002;
t831 = t837 * t979 + t843 * t976;
t985 = -m(3) * t953 + t964 * mrSges(3,1) - mrSges(3,2) * t963 + t968 * t1009 - t1020 * t967 - t831;
t826 = m(2) * t969 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t982 + t985;
t1023 = t978 * t822 + t980 * t826;
t824 = t1031 * t836 + t977 * t830;
t1003 = t980 * t822 - t826 * t978;
t950 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t977 + Ifges(3,6) * t1031) * qJD(1);
t815 = mrSges(3,2) * t953 - mrSges(3,3) * t943 + Ifges(3,1) * t963 + Ifges(3,4) * t964 + Ifges(3,5) * qJDD(2) - pkin(9) * t831 - qJD(2) * t951 + t1009 * t950 - t976 * t818 + t979 * t819;
t817 = -mrSges(3,1) * t953 + mrSges(3,3) * t944 + Ifges(3,4) * t963 + Ifges(3,2) * t964 + Ifges(3,6) * qJDD(2) - pkin(2) * t831 + qJD(2) * t952 - t1020 * t950 - t1035;
t990 = mrSges(2,1) * t969 - mrSges(2,2) * t970 + Ifges(2,3) * qJDD(1) + pkin(1) * t985 + pkin(8) * t1002 + t1031 * t817 + t977 * t815;
t813 = mrSges(2,1) * g(3) + mrSges(2,3) * t970 + t982 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t824 - t1034;
t812 = -mrSges(2,2) * g(3) - mrSges(2,3) * t969 + Ifges(2,5) * qJDD(1) - t982 * Ifges(2,6) - pkin(8) * t824 + t1031 * t815 - t977 * t817;
t1 = [-m(1) * g(1) + t1003; -m(1) * g(2) + t1023; (-m(1) - m(2)) * g(3) + t824; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1023 + t980 * t812 - t978 * t813; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1003 + t978 * t812 + t980 * t813; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t990; t990; t1034; t1035; t850; t851; t854;];
tauJB  = t1;
