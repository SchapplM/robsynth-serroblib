% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 10:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:02:17
% EndTime: 2019-05-06 10:03:12
% DurationCPUTime: 52.73s
% Computational Cost: add. (795501->398), mult. (2149006->518), div. (0->0), fcn. (1713502->14), ass. (0->163)
t1025 = -2 * qJD(3);
t1023 = cos(pkin(11));
t991 = cos(pkin(6));
t998 = cos(qJ(2));
t1019 = t991 * t998;
t1000 = qJD(1) ^ 2;
t989 = sin(pkin(6));
t1024 = pkin(8) * t989;
t995 = sin(qJ(1));
t999 = cos(qJ(1));
t978 = t995 * g(1) - t999 * g(2);
t971 = qJDD(1) * pkin(1) + t1000 * t1024 + t978;
t979 = -t999 * g(1) - t995 * g(2);
t972 = -t1000 * pkin(1) + qJDD(1) * t1024 + t979;
t994 = sin(qJ(2));
t1006 = t971 * t1019 - t994 * t972;
t1015 = t1000 * t989 ^ 2;
t1014 = qJD(1) * qJD(2);
t974 = (qJDD(1) * t994 + t1014 * t998) * t989;
t981 = t991 * qJDD(1) + qJDD(2);
t982 = t991 * qJD(1) + qJD(2);
t912 = t981 * pkin(2) - t974 * qJ(3) + (pkin(2) * t994 * t1015 + (qJ(3) * qJD(1) * t982 - g(3)) * t989) * t998 + t1006;
t1011 = t998 ^ 2 * t1015;
t1020 = t991 * t994;
t1022 = t989 * t994;
t939 = -g(3) * t1022 + t971 * t1020 + t998 * t972;
t1017 = qJD(1) * t989;
t1013 = t994 * t1017;
t968 = t982 * pkin(2) - qJ(3) * t1013;
t975 = (qJDD(1) * t998 - t1014 * t994) * t989;
t915 = -pkin(2) * t1011 + t975 * qJ(3) - t982 * t968 + t939;
t988 = sin(pkin(11));
t965 = (t1023 * t994 + t988 * t998) * t1017;
t889 = t1023 * t912 + t965 * t1025 - t988 * t915;
t1021 = t989 * t998;
t1012 = t998 * t1017;
t964 = -t1023 * t1012 + t1013 * t988;
t890 = t1023 * t915 + t964 * t1025 + t988 * t912;
t940 = t964 * pkin(3) - t965 * qJ(4);
t980 = t982 ^ 2;
t884 = -t980 * pkin(3) + t981 * qJ(4) - t964 * t940 + t890;
t956 = -t991 * g(3) - t989 * t971;
t923 = -t975 * pkin(2) - qJ(3) * t1011 + t968 * t1013 + qJDD(3) + t956;
t944 = -t1023 * t975 + t988 * t974;
t945 = t1023 * t974 + t988 * t975;
t893 = (t964 * t982 - t945) * qJ(4) + (t965 * t982 + t944) * pkin(3) + t923;
t987 = sin(pkin(12));
t990 = cos(pkin(12));
t950 = t990 * t965 + t987 * t982;
t876 = -0.2e1 * qJD(4) * t950 - t987 * t884 + t990 * t893;
t933 = t990 * t945 + t987 * t981;
t949 = -t987 * t965 + t990 * t982;
t873 = (t949 * t964 - t933) * pkin(9) + (t949 * t950 + t944) * pkin(4) + t876;
t877 = 0.2e1 * qJD(4) * t949 + t990 * t884 + t987 * t893;
t930 = t964 * pkin(4) - t950 * pkin(9);
t932 = -t987 * t945 + t990 * t981;
t948 = t949 ^ 2;
t875 = -t948 * pkin(4) + t932 * pkin(9) - t964 * t930 + t877;
t993 = sin(qJ(5));
t997 = cos(qJ(5));
t870 = t993 * t873 + t997 * t875;
t924 = t997 * t949 - t993 * t950;
t925 = t993 * t949 + t997 * t950;
t906 = -t924 * pkin(5) - t925 * pkin(10);
t943 = qJDD(5) + t944;
t963 = qJD(5) + t964;
t962 = t963 ^ 2;
t868 = -t962 * pkin(5) + t943 * pkin(10) + t924 * t906 + t870;
t883 = -t981 * pkin(3) - t980 * qJ(4) + t965 * t940 + qJDD(4) - t889;
t878 = -t932 * pkin(4) - t948 * pkin(9) + t950 * t930 + t883;
t897 = -t925 * qJD(5) + t997 * t932 - t993 * t933;
t898 = t924 * qJD(5) + t993 * t932 + t997 * t933;
t871 = (-t924 * t963 - t898) * pkin(10) + (t925 * t963 - t897) * pkin(5) + t878;
t992 = sin(qJ(6));
t996 = cos(qJ(6));
t865 = -t992 * t868 + t996 * t871;
t908 = -t992 * t925 + t996 * t963;
t881 = t908 * qJD(6) + t996 * t898 + t992 * t943;
t909 = t996 * t925 + t992 * t963;
t894 = -t908 * mrSges(7,1) + t909 * mrSges(7,2);
t896 = qJDD(6) - t897;
t922 = qJD(6) - t924;
t899 = -t922 * mrSges(7,2) + t908 * mrSges(7,3);
t862 = m(7) * t865 + t896 * mrSges(7,1) - t881 * mrSges(7,3) - t909 * t894 + t922 * t899;
t866 = t996 * t868 + t992 * t871;
t880 = -t909 * qJD(6) - t992 * t898 + t996 * t943;
t900 = t922 * mrSges(7,1) - t909 * mrSges(7,3);
t863 = m(7) * t866 - t896 * mrSges(7,2) + t880 * mrSges(7,3) + t908 * t894 - t922 * t900;
t854 = -t992 * t862 + t996 * t863;
t905 = -t924 * mrSges(6,1) + t925 * mrSges(6,2);
t914 = t963 * mrSges(6,1) - t925 * mrSges(6,3);
t851 = m(6) * t870 - t943 * mrSges(6,2) + t897 * mrSges(6,3) + t924 * t905 - t963 * t914 + t854;
t869 = t997 * t873 - t993 * t875;
t867 = -t943 * pkin(5) - t962 * pkin(10) + t925 * t906 - t869;
t864 = -m(7) * t867 + t880 * mrSges(7,1) - t881 * mrSges(7,2) + t908 * t899 - t909 * t900;
t913 = -t963 * mrSges(6,2) + t924 * mrSges(6,3);
t858 = m(6) * t869 + t943 * mrSges(6,1) - t898 * mrSges(6,3) - t925 * t905 + t963 * t913 + t864;
t846 = t993 * t851 + t997 * t858;
t926 = -t949 * mrSges(5,1) + t950 * mrSges(5,2);
t928 = -t964 * mrSges(5,2) + t949 * mrSges(5,3);
t844 = m(5) * t876 + t944 * mrSges(5,1) - t933 * mrSges(5,3) - t950 * t926 + t964 * t928 + t846;
t1007 = t997 * t851 - t993 * t858;
t929 = t964 * mrSges(5,1) - t950 * mrSges(5,3);
t845 = m(5) * t877 - t944 * mrSges(5,2) + t932 * mrSges(5,3) + t949 * t926 - t964 * t929 + t1007;
t1008 = -t987 * t844 + t990 * t845;
t941 = t964 * mrSges(4,1) + t965 * mrSges(4,2);
t952 = t982 * mrSges(4,1) - t965 * mrSges(4,3);
t834 = m(4) * t890 - t981 * mrSges(4,2) - t944 * mrSges(4,3) - t964 * t941 - t982 * t952 + t1008;
t853 = t996 * t862 + t992 * t863;
t1003 = m(6) * t878 - t897 * mrSges(6,1) + t898 * mrSges(6,2) - t924 * t913 + t925 * t914 + t853;
t852 = m(5) * t883 - t932 * mrSges(5,1) + t933 * mrSges(5,2) - t949 * t928 + t950 * t929 + t1003;
t951 = -t982 * mrSges(4,2) - t964 * mrSges(4,3);
t848 = m(4) * t889 + t981 * mrSges(4,1) - t945 * mrSges(4,3) - t965 * t941 + t982 * t951 - t852;
t831 = t1023 * t848 + t988 * t834;
t938 = -g(3) * t1021 + t1006;
t970 = -t982 * mrSges(3,2) + mrSges(3,3) * t1012;
t973 = (-mrSges(3,1) * t998 + mrSges(3,2) * t994) * t1017;
t829 = m(3) * t938 + t981 * mrSges(3,1) - t974 * mrSges(3,3) - t1013 * t973 + t982 * t970 + t831;
t1009 = t1023 * t834 - t988 * t848;
t969 = t982 * mrSges(3,1) - mrSges(3,3) * t1013;
t830 = m(3) * t939 - t981 * mrSges(3,2) + t975 * mrSges(3,3) + t1012 * t973 - t982 * t969 + t1009;
t838 = t990 * t844 + t987 * t845;
t837 = m(4) * t923 + t944 * mrSges(4,1) + t945 * mrSges(4,2) + t964 * t951 + t965 * t952 + t838;
t836 = m(3) * t956 - t975 * mrSges(3,1) + t974 * mrSges(3,2) + (t969 * t994 - t970 * t998) * t1017 + t837;
t815 = t829 * t1019 + t830 * t1020 - t989 * t836;
t812 = m(2) * t978 + qJDD(1) * mrSges(2,1) - t1000 * mrSges(2,2) + t815;
t820 = -t994 * t829 + t998 * t830;
t818 = m(2) * t979 - t1000 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t820;
t1018 = t999 * t812 + t995 * t818;
t814 = t829 * t1021 + t830 * t1022 + t991 * t836;
t1010 = -t995 * t812 + t999 * t818;
t885 = Ifges(7,5) * t909 + Ifges(7,6) * t908 + Ifges(7,3) * t922;
t887 = Ifges(7,1) * t909 + Ifges(7,4) * t908 + Ifges(7,5) * t922;
t855 = -mrSges(7,1) * t867 + mrSges(7,3) * t866 + Ifges(7,4) * t881 + Ifges(7,2) * t880 + Ifges(7,6) * t896 - t909 * t885 + t922 * t887;
t886 = Ifges(7,4) * t909 + Ifges(7,2) * t908 + Ifges(7,6) * t922;
t856 = mrSges(7,2) * t867 - mrSges(7,3) * t865 + Ifges(7,1) * t881 + Ifges(7,4) * t880 + Ifges(7,5) * t896 + t908 * t885 - t922 * t886;
t901 = Ifges(6,5) * t925 + Ifges(6,6) * t924 + Ifges(6,3) * t963;
t902 = Ifges(6,4) * t925 + Ifges(6,2) * t924 + Ifges(6,6) * t963;
t839 = mrSges(6,2) * t878 - mrSges(6,3) * t869 + Ifges(6,1) * t898 + Ifges(6,4) * t897 + Ifges(6,5) * t943 - pkin(10) * t853 - t992 * t855 + t996 * t856 + t924 * t901 - t963 * t902;
t1002 = mrSges(7,1) * t865 - mrSges(7,2) * t866 + Ifges(7,5) * t881 + Ifges(7,6) * t880 + Ifges(7,3) * t896 + t909 * t886 - t908 * t887;
t903 = Ifges(6,1) * t925 + Ifges(6,4) * t924 + Ifges(6,5) * t963;
t840 = -mrSges(6,1) * t878 + mrSges(6,3) * t870 + Ifges(6,4) * t898 + Ifges(6,2) * t897 + Ifges(6,6) * t943 - pkin(5) * t853 - t925 * t901 + t963 * t903 - t1002;
t916 = Ifges(5,5) * t950 + Ifges(5,6) * t949 + Ifges(5,3) * t964;
t918 = Ifges(5,1) * t950 + Ifges(5,4) * t949 + Ifges(5,5) * t964;
t822 = -mrSges(5,1) * t883 + mrSges(5,3) * t877 + Ifges(5,4) * t933 + Ifges(5,2) * t932 + Ifges(5,6) * t944 - pkin(4) * t1003 + pkin(9) * t1007 + t993 * t839 + t997 * t840 - t950 * t916 + t964 * t918;
t917 = Ifges(5,4) * t950 + Ifges(5,2) * t949 + Ifges(5,6) * t964;
t823 = mrSges(5,2) * t883 - mrSges(5,3) * t876 + Ifges(5,1) * t933 + Ifges(5,4) * t932 + Ifges(5,5) * t944 - pkin(9) * t846 + t997 * t839 - t993 * t840 + t949 * t916 - t964 * t917;
t934 = Ifges(4,5) * t965 - Ifges(4,6) * t964 + Ifges(4,3) * t982;
t935 = Ifges(4,4) * t965 - Ifges(4,2) * t964 + Ifges(4,6) * t982;
t810 = mrSges(4,2) * t923 - mrSges(4,3) * t889 + Ifges(4,1) * t945 - Ifges(4,4) * t944 + Ifges(4,5) * t981 - qJ(4) * t838 - t987 * t822 + t990 * t823 - t964 * t934 - t982 * t935;
t1001 = mrSges(6,1) * t869 - mrSges(6,2) * t870 + Ifges(6,5) * t898 + Ifges(6,6) * t897 + Ifges(6,3) * t943 + pkin(5) * t864 + pkin(10) * t854 + t996 * t855 + t992 * t856 + t925 * t902 - t924 * t903;
t936 = Ifges(4,1) * t965 - Ifges(4,4) * t964 + Ifges(4,5) * t982;
t821 = -t1001 + (-Ifges(5,3) - Ifges(4,2)) * t944 - mrSges(5,1) * t876 + Ifges(4,4) * t945 - pkin(4) * t846 + mrSges(5,2) * t877 - mrSges(4,1) * t923 + t949 * t918 - t950 * t917 - Ifges(5,6) * t932 - Ifges(5,5) * t933 - t965 * t934 - pkin(3) * t838 + mrSges(4,3) * t890 + Ifges(4,6) * t981 + t982 * t936;
t953 = Ifges(3,3) * t982 + (Ifges(3,5) * t994 + Ifges(3,6) * t998) * t1017;
t955 = Ifges(3,5) * t982 + (Ifges(3,1) * t994 + Ifges(3,4) * t998) * t1017;
t805 = -mrSges(3,1) * t956 + mrSges(3,3) * t939 + Ifges(3,4) * t974 + Ifges(3,2) * t975 + Ifges(3,6) * t981 - pkin(2) * t837 + qJ(3) * t1009 - t953 * t1013 + t1023 * t821 + t988 * t810 + t982 * t955;
t954 = Ifges(3,6) * t982 + (Ifges(3,4) * t994 + Ifges(3,2) * t998) * t1017;
t807 = mrSges(3,2) * t956 - mrSges(3,3) * t938 + Ifges(3,1) * t974 + Ifges(3,4) * t975 + Ifges(3,5) * t981 - qJ(3) * t831 + t953 * t1012 + t1023 * t810 - t988 * t821 - t982 * t954;
t809 = Ifges(3,5) * t974 + Ifges(3,6) * t975 + mrSges(3,1) * t938 - mrSges(3,2) * t939 + Ifges(4,5) * t945 - Ifges(4,6) * t944 + t965 * t935 + t964 * t936 + mrSges(4,1) * t889 - mrSges(4,2) * t890 + t987 * t823 + t990 * t822 - pkin(3) * t852 + qJ(4) * t1008 + pkin(2) * t831 + (Ifges(3,3) + Ifges(4,3)) * t981 + (t954 * t994 - t955 * t998) * t1017;
t1004 = mrSges(2,1) * t978 - mrSges(2,2) * t979 + Ifges(2,3) * qJDD(1) + pkin(1) * t815 + t805 * t1021 + t807 * t1022 + t1024 * t820 + t991 * t809;
t803 = -mrSges(2,2) * g(3) - mrSges(2,3) * t978 + Ifges(2,5) * qJDD(1) - t1000 * Ifges(2,6) - t994 * t805 + t998 * t807 + (-t814 * t989 - t815 * t991) * pkin(8);
t802 = mrSges(2,1) * g(3) + mrSges(2,3) * t979 + t1000 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t814 - t989 * t809 + (pkin(8) * t820 + t805 * t998 + t807 * t994) * t991;
t1 = [-m(1) * g(1) + t1010; -m(1) * g(2) + t1018; (-m(1) - m(2)) * g(3) + t814; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1018 - t995 * t802 + t999 * t803; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1010 + t999 * t802 + t995 * t803; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1004; t1004; t809; t837; t852; t1001; t1002;];
tauJB  = t1;
