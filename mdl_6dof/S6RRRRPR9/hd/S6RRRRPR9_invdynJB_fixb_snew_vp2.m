% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 22:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:16:34
% EndTime: 2019-05-07 22:18:15
% DurationCPUTime: 72.74s
% Computational Cost: add. (1235337->401), mult. (2644091->518), div. (0->0), fcn. (2161265->14), ass. (0->166)
t1029 = cos(qJ(4));
t990 = sin(pkin(6));
t1028 = pkin(8) * t990;
t992 = cos(pkin(6));
t1027 = t992 * g(3);
t996 = sin(qJ(2));
t1026 = t990 * t996;
t1025 = t992 * t996;
t1001 = cos(qJ(1));
t1002 = qJD(1) ^ 2;
t1000 = cos(qJ(2));
t1020 = t1000 * t992;
t1022 = qJD(1) * t990;
t1017 = t996 * t1022;
t985 = qJD(1) * t992 + qJD(2);
t995 = sin(qJ(3));
t999 = cos(qJ(3));
t962 = -t1017 * t995 + t985 * t999;
t963 = t1017 * t999 + t985 * t995;
t994 = sin(qJ(4));
t948 = t1029 * t963 + t994 * t962;
t1019 = qJD(1) * t1000;
t1016 = t990 * t1019;
t979 = qJD(3) - t1016;
t977 = qJD(4) + t979;
t989 = sin(pkin(12));
t991 = cos(pkin(12));
t931 = -t948 * t989 + t977 * t991;
t947 = -t1029 * t962 + t963 * t994;
t1010 = -mrSges(6,2) * t947 + mrSges(6,3) * t931;
t997 = sin(qJ(1));
t980 = t997 * g(1) - g(2) * t1001;
t970 = qJDD(1) * pkin(1) + t1002 * t1028 + t980;
t1018 = qJDD(1) * t990;
t981 = -g(1) * t1001 - g(2) * t997;
t971 = -pkin(1) * t1002 + pkin(8) * t1018 + t981;
t1023 = t1000 * t971 + t970 * t1025;
t973 = (-pkin(2) * t1000 - pkin(9) * t996) * t1022;
t983 = t985 ^ 2;
t984 = qJDD(1) * t992 + qJDD(2);
t926 = -t983 * pkin(2) + t984 * pkin(9) + (-g(3) * t996 + t1019 * t973) * t990 + t1023;
t974 = (qJD(2) * t1019 + qJDD(1) * t996) * t990;
t975 = -qJD(2) * t1017 + t1000 * t1018;
t927 = -t975 * pkin(2) - t974 * pkin(9) - t1027 + (-t970 + (pkin(2) * t996 - pkin(9) * t1000) * t985 * qJD(1)) * t990;
t895 = -t995 * t926 + t999 * t927;
t943 = qJD(3) * t962 + t974 * t999 + t984 * t995;
t967 = qJDD(3) - t975;
t882 = (t962 * t979 - t943) * pkin(10) + (t962 * t963 + t967) * pkin(3) + t895;
t896 = t999 * t926 + t995 * t927;
t942 = -qJD(3) * t963 - t974 * t995 + t984 * t999;
t952 = pkin(3) * t979 - pkin(10) * t963;
t961 = t962 ^ 2;
t889 = -pkin(3) * t961 + pkin(10) * t942 - t952 * t979 + t896;
t877 = t1029 * t889 + t994 * t882;
t920 = pkin(4) * t947 - qJ(5) * t948;
t966 = qJDD(4) + t967;
t976 = t977 ^ 2;
t871 = -pkin(4) * t976 + qJ(5) * t966 - t920 * t947 + t877;
t1021 = t1000 * t990;
t944 = -g(3) * t1021 + t1020 * t970 - t996 * t971;
t925 = -t984 * pkin(2) - t983 * pkin(9) + t973 * t1017 - t944;
t891 = -t942 * pkin(3) - t961 * pkin(10) + t963 * t952 + t925;
t907 = qJD(4) * t948 - t1029 * t942 + t943 * t994;
t908 = -t947 * qJD(4) + t1029 * t943 + t994 * t942;
t874 = (t947 * t977 - t908) * qJ(5) + (t948 * t977 + t907) * pkin(4) + t891;
t932 = t948 * t991 + t977 * t989;
t866 = -0.2e1 * qJD(5) * t932 - t989 * t871 + t991 * t874;
t898 = t908 * t991 + t966 * t989;
t864 = (t931 * t947 - t898) * pkin(11) + (t931 * t932 + t907) * pkin(5) + t866;
t867 = 0.2e1 * qJD(5) * t931 + t991 * t871 + t989 * t874;
t897 = -t908 * t989 + t966 * t991;
t915 = pkin(5) * t947 - pkin(11) * t932;
t930 = t931 ^ 2;
t865 = -pkin(5) * t930 + pkin(11) * t897 - t915 * t947 + t867;
t993 = sin(qJ(6));
t998 = cos(qJ(6));
t862 = t864 * t998 - t865 * t993;
t911 = t931 * t998 - t932 * t993;
t880 = qJD(6) * t911 + t897 * t993 + t898 * t998;
t912 = t931 * t993 + t932 * t998;
t890 = -mrSges(7,1) * t911 + mrSges(7,2) * t912;
t946 = qJD(6) + t947;
t893 = -mrSges(7,2) * t946 + mrSges(7,3) * t911;
t906 = qJDD(6) + t907;
t858 = m(7) * t862 + mrSges(7,1) * t906 - mrSges(7,3) * t880 - t890 * t912 + t893 * t946;
t863 = t864 * t993 + t865 * t998;
t879 = -qJD(6) * t912 + t897 * t998 - t898 * t993;
t894 = mrSges(7,1) * t946 - mrSges(7,3) * t912;
t859 = m(7) * t863 - mrSges(7,2) * t906 + mrSges(7,3) * t879 + t890 * t911 - t894 * t946;
t850 = t998 * t858 + t993 * t859;
t913 = -mrSges(6,1) * t931 + mrSges(6,2) * t932;
t848 = m(6) * t866 + t907 * mrSges(6,1) - t898 * mrSges(6,3) + t1010 * t947 - t932 * t913 + t850;
t1011 = -t858 * t993 + t998 * t859;
t914 = mrSges(6,1) * t947 - mrSges(6,3) * t932;
t849 = m(6) * t867 - mrSges(6,2) * t907 + mrSges(6,3) * t897 + t913 * t931 - t914 * t947 + t1011;
t1012 = -t848 * t989 + t991 * t849;
t921 = mrSges(5,1) * t947 + mrSges(5,2) * t948;
t934 = mrSges(5,1) * t977 - mrSges(5,3) * t948;
t841 = m(5) * t877 - mrSges(5,2) * t966 - mrSges(5,3) * t907 - t921 * t947 - t934 * t977 + t1012;
t876 = t1029 * t882 - t994 * t889;
t870 = -t966 * pkin(4) - t976 * qJ(5) + t948 * t920 + qJDD(5) - t876;
t868 = -t897 * pkin(5) - t930 * pkin(11) + t932 * t915 + t870;
t1009 = m(7) * t868 - t879 * mrSges(7,1) + mrSges(7,2) * t880 - t911 * t893 + t894 * t912;
t861 = m(6) * t870 - t897 * mrSges(6,1) + mrSges(6,2) * t898 - t931 * t1010 + t914 * t932 + t1009;
t933 = -mrSges(5,2) * t977 - mrSges(5,3) * t947;
t854 = m(5) * t876 + mrSges(5,1) * t966 - mrSges(5,3) * t908 - t921 * t948 + t933 * t977 - t861;
t831 = t1029 * t854 + t994 * t841;
t949 = -mrSges(4,1) * t962 + mrSges(4,2) * t963;
t950 = -mrSges(4,2) * t979 + mrSges(4,3) * t962;
t829 = m(4) * t895 + mrSges(4,1) * t967 - mrSges(4,3) * t943 - t949 * t963 + t950 * t979 + t831;
t1013 = t1029 * t841 - t854 * t994;
t951 = mrSges(4,1) * t979 - mrSges(4,3) * t963;
t830 = m(4) * t896 - mrSges(4,2) * t967 + mrSges(4,3) * t942 + t949 * t962 - t951 * t979 + t1013;
t1014 = -t829 * t995 + t999 * t830;
t945 = -g(3) * t1026 + t1023;
t968 = mrSges(3,1) * t985 - mrSges(3,3) * t1017;
t972 = (-mrSges(3,1) * t1000 + mrSges(3,2) * t996) * t1022;
t821 = m(3) * t945 - mrSges(3,2) * t984 + mrSges(3,3) * t975 + t1016 * t972 - t968 * t985 + t1014;
t824 = t999 * t829 + t995 * t830;
t956 = -t990 * t970 - t1027;
t969 = -mrSges(3,2) * t985 + mrSges(3,3) * t1016;
t823 = m(3) * t956 - t975 * mrSges(3,1) + t974 * mrSges(3,2) + (-t1000 * t969 + t968 * t996) * t1022 + t824;
t843 = t991 * t848 + t989 * t849;
t1007 = m(5) * t891 + t907 * mrSges(5,1) + mrSges(5,2) * t908 + t947 * t933 + t934 * t948 + t843;
t1004 = -m(4) * t925 + t942 * mrSges(4,1) - mrSges(4,2) * t943 + t962 * t950 - t951 * t963 - t1007;
t838 = m(3) * t944 + mrSges(3,1) * t984 - mrSges(3,3) * t974 - t1017 * t972 + t969 * t985 + t1004;
t811 = t838 * t1020 + t821 * t1025 - t823 * t990;
t808 = m(2) * t980 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1002 + t811;
t816 = t1000 * t821 - t838 * t996;
t814 = m(2) * t981 - mrSges(2,1) * t1002 - qJDD(1) * mrSges(2,2) + t816;
t1024 = t1001 * t808 + t997 * t814;
t810 = t838 * t1021 + t821 * t1026 + t992 * t823;
t1015 = t1001 * t814 - t808 * t997;
t883 = Ifges(7,5) * t912 + Ifges(7,6) * t911 + Ifges(7,3) * t946;
t885 = Ifges(7,1) * t912 + Ifges(7,4) * t911 + Ifges(7,5) * t946;
t851 = -mrSges(7,1) * t868 + mrSges(7,3) * t863 + Ifges(7,4) * t880 + Ifges(7,2) * t879 + Ifges(7,6) * t906 - t883 * t912 + t885 * t946;
t884 = Ifges(7,4) * t912 + Ifges(7,2) * t911 + Ifges(7,6) * t946;
t852 = mrSges(7,2) * t868 - mrSges(7,3) * t862 + Ifges(7,1) * t880 + Ifges(7,4) * t879 + Ifges(7,5) * t906 + t883 * t911 - t884 * t946;
t899 = Ifges(6,5) * t932 + Ifges(6,6) * t931 + Ifges(6,3) * t947;
t901 = Ifges(6,1) * t932 + Ifges(6,4) * t931 + Ifges(6,5) * t947;
t833 = -mrSges(6,1) * t870 + mrSges(6,3) * t867 + Ifges(6,4) * t898 + Ifges(6,2) * t897 + Ifges(6,6) * t907 - pkin(5) * t1009 + pkin(11) * t1011 + t998 * t851 + t993 * t852 - t932 * t899 + t947 * t901;
t900 = Ifges(6,4) * t932 + Ifges(6,2) * t931 + Ifges(6,6) * t947;
t835 = mrSges(6,2) * t870 - mrSges(6,3) * t866 + Ifges(6,1) * t898 + Ifges(6,4) * t897 + Ifges(6,5) * t907 - pkin(11) * t850 - t851 * t993 + t852 * t998 + t899 * t931 - t900 * t947;
t916 = Ifges(5,5) * t948 - Ifges(5,6) * t947 + Ifges(5,3) * t977;
t917 = Ifges(5,4) * t948 - Ifges(5,2) * t947 + Ifges(5,6) * t977;
t817 = mrSges(5,2) * t891 - mrSges(5,3) * t876 + Ifges(5,1) * t908 - Ifges(5,4) * t907 + Ifges(5,5) * t966 - qJ(5) * t843 - t833 * t989 + t835 * t991 - t916 * t947 - t917 * t977;
t1005 = mrSges(7,1) * t862 - mrSges(7,2) * t863 + Ifges(7,5) * t880 + Ifges(7,6) * t879 + Ifges(7,3) * t906 + t912 * t884 - t911 * t885;
t918 = Ifges(5,1) * t948 - Ifges(5,4) * t947 + Ifges(5,5) * t977;
t825 = -t1005 - pkin(4) * t843 + t977 * t918 + Ifges(5,6) * t966 + Ifges(5,4) * t908 - Ifges(6,6) * t897 - Ifges(6,5) * t898 - mrSges(6,1) * t866 + mrSges(6,2) * t867 + mrSges(5,3) * t877 - t948 * t916 - pkin(5) * t850 + t931 * t901 - t932 * t900 - mrSges(5,1) * t891 + (-Ifges(5,2) - Ifges(6,3)) * t907;
t936 = Ifges(4,5) * t963 + Ifges(4,6) * t962 + Ifges(4,3) * t979;
t938 = Ifges(4,1) * t963 + Ifges(4,4) * t962 + Ifges(4,5) * t979;
t803 = -mrSges(4,1) * t925 + mrSges(4,3) * t896 + Ifges(4,4) * t943 + Ifges(4,2) * t942 + Ifges(4,6) * t967 - pkin(3) * t1007 + pkin(10) * t1013 + t1029 * t825 + t994 * t817 - t963 * t936 + t979 * t938;
t937 = Ifges(4,4) * t963 + Ifges(4,2) * t962 + Ifges(4,6) * t979;
t806 = mrSges(4,2) * t925 - mrSges(4,3) * t895 + Ifges(4,1) * t943 + Ifges(4,4) * t942 + Ifges(4,5) * t967 - pkin(10) * t831 + t1029 * t817 - t994 * t825 + t962 * t936 - t979 * t937;
t954 = Ifges(3,6) * t985 + (Ifges(3,4) * t996 + Ifges(3,2) * t1000) * t1022;
t955 = Ifges(3,5) * t985 + (Ifges(3,1) * t996 + Ifges(3,4) * t1000) * t1022;
t800 = Ifges(3,5) * t974 + Ifges(3,6) * t975 + Ifges(3,3) * t984 + mrSges(3,1) * t944 - mrSges(3,2) * t945 + t995 * t806 + t999 * t803 + pkin(2) * t1004 + pkin(9) * t1014 + (-t1000 * t955 + t954 * t996) * t1022;
t953 = Ifges(3,3) * t985 + (Ifges(3,5) * t996 + Ifges(3,6) * t1000) * t1022;
t802 = mrSges(3,2) * t956 - mrSges(3,3) * t944 + Ifges(3,1) * t974 + Ifges(3,4) * t975 + Ifges(3,5) * t984 - pkin(9) * t824 + t1016 * t953 - t803 * t995 + t806 * t999 - t954 * t985;
t1006 = -mrSges(5,1) * t876 + mrSges(5,2) * t877 - Ifges(5,5) * t908 + Ifges(5,6) * t907 - Ifges(5,3) * t966 + pkin(4) * t861 - qJ(5) * t1012 - t991 * t833 - t989 * t835 - t948 * t917 - t947 * t918;
t1003 = mrSges(4,1) * t895 - mrSges(4,2) * t896 + Ifges(4,5) * t943 + Ifges(4,6) * t942 + Ifges(4,3) * t967 + pkin(3) * t831 + t963 * t937 - t962 * t938 - t1006;
t805 = -mrSges(3,1) * t956 + mrSges(3,3) * t945 + Ifges(3,4) * t974 + Ifges(3,2) * t975 + Ifges(3,6) * t984 - pkin(2) * t824 - t1017 * t953 + t985 * t955 - t1003;
t1008 = mrSges(2,1) * t980 - mrSges(2,2) * t981 + Ifges(2,3) * qJDD(1) + pkin(1) * t811 + t805 * t1021 + t802 * t1026 + t816 * t1028 + t992 * t800;
t798 = -mrSges(2,2) * g(3) - mrSges(2,3) * t980 + Ifges(2,5) * qJDD(1) - t1002 * Ifges(2,6) + t1000 * t802 - t996 * t805 + (-t810 * t990 - t811 * t992) * pkin(8);
t797 = mrSges(2,1) * g(3) + mrSges(2,3) * t981 + t1002 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t810 - t990 * t800 + (pkin(8) * t816 + t1000 * t805 + t802 * t996) * t992;
t1 = [-m(1) * g(1) + t1015; -m(1) * g(2) + t1024; (-m(1) - m(2)) * g(3) + t810; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1024 + t1001 * t798 - t997 * t797; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1015 + t1001 * t797 + t997 * t798; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1008; t1008; t800; t1003; -t1006; t861; t1005;];
tauJB  = t1;
