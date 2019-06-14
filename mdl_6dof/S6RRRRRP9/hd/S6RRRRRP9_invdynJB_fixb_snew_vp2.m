% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP9
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
% Datum: 2019-05-08 06:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:00:26
% EndTime: 2019-05-08 06:01:13
% DurationCPUTime: 36.40s
% Computational Cost: add. (611777->379), mult. (1301353->476), div. (0->0), fcn. (1043260->12), ass. (0->158)
t1037 = Ifges(6,4) + Ifges(7,4);
t1047 = Ifges(6,2) + Ifges(7,2);
t1043 = Ifges(6,6) + Ifges(7,6);
t1044 = Ifges(6,5) + Ifges(7,5);
t1045 = Ifges(6,1) + Ifges(7,1);
t1000 = sin(qJ(5));
t1005 = cos(qJ(5));
t1001 = sin(qJ(4));
t1006 = cos(qJ(4));
t1002 = sin(qJ(3));
t1007 = cos(qJ(3));
t1003 = sin(qJ(2));
t1027 = qJD(1) * t1003;
t998 = sin(pkin(6));
t1023 = t998 * t1027;
t999 = cos(pkin(6));
t994 = qJD(1) * t999 + qJD(2);
t973 = t1002 * t994 + t1007 * t1023;
t1008 = cos(qJ(2));
t1026 = qJD(1) * t1008;
t1022 = t998 * t1026;
t989 = qJD(3) - t1022;
t959 = -t1001 * t973 + t1006 * t989;
t960 = t1001 * t989 + t1006 * t973;
t935 = -t1000 * t960 + t1005 * t959;
t936 = t1000 * t959 + t1005 * t960;
t972 = -t1002 * t1023 + t1007 * t994;
t971 = qJD(4) - t972;
t969 = qJD(5) + t971;
t1034 = -t1037 * t935 - t1044 * t969 - t1045 * t936;
t1041 = t1037 * t936 + t1043 * t969 + t1047 * t935;
t1042 = Ifges(6,3) + Ifges(7,3);
t1030 = t1003 * t999;
t1010 = qJD(1) ^ 2;
t1039 = pkin(8) * t998;
t1004 = sin(qJ(1));
t1009 = cos(qJ(1));
t990 = t1004 * g(1) - g(2) * t1009;
t980 = qJDD(1) * pkin(1) + t1010 * t1039 + t990;
t991 = -g(1) * t1009 - g(2) * t1004;
t981 = -pkin(1) * t1010 + qJDD(1) * t1039 + t991;
t1033 = t1008 * t981 + t980 * t1030;
t1032 = qJD(1) * t998;
t983 = (-pkin(2) * t1008 - pkin(9) * t1003) * t1032;
t992 = t994 ^ 2;
t993 = qJDD(1) * t999 + qJDD(2);
t932 = -t992 * pkin(2) + t993 * pkin(9) + (-g(3) * t1003 + t1026 * t983) * t998 + t1033;
t1038 = t999 * g(3);
t984 = (qJD(2) * t1026 + qJDD(1) * t1003) * t998;
t985 = (qJD(2) * t1027 - qJDD(1) * t1008) * t998;
t933 = t985 * pkin(2) - t984 * pkin(9) - t1038 + (-t980 + (pkin(2) * t1003 - pkin(9) * t1008) * t994 * qJD(1)) * t998;
t902 = t1002 * t933 + t1007 * t932;
t957 = -pkin(3) * t972 - pkin(10) * t973;
t977 = qJDD(3) + t985;
t987 = t989 ^ 2;
t894 = -pkin(3) * t987 + pkin(10) * t977 + t957 * t972 + t902;
t1028 = t1008 * t999;
t1029 = t1008 * t998;
t954 = -g(3) * t1029 - t1003 * t981 + t1028 * t980;
t931 = -t993 * pkin(2) - t992 * pkin(9) + t983 * t1023 - t954;
t952 = -t973 * qJD(3) - t1002 * t984 + t1007 * t993;
t953 = qJD(3) * t972 + t1002 * t993 + t1007 * t984;
t899 = (-t972 * t989 - t953) * pkin(10) + (t973 * t989 - t952) * pkin(3) + t931;
t879 = -t1001 * t894 + t1006 * t899;
t917 = qJD(4) * t959 + t1001 * t977 + t1006 * t953;
t950 = qJDD(4) - t952;
t876 = (t959 * t971 - t917) * pkin(11) + (t959 * t960 + t950) * pkin(4) + t879;
t880 = t1001 * t899 + t1006 * t894;
t916 = -qJD(4) * t960 - t1001 * t953 + t1006 * t977;
t941 = pkin(4) * t971 - pkin(11) * t960;
t958 = t959 ^ 2;
t878 = -pkin(4) * t958 + pkin(11) * t916 - t941 * t971 + t880;
t870 = -t1000 * t878 + t1005 * t876;
t891 = qJD(5) * t935 + t1000 * t916 + t1005 * t917;
t945 = qJDD(5) + t950;
t865 = -0.2e1 * qJD(6) * t936 + (t935 * t969 - t891) * qJ(6) + (t935 * t936 + t945) * pkin(5) + t870;
t918 = -mrSges(7,2) * t969 + mrSges(7,3) * t935;
t1025 = m(7) * t865 + t945 * mrSges(7,1) + t969 * t918;
t912 = -mrSges(7,1) * t935 + mrSges(7,2) * t936;
t862 = -t891 * mrSges(7,3) - t936 * t912 + t1025;
t871 = t1000 * t876 + t1005 * t878;
t890 = -qJD(5) * t936 - t1000 * t917 + t1005 * t916;
t920 = pkin(5) * t969 - qJ(6) * t936;
t934 = t935 ^ 2;
t867 = -pkin(5) * t934 + qJ(6) * t890 + 0.2e1 * qJD(6) * t935 - t920 * t969 + t871;
t1046 = mrSges(6,1) * t870 + mrSges(7,1) * t865 - mrSges(6,2) * t871 - mrSges(7,2) * t867 + pkin(5) * t862 + t1034 * t935 + t1041 * t936 + t1042 * t945 + t1043 * t890 + t1044 * t891;
t913 = -mrSges(6,1) * t935 + mrSges(6,2) * t936;
t919 = -mrSges(6,2) * t969 + mrSges(6,3) * t935;
t854 = m(6) * t870 + t945 * mrSges(6,1) + t969 * t919 + (-t912 - t913) * t936 + (-mrSges(6,3) - mrSges(7,3)) * t891 + t1025;
t1024 = m(7) * t867 + t890 * mrSges(7,3) + t935 * t912;
t921 = mrSges(7,1) * t969 - mrSges(7,3) * t936;
t922 = mrSges(6,1) * t969 - mrSges(6,3) * t936;
t857 = m(6) * t871 + t890 * mrSges(6,3) + t935 * t913 + (-t921 - t922) * t969 + (-mrSges(6,2) - mrSges(7,2)) * t945 + t1024;
t852 = t1000 * t857 + t1005 * t854;
t924 = Ifges(5,4) * t960 + Ifges(5,2) * t959 + Ifges(5,6) * t971;
t925 = Ifges(5,1) * t960 + Ifges(5,4) * t959 + Ifges(5,5) * t971;
t1040 = mrSges(5,1) * t879 - mrSges(5,2) * t880 + Ifges(5,5) * t917 + Ifges(5,6) * t916 + Ifges(5,3) * t950 + pkin(4) * t852 + t960 * t924 - t959 * t925 + t1046;
t937 = -mrSges(5,1) * t959 + mrSges(5,2) * t960;
t939 = -mrSges(5,2) * t971 + mrSges(5,3) * t959;
t849 = m(5) * t879 + mrSges(5,1) * t950 - mrSges(5,3) * t917 - t937 * t960 + t939 * t971 + t852;
t1020 = -t1000 * t854 + t1005 * t857;
t940 = mrSges(5,1) * t971 - mrSges(5,3) * t960;
t850 = m(5) * t880 - mrSges(5,2) * t950 + mrSges(5,3) * t916 + t937 * t959 - t940 * t971 + t1020;
t846 = -t1001 * t849 + t1006 * t850;
t956 = -mrSges(4,1) * t972 + mrSges(4,2) * t973;
t962 = mrSges(4,1) * t989 - mrSges(4,3) * t973;
t844 = m(4) * t902 - mrSges(4,2) * t977 + mrSges(4,3) * t952 + t956 * t972 - t962 * t989 + t846;
t901 = -t1002 * t932 + t1007 * t933;
t893 = -pkin(3) * t977 - pkin(10) * t987 + t973 * t957 - t901;
t881 = -pkin(4) * t916 - pkin(11) * t958 + t960 * t941 + t893;
t873 = -pkin(5) * t890 - qJ(6) * t934 + t920 * t936 + qJDD(6) + t881;
t868 = m(7) * t873 - t890 * mrSges(7,1) + t891 * mrSges(7,2) - t935 * t918 + t936 * t921;
t1015 = m(6) * t881 - t890 * mrSges(6,1) + mrSges(6,2) * t891 - t935 * t919 + t922 * t936 + t868;
t860 = -m(5) * t893 + t916 * mrSges(5,1) - mrSges(5,2) * t917 + t959 * t939 - t940 * t960 - t1015;
t961 = -mrSges(4,2) * t989 + mrSges(4,3) * t972;
t859 = m(4) * t901 + mrSges(4,1) * t977 - mrSges(4,3) * t953 - t956 * t973 + t961 * t989 + t860;
t1019 = -t1002 * t859 + t1007 * t844;
t1031 = t1003 * t998;
t955 = -g(3) * t1031 + t1033;
t978 = mrSges(3,1) * t994 - mrSges(3,3) * t1023;
t982 = (-mrSges(3,1) * t1008 + mrSges(3,2) * t1003) * t1032;
t835 = m(3) * t955 - mrSges(3,2) * t993 - mrSges(3,3) * t985 + t1022 * t982 - t978 * t994 + t1019;
t838 = t1002 * t844 + t1007 * t859;
t966 = -t998 * t980 - t1038;
t979 = -mrSges(3,2) * t994 + mrSges(3,3) * t1022;
t837 = m(3) * t966 + t985 * mrSges(3,1) + t984 * mrSges(3,2) + (t1003 * t978 - t1008 * t979) * t1032 + t838;
t845 = t1001 * t850 + t1006 * t849;
t1014 = -m(4) * t931 + t952 * mrSges(4,1) - mrSges(4,2) * t953 + t972 * t961 - t962 * t973 - t845;
t841 = m(3) * t954 + mrSges(3,1) * t993 - mrSges(3,3) * t984 - t1023 * t982 + t979 * t994 + t1014;
t823 = t841 * t1028 + t835 * t1030 - t837 * t998;
t820 = m(2) * t990 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1010 + t823;
t828 = -t1003 * t841 + t1008 * t835;
t826 = m(2) * t991 - mrSges(2,1) * t1010 - qJDD(1) * mrSges(2,2) + t828;
t1036 = t1004 * t826 + t1009 * t820;
t1035 = -t1042 * t969 - t1043 * t935 - t1044 * t936;
t822 = t841 * t1029 + t835 * t1031 + t999 * t837;
t1018 = -t1004 * t820 + t1009 * t826;
t847 = -mrSges(6,1) * t881 + mrSges(6,3) * t871 - mrSges(7,1) * t873 + mrSges(7,3) * t867 - pkin(5) * t868 + qJ(6) * t1024 + (-qJ(6) * t921 - t1034) * t969 + (-mrSges(7,2) * qJ(6) + t1043) * t945 + t1035 * t936 + t1037 * t891 + t1047 * t890;
t851 = mrSges(6,2) * t881 + mrSges(7,2) * t873 - mrSges(6,3) * t870 - mrSges(7,3) * t865 - qJ(6) * t862 - t1035 * t935 + t1037 * t890 - t1041 * t969 + t1044 * t945 + t1045 * t891;
t923 = Ifges(5,5) * t960 + Ifges(5,6) * t959 + Ifges(5,3) * t971;
t830 = -mrSges(5,1) * t893 + mrSges(5,3) * t880 + Ifges(5,4) * t917 + Ifges(5,2) * t916 + Ifges(5,6) * t950 - pkin(4) * t1015 + pkin(11) * t1020 + t1000 * t851 + t1005 * t847 - t960 * t923 + t971 * t925;
t831 = mrSges(5,2) * t893 - mrSges(5,3) * t879 + Ifges(5,1) * t917 + Ifges(5,4) * t916 + Ifges(5,5) * t950 - pkin(11) * t852 - t1000 * t847 + t1005 * t851 + t923 * t959 - t924 * t971;
t946 = Ifges(4,5) * t973 + Ifges(4,6) * t972 + Ifges(4,3) * t989;
t947 = Ifges(4,4) * t973 + Ifges(4,2) * t972 + Ifges(4,6) * t989;
t818 = mrSges(4,2) * t931 - mrSges(4,3) * t901 + Ifges(4,1) * t953 + Ifges(4,4) * t952 + Ifges(4,5) * t977 - pkin(10) * t845 - t1001 * t830 + t1006 * t831 + t946 * t972 - t947 * t989;
t948 = Ifges(4,1) * t973 + Ifges(4,4) * t972 + Ifges(4,5) * t989;
t829 = -mrSges(4,1) * t931 + mrSges(4,3) * t902 + Ifges(4,4) * t953 + Ifges(4,2) * t952 + Ifges(4,6) * t977 - pkin(3) * t845 - t973 * t946 + t989 * t948 - t1040;
t964 = Ifges(3,6) * t994 + (Ifges(3,4) * t1003 + Ifges(3,2) * t1008) * t1032;
t965 = Ifges(3,5) * t994 + (Ifges(3,1) * t1003 + Ifges(3,4) * t1008) * t1032;
t813 = Ifges(3,5) * t984 - Ifges(3,6) * t985 + Ifges(3,3) * t993 + mrSges(3,1) * t954 - mrSges(3,2) * t955 + t1002 * t818 + t1007 * t829 + pkin(2) * t1014 + pkin(9) * t1019 + (t1003 * t964 - t1008 * t965) * t1032;
t963 = Ifges(3,3) * t994 + (Ifges(3,5) * t1003 + Ifges(3,6) * t1008) * t1032;
t815 = mrSges(3,2) * t966 - mrSges(3,3) * t954 + Ifges(3,1) * t984 - Ifges(3,4) * t985 + Ifges(3,5) * t993 - pkin(9) * t838 - t1002 * t829 + t1007 * t818 + t1022 * t963 - t964 * t994;
t1012 = mrSges(4,1) * t901 - mrSges(4,2) * t902 + Ifges(4,5) * t953 + Ifges(4,6) * t952 + Ifges(4,3) * t977 + pkin(3) * t860 + pkin(10) * t846 + t1001 * t831 + t1006 * t830 + t973 * t947 - t972 * t948;
t817 = -mrSges(3,1) * t966 + mrSges(3,3) * t955 + Ifges(3,4) * t984 - Ifges(3,2) * t985 + Ifges(3,6) * t993 - pkin(2) * t838 - t1023 * t963 + t994 * t965 - t1012;
t1016 = mrSges(2,1) * t990 - mrSges(2,2) * t991 + Ifges(2,3) * qJDD(1) + pkin(1) * t823 + t817 * t1029 + t815 * t1031 + t1039 * t828 + t999 * t813;
t811 = -mrSges(2,2) * g(3) - mrSges(2,3) * t990 + Ifges(2,5) * qJDD(1) - t1010 * Ifges(2,6) - t1003 * t817 + t1008 * t815 + (-t822 * t998 - t823 * t999) * pkin(8);
t810 = mrSges(2,1) * g(3) + mrSges(2,3) * t991 + t1010 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t822 - t998 * t813 + (pkin(8) * t828 + t1003 * t815 + t1008 * t817) * t999;
t1 = [-m(1) * g(1) + t1018; -m(1) * g(2) + t1036; (-m(1) - m(2)) * g(3) + t822; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1036 - t1004 * t810 + t1009 * t811; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1018 + t1004 * t811 + t1009 * t810; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1016; t1016; t813; t1012; t1040; t1046; t868;];
tauJB  = t1;
