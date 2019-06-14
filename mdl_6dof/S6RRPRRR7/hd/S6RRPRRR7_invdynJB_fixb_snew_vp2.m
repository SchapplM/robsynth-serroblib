% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-06 22:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:17:23
% EndTime: 2019-05-06 22:17:38
% DurationCPUTime: 11.07s
% Computational Cost: add. (174446->366), mult. (354002->444), div. (0->0), fcn. (231469->10), ass. (0->147)
t1041 = Ifges(3,1) + Ifges(4,1);
t1033 = Ifges(3,4) - Ifges(4,5);
t1032 = (Ifges(3,5) + Ifges(4,4));
t1040 = Ifges(3,2) + Ifges(4,3);
t1031 = (Ifges(3,6) - Ifges(4,6));
t1039 = Ifges(3,3) + Ifges(4,2);
t1000 = cos(qJ(2));
t1003 = qJD(1) ^ 2;
t1024 = t1003 * t1000 ^ 2;
t1002 = qJD(2) ^ 2;
t1023 = qJD(1) * t1000;
t1035 = 2 * qJD(3);
t1001 = cos(qJ(1));
t996 = sin(qJ(1));
t971 = -t1001 * g(1) - t996 * g(2);
t947 = -t1003 * pkin(1) + qJDD(1) * pkin(7) + t971;
t995 = sin(qJ(2));
t927 = -t995 * g(3) + t1000 * t947;
t957 = (-pkin(2) * t1000 - qJ(3) * t995) * qJD(1);
t904 = -t1002 * pkin(2) + qJDD(2) * qJ(3) + (qJD(2) * t1035) + t957 * t1023 + t927;
t1025 = qJD(1) * t995;
t1022 = qJD(2) * t1025;
t961 = t1000 * qJDD(1) - t1022;
t969 = -(qJD(2) * pkin(3)) - pkin(8) * t1025;
t889 = -pkin(3) * t1024 - t961 * pkin(8) + qJD(2) * t969 + t904;
t1021 = qJD(2) * t1023;
t926 = -t1000 * g(3) - t995 * t947;
t910 = -qJDD(2) * pkin(2) - t1002 * qJ(3) + t957 * t1025 + qJDD(3) - t926;
t960 = t995 * qJDD(1) + t1021;
t890 = (-t960 + t1021) * pkin(8) + (-t1000 * t1003 * t995 - qJDD(2)) * pkin(3) + t910;
t994 = sin(qJ(4));
t999 = cos(qJ(4));
t872 = -t994 * t889 + t999 * t890;
t944 = -t999 * t1023 - t994 * t1025;
t945 = (-t1000 * t994 + t995 * t999) * qJD(1);
t922 = -t944 * pkin(4) - t945 * pkin(9);
t984 = -qJD(2) + qJD(4);
t982 = t984 ^ 2;
t983 = -qJDD(2) + qJDD(4);
t868 = -t983 * pkin(4) - t982 * pkin(9) + t945 * t922 - t872;
t912 = t944 * qJD(4) + t999 * t960 - t994 * t961;
t993 = sin(qJ(5));
t998 = cos(qJ(5));
t925 = t998 * t945 + t993 * t984;
t882 = -t925 * qJD(5) - t993 * t912 + t998 * t983;
t937 = qJD(5) - t944;
t915 = t937 * pkin(5) - t925 * pkin(10);
t924 = -t993 * t945 + t998 * t984;
t923 = t924 ^ 2;
t858 = -t882 * pkin(5) - t923 * pkin(10) + t925 * t915 + t868;
t883 = t924 * qJD(5) + t998 * t912 + t993 * t983;
t992 = sin(qJ(6));
t997 = cos(qJ(6));
t897 = t992 * t924 + t997 * t925;
t862 = -t897 * qJD(6) + t997 * t882 - t992 * t883;
t896 = t997 * t924 - t992 * t925;
t863 = t896 * qJD(6) + t992 * t882 + t997 * t883;
t932 = qJD(6) + t937;
t887 = -t932 * mrSges(7,2) + t896 * mrSges(7,3);
t888 = t932 * mrSges(7,1) - t897 * mrSges(7,3);
t1011 = m(7) * t858 - t862 * mrSges(7,1) + t863 * mrSges(7,2) - t896 * t887 + t897 * t888;
t913 = -t937 * mrSges(6,2) + t924 * mrSges(6,3);
t914 = t937 * mrSges(6,1) - t925 * mrSges(6,3);
t1007 = -m(6) * t868 + t882 * mrSges(6,1) - t883 * mrSges(6,2) + t924 * t913 - t925 * t914 - t1011;
t970 = t996 * g(1) - t1001 * g(2);
t946 = -qJDD(1) * pkin(1) - t1003 * pkin(7) - t970;
t1014 = -t961 * pkin(2) + t946 + (-t1021 - t960) * qJ(3);
t880 = -pkin(2) * t1022 + t961 * pkin(3) - pkin(8) * t1024 - t1014 + (t1035 + t969) * t1025;
t911 = -t945 * qJD(4) - t994 * t960 - t999 * t961;
t866 = t880 + (-t944 * t984 - t912) * pkin(9) + (t945 * t984 - t911) * pkin(4);
t873 = t999 * t889 + t994 * t890;
t869 = -t982 * pkin(4) + t983 * pkin(9) + t944 * t922 + t873;
t856 = t998 * t866 - t993 * t869;
t909 = qJDD(5) - t911;
t854 = (t924 * t937 - t883) * pkin(10) + (t924 * t925 + t909) * pkin(5) + t856;
t857 = t993 * t866 + t998 * t869;
t855 = -t923 * pkin(5) + t882 * pkin(10) - t937 * t915 + t857;
t852 = t997 * t854 - t992 * t855;
t878 = -t896 * mrSges(7,1) + t897 * mrSges(7,2);
t902 = qJDD(6) + t909;
t847 = m(7) * t852 + t902 * mrSges(7,1) - t863 * mrSges(7,3) - t897 * t878 + t932 * t887;
t853 = t992 * t854 + t997 * t855;
t848 = m(7) * t853 - t902 * mrSges(7,2) + t862 * mrSges(7,3) + t896 * t878 - t932 * t888;
t840 = t997 * t847 + t992 * t848;
t898 = -t924 * mrSges(6,1) + t925 * mrSges(6,2);
t838 = m(6) * t856 + t909 * mrSges(6,1) - t883 * mrSges(6,3) - t925 * t898 + t937 * t913 + t840;
t1016 = -t992 * t847 + t997 * t848;
t839 = m(6) * t857 - t909 * mrSges(6,2) + t882 * mrSges(6,3) + t924 * t898 - t937 * t914 + t1016;
t1017 = -t993 * t838 + t998 * t839;
t874 = Ifges(7,5) * t897 + Ifges(7,6) * t896 + Ifges(7,3) * t932;
t876 = Ifges(7,1) * t897 + Ifges(7,4) * t896 + Ifges(7,5) * t932;
t841 = -mrSges(7,1) * t858 + mrSges(7,3) * t853 + Ifges(7,4) * t863 + Ifges(7,2) * t862 + Ifges(7,6) * t902 - t897 * t874 + t932 * t876;
t875 = Ifges(7,4) * t897 + Ifges(7,2) * t896 + Ifges(7,6) * t932;
t842 = mrSges(7,2) * t858 - mrSges(7,3) * t852 + Ifges(7,1) * t863 + Ifges(7,4) * t862 + Ifges(7,5) * t902 + t896 * t874 - t932 * t875;
t891 = Ifges(6,5) * t925 + Ifges(6,6) * t924 + Ifges(6,3) * t937;
t893 = Ifges(6,1) * t925 + Ifges(6,4) * t924 + Ifges(6,5) * t937;
t823 = -mrSges(6,1) * t868 + mrSges(6,3) * t857 + Ifges(6,4) * t883 + Ifges(6,2) * t882 + Ifges(6,6) * t909 - pkin(5) * t1011 + pkin(10) * t1016 + t997 * t841 + t992 * t842 - t925 * t891 + t937 * t893;
t892 = Ifges(6,4) * t925 + Ifges(6,2) * t924 + Ifges(6,6) * t937;
t825 = mrSges(6,2) * t868 - mrSges(6,3) * t856 + Ifges(6,1) * t883 + Ifges(6,4) * t882 + Ifges(6,5) * t909 - pkin(10) * t840 - t992 * t841 + t997 * t842 + t924 * t891 - t937 * t892;
t917 = Ifges(5,4) * t945 + Ifges(5,2) * t944 + (Ifges(5,6) * t984);
t918 = Ifges(5,1) * t945 + Ifges(5,4) * t944 + (Ifges(5,5) * t984);
t1008 = -mrSges(5,1) * t872 + mrSges(5,2) * t873 - Ifges(5,5) * t912 - Ifges(5,6) * t911 - Ifges(5,3) * t983 - pkin(4) * t1007 - pkin(9) * t1017 - t998 * t823 - t993 * t825 - t945 * t917 + t944 * t918;
t921 = -t944 * mrSges(5,1) + t945 * mrSges(5,2);
t929 = t984 * mrSges(5,1) - t945 * mrSges(5,3);
t831 = m(5) * t873 - t983 * mrSges(5,2) + t911 * mrSges(5,3) + t944 * t921 - t984 * t929 + t1017;
t928 = -t984 * mrSges(5,2) + t944 * mrSges(5,3);
t843 = m(5) * t872 + t983 * mrSges(5,1) - t912 * mrSges(5,3) - t945 * t921 + t984 * t928 + t1007;
t1018 = t999 * t831 - t994 * t843;
t958 = (-mrSges(4,1) * t1000 - mrSges(4,3) * t995) * qJD(1);
t966 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t1025;
t1013 = m(4) * t904 + qJDD(2) * mrSges(4,3) + qJD(2) * t966 + t958 * t1023 + t1018;
t1026 = (t1032 * qJD(2)) + (t1033 * t1000 + t1041 * t995) * qJD(1);
t1027 = -(t1031 * qJD(2)) + (-t1040 * t1000 - t1033 * t995) * qJD(1);
t821 = t994 * t831 + t999 * t843;
t968 = mrSges(4,2) * t1023 + qJD(2) * mrSges(4,3);
t1010 = -m(4) * t910 + qJDD(2) * mrSges(4,1) + qJD(2) * t968 - t821;
t820 = t960 * mrSges(4,2) + t1025 * t958 - t1010;
t1038 = -(t1000 * t1026 + t1027 * t995) * qJD(1) + t1039 * qJDD(2) + t1031 * t961 + t1032 * t960 + mrSges(3,1) * t926 - mrSges(4,1) * t910 - mrSges(3,2) * t927 + mrSges(4,3) * t904 - pkin(2) * t820 - pkin(3) * t821 + qJ(3) * (t961 * mrSges(4,2) + t1013) + t1008;
t1034 = mrSges(3,3) + mrSges(4,2);
t959 = (-mrSges(3,1) * t1000 + mrSges(3,2) * t995) * qJD(1);
t965 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1025;
t817 = m(3) * t927 - qJDD(2) * mrSges(3,2) - qJD(2) * t965 + t1023 * t959 + t1034 * t961 + t1013;
t967 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1023;
t818 = m(3) * t926 + qJDD(2) * mrSges(3,1) + qJD(2) * t967 - t1034 * t960 + (-t958 - t959) * t1025 + t1010;
t1019 = t1000 * t817 - t995 * t818;
t810 = m(2) * t971 - t1003 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t1019;
t833 = t998 * t838 + t993 * t839;
t1015 = -m(5) * t880 + t911 * mrSges(5,1) - t912 * mrSges(5,2) + t944 * t928 - t945 * t929 - t833;
t895 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t1025 + t1014;
t829 = m(4) * t895 - t961 * mrSges(4,1) - t960 * mrSges(4,3) - t968 * t1023 - t966 * t1025 + t1015;
t1006 = -m(3) * t946 + t961 * mrSges(3,1) - t960 * mrSges(3,2) + t967 * t1023 - t1025 * t965 - t829;
t827 = m(2) * t970 + qJDD(1) * mrSges(2,1) - t1003 * mrSges(2,2) + t1006;
t1029 = t1001 * t827 + t996 * t810;
t812 = t1000 * t818 + t995 * t817;
t1028 = t1039 * qJD(2) + (t1031 * t1000 + t1032 * t995) * qJD(1);
t1020 = t1001 * t810 - t996 * t827;
t916 = Ifges(5,5) * t945 + Ifges(5,6) * t944 + Ifges(5,3) * t984;
t807 = mrSges(5,2) * t880 - mrSges(5,3) * t872 + Ifges(5,1) * t912 + Ifges(5,4) * t911 + Ifges(5,5) * t983 - pkin(9) * t833 - t993 * t823 + t998 * t825 + t944 * t916 - t984 * t917;
t1009 = -mrSges(7,1) * t852 + mrSges(7,2) * t853 - Ifges(7,5) * t863 - Ifges(7,6) * t862 - Ifges(7,3) * t902 - t897 * t875 + t896 * t876;
t1005 = mrSges(6,1) * t856 - mrSges(6,2) * t857 + Ifges(6,5) * t883 + Ifges(6,6) * t882 + Ifges(6,3) * t909 + pkin(5) * t840 + t925 * t892 - t924 * t893 - t1009;
t813 = -mrSges(5,1) * t880 + mrSges(5,3) * t873 + Ifges(5,4) * t912 + Ifges(5,2) * t911 + Ifges(5,6) * t983 - pkin(4) * t833 - t945 * t916 + t984 * t918 - t1005;
t804 = -mrSges(3,1) * t946 - mrSges(4,1) * t895 + mrSges(4,2) * t904 + mrSges(3,3) * t927 - pkin(2) * t829 - pkin(3) * t1015 - pkin(8) * t1018 + t1026 * qJD(2) + t1031 * qJDD(2) - t1028 * t1025 + t1033 * t960 + t1040 * t961 - t994 * t807 - t999 * t813;
t806 = mrSges(3,2) * t946 + mrSges(4,2) * t910 - mrSges(3,3) * t926 - mrSges(4,3) * t895 - pkin(8) * t821 - qJ(3) * t829 + t1027 * qJD(2) + t1032 * qJDD(2) + t1028 * t1023 + t1033 * t961 + t1041 * t960 + t999 * t807 - t994 * t813;
t1012 = mrSges(2,1) * t970 - mrSges(2,2) * t971 + Ifges(2,3) * qJDD(1) + pkin(1) * t1006 + pkin(7) * t1019 + t1000 * t804 + t995 * t806;
t802 = mrSges(2,1) * g(3) + mrSges(2,3) * t971 + t1003 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t812 - t1038;
t801 = -mrSges(2,2) * g(3) - mrSges(2,3) * t970 + Ifges(2,5) * qJDD(1) - t1003 * Ifges(2,6) - pkin(7) * t812 + t1000 * t806 - t995 * t804;
t1 = [-m(1) * g(1) + t1020; -m(1) * g(2) + t1029; (-m(1) - m(2)) * g(3) + t812; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1029 + t1001 * t801 - t996 * t802; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1020 + t1001 * t802 + t996 * t801; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1012; t1012; t1038; t820; -t1008; t1005; -t1009;];
tauJB  = t1;
