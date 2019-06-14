% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 07:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:00:08
% EndTime: 2019-05-07 07:00:33
% DurationCPUTime: 23.86s
% Computational Cost: add. (399654->381), mult. (872689->476), div. (0->0), fcn. (676558->12), ass. (0->157)
t1051 = Ifges(4,1) + Ifges(5,2);
t1046 = Ifges(4,5) - Ifges(5,4);
t1050 = -Ifges(4,2) - Ifges(5,3);
t1045 = Ifges(4,6) - Ifges(5,5);
t1044 = -Ifges(5,6) - Ifges(4,4);
t1049 = Ifges(4,3) + Ifges(5,1);
t1000 = sin(pkin(6));
t1005 = sin(qJ(2));
t1008 = cos(qJ(2));
t983 = (qJD(1) * qJD(2) * t1005 - qJDD(1) * t1008) * t1000;
t1004 = sin(qJ(3));
t1047 = cos(qJ(3));
t1031 = qJD(1) * t1008;
t1002 = cos(pkin(6));
t1028 = t1002 * t1005;
t1010 = qJD(1) ^ 2;
t1042 = pkin(8) * t1000;
t1006 = sin(qJ(1));
t1009 = cos(qJ(1));
t991 = t1006 * g(1) - g(2) * t1009;
t978 = qJDD(1) * pkin(1) + t1010 * t1042 + t991;
t992 = -g(1) * t1009 - g(2) * t1006;
t979 = -pkin(1) * t1010 + qJDD(1) * t1042 + t992;
t1034 = t1008 * t979 + t978 * t1028;
t1032 = qJD(1) * t1000;
t981 = (-pkin(2) * t1008 - pkin(9) * t1005) * t1032;
t995 = qJD(1) * t1002 + qJD(2);
t993 = t995 ^ 2;
t994 = qJDD(1) * t1002 + qJDD(2);
t913 = -t993 * pkin(2) + t994 * pkin(9) + (-g(3) * t1005 + t981 * t1031) * t1000 + t1034;
t1041 = t1002 * g(3);
t982 = (qJD(2) * t1031 + qJDD(1) * t1005) * t1000;
t914 = t983 * pkin(2) - t982 * pkin(9) - t1041 + (-t978 + (pkin(2) * t1005 - pkin(9) * t1008) * t995 * qJD(1)) * t1000;
t897 = t1004 * t914 + t1047 * t913;
t1030 = t1000 * t1005;
t1026 = qJD(1) * t1030;
t969 = t1004 * t1026 - t1047 * t995;
t970 = t1004 * t995 + t1047 * t1026;
t944 = pkin(3) * t969 - qJ(4) * t970;
t975 = qJDD(3) + t983;
t1029 = t1000 * t1008;
t1025 = qJD(1) * t1029;
t989 = -qJD(3) + t1025;
t988 = t989 ^ 2;
t887 = pkin(3) * t988 - t975 * qJ(4) + 0.2e1 * qJD(4) * t989 + t969 * t944 - t897;
t1001 = cos(pkin(11));
t940 = qJD(3) * t970 + t1004 * t982 - t1047 * t994;
t956 = pkin(4) * t970 + qJ(5) * t989;
t968 = t969 ^ 2;
t885 = -pkin(4) * t940 - qJ(5) * t968 - t989 * t956 + qJDD(5) - t887;
t999 = sin(pkin(11));
t918 = t1001 * t940 - t975 * t999;
t953 = -t1001 * t989 + t969 * t999;
t926 = pkin(5) * t970 - pkin(10) * t953;
t952 = t1001 * t969 + t989 * t999;
t951 = t952 ^ 2;
t879 = -pkin(5) * t918 - pkin(10) * t951 + t926 * t953 + t885;
t1003 = sin(qJ(6));
t1007 = cos(qJ(6));
t916 = t1003 * t952 + t1007 * t953;
t919 = t1001 * t975 + t940 * t999;
t894 = -qJD(6) * t916 - t1003 * t919 + t1007 * t918;
t915 = -t1003 * t953 + t1007 * t952;
t895 = qJD(6) * t915 + t1003 * t918 + t1007 * t919;
t967 = qJD(6) + t970;
t903 = -mrSges(7,2) * t967 + mrSges(7,3) * t915;
t904 = mrSges(7,1) * t967 - mrSges(7,3) * t916;
t1018 = m(7) * t879 - t894 * mrSges(7,1) + t895 * mrSges(7,2) - t915 * t903 + t916 * t904;
t924 = -mrSges(6,2) * t970 + mrSges(6,3) * t952;
t925 = mrSges(6,1) * t970 - mrSges(6,3) * t953;
t870 = m(6) * t885 - t918 * mrSges(6,1) + t919 * mrSges(6,2) - t952 * t924 + t953 * t925 + t1018;
t958 = mrSges(5,1) * t970 - mrSges(5,2) * t989;
t1012 = -m(5) * t887 + t975 * mrSges(5,3) - t989 * t958 + t870;
t1035 = t1044 * t969 - t1046 * t989 + t1051 * t970;
t1036 = -t1044 * t970 - t1045 * t989 + t1050 * t969;
t1040 = t969 * t989;
t896 = -t1004 * t913 + t1047 * t914;
t888 = -t975 * pkin(3) - t988 * qJ(4) + t970 * t944 + qJDD(4) - t896;
t941 = -t969 * qJD(3) + t1004 * t994 + t1047 * t982;
t882 = (t969 * t970 - t975) * qJ(5) + (t941 - t1040) * pkin(4) + t888;
t1027 = t1002 * t1008;
t942 = -g(3) * t1029 - t1005 * t979 + t978 * t1027;
t912 = -t994 * pkin(2) - t993 * pkin(9) + t981 * t1026 - t942;
t1013 = (-t941 - t1040) * qJ(4) + t912 + (-t989 * pkin(3) - 0.2e1 * qJD(4)) * t970;
t886 = -t968 * pkin(4) - t970 * t956 + (pkin(3) + qJ(5)) * t940 + t1013;
t876 = -0.2e1 * qJD(5) * t953 + t1001 * t882 - t999 * t886;
t874 = (t952 * t970 - t919) * pkin(10) + (t952 * t953 + t941) * pkin(5) + t876;
t877 = 0.2e1 * qJD(5) * t952 + t1001 * t886 + t999 * t882;
t875 = -pkin(5) * t951 + pkin(10) * t918 - t926 * t970 + t877;
t872 = -t1003 * t875 + t1007 * t874;
t902 = -mrSges(7,1) * t915 + mrSges(7,2) * t916;
t937 = qJDD(6) + t941;
t867 = m(7) * t872 + mrSges(7,1) * t937 - mrSges(7,3) * t895 - t902 * t916 + t903 * t967;
t873 = t1003 * t874 + t1007 * t875;
t868 = m(7) * t873 - mrSges(7,2) * t937 + mrSges(7,3) * t894 + t902 * t915 - t904 * t967;
t1024 = -t1003 * t867 + t1007 * t868;
t898 = Ifges(7,5) * t916 + Ifges(7,6) * t915 + Ifges(7,3) * t967;
t900 = Ifges(7,1) * t916 + Ifges(7,4) * t915 + Ifges(7,5) * t967;
t859 = -mrSges(7,1) * t879 + mrSges(7,3) * t873 + Ifges(7,4) * t895 + Ifges(7,2) * t894 + Ifges(7,6) * t937 - t898 * t916 + t900 * t967;
t899 = Ifges(7,4) * t916 + Ifges(7,2) * t915 + Ifges(7,6) * t967;
t860 = mrSges(7,2) * t879 - mrSges(7,3) * t872 + Ifges(7,1) * t895 + Ifges(7,4) * t894 + Ifges(7,5) * t937 + t898 * t915 - t899 * t967;
t905 = Ifges(6,5) * t953 + Ifges(6,6) * t952 + Ifges(6,3) * t970;
t907 = Ifges(6,1) * t953 + Ifges(6,4) * t952 + Ifges(6,5) * t970;
t844 = -mrSges(6,1) * t885 + mrSges(6,3) * t877 + Ifges(6,4) * t919 + Ifges(6,2) * t918 + Ifges(6,6) * t941 - pkin(5) * t1018 + pkin(10) * t1024 + t1003 * t860 + t1007 * t859 - t953 * t905 + t970 * t907;
t858 = t1003 * t868 + t1007 * t867;
t906 = Ifges(6,4) * t953 + Ifges(6,2) * t952 + Ifges(6,6) * t970;
t845 = mrSges(6,2) * t885 - mrSges(6,3) * t876 + Ifges(6,1) * t919 + Ifges(6,4) * t918 + Ifges(6,5) * t941 - pkin(10) * t858 - t1003 * t859 + t1007 * t860 + t905 * t952 - t906 * t970;
t920 = -mrSges(6,1) * t952 + mrSges(6,2) * t953;
t856 = m(6) * t876 + mrSges(6,1) * t941 - mrSges(6,3) * t919 - t920 * t953 + t924 * t970 + t858;
t857 = m(6) * t877 - mrSges(6,2) * t941 + mrSges(6,3) * t918 + t920 * t952 - t925 * t970 + t1024;
t853 = t1001 * t856 + t999 * t857;
t946 = -mrSges(5,2) * t969 - mrSges(5,3) * t970;
t1016 = -m(5) * t888 - t941 * mrSges(5,1) - t970 * t946 - t853;
t957 = mrSges(5,1) * t969 + mrSges(5,3) * t989;
t852 = t975 * mrSges(5,2) - t989 * t957 - t1016;
t1048 = t1035 * t969 + t1036 * t970 + t1049 * t975 - t1045 * t940 + t1046 * t941 + mrSges(4,1) * t896 - mrSges(4,2) * t897 + mrSges(5,2) * t888 - mrSges(5,3) * t887 - pkin(3) * t852 + qJ(4) * (-t940 * mrSges(5,1) - t969 * t946 + t1012) - qJ(5) * t853 + t1001 * t845 - t999 * t844;
t945 = mrSges(4,1) * t969 + mrSges(4,2) * t970;
t954 = mrSges(4,2) * t989 - mrSges(4,3) * t969;
t850 = m(4) * t896 - t941 * mrSges(4,3) - t970 * t945 + (-t954 + t957) * t989 + (mrSges(4,1) - mrSges(5,2)) * t975 + t1016;
t955 = -mrSges(4,1) * t989 - mrSges(4,3) * t970;
t863 = t1012 + t989 * t955 - t975 * mrSges(4,2) + m(4) * t897 + (-t945 - t946) * t969 + (-mrSges(4,3) - mrSges(5,1)) * t940;
t1023 = -t1004 * t850 + t1047 * t863;
t943 = -g(3) * t1030 + t1034;
t976 = mrSges(3,1) * t995 - mrSges(3,3) * t1026;
t980 = (-mrSges(3,1) * t1008 + mrSges(3,2) * t1005) * t1032;
t840 = m(3) * t943 - mrSges(3,2) * t994 - mrSges(3,3) * t983 + t980 * t1025 - t976 * t995 + t1023;
t843 = t1004 * t863 + t1047 * t850;
t963 = -t1000 * t978 - t1041;
t977 = -mrSges(3,2) * t995 + mrSges(3,3) * t1025;
t842 = m(3) * t963 + t983 * mrSges(3,1) + t982 * mrSges(3,2) + (t1005 * t976 - t1008 * t977) * t1032 + t843;
t1038 = t1001 * t857 - t999 * t856;
t889 = t940 * pkin(3) + t1013;
t1021 = -m(5) * t889 + t940 * mrSges(5,2) + t969 * t957 - t1038;
t1014 = -m(4) * t912 - t940 * mrSges(4,1) - t969 * t954 + (-t955 + t958) * t970 + (-mrSges(4,2) + mrSges(5,3)) * t941 + t1021;
t848 = m(3) * t942 + t994 * mrSges(3,1) - t982 * mrSges(3,3) - t980 * t1026 + t995 * t977 + t1014;
t830 = -t1000 * t842 + t848 * t1027 + t840 * t1028;
t827 = m(2) * t991 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1010 + t830;
t836 = -t1005 * t848 + t1008 * t840;
t834 = m(2) * t992 - mrSges(2,1) * t1010 - qJDD(1) * mrSges(2,2) + t836;
t1039 = t1006 * t834 + t1009 * t827;
t1037 = t1045 * t969 - t1046 * t970 + t1049 * t989;
t829 = t1002 * t842 + t848 * t1029 + t840 * t1030;
t1022 = -t1006 * t827 + t1009 * t834;
t851 = -t941 * mrSges(5,3) - t970 * t958 - t1021;
t825 = -mrSges(4,1) * t912 - mrSges(5,1) * t887 + mrSges(5,2) * t889 + mrSges(4,3) * t897 - pkin(3) * t851 + pkin(4) * t870 - qJ(5) * t1038 - t1001 * t844 - t1035 * t989 + t1037 * t970 - t1044 * t941 + t1045 * t975 + t1050 * t940 - t999 * t845;
t1015 = mrSges(7,1) * t872 - mrSges(7,2) * t873 + Ifges(7,5) * t895 + Ifges(7,6) * t894 + Ifges(7,3) * t937 + t916 * t899 - t915 * t900;
t831 = (Ifges(6,3) + t1051) * t941 + t1046 * t975 + t1037 * t969 + t1044 * t940 + t1015 + t1036 * t989 - t952 * t907 + t953 * t906 + Ifges(6,6) * t918 + Ifges(6,5) * t919 + mrSges(4,2) * t912 - mrSges(4,3) * t896 - mrSges(5,3) * t889 + mrSges(5,1) * t888 - mrSges(6,2) * t877 + mrSges(6,1) * t876 + pkin(5) * t858 + pkin(4) * t853 - qJ(4) * t851;
t961 = Ifges(3,6) * t995 + (Ifges(3,4) * t1005 + Ifges(3,2) * t1008) * t1032;
t962 = Ifges(3,5) * t995 + (Ifges(3,1) * t1005 + Ifges(3,4) * t1008) * t1032;
t820 = Ifges(3,5) * t982 - Ifges(3,6) * t983 + Ifges(3,3) * t994 + mrSges(3,1) * t942 - mrSges(3,2) * t943 + t1004 * t831 + t1047 * t825 + pkin(2) * t1014 + pkin(9) * t1023 + (t1005 * t961 - t1008 * t962) * t1032;
t960 = Ifges(3,3) * t995 + (Ifges(3,5) * t1005 + Ifges(3,6) * t1008) * t1032;
t822 = mrSges(3,2) * t963 - mrSges(3,3) * t942 + Ifges(3,1) * t982 - Ifges(3,4) * t983 + Ifges(3,5) * t994 - pkin(9) * t843 - t1004 * t825 + t960 * t1025 + t1047 * t831 - t995 * t961;
t824 = -mrSges(3,1) * t963 + mrSges(3,3) * t943 + Ifges(3,4) * t982 - Ifges(3,2) * t983 + Ifges(3,6) * t994 - pkin(2) * t843 - t960 * t1026 + t995 * t962 - t1048;
t1017 = mrSges(2,1) * t991 - mrSges(2,2) * t992 + Ifges(2,3) * qJDD(1) + pkin(1) * t830 + t1002 * t820 + t824 * t1029 + t822 * t1030 + t836 * t1042;
t818 = -mrSges(2,2) * g(3) - mrSges(2,3) * t991 + Ifges(2,5) * qJDD(1) - t1010 * Ifges(2,6) - t1005 * t824 + t1008 * t822 + (-t1000 * t829 - t1002 * t830) * pkin(8);
t817 = mrSges(2,1) * g(3) + mrSges(2,3) * t992 + t1010 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t829 - t1000 * t820 + (pkin(8) * t836 + t1005 * t822 + t1008 * t824) * t1002;
t1 = [-m(1) * g(1) + t1022; -m(1) * g(2) + t1039; (-m(1) - m(2)) * g(3) + t829; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1039 - t1006 * t817 + t1009 * t818; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1022 + t1006 * t818 + t1009 * t817; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1017; t1017; t820; t1048; t852; t870; t1015;];
tauJB  = t1;
