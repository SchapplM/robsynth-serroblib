% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:20:48
% EndTime: 2019-05-07 04:21:03
% DurationCPUTime: 14.18s
% Computational Cost: add. (201135->370), mult. (458068->448), div. (0->0), fcn. (331576->10), ass. (0->150)
t1031 = Ifges(5,4) + Ifges(6,6);
t1043 = -Ifges(5,2) - Ifges(6,3);
t1037 = Ifges(5,6) - Ifges(6,5);
t1030 = cos(pkin(10));
t1019 = qJD(1) * qJD(2);
t994 = sin(qJ(1));
t998 = cos(qJ(1));
t979 = -g(1) * t998 - g(2) * t994;
t999 = qJD(1) ^ 2;
t967 = -pkin(1) * t999 + qJDD(1) * pkin(7) + t979;
t993 = sin(qJ(2));
t1028 = t993 * t967;
t1032 = pkin(2) * t999;
t997 = cos(qJ(2));
t972 = qJDD(1) * t993 + t997 * t1019;
t922 = qJDD(2) * pkin(2) - t972 * pkin(8) - t1028 + (pkin(8) * t1019 + t993 * t1032 - g(3)) * t997;
t954 = -g(3) * t993 + t997 * t967;
t973 = qJDD(1) * t997 - t993 * t1019;
t1022 = qJD(1) * t993;
t977 = qJD(2) * pkin(2) - pkin(8) * t1022;
t989 = t997 ^ 2;
t923 = pkin(8) * t973 - qJD(2) * t977 - t989 * t1032 + t954;
t992 = sin(qJ(3));
t996 = cos(qJ(3));
t887 = t996 * t922 - t992 * t923;
t964 = (-t992 * t993 + t996 * t997) * qJD(1);
t932 = qJD(3) * t964 + t972 * t996 + t973 * t992;
t965 = (t992 * t997 + t993 * t996) * qJD(1);
t986 = qJDD(2) + qJDD(3);
t987 = qJD(2) + qJD(3);
t866 = (t964 * t987 - t932) * qJ(4) + (t964 * t965 + t986) * pkin(3) + t887;
t888 = t992 * t922 + t996 * t923;
t931 = -qJD(3) * t965 - t972 * t992 + t973 * t996;
t956 = pkin(3) * t987 - qJ(4) * t965;
t960 = t964 ^ 2;
t869 = -pkin(3) * t960 + qJ(4) * t931 - t956 * t987 + t888;
t990 = sin(pkin(10));
t1025 = t1030 * t869 + t990 * t866;
t985 = t987 ^ 2;
t1010 = t985 * pkin(4) - t986 * qJ(5) - t1025;
t900 = -t1030 * t931 + t932 * t990;
t950 = -t1030 * t964 + t965 * t990;
t951 = t1030 * t965 + t990 * t964;
t915 = pkin(4) * t950 - qJ(5) * t951;
t938 = pkin(5) * t951 - pkin(9) * t987;
t1041 = -2 * qJD(4);
t941 = t950 * t1041;
t947 = t950 ^ 2;
t855 = -t900 * pkin(5) - t947 * pkin(9) - t950 * t915 + t941 + ((2 * qJD(5)) + t938) * t987 - t1010;
t991 = sin(qJ(6));
t995 = cos(qJ(6));
t930 = t950 * t991 + t987 * t995;
t874 = -qJD(6) * t930 + t900 * t995 - t986 * t991;
t929 = t950 * t995 - t987 * t991;
t875 = qJD(6) * t929 + t900 * t991 + t986 * t995;
t943 = qJD(6) + t951;
t903 = -mrSges(7,2) * t943 + mrSges(7,3) * t929;
t904 = mrSges(7,1) * t943 - mrSges(7,3) * t930;
t1009 = -m(7) * t855 + t874 * mrSges(7,1) - t875 * mrSges(7,2) + t929 * t903 - t930 * t904;
t1033 = -2 * qJD(5);
t857 = t987 * t1033 + ((2 * qJD(4)) + t915) * t950 + t1010;
t937 = mrSges(6,1) * t951 + mrSges(6,2) * t987;
t1005 = -m(6) * t857 + t986 * mrSges(6,3) + t987 * t937 - t1009;
t1038 = Ifges(5,5) - Ifges(6,4);
t1040 = Ifges(5,1) + Ifges(6,2);
t1023 = -t1031 * t950 + t1038 * t987 + t1040 * t951;
t1036 = t1031 * t951 + t1037 * t987 + t1043 * t950;
t1039 = -Ifges(6,1) - Ifges(5,3);
t1029 = t950 * t987;
t863 = t1030 * t866 + t1041 * t951 - t990 * t869;
t859 = -t986 * pkin(4) - t985 * qJ(5) + t951 * t915 + qJDD(5) - t863;
t901 = t1030 * t932 + t990 * t931;
t853 = (t950 * t951 - t986) * pkin(9) + (t901 + t1029) * pkin(5) + t859;
t978 = t994 * g(1) - t998 * g(2);
t1011 = -qJDD(1) * pkin(1) - t978;
t933 = -t973 * pkin(2) + t977 * t1022 + (-pkin(8) * t989 - pkin(7)) * t999 + t1011;
t877 = -t931 * pkin(3) - t960 * qJ(4) + t965 * t956 + qJDD(4) + t933;
t1003 = (-t901 + t1029) * qJ(5) + t877 + (t987 * pkin(4) + t1033) * t951;
t856 = (pkin(4) + pkin(9)) * t900 + t1003 - t951 * t938 - t947 * pkin(5);
t851 = t853 * t995 - t856 * t991;
t899 = qJDD(6) + t901;
t902 = -mrSges(7,1) * t929 + mrSges(7,2) * t930;
t848 = m(7) * t851 + mrSges(7,1) * t899 - mrSges(7,3) * t875 - t902 * t930 + t903 * t943;
t852 = t853 * t991 + t856 * t995;
t849 = m(7) * t852 - mrSges(7,2) * t899 + mrSges(7,3) * t874 + t902 * t929 - t904 * t943;
t837 = t995 * t848 + t991 * t849;
t917 = -mrSges(6,2) * t950 - mrSges(6,3) * t951;
t1007 = -m(6) * t859 - t901 * mrSges(6,1) - t951 * t917 - t837;
t916 = mrSges(5,1) * t950 + mrSges(5,2) * t951;
t934 = -mrSges(5,2) * t987 - mrSges(5,3) * t950;
t936 = mrSges(6,1) * t950 - mrSges(6,3) * t987;
t832 = m(5) * t863 - t901 * mrSges(5,3) - t951 * t916 + (t934 - t936) * t987 + (mrSges(5,1) - mrSges(6,2)) * t986 + t1007;
t864 = t941 + t1025;
t935 = mrSges(5,1) * t987 - mrSges(5,3) * t951;
t843 = m(5) * t864 - t986 * mrSges(5,2) - t987 * t935 + (-t916 - t917) * t950 + (-mrSges(5,3) - mrSges(6,1)) * t900 + t1005;
t827 = t1030 * t832 + t990 * t843;
t835 = t986 * mrSges(6,2) + t987 * t936 - t1007;
t878 = Ifges(7,5) * t930 + Ifges(7,6) * t929 + Ifges(7,3) * t943;
t880 = Ifges(7,1) * t930 + Ifges(7,4) * t929 + Ifges(7,5) * t943;
t839 = -mrSges(7,1) * t855 + mrSges(7,3) * t852 + Ifges(7,4) * t875 + Ifges(7,2) * t874 + Ifges(7,6) * t899 - t878 * t930 + t880 * t943;
t879 = Ifges(7,4) * t930 + Ifges(7,2) * t929 + Ifges(7,6) * t943;
t840 = mrSges(7,2) * t855 - mrSges(7,3) * t851 + Ifges(7,1) * t875 + Ifges(7,4) * t874 + Ifges(7,5) * t899 + t878 * t929 - t879 * t943;
t945 = Ifges(4,4) * t965 + Ifges(4,2) * t964 + Ifges(4,6) * t987;
t946 = Ifges(4,1) * t965 + Ifges(4,4) * t964 + Ifges(4,5) * t987;
t1042 = pkin(3) * t827 + qJ(5) * (-t950 * t917 + t1005) + t995 * t840 - t991 * t839 + t965 * t945 - t964 * t946 + Ifges(4,6) * t931 + Ifges(4,5) * t932 + mrSges(4,1) * t887 - mrSges(4,2) * t888 + mrSges(5,1) * t863 - mrSges(5,2) * t864 + mrSges(6,2) * t859 - mrSges(6,3) * t857 - pkin(4) * t835 - pkin(9) * t837 + t1036 * t951 + t1038 * t901 + (-qJ(5) * mrSges(6,1) - t1037) * t900 + (Ifges(4,3) - t1039) * t986 + t1023 * t950;
t952 = -mrSges(4,1) * t964 + mrSges(4,2) * t965;
t955 = -mrSges(4,2) * t987 + mrSges(4,3) * t964;
t824 = m(4) * t887 + mrSges(4,1) * t986 - mrSges(4,3) * t932 - t952 * t965 + t955 * t987 + t827;
t1014 = t1030 * t843 - t832 * t990;
t957 = mrSges(4,1) * t987 - mrSges(4,3) * t965;
t825 = m(4) * t888 - mrSges(4,2) * t986 + mrSges(4,3) * t931 + t952 * t964 - t957 * t987 + t1014;
t819 = t996 * t824 + t992 * t825;
t953 = -t997 * g(3) - t1028;
t962 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t993 + Ifges(3,2) * t997) * qJD(1);
t963 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t993 + Ifges(3,4) * t997) * qJD(1);
t1035 = mrSges(3,1) * t953 - mrSges(3,2) * t954 + Ifges(3,5) * t972 + Ifges(3,6) * t973 + Ifges(3,3) * qJDD(2) + pkin(2) * t819 + (t962 * t993 - t963 * t997) * qJD(1) + t1042;
t971 = (-mrSges(3,1) * t997 + mrSges(3,2) * t993) * qJD(1);
t1021 = qJD(1) * t997;
t976 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1021;
t817 = m(3) * t953 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t972 + qJD(2) * t976 - t971 * t1022 + t819;
t1015 = -t824 * t992 + t996 * t825;
t975 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1022;
t818 = m(3) * t954 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t973 - qJD(2) * t975 + t971 * t1021 + t1015;
t1016 = -t817 * t993 + t997 * t818;
t810 = m(2) * t979 - mrSges(2,1) * t999 - qJDD(1) * mrSges(2,2) + t1016;
t1026 = -t991 * t848 + t995 * t849;
t861 = t900 * pkin(4) + t1003;
t836 = m(6) * t861 - t900 * mrSges(6,2) - t901 * mrSges(6,3) - t950 * t936 - t951 * t937 + t1026;
t833 = m(5) * t877 + t900 * mrSges(5,1) + t901 * mrSges(5,2) + t950 * t934 + t951 * t935 + t836;
t1004 = m(4) * t933 - t931 * mrSges(4,1) + t932 * mrSges(4,2) - t964 * t955 + t965 * t957 + t833;
t966 = -t999 * pkin(7) + t1011;
t1002 = -m(3) * t966 + t973 * mrSges(3,1) - t972 * mrSges(3,2) + t976 * t1021 - t975 * t1022 - t1004;
t829 = m(2) * t978 + qJDD(1) * mrSges(2,1) - t999 * mrSges(2,2) + t1002;
t1027 = t994 * t810 + t998 * t829;
t812 = t997 * t817 + t993 * t818;
t1024 = t1037 * t950 - t1038 * t951 + t1039 * t987;
t1017 = t998 * t810 - t829 * t994;
t813 = -mrSges(5,1) * t877 - mrSges(6,1) * t857 + mrSges(6,2) * t861 + mrSges(5,3) * t864 - pkin(4) * t836 - pkin(5) * t1009 - pkin(9) * t1026 + t1023 * t987 + t1024 * t951 + t1031 * t901 + t1037 * t986 + t1043 * t900 - t995 * t839 - t991 * t840;
t1006 = mrSges(7,1) * t851 - mrSges(7,2) * t852 + Ifges(7,5) * t875 + Ifges(7,6) * t874 + Ifges(7,3) * t899 + t930 * t879 - t929 * t880;
t820 = mrSges(6,1) * t859 + mrSges(5,2) * t877 - mrSges(5,3) * t863 - mrSges(6,3) * t861 + pkin(5) * t837 - qJ(5) * t836 + t1024 * t950 - t1031 * t900 - t1036 * t987 + t1038 * t986 + t1040 * t901 + t1006;
t944 = Ifges(4,5) * t965 + Ifges(4,6) * t964 + Ifges(4,3) * t987;
t806 = -mrSges(4,1) * t933 + mrSges(4,3) * t888 + Ifges(4,4) * t932 + Ifges(4,2) * t931 + Ifges(4,6) * t986 - pkin(3) * t833 + qJ(4) * t1014 + t1030 * t813 + t990 * t820 - t965 * t944 + t987 * t946;
t807 = mrSges(4,2) * t933 - mrSges(4,3) * t887 + Ifges(4,1) * t932 + Ifges(4,4) * t931 + Ifges(4,5) * t986 - qJ(4) * t827 + t1030 * t820 - t990 * t813 + t964 * t944 - t987 * t945;
t961 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t993 + Ifges(3,6) * t997) * qJD(1);
t802 = -mrSges(3,1) * t966 + mrSges(3,3) * t954 + Ifges(3,4) * t972 + Ifges(3,2) * t973 + Ifges(3,6) * qJDD(2) - pkin(2) * t1004 + pkin(8) * t1015 + qJD(2) * t963 - t961 * t1022 + t996 * t806 + t992 * t807;
t804 = mrSges(3,2) * t966 - mrSges(3,3) * t953 + Ifges(3,1) * t972 + Ifges(3,4) * t973 + Ifges(3,5) * qJDD(2) - pkin(8) * t819 - qJD(2) * t962 + t961 * t1021 - t806 * t992 + t807 * t996;
t1008 = mrSges(2,1) * t978 - mrSges(2,2) * t979 + Ifges(2,3) * qJDD(1) + pkin(1) * t1002 + pkin(7) * t1016 + t997 * t802 + t993 * t804;
t805 = mrSges(2,1) * g(3) + mrSges(2,3) * t979 + t999 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t812 - t1035;
t800 = -mrSges(2,2) * g(3) - mrSges(2,3) * t978 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t999 - pkin(7) * t812 - t802 * t993 + t804 * t997;
t1 = [-m(1) * g(1) + t1017; -m(1) * g(2) + t1027; (-m(1) - m(2)) * g(3) + t812; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1027 + t998 * t800 - t994 * t805; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1017 + t994 * t800 + t998 * t805; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1008; t1008; t1035; t1042; t833; t835; t1006;];
tauJB  = t1;
