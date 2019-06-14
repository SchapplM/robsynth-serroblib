% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:52:59
% EndTime: 2019-05-05 17:53:07
% DurationCPUTime: 6.91s
% Computational Cost: add. (65447->325), mult. (156992->374), div. (0->0), fcn. (109831->8), ass. (0->140)
t1023 = Ifges(6,4) + Ifges(7,4);
t1038 = Ifges(6,2) + Ifges(7,2);
t1032 = Ifges(6,6) + Ifges(7,6);
t1037 = -2 * qJD(4);
t1036 = Ifges(4,1) + Ifges(5,2);
t1035 = Ifges(6,1) + Ifges(7,1);
t1022 = Ifges(4,5) - Ifges(5,4);
t1034 = Ifges(6,5) + Ifges(7,5);
t1033 = -Ifges(4,2) - Ifges(5,3);
t1021 = Ifges(4,6) - Ifges(5,5);
t1020 = -Ifges(5,6) - Ifges(4,4);
t1031 = Ifges(4,3) + Ifges(5,1);
t1030 = Ifges(6,3) + Ifges(7,3);
t969 = sin(pkin(9));
t1006 = t969 * qJD(1);
t972 = sin(qJ(3));
t1026 = cos(qJ(3));
t970 = cos(pkin(9));
t998 = t970 * t1026;
t947 = -qJD(1) * t998 + t1006 * t972;
t971 = sin(qJ(5));
t974 = cos(qJ(5));
t932 = -qJD(3) * t971 + t947 * t974;
t933 = qJD(3) * t974 + t947 * t971;
t987 = t1026 * t969 + t970 * t972;
t948 = t987 * qJD(1);
t944 = qJD(5) + t948;
t1029 = t1023 * t933 + t1032 * t944 + t1038 * t932;
t977 = qJD(1) ^ 2;
t1025 = pkin(2) * t970;
t973 = sin(qJ(1));
t975 = cos(qJ(1));
t954 = -g(1) * t975 - g(2) * t973;
t949 = -pkin(1) * t977 + qJDD(1) * qJ(2) + t954;
t1005 = qJD(1) * qJD(2);
t997 = -t970 * g(3) - 0.2e1 * t1005 * t969;
t904 = (-pkin(7) * qJDD(1) + t1025 * t977 - t949) * t969 + t997;
t1002 = qJDD(1) * t970;
t965 = t970 ^ 2;
t1017 = t965 * t977;
t930 = -t969 * g(3) + (0.2e1 * t1005 + t949) * t970;
t905 = -pkin(2) * t1017 + pkin(7) * t1002 + t930;
t871 = t1026 * t905 + t904 * t972;
t918 = pkin(3) * t947 - qJ(4) * t948;
t976 = qJD(3) ^ 2;
t866 = pkin(3) * t976 - qJDD(3) * qJ(4) + qJD(3) * t1037 + t918 * t947 - t871;
t1003 = qJDD(1) * t969;
t1007 = t948 * qJD(3);
t927 = -qJDD(1) * t998 + t1003 * t972 + t1007;
t940 = pkin(4) * t948 - qJD(3) * pkin(8);
t946 = t947 ^ 2;
t863 = -pkin(4) * t927 - pkin(8) * t946 + qJD(3) * t940 - t866;
t888 = -qJD(5) * t933 - qJDD(3) * t971 + t927 * t974;
t889 = qJD(5) * t932 + qJDD(3) * t974 + t927 * t971;
t897 = -mrSges(7,2) * t944 + mrSges(7,3) * t932;
t898 = -mrSges(6,2) * t944 + mrSges(6,3) * t932;
t901 = mrSges(6,1) * t944 - mrSges(6,3) * t933;
t899 = pkin(5) * t944 - qJ(6) * t933;
t931 = t932 ^ 2;
t856 = -pkin(5) * t888 - qJ(6) * t931 + t899 * t933 + qJDD(6) + t863;
t900 = mrSges(7,1) * t944 - mrSges(7,3) * t933;
t999 = m(7) * t856 + mrSges(7,2) * t889 + t900 * t933;
t1028 = -m(6) * t863 - t889 * mrSges(6,2) + (t897 + t898) * t932 + (mrSges(6,1) + mrSges(7,1)) * t888 - t933 * t901 - t999;
t1009 = qJD(3) * t1022 + t1020 * t947 + t1036 * t948;
t1010 = qJD(3) * t1021 - t1020 * t948 + t1033 * t947;
t1008 = t947 * qJD(3);
t964 = t969 ^ 2;
t953 = t973 * g(1) - g(2) * t975;
t993 = qJDD(2) - t953;
t926 = (-pkin(1) - t1025) * qJDD(1) + (-qJ(2) + (-t964 - t965) * pkin(7)) * t977 + t993;
t928 = qJDD(1) * t987 - t1008;
t978 = pkin(3) * t1007 + t948 * t1037 + (-t928 + t1008) * qJ(4) + t926;
t858 = -t946 * pkin(4) - t948 * t940 + (pkin(3) + pkin(8)) * t927 + t978;
t870 = t1026 * t904 - t905 * t972;
t867 = -qJDD(3) * pkin(3) - t976 * qJ(4) + t918 * t948 + qJDD(4) - t870;
t861 = (t947 * t948 - qJDD(3)) * pkin(8) + (t928 + t1008) * pkin(4) + t867;
t854 = t858 * t974 + t861 * t971;
t850 = -pkin(5) * t931 + qJ(6) * t888 + 0.2e1 * qJD(6) * t932 - t899 * t944 + t854;
t891 = -mrSges(7,1) * t932 + mrSges(7,2) * t933;
t1000 = m(7) * t850 + mrSges(7,3) * t888 + t891 * t932;
t1013 = -t1023 * t932 - t1034 * t944 - t1035 * t933;
t1014 = -t1030 * t944 - t1032 * t932 - t1034 * t933;
t851 = -t888 * mrSges(7,1) - t932 * t897 + t999;
t925 = qJDD(5) + t928;
t825 = -mrSges(6,1) * t863 + mrSges(6,3) * t854 - mrSges(7,1) * t856 + mrSges(7,3) * t850 - pkin(5) * t851 + qJ(6) * t1000 + (-qJ(6) * t900 - t1013) * t944 + t1014 * t933 + (-mrSges(7,2) * qJ(6) + t1032) * t925 + t1023 * t889 + t1038 * t888;
t938 = mrSges(5,1) * t947 - qJD(3) * mrSges(5,3);
t853 = -t971 * t858 + t861 * t974;
t848 = -0.2e1 * qJD(6) * t933 + (t932 * t944 - t889) * qJ(6) + (t932 * t933 + t925) * pkin(5) + t853;
t1001 = m(7) * t848 + mrSges(7,1) * t925 + t897 * t944;
t892 = -mrSges(6,1) * t932 + mrSges(6,2) * t933;
t837 = m(6) * t853 + t925 * mrSges(6,1) + t944 * t898 + (-t891 - t892) * t933 + (-mrSges(6,3) - mrSges(7,3)) * t889 + t1001;
t839 = m(6) * t854 + t888 * mrSges(6,3) + t932 * t892 + (-t900 - t901) * t944 + (-mrSges(6,2) - mrSges(7,2)) * t925 + t1000;
t835 = t974 * t837 + t971 * t839;
t920 = -mrSges(5,2) * t947 - mrSges(5,3) * t948;
t984 = -m(5) * t867 - mrSges(5,1) * t928 - t920 * t948 - t835;
t833 = qJDD(3) * mrSges(5,2) + qJD(3) * t938 - t984;
t845 = -t889 * mrSges(7,3) - t933 * t891 + t1001;
t834 = mrSges(6,2) * t863 + mrSges(7,2) * t856 - mrSges(6,3) * t853 - mrSges(7,3) * t848 - qJ(6) * t845 - t1014 * t932 + t1023 * t888 - t1029 * t944 + t1034 * t925 + t1035 * t889;
t939 = mrSges(5,1) * t948 + qJD(3) * mrSges(5,2);
t982 = -m(5) * t866 + qJDD(3) * mrSges(5,3) + qJD(3) * t939 - t1028;
t1027 = t1031 * qJDD(3) + t1009 * t947 + t1010 * t948 - t1021 * t927 + t1022 * t928 + mrSges(4,1) * t870 - mrSges(4,2) * t871 + mrSges(5,2) * t867 - mrSges(5,3) * t866 - pkin(3) * t833 - pkin(8) * t835 + qJ(4) * (-t927 * mrSges(5,1) - t947 * t920 + t982) - t971 * t825 + t974 * t834;
t1018 = mrSges(3,2) * t969;
t919 = mrSges(4,1) * t947 + mrSges(4,2) * t948;
t936 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t947;
t831 = m(4) * t870 - t928 * mrSges(4,3) - t948 * t919 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t936 - t938) * qJD(3) + t984;
t937 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t948;
t842 = (-t919 - t920) * t947 + (-mrSges(4,3) - mrSges(5,1)) * t927 - qJD(3) * t937 + m(4) * t871 + t982 - qJDD(3) * mrSges(4,2);
t824 = t1026 * t831 + t842 * t972;
t929 = -t969 * t949 + t997;
t989 = mrSges(3,3) * qJDD(1) + t977 * (-mrSges(3,1) * t970 + t1018);
t822 = m(3) * t929 - t969 * t989 + t824;
t994 = t1026 * t842 - t972 * t831;
t823 = m(3) * t930 + t970 * t989 + t994;
t995 = -t822 * t969 + t823 * t970;
t816 = m(2) * t954 - mrSges(2,1) * t977 - qJDD(1) * mrSges(2,2) + t995;
t945 = -qJDD(1) * pkin(1) - t977 * qJ(2) + t993;
t1015 = -t837 * t971 + t839 * t974;
t865 = t927 * pkin(3) + t978;
t832 = m(5) * t865 - mrSges(5,2) * t927 - mrSges(5,3) * t928 - t938 * t947 - t939 * t948 + t1015;
t983 = m(4) * t926 + mrSges(4,1) * t927 + t928 * mrSges(4,2) + t936 * t947 + t948 * t937 + t832;
t980 = -m(3) * t945 + mrSges(3,1) * t1002 - t983 + (t964 * t977 + t1017) * mrSges(3,3);
t827 = (mrSges(2,1) - t1018) * qJDD(1) - t977 * mrSges(2,2) + m(2) * t953 + t980;
t1016 = t816 * t973 + t827 * t975;
t818 = t822 * t970 + t823 * t969;
t1011 = -qJD(3) * t1031 + t1021 * t947 - t1022 * t948;
t996 = t816 * t975 - t827 * t973;
t992 = Ifges(3,1) * t969 + Ifges(3,4) * t970;
t991 = Ifges(3,4) * t969 + Ifges(3,2) * t970;
t990 = Ifges(3,5) * t969 + Ifges(3,6) * t970;
t812 = -mrSges(4,1) * t926 - mrSges(5,1) * t866 + mrSges(5,2) * t865 + mrSges(4,3) * t871 - pkin(3) * t832 - pkin(4) * t1028 - pkin(8) * t1015 + qJD(3) * t1009 + qJDD(3) * t1021 + t1011 * t948 - t1020 * t928 + t1033 * t927 - t974 * t825 - t971 * t834;
t981 = mrSges(6,1) * t853 + mrSges(7,1) * t848 - mrSges(6,2) * t854 - mrSges(7,2) * t850 + pkin(5) * t845 + t1013 * t932 + t1029 * t933 + t1030 * t925 + t1032 * t888 + t1034 * t889;
t813 = mrSges(5,1) * t867 + mrSges(4,2) * t926 - mrSges(4,3) * t870 - mrSges(5,3) * t865 + pkin(4) * t835 - qJ(4) * t832 - qJD(3) * t1010 + qJDD(3) * t1022 + t1011 * t947 + t1020 * t927 + t1036 * t928 + t981;
t951 = t990 * qJD(1);
t808 = -mrSges(3,1) * t945 + mrSges(3,3) * t930 - pkin(2) * t983 + pkin(7) * t994 + qJDD(1) * t991 - t1006 * t951 + t1026 * t812 + t972 * t813;
t810 = t970 * qJD(1) * t951 + mrSges(3,2) * t945 - mrSges(3,3) * t929 - pkin(7) * t824 + qJDD(1) * t992 + t1026 * t813 - t972 * t812;
t829 = mrSges(3,2) * t1003 - t980;
t985 = mrSges(2,1) * t953 - mrSges(2,2) * t954 + Ifges(2,3) * qJDD(1) - pkin(1) * t829 + qJ(2) * t995 + t808 * t970 + t810 * t969;
t811 = mrSges(2,3) * t954 - mrSges(3,1) * t929 + mrSges(3,2) * t930 - pkin(2) * t824 + mrSges(2,1) * g(3) - pkin(1) * t818 + (Ifges(2,6) - t990) * qJDD(1) + (-t969 * t991 + t970 * t992 + Ifges(2,5)) * t977 - t1027;
t806 = -mrSges(2,2) * g(3) - mrSges(2,3) * t953 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t977 - qJ(2) * t818 - t808 * t969 + t810 * t970;
t1 = [-m(1) * g(1) + t996; -m(1) * g(2) + t1016; (-m(1) - m(2)) * g(3) + t818; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1016 + t806 * t975 - t811 * t973; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t996 + t973 * t806 + t975 * t811; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t985; t985; t829; t1027; t833; t981; t851;];
tauJB  = t1;
