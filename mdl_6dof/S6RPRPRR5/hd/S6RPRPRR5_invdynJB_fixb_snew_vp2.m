% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:52:42
% EndTime: 2019-05-05 18:52:53
% DurationCPUTime: 9.85s
% Computational Cost: add. (137081->342), mult. (331942->416), div. (0->0), fcn. (245716->10), ass. (0->148)
t1037 = Ifges(4,1) + Ifges(5,1);
t1029 = Ifges(4,4) - Ifges(5,5);
t1028 = Ifges(4,5) + Ifges(5,4);
t1036 = Ifges(4,2) + Ifges(5,3);
t1027 = Ifges(4,6) - Ifges(5,6);
t1035 = Ifges(4,3) + Ifges(5,2);
t979 = cos(pkin(10));
t988 = qJD(1) ^ 2;
t1023 = t979 ^ 2 * t988;
t978 = sin(pkin(10));
t1034 = t978 ^ 2 * t988 + t1023;
t1031 = cos(qJ(3));
t1011 = t979 * t1031;
t1017 = qJD(1) * t978;
t982 = sin(qJ(3));
t945 = -qJD(1) * t1011 + t982 * t1017;
t1000 = t1031 * t978 + t979 * t982;
t946 = t1000 * qJD(1);
t1019 = t1028 * qJD(3) - t1029 * t945 + t1037 * t946;
t1021 = -qJD(3) * t1027 - t1029 * t946 + t1036 * t945;
t922 = mrSges(5,1) * t945 - mrSges(5,3) * t946;
t1015 = t945 * qJD(3);
t930 = t1000 * qJDD(1) - t1015;
t1014 = qJD(1) * qJD(2);
t1010 = -t979 * g(3) - 0.2e1 * t978 * t1014;
t983 = sin(qJ(1));
t986 = cos(qJ(1));
t952 = -g(1) * t986 - g(2) * t983;
t947 = -pkin(1) * t988 + qJDD(1) * qJ(2) + t952;
t905 = (pkin(2) * t979 * t988 - pkin(7) * qJDD(1) - t947) * t978 + t1010;
t1012 = qJDD(1) * t979;
t932 = -g(3) * t978 + (0.2e1 * t1014 + t947) * t979;
t906 = -pkin(2) * t1023 + pkin(7) * t1012 + t932;
t885 = t1031 * t905 - t982 * t906;
t921 = pkin(3) * t945 - qJ(4) * t946;
t987 = qJD(3) ^ 2;
t871 = -qJDD(3) * pkin(3) - t987 * qJ(4) + t946 * t921 + qJDD(4) - t885;
t861 = (-t930 - t1015) * pkin(8) + (t945 * t946 - qJDD(3)) * pkin(4) + t871;
t1032 = 2 * qJD(4);
t886 = t1031 * t906 + t982 * t905;
t870 = -pkin(3) * t987 + qJDD(3) * qJ(4) + qJD(3) * t1032 - t945 * t921 + t886;
t1013 = qJDD(1) * t978;
t1016 = qJD(3) * t946;
t929 = -qJDD(1) * t1011 + t982 * t1013 + t1016;
t939 = -qJD(3) * pkin(4) - pkin(8) * t946;
t944 = t945 ^ 2;
t863 = -pkin(4) * t944 + pkin(8) * t929 + qJD(3) * t939 + t870;
t981 = sin(qJ(5));
t985 = cos(qJ(5));
t857 = t981 * t861 + t985 * t863;
t914 = t945 * t985 - t946 * t981;
t915 = t945 * t981 + t946 * t985;
t895 = -pkin(5) * t914 - pkin(9) * t915;
t971 = -qJD(3) + qJD(5);
t967 = t971 ^ 2;
t968 = -qJDD(3) + qJDD(5);
t853 = -pkin(5) * t967 + pkin(9) * t968 + t895 * t914 + t857;
t951 = t983 * g(1) - t986 * g(2);
t943 = -qJDD(1) * pkin(1) - t988 * qJ(2) + qJDD(2) - t951;
t928 = -pkin(2) * t1012 - pkin(7) * t1034 + t943;
t994 = t929 * pkin(3) + t928 + (t1015 - t930) * qJ(4);
t859 = -pkin(3) * t1016 - t929 * pkin(4) - t944 * pkin(8) - t994 + (t1032 + t939) * t946;
t882 = -qJD(5) * t915 + t929 * t985 - t930 * t981;
t883 = qJD(5) * t914 + t929 * t981 + t930 * t985;
t854 = t859 + (t915 * t971 - t882) * pkin(5) + (-t914 * t971 - t883) * pkin(9);
t980 = sin(qJ(6));
t984 = cos(qJ(6));
t850 = -t853 * t980 + t854 * t984;
t898 = -t915 * t980 + t971 * t984;
t866 = qJD(6) * t898 + t883 * t984 + t968 * t980;
t881 = qJDD(6) - t882;
t899 = t915 * t984 + t971 * t980;
t884 = -mrSges(7,1) * t898 + mrSges(7,2) * t899;
t907 = qJD(6) - t914;
t887 = -mrSges(7,2) * t907 + mrSges(7,3) * t898;
t846 = m(7) * t850 + mrSges(7,1) * t881 - mrSges(7,3) * t866 - t884 * t899 + t887 * t907;
t851 = t853 * t984 + t854 * t980;
t865 = -qJD(6) * t899 - t883 * t980 + t968 * t984;
t888 = mrSges(7,1) * t907 - mrSges(7,3) * t899;
t847 = m(7) * t851 - mrSges(7,2) * t881 + mrSges(7,3) * t865 + t884 * t898 - t888 * t907;
t1005 = -t846 * t980 + t984 * t847;
t894 = -mrSges(6,1) * t914 + mrSges(6,2) * t915;
t902 = mrSges(6,1) * t971 - mrSges(6,3) * t915;
t834 = m(6) * t857 - mrSges(6,2) * t968 + mrSges(6,3) * t882 + t894 * t914 - t902 * t971 + t1005;
t856 = t861 * t985 - t863 * t981;
t901 = -mrSges(6,2) * t971 + mrSges(6,3) * t914;
t852 = -pkin(5) * t968 - pkin(9) * t967 + t895 * t915 - t856;
t996 = -m(7) * t852 + t865 * mrSges(7,1) - mrSges(7,2) * t866 + t898 * t887 - t888 * t899;
t842 = m(6) * t856 + mrSges(6,1) * t968 - mrSges(6,3) * t883 - t894 * t915 + t901 * t971 + t996;
t827 = t981 * t834 + t985 * t842;
t938 = -mrSges(5,2) * t945 + qJD(3) * mrSges(5,3);
t995 = -m(5) * t871 + qJDD(3) * mrSges(5,1) + qJD(3) * t938 - t827;
t826 = t930 * mrSges(5,2) + t946 * t922 - t995;
t872 = Ifges(7,5) * t899 + Ifges(7,6) * t898 + Ifges(7,3) * t907;
t874 = Ifges(7,1) * t899 + Ifges(7,4) * t898 + Ifges(7,5) * t907;
t840 = -mrSges(7,1) * t852 + mrSges(7,3) * t851 + Ifges(7,4) * t866 + Ifges(7,2) * t865 + Ifges(7,6) * t881 - t872 * t899 + t874 * t907;
t873 = Ifges(7,4) * t899 + Ifges(7,2) * t898 + Ifges(7,6) * t907;
t841 = mrSges(7,2) * t852 - mrSges(7,3) * t850 + Ifges(7,1) * t866 + Ifges(7,4) * t865 + Ifges(7,5) * t881 + t872 * t898 - t873 * t907;
t890 = Ifges(6,4) * t915 + Ifges(6,2) * t914 + Ifges(6,6) * t971;
t891 = Ifges(6,1) * t915 + Ifges(6,4) * t914 + Ifges(6,5) * t971;
t993 = -mrSges(6,1) * t856 + mrSges(6,2) * t857 - Ifges(6,5) * t883 - Ifges(6,6) * t882 - Ifges(6,3) * t968 - pkin(5) * t996 - pkin(9) * t1005 - t984 * t840 - t980 * t841 - t915 * t890 + t914 * t891;
t1006 = t985 * t834 - t981 * t842;
t937 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t946;
t999 = m(5) * t870 + qJDD(3) * mrSges(5,3) + qJD(3) * t937 + t1006;
t1033 = t1035 * qJDD(3) + t1019 * t945 - t1021 * t946 - t1027 * t929 + t1028 * t930 + mrSges(4,1) * t885 - mrSges(5,1) * t871 - mrSges(4,2) * t886 + mrSges(5,3) * t870 - pkin(3) * t826 - pkin(4) * t827 + qJ(4) * (-t929 * mrSges(5,2) - t945 * t922 + t999) + t993;
t1030 = -mrSges(4,3) - mrSges(5,2);
t1025 = mrSges(3,2) * t978;
t1001 = mrSges(3,3) * qJDD(1) + t988 * (-mrSges(3,1) * t979 + t1025);
t1018 = -mrSges(4,1) * t945 - mrSges(4,2) * t946 - t922;
t936 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t946;
t823 = m(4) * t886 - qJDD(3) * mrSges(4,2) - qJD(3) * t936 + t1018 * t945 + t1030 * t929 + t999;
t935 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t945;
t824 = m(4) * t885 + qJDD(3) * mrSges(4,1) + qJD(3) * t935 + t1018 * t946 + t1030 * t930 + t995;
t817 = t1031 * t824 + t982 * t823;
t931 = -t978 * t947 + t1010;
t815 = m(3) * t931 - t1001 * t978 + t817;
t1007 = t1031 * t823 - t982 * t824;
t816 = m(3) * t932 + t1001 * t979 + t1007;
t1008 = -t815 * t978 + t979 * t816;
t809 = m(2) * t952 - mrSges(2,1) * t988 - qJDD(1) * mrSges(2,2) + t1008;
t868 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t946 + t994;
t836 = t984 * t846 + t980 * t847;
t998 = -m(6) * t859 + t882 * mrSges(6,1) - t883 * mrSges(6,2) + t914 * t901 - t915 * t902 - t836;
t832 = m(5) * t868 + t929 * mrSges(5,1) - t930 * mrSges(5,3) - t946 * t937 + t945 * t938 + t998;
t991 = m(4) * t928 + t929 * mrSges(4,1) + t930 * mrSges(4,2) + t945 * t935 + t946 * t936 + t832;
t990 = -m(3) * t943 + mrSges(3,1) * t1012 + mrSges(3,3) * t1034 - t991;
t829 = (mrSges(2,1) - t1025) * qJDD(1) + t990 - t988 * mrSges(2,2) + m(2) * t951;
t1022 = t983 * t809 + t986 * t829;
t811 = t979 * t815 + t978 * t816;
t1020 = -qJD(3) * t1035 + t1027 * t945 - t1028 * t946;
t1009 = t986 * t809 - t829 * t983;
t1004 = Ifges(3,1) * t978 + Ifges(3,4) * t979;
t1003 = Ifges(3,4) * t978 + Ifges(3,2) * t979;
t1002 = Ifges(3,5) * t978 + Ifges(3,6) * t979;
t889 = Ifges(6,5) * t915 + Ifges(6,6) * t914 + Ifges(6,3) * t971;
t818 = mrSges(6,2) * t859 - mrSges(6,3) * t856 + Ifges(6,1) * t883 + Ifges(6,4) * t882 + Ifges(6,5) * t968 - pkin(9) * t836 - t840 * t980 + t841 * t984 + t889 * t914 - t890 * t971;
t992 = mrSges(7,1) * t850 - mrSges(7,2) * t851 + Ifges(7,5) * t866 + Ifges(7,6) * t865 + Ifges(7,3) * t881 + t873 * t899 - t874 * t898;
t819 = -mrSges(6,1) * t859 + mrSges(6,3) * t857 + Ifges(6,4) * t883 + Ifges(6,2) * t882 + Ifges(6,6) * t968 - pkin(5) * t836 - t889 * t915 + t891 * t971 - t992;
t805 = -mrSges(4,1) * t928 - mrSges(5,1) * t868 + mrSges(5,2) * t870 + mrSges(4,3) * t886 - pkin(3) * t832 - pkin(4) * t998 - pkin(8) * t1006 + t1019 * qJD(3) + t1027 * qJDD(3) + t1020 * t946 + t1029 * t930 - t1036 * t929 - t981 * t818 - t985 * t819;
t806 = mrSges(4,2) * t928 + mrSges(5,2) * t871 - mrSges(4,3) * t885 - mrSges(5,3) * t868 - pkin(8) * t827 - qJ(4) * t832 + t1021 * qJD(3) + t1028 * qJDD(3) + t1020 * t945 - t1029 * t929 + t1037 * t930 + t985 * t818 - t981 * t819;
t949 = t1002 * qJD(1);
t801 = -mrSges(3,1) * t943 + mrSges(3,3) * t932 - pkin(2) * t991 + pkin(7) * t1007 + t1003 * qJDD(1) - t949 * t1017 + t1031 * t805 + t982 * t806;
t803 = t979 * qJD(1) * t949 + mrSges(3,2) * t943 - mrSges(3,3) * t931 - pkin(7) * t817 + t1004 * qJDD(1) + t1031 * t806 - t982 * t805;
t831 = mrSges(3,2) * t1013 - t990;
t997 = mrSges(2,1) * t951 - mrSges(2,2) * t952 + Ifges(2,3) * qJDD(1) - pkin(1) * t831 + qJ(2) * t1008 + t979 * t801 + t978 * t803;
t804 = -pkin(1) * t811 + mrSges(2,1) * g(3) + (Ifges(2,6) - t1002) * qJDD(1) - pkin(2) * t817 + mrSges(2,3) * t952 - mrSges(3,1) * t931 + mrSges(3,2) * t932 + (-t978 * t1003 + t979 * t1004 + Ifges(2,5)) * t988 - t1033;
t799 = -mrSges(2,2) * g(3) - mrSges(2,3) * t951 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t988 - qJ(2) * t811 - t801 * t978 + t803 * t979;
t1 = [-m(1) * g(1) + t1009; -m(1) * g(2) + t1022; (-m(1) - m(2)) * g(3) + t811; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1022 + t986 * t799 - t983 * t804; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1009 + t983 * t799 + t986 * t804; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t997; t997; t831; t1033; t826; -t993; t992;];
tauJB  = t1;
