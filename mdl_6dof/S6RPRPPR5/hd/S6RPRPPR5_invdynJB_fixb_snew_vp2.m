% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-05-05 17:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:58:23
% EndTime: 2019-05-05 16:58:33
% DurationCPUTime: 9.93s
% Computational Cost: add. (137996->344), mult. (342827->417), div. (0->0), fcn. (247582->10), ass. (0->147)
t1003 = -2 * qJD(4);
t1002 = Ifges(4,1) + Ifges(5,2);
t996 = (Ifges(4,5) - Ifges(5,4));
t1001 = -Ifges(4,2) - Ifges(5,3);
t995 = (Ifges(4,6) - Ifges(5,5));
t994 = -Ifges(5,6) - Ifges(4,4);
t1000 = Ifges(4,3) + Ifges(5,1);
t956 = qJD(1) ^ 2;
t952 = sin(qJ(1));
t954 = cos(qJ(1));
t931 = -g(1) * t954 - g(2) * t952;
t926 = -pkin(1) * t956 + qJDD(1) * qJ(2) + t931;
t947 = sin(pkin(9));
t949 = cos(pkin(9));
t982 = qJD(1) * qJD(2);
t977 = -t949 * g(3) - 0.2e1 * t947 * t982;
t997 = pkin(2) * t949;
t875 = (-pkin(7) * qJDD(1) + t956 * t997 - t926) * t947 + t977;
t905 = -g(3) * t947 + (t926 + 0.2e1 * t982) * t949;
t979 = qJDD(1) * t949;
t942 = t949 ^ 2;
t991 = t942 * t956;
t882 = -pkin(2) * t991 + pkin(7) * t979 + t905;
t951 = sin(qJ(3));
t998 = cos(qJ(3));
t859 = t951 * t875 + t998 * t882;
t978 = t949 * t998;
t985 = qJD(1) * t947;
t924 = -qJD(1) * t978 + t951 * t985;
t967 = t947 * t998 + t949 * t951;
t925 = t967 * qJD(1);
t893 = pkin(3) * t924 - qJ(4) * t925;
t955 = qJD(3) ^ 2;
t846 = pkin(3) * t955 - qJDD(3) * qJ(4) + (qJD(3) * t1003) + t924 * t893 - t859;
t980 = qJDD(1) * t947;
t983 = t925 * qJD(3);
t902 = -qJDD(1) * t978 + t951 * t980 + t983;
t915 = pkin(4) * t925 - (qJD(3) * qJ(5));
t923 = t924 ^ 2;
t941 = t947 ^ 2;
t930 = t952 * g(1) - t954 * g(2);
t972 = qJDD(2) - t930;
t901 = (-pkin(1) - t997) * qJDD(1) + (-qJ(2) + (-t941 - t942) * pkin(7)) * t956 + t972;
t984 = t924 * qJD(3);
t903 = qJDD(1) * t967 - t984;
t957 = pkin(3) * t983 + t925 * t1003 + (-t903 + t984) * qJ(4) + t901;
t837 = -t923 * pkin(4) - t925 * t915 + (pkin(3) + qJ(5)) * t902 + t957;
t858 = t875 * t998 - t951 * t882;
t850 = -qJDD(3) * pkin(3) - t955 * qJ(4) + t925 * t893 + qJDD(4) - t858;
t840 = (t924 * t925 - qJDD(3)) * qJ(5) + (t903 + t984) * pkin(4) + t850;
t946 = sin(pkin(10));
t948 = cos(pkin(10));
t910 = qJD(3) * t948 + t924 * t946;
t832 = -0.2e1 * qJD(5) * t910 - t946 * t837 + t948 * t840;
t881 = qJDD(3) * t948 + t902 * t946;
t909 = -qJD(3) * t946 + t924 * t948;
t830 = (t909 * t925 - t881) * pkin(8) + (t909 * t910 + t903) * pkin(5) + t832;
t833 = 0.2e1 * qJD(5) * t909 + t948 * t837 + t946 * t840;
t879 = pkin(5) * t925 - pkin(8) * t910;
t880 = -qJDD(3) * t946 + t902 * t948;
t908 = t909 ^ 2;
t831 = -pkin(5) * t908 + pkin(8) * t880 - t879 * t925 + t833;
t950 = sin(qJ(6));
t953 = cos(qJ(6));
t829 = t830 * t950 + t831 * t953;
t842 = -pkin(4) * t902 - qJ(5) * t923 + qJD(3) * t915 + qJDD(5) - t846;
t835 = -pkin(5) * t880 - pkin(8) * t908 + t879 * t910 + t842;
t867 = t909 * t950 + t910 * t953;
t851 = -t867 * qJD(6) + t880 * t953 - t881 * t950;
t866 = t909 * t953 - t910 * t950;
t852 = t866 * qJD(6) + t880 * t950 + t881 * t953;
t921 = qJD(6) + t925;
t853 = Ifges(7,5) * t867 + Ifges(7,6) * t866 + Ifges(7,3) * t921;
t855 = Ifges(7,1) * t867 + Ifges(7,4) * t866 + Ifges(7,5) * t921;
t900 = qJDD(6) + t903;
t815 = -mrSges(7,1) * t835 + mrSges(7,3) * t829 + Ifges(7,4) * t852 + Ifges(7,2) * t851 + Ifges(7,6) * t900 - t867 * t853 + t855 * t921;
t828 = t830 * t953 - t831 * t950;
t854 = Ifges(7,4) * t867 + Ifges(7,2) * t866 + Ifges(7,6) * t921;
t816 = mrSges(7,2) * t835 - mrSges(7,3) * t828 + Ifges(7,1) * t852 + Ifges(7,4) * t851 + Ifges(7,5) * t900 + t866 * t853 - t854 * t921;
t862 = Ifges(6,5) * t910 + Ifges(6,6) * t909 + Ifges(6,3) * t925;
t864 = Ifges(6,1) * t910 + Ifges(6,4) * t909 + Ifges(6,5) * t925;
t860 = -mrSges(7,2) * t921 + t866 * mrSges(7,3);
t861 = mrSges(7,1) * t921 - t867 * mrSges(7,3);
t965 = m(7) * t835 - t851 * mrSges(7,1) + t852 * mrSges(7,2) - t866 * t860 + t867 * t861;
t857 = -mrSges(7,1) * t866 + mrSges(7,2) * t867;
t823 = m(7) * t828 + mrSges(7,1) * t900 - t852 * mrSges(7,3) - t867 * t857 + t860 * t921;
t824 = m(7) * t829 - mrSges(7,2) * t900 + t851 * mrSges(7,3) + t866 * t857 - t861 * t921;
t973 = -t823 * t950 + t953 * t824;
t799 = -mrSges(6,1) * t842 + mrSges(6,3) * t833 + Ifges(6,4) * t881 + Ifges(6,2) * t880 + Ifges(6,6) * t903 - pkin(5) * t965 + pkin(8) * t973 + t953 * t815 + t950 * t816 - t910 * t862 + t925 * t864;
t814 = t953 * t823 + t950 * t824;
t863 = Ifges(6,4) * t910 + Ifges(6,2) * t909 + Ifges(6,6) * t925;
t800 = mrSges(6,2) * t842 - mrSges(6,3) * t832 + Ifges(6,1) * t881 + Ifges(6,4) * t880 + Ifges(6,5) * t903 - pkin(8) * t814 - t815 * t950 + t816 * t953 + t862 * t909 - t863 * t925;
t916 = mrSges(5,1) * t924 - (qJD(3) * mrSges(5,3));
t868 = -mrSges(6,1) * t909 + mrSges(6,2) * t910;
t877 = -mrSges(6,2) * t925 + mrSges(6,3) * t909;
t812 = m(6) * t832 + mrSges(6,1) * t903 - mrSges(6,3) * t881 - t868 * t910 + t877 * t925 + t814;
t878 = mrSges(6,1) * t925 - mrSges(6,3) * t910;
t813 = m(6) * t833 - mrSges(6,2) * t903 + mrSges(6,3) * t880 + t868 * t909 - t878 * t925 + t973;
t809 = t948 * t812 + t946 * t813;
t895 = -mrSges(5,2) * t924 - mrSges(5,3) * t925;
t963 = -m(5) * t850 - t903 * mrSges(5,1) - t925 * t895 - t809;
t808 = qJDD(3) * mrSges(5,2) + qJD(3) * t916 - t963;
t826 = m(6) * t842 - t880 * mrSges(6,1) + t881 * mrSges(6,2) - t909 * t877 + t910 * t878 + t965;
t917 = mrSges(5,1) * t925 + qJD(3) * mrSges(5,2);
t960 = -m(5) * t846 + qJDD(3) * mrSges(5,3) + qJD(3) * t917 + t826;
t986 = (qJD(3) * t996) + t1002 * t925 + t924 * t994;
t987 = (qJD(3) * t995) + t1001 * t924 - t925 * t994;
t999 = t1000 * qJDD(3) - t995 * t902 + t996 * t903 + t986 * t924 + t987 * t925 + mrSges(4,1) * t858 - mrSges(4,2) * t859 + mrSges(5,2) * t850 - mrSges(5,3) * t846 - pkin(3) * t808 + qJ(4) * (-t902 * mrSges(5,1) - t924 * t895 + t960) - qJ(5) * t809 - t946 * t799 + t948 * t800;
t992 = mrSges(3,2) * t947;
t894 = mrSges(4,1) * t924 + mrSges(4,2) * t925;
t913 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t924;
t806 = m(4) * t858 - t903 * mrSges(4,3) - t925 * t894 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t913 - t916) * qJD(3) + t963;
t914 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t925;
t819 = (-mrSges(4,3) - mrSges(5,1)) * t902 + (-t894 - t895) * t924 - qJDD(3) * mrSges(4,2) + t960 - qJD(3) * t914 + m(4) * t859;
t798 = t998 * t806 + t951 * t819;
t904 = -t947 * t926 + t977;
t968 = qJDD(1) * mrSges(3,3) + t956 * (-mrSges(3,1) * t949 + t992);
t796 = m(3) * t904 - t947 * t968 + t798;
t974 = -t951 * t806 + t998 * t819;
t797 = m(3) * t905 + t949 * t968 + t974;
t975 = -t796 * t947 + t949 * t797;
t790 = m(2) * t931 - mrSges(2,1) * t956 - qJDD(1) * mrSges(2,2) + t975;
t922 = -qJDD(1) * pkin(1) - t956 * qJ(2) + t972;
t845 = t902 * pkin(3) + t957;
t989 = -t946 * t812 + t948 * t813;
t807 = m(5) * t845 - t902 * mrSges(5,2) - t903 * mrSges(5,3) - t924 * t916 - t925 * t917 + t989;
t961 = m(4) * t901 + t902 * mrSges(4,1) + t903 * mrSges(4,2) + t924 * t913 + t925 * t914 + t807;
t959 = -m(3) * t922 + mrSges(3,1) * t979 - t961 + (t941 * t956 + t991) * mrSges(3,3);
t802 = (mrSges(2,1) - t992) * qJDD(1) + t959 - t956 * mrSges(2,2) + m(2) * t930;
t990 = t952 * t790 + t954 * t802;
t792 = t949 * t796 + t947 * t797;
t988 = -qJD(3) * t1000 + t924 * t995 - t925 * t996;
t976 = t954 * t790 - t802 * t952;
t971 = Ifges(3,1) * t947 + Ifges(3,4) * t949;
t970 = Ifges(3,4) * t947 + Ifges(3,2) * t949;
t969 = Ifges(3,5) * t947 + Ifges(3,6) * t949;
t786 = -mrSges(4,1) * t901 - mrSges(5,1) * t846 + mrSges(5,2) * t845 + mrSges(4,3) * t859 - pkin(3) * t807 + pkin(4) * t826 - qJ(5) * t989 + t986 * qJD(3) + t995 * qJDD(3) + t1001 * t902 - t948 * t799 - t946 * t800 - t994 * t903 + t988 * t925;
t962 = mrSges(7,1) * t828 - mrSges(7,2) * t829 + Ifges(7,5) * t852 + Ifges(7,6) * t851 + Ifges(7,3) * t900 + t867 * t854 - t866 * t855;
t787 = t988 * t924 + (Ifges(6,3) + t1002) * t903 + t994 * t902 - t987 * qJD(3) + t962 - t909 * t864 + t910 * t863 + mrSges(4,2) * t901 + Ifges(6,6) * t880 + Ifges(6,5) * t881 - mrSges(4,3) * t858 + mrSges(5,1) * t850 - mrSges(5,3) * t845 - mrSges(6,2) * t833 + mrSges(6,1) * t832 + pkin(5) * t814 + pkin(4) * t809 - qJ(4) * t807 + t996 * qJDD(3);
t928 = t969 * qJD(1);
t782 = -mrSges(3,1) * t922 + mrSges(3,3) * t905 - pkin(2) * t961 + pkin(7) * t974 + qJDD(1) * t970 + t786 * t998 + t951 * t787 - t928 * t985;
t784 = t949 * qJD(1) * t928 + mrSges(3,2) * t922 - mrSges(3,3) * t904 - pkin(7) * t798 + qJDD(1) * t971 - t951 * t786 + t787 * t998;
t804 = mrSges(3,2) * t980 - t959;
t964 = mrSges(2,1) * t930 - mrSges(2,2) * t931 + Ifges(2,3) * qJDD(1) - pkin(1) * t804 + qJ(2) * t975 + t949 * t782 + t947 * t784;
t785 = -pkin(1) * t792 + (Ifges(2,6) - t969) * qJDD(1) + mrSges(2,1) * g(3) + mrSges(2,3) * t931 - mrSges(3,1) * t904 + mrSges(3,2) * t905 - pkin(2) * t798 + (-t947 * t970 + t949 * t971 + Ifges(2,5)) * t956 - t999;
t780 = -mrSges(2,2) * g(3) - mrSges(2,3) * t930 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t956 - qJ(2) * t792 - t782 * t947 + t784 * t949;
t1 = [-m(1) * g(1) + t976; -m(1) * g(2) + t990; (-m(1) - m(2)) * g(3) + t792; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t990 + t954 * t780 - t952 * t785; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t976 + t952 * t780 + t954 * t785; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t964; t964; t804; t999; t808; t826; t962;];
tauJB  = t1;
