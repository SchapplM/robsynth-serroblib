% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR4
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
% Datum: 2019-05-07 04:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:42:15
% EndTime: 2019-05-07 04:42:36
% DurationCPUTime: 13.18s
% Computational Cost: add. (202325->365), mult. (420592->441), div. (0->0), fcn. (294938->10), ass. (0->147)
t1013 = -2 * qJD(4);
t1012 = Ifges(5,1) + Ifges(6,1);
t1004 = Ifges(5,4) - Ifges(6,5);
t1003 = Ifges(5,5) + Ifges(6,4);
t1011 = Ifges(5,2) + Ifges(6,3);
t1010 = -Ifges(6,2) - Ifges(5,3);
t1002 = Ifges(6,6) - Ifges(5,6);
t1001 = cos(pkin(10));
t964 = sin(qJ(1));
t968 = cos(qJ(1));
t951 = t964 * g(1) - t968 * g(2);
t970 = qJD(1) ^ 2;
t931 = -qJDD(1) * pkin(1) - t970 * pkin(7) - t951;
t963 = sin(qJ(2));
t967 = cos(qJ(2));
t991 = qJD(1) * qJD(2);
t988 = t967 * t991;
t945 = t963 * qJDD(1) + t988;
t989 = t963 * t991;
t946 = t967 * qJDD(1) - t989;
t892 = (-t945 - t988) * pkin(8) + (-t946 + t989) * pkin(2) + t931;
t952 = -t968 * g(1) - t964 * g(2);
t932 = -t970 * pkin(1) + qJDD(1) * pkin(7) + t952;
t922 = -t963 * g(3) + t967 * t932;
t944 = (-pkin(2) * t967 - pkin(8) * t963) * qJD(1);
t969 = qJD(2) ^ 2;
t992 = t967 * qJD(1);
t897 = -t969 * pkin(2) + qJDD(2) * pkin(8) + t944 * t992 + t922;
t962 = sin(qJ(3));
t966 = cos(qJ(3));
t865 = t966 * t892 - t962 * t897;
t994 = qJD(1) * t963;
t941 = t966 * qJD(2) - t962 * t994;
t913 = t941 * qJD(3) + t962 * qJDD(2) + t966 * t945;
t940 = qJDD(3) - t946;
t942 = t962 * qJD(2) + t966 * t994;
t955 = qJD(3) - t992;
t854 = (t941 * t955 - t913) * qJ(4) + (t941 * t942 + t940) * pkin(3) + t865;
t866 = t962 * t892 + t966 * t897;
t912 = -t942 * qJD(3) + t966 * qJDD(2) - t962 * t945;
t919 = t955 * pkin(3) - t942 * qJ(4);
t939 = t941 ^ 2;
t857 = -t939 * pkin(3) + t912 * qJ(4) - t955 * t919 + t866;
t960 = sin(pkin(10));
t916 = t1001 * t942 + t960 * t941;
t842 = t1001 * t854 + t916 * t1013 - t960 * t857;
t915 = -t1001 * t941 + t960 * t942;
t1000 = t915 * t955;
t881 = t1001 * t913 + t960 * t912;
t921 = -t967 * g(3) - t963 * t932;
t978 = qJDD(2) * pkin(2) + t969 * pkin(8) - t944 * t994 + t921;
t974 = t912 * pkin(3) + t939 * qJ(4) - t942 * t919 - qJDD(4) + t978;
t1009 = (-t881 + t1000) * qJ(5) - t974;
t887 = t915 * pkin(4) - t916 * qJ(5);
t954 = t955 ^ 2;
t841 = -t940 * pkin(4) - t954 * qJ(5) + t916 * t887 + qJDD(5) - t842;
t835 = (-t881 - t1000) * pkin(9) + (t915 * t916 - t940) * pkin(5) + t841;
t1006 = 2 * qJD(5);
t843 = t1001 * t857 + t915 * t1013 + t960 * t854;
t840 = -t954 * pkin(4) + t940 * qJ(5) + t955 * t1006 - t915 * t887 + t843;
t880 = -t1001 * t912 + t960 * t913;
t902 = -t955 * pkin(5) - t916 * pkin(9);
t914 = t915 ^ 2;
t836 = -t914 * pkin(5) + t880 * pkin(9) + t955 * t902 + t840;
t961 = sin(qJ(6));
t965 = cos(qJ(6));
t834 = t961 * t835 + t965 * t836;
t838 = -t914 * pkin(9) + (-pkin(4) - pkin(5)) * t880 + (-pkin(4) * t955 + t1006 + t902) * t916 - t1009;
t886 = t961 * t915 + t965 * t916;
t850 = -t886 * qJD(6) + t965 * t880 - t961 * t881;
t885 = t965 * t915 - t961 * t916;
t851 = t885 * qJD(6) + t961 * t880 + t965 * t881;
t953 = qJD(6) - t955;
t858 = Ifges(7,5) * t886 + Ifges(7,6) * t885 + Ifges(7,3) * t953;
t860 = Ifges(7,1) * t886 + Ifges(7,4) * t885 + Ifges(7,5) * t953;
t936 = qJDD(6) - t940;
t822 = -mrSges(7,1) * t838 + mrSges(7,3) * t834 + Ifges(7,4) * t851 + Ifges(7,2) * t850 + Ifges(7,6) * t936 - t886 * t858 + t953 * t860;
t833 = t965 * t835 - t961 * t836;
t859 = Ifges(7,4) * t886 + Ifges(7,2) * t885 + Ifges(7,6) * t953;
t823 = mrSges(7,2) * t838 - mrSges(7,3) * t833 + Ifges(7,1) * t851 + Ifges(7,4) * t850 + Ifges(7,5) * t936 + t885 * t858 - t953 * t859;
t845 = -0.2e1 * qJD(5) * t916 + (t916 * t955 + t880) * pkin(4) + t1009;
t900 = -t955 * mrSges(6,1) + t916 * mrSges(6,2);
t901 = -t915 * mrSges(6,2) + t955 * mrSges(6,3);
t869 = -t953 * mrSges(7,2) + t885 * mrSges(7,3);
t870 = t953 * mrSges(7,1) - t886 * mrSges(7,3);
t981 = -m(7) * t838 + t850 * mrSges(7,1) - t851 * mrSges(7,2) + t885 * t869 - t886 * t870;
t828 = m(6) * t845 + t880 * mrSges(6,1) - t881 * mrSges(6,3) - t916 * t900 + t915 * t901 + t981;
t863 = -t885 * mrSges(7,1) + t886 * mrSges(7,2);
t830 = m(7) * t833 + t936 * mrSges(7,1) - t851 * mrSges(7,3) - t886 * t863 + t953 * t869;
t831 = m(7) * t834 - t936 * mrSges(7,2) + t850 * mrSges(7,3) + t885 * t863 - t953 * t870;
t984 = -t961 * t830 + t965 * t831;
t996 = t1003 * t955 - t1004 * t915 + t1012 * t916;
t997 = -t1002 * t915 - t1003 * t916 + t1010 * t955;
t808 = mrSges(5,1) * t974 - mrSges(6,1) * t845 + mrSges(6,2) * t840 + mrSges(5,3) * t843 - pkin(4) * t828 - pkin(5) * t981 - pkin(9) * t984 - t1002 * t940 + t1004 * t881 - t1011 * t880 - t965 * t822 - t961 * t823 + t997 * t916 + t996 * t955;
t821 = t965 * t830 + t961 * t831;
t998 = t1002 * t955 - t1004 * t916 + t1011 * t915;
t809 = -mrSges(5,2) * t974 + mrSges(6,2) * t841 - mrSges(5,3) * t842 - mrSges(6,3) * t845 - pkin(9) * t821 - qJ(5) * t828 + t1003 * t940 - t1004 * t880 + t1012 * t881 - t961 * t822 + t965 * t823 + t997 * t915 + t998 * t955;
t898 = -t955 * mrSges(5,2) - t915 * mrSges(5,3);
t899 = t955 * mrSges(5,1) - t916 * mrSges(5,3);
t827 = -m(5) * t974 + t880 * mrSges(5,1) + t881 * mrSges(5,2) + t915 * t898 + t916 * t899 + t828;
t904 = Ifges(4,5) * t942 + Ifges(4,6) * t941 + Ifges(4,3) * t955;
t906 = Ifges(4,1) * t942 + Ifges(4,4) * t941 + Ifges(4,5) * t955;
t1005 = -mrSges(5,3) - mrSges(6,2);
t979 = m(6) * t840 + t940 * mrSges(6,3) + t955 * t900 + t984;
t888 = t915 * mrSges(6,1) - t916 * mrSges(6,3);
t995 = -t915 * mrSges(5,1) - t916 * mrSges(5,2) - t888;
t816 = m(5) * t843 - t940 * mrSges(5,2) + t1005 * t880 - t955 * t899 + t995 * t915 + t979;
t976 = -m(6) * t841 + t940 * mrSges(6,1) + t955 * t901 - t821;
t818 = m(5) * t842 + t940 * mrSges(5,1) + t1005 * t881 + t955 * t898 + t995 * t916 + t976;
t985 = t1001 * t816 - t960 * t818;
t793 = mrSges(4,1) * t978 + mrSges(4,3) * t866 + Ifges(4,4) * t913 + Ifges(4,2) * t912 + Ifges(4,6) * t940 - pkin(3) * t827 + qJ(4) * t985 + t1001 * t808 + t960 * t809 - t942 * t904 + t955 * t906;
t813 = t1001 * t818 + t960 * t816;
t905 = Ifges(4,4) * t942 + Ifges(4,2) * t941 + Ifges(4,6) * t955;
t794 = -mrSges(4,2) * t978 - mrSges(4,3) * t865 + Ifges(4,1) * t913 + Ifges(4,4) * t912 + Ifges(4,5) * t940 - qJ(4) * t813 + t1001 * t809 - t960 * t808 + t941 * t904 - t955 * t905;
t917 = -t941 * mrSges(4,1) + t942 * mrSges(4,2);
t918 = -t955 * mrSges(4,2) + t941 * mrSges(4,3);
t811 = m(4) * t865 + t940 * mrSges(4,1) - t913 * mrSges(4,3) - t942 * t917 + t955 * t918 + t813;
t920 = t955 * mrSges(4,1) - t942 * mrSges(4,3);
t812 = m(4) * t866 - t940 * mrSges(4,2) + t912 * mrSges(4,3) + t941 * t917 - t955 * t920 + t985;
t807 = -t962 * t811 + t966 * t812;
t826 = m(4) * t978 + t912 * mrSges(4,1) - t913 * mrSges(4,2) + t941 * t918 - t942 * t920 - t827;
t929 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t963 + Ifges(3,2) * t967) * qJD(1);
t930 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t963 + Ifges(3,4) * t967) * qJD(1);
t1008 = mrSges(3,1) * t921 - mrSges(3,2) * t922 + Ifges(3,5) * t945 + Ifges(3,6) * t946 + Ifges(3,3) * qJDD(2) + pkin(2) * t826 + pkin(8) * t807 + t966 * t793 + t962 * t794 + (t963 * t929 - t967 * t930) * qJD(1);
t820 = t881 * mrSges(6,2) + t916 * t888 - t976;
t975 = -mrSges(7,1) * t833 + mrSges(7,2) * t834 - Ifges(7,5) * t851 - Ifges(7,6) * t850 - Ifges(7,3) * t936 - t886 * t859 + t885 * t860;
t1007 = t1002 * t880 + t1003 * t881 + t996 * t915 - t998 * t916 + (Ifges(4,3) - t1010) * t940 + mrSges(4,1) * t865 + mrSges(5,1) * t842 - mrSges(6,1) * t841 - mrSges(4,2) * t866 - mrSges(5,2) * t843 + mrSges(6,3) * t840 + Ifges(4,5) * t913 + Ifges(4,6) * t912 + pkin(3) * t813 - pkin(4) * t820 - pkin(5) * t821 + qJ(5) * (-t880 * mrSges(6,2) - t915 * t888 + t979) + t942 * t905 - t941 * t906 + t975;
t943 = (-mrSges(3,1) * t967 + mrSges(3,2) * t963) * qJD(1);
t948 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t994;
t805 = m(3) * t922 - qJDD(2) * mrSges(3,2) + t946 * mrSges(3,3) - qJD(2) * t948 + t943 * t992 + t807;
t949 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t992;
t825 = m(3) * t921 + qJDD(2) * mrSges(3,1) - t945 * mrSges(3,3) + qJD(2) * t949 - t943 * t994 + t826;
t986 = t967 * t805 - t963 * t825;
t797 = m(2) * t952 - t970 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t986;
t806 = t966 * t811 + t962 * t812;
t973 = -m(3) * t931 + t946 * mrSges(3,1) - t945 * mrSges(3,2) - t948 * t994 + t949 * t992 - t806;
t801 = m(2) * t951 + qJDD(1) * mrSges(2,1) - t970 * mrSges(2,2) + t973;
t999 = t964 * t797 + t968 * t801;
t799 = t963 * t805 + t967 * t825;
t987 = t968 * t797 - t964 * t801;
t928 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t963 + Ifges(3,6) * t967) * qJD(1);
t790 = mrSges(3,2) * t931 - mrSges(3,3) * t921 + Ifges(3,1) * t945 + Ifges(3,4) * t946 + Ifges(3,5) * qJDD(2) - pkin(8) * t806 - qJD(2) * t929 - t962 * t793 + t966 * t794 + t928 * t992;
t792 = -mrSges(3,1) * t931 + mrSges(3,3) * t922 + Ifges(3,4) * t945 + Ifges(3,2) * t946 + Ifges(3,6) * qJDD(2) - pkin(2) * t806 + qJD(2) * t930 - t928 * t994 - t1007;
t977 = mrSges(2,1) * t951 - mrSges(2,2) * t952 + Ifges(2,3) * qJDD(1) + pkin(1) * t973 + pkin(7) * t986 + t963 * t790 + t967 * t792;
t788 = mrSges(2,1) * g(3) + mrSges(2,3) * t952 + t970 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t799 - t1008;
t787 = -mrSges(2,2) * g(3) - mrSges(2,3) * t951 + Ifges(2,5) * qJDD(1) - t970 * Ifges(2,6) - pkin(7) * t799 + t967 * t790 - t963 * t792;
t1 = [-m(1) * g(1) + t987; -m(1) * g(2) + t999; (-m(1) - m(2)) * g(3) + t799; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t999 + t968 * t787 - t964 * t788; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t987 + t964 * t787 + t968 * t788; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t977; t977; t1008; t1007; t827; t820; -t975;];
tauJB  = t1;
