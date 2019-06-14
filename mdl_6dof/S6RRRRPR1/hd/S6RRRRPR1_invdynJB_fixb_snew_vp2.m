% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 19:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:28:23
% EndTime: 2019-05-07 19:29:01
% DurationCPUTime: 39.55s
% Computational Cost: add. (667207->389), mult. (1544735->490), div. (0->0), fcn. (1180691->12), ass. (0->156)
t963 = sin(qJ(2));
t968 = cos(qJ(2));
t989 = qJD(1) * qJD(2);
t939 = qJDD(1) * t963 + t968 * t989;
t964 = sin(qJ(1));
t969 = cos(qJ(1));
t946 = -g(1) * t969 - g(2) * t964;
t970 = qJD(1) ^ 2;
t934 = -pkin(1) * t970 + qJDD(1) * pkin(7) + t946;
t993 = t963 * t934;
t994 = pkin(2) * t970;
t901 = qJDD(2) * pkin(2) - t939 * pkin(8) - t993 + (pkin(8) * t989 + t963 * t994 - g(3)) * t968;
t922 = -g(3) * t963 + t968 * t934;
t940 = qJDD(1) * t968 - t963 * t989;
t991 = qJD(1) * t963;
t944 = qJD(2) * pkin(2) - pkin(8) * t991;
t957 = t968 ^ 2;
t902 = pkin(8) * t940 - qJD(2) * t944 - t957 * t994 + t922;
t962 = sin(qJ(3));
t967 = cos(qJ(3));
t880 = t967 * t901 - t962 * t902;
t931 = (-t962 * t963 + t967 * t968) * qJD(1);
t907 = qJD(3) * t931 + t939 * t967 + t940 * t962;
t932 = (t962 * t968 + t963 * t967) * qJD(1);
t954 = qJDD(2) + qJDD(3);
t955 = qJD(2) + qJD(3);
t857 = (t931 * t955 - t907) * pkin(9) + (t931 * t932 + t954) * pkin(3) + t880;
t881 = t962 * t901 + t967 * t902;
t906 = -qJD(3) * t932 - t939 * t962 + t940 * t967;
t925 = pkin(3) * t955 - pkin(9) * t932;
t927 = t931 ^ 2;
t859 = -pkin(3) * t927 + pkin(9) * t906 - t925 * t955 + t881;
t961 = sin(qJ(4));
t966 = cos(qJ(4));
t836 = t966 * t857 - t961 * t859;
t918 = t931 * t966 - t932 * t961;
t877 = qJD(4) * t918 + t906 * t961 + t907 * t966;
t919 = t931 * t961 + t932 * t966;
t951 = qJDD(4) + t954;
t952 = qJD(4) + t955;
t832 = (t918 * t952 - t877) * qJ(5) + (t918 * t919 + t951) * pkin(4) + t836;
t837 = t961 * t857 + t966 * t859;
t876 = -qJD(4) * t919 + t906 * t966 - t907 * t961;
t910 = pkin(4) * t952 - qJ(5) * t919;
t917 = t918 ^ 2;
t834 = -pkin(4) * t917 + qJ(5) * t876 - t910 * t952 + t837;
t958 = sin(pkin(11));
t959 = cos(pkin(11));
t894 = t918 * t959 - t919 * t958;
t995 = 2 * qJD(5);
t829 = t958 * t832 + t959 * t834 + t894 * t995;
t848 = t876 * t959 - t877 * t958;
t895 = t918 * t958 + t919 * t959;
t868 = -mrSges(6,1) * t894 + mrSges(6,2) * t895;
t885 = mrSges(6,1) * t952 - mrSges(6,3) * t895;
t869 = -pkin(5) * t894 - pkin(10) * t895;
t950 = t952 ^ 2;
t826 = -pkin(5) * t950 + pkin(10) * t951 + t869 * t894 + t829;
t945 = t964 * g(1) - t969 * g(2);
t980 = -qJDD(1) * pkin(1) - t945;
t908 = -t940 * pkin(2) + t944 * t991 + (-pkin(8) * t957 - pkin(7)) * t970 + t980;
t871 = -t906 * pkin(3) - t927 * pkin(9) + t932 * t925 + t908;
t839 = -t876 * pkin(4) - t917 * qJ(5) + t919 * t910 + qJDD(5) + t871;
t849 = t876 * t958 + t877 * t959;
t830 = t839 + (-t894 * t952 - t849) * pkin(10) + (t895 * t952 - t848) * pkin(5);
t960 = sin(qJ(6));
t965 = cos(qJ(6));
t823 = -t826 * t960 + t830 * t965;
t882 = -t895 * t960 + t952 * t965;
t842 = qJD(6) * t882 + t849 * t965 + t951 * t960;
t847 = qJDD(6) - t848;
t883 = t895 * t965 + t952 * t960;
t860 = -mrSges(7,1) * t882 + mrSges(7,2) * t883;
t893 = qJD(6) - t894;
t861 = -mrSges(7,2) * t893 + mrSges(7,3) * t882;
t819 = m(7) * t823 + mrSges(7,1) * t847 - mrSges(7,3) * t842 - t860 * t883 + t861 * t893;
t824 = t826 * t965 + t830 * t960;
t841 = -qJD(6) * t883 - t849 * t960 + t951 * t965;
t862 = mrSges(7,1) * t893 - mrSges(7,3) * t883;
t820 = m(7) * t824 - mrSges(7,2) * t847 + mrSges(7,3) * t841 + t860 * t882 - t862 * t893;
t983 = -t819 * t960 + t965 * t820;
t805 = m(6) * t829 - mrSges(6,2) * t951 + mrSges(6,3) * t848 + t868 * t894 - t885 * t952 + t983;
t982 = -t959 * t832 + t958 * t834;
t828 = -0.2e1 * qJD(5) * t895 - t982;
t884 = -mrSges(6,2) * t952 + mrSges(6,3) * t894;
t825 = -t951 * pkin(5) - t950 * pkin(10) + (t995 + t869) * t895 + t982;
t978 = -m(7) * t825 + t841 * mrSges(7,1) - mrSges(7,2) * t842 + t882 * t861 - t862 * t883;
t815 = m(6) * t828 + mrSges(6,1) * t951 - mrSges(6,3) * t849 - t868 * t895 + t884 * t952 + t978;
t799 = t958 * t805 + t959 * t815;
t896 = -mrSges(5,1) * t918 + mrSges(5,2) * t919;
t909 = -mrSges(5,2) * t952 + mrSges(5,3) * t918;
t796 = m(5) * t836 + mrSges(5,1) * t951 - mrSges(5,3) * t877 - t896 * t919 + t909 * t952 + t799;
t911 = mrSges(5,1) * t952 - mrSges(5,3) * t919;
t984 = t959 * t805 - t815 * t958;
t797 = m(5) * t837 - mrSges(5,2) * t951 + mrSges(5,3) * t876 + t896 * t918 - t911 * t952 + t984;
t790 = t966 * t796 + t961 * t797;
t920 = -mrSges(4,1) * t931 + mrSges(4,2) * t932;
t923 = -mrSges(4,2) * t955 + mrSges(4,3) * t931;
t787 = m(4) * t880 + mrSges(4,1) * t954 - mrSges(4,3) * t907 - t920 * t932 + t923 * t955 + t790;
t924 = mrSges(4,1) * t955 - mrSges(4,3) * t932;
t985 = -t796 * t961 + t966 * t797;
t788 = m(4) * t881 - mrSges(4,2) * t954 + mrSges(4,3) * t906 + t920 * t931 - t924 * t955 + t985;
t782 = t967 * t787 + t962 * t788;
t921 = -t968 * g(3) - t993;
t929 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t963 + Ifges(3,2) * t968) * qJD(1);
t930 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t963 + Ifges(3,4) * t968) * qJD(1);
t913 = Ifges(4,4) * t932 + Ifges(4,2) * t931 + Ifges(4,6) * t955;
t914 = Ifges(4,1) * t932 + Ifges(4,4) * t931 + Ifges(4,5) * t955;
t850 = Ifges(7,5) * t883 + Ifges(7,6) * t882 + Ifges(7,3) * t893;
t852 = Ifges(7,1) * t883 + Ifges(7,4) * t882 + Ifges(7,5) * t893;
t812 = -mrSges(7,1) * t825 + mrSges(7,3) * t824 + Ifges(7,4) * t842 + Ifges(7,2) * t841 + Ifges(7,6) * t847 - t850 * t883 + t852 * t893;
t851 = Ifges(7,4) * t883 + Ifges(7,2) * t882 + Ifges(7,6) * t893;
t813 = mrSges(7,2) * t825 - mrSges(7,3) * t823 + Ifges(7,1) * t842 + Ifges(7,4) * t841 + Ifges(7,5) * t847 + t850 * t882 - t851 * t893;
t864 = Ifges(6,4) * t895 + Ifges(6,2) * t894 + Ifges(6,6) * t952;
t865 = Ifges(6,1) * t895 + Ifges(6,4) * t894 + Ifges(6,5) * t952;
t887 = Ifges(5,4) * t919 + Ifges(5,2) * t918 + Ifges(5,6) * t952;
t888 = Ifges(5,1) * t919 + Ifges(5,4) * t918 + Ifges(5,5) * t952;
t974 = -mrSges(5,1) * t836 - mrSges(6,1) * t828 + mrSges(5,2) * t837 + mrSges(6,2) * t829 - Ifges(6,6) * t848 - pkin(4) * t799 - pkin(5) * t978 - pkin(10) * t983 - t965 * t812 - t960 * t813 + t894 * t865 + t918 * t888 - Ifges(6,5) * t849 - t895 * t864 - Ifges(5,6) * t876 - Ifges(5,5) * t877 - t919 * t887 + (-Ifges(5,3) - Ifges(6,3)) * t951;
t972 = mrSges(4,1) * t880 - mrSges(4,2) * t881 + Ifges(4,5) * t907 + Ifges(4,6) * t906 + Ifges(4,3) * t954 + pkin(3) * t790 + t932 * t913 - t931 * t914 - t974;
t996 = mrSges(3,1) * t921 - mrSges(3,2) * t922 + Ifges(3,5) * t939 + Ifges(3,6) * t940 + Ifges(3,3) * qJDD(2) + pkin(2) * t782 + (t929 * t963 - t930 * t968) * qJD(1) + t972;
t938 = (-mrSges(3,1) * t968 + mrSges(3,2) * t963) * qJD(1);
t990 = qJD(1) * t968;
t943 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t990;
t780 = m(3) * t921 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t939 + qJD(2) * t943 - t938 * t991 + t782;
t942 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t991;
t986 = -t787 * t962 + t967 * t788;
t781 = m(3) * t922 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t940 - qJD(2) * t942 + t938 * t990 + t986;
t987 = -t780 * t963 + t968 * t781;
t773 = m(2) * t946 - mrSges(2,1) * t970 - qJDD(1) * mrSges(2,2) + t987;
t933 = -t970 * pkin(7) + t980;
t808 = t965 * t819 + t960 * t820;
t806 = m(6) * t839 - t848 * mrSges(6,1) + t849 * mrSges(6,2) - t894 * t884 + t895 * t885 + t808;
t977 = m(5) * t871 - t876 * mrSges(5,1) + t877 * mrSges(5,2) - t918 * t909 + t919 * t911 + t806;
t975 = m(4) * t908 - t906 * mrSges(4,1) + t907 * mrSges(4,2) - t931 * t923 + t932 * t924 + t977;
t973 = -m(3) * t933 + t940 * mrSges(3,1) - t939 * mrSges(3,2) - t942 * t991 + t943 * t990 - t975;
t801 = m(2) * t945 + qJDD(1) * mrSges(2,1) - t970 * mrSges(2,2) + t973;
t992 = t964 * t773 + t969 * t801;
t775 = t968 * t780 + t963 * t781;
t988 = t969 * t773 - t801 * t964;
t863 = Ifges(6,5) * t895 + Ifges(6,6) * t894 + Ifges(6,3) * t952;
t791 = mrSges(6,2) * t839 - mrSges(6,3) * t828 + Ifges(6,1) * t849 + Ifges(6,4) * t848 + Ifges(6,5) * t951 - pkin(10) * t808 - t812 * t960 + t813 * t965 + t863 * t894 - t864 * t952;
t976 = mrSges(7,1) * t823 - mrSges(7,2) * t824 + Ifges(7,5) * t842 + Ifges(7,6) * t841 + Ifges(7,3) * t847 + t851 * t883 - t852 * t882;
t792 = -mrSges(6,1) * t839 + mrSges(6,3) * t829 + Ifges(6,4) * t849 + Ifges(6,2) * t848 + Ifges(6,6) * t951 - pkin(5) * t808 - t863 * t895 + t865 * t952 - t976;
t886 = Ifges(5,5) * t919 + Ifges(5,6) * t918 + Ifges(5,3) * t952;
t776 = -mrSges(5,1) * t871 + mrSges(5,3) * t837 + Ifges(5,4) * t877 + Ifges(5,2) * t876 + Ifges(5,6) * t951 - pkin(4) * t806 + qJ(5) * t984 + t958 * t791 + t959 * t792 - t919 * t886 + t952 * t888;
t783 = mrSges(5,2) * t871 - mrSges(5,3) * t836 + Ifges(5,1) * t877 + Ifges(5,4) * t876 + Ifges(5,5) * t951 - qJ(5) * t799 + t791 * t959 - t792 * t958 + t886 * t918 - t887 * t952;
t912 = Ifges(4,5) * t932 + Ifges(4,6) * t931 + Ifges(4,3) * t955;
t769 = -mrSges(4,1) * t908 + mrSges(4,3) * t881 + Ifges(4,4) * t907 + Ifges(4,2) * t906 + Ifges(4,6) * t954 - pkin(3) * t977 + pkin(9) * t985 + t966 * t776 + t961 * t783 - t932 * t912 + t955 * t914;
t770 = mrSges(4,2) * t908 - mrSges(4,3) * t880 + Ifges(4,1) * t907 + Ifges(4,4) * t906 + Ifges(4,5) * t954 - pkin(9) * t790 - t776 * t961 + t783 * t966 + t912 * t931 - t913 * t955;
t928 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t963 + Ifges(3,6) * t968) * qJD(1);
t765 = -mrSges(3,1) * t933 + mrSges(3,3) * t922 + Ifges(3,4) * t939 + Ifges(3,2) * t940 + Ifges(3,6) * qJDD(2) - pkin(2) * t975 + pkin(8) * t986 + qJD(2) * t930 + t967 * t769 + t962 * t770 - t928 * t991;
t767 = mrSges(3,2) * t933 - mrSges(3,3) * t921 + Ifges(3,1) * t939 + Ifges(3,4) * t940 + Ifges(3,5) * qJDD(2) - pkin(8) * t782 - qJD(2) * t929 - t769 * t962 + t770 * t967 + t928 * t990;
t979 = mrSges(2,1) * t945 - mrSges(2,2) * t946 + Ifges(2,3) * qJDD(1) + pkin(1) * t973 + pkin(7) * t987 + t968 * t765 + t963 * t767;
t768 = mrSges(2,1) * g(3) + mrSges(2,3) * t946 + t970 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t775 - t996;
t763 = -mrSges(2,2) * g(3) - mrSges(2,3) * t945 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t970 - pkin(7) * t775 - t765 * t963 + t767 * t968;
t1 = [-m(1) * g(1) + t988; -m(1) * g(2) + t992; (-m(1) - m(2)) * g(3) + t775; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t992 + t969 * t763 - t964 * t768; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t988 + t964 * t763 + t969 * t768; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t979; t979; t996; t972; -t974; t806; t976;];
tauJB  = t1;
