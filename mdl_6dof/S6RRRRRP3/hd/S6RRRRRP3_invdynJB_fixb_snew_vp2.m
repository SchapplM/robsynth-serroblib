% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:39:32
% EndTime: 2019-05-08 04:40:02
% DurationCPUTime: 18.29s
% Computational Cost: add. (295277->367), mult. (588098->446), div. (0->0), fcn. (428180->10), ass. (0->146)
t1001 = Ifges(6,4) + Ifges(7,4);
t1011 = Ifges(6,2) + Ifges(7,2);
t1007 = Ifges(6,6) + Ifges(7,6);
t966 = sin(qJ(3));
t967 = sin(qJ(2));
t971 = cos(qJ(3));
t972 = cos(qJ(2));
t940 = (t966 * t972 + t967 * t971) * qJD(1);
t961 = qJD(2) + qJD(3);
t965 = sin(qJ(4));
t970 = cos(qJ(4));
t926 = -t940 * t965 + t961 * t970;
t927 = t940 * t970 + t961 * t965;
t964 = sin(qJ(5));
t969 = cos(qJ(5));
t895 = t926 * t969 - t927 * t964;
t896 = t926 * t964 + t927 * t969;
t995 = qJD(1) * t972;
t996 = qJD(1) * t967;
t939 = -t966 * t996 + t971 * t995;
t935 = qJD(4) - t939;
t933 = qJD(5) + t935;
t1005 = t1001 * t896 + t1007 * t933 + t1011 * t895;
t1006 = Ifges(6,3) + Ifges(7,3);
t1008 = Ifges(6,5) + Ifges(7,5);
t994 = qJD(1) * qJD(2);
t948 = qJDD(1) * t967 + t972 * t994;
t949 = qJDD(1) * t972 - t967 * t994;
t914 = qJD(3) * t939 + t948 * t971 + t949 * t966;
t960 = qJDD(2) + qJDD(3);
t882 = -qJD(4) * t927 - t914 * t965 + t960 * t970;
t883 = qJD(4) * t926 + t914 * t970 + t960 * t965;
t856 = qJD(5) * t895 + t882 * t964 + t883 * t969;
t875 = -mrSges(7,1) * t895 + mrSges(7,2) * t896;
t913 = -t940 * qJD(3) - t966 * t948 + t949 * t971;
t953 = qJD(2) * pkin(2) - pkin(8) * t996;
t963 = t972 ^ 2;
t974 = qJD(1) ^ 2;
t968 = sin(qJ(1));
t973 = cos(qJ(1));
t954 = t968 * g(1) - t973 * g(2);
t984 = -qJDD(1) * pkin(1) - t954;
t915 = -t949 * pkin(2) + t953 * t996 + (-pkin(8) * t963 - pkin(7)) * t974 + t984;
t861 = (-t939 * t961 - t914) * pkin(9) + (t940 * t961 - t913) * pkin(3) + t915;
t955 = -g(1) * t973 - g(2) * t968;
t942 = -pkin(1) * t974 + qJDD(1) * pkin(7) + t955;
t1000 = t967 * t942;
t1002 = pkin(2) * t974;
t902 = qJDD(2) * pkin(2) - t948 * pkin(8) - t1000 + (pkin(8) * t994 + t1002 * t967 - g(3)) * t972;
t929 = -g(3) * t967 + t972 * t942;
t903 = pkin(8) * t949 - qJD(2) * t953 - t1002 * t963 + t929;
t879 = t966 * t902 + t971 * t903;
t924 = -pkin(3) * t939 - pkin(9) * t940;
t959 = t961 ^ 2;
t865 = -pkin(3) * t959 + pkin(9) * t960 + t924 * t939 + t879;
t844 = t970 * t861 - t965 * t865;
t912 = qJDD(4) - t913;
t841 = (t926 * t935 - t883) * pkin(10) + (t926 * t927 + t912) * pkin(4) + t844;
t845 = t965 * t861 + t970 * t865;
t918 = pkin(4) * t935 - pkin(10) * t927;
t925 = t926 ^ 2;
t843 = -pkin(4) * t925 + pkin(10) * t882 - t918 * t935 + t845;
t835 = t969 * t841 - t964 * t843;
t907 = qJDD(5) + t912;
t830 = -0.2e1 * qJD(6) * t896 + (t895 * t933 - t856) * qJ(6) + (t895 * t896 + t907) * pkin(5) + t835;
t884 = -mrSges(7,2) * t933 + mrSges(7,3) * t895;
t993 = m(7) * t830 + t907 * mrSges(7,1) + t933 * t884;
t827 = -t856 * mrSges(7,3) - t896 * t875 + t993;
t836 = t964 * t841 + t969 * t843;
t855 = -qJD(5) * t896 + t882 * t969 - t883 * t964;
t886 = pkin(5) * t933 - qJ(6) * t896;
t894 = t895 ^ 2;
t832 = -pkin(5) * t894 + qJ(6) * t855 + 0.2e1 * qJD(6) * t895 - t886 * t933 + t836;
t1009 = Ifges(6,1) + Ifges(7,1);
t997 = -t1001 * t895 - t1008 * t933 - t1009 * t896;
t1010 = mrSges(6,1) * t835 + mrSges(7,1) * t830 - mrSges(6,2) * t836 - mrSges(7,2) * t832 + pkin(5) * t827 + t1005 * t896 + t1006 * t907 + t1007 * t855 + t1008 * t856 + t997 * t895;
t876 = -mrSges(6,1) * t895 + mrSges(6,2) * t896;
t885 = -mrSges(6,2) * t933 + mrSges(6,3) * t895;
t818 = m(6) * t835 + t907 * mrSges(6,1) + t933 * t885 + (-t875 - t876) * t896 + (-mrSges(6,3) - mrSges(7,3)) * t856 + t993;
t887 = mrSges(7,1) * t933 - mrSges(7,3) * t896;
t888 = mrSges(6,1) * t933 - mrSges(6,3) * t896;
t992 = m(7) * t832 + t855 * mrSges(7,3) + t895 * t875;
t821 = m(6) * t836 + t855 * mrSges(6,3) + t895 * t876 + (-t887 - t888) * t933 + (-mrSges(6,2) - mrSges(7,2)) * t907 + t992;
t816 = t969 * t818 + t964 * t821;
t890 = Ifges(5,4) * t927 + Ifges(5,2) * t926 + Ifges(5,6) * t935;
t891 = Ifges(5,1) * t927 + Ifges(5,4) * t926 + Ifges(5,5) * t935;
t1004 = mrSges(5,1) * t844 - mrSges(5,2) * t845 + Ifges(5,5) * t883 + Ifges(5,6) * t882 + Ifges(5,3) * t912 + pkin(4) * t816 + t927 * t890 - t926 * t891 + t1010;
t923 = -mrSges(4,1) * t939 + mrSges(4,2) * t940;
t931 = mrSges(4,1) * t961 - mrSges(4,3) * t940;
t900 = -mrSges(5,1) * t926 + mrSges(5,2) * t927;
t916 = -mrSges(5,2) * t935 + mrSges(5,3) * t926;
t813 = m(5) * t844 + mrSges(5,1) * t912 - mrSges(5,3) * t883 - t900 * t927 + t916 * t935 + t816;
t917 = mrSges(5,1) * t935 - mrSges(5,3) * t927;
t986 = -t818 * t964 + t969 * t821;
t814 = m(5) * t845 - mrSges(5,2) * t912 + mrSges(5,3) * t882 + t900 * t926 - t917 * t935 + t986;
t987 = -t813 * t965 + t970 * t814;
t805 = m(4) * t879 - mrSges(4,2) * t960 + mrSges(4,3) * t913 + t923 * t939 - t931 * t961 + t987;
t878 = t902 * t971 - t966 * t903;
t930 = -mrSges(4,2) * t961 + mrSges(4,3) * t939;
t864 = -pkin(3) * t960 - pkin(9) * t959 + t940 * t924 - t878;
t846 = -pkin(4) * t882 - pkin(10) * t925 + t927 * t918 + t864;
t838 = -pkin(5) * t855 - qJ(6) * t894 + t886 * t896 + qJDD(6) + t846;
t833 = m(7) * t838 - t855 * mrSges(7,1) + t856 * mrSges(7,2) - t895 * t884 + t896 * t887;
t980 = m(6) * t846 - t855 * mrSges(6,1) + mrSges(6,2) * t856 - t895 * t885 + t888 * t896 + t833;
t977 = -m(5) * t864 + t882 * mrSges(5,1) - mrSges(5,2) * t883 + t926 * t916 - t917 * t927 - t980;
t823 = m(4) * t878 + mrSges(4,1) * t960 - mrSges(4,3) * t914 - t923 * t940 + t930 * t961 + t977;
t799 = t966 * t805 + t971 * t823;
t928 = -t972 * g(3) - t1000;
t937 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t967 + Ifges(3,2) * t972) * qJD(1);
t938 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t967 + Ifges(3,4) * t972) * qJD(1);
t998 = -t1006 * t933 - t1007 * t895 - t1008 * t896;
t809 = -mrSges(6,1) * t846 + mrSges(6,3) * t836 - mrSges(7,1) * t838 + mrSges(7,3) * t832 - pkin(5) * t833 + qJ(6) * t992 + (-qJ(6) * t887 - t997) * t933 + (-mrSges(7,2) * qJ(6) + t1007) * t907 + t998 * t896 + t1001 * t856 + t1011 * t855;
t815 = mrSges(6,2) * t846 + mrSges(7,2) * t838 - mrSges(6,3) * t835 - mrSges(7,3) * t830 - qJ(6) * t827 + t1001 * t855 - t1005 * t933 + t1008 * t907 + t1009 * t856 - t998 * t895;
t889 = Ifges(5,5) * t927 + Ifges(5,6) * t926 + Ifges(5,3) * t935;
t791 = -mrSges(5,1) * t864 + mrSges(5,3) * t845 + Ifges(5,4) * t883 + Ifges(5,2) * t882 + Ifges(5,6) * t912 - pkin(4) * t980 + pkin(10) * t986 + t969 * t809 + t964 * t815 - t927 * t889 + t935 * t891;
t793 = mrSges(5,2) * t864 - mrSges(5,3) * t844 + Ifges(5,1) * t883 + Ifges(5,4) * t882 + Ifges(5,5) * t912 - pkin(10) * t816 - t809 * t964 + t815 * t969 + t889 * t926 - t890 * t935;
t920 = Ifges(4,4) * t940 + Ifges(4,2) * t939 + Ifges(4,6) * t961;
t921 = Ifges(4,1) * t940 + Ifges(4,4) * t939 + Ifges(4,5) * t961;
t981 = -mrSges(4,1) * t878 + mrSges(4,2) * t879 - Ifges(4,5) * t914 - Ifges(4,6) * t913 - Ifges(4,3) * t960 - pkin(3) * t977 - pkin(9) * t987 - t970 * t791 - t965 * t793 - t940 * t920 + t939 * t921;
t1003 = mrSges(3,1) * t928 - mrSges(3,2) * t929 + Ifges(3,5) * t948 + Ifges(3,6) * t949 + Ifges(3,3) * qJDD(2) + pkin(2) * t799 + (t937 * t967 - t938 * t972) * qJD(1) - t981;
t947 = (-mrSges(3,1) * t972 + mrSges(3,2) * t967) * qJD(1);
t952 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t995;
t797 = m(3) * t928 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t948 + qJD(2) * t952 - t947 * t996 + t799;
t951 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t996;
t988 = t971 * t805 - t823 * t966;
t798 = m(3) * t929 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t949 - qJD(2) * t951 + t947 * t995 + t988;
t989 = -t797 * t967 + t972 * t798;
t786 = m(2) * t955 - mrSges(2,1) * t974 - qJDD(1) * mrSges(2,2) + t989;
t941 = -t974 * pkin(7) + t984;
t807 = t970 * t813 + t965 * t814;
t982 = m(4) * t915 - t913 * mrSges(4,1) + mrSges(4,2) * t914 - t939 * t930 + t931 * t940 + t807;
t978 = -m(3) * t941 + t949 * mrSges(3,1) - mrSges(3,2) * t948 - t951 * t996 + t952 * t995 - t982;
t801 = m(2) * t954 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t974 + t978;
t999 = t968 * t786 + t973 * t801;
t788 = t972 * t797 + t967 * t798;
t990 = t973 * t786 - t801 * t968;
t919 = Ifges(4,5) * t940 + Ifges(4,6) * t939 + Ifges(4,3) * t961;
t783 = mrSges(4,2) * t915 - mrSges(4,3) * t878 + Ifges(4,1) * t914 + Ifges(4,4) * t913 + Ifges(4,5) * t960 - pkin(9) * t807 - t791 * t965 + t793 * t970 + t919 * t939 - t920 * t961;
t789 = -mrSges(4,1) * t915 + mrSges(4,3) * t879 + Ifges(4,4) * t914 + Ifges(4,2) * t913 + Ifges(4,6) * t960 - pkin(3) * t807 - t940 * t919 + t961 * t921 - t1004;
t936 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t967 + Ifges(3,6) * t972) * qJD(1);
t779 = -mrSges(3,1) * t941 + mrSges(3,3) * t929 + Ifges(3,4) * t948 + Ifges(3,2) * t949 + Ifges(3,6) * qJDD(2) - pkin(2) * t982 + pkin(8) * t988 + qJD(2) * t938 + t966 * t783 + t971 * t789 - t936 * t996;
t782 = mrSges(3,2) * t941 - mrSges(3,3) * t928 + Ifges(3,1) * t948 + Ifges(3,4) * t949 + Ifges(3,5) * qJDD(2) - pkin(8) * t799 - qJD(2) * t937 + t783 * t971 - t789 * t966 + t936 * t995;
t983 = mrSges(2,1) * t954 - mrSges(2,2) * t955 + Ifges(2,3) * qJDD(1) + pkin(1) * t978 + pkin(7) * t989 + t972 * t779 + t967 * t782;
t780 = mrSges(2,1) * g(3) + mrSges(2,3) * t955 + t974 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t788 - t1003;
t777 = -mrSges(2,2) * g(3) - mrSges(2,3) * t954 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t974 - pkin(7) * t788 - t779 * t967 + t782 * t972;
t1 = [-m(1) * g(1) + t990; -m(1) * g(2) + t999; (-m(1) - m(2)) * g(3) + t788; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t999 + t973 * t777 - t968 * t780; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t990 + t968 * t777 + t973 * t780; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t983; t983; t1003; -t981; t1004; t1010; t833;];
tauJB  = t1;
