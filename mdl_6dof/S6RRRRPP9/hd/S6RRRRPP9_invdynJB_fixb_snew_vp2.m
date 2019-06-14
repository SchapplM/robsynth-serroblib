% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPP9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:13:10
% EndTime: 2019-05-07 19:13:33
% DurationCPUTime: 13.40s
% Computational Cost: add. (209836->356), mult. (444882->433), div. (0->0), fcn. (343766->10), ass. (0->149)
t1000 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t978 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t977 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t999 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t976 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t998 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t946 = sin(pkin(6));
t950 = sin(qJ(2));
t953 = cos(qJ(2));
t979 = qJD(1) * qJD(2);
t932 = (-qJDD(1) * t953 + t950 * t979) * t946;
t981 = qJD(1) * t946;
t930 = (-pkin(2) * t953 - pkin(9) * t950) * t981;
t947 = cos(pkin(6));
t942 = t947 * qJD(1) + qJD(2);
t940 = t942 ^ 2;
t941 = t947 * qJDD(1) + qJDD(2);
t980 = qJD(1) * t953;
t951 = sin(qJ(1));
t954 = cos(qJ(1));
t938 = t951 * g(1) - t954 * g(2);
t955 = qJD(1) ^ 2;
t992 = pkin(8) * t946;
t927 = qJDD(1) * pkin(1) + t955 * t992 + t938;
t939 = -t954 * g(1) - t951 * g(2);
t928 = -t955 * pkin(1) + qJDD(1) * t992 + t939;
t986 = t947 * t950;
t982 = t927 * t986 + t953 * t928;
t872 = -t940 * pkin(2) + t941 * pkin(9) + (-g(3) * t950 + t930 * t980) * t946 + t982;
t931 = (qJDD(1) * t950 + t953 * t979) * t946;
t991 = t947 * g(3);
t873 = t932 * pkin(2) - t931 * pkin(9) - t991 + (-t927 + (pkin(2) * t950 - pkin(9) * t953) * t942 * qJD(1)) * t946;
t949 = sin(qJ(3));
t952 = cos(qJ(3));
t845 = t952 * t872 + t949 * t873;
t970 = t950 * t981;
t920 = t952 * t942 - t949 * t970;
t921 = t949 * t942 + t952 * t970;
t904 = -t920 * pkin(3) - t921 * pkin(10);
t924 = qJDD(3) + t932;
t969 = t946 * t980;
t937 = qJD(3) - t969;
t936 = t937 ^ 2;
t841 = -t936 * pkin(3) + t924 * pkin(10) + t920 * t904 + t845;
t985 = t947 * t953;
t987 = t946 * t953;
t901 = -g(3) * t987 + t927 * t985 - t950 * t928;
t871 = -t941 * pkin(2) - t940 * pkin(9) + t930 * t970 - t901;
t899 = -t921 * qJD(3) - t949 * t931 + t952 * t941;
t900 = t920 * qJD(3) + t952 * t931 + t949 * t941;
t843 = (-t920 * t937 - t900) * pkin(10) + (t921 * t937 - t899) * pkin(3) + t871;
t948 = sin(qJ(4));
t994 = cos(qJ(4));
t836 = -t948 * t841 + t843 * t994;
t907 = t948 * t921 - t937 * t994;
t908 = t921 * t994 + t948 * t937;
t878 = t907 * pkin(4) - t908 * qJ(5);
t897 = qJDD(4) - t899;
t918 = qJD(4) - t920;
t917 = t918 ^ 2;
t834 = -t897 * pkin(4) - t917 * qJ(5) + t908 * t878 + qJDD(5) - t836;
t852 = -t907 * qJD(4) + t900 * t994 + t948 * t924;
t880 = -t907 * mrSges(6,2) - t908 * mrSges(6,3);
t997 = -m(6) * t834 - t852 * mrSges(6,1) - t908 * t880;
t877 = -t908 * mrSges(7,2) + t907 * mrSges(7,3);
t989 = t907 * t918;
t828 = -0.2e1 * qJD(6) * t918 + (t907 * t908 - t897) * qJ(6) + (t852 + t989) * pkin(5) + t834;
t885 = -t907 * mrSges(7,1) + t918 * mrSges(7,2);
t966 = -m(7) * t828 + t897 * mrSges(7,3) + t918 * t885;
t826 = t852 * mrSges(7,1) + t908 * t877 - t966;
t884 = t907 * mrSges(6,1) - t918 * mrSges(6,3);
t823 = t897 * mrSges(6,2) + t918 * t884 + t826 - t997;
t851 = t908 * qJD(4) + t948 * t900 - t924 * t994;
t882 = t908 * pkin(5) - t918 * qJ(6);
t906 = t907 ^ 2;
t837 = t994 * t841 + t948 * t843;
t961 = -t917 * pkin(4) + t897 * qJ(5) - t907 * t878 + t837;
t830 = -t851 * pkin(5) - t906 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t882) * t918 + t961;
t995 = -2 * qJD(5);
t833 = t918 * t995 - t961;
t886 = t908 * mrSges(6,1) + t918 * mrSges(6,2);
t883 = t908 * mrSges(7,1) - t918 * mrSges(7,3);
t974 = m(7) * t830 + t897 * mrSges(7,2) + t918 * t883;
t963 = -m(6) * t833 + t897 * mrSges(6,3) + t918 * t886 + t974;
t971 = t1000 * t908 - t978 * t907 + t977 * t918;
t972 = t907 * t999 + t908 * t978 + t918 * t976;
t983 = -t877 - t880;
t996 = -t851 * t976 + t852 * t977 + t998 * t897 + t907 * t971 + t908 * t972 + mrSges(5,1) * t836 - mrSges(5,2) * t837 + mrSges(6,2) * t834 + mrSges(7,2) * t830 - mrSges(6,3) * t833 - mrSges(7,3) * t828 - pkin(4) * t823 + qJ(5) * (t983 * t907 + (-mrSges(6,1) - mrSges(7,1)) * t851 + t963) - qJ(6) * t826;
t990 = -mrSges(7,1) - mrSges(5,3);
t988 = t946 * t950;
t902 = -g(3) * t988 + t982;
t925 = t942 * mrSges(3,1) - mrSges(3,3) * t970;
t929 = (-mrSges(3,1) * t953 + mrSges(3,2) * t950) * t981;
t879 = t907 * mrSges(5,1) + t908 * mrSges(5,2);
t887 = -t918 * mrSges(5,2) - t907 * mrSges(5,3);
t819 = m(5) * t836 + (-t884 + t887) * t918 + (-t877 - t879) * t908 + (mrSges(5,1) - mrSges(6,2)) * t897 + t990 * t852 + t966 + t997;
t888 = t918 * mrSges(5,1) - t908 * mrSges(5,3);
t821 = m(5) * t837 - t897 * mrSges(5,2) - t918 * t888 + (-t879 + t983) * t907 + (-mrSges(6,1) + t990) * t851 + t963;
t816 = -t948 * t819 + t994 * t821;
t903 = -t920 * mrSges(4,1) + t921 * mrSges(4,2);
t910 = t937 * mrSges(4,1) - t921 * mrSges(4,3);
t814 = m(4) * t845 - t924 * mrSges(4,2) + t899 * mrSges(4,3) + t920 * t903 - t937 * t910 + t816;
t844 = -t949 * t872 + t952 * t873;
t840 = -t924 * pkin(3) - t936 * pkin(10) + t921 * t904 - t844;
t958 = (-t852 + t989) * qJ(5) + t840 + (pkin(4) * t918 + t995) * t908;
t835 = t851 * pkin(4) + t958;
t832 = -t906 * pkin(5) + 0.2e1 * qJD(6) * t907 - t908 * t882 + (pkin(4) + qJ(6)) * t851 + t958;
t965 = m(7) * t832 - t852 * mrSges(7,2) + t851 * mrSges(7,3) - t908 * t883 + t907 * t885;
t962 = -m(6) * t835 + t851 * mrSges(6,2) + t907 * t884 - t965;
t822 = -m(5) * t840 - t851 * mrSges(5,1) - t907 * t887 + (t886 - t888) * t908 + (-mrSges(5,2) + mrSges(6,3)) * t852 + t962;
t909 = -t937 * mrSges(4,2) + t920 * mrSges(4,3);
t818 = m(4) * t844 + t924 * mrSges(4,1) - t900 * mrSges(4,3) - t921 * t903 + t937 * t909 + t822;
t967 = t952 * t814 - t949 * t818;
t804 = m(3) * t902 - t941 * mrSges(3,2) - t932 * mrSges(3,3) - t942 * t925 + t929 * t969 + t967;
t807 = t949 * t814 + t952 * t818;
t914 = -t946 * t927 - t991;
t926 = -t942 * mrSges(3,2) + mrSges(3,3) * t969;
t806 = m(3) * t914 + t932 * mrSges(3,1) + t931 * mrSges(3,2) + (t925 * t950 - t926 * t953) * t981 + t807;
t815 = t819 * t994 + t948 * t821;
t959 = -m(4) * t871 + t899 * mrSges(4,1) - t900 * mrSges(4,2) + t920 * t909 - t921 * t910 - t815;
t811 = m(3) * t901 + t941 * mrSges(3,1) - t931 * mrSges(3,3) + t942 * t926 - t929 * t970 + t959;
t792 = t804 * t986 - t946 * t806 + t811 * t985;
t789 = m(2) * t938 + qJDD(1) * mrSges(2,1) - t955 * mrSges(2,2) + t792;
t799 = t953 * t804 - t950 * t811;
t797 = m(2) * t939 - t955 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t799;
t984 = t954 * t789 + t951 * t797;
t791 = t804 * t988 + t947 * t806 + t811 * t987;
t973 = t907 * t976 - t908 * t977 - t918 * t998;
t968 = -t951 * t789 + t954 * t797;
t825 = -t852 * mrSges(6,3) - t908 * t886 - t962;
t827 = -t851 * mrSges(7,1) - t907 * t877 + t974;
t800 = -mrSges(5,1) * t840 - mrSges(6,1) * t833 + mrSges(7,1) * t830 + mrSges(6,2) * t835 + mrSges(5,3) * t837 - mrSges(7,3) * t832 - pkin(4) * t825 + pkin(5) * t827 - qJ(6) * t965 + t851 * t999 + t978 * t852 + t976 * t897 + t973 * t908 + t971 * t918;
t808 = mrSges(6,1) * t834 + mrSges(7,1) * t828 + mrSges(5,2) * t840 - mrSges(7,2) * t832 - mrSges(5,3) * t836 - mrSges(6,3) * t835 + pkin(5) * t826 - qJ(5) * t825 + t1000 * t852 - t978 * t851 + t977 * t897 + t973 * t907 - t972 * t918;
t893 = Ifges(4,5) * t921 + Ifges(4,6) * t920 + Ifges(4,3) * t937;
t894 = Ifges(4,4) * t921 + Ifges(4,2) * t920 + Ifges(4,6) * t937;
t793 = mrSges(4,2) * t871 - mrSges(4,3) * t844 + Ifges(4,1) * t900 + Ifges(4,4) * t899 + Ifges(4,5) * t924 - pkin(10) * t815 - t948 * t800 + t808 * t994 + t920 * t893 - t937 * t894;
t895 = Ifges(4,1) * t921 + Ifges(4,4) * t920 + Ifges(4,5) * t937;
t794 = -mrSges(4,1) * t871 + mrSges(4,3) * t845 + Ifges(4,4) * t900 + Ifges(4,2) * t899 + Ifges(4,6) * t924 - pkin(3) * t815 - t921 * t893 + t937 * t895 - t996;
t912 = Ifges(3,6) * t942 + (Ifges(3,4) * t950 + Ifges(3,2) * t953) * t981;
t913 = Ifges(3,5) * t942 + (Ifges(3,1) * t950 + Ifges(3,4) * t953) * t981;
t783 = Ifges(3,5) * t931 - Ifges(3,6) * t932 + Ifges(3,3) * t941 + mrSges(3,1) * t901 - mrSges(3,2) * t902 + t949 * t793 + t952 * t794 + pkin(2) * t959 + pkin(9) * t967 + (t912 * t950 - t913 * t953) * t981;
t911 = Ifges(3,3) * t942 + (Ifges(3,5) * t950 + Ifges(3,6) * t953) * t981;
t785 = mrSges(3,2) * t914 - mrSges(3,3) * t901 + Ifges(3,1) * t931 - Ifges(3,4) * t932 + Ifges(3,5) * t941 - pkin(9) * t807 + t952 * t793 - t949 * t794 + t911 * t969 - t942 * t912;
t956 = mrSges(4,1) * t844 - mrSges(4,2) * t845 + Ifges(4,5) * t900 + Ifges(4,6) * t899 + Ifges(4,3) * t924 + pkin(3) * t822 + pkin(10) * t816 + t800 * t994 + t948 * t808 + t921 * t894 - t920 * t895;
t787 = -mrSges(3,1) * t914 + mrSges(3,3) * t902 + Ifges(3,4) * t931 - Ifges(3,2) * t932 + Ifges(3,6) * t941 - pkin(2) * t807 - t911 * t970 + t942 * t913 - t956;
t960 = mrSges(2,1) * t938 - mrSges(2,2) * t939 + Ifges(2,3) * qJDD(1) + pkin(1) * t792 + t947 * t783 + t785 * t988 + t787 * t987 + t799 * t992;
t781 = -mrSges(2,2) * g(3) - mrSges(2,3) * t938 + Ifges(2,5) * qJDD(1) - t955 * Ifges(2,6) + t953 * t785 - t950 * t787 + (-t791 * t946 - t792 * t947) * pkin(8);
t780 = mrSges(2,1) * g(3) + mrSges(2,3) * t939 + t955 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t791 - t946 * t783 + (pkin(8) * t799 + t785 * t950 + t787 * t953) * t947;
t1 = [-m(1) * g(1) + t968; -m(1) * g(2) + t984; (-m(1) - m(2)) * g(3) + t791; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t984 - t951 * t780 + t954 * t781; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t968 + t954 * t780 + t951 * t781; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t960; t960; t783; t956; t996; t823; t827;];
tauJB  = t1;
