% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:31:37
% EndTime: 2019-05-07 07:31:54
% DurationCPUTime: 16.44s
% Computational Cost: add. (253406->364), mult. (569165->448), div. (0->0), fcn. (416384->10), ass. (0->144)
t993 = Ifges(6,1) + Ifges(7,1);
t985 = Ifges(6,4) - Ifges(7,5);
t984 = -Ifges(6,5) - Ifges(7,4);
t992 = Ifges(6,2) + Ifges(7,3);
t983 = Ifges(6,6) - Ifges(7,6);
t991 = -Ifges(6,3) - Ifges(7,2);
t949 = sin(qJ(2));
t952 = cos(qJ(2));
t972 = qJD(1) * qJD(2);
t927 = qJDD(1) * t949 + t952 * t972;
t950 = sin(qJ(1));
t953 = cos(qJ(1));
t934 = -g(1) * t953 - g(2) * t950;
t954 = qJD(1) ^ 2;
t922 = -pkin(1) * t954 + qJDD(1) * pkin(7) + t934;
t981 = t949 * t922;
t987 = pkin(2) * t954;
t885 = qJDD(2) * pkin(2) - t927 * pkin(8) - t981 + (pkin(8) * t972 + t949 * t987 - g(3)) * t952;
t910 = -g(3) * t949 + t952 * t922;
t928 = qJDD(1) * t952 - t949 * t972;
t975 = qJD(1) * t949;
t932 = qJD(2) * pkin(2) - pkin(8) * t975;
t944 = t952 ^ 2;
t886 = pkin(8) * t928 - qJD(2) * t932 - t944 * t987 + t910;
t948 = sin(qJ(3));
t951 = cos(qJ(3));
t856 = t951 * t885 - t948 * t886;
t919 = (-t948 * t949 + t951 * t952) * qJD(1);
t893 = qJD(3) * t919 + t927 * t951 + t928 * t948;
t920 = (t948 * t952 + t949 * t951) * qJD(1);
t941 = qJDD(2) + qJDD(3);
t942 = qJD(2) + qJD(3);
t833 = (t919 * t942 - t893) * qJ(4) + (t919 * t920 + t941) * pkin(3) + t856;
t857 = t948 * t885 + t951 * t886;
t892 = -qJD(3) * t920 - t927 * t948 + t928 * t951;
t912 = pkin(3) * t942 - qJ(4) * t920;
t915 = t919 ^ 2;
t836 = -pkin(3) * t915 + qJ(4) * t892 - t912 * t942 + t857;
t945 = sin(pkin(10));
t946 = cos(pkin(10));
t907 = t919 * t945 + t920 * t946;
t828 = -0.2e1 * qJD(4) * t907 + t946 * t833 - t945 * t836;
t906 = t919 * t946 - t920 * t945;
t829 = 0.2e1 * qJD(4) * t906 + t945 * t833 + t946 * t836;
t865 = t892 * t946 - t893 * t945;
t879 = -mrSges(5,1) * t906 + mrSges(5,2) * t907;
t896 = mrSges(5,1) * t942 - mrSges(5,3) * t907;
t880 = -pkin(4) * t906 - pkin(9) * t907;
t940 = t942 ^ 2;
t826 = -pkin(4) * t940 + pkin(9) * t941 + t880 * t906 + t829;
t933 = t950 * g(1) - t953 * g(2);
t962 = -qJDD(1) * pkin(1) - t933;
t894 = -t928 * pkin(2) + t932 * t975 + (-pkin(8) * t944 - pkin(7)) * t954 + t962;
t842 = -t892 * pkin(3) - t915 * qJ(4) + t920 * t912 + qJDD(4) + t894;
t866 = t892 * t945 + t893 * t946;
t831 = (-t906 * t942 - t866) * pkin(9) + (t907 * t942 - t865) * pkin(4) + t842;
t947 = sin(qJ(5));
t988 = cos(qJ(5));
t823 = t988 * t826 + t947 * t831;
t891 = t907 * t988 + t947 * t942;
t839 = qJD(5) * t891 + t866 * t947 - t941 * t988;
t864 = qJDD(5) - t865;
t900 = qJD(5) - t906;
t872 = mrSges(6,1) * t900 - mrSges(6,3) * t891;
t890 = t907 * t947 - t942 * t988;
t867 = pkin(5) * t890 - qJ(6) * t891;
t899 = t900 ^ 2;
t819 = -pkin(5) * t899 + qJ(6) * t864 + 0.2e1 * qJD(6) * t900 - t867 * t890 + t823;
t873 = -mrSges(7,1) * t900 + mrSges(7,2) * t891;
t971 = m(7) * t819 + t864 * mrSges(7,3) + t900 * t873;
t868 = mrSges(7,1) * t890 - mrSges(7,3) * t891;
t976 = -mrSges(6,1) * t890 - mrSges(6,2) * t891 - t868;
t986 = -mrSges(6,3) - mrSges(7,2);
t810 = m(6) * t823 - t864 * mrSges(6,2) + t839 * t986 - t900 * t872 + t890 * t976 + t971;
t822 = -t947 * t826 + t831 * t988;
t840 = -t890 * qJD(5) + t866 * t988 + t947 * t941;
t871 = -mrSges(6,2) * t900 - mrSges(6,3) * t890;
t820 = -t864 * pkin(5) - t899 * qJ(6) + t891 * t867 + qJDD(6) - t822;
t870 = -mrSges(7,2) * t890 + mrSges(7,3) * t900;
t964 = -m(7) * t820 + t864 * mrSges(7,1) + t900 * t870;
t812 = m(6) * t822 + t864 * mrSges(6,1) + t840 * t986 + t900 * t871 + t891 * t976 + t964;
t966 = t988 * t810 - t812 * t947;
t797 = m(5) * t829 - mrSges(5,2) * t941 + mrSges(5,3) * t865 + t879 * t906 - t896 * t942 + t966;
t895 = -mrSges(5,2) * t942 + mrSges(5,3) * t906;
t825 = -t941 * pkin(4) - t940 * pkin(9) + t907 * t880 - t828;
t821 = -0.2e1 * qJD(6) * t891 + (t890 * t900 - t840) * qJ(6) + (t891 * t900 + t839) * pkin(5) + t825;
t817 = m(7) * t821 + mrSges(7,1) * t839 - t840 * mrSges(7,3) + t870 * t890 - t891 * t873;
t958 = -m(6) * t825 - t839 * mrSges(6,1) - mrSges(6,2) * t840 - t890 * t871 - t872 * t891 - t817;
t807 = m(5) * t828 + mrSges(5,1) * t941 - mrSges(5,3) * t866 - t879 * t907 + t895 * t942 + t958;
t791 = t945 * t797 + t946 * t807;
t908 = -mrSges(4,1) * t919 + mrSges(4,2) * t920;
t911 = -mrSges(4,2) * t942 + mrSges(4,3) * t919;
t788 = m(4) * t856 + mrSges(4,1) * t941 - mrSges(4,3) * t893 - t908 * t920 + t911 * t942 + t791;
t913 = mrSges(4,1) * t942 - mrSges(4,3) * t920;
t967 = t946 * t797 - t807 * t945;
t789 = m(4) * t857 - mrSges(4,2) * t941 + mrSges(4,3) * t892 + t908 * t919 - t913 * t942 + t967;
t782 = t951 * t788 + t948 * t789;
t909 = -t952 * g(3) - t981;
t917 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t949 + Ifges(3,2) * t952) * qJD(1);
t918 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t949 + Ifges(3,4) * t952) * qJD(1);
t977 = t985 * t890 - t891 * t993 + t984 * t900;
t979 = t890 * t983 + t891 * t984 + t900 * t991;
t799 = -mrSges(6,1) * t825 - mrSges(7,1) * t821 + mrSges(7,2) * t819 + mrSges(6,3) * t823 - pkin(5) * t817 - t839 * t992 + t985 * t840 + t983 * t864 + t979 * t891 - t977 * t900;
t978 = t890 * t992 - t891 * t985 - t900 * t983;
t802 = mrSges(6,2) * t825 + mrSges(7,2) * t820 - mrSges(6,3) * t822 - mrSges(7,3) * t821 - qJ(6) * t817 - t985 * t839 + t840 * t993 - t984 * t864 + t979 * t890 + t978 * t900;
t875 = Ifges(5,4) * t907 + Ifges(5,2) * t906 + Ifges(5,6) * t942;
t876 = Ifges(5,1) * t907 + Ifges(5,4) * t906 + Ifges(5,5) * t942;
t902 = Ifges(4,4) * t920 + Ifges(4,2) * t919 + Ifges(4,6) * t942;
t903 = Ifges(4,1) * t920 + Ifges(4,4) * t919 + Ifges(4,5) * t942;
t957 = -mrSges(4,1) * t856 - mrSges(5,1) * t828 + mrSges(4,2) * t857 + mrSges(5,2) * t829 - pkin(3) * t791 - pkin(4) * t958 - pkin(9) * t966 - t988 * t799 - t947 * t802 - t907 * t875 + t906 * t876 + t919 * t903 - Ifges(5,6) * t865 - Ifges(5,5) * t866 - t920 * t902 - Ifges(4,6) * t892 - Ifges(4,5) * t893 + (-Ifges(4,3) - Ifges(5,3)) * t941;
t990 = mrSges(3,1) * t909 - mrSges(3,2) * t910 + Ifges(3,5) * t927 + Ifges(3,6) * t928 + Ifges(3,3) * qJDD(2) + pkin(2) * t782 + (t917 * t949 - t918 * t952) * qJD(1) - t957;
t816 = t840 * mrSges(7,2) + t891 * t868 - t964;
t989 = -t839 * t983 - t840 * t984 - t991 * t864 - t890 * t977 - t978 * t891 + mrSges(6,1) * t822 - mrSges(7,1) * t820 - mrSges(6,2) * t823 + mrSges(7,3) * t819 - pkin(5) * t816 + qJ(6) * (-t839 * mrSges(7,2) - t890 * t868 + t971);
t926 = (-mrSges(3,1) * t952 + mrSges(3,2) * t949) * qJD(1);
t974 = qJD(1) * t952;
t931 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t974;
t780 = m(3) * t909 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t927 + qJD(2) * t931 - t926 * t975 + t782;
t930 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t975;
t968 = -t788 * t948 + t951 * t789;
t781 = m(3) * t910 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t928 - qJD(2) * t930 + t926 * t974 + t968;
t969 = -t780 * t949 + t952 * t781;
t774 = m(2) * t934 - mrSges(2,1) * t954 - qJDD(1) * mrSges(2,2) + t969;
t921 = -t954 * pkin(7) + t962;
t804 = t947 * t810 + t988 * t812;
t800 = m(5) * t842 - t865 * mrSges(5,1) + t866 * mrSges(5,2) - t906 * t895 + t907 * t896 + t804;
t959 = m(4) * t894 - t892 * mrSges(4,1) + mrSges(4,2) * t893 - t919 * t911 + t913 * t920 + t800;
t956 = -m(3) * t921 + t928 * mrSges(3,1) - mrSges(3,2) * t927 - t930 * t975 + t931 * t974 - t959;
t793 = m(2) * t933 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t954 + t956;
t980 = t950 * t774 + t953 * t793;
t776 = t952 * t780 + t949 * t781;
t970 = t953 * t774 - t793 * t950;
t874 = Ifges(5,5) * t907 + Ifges(5,6) * t906 + Ifges(5,3) * t942;
t783 = mrSges(5,2) * t842 - mrSges(5,3) * t828 + Ifges(5,1) * t866 + Ifges(5,4) * t865 + Ifges(5,5) * t941 - pkin(9) * t804 - t947 * t799 + t802 * t988 + t906 * t874 - t942 * t875;
t784 = -mrSges(5,1) * t842 + mrSges(5,3) * t829 + Ifges(5,4) * t866 + Ifges(5,2) * t865 + Ifges(5,6) * t941 - pkin(4) * t804 - t907 * t874 + t942 * t876 - t989;
t901 = Ifges(4,5) * t920 + Ifges(4,6) * t919 + Ifges(4,3) * t942;
t770 = -mrSges(4,1) * t894 + mrSges(4,3) * t857 + Ifges(4,4) * t893 + Ifges(4,2) * t892 + Ifges(4,6) * t941 - pkin(3) * t800 + qJ(4) * t967 + t945 * t783 + t946 * t784 - t920 * t901 + t942 * t903;
t771 = mrSges(4,2) * t894 - mrSges(4,3) * t856 + Ifges(4,1) * t893 + Ifges(4,4) * t892 + Ifges(4,5) * t941 - qJ(4) * t791 + t783 * t946 - t784 * t945 + t901 * t919 - t902 * t942;
t916 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t949 + Ifges(3,6) * t952) * qJD(1);
t766 = -mrSges(3,1) * t921 + mrSges(3,3) * t910 + Ifges(3,4) * t927 + Ifges(3,2) * t928 + Ifges(3,6) * qJDD(2) - pkin(2) * t959 + pkin(8) * t968 + qJD(2) * t918 + t951 * t770 + t948 * t771 - t916 * t975;
t768 = mrSges(3,2) * t921 - mrSges(3,3) * t909 + Ifges(3,1) * t927 + Ifges(3,4) * t928 + Ifges(3,5) * qJDD(2) - pkin(8) * t782 - qJD(2) * t917 - t770 * t948 + t771 * t951 + t916 * t974;
t961 = mrSges(2,1) * t933 - mrSges(2,2) * t934 + Ifges(2,3) * qJDD(1) + pkin(1) * t956 + pkin(7) * t969 + t952 * t766 + t949 * t768;
t769 = mrSges(2,1) * g(3) + mrSges(2,3) * t934 + t954 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t776 - t990;
t764 = -mrSges(2,2) * g(3) - mrSges(2,3) * t933 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t954 - pkin(7) * t776 - t766 * t949 + t768 * t952;
t1 = [-m(1) * g(1) + t970; -m(1) * g(2) + t980; (-m(1) - m(2)) * g(3) + t776; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t980 + t953 * t764 - t950 * t769; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t970 + t950 * t764 + t953 * t769; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t961; t961; t990; -t957; t800; t989; t816;];
tauJB  = t1;
