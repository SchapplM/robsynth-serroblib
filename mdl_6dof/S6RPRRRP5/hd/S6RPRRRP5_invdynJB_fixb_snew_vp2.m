% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:29:15
% EndTime: 2019-05-06 01:29:30
% DurationCPUTime: 14.85s
% Computational Cost: add. (219127->342), mult. (525738->414), div. (0->0), fcn. (409379->10), ass. (0->146)
t979 = Ifges(6,1) + Ifges(7,1);
t972 = Ifges(6,4) - Ifges(7,5);
t971 = -Ifges(6,5) - Ifges(7,4);
t978 = Ifges(6,2) + Ifges(7,3);
t970 = Ifges(6,6) - Ifges(7,6);
t977 = -Ifges(6,3) - Ifges(7,2);
t934 = qJD(1) ^ 2;
t925 = sin(pkin(10));
t926 = cos(pkin(10));
t929 = sin(qJ(3));
t932 = cos(qJ(3));
t944 = -t925 * t929 + t926 * t932;
t902 = t944 * qJD(1);
t945 = t925 * t932 + t926 * t929;
t903 = t945 * qJD(1);
t928 = sin(qJ(4));
t931 = cos(qJ(4));
t885 = t902 * t931 - t903 * t928;
t892 = -t903 * qJD(3) + qJDD(1) * t944;
t960 = t902 * qJD(3);
t893 = qJDD(1) * t945 + t960;
t854 = qJD(4) * t885 + t892 * t928 + t893 * t931;
t886 = t902 * t928 + t903 * t931;
t923 = qJD(3) + qJD(4);
t927 = sin(qJ(5));
t975 = cos(qJ(5));
t871 = t886 * t927 - t923 * t975;
t920 = qJDD(3) + qJDD(4);
t828 = -qJD(5) * t871 + t854 * t975 + t920 * t927;
t872 = t886 * t975 + t923 * t927;
t856 = mrSges(7,1) * t871 - mrSges(7,3) * t872;
t930 = sin(qJ(1));
t933 = cos(qJ(1));
t909 = -g(1) * t933 - g(2) * t930;
t904 = -pkin(1) * t934 + qJDD(1) * qJ(2) + t909;
t959 = qJD(1) * qJD(2);
t956 = -t926 * g(3) - 0.2e1 * t925 * t959;
t974 = pkin(2) * t926;
t878 = (-pkin(7) * qJDD(1) + t934 * t974 - t904) * t925 + t956;
t895 = -g(3) * t925 + (t904 + 0.2e1 * t959) * t926;
t958 = qJDD(1) * t926;
t922 = t926 ^ 2;
t967 = t922 * t934;
t879 = -pkin(2) * t967 + pkin(7) * t958 + t895;
t858 = t878 * t932 - t929 * t879;
t824 = (-t893 + t960) * pkin(8) + (t902 * t903 + qJDD(3)) * pkin(3) + t858;
t859 = t878 * t929 + t879 * t932;
t898 = qJD(3) * pkin(3) - pkin(8) * t903;
t901 = t902 ^ 2;
t831 = -pkin(3) * t901 + pkin(8) * t892 - qJD(3) * t898 + t859;
t822 = t824 * t928 + t831 * t931;
t870 = -pkin(4) * t885 - pkin(9) * t886;
t919 = t923 ^ 2;
t817 = -pkin(4) * t919 + pkin(9) * t920 + t870 * t885 + t822;
t921 = t925 ^ 2;
t908 = t930 * g(1) - g(2) * t933;
t949 = qJDD(2) - t908;
t891 = (-pkin(1) - t974) * qJDD(1) + (-qJ(2) + (-t921 - t922) * pkin(7)) * t934 + t949;
t844 = -pkin(3) * t892 - pkin(8) * t901 + t898 * t903 + t891;
t853 = -qJD(4) * t886 + t892 * t931 - t893 * t928;
t819 = (-t885 * t923 - t854) * pkin(9) + (t886 * t923 - t853) * pkin(4) + t844;
t813 = -t817 * t927 + t819 * t975;
t852 = qJDD(5) - t853;
t855 = pkin(5) * t871 - qJ(6) * t872;
t881 = qJD(5) - t885;
t880 = t881 ^ 2;
t811 = -pkin(5) * t852 - qJ(6) * t880 + t855 * t872 + qJDD(6) - t813;
t860 = -mrSges(7,2) * t871 + mrSges(7,3) * t881;
t950 = -m(7) * t811 + mrSges(7,1) * t852 + t860 * t881;
t807 = mrSges(7,2) * t828 + t856 * t872 - t950;
t814 = t817 * t975 + t819 * t927;
t810 = -pkin(5) * t880 + qJ(6) * t852 + 0.2e1 * qJD(6) * t881 - t855 * t871 + t814;
t827 = qJD(5) * t872 + t854 * t927 - t920 * t975;
t863 = -mrSges(7,1) * t881 + mrSges(7,2) * t872;
t957 = m(7) * t810 + mrSges(7,3) * t852 + t863 * t881;
t963 = t871 * t972 - t872 * t979 + t881 * t971;
t964 = t871 * t978 - t872 * t972 - t881 * t970;
t976 = -t970 * t827 - t971 * t828 - t977 * t852 - t963 * t871 - t964 * t872 + mrSges(6,1) * t813 - mrSges(7,1) * t811 - mrSges(6,2) * t814 + mrSges(7,3) * t810 - pkin(5) * t807 + qJ(6) * (-t827 * mrSges(7,2) - t856 * t871 + t957);
t973 = -mrSges(6,3) - mrSges(7,2);
t968 = mrSges(3,2) * t925;
t869 = -mrSges(5,1) * t885 + mrSges(5,2) * t886;
t876 = mrSges(5,1) * t923 - mrSges(5,3) * t886;
t862 = mrSges(6,1) * t881 - mrSges(6,3) * t872;
t962 = -mrSges(6,1) * t871 - mrSges(6,2) * t872 - t856;
t801 = m(6) * t814 - mrSges(6,2) * t852 + t827 * t973 - t862 * t881 + t871 * t962 + t957;
t861 = -mrSges(6,2) * t881 - mrSges(6,3) * t871;
t803 = m(6) * t813 + mrSges(6,1) * t852 + t828 * t973 + t861 * t881 + t872 * t962 + t950;
t951 = t801 * t975 - t803 * t927;
t789 = m(5) * t822 - mrSges(5,2) * t920 + mrSges(5,3) * t853 + t869 * t885 - t876 * t923 + t951;
t821 = t824 * t931 - t831 * t928;
t875 = -mrSges(5,2) * t923 + mrSges(5,3) * t885;
t816 = -pkin(4) * t920 - pkin(9) * t919 + t870 * t886 - t821;
t812 = -0.2e1 * qJD(6) * t872 + (t871 * t881 - t828) * qJ(6) + (t872 * t881 + t827) * pkin(5) + t816;
t808 = m(7) * t812 + mrSges(7,1) * t827 - mrSges(7,3) * t828 + t860 * t871 - t863 * t872;
t937 = -m(6) * t816 - mrSges(6,1) * t827 - mrSges(6,2) * t828 - t861 * t871 - t862 * t872 - t808;
t798 = m(5) * t821 + mrSges(5,1) * t920 - mrSges(5,3) * t854 - t869 * t886 + t875 * t923 + t937;
t782 = t789 * t928 + t798 * t931;
t889 = -mrSges(4,1) * t902 + mrSges(4,2) * t903;
t896 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t902;
t780 = m(4) * t858 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t893 + qJD(3) * t896 - t889 * t903 + t782;
t897 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t903;
t952 = t789 * t931 - t798 * t928;
t781 = m(4) * t859 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t892 - qJD(3) * t897 + t889 * t902 + t952;
t774 = t780 * t932 + t781 * t929;
t894 = -t925 * t904 + t956;
t943 = mrSges(3,3) * qJDD(1) + t934 * (-mrSges(3,1) * t926 + t968);
t772 = m(3) * t894 - t925 * t943 + t774;
t953 = -t929 * t780 + t781 * t932;
t773 = m(3) * t895 + t926 * t943 + t953;
t954 = -t772 * t925 + t773 * t926;
t766 = m(2) * t909 - mrSges(2,1) * t934 - qJDD(1) * mrSges(2,2) + t954;
t900 = -qJDD(1) * pkin(1) - t934 * qJ(2) + t949;
t795 = t801 * t927 + t803 * t975;
t942 = m(5) * t844 - mrSges(5,1) * t853 + mrSges(5,2) * t854 - t875 * t885 + t876 * t886 + t795;
t938 = m(4) * t891 - mrSges(4,1) * t892 + mrSges(4,2) * t893 - t896 * t902 + t897 * t903 + t942;
t936 = -m(3) * t900 + mrSges(3,1) * t958 - t938 + (t921 * t934 + t967) * mrSges(3,3);
t784 = t936 + (mrSges(2,1) - t968) * qJDD(1) + m(2) * t908 - mrSges(2,2) * t934;
t966 = t766 * t930 + t784 * t933;
t768 = t772 * t926 + t773 * t925;
t965 = t871 * t970 + t872 * t971 + t881 * t977;
t946 = Ifges(3,5) * t925 + Ifges(3,6) * t926;
t961 = t934 * t946;
t955 = t766 * t933 - t784 * t930;
t948 = Ifges(3,1) * t925 + Ifges(3,4) * t926;
t947 = Ifges(3,4) * t925 + Ifges(3,2) * t926;
t791 = -mrSges(6,1) * t816 - mrSges(7,1) * t812 + mrSges(7,2) * t810 + mrSges(6,3) * t814 - pkin(5) * t808 - t827 * t978 + t828 * t972 + t852 * t970 + t872 * t965 - t881 * t963;
t793 = mrSges(6,2) * t816 + mrSges(7,2) * t811 - mrSges(6,3) * t813 - mrSges(7,3) * t812 - qJ(6) * t808 - t827 * t972 + t828 * t979 - t852 * t971 + t871 * t965 + t881 * t964;
t864 = Ifges(5,5) * t886 + Ifges(5,6) * t885 + Ifges(5,3) * t923;
t865 = Ifges(5,4) * t886 + Ifges(5,2) * t885 + Ifges(5,6) * t923;
t775 = mrSges(5,2) * t844 - mrSges(5,3) * t821 + Ifges(5,1) * t854 + Ifges(5,4) * t853 + Ifges(5,5) * t920 - pkin(9) * t795 - t791 * t927 + t793 * t975 + t864 * t885 - t865 * t923;
t866 = Ifges(5,1) * t886 + Ifges(5,4) * t885 + Ifges(5,5) * t923;
t776 = -mrSges(5,1) * t844 + mrSges(5,3) * t822 + Ifges(5,4) * t854 + Ifges(5,2) * t853 + Ifges(5,6) * t920 - pkin(4) * t795 - t886 * t864 + t923 * t866 - t976;
t882 = Ifges(4,5) * t903 + Ifges(4,6) * t902 + Ifges(4,3) * qJD(3);
t884 = Ifges(4,1) * t903 + Ifges(4,4) * t902 + Ifges(4,5) * qJD(3);
t762 = -mrSges(4,1) * t891 + mrSges(4,3) * t859 + Ifges(4,4) * t893 + Ifges(4,2) * t892 + Ifges(4,6) * qJDD(3) - pkin(3) * t942 + pkin(8) * t952 + qJD(3) * t884 + t928 * t775 + t931 * t776 - t903 * t882;
t883 = Ifges(4,4) * t903 + Ifges(4,2) * t902 + Ifges(4,6) * qJD(3);
t763 = mrSges(4,2) * t891 - mrSges(4,3) * t858 + Ifges(4,1) * t893 + Ifges(4,4) * t892 + Ifges(4,5) * qJDD(3) - pkin(8) * t782 - qJD(3) * t883 + t775 * t931 - t776 * t928 + t882 * t902;
t758 = -mrSges(3,1) * t900 + mrSges(3,3) * t895 - pkin(2) * t938 + pkin(7) * t953 + qJDD(1) * t947 + t932 * t762 + t929 * t763 - t925 * t961;
t760 = mrSges(3,2) * t900 - mrSges(3,3) * t894 - pkin(7) * t774 + qJDD(1) * t948 - t762 * t929 + t763 * t932 + t926 * t961;
t786 = qJDD(1) * t968 - t936;
t941 = mrSges(2,1) * t908 - mrSges(2,2) * t909 + Ifges(2,3) * qJDD(1) - pkin(1) * t786 + qJ(2) * t954 + t758 * t926 + t760 * t925;
t940 = -mrSges(5,1) * t821 + mrSges(5,2) * t822 - Ifges(5,5) * t854 - Ifges(5,6) * t853 - Ifges(5,3) * t920 - pkin(4) * t937 - pkin(9) * t951 - t791 * t975 - t793 * t927 - t865 * t886 + t885 * t866;
t935 = mrSges(4,1) * t858 - mrSges(4,2) * t859 + Ifges(4,5) * t893 + Ifges(4,6) * t892 + Ifges(4,3) * qJDD(3) + pkin(3) * t782 + t903 * t883 - t902 * t884 - t940;
t761 = mrSges(2,3) * t909 - mrSges(3,1) * t894 + mrSges(3,2) * t895 + mrSges(2,1) * g(3) - t935 - pkin(2) * t774 - pkin(1) * t768 + (Ifges(2,6) - t946) * qJDD(1) + (-t925 * t947 + t926 * t948 + Ifges(2,5)) * t934;
t756 = -mrSges(2,2) * g(3) - mrSges(2,3) * t908 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t934 - qJ(2) * t768 - t758 * t925 + t760 * t926;
t1 = [-m(1) * g(1) + t955; -m(1) * g(2) + t966; (-m(1) - m(2)) * g(3) + t768; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t966 + t756 * t933 - t761 * t930; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t955 + t930 * t756 + t933 * t761; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t941; t941; t786; t935; -t940; t976; t807;];
tauJB  = t1;
