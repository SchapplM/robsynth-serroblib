% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 12:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:49:06
% EndTime: 2019-05-06 12:49:37
% DurationCPUTime: 31.99s
% Computational Cost: add. (515915->387), mult. (1212060->491), div. (0->0), fcn. (910372->12), ass. (0->151)
t938 = sin(qJ(2));
t941 = cos(qJ(2));
t960 = qJD(1) * qJD(2);
t917 = qJDD(1) * t938 + t941 * t960;
t939 = sin(qJ(1));
t942 = cos(qJ(1));
t924 = -g(1) * t942 - g(2) * t939;
t943 = qJD(1) ^ 2;
t912 = -pkin(1) * t943 + qJDD(1) * pkin(7) + t924;
t964 = t938 * t912;
t966 = pkin(2) * t943;
t873 = qJDD(2) * pkin(2) - t917 * qJ(3) - t964 + (qJ(3) * t960 + t938 * t966 - g(3)) * t941;
t897 = -g(3) * t938 + t912 * t941;
t918 = qJDD(1) * t941 - t938 * t960;
t962 = qJD(1) * t938;
t920 = qJD(2) * pkin(2) - qJ(3) * t962;
t931 = t941 ^ 2;
t874 = qJ(3) * t918 - qJD(2) * t920 - t931 * t966 + t897;
t933 = sin(pkin(10));
t935 = cos(pkin(10));
t907 = (t933 * t941 + t935 * t938) * qJD(1);
t844 = -0.2e1 * qJD(3) * t907 + t873 * t935 - t933 * t874;
t895 = t917 * t935 + t918 * t933;
t906 = (-t933 * t938 + t935 * t941) * qJD(1);
t831 = (qJD(2) * t906 - t895) * pkin(8) + (t906 * t907 + qJDD(2)) * pkin(3) + t844;
t845 = 0.2e1 * qJD(3) * t906 + t873 * t933 + t874 * t935;
t894 = -t917 * t933 + t918 * t935;
t900 = qJD(2) * pkin(3) - pkin(8) * t907;
t905 = t906 ^ 2;
t835 = -pkin(3) * t905 + pkin(8) * t894 - qJD(2) * t900 + t845;
t937 = sin(qJ(4));
t967 = cos(qJ(4));
t822 = t831 * t937 + t835 * t967;
t888 = t906 * t937 + t907 * t967;
t856 = qJD(4) * t888 - t894 * t967 + t895 * t937;
t887 = -t906 * t967 + t907 * t937;
t869 = mrSges(5,1) * t887 + mrSges(5,2) * t888;
t929 = qJD(2) + qJD(4);
t882 = mrSges(5,1) * t929 - mrSges(5,3) * t888;
t928 = qJDD(2) + qJDD(4);
t868 = pkin(4) * t887 - qJ(5) * t888;
t927 = t929 ^ 2;
t816 = -pkin(4) * t927 + qJ(5) * t928 - t868 * t887 + t822;
t923 = t939 * g(1) - g(2) * t942;
t951 = -qJDD(1) * pkin(1) - t923;
t880 = -t918 * pkin(2) + qJDD(3) + t920 * t962 + (-qJ(3) * t931 - pkin(7)) * t943 + t951;
t843 = -t894 * pkin(3) - t905 * pkin(8) + t900 * t907 + t880;
t857 = -qJD(4) * t887 + t894 * t937 + t895 * t967;
t819 = (t887 * t929 - t857) * qJ(5) + (t888 * t929 + t856) * pkin(4) + t843;
t932 = sin(pkin(11));
t934 = cos(pkin(11));
t879 = t888 * t934 + t929 * t932;
t811 = -0.2e1 * qJD(5) * t879 - t932 * t816 + t819 * t934;
t849 = t857 * t934 + t928 * t932;
t878 = -t888 * t932 + t929 * t934;
t809 = (t878 * t887 - t849) * pkin(9) + (t878 * t879 + t856) * pkin(5) + t811;
t812 = 0.2e1 * qJD(5) * t878 + t816 * t934 + t819 * t932;
t848 = -t857 * t932 + t928 * t934;
t862 = pkin(5) * t887 - pkin(9) * t879;
t877 = t878 ^ 2;
t810 = -pkin(5) * t877 + pkin(9) * t848 - t862 * t887 + t812;
t936 = sin(qJ(6));
t940 = cos(qJ(6));
t807 = t809 * t940 - t810 * t936;
t858 = t878 * t940 - t879 * t936;
t825 = qJD(6) * t858 + t848 * t936 + t849 * t940;
t859 = t878 * t936 + t879 * t940;
t832 = -mrSges(7,1) * t858 + mrSges(7,2) * t859;
t883 = qJD(6) + t887;
t836 = -mrSges(7,2) * t883 + mrSges(7,3) * t858;
t855 = qJDD(6) + t856;
t803 = m(7) * t807 + mrSges(7,1) * t855 - mrSges(7,3) * t825 - t832 * t859 + t836 * t883;
t808 = t809 * t936 + t810 * t940;
t824 = -qJD(6) * t859 + t848 * t940 - t849 * t936;
t837 = mrSges(7,1) * t883 - mrSges(7,3) * t859;
t804 = m(7) * t808 - mrSges(7,2) * t855 + mrSges(7,3) * t824 + t832 * t858 - t837 * t883;
t795 = t803 * t940 + t804 * t936;
t860 = -mrSges(6,1) * t878 + mrSges(6,2) * t879;
t953 = -mrSges(6,2) * t887 + mrSges(6,3) * t878;
t793 = m(6) * t811 + t856 * mrSges(6,1) - t849 * mrSges(6,3) - t879 * t860 + t887 * t953 + t795;
t861 = mrSges(6,1) * t887 - mrSges(6,3) * t879;
t954 = -t803 * t936 + t804 * t940;
t794 = m(6) * t812 - mrSges(6,2) * t856 + mrSges(6,3) * t848 + t860 * t878 - t861 * t887 + t954;
t955 = -t793 * t932 + t794 * t934;
t785 = m(5) * t822 - mrSges(5,2) * t928 - mrSges(5,3) * t856 - t869 * t887 - t882 * t929 + t955;
t821 = t831 * t967 - t835 * t937;
t815 = -pkin(4) * t928 - qJ(5) * t927 + t868 * t888 + qJDD(5) - t821;
t813 = -pkin(5) * t848 - pkin(9) * t877 + t862 * t879 + t815;
t948 = m(7) * t813 - mrSges(7,1) * t824 + mrSges(7,2) * t825 - t836 * t858 + t837 * t859;
t806 = m(6) * t815 - mrSges(6,1) * t848 + mrSges(6,2) * t849 + t861 * t879 - t878 * t953 + t948;
t881 = -mrSges(5,2) * t929 - mrSges(5,3) * t887;
t799 = m(5) * t821 + mrSges(5,1) * t928 - mrSges(5,3) * t857 - t869 * t888 + t881 * t929 - t806;
t775 = t785 * t937 + t799 * t967;
t891 = -mrSges(4,1) * t906 + mrSges(4,2) * t907;
t898 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t906;
t773 = m(4) * t844 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t895 + qJD(2) * t898 - t891 * t907 + t775;
t899 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t907;
t956 = t785 * t967 - t799 * t937;
t774 = m(4) * t845 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t894 - qJD(2) * t899 + t891 * t906 + t956;
t768 = t773 * t935 + t774 * t933;
t885 = Ifges(4,4) * t907 + Ifges(4,2) * t906 + Ifges(4,6) * qJD(2);
t886 = Ifges(4,1) * t907 + Ifges(4,4) * t906 + Ifges(4,5) * qJD(2);
t896 = -t941 * g(3) - t964;
t909 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t938 + Ifges(3,2) * t941) * qJD(1);
t910 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t938 + Ifges(3,4) * t941) * qJD(1);
t826 = Ifges(7,5) * t859 + Ifges(7,6) * t858 + Ifges(7,3) * t883;
t828 = Ifges(7,1) * t859 + Ifges(7,4) * t858 + Ifges(7,5) * t883;
t796 = -mrSges(7,1) * t813 + mrSges(7,3) * t808 + Ifges(7,4) * t825 + Ifges(7,2) * t824 + Ifges(7,6) * t855 - t826 * t859 + t828 * t883;
t827 = Ifges(7,4) * t859 + Ifges(7,2) * t858 + Ifges(7,6) * t883;
t797 = mrSges(7,2) * t813 - mrSges(7,3) * t807 + Ifges(7,1) * t825 + Ifges(7,4) * t824 + Ifges(7,5) * t855 + t826 * t858 - t827 * t883;
t838 = Ifges(6,5) * t879 + Ifges(6,6) * t878 + Ifges(6,3) * t887;
t840 = Ifges(6,1) * t879 + Ifges(6,4) * t878 + Ifges(6,5) * t887;
t777 = -mrSges(6,1) * t815 + mrSges(6,3) * t812 + Ifges(6,4) * t849 + Ifges(6,2) * t848 + Ifges(6,6) * t856 - pkin(5) * t948 + pkin(9) * t954 + t940 * t796 + t936 * t797 - t879 * t838 + t887 * t840;
t839 = Ifges(6,4) * t879 + Ifges(6,2) * t878 + Ifges(6,6) * t887;
t779 = mrSges(6,2) * t815 - mrSges(6,3) * t811 + Ifges(6,1) * t849 + Ifges(6,4) * t848 + Ifges(6,5) * t856 - pkin(9) * t795 - t796 * t936 + t797 * t940 + t838 * t878 - t839 * t887;
t864 = Ifges(5,4) * t888 - Ifges(5,2) * t887 + Ifges(5,6) * t929;
t865 = Ifges(5,1) * t888 - Ifges(5,4) * t887 + Ifges(5,5) * t929;
t947 = -mrSges(5,1) * t821 + mrSges(5,2) * t822 - Ifges(5,5) * t857 + Ifges(5,6) * t856 - Ifges(5,3) * t928 + pkin(4) * t806 - qJ(5) * t955 - t777 * t934 - t779 * t932 - t864 * t888 - t887 * t865;
t968 = mrSges(3,1) * t896 + mrSges(4,1) * t844 - mrSges(3,2) * t897 - mrSges(4,2) * t845 + Ifges(3,5) * t917 + Ifges(4,5) * t895 + Ifges(3,6) * t918 + Ifges(4,6) * t894 + pkin(2) * t768 + pkin(3) * t775 + (t909 * t938 - t910 * t941) * qJD(1) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t907 * t885 - t906 * t886 - t947;
t916 = (-mrSges(3,1) * t941 + mrSges(3,2) * t938) * qJD(1);
t961 = qJD(1) * t941;
t922 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t961;
t766 = m(3) * t896 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t917 + qJD(2) * t922 - t916 * t962 + t768;
t921 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t962;
t957 = -t773 * t933 + t774 * t935;
t767 = m(3) * t897 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t918 - qJD(2) * t921 + t916 * t961 + t957;
t958 = -t766 * t938 + t767 * t941;
t759 = m(2) * t924 - mrSges(2,1) * t943 - qJDD(1) * mrSges(2,2) + t958;
t788 = t793 * t934 + t794 * t932;
t950 = m(5) * t843 + mrSges(5,1) * t856 + mrSges(5,2) * t857 + t881 * t887 + t882 * t888 + t788;
t786 = m(4) * t880 - mrSges(4,1) * t894 + mrSges(4,2) * t895 - t898 * t906 + t899 * t907 + t950;
t911 = -t943 * pkin(7) + t951;
t945 = -m(3) * t911 + mrSges(3,1) * t918 - mrSges(3,2) * t917 - t921 * t962 + t922 * t961 - t786;
t781 = m(2) * t923 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t943 + t945;
t963 = t759 * t939 + t781 * t942;
t761 = t766 * t941 + t767 * t938;
t959 = t759 * t942 - t781 * t939;
t863 = Ifges(5,5) * t888 - Ifges(5,6) * t887 + Ifges(5,3) * t929;
t762 = mrSges(5,2) * t843 - mrSges(5,3) * t821 + Ifges(5,1) * t857 - Ifges(5,4) * t856 + Ifges(5,5) * t928 - qJ(5) * t788 - t777 * t932 + t779 * t934 - t863 * t887 - t864 * t929;
t946 = mrSges(7,1) * t807 - mrSges(7,2) * t808 + Ifges(7,5) * t825 + Ifges(7,6) * t824 + Ifges(7,3) * t855 + t859 * t827 - t858 * t828;
t769 = -t946 + (-Ifges(5,2) - Ifges(6,3)) * t856 + Ifges(5,6) * t928 + t929 * t865 - t888 * t863 + t878 * t840 - t879 * t839 + Ifges(5,4) * t857 - Ifges(6,5) * t849 - mrSges(5,1) * t843 - Ifges(6,6) * t848 + mrSges(5,3) * t822 - mrSges(6,1) * t811 + mrSges(6,2) * t812 - pkin(5) * t795 - pkin(4) * t788;
t884 = Ifges(4,5) * t907 + Ifges(4,6) * t906 + Ifges(4,3) * qJD(2);
t755 = -mrSges(4,1) * t880 + mrSges(4,3) * t845 + Ifges(4,4) * t895 + Ifges(4,2) * t894 + Ifges(4,6) * qJDD(2) - pkin(3) * t950 + pkin(8) * t956 + qJD(2) * t886 + t937 * t762 + t769 * t967 - t907 * t884;
t756 = mrSges(4,2) * t880 - mrSges(4,3) * t844 + Ifges(4,1) * t895 + Ifges(4,4) * t894 + Ifges(4,5) * qJDD(2) - pkin(8) * t775 - qJD(2) * t885 + t762 * t967 - t769 * t937 + t884 * t906;
t908 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t938 + Ifges(3,6) * t941) * qJD(1);
t751 = -mrSges(3,1) * t911 + mrSges(3,3) * t897 + Ifges(3,4) * t917 + Ifges(3,2) * t918 + Ifges(3,6) * qJDD(2) - pkin(2) * t786 + qJ(3) * t957 + qJD(2) * t910 + t935 * t755 + t933 * t756 - t908 * t962;
t753 = mrSges(3,2) * t911 - mrSges(3,3) * t896 + Ifges(3,1) * t917 + Ifges(3,4) * t918 + Ifges(3,5) * qJDD(2) - qJ(3) * t768 - qJD(2) * t909 - t755 * t933 + t756 * t935 + t908 * t961;
t949 = mrSges(2,1) * t923 - mrSges(2,2) * t924 + Ifges(2,3) * qJDD(1) + pkin(1) * t945 + pkin(7) * t958 + t751 * t941 + t753 * t938;
t754 = mrSges(2,1) * g(3) + mrSges(2,3) * t924 + t943 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t761 - t968;
t749 = -mrSges(2,2) * g(3) - mrSges(2,3) * t923 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t943 - pkin(7) * t761 - t751 * t938 + t753 * t941;
t1 = [-m(1) * g(1) + t959; -m(1) * g(2) + t963; (-m(1) - m(2)) * g(3) + t761; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t963 + t749 * t942 - t754 * t939; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t959 + t939 * t749 + t942 * t754; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t949; t949; t968; t786; -t947; t806; t946;];
tauJB  = t1;
