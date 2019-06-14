% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 22:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:39:38
% EndTime: 2019-05-06 22:40:23
% DurationCPUTime: 44.73s
% Computational Cost: add. (743147->387), mult. (1637598->484), div. (0->0), fcn. (1220287->12), ass. (0->153)
t937 = sin(qJ(1));
t942 = cos(qJ(1));
t923 = t937 * g(1) - t942 * g(2);
t944 = qJD(1) ^ 2;
t906 = -qJDD(1) * pkin(1) - t944 * pkin(7) - t923;
t936 = sin(qJ(2));
t941 = cos(qJ(2));
t962 = qJD(1) * qJD(2);
t961 = t941 * t962;
t917 = qJDD(1) * t936 + t961;
t927 = t936 * t962;
t918 = qJDD(1) * t941 - t927;
t872 = (-t917 - t961) * qJ(3) + (-t918 + t927) * pkin(2) + t906;
t924 = -g(1) * t942 - g(2) * t937;
t907 = -pkin(1) * t944 + qJDD(1) * pkin(7) + t924;
t889 = -g(3) * t936 + t941 * t907;
t915 = (-pkin(2) * t941 - qJ(3) * t936) * qJD(1);
t943 = qJD(2) ^ 2;
t963 = qJD(1) * t941;
t875 = -pkin(2) * t943 + qJDD(2) * qJ(3) + t915 * t963 + t889;
t931 = sin(pkin(11));
t932 = cos(pkin(11));
t964 = qJD(1) * t936;
t912 = qJD(2) * t931 + t932 * t964;
t852 = -0.2e1 * qJD(3) * t912 + t932 * t872 - t931 * t875;
t893 = qJDD(2) * t931 + t917 * t932;
t911 = qJD(2) * t932 - t931 * t964;
t840 = (-t911 * t963 - t893) * pkin(8) + (t911 * t912 - t918) * pkin(3) + t852;
t853 = 0.2e1 * qJD(3) * t911 + t931 * t872 + t932 * t875;
t892 = qJDD(2) * t932 - t917 * t931;
t894 = -pkin(3) * t963 - pkin(8) * t912;
t910 = t911 ^ 2;
t842 = -pkin(3) * t910 + pkin(8) * t892 + t894 * t963 + t853;
t935 = sin(qJ(4));
t940 = cos(qJ(4));
t823 = t940 * t840 - t935 * t842;
t885 = t911 * t940 - t912 * t935;
t860 = qJD(4) * t885 + t892 * t935 + t893 * t940;
t886 = t911 * t935 + t912 * t940;
t914 = qJDD(4) - t918;
t926 = qJD(4) - t963;
t819 = (t885 * t926 - t860) * pkin(9) + (t885 * t886 + t914) * pkin(4) + t823;
t824 = t935 * t840 + t940 * t842;
t859 = -qJD(4) * t886 + t892 * t940 - t893 * t935;
t878 = pkin(4) * t926 - pkin(9) * t886;
t884 = t885 ^ 2;
t821 = -pkin(4) * t884 + pkin(9) * t859 - t878 * t926 + t824;
t934 = sin(qJ(5));
t939 = cos(qJ(5));
t807 = t939 * t819 - t934 * t821;
t867 = t885 * t939 - t886 * t934;
t836 = qJD(5) * t867 + t859 * t934 + t860 * t939;
t868 = t885 * t934 + t886 * t939;
t908 = qJDD(5) + t914;
t925 = qJD(5) + t926;
t804 = (t867 * t925 - t836) * pkin(10) + (t867 * t868 + t908) * pkin(5) + t807;
t808 = t934 * t819 + t939 * t821;
t835 = -qJD(5) * t868 + t859 * t939 - t860 * t934;
t858 = pkin(5) * t925 - pkin(10) * t868;
t866 = t867 ^ 2;
t805 = -pkin(5) * t866 + pkin(10) * t835 - t858 * t925 + t808;
t933 = sin(qJ(6));
t938 = cos(qJ(6));
t803 = t804 * t933 + t805 * t938;
t888 = -t941 * g(3) - t936 * t907;
t874 = -qJDD(2) * pkin(2) - qJ(3) * t943 + t915 * t964 + qJDD(3) - t888;
t854 = -pkin(3) * t892 - pkin(8) * t910 + t912 * t894 + t874;
t830 = -pkin(4) * t859 - pkin(9) * t884 + t886 * t878 + t854;
t810 = -pkin(5) * t835 - pkin(10) * t866 + t858 * t868 + t830;
t850 = t867 * t933 + t868 * t938;
t815 = -qJD(6) * t850 + t835 * t938 - t836 * t933;
t849 = t867 * t938 - t868 * t933;
t816 = qJD(6) * t849 + t835 * t933 + t836 * t938;
t920 = qJD(6) + t925;
t825 = Ifges(7,5) * t850 + Ifges(7,6) * t849 + Ifges(7,3) * t920;
t827 = Ifges(7,1) * t850 + Ifges(7,4) * t849 + Ifges(7,5) * t920;
t900 = qJDD(6) + t908;
t791 = -mrSges(7,1) * t810 + mrSges(7,3) * t803 + Ifges(7,4) * t816 + Ifges(7,2) * t815 + Ifges(7,6) * t900 - t825 * t850 + t827 * t920;
t802 = t804 * t938 - t805 * t933;
t826 = Ifges(7,4) * t850 + Ifges(7,2) * t849 + Ifges(7,6) * t920;
t792 = mrSges(7,2) * t810 - mrSges(7,3) * t802 + Ifges(7,1) * t816 + Ifges(7,4) * t815 + Ifges(7,5) * t900 + t825 * t849 - t826 * t920;
t845 = Ifges(6,5) * t868 + Ifges(6,6) * t867 + Ifges(6,3) * t925;
t847 = Ifges(6,1) * t868 + Ifges(6,4) * t867 + Ifges(6,5) * t925;
t843 = -mrSges(7,2) * t920 + mrSges(7,3) * t849;
t844 = mrSges(7,1) * t920 - mrSges(7,3) * t850;
t955 = m(7) * t810 - t815 * mrSges(7,1) + t816 * mrSges(7,2) - t849 * t843 + t850 * t844;
t831 = -mrSges(7,1) * t849 + mrSges(7,2) * t850;
t796 = m(7) * t802 + mrSges(7,1) * t900 - mrSges(7,3) * t816 - t831 * t850 + t843 * t920;
t797 = m(7) * t803 - mrSges(7,2) * t900 + mrSges(7,3) * t815 + t831 * t849 - t844 * t920;
t956 = -t796 * t933 + t938 * t797;
t776 = -mrSges(6,1) * t830 + mrSges(6,3) * t808 + Ifges(6,4) * t836 + Ifges(6,2) * t835 + Ifges(6,6) * t908 - pkin(5) * t955 + pkin(10) * t956 + t938 * t791 + t933 * t792 - t868 * t845 + t925 * t847;
t790 = t938 * t796 + t933 * t797;
t846 = Ifges(6,4) * t868 + Ifges(6,2) * t867 + Ifges(6,6) * t925;
t777 = mrSges(6,2) * t830 - mrSges(6,3) * t807 + Ifges(6,1) * t836 + Ifges(6,4) * t835 + Ifges(6,5) * t908 - pkin(10) * t790 - t791 * t933 + t792 * t938 + t845 * t867 - t846 * t925;
t861 = Ifges(5,5) * t886 + Ifges(5,6) * t885 + Ifges(5,3) * t926;
t863 = Ifges(5,1) * t886 + Ifges(5,4) * t885 + Ifges(5,5) * t926;
t855 = -mrSges(6,2) * t925 + mrSges(6,3) * t867;
t856 = mrSges(6,1) * t925 - mrSges(6,3) * t868;
t950 = m(6) * t830 - t835 * mrSges(6,1) + t836 * mrSges(6,2) - t867 * t855 + t868 * t856 + t955;
t851 = -mrSges(6,1) * t867 + mrSges(6,2) * t868;
t787 = m(6) * t807 + mrSges(6,1) * t908 - mrSges(6,3) * t836 - t851 * t868 + t855 * t925 + t790;
t788 = m(6) * t808 - mrSges(6,2) * t908 + mrSges(6,3) * t835 + t851 * t867 - t856 * t925 + t956;
t957 = -t787 * t934 + t939 * t788;
t770 = -mrSges(5,1) * t854 + mrSges(5,3) * t824 + Ifges(5,4) * t860 + Ifges(5,2) * t859 + Ifges(5,6) * t914 - pkin(4) * t950 + pkin(9) * t957 + t939 * t776 + t934 * t777 - t886 * t861 + t926 * t863;
t783 = t939 * t787 + t934 * t788;
t862 = Ifges(5,4) * t886 + Ifges(5,2) * t885 + Ifges(5,6) * t926;
t771 = mrSges(5,2) * t854 - mrSges(5,3) * t823 + Ifges(5,1) * t860 + Ifges(5,4) * t859 + Ifges(5,5) * t914 - pkin(9) * t783 - t776 * t934 + t777 * t939 + t861 * t885 - t862 * t926;
t879 = Ifges(4,5) * t912 + Ifges(4,6) * t911 - Ifges(4,3) * t963;
t881 = Ifges(4,1) * t912 + Ifges(4,4) * t911 - Ifges(4,5) * t963;
t876 = -mrSges(5,2) * t926 + mrSges(5,3) * t885;
t877 = mrSges(5,1) * t926 - mrSges(5,3) * t886;
t946 = m(5) * t854 - t859 * mrSges(5,1) + t860 * mrSges(5,2) - t885 * t876 + t886 * t877 + t950;
t869 = -mrSges(5,1) * t885 + mrSges(5,2) * t886;
t781 = m(5) * t823 + mrSges(5,1) * t914 - mrSges(5,3) * t860 - t869 * t886 + t876 * t926 + t783;
t782 = m(5) * t824 - mrSges(5,2) * t914 + mrSges(5,3) * t859 + t869 * t885 - t877 * t926 + t957;
t958 = -t781 * t935 + t940 * t782;
t753 = -mrSges(4,1) * t874 + mrSges(4,3) * t853 + Ifges(4,4) * t893 + Ifges(4,2) * t892 - Ifges(4,6) * t918 - pkin(3) * t946 + pkin(8) * t958 + t940 * t770 + t935 * t771 - t912 * t879 - t881 * t963;
t775 = t940 * t781 + t935 * t782;
t880 = Ifges(4,4) * t912 + Ifges(4,2) * t911 - Ifges(4,6) * t963;
t754 = mrSges(4,2) * t874 - mrSges(4,3) * t852 + Ifges(4,1) * t893 + Ifges(4,4) * t892 - Ifges(4,5) * t918 - pkin(8) * t775 - t770 * t935 + t771 * t940 + t879 * t911 + t880 * t963;
t887 = -mrSges(4,1) * t911 + mrSges(4,2) * t912;
t953 = mrSges(4,2) * t963 + mrSges(4,3) * t911;
t773 = m(4) * t852 - t918 * mrSges(4,1) - t893 * mrSges(4,3) - t912 * t887 - t953 * t963 + t775;
t891 = -mrSges(4,1) * t963 - mrSges(4,3) * t912;
t774 = m(4) * t853 + mrSges(4,2) * t918 + mrSges(4,3) * t892 + t887 * t911 + t891 * t963 + t958;
t769 = -t773 * t931 + t932 * t774;
t800 = m(4) * t874 - t892 * mrSges(4,1) + t893 * mrSges(4,2) + t912 * t891 - t911 * t953 + t946;
t904 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t936 + Ifges(3,2) * t941) * qJD(1);
t905 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t936 + Ifges(3,4) * t941) * qJD(1);
t966 = mrSges(3,1) * t888 - mrSges(3,2) * t889 + Ifges(3,5) * t917 + Ifges(3,6) * t918 + Ifges(3,3) * qJDD(2) - pkin(2) * t800 + qJ(3) * t769 + t932 * t753 + t931 * t754 + (t904 * t936 - t905 * t941) * qJD(1);
t916 = (-mrSges(3,1) * t941 + mrSges(3,2) * t936) * qJD(1);
t921 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t964;
t767 = m(3) * t889 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t918 - qJD(2) * t921 + t916 * t963 + t769;
t922 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t963;
t799 = m(3) * t888 + qJDD(2) * mrSges(3,1) - t917 * mrSges(3,3) + qJD(2) * t922 - t916 * t964 - t800;
t959 = t941 * t767 - t799 * t936;
t759 = m(2) * t924 - mrSges(2,1) * t944 - qJDD(1) * mrSges(2,2) + t959;
t768 = t773 * t932 + t774 * t931;
t949 = -m(3) * t906 + t918 * mrSges(3,1) - mrSges(3,2) * t917 - t921 * t964 + t922 * t963 - t768;
t763 = m(2) * t923 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t944 + t949;
t965 = t937 * t759 + t942 * t763;
t761 = t936 * t767 + t941 * t799;
t960 = t942 * t759 - t763 * t937;
t903 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t936 + Ifges(3,6) * t941) * qJD(1);
t752 = mrSges(3,2) * t906 - mrSges(3,3) * t888 + Ifges(3,1) * t917 + Ifges(3,4) * t918 + Ifges(3,5) * qJDD(2) - qJ(3) * t768 - qJD(2) * t904 - t753 * t931 + t754 * t932 + t903 * t963;
t951 = -mrSges(7,1) * t802 + mrSges(7,2) * t803 - Ifges(7,5) * t816 - Ifges(7,6) * t815 - Ifges(7,3) * t900 - t850 * t826 + t849 * t827;
t948 = -mrSges(6,1) * t807 + mrSges(6,2) * t808 - Ifges(6,5) * t836 - Ifges(6,6) * t835 - Ifges(6,3) * t908 - pkin(5) * t790 - t868 * t846 + t867 * t847 + t951;
t945 = mrSges(5,1) * t823 - mrSges(5,2) * t824 + Ifges(5,5) * t860 + Ifges(5,6) * t859 + Ifges(5,3) * t914 + pkin(4) * t783 + t886 * t862 - t885 * t863 - t948;
t756 = -t945 - pkin(3) * t775 - t903 * t964 + (Ifges(3,2) + Ifges(4,3)) * t918 + Ifges(3,6) * qJDD(2) - t912 * t880 + Ifges(3,4) * t917 + qJD(2) * t905 - mrSges(3,1) * t906 + t911 * t881 + mrSges(3,3) * t889 - Ifges(4,6) * t892 - Ifges(4,5) * t893 - pkin(2) * t768 - mrSges(4,1) * t852 + mrSges(4,2) * t853;
t952 = mrSges(2,1) * t923 - mrSges(2,2) * t924 + Ifges(2,3) * qJDD(1) + pkin(1) * t949 + pkin(7) * t959 + t936 * t752 + t941 * t756;
t750 = mrSges(2,1) * g(3) + mrSges(2,3) * t924 + t944 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t761 - t966;
t749 = -mrSges(2,2) * g(3) - mrSges(2,3) * t923 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t944 - pkin(7) * t761 + t752 * t941 - t756 * t936;
t1 = [-m(1) * g(1) + t960; -m(1) * g(2) + t965; (-m(1) - m(2)) * g(3) + t761; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t965 + t942 * t749 - t937 * t750; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t960 + t937 * t749 + t942 * t750; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t952; t952; t966; t800; t945; -t948; -t951;];
tauJB  = t1;
