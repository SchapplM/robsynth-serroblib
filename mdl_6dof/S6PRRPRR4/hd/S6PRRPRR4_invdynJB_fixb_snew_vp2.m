% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 05:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:07:04
% EndTime: 2019-05-05 05:07:13
% DurationCPUTime: 7.38s
% Computational Cost: add. (117172->323), mult. (231869->402), div. (0->0), fcn. (149164->12), ass. (0->144)
t977 = Ifges(4,1) + Ifges(5,1);
t969 = Ifges(4,4) - Ifges(5,5);
t968 = Ifges(4,5) + Ifges(5,4);
t976 = Ifges(4,2) + Ifges(5,3);
t967 = Ifges(4,6) - Ifges(5,6);
t975 = Ifges(4,3) + Ifges(5,2);
t926 = sin(qJ(3));
t930 = cos(qJ(3));
t889 = (-mrSges(5,1) * t930 - mrSges(5,3) * t926) * qJD(2);
t953 = qJD(2) * qJD(3);
t951 = t930 * t953;
t891 = t926 * qJDD(2) + t951;
t919 = sin(pkin(11));
t921 = cos(pkin(11));
t895 = t919 * g(1) - t921 * g(2);
t896 = -t921 * g(1) - t919 * g(2);
t917 = -g(3) + qJDD(1);
t931 = cos(qJ(2));
t922 = cos(pkin(6));
t927 = sin(qJ(2));
t962 = t922 * t927;
t920 = sin(pkin(6));
t964 = t920 * t927;
t853 = t895 * t962 + t931 * t896 + t917 * t964;
t933 = qJD(2) ^ 2;
t848 = -t933 * pkin(2) + qJDD(2) * pkin(8) + t853;
t863 = -t920 * t895 + t922 * t917;
t833 = t930 * t848 + t926 * t863;
t888 = (-pkin(3) * t930 - qJ(4) * t926) * qJD(2);
t932 = qJD(3) ^ 2;
t954 = qJD(2) * t930;
t971 = 2 * qJD(4);
t822 = -t932 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t971 + t888 * t954 + t833;
t952 = t926 * t953;
t892 = t930 * qJDD(2) - t952;
t955 = qJD(2) * t926;
t903 = -qJD(3) * pkin(4) - pkin(9) * t955;
t965 = t930 ^ 2 * t933;
t817 = -pkin(4) * t965 - t892 * pkin(9) + qJD(3) * t903 + t822;
t841 = t926 * t848;
t944 = -t932 * qJ(4) + t888 * t955 + qJDD(4) + t841;
t818 = -t891 * pkin(9) + (-pkin(3) - pkin(4)) * qJDD(3) + (-pkin(4) * t926 * t933 + pkin(9) * t953 - t863) * t930 + t944;
t925 = sin(qJ(5));
t929 = cos(qJ(5));
t813 = t929 * t817 + t925 * t818;
t877 = (-t925 * t930 + t926 * t929) * qJD(2);
t843 = -t877 * qJD(5) - t925 * t891 - t929 * t892;
t876 = (t925 * t926 + t929 * t930) * qJD(2);
t856 = t876 * mrSges(6,1) + t877 * mrSges(6,2);
t912 = -qJD(3) + qJD(5);
t862 = t912 * mrSges(6,1) - t877 * mrSges(6,3);
t911 = -qJDD(3) + qJDD(5);
t857 = t876 * pkin(5) - t877 * pkin(10);
t910 = t912 ^ 2;
t810 = -t910 * pkin(5) + t911 * pkin(10) - t876 * t857 + t813;
t961 = t922 * t931;
t963 = t920 * t931;
t852 = t895 * t961 - t927 * t896 + t917 * t963;
t847 = -qJDD(2) * pkin(2) - t933 * pkin(8) - t852;
t940 = -t892 * pkin(3) + t847 + (-t891 - t951) * qJ(4);
t820 = -pkin(3) * t952 + t892 * pkin(4) - pkin(9) * t965 - t940 + (t903 + t971) * t955;
t844 = -t876 * qJD(5) + t929 * t891 - t925 * t892;
t814 = (t877 * t912 - t843) * pkin(5) + (t876 * t912 - t844) * pkin(10) + t820;
t924 = sin(qJ(6));
t928 = cos(qJ(6));
t807 = -t924 * t810 + t928 * t814;
t858 = -t924 * t877 + t928 * t912;
t827 = t858 * qJD(6) + t928 * t844 + t924 * t911;
t859 = t928 * t877 + t924 * t912;
t834 = -t858 * mrSges(7,1) + t859 * mrSges(7,2);
t840 = qJDD(6) - t843;
t867 = qJD(6) + t876;
t845 = -t867 * mrSges(7,2) + t858 * mrSges(7,3);
t803 = m(7) * t807 + t840 * mrSges(7,1) - t827 * mrSges(7,3) - t859 * t834 + t867 * t845;
t808 = t928 * t810 + t924 * t814;
t826 = -t859 * qJD(6) - t924 * t844 + t928 * t911;
t846 = t867 * mrSges(7,1) - t859 * mrSges(7,3);
t804 = m(7) * t808 - t840 * mrSges(7,2) + t826 * mrSges(7,3) + t858 * t834 - t867 * t846;
t947 = -t924 * t803 + t928 * t804;
t791 = m(6) * t813 - t911 * mrSges(6,2) + t843 * mrSges(6,3) - t876 * t856 - t912 * t862 + t947;
t812 = -t925 * t817 + t929 * t818;
t861 = -t912 * mrSges(6,2) - t876 * mrSges(6,3);
t809 = -t911 * pkin(5) - t910 * pkin(10) + t877 * t857 - t812;
t939 = -m(7) * t809 + t826 * mrSges(7,1) - t827 * mrSges(7,2) + t858 * t845 - t859 * t846;
t799 = m(6) * t812 + t911 * mrSges(6,1) - t844 * mrSges(6,3) - t877 * t856 + t912 * t861 + t939;
t785 = t925 * t791 + t929 * t799;
t960 = t930 * t863;
t823 = -qJDD(3) * pkin(3) + t944 - t960;
t900 = mrSges(5,2) * t954 + qJD(3) * mrSges(5,3);
t938 = -m(5) * t823 + qJDD(3) * mrSges(5,1) + qJD(3) * t900 - t785;
t784 = t891 * mrSges(5,2) + t889 * t955 - t938;
t832 = -t841 + t960;
t828 = Ifges(7,5) * t859 + Ifges(7,6) * t858 + Ifges(7,3) * t867;
t830 = Ifges(7,1) * t859 + Ifges(7,4) * t858 + Ifges(7,5) * t867;
t797 = -mrSges(7,1) * t809 + mrSges(7,3) * t808 + Ifges(7,4) * t827 + Ifges(7,2) * t826 + Ifges(7,6) * t840 - t859 * t828 + t867 * t830;
t829 = Ifges(7,4) * t859 + Ifges(7,2) * t858 + Ifges(7,6) * t867;
t798 = mrSges(7,2) * t809 - mrSges(7,3) * t807 + Ifges(7,1) * t827 + Ifges(7,4) * t826 + Ifges(7,5) * t840 + t858 * t828 - t867 * t829;
t850 = Ifges(6,4) * t877 - Ifges(6,2) * t876 + Ifges(6,6) * t912;
t851 = Ifges(6,1) * t877 - Ifges(6,4) * t876 + Ifges(6,5) * t912;
t937 = -mrSges(6,1) * t812 + mrSges(6,2) * t813 - Ifges(6,5) * t844 - Ifges(6,6) * t843 - Ifges(6,3) * t911 - pkin(5) * t939 - pkin(10) * t947 - t928 * t797 - t924 * t798 - t877 * t850 - t876 * t851;
t898 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t955;
t948 = t929 * t791 - t925 * t799;
t941 = m(5) * t822 + qJDD(3) * mrSges(5,3) + qJD(3) * t898 + t889 * t954 + t948;
t956 = t968 * qJD(3) + (t977 * t926 + t969 * t930) * qJD(2);
t957 = -t967 * qJD(3) + (-t969 * t926 - t976 * t930) * qJD(2);
t974 = -(t957 * t926 + t956 * t930) * qJD(2) + t975 * qJDD(3) + t968 * t891 + t967 * t892 + mrSges(4,1) * t832 - mrSges(5,1) * t823 - mrSges(4,2) * t833 + mrSges(5,3) * t822 - pkin(3) * t784 - pkin(4) * t785 + qJ(4) * (t892 * mrSges(5,2) + t941) + t937;
t970 = mrSges(4,3) + mrSges(5,2);
t890 = (-mrSges(4,1) * t930 + mrSges(4,2) * t926) * qJD(2);
t897 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t955;
t781 = m(4) * t833 - qJDD(3) * mrSges(4,2) - qJD(3) * t897 + t890 * t954 + t970 * t892 + t941;
t899 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t954;
t782 = m(4) * t832 + qJDD(3) * mrSges(4,1) + qJD(3) * t899 - t970 * t891 + (-t889 - t890) * t955 + t938;
t949 = t930 * t781 - t926 * t782;
t772 = m(3) * t853 - t933 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t949;
t775 = t926 * t781 + t930 * t782;
t774 = m(3) * t863 + t775;
t824 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t955 + t940;
t793 = t928 * t803 + t924 * t804;
t943 = -m(6) * t820 + t843 * mrSges(6,1) - t844 * mrSges(6,2) - t876 * t861 - t877 * t862 - t793;
t789 = m(5) * t824 - t892 * mrSges(5,1) - t891 * mrSges(5,3) - t898 * t955 - t900 * t954 + t943;
t935 = -m(4) * t847 + t892 * mrSges(4,1) - t891 * mrSges(4,2) - t897 * t955 + t899 * t954 - t789;
t788 = m(3) * t852 + qJDD(2) * mrSges(3,1) - t933 * mrSges(3,2) + t935;
t763 = t772 * t962 - t920 * t774 + t788 * t961;
t761 = m(2) * t895 + t763;
t768 = t931 * t772 - t927 * t788;
t767 = m(2) * t896 + t768;
t959 = t921 * t761 + t919 * t767;
t958 = t975 * qJD(3) + (t968 * t926 + t967 * t930) * qJD(2);
t762 = t772 * t964 + t922 * t774 + t788 * t963;
t950 = -t919 * t761 + t921 * t767;
t946 = m(2) * t917 + t762;
t849 = Ifges(6,5) * t877 - Ifges(6,6) * t876 + Ifges(6,3) * t912;
t776 = mrSges(6,2) * t820 - mrSges(6,3) * t812 + Ifges(6,1) * t844 + Ifges(6,4) * t843 + Ifges(6,5) * t911 - pkin(10) * t793 - t924 * t797 + t928 * t798 - t876 * t849 - t912 * t850;
t936 = mrSges(7,1) * t807 - mrSges(7,2) * t808 + Ifges(7,5) * t827 + Ifges(7,6) * t826 + Ifges(7,3) * t840 + t859 * t829 - t858 * t830;
t777 = -mrSges(6,1) * t820 + mrSges(6,3) * t813 + Ifges(6,4) * t844 + Ifges(6,2) * t843 + Ifges(6,6) * t911 - pkin(5) * t793 - t877 * t849 + t912 * t851 - t936;
t759 = -mrSges(4,1) * t847 - mrSges(5,1) * t824 + mrSges(5,2) * t822 + mrSges(4,3) * t833 - pkin(3) * t789 - pkin(4) * t943 - pkin(9) * t948 + t956 * qJD(3) + t967 * qJDD(3) - t925 * t776 - t929 * t777 + t969 * t891 + t976 * t892 - t958 * t955;
t764 = mrSges(4,2) * t847 + mrSges(5,2) * t823 - mrSges(4,3) * t832 - mrSges(5,3) * t824 - pkin(9) * t785 - qJ(4) * t789 + t957 * qJD(3) + t968 * qJDD(3) + t929 * t776 - t925 * t777 + t977 * t891 + t969 * t892 + t958 * t954;
t757 = mrSges(3,2) * t863 - mrSges(3,3) * t852 + Ifges(3,5) * qJDD(2) - t933 * Ifges(3,6) - pkin(8) * t775 - t926 * t759 + t930 * t764;
t758 = -mrSges(3,1) * t863 + mrSges(3,3) * t853 + t933 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t775 - t974;
t942 = pkin(7) * t768 + t757 * t927 + t758 * t931;
t756 = mrSges(3,1) * t852 - mrSges(3,2) * t853 + Ifges(3,3) * qJDD(2) + pkin(2) * t935 + pkin(8) * t949 + t930 * t759 + t926 * t764;
t755 = mrSges(2,2) * t917 - mrSges(2,3) * t895 + t931 * t757 - t927 * t758 + (-t762 * t920 - t763 * t922) * pkin(7);
t754 = -mrSges(2,1) * t917 + mrSges(2,3) * t896 - pkin(1) * t762 - t920 * t756 + t922 * t942;
t1 = [-m(1) * g(1) + t950; -m(1) * g(2) + t959; -m(1) * g(3) + t946; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t959 - t919 * t754 + t921 * t755; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t950 + t921 * t754 + t919 * t755; -mrSges(1,1) * g(2) + mrSges(2,1) * t895 + mrSges(1,2) * g(1) - mrSges(2,2) * t896 + pkin(1) * t763 + t922 * t756 + t920 * t942; t946; t756; t974; t784; -t937; t936;];
tauJB  = t1;
