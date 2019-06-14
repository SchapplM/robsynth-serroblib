% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:23:44
% EndTime: 2019-05-06 17:24:01
% DurationCPUTime: 16.24s
% Computational Cost: add. (239141->364), mult. (553607->448), div. (0->0), fcn. (406963->10), ass. (0->144)
t971 = Ifges(6,1) + Ifges(7,1);
t963 = Ifges(6,4) - Ifges(7,5);
t962 = -Ifges(6,5) - Ifges(7,4);
t970 = Ifges(6,2) + Ifges(7,3);
t961 = Ifges(6,6) - Ifges(7,6);
t969 = -Ifges(6,3) - Ifges(7,2);
t928 = sin(qJ(2));
t931 = cos(qJ(2));
t950 = qJD(1) * qJD(2);
t907 = qJDD(1) * t928 + t931 * t950;
t929 = sin(qJ(1));
t932 = cos(qJ(1));
t914 = -g(1) * t932 - g(2) * t929;
t933 = qJD(1) ^ 2;
t902 = -pkin(1) * t933 + qJDD(1) * pkin(7) + t914;
t958 = t928 * t902;
t965 = pkin(2) * t933;
t867 = qJDD(2) * pkin(2) - t907 * qJ(3) - t958 + (qJ(3) * t950 + t928 * t965 - g(3)) * t931;
t888 = -g(3) * t928 + t931 * t902;
t908 = qJDD(1) * t931 - t928 * t950;
t952 = qJD(1) * t928;
t910 = qJD(2) * pkin(2) - qJ(3) * t952;
t923 = t931 ^ 2;
t868 = qJ(3) * t908 - qJD(2) * t910 - t923 * t965 + t888;
t924 = sin(pkin(10));
t925 = cos(pkin(10));
t897 = (t924 * t931 + t925 * t928) * qJD(1);
t837 = -0.2e1 * qJD(3) * t897 + t925 * t867 - t924 * t868;
t886 = t907 * t925 + t908 * t924;
t896 = (-t924 * t928 + t925 * t931) * qJD(1);
t817 = (qJD(2) * t896 - t886) * pkin(8) + (t896 * t897 + qJDD(2)) * pkin(3) + t837;
t838 = 0.2e1 * qJD(3) * t896 + t924 * t867 + t925 * t868;
t885 = -t907 * t924 + t908 * t925;
t891 = qJD(2) * pkin(3) - pkin(8) * t897;
t895 = t896 ^ 2;
t820 = -pkin(3) * t895 + pkin(8) * t885 - qJD(2) * t891 + t838;
t927 = sin(qJ(4));
t930 = cos(qJ(4));
t815 = t927 * t817 + t930 * t820;
t880 = t896 * t927 + t897 * t930;
t848 = -qJD(4) * t880 + t885 * t930 - t886 * t927;
t879 = t896 * t930 - t897 * t927;
t862 = -mrSges(5,1) * t879 + mrSges(5,2) * t880;
t921 = qJD(2) + qJD(4);
t873 = mrSges(5,1) * t921 - mrSges(5,3) * t880;
t920 = qJDD(2) + qJDD(4);
t863 = -pkin(4) * t879 - pkin(9) * t880;
t919 = t921 ^ 2;
t810 = -pkin(4) * t919 + pkin(9) * t920 + t863 * t879 + t815;
t913 = t929 * g(1) - t932 * g(2);
t941 = -qJDD(1) * pkin(1) - t913;
t869 = -t908 * pkin(2) + qJDD(3) + t910 * t952 + (-qJ(3) * t923 - pkin(7)) * t933 + t941;
t834 = -t885 * pkin(3) - t895 * pkin(8) + t897 * t891 + t869;
t849 = qJD(4) * t879 + t885 * t927 + t886 * t930;
t812 = (-t879 * t921 - t849) * pkin(9) + (t880 * t921 - t848) * pkin(4) + t834;
t926 = sin(qJ(5));
t966 = cos(qJ(5));
t807 = t966 * t810 + t926 * t812;
t871 = t966 * t880 + t926 * t921;
t823 = qJD(5) * t871 + t849 * t926 - t966 * t920;
t847 = qJDD(5) - t848;
t875 = qJD(5) - t879;
t855 = mrSges(6,1) * t875 - mrSges(6,3) * t871;
t870 = t880 * t926 - t966 * t921;
t850 = pkin(5) * t870 - qJ(6) * t871;
t874 = t875 ^ 2;
t803 = -pkin(5) * t874 + qJ(6) * t847 + 0.2e1 * qJD(6) * t875 - t850 * t870 + t807;
t856 = -mrSges(7,1) * t875 + mrSges(7,2) * t871;
t949 = m(7) * t803 + t847 * mrSges(7,3) + t875 * t856;
t851 = mrSges(7,1) * t870 - mrSges(7,3) * t871;
t953 = -mrSges(6,1) * t870 - mrSges(6,2) * t871 - t851;
t964 = -mrSges(6,3) - mrSges(7,2);
t794 = m(6) * t807 - t847 * mrSges(6,2) + t964 * t823 - t875 * t855 + t953 * t870 + t949;
t806 = -t926 * t810 + t966 * t812;
t824 = -t870 * qJD(5) + t966 * t849 + t926 * t920;
t854 = -mrSges(6,2) * t875 - mrSges(6,3) * t870;
t804 = -t847 * pkin(5) - t874 * qJ(6) + t871 * t850 + qJDD(6) - t806;
t853 = -mrSges(7,2) * t870 + mrSges(7,3) * t875;
t943 = -m(7) * t804 + t847 * mrSges(7,1) + t875 * t853;
t796 = m(6) * t806 + t847 * mrSges(6,1) + t964 * t824 + t875 * t854 + t953 * t871 + t943;
t944 = t966 * t794 - t796 * t926;
t781 = m(5) * t815 - mrSges(5,2) * t920 + mrSges(5,3) * t848 + t862 * t879 - t873 * t921 + t944;
t814 = t930 * t817 - t927 * t820;
t872 = -mrSges(5,2) * t921 + mrSges(5,3) * t879;
t809 = -t920 * pkin(4) - t919 * pkin(9) + t880 * t863 - t814;
t805 = -0.2e1 * qJD(6) * t871 + (t870 * t875 - t824) * qJ(6) + (t871 * t875 + t823) * pkin(5) + t809;
t801 = m(7) * t805 + mrSges(7,1) * t823 - t824 * mrSges(7,3) + t853 * t870 - t871 * t856;
t936 = -m(6) * t809 - t823 * mrSges(6,1) - mrSges(6,2) * t824 - t870 * t854 - t855 * t871 - t801;
t791 = m(5) * t814 + mrSges(5,1) * t920 - mrSges(5,3) * t849 - t862 * t880 + t872 * t921 + t936;
t775 = t927 * t781 + t930 * t791;
t883 = -mrSges(4,1) * t896 + mrSges(4,2) * t897;
t889 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t896;
t773 = m(4) * t837 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t886 + qJD(2) * t889 - t883 * t897 + t775;
t890 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t897;
t945 = t930 * t781 - t791 * t927;
t774 = m(4) * t838 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t885 - qJD(2) * t890 + t883 * t896 + t945;
t767 = t925 * t773 + t924 * t774;
t877 = Ifges(4,4) * t897 + Ifges(4,2) * t896 + Ifges(4,6) * qJD(2);
t878 = Ifges(4,1) * t897 + Ifges(4,4) * t896 + Ifges(4,5) * qJD(2);
t887 = -t931 * g(3) - t958;
t899 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t928 + Ifges(3,2) * t931) * qJD(1);
t900 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t928 + Ifges(3,4) * t931) * qJD(1);
t954 = t963 * t870 - t871 * t971 + t962 * t875;
t956 = t870 * t961 + t871 * t962 + t875 * t969;
t784 = -mrSges(6,1) * t809 - mrSges(7,1) * t805 + mrSges(7,2) * t803 + mrSges(6,3) * t807 - pkin(5) * t801 - t823 * t970 + t963 * t824 + t961 * t847 + t956 * t871 - t954 * t875;
t955 = t870 * t970 - t871 * t963 - t875 * t961;
t786 = mrSges(6,2) * t809 + mrSges(7,2) * t804 - mrSges(6,3) * t806 - mrSges(7,3) * t805 - qJ(6) * t801 - t963 * t823 + t824 * t971 - t962 * t847 + t956 * t870 + t955 * t875;
t858 = Ifges(5,4) * t880 + Ifges(5,2) * t879 + Ifges(5,6) * t921;
t859 = Ifges(5,1) * t880 + Ifges(5,4) * t879 + Ifges(5,5) * t921;
t938 = -mrSges(5,1) * t814 + mrSges(5,2) * t815 - Ifges(5,5) * t849 - Ifges(5,6) * t848 - Ifges(5,3) * t920 - pkin(4) * t936 - pkin(9) * t944 - t966 * t784 - t926 * t786 - t880 * t858 + t879 * t859;
t968 = mrSges(3,1) * t887 + mrSges(4,1) * t837 - mrSges(3,2) * t888 - mrSges(4,2) * t838 + Ifges(3,5) * t907 + Ifges(4,5) * t886 + Ifges(3,6) * t908 + Ifges(4,6) * t885 + pkin(2) * t767 + pkin(3) * t775 + (t899 * t928 - t900 * t931) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t897 * t877 - t896 * t878 - t938;
t800 = t824 * mrSges(7,2) + t871 * t851 - t943;
t967 = -t961 * t823 - t962 * t824 - t969 * t847 - t954 * t870 - t955 * t871 + mrSges(6,1) * t806 - mrSges(7,1) * t804 - mrSges(6,2) * t807 + mrSges(7,3) * t803 - pkin(5) * t800 + qJ(6) * (-t823 * mrSges(7,2) - t870 * t851 + t949);
t906 = (-mrSges(3,1) * t931 + mrSges(3,2) * t928) * qJD(1);
t951 = qJD(1) * t931;
t912 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t951;
t765 = m(3) * t887 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t907 + qJD(2) * t912 - t906 * t952 + t767;
t911 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t952;
t946 = -t773 * t924 + t925 * t774;
t766 = m(3) * t888 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t908 - qJD(2) * t911 + t906 * t951 + t946;
t947 = -t765 * t928 + t931 * t766;
t759 = m(2) * t914 - mrSges(2,1) * t933 - qJDD(1) * mrSges(2,2) + t947;
t788 = t926 * t794 + t966 * t796;
t940 = m(5) * t834 - t848 * mrSges(5,1) + t849 * mrSges(5,2) - t879 * t872 + t880 * t873 + t788;
t782 = m(4) * t869 - t885 * mrSges(4,1) + mrSges(4,2) * t886 - t896 * t889 + t890 * t897 + t940;
t901 = -pkin(7) * t933 + t941;
t935 = -m(3) * t901 + t908 * mrSges(3,1) - mrSges(3,2) * t907 - t911 * t952 + t912 * t951 - t782;
t777 = m(2) * t913 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t933 + t935;
t957 = t929 * t759 + t932 * t777;
t761 = t931 * t765 + t928 * t766;
t948 = t932 * t759 - t777 * t929;
t857 = Ifges(5,5) * t880 + Ifges(5,6) * t879 + Ifges(5,3) * t921;
t768 = mrSges(5,2) * t834 - mrSges(5,3) * t814 + Ifges(5,1) * t849 + Ifges(5,4) * t848 + Ifges(5,5) * t920 - pkin(9) * t788 - t926 * t784 + t966 * t786 + t879 * t857 - t921 * t858;
t769 = -mrSges(5,1) * t834 + mrSges(5,3) * t815 + Ifges(5,4) * t849 + Ifges(5,2) * t848 + Ifges(5,6) * t920 - pkin(4) * t788 - t880 * t857 + t921 * t859 - t967;
t876 = Ifges(4,5) * t897 + Ifges(4,6) * t896 + Ifges(4,3) * qJD(2);
t755 = -mrSges(4,1) * t869 + mrSges(4,3) * t838 + Ifges(4,4) * t886 + Ifges(4,2) * t885 + Ifges(4,6) * qJDD(2) - pkin(3) * t940 + pkin(8) * t945 + qJD(2) * t878 + t927 * t768 + t930 * t769 - t897 * t876;
t756 = mrSges(4,2) * t869 - mrSges(4,3) * t837 + Ifges(4,1) * t886 + Ifges(4,4) * t885 + Ifges(4,5) * qJDD(2) - pkin(8) * t775 - qJD(2) * t877 + t768 * t930 - t769 * t927 + t876 * t896;
t898 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t928 + Ifges(3,6) * t931) * qJD(1);
t751 = -mrSges(3,1) * t901 + mrSges(3,3) * t888 + Ifges(3,4) * t907 + Ifges(3,2) * t908 + Ifges(3,6) * qJDD(2) - pkin(2) * t782 + qJ(3) * t946 + qJD(2) * t900 + t925 * t755 + t924 * t756 - t898 * t952;
t753 = mrSges(3,2) * t901 - mrSges(3,3) * t887 + Ifges(3,1) * t907 + Ifges(3,4) * t908 + Ifges(3,5) * qJDD(2) - qJ(3) * t767 - qJD(2) * t899 - t755 * t924 + t756 * t925 + t898 * t951;
t939 = mrSges(2,1) * t913 - mrSges(2,2) * t914 + Ifges(2,3) * qJDD(1) + pkin(1) * t935 + pkin(7) * t947 + t931 * t751 + t928 * t753;
t754 = mrSges(2,1) * g(3) + mrSges(2,3) * t914 + t933 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t761 - t968;
t749 = -mrSges(2,2) * g(3) - mrSges(2,3) * t913 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t933 - pkin(7) * t761 - t751 * t928 + t753 * t931;
t1 = [-m(1) * g(1) + t948; -m(1) * g(2) + t957; (-m(1) - m(2)) * g(3) + t761; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t957 + t932 * t749 - t929 * t754; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t948 + t929 * t749 + t932 * t754; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t939; t939; t968; t782; -t938; t967; t800;];
tauJB  = t1;
