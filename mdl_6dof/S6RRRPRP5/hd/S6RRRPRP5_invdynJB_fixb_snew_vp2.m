% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP5
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
% Datum: 2019-05-07 07:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:52:32
% EndTime: 2019-05-07 07:52:51
% DurationCPUTime: 18.73s
% Computational Cost: add. (298764->362), mult. (624788->442), div. (0->0), fcn. (447823->10), ass. (0->143)
t976 = Ifges(6,1) + Ifges(7,1);
t965 = Ifges(6,4) - Ifges(7,5);
t974 = Ifges(7,4) + Ifges(6,5);
t975 = Ifges(6,2) + Ifges(7,3);
t973 = Ifges(6,6) - Ifges(7,6);
t972 = -Ifges(6,3) - Ifges(7,2);
t936 = sin(qJ(3));
t939 = cos(qJ(3));
t937 = sin(qJ(2));
t960 = qJD(1) * t937;
t915 = t939 * qJD(2) - t936 * t960;
t916 = t936 * qJD(2) + t939 * t960;
t933 = sin(pkin(10));
t934 = cos(pkin(10));
t890 = t934 * t915 - t933 * t916;
t891 = t933 * t915 + t934 * t916;
t935 = sin(qJ(5));
t967 = cos(qJ(5));
t866 = -t967 * t890 + t935 * t891;
t867 = t935 * t890 + t967 * t891;
t940 = cos(qJ(2));
t959 = t940 * qJD(1);
t928 = qJD(3) - t959;
t927 = qJD(5) + t928;
t971 = t866 * t975 - t867 * t965 - t927 * t973;
t970 = -t965 * t866 + t867 * t976 + t974 * t927;
t938 = sin(qJ(1));
t941 = cos(qJ(1));
t924 = t938 * g(1) - t941 * g(2);
t943 = qJD(1) ^ 2;
t908 = -qJDD(1) * pkin(1) - t943 * pkin(7) - t924;
t958 = qJD(1) * qJD(2);
t956 = t940 * t958;
t919 = t937 * qJDD(1) + t956;
t929 = t937 * t958;
t920 = t940 * qJDD(1) - t929;
t871 = (-t919 - t956) * pkin(8) + (-t920 + t929) * pkin(2) + t908;
t925 = -t941 * g(1) - t938 * g(2);
t909 = -t943 * pkin(1) + qJDD(1) * pkin(7) + t925;
t897 = -t937 * g(3) + t940 * t909;
t918 = (-pkin(2) * t940 - pkin(8) * t937) * qJD(1);
t942 = qJD(2) ^ 2;
t874 = -t942 * pkin(2) + qJDD(2) * pkin(8) + t918 * t959 + t897;
t849 = t939 * t871 - t936 * t874;
t888 = t915 * qJD(3) + t936 * qJDD(2) + t939 * t919;
t914 = qJDD(3) - t920;
t832 = (t915 * t928 - t888) * qJ(4) + (t915 * t916 + t914) * pkin(3) + t849;
t850 = t936 * t871 + t939 * t874;
t887 = -t916 * qJD(3) + t939 * qJDD(2) - t936 * t919;
t894 = t928 * pkin(3) - t916 * qJ(4);
t913 = t915 ^ 2;
t834 = -t913 * pkin(3) + t887 * qJ(4) - t928 * t894 + t850;
t814 = -0.2e1 * qJD(4) * t891 + t934 * t832 - t933 * t834;
t863 = t933 * t887 + t934 * t888;
t811 = (t890 * t928 - t863) * pkin(9) + (t890 * t891 + t914) * pkin(4) + t814;
t815 = 0.2e1 * qJD(4) * t890 + t933 * t832 + t934 * t834;
t862 = t934 * t887 - t933 * t888;
t877 = t928 * pkin(4) - t891 * pkin(9);
t889 = t890 ^ 2;
t813 = -t889 * pkin(4) + t862 * pkin(9) - t928 * t877 + t815;
t807 = t935 * t811 + t967 * t813;
t826 = t867 * qJD(5) - t967 * t862 + t935 * t863;
t855 = t927 * mrSges(6,1) - t867 * mrSges(6,3);
t910 = qJDD(5) + t914;
t845 = t866 * pkin(5) - t867 * qJ(6);
t926 = t927 ^ 2;
t803 = -t926 * pkin(5) + t910 * qJ(6) + 0.2e1 * qJD(6) * t927 - t866 * t845 + t807;
t856 = -t927 * mrSges(7,1) + t867 * mrSges(7,2);
t957 = m(7) * t803 + t910 * mrSges(7,3) + t927 * t856;
t846 = t866 * mrSges(7,1) - t867 * mrSges(7,3);
t961 = -t866 * mrSges(6,1) - t867 * mrSges(6,2) - t846;
t966 = -mrSges(6,3) - mrSges(7,2);
t789 = m(6) * t807 - t910 * mrSges(6,2) + t966 * t826 - t927 * t855 + t961 * t866 + t957;
t806 = t967 * t811 - t935 * t813;
t827 = -t866 * qJD(5) + t935 * t862 + t967 * t863;
t854 = -t927 * mrSges(6,2) - t866 * mrSges(6,3);
t804 = -t910 * pkin(5) - t926 * qJ(6) + t867 * t845 + qJDD(6) - t806;
t853 = -t866 * mrSges(7,2) + t927 * mrSges(7,3);
t951 = -m(7) * t804 + t910 * mrSges(7,1) + t927 * t853;
t791 = m(6) * t806 + t910 * mrSges(6,1) + t966 * t827 + t927 * t854 + t961 * t867 + t951;
t784 = t935 * t789 + t967 * t791;
t868 = -t890 * mrSges(5,1) + t891 * mrSges(5,2);
t875 = -t928 * mrSges(5,2) + t890 * mrSges(5,3);
t782 = m(5) * t814 + t914 * mrSges(5,1) - t863 * mrSges(5,3) - t891 * t868 + t928 * t875 + t784;
t876 = t928 * mrSges(5,1) - t891 * mrSges(5,3);
t952 = t967 * t789 - t935 * t791;
t783 = m(5) * t815 - t914 * mrSges(5,2) + t862 * mrSges(5,3) + t890 * t868 - t928 * t876 + t952;
t778 = t934 * t782 + t933 * t783;
t860 = Ifges(5,4) * t891 + Ifges(5,2) * t890 + Ifges(5,6) * t928;
t861 = Ifges(5,1) * t891 + Ifges(5,4) * t890 + Ifges(5,5) * t928;
t879 = Ifges(4,4) * t916 + Ifges(4,2) * t915 + Ifges(4,6) * t928;
t880 = Ifges(4,1) * t916 + Ifges(4,4) * t915 + Ifges(4,5) * t928;
t799 = t827 * mrSges(7,2) + t867 * t846 - t951;
t946 = -mrSges(6,1) * t806 + mrSges(7,1) * t804 + mrSges(6,2) * t807 - mrSges(7,3) * t803 + pkin(5) * t799 - qJ(6) * t957 + t972 * t910 + t971 * t867 + (qJ(6) * t846 - t970) * t866 - t974 * t827 + (qJ(6) * mrSges(7,2) + t973) * t826;
t969 = mrSges(4,1) * t849 + mrSges(5,1) * t814 - mrSges(4,2) * t850 - mrSges(5,2) * t815 + Ifges(4,5) * t888 + Ifges(5,5) * t863 + Ifges(4,6) * t887 + Ifges(5,6) * t862 + pkin(3) * t778 + pkin(4) * t784 + t891 * t860 - t890 * t861 + t916 * t879 - t915 * t880 + (Ifges(4,3) + Ifges(5,3)) * t914 - t946;
t896 = -t940 * g(3) - t937 * t909;
t873 = -qJDD(2) * pkin(2) - t942 * pkin(8) + t918 * t960 - t896;
t848 = -t887 * pkin(3) - t913 * qJ(4) + t916 * t894 + qJDD(4) + t873;
t817 = -t862 * pkin(4) - t889 * pkin(9) + t891 * t877 + t848;
t809 = t817 - 0.2e1 * qJD(6) * t867 + (t867 * t927 + t826) * pkin(5) + (t866 * t927 - t827) * qJ(6);
t800 = m(7) * t809 + t826 * mrSges(7,1) - t827 * mrSges(7,3) + t866 * t853 - t867 * t856;
t962 = t973 * t866 - t974 * t867 + t972 * t927;
t785 = -mrSges(6,1) * t817 - mrSges(7,1) * t809 + mrSges(7,2) * t803 + mrSges(6,3) * t807 - pkin(5) * t800 - t826 * t975 + t965 * t827 + t962 * t867 + t973 * t910 + t970 * t927;
t786 = mrSges(6,2) * t817 + mrSges(7,2) * t804 - mrSges(6,3) * t806 - mrSges(7,3) * t809 - qJ(6) * t800 - t965 * t826 + t827 * t976 + t962 * t866 + t974 * t910 + t971 * t927;
t859 = Ifges(5,5) * t891 + Ifges(5,6) * t890 + Ifges(5,3) * t928;
t948 = m(6) * t817 + t826 * mrSges(6,1) + t827 * mrSges(6,2) + t866 * t854 + t867 * t855 + t800;
t773 = -mrSges(5,1) * t848 + mrSges(5,3) * t815 + Ifges(5,4) * t863 + Ifges(5,2) * t862 + Ifges(5,6) * t914 - pkin(4) * t948 + pkin(9) * t952 + t967 * t785 + t935 * t786 - t891 * t859 + t928 * t861;
t774 = mrSges(5,2) * t848 - mrSges(5,3) * t814 + Ifges(5,1) * t863 + Ifges(5,4) * t862 + Ifges(5,5) * t914 - pkin(9) * t784 - t935 * t785 + t967 * t786 + t890 * t859 - t928 * t860;
t795 = m(5) * t848 - t862 * mrSges(5,1) + t863 * mrSges(5,2) - t890 * t875 + t891 * t876 + t948;
t878 = Ifges(4,5) * t916 + Ifges(4,6) * t915 + Ifges(4,3) * t928;
t953 = -t933 * t782 + t934 * t783;
t758 = -mrSges(4,1) * t873 + mrSges(4,3) * t850 + Ifges(4,4) * t888 + Ifges(4,2) * t887 + Ifges(4,6) * t914 - pkin(3) * t795 + qJ(4) * t953 + t934 * t773 + t933 * t774 - t916 * t878 + t928 * t880;
t759 = mrSges(4,2) * t873 - mrSges(4,3) * t849 + Ifges(4,1) * t888 + Ifges(4,4) * t887 + Ifges(4,5) * t914 - qJ(4) * t778 - t933 * t773 + t934 * t774 + t915 * t878 - t928 * t879;
t892 = -t915 * mrSges(4,1) + t916 * mrSges(4,2);
t893 = -t928 * mrSges(4,2) + t915 * mrSges(4,3);
t776 = m(4) * t849 + t914 * mrSges(4,1) - t888 * mrSges(4,3) - t916 * t892 + t928 * t893 + t778;
t895 = t928 * mrSges(4,1) - t916 * mrSges(4,3);
t777 = m(4) * t850 - t914 * mrSges(4,2) + t887 * mrSges(4,3) + t915 * t892 - t928 * t895 + t953;
t772 = -t936 * t776 + t939 * t777;
t794 = -m(4) * t873 + t887 * mrSges(4,1) - t888 * mrSges(4,2) + t915 * t893 - t916 * t895 - t795;
t906 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t937 + Ifges(3,2) * t940) * qJD(1);
t907 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t937 + Ifges(3,4) * t940) * qJD(1);
t968 = mrSges(3,1) * t896 - mrSges(3,2) * t897 + Ifges(3,5) * t919 + Ifges(3,6) * t920 + Ifges(3,3) * qJDD(2) + pkin(2) * t794 + pkin(8) * t772 + t939 * t758 + t936 * t759 + (t937 * t906 - t940 * t907) * qJD(1);
t917 = (-mrSges(3,1) * t940 + mrSges(3,2) * t937) * qJD(1);
t922 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t960;
t770 = m(3) * t897 - qJDD(2) * mrSges(3,2) + t920 * mrSges(3,3) - qJD(2) * t922 + t917 * t959 + t772;
t923 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t959;
t793 = m(3) * t896 + qJDD(2) * mrSges(3,1) - t919 * mrSges(3,3) + qJD(2) * t923 - t917 * t960 + t794;
t954 = t940 * t770 - t937 * t793;
t762 = m(2) * t925 - t943 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t954;
t771 = t939 * t776 + t936 * t777;
t947 = -m(3) * t908 + t920 * mrSges(3,1) - t919 * mrSges(3,2) - t922 * t960 + t923 * t959 - t771;
t766 = m(2) * t924 + qJDD(1) * mrSges(2,1) - t943 * mrSges(2,2) + t947;
t963 = t938 * t762 + t941 * t766;
t764 = t937 * t770 + t940 * t793;
t955 = t941 * t762 - t938 * t766;
t905 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t937 + Ifges(3,6) * t940) * qJD(1);
t755 = mrSges(3,2) * t908 - mrSges(3,3) * t896 + Ifges(3,1) * t919 + Ifges(3,4) * t920 + Ifges(3,5) * qJDD(2) - pkin(8) * t771 - qJD(2) * t906 - t936 * t758 + t939 * t759 + t905 * t959;
t757 = -mrSges(3,1) * t908 + mrSges(3,3) * t897 + Ifges(3,4) * t919 + Ifges(3,2) * t920 + Ifges(3,6) * qJDD(2) - pkin(2) * t771 + qJD(2) * t907 - t905 * t960 - t969;
t949 = mrSges(2,1) * t924 - mrSges(2,2) * t925 + Ifges(2,3) * qJDD(1) + pkin(1) * t947 + pkin(7) * t954 + t755 * t937 + t757 * t940;
t753 = mrSges(2,1) * g(3) + mrSges(2,3) * t925 + t943 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t764 - t968;
t752 = -mrSges(2,2) * g(3) - mrSges(2,3) * t924 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t943 - pkin(7) * t764 + t755 * t940 - t757 * t937;
t1 = [-m(1) * g(1) + t955; -m(1) * g(2) + t963; (-m(1) - m(2)) * g(3) + t764; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t963 + t752 * t941 - t753 * t938; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t955 + t938 * t752 + t941 * t753; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t949; t949; t968; t969; t795; -t946; t799;];
tauJB  = t1;
