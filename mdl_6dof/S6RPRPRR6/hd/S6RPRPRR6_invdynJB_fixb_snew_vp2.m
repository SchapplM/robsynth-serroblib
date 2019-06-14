% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:02:53
% EndTime: 2019-05-05 19:03:21
% DurationCPUTime: 28.61s
% Computational Cost: add. (467551->366), mult. (1154958->460), div. (0->0), fcn. (888919->12), ass. (0->156)
t942 = qJD(1) ^ 2;
t975 = cos(qJ(3));
t933 = cos(pkin(10));
t974 = pkin(2) * t933;
t931 = sin(pkin(10));
t973 = mrSges(3,2) * t931;
t928 = t933 ^ 2;
t972 = t928 * t942;
t937 = sin(qJ(1));
t940 = cos(qJ(1));
t918 = -t940 * g(1) - t937 * g(2);
t913 = -t942 * pkin(1) + qJDD(1) * qJ(2) + t918;
t967 = qJD(1) * qJD(2);
t963 = -t933 * g(3) - 0.2e1 * t931 * t967;
t878 = (-pkin(7) * qJDD(1) + t942 * t974 - t913) * t931 + t963;
t899 = -t931 * g(3) + (t913 + 0.2e1 * t967) * t933;
t965 = qJDD(1) * t933;
t884 = -pkin(2) * t972 + pkin(7) * t965 + t899;
t936 = sin(qJ(3));
t861 = t936 * t878 + t975 * t884;
t964 = t933 * t975;
t968 = t931 * qJD(1);
t911 = -qJD(1) * t964 + t936 * t968;
t951 = t975 * t931 + t933 * t936;
t912 = t951 * qJD(1);
t891 = t911 * pkin(3) - t912 * qJ(4);
t941 = qJD(3) ^ 2;
t847 = -t941 * pkin(3) + qJDD(3) * qJ(4) - t911 * t891 + t861;
t927 = t931 ^ 2;
t917 = t937 * g(1) - t940 * g(2);
t957 = qJDD(2) - t917;
t895 = (-pkin(1) - t974) * qJDD(1) + (-qJ(2) + (-t927 - t928) * pkin(7)) * t942 + t957;
t966 = qJDD(1) * t931;
t969 = t912 * qJD(3);
t896 = -qJDD(1) * t964 + t936 * t966 + t969;
t970 = t911 * qJD(3);
t897 = t951 * qJDD(1) - t970;
t852 = (-t897 + t970) * qJ(4) + (t896 + t969) * pkin(3) + t895;
t930 = sin(pkin(11));
t932 = cos(pkin(11));
t904 = t930 * qJD(3) + t932 * t912;
t834 = -0.2e1 * qJD(4) * t904 - t930 * t847 + t932 * t852;
t883 = t930 * qJDD(3) + t932 * t897;
t903 = t932 * qJD(3) - t930 * t912;
t825 = (t911 * t903 - t883) * pkin(8) + (t903 * t904 + t896) * pkin(4) + t834;
t835 = 0.2e1 * qJD(4) * t903 + t932 * t847 + t930 * t852;
t881 = t911 * pkin(4) - t904 * pkin(8);
t882 = t932 * qJDD(3) - t930 * t897;
t902 = t903 ^ 2;
t827 = -t902 * pkin(4) + t882 * pkin(8) - t911 * t881 + t835;
t935 = sin(qJ(5));
t939 = cos(qJ(5));
t820 = t939 * t825 - t935 * t827;
t871 = t939 * t903 - t935 * t904;
t846 = t871 * qJD(5) + t935 * t882 + t939 * t883;
t872 = t935 * t903 + t939 * t904;
t894 = qJDD(5) + t896;
t909 = qJD(5) + t911;
t818 = (t871 * t909 - t846) * pkin(9) + (t871 * t872 + t894) * pkin(5) + t820;
t821 = t935 * t825 + t939 * t827;
t845 = -t872 * qJD(5) + t939 * t882 - t935 * t883;
t864 = t909 * pkin(5) - t872 * pkin(9);
t870 = t871 ^ 2;
t819 = -t870 * pkin(5) + t845 * pkin(9) - t909 * t864 + t821;
t934 = sin(qJ(6));
t938 = cos(qJ(6));
t816 = t938 * t818 - t934 * t819;
t857 = t938 * t871 - t934 * t872;
t832 = t857 * qJD(6) + t934 * t845 + t938 * t846;
t858 = t934 * t871 + t938 * t872;
t841 = -t857 * mrSges(7,1) + t858 * mrSges(7,2);
t908 = qJD(6) + t909;
t850 = -t908 * mrSges(7,2) + t857 * mrSges(7,3);
t890 = qJDD(6) + t894;
t810 = m(7) * t816 + t890 * mrSges(7,1) - t832 * mrSges(7,3) - t858 * t841 + t908 * t850;
t817 = t934 * t818 + t938 * t819;
t831 = -t858 * qJD(6) + t938 * t845 - t934 * t846;
t851 = t908 * mrSges(7,1) - t858 * mrSges(7,3);
t811 = m(7) * t817 - t890 * mrSges(7,2) + t831 * mrSges(7,3) + t857 * t841 - t908 * t851;
t804 = t938 * t810 + t934 * t811;
t859 = -t871 * mrSges(6,1) + t872 * mrSges(6,2);
t862 = -t909 * mrSges(6,2) + t871 * mrSges(6,3);
t802 = m(6) * t820 + t894 * mrSges(6,1) - t846 * mrSges(6,3) - t872 * t859 + t909 * t862 + t804;
t863 = t909 * mrSges(6,1) - t872 * mrSges(6,3);
t958 = -t934 * t810 + t938 * t811;
t803 = m(6) * t821 - t894 * mrSges(6,2) + t845 * mrSges(6,3) + t871 * t859 - t909 * t863 + t958;
t798 = t939 * t802 + t935 * t803;
t873 = -t903 * mrSges(5,1) + t904 * mrSges(5,2);
t956 = -t911 * mrSges(5,2) + t903 * mrSges(5,3);
t796 = m(5) * t834 + t896 * mrSges(5,1) - t883 * mrSges(5,3) - t904 * t873 + t911 * t956 + t798;
t880 = t911 * mrSges(5,1) - t904 * mrSges(5,3);
t959 = -t935 * t802 + t939 * t803;
t797 = m(5) * t835 - t896 * mrSges(5,2) + t882 * mrSges(5,3) + t903 * t873 - t911 * t880 + t959;
t790 = -t930 * t796 + t932 * t797;
t892 = t911 * mrSges(4,1) + t912 * mrSges(4,2);
t906 = qJD(3) * mrSges(4,1) - t912 * mrSges(4,3);
t788 = m(4) * t861 - qJDD(3) * mrSges(4,2) - t896 * mrSges(4,3) - qJD(3) * t906 - t911 * t892 + t790;
t860 = t975 * t878 - t936 * t884;
t844 = -qJDD(3) * pkin(3) - t941 * qJ(4) + t912 * t891 + qJDD(4) - t860;
t836 = -t882 * pkin(4) - t902 * pkin(8) + t904 * t881 + t844;
t822 = -t845 * pkin(5) - t870 * pkin(9) + t872 * t864 + t836;
t949 = m(7) * t822 - t831 * mrSges(7,1) + t832 * mrSges(7,2) - t857 * t850 + t858 * t851;
t945 = m(6) * t836 - t845 * mrSges(6,1) + t846 * mrSges(6,2) - t871 * t862 + t872 * t863 + t949;
t814 = m(5) * t844 - t882 * mrSges(5,1) + t883 * mrSges(5,2) + t904 * t880 - t903 * t956 + t945;
t905 = -qJD(3) * mrSges(4,2) - t911 * mrSges(4,3);
t813 = m(4) * t860 + qJDD(3) * mrSges(4,1) - t897 * mrSges(4,3) + qJD(3) * t905 - t912 * t892 - t814;
t781 = t936 * t788 + t975 * t813;
t898 = -t931 * t913 + t963;
t952 = mrSges(3,3) * qJDD(1) + t942 * (-mrSges(3,1) * t933 + t973);
t779 = m(3) * t898 - t952 * t931 + t781;
t960 = t975 * t788 - t936 * t813;
t780 = m(3) * t899 + t952 * t933 + t960;
t961 = -t931 * t779 + t933 * t780;
t770 = m(2) * t918 - t942 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t961;
t910 = -qJDD(1) * pkin(1) - t942 * qJ(2) + t957;
t789 = t932 * t796 + t930 * t797;
t947 = m(4) * t895 + t896 * mrSges(4,1) + t897 * mrSges(4,2) + t911 * t905 + t912 * t906 + t789;
t946 = -m(3) * t910 + mrSges(3,1) * t965 - t947 + (t927 * t942 + t972) * mrSges(3,3);
t783 = -t942 * mrSges(2,2) + m(2) * t917 + t946 + (mrSges(2,1) - t973) * qJDD(1);
t971 = t937 * t770 + t940 * t783;
t772 = t933 * t779 + t931 * t780;
t962 = t940 * t770 - t937 * t783;
t955 = Ifges(3,1) * t931 + Ifges(3,4) * t933;
t954 = Ifges(3,4) * t931 + Ifges(3,2) * t933;
t953 = Ifges(3,5) * t931 + Ifges(3,6) * t933;
t837 = Ifges(7,5) * t858 + Ifges(7,6) * t857 + Ifges(7,3) * t908;
t839 = Ifges(7,1) * t858 + Ifges(7,4) * t857 + Ifges(7,5) * t908;
t805 = -mrSges(7,1) * t822 + mrSges(7,3) * t817 + Ifges(7,4) * t832 + Ifges(7,2) * t831 + Ifges(7,6) * t890 - t858 * t837 + t908 * t839;
t838 = Ifges(7,4) * t858 + Ifges(7,2) * t857 + Ifges(7,6) * t908;
t806 = mrSges(7,2) * t822 - mrSges(7,3) * t816 + Ifges(7,1) * t832 + Ifges(7,4) * t831 + Ifges(7,5) * t890 + t857 * t837 - t908 * t838;
t853 = Ifges(6,5) * t872 + Ifges(6,6) * t871 + Ifges(6,3) * t909;
t855 = Ifges(6,1) * t872 + Ifges(6,4) * t871 + Ifges(6,5) * t909;
t791 = -mrSges(6,1) * t836 + mrSges(6,3) * t821 + Ifges(6,4) * t846 + Ifges(6,2) * t845 + Ifges(6,6) * t894 - pkin(5) * t949 + pkin(9) * t958 + t938 * t805 + t934 * t806 - t872 * t853 + t909 * t855;
t854 = Ifges(6,4) * t872 + Ifges(6,2) * t871 + Ifges(6,6) * t909;
t792 = mrSges(6,2) * t836 - mrSges(6,3) * t820 + Ifges(6,1) * t846 + Ifges(6,4) * t845 + Ifges(6,5) * t894 - pkin(9) * t804 - t934 * t805 + t938 * t806 + t871 * t853 - t909 * t854;
t865 = Ifges(5,5) * t904 + Ifges(5,6) * t903 + Ifges(5,3) * t911;
t867 = Ifges(5,1) * t904 + Ifges(5,4) * t903 + Ifges(5,5) * t911;
t774 = -mrSges(5,1) * t844 + mrSges(5,3) * t835 + Ifges(5,4) * t883 + Ifges(5,2) * t882 + Ifges(5,6) * t896 - pkin(4) * t945 + pkin(8) * t959 + t939 * t791 + t935 * t792 - t904 * t865 + t911 * t867;
t866 = Ifges(5,4) * t904 + Ifges(5,2) * t903 + Ifges(5,6) * t911;
t775 = mrSges(5,2) * t844 - mrSges(5,3) * t834 + Ifges(5,1) * t883 + Ifges(5,4) * t882 + Ifges(5,5) * t896 - pkin(8) * t798 - t935 * t791 + t939 * t792 + t903 * t865 - t911 * t866;
t885 = Ifges(4,5) * t912 - Ifges(4,6) * t911 + Ifges(4,3) * qJD(3);
t886 = Ifges(4,4) * t912 - Ifges(4,2) * t911 + Ifges(4,6) * qJD(3);
t767 = mrSges(4,2) * t895 - mrSges(4,3) * t860 + Ifges(4,1) * t897 - Ifges(4,4) * t896 + Ifges(4,5) * qJDD(3) - qJ(4) * t789 - qJD(3) * t886 - t930 * t774 + t932 * t775 - t911 * t885;
t887 = Ifges(4,1) * t912 - Ifges(4,4) * t911 + Ifges(4,5) * qJD(3);
t948 = -mrSges(7,1) * t816 + mrSges(7,2) * t817 - Ifges(7,5) * t832 - Ifges(7,6) * t831 - Ifges(7,3) * t890 - t858 * t838 + t857 * t839;
t943 = mrSges(6,1) * t820 - mrSges(6,2) * t821 + Ifges(6,5) * t846 + Ifges(6,6) * t845 + Ifges(6,3) * t894 + pkin(5) * t804 + t872 * t854 - t871 * t855 - t948;
t773 = -t912 * t885 + t903 * t867 - t904 * t866 - mrSges(4,1) * t895 + Ifges(4,4) * t897 - Ifges(5,6) * t882 - Ifges(5,5) * t883 + qJD(3) * t887 + mrSges(4,3) * t861 + mrSges(5,2) * t835 - mrSges(5,1) * t834 - pkin(4) * t798 - pkin(3) * t789 + Ifges(4,6) * qJDD(3) + (-Ifges(5,3) - Ifges(4,2)) * t896 - t943;
t915 = t953 * qJD(1);
t763 = -mrSges(3,1) * t910 + mrSges(3,3) * t899 - pkin(2) * t947 + pkin(7) * t960 + t954 * qJDD(1) + t936 * t767 + t975 * t773 - t915 * t968;
t766 = t933 * qJD(1) * t915 + mrSges(3,2) * t910 - mrSges(3,3) * t898 - pkin(7) * t781 + t955 * qJDD(1) + t975 * t767 - t936 * t773;
t785 = mrSges(3,2) * t966 - t946;
t950 = mrSges(2,1) * t917 - mrSges(2,2) * t918 + Ifges(2,3) * qJDD(1) - pkin(1) * t785 + qJ(2) * t961 + t933 * t763 + t931 * t766;
t944 = mrSges(4,1) * t860 - mrSges(4,2) * t861 + Ifges(4,5) * t897 - Ifges(4,6) * t896 + Ifges(4,3) * qJDD(3) - pkin(3) * t814 + qJ(4) * t790 + t932 * t774 + t930 * t775 + t912 * t886 + t911 * t887;
t764 = mrSges(2,1) * g(3) + mrSges(2,3) * t918 - mrSges(3,1) * t898 + mrSges(3,2) * t899 - pkin(1) * t772 + (Ifges(2,6) - t953) * qJDD(1) - pkin(2) * t781 - t944 + (-t931 * t954 + t933 * t955 + Ifges(2,5)) * t942;
t761 = -mrSges(2,2) * g(3) - mrSges(2,3) * t917 + Ifges(2,5) * qJDD(1) - t942 * Ifges(2,6) - qJ(2) * t772 - t931 * t763 + t933 * t766;
t1 = [-m(1) * g(1) + t962; -m(1) * g(2) + t971; (-m(1) - m(2)) * g(3) + t772; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t971 + t940 * t761 - t937 * t764; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t962 + t937 * t761 + t940 * t764; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t950; t950; t785; t944; t814; t943; -t948;];
tauJB  = t1;
