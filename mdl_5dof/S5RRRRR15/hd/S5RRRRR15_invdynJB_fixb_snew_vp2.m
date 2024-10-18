% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR15_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR15_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:25:47
% EndTime: 2024-09-27 22:26:16
% DurationCPUTime: 13.78s
% Computational Cost: add. (700306->342), mult. (1969467->460), div. (0->0), fcn. (1595511->14), ass. (0->151)
t898 = sin(pkin(5));
t939 = pkin(8) * t898;
t897 = sin(pkin(6));
t938 = pkin(11) * t897;
t899 = cos(pkin(6));
t937 = pkin(11) * t899;
t911 = qJD(1) ^ 2;
t936 = t898 ^ 2 * t911;
t901 = sin(qJ(5));
t935 = t897 * t901;
t906 = cos(qJ(5));
t934 = t897 * t906;
t904 = sin(qJ(2));
t933 = t898 * t904;
t909 = cos(qJ(2));
t932 = t898 * t909;
t900 = cos(pkin(5));
t931 = t900 * t904;
t930 = t900 * t909;
t927 = qJD(1) * qJD(2);
t878 = (qJDD(1) * t904 + t909 * t927) * t898;
t889 = t900 * qJDD(1) + qJDD(2);
t890 = t900 * qJD(1) + qJD(2);
t905 = sin(qJ(1));
t910 = cos(qJ(1));
t886 = t905 * g(1) - t910 * g(2);
t874 = qJDD(1) * pkin(1) + t911 * t939 + t886;
t887 = -t910 * g(1) - t905 * g(2);
t875 = -t911 * pkin(1) + qJDD(1) * t939 + t887;
t920 = t874 * t930 - t904 * t875;
t831 = t889 * pkin(2) - t878 * pkin(9) + (pkin(2) * t904 * t936 + (pkin(9) * qJD(1) * t890 - g(3)) * t898) * t909 + t920;
t854 = -g(3) * t933 + t874 * t931 + t909 * t875;
t928 = qJD(1) * t898;
t925 = t904 * t928;
t877 = t890 * pkin(2) - pkin(9) * t925;
t879 = (qJDD(1) * t909 - t904 * t927) * t898;
t926 = t909 ^ 2 * t936;
t832 = -pkin(2) * t926 + t879 * pkin(9) - t890 * t877 + t854;
t903 = sin(qJ(3));
t908 = cos(qJ(3));
t811 = t908 * t831 - t903 * t832;
t870 = (-t903 * t904 + t908 * t909) * t928;
t845 = t870 * qJD(3) + t908 * t878 + t903 * t879;
t871 = (t903 * t909 + t904 * t908) * t928;
t885 = qJDD(3) + t889;
t888 = qJD(3) + t890;
t801 = (t870 * t888 - t845) * pkin(10) + (t870 * t871 + t885) * pkin(3) + t811;
t812 = t903 * t831 + t908 * t832;
t844 = -t871 * qJD(3) - t903 * t878 + t908 * t879;
t860 = t888 * pkin(3) - t871 * pkin(10);
t868 = t870 ^ 2;
t803 = -t868 * pkin(3) + t844 * pkin(10) - t888 * t860 + t812;
t902 = sin(qJ(4));
t907 = cos(qJ(4));
t796 = t902 * t801 + t907 * t803;
t855 = t907 * t870 - t902 * t871;
t856 = t902 * t870 + t907 * t871;
t833 = -t855 * pkin(4) - t856 * t938;
t884 = qJD(4) + t888;
t843 = t884 * pkin(4) - t856 * t937;
t819 = -t856 * qJD(4) + t907 * t844 - t902 * t845;
t883 = qJDD(4) + t885;
t918 = t819 * t899 + t883 * t897;
t792 = t918 * pkin(11) + t855 * t833 - t884 * t843 + t796;
t795 = t907 * t801 - t902 * t803;
t820 = t855 * qJD(4) + t902 * t844 + t907 * t845;
t917 = t855 * t899 + t884 * t897;
t838 = t917 * pkin(11);
t791 = t883 * pkin(4) - t820 * t937 - t856 * t833 + t884 * t838 + t795;
t864 = -t900 * g(3) - t898 * t874;
t836 = -t879 * pkin(2) - pkin(9) * t926 + t877 * t925 + t864;
t809 = -t844 * pkin(3) - t868 * pkin(10) + t871 * t860 + t836;
t793 = -t819 * pkin(4) - t820 * t938 - t855 * t838 + t856 * t843 + t809;
t919 = t791 * t899 + t793 * t897;
t788 = -t901 * t792 + t919 * t906;
t822 = -t901 * t856 + t917 * t906;
t798 = t822 * qJD(5) + t906 * t820 + t918 * t901;
t823 = t906 * t856 + t917 * t901;
t807 = -t822 * mrSges(6,1) + t823 * mrSges(6,2);
t813 = -t897 * t819 + t899 * t883 + qJDD(5);
t839 = -t897 * t855 + t899 * t884 + qJD(5);
t814 = -t839 * mrSges(6,2) + t822 * mrSges(6,3);
t784 = m(6) * t788 + t813 * mrSges(6,1) - t798 * mrSges(6,3) - t823 * t807 + t839 * t814;
t789 = t906 * t792 + t919 * t901;
t797 = -t823 * qJD(5) - t901 * t820 + t918 * t906;
t815 = t839 * mrSges(6,1) - t823 * mrSges(6,3);
t785 = m(6) * t789 - t813 * mrSges(6,2) + t797 * mrSges(6,3) + t822 * t807 - t839 * t815;
t790 = -t897 * t791 + t899 * t793;
t787 = m(6) * t790 - t797 * mrSges(6,1) + t798 * mrSges(6,2) - t822 * t814 + t823 * t815;
t767 = -t897 * t787 + (t784 * t906 + t785 * t901) * t899;
t834 = -t855 * mrSges(5,1) + t856 * mrSges(5,2);
t846 = -t884 * mrSges(5,2) + t855 * mrSges(5,3);
t764 = m(5) * t795 + t883 * mrSges(5,1) - t820 * mrSges(5,3) - t856 * t834 + t884 * t846 + t767;
t772 = -t901 * t784 + t906 * t785;
t847 = t884 * mrSges(5,1) - t856 * mrSges(5,3);
t770 = m(5) * t796 - t883 * mrSges(5,2) + t819 * mrSges(5,3) + t855 * t834 - t884 * t847 + t772;
t760 = t907 * t764 + t902 * t770;
t857 = -t870 * mrSges(4,1) + t871 * mrSges(4,2);
t858 = -t888 * mrSges(4,2) + t870 * mrSges(4,3);
t757 = m(4) * t811 + t885 * mrSges(4,1) - t845 * mrSges(4,3) - t871 * t857 + t888 * t858 + t760;
t859 = t888 * mrSges(4,1) - t871 * mrSges(4,3);
t921 = -t902 * t764 + t907 * t770;
t758 = m(4) * t812 - t885 * mrSges(4,2) + t844 * mrSges(4,3) + t870 * t857 - t888 * t859 + t921;
t751 = t908 * t757 + t903 * t758;
t853 = -g(3) * t932 + t920;
t924 = t909 * t928;
t873 = -t890 * mrSges(3,2) + mrSges(3,3) * t924;
t876 = (-mrSges(3,1) * t909 + mrSges(3,2) * t904) * t928;
t749 = m(3) * t853 + t889 * mrSges(3,1) - t878 * mrSges(3,3) + t890 * t873 - t876 * t925 + t751;
t872 = t890 * mrSges(3,1) - mrSges(3,3) * t925;
t922 = -t903 * t757 + t908 * t758;
t750 = m(3) * t854 - t889 * mrSges(3,2) + t879 * mrSges(3,3) - t890 * t872 + t876 * t924 + t922;
t766 = t784 * t934 + t785 * t935 + t899 * t787;
t915 = m(5) * t809 - t819 * mrSges(5,1) + t820 * mrSges(5,2) - t855 * t846 + t856 * t847 + t766;
t913 = m(4) * t836 - t844 * mrSges(4,1) + t845 * mrSges(4,2) - t870 * t858 + t871 * t859 + t915;
t762 = t878 * mrSges(3,2) - t879 * mrSges(3,1) + m(3) * t864 + (t872 * t904 - t873 * t909) * t928 + t913;
t738 = t749 * t930 + t750 * t931 - t898 * t762;
t735 = m(2) * t886 + qJDD(1) * mrSges(2,1) - t911 * mrSges(2,2) + t738;
t743 = -t904 * t749 + t909 * t750;
t741 = m(2) * t887 - t911 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t743;
t929 = t910 * t735 + t905 * t741;
t737 = t749 * t932 + t750 * t933 + t900 * t762;
t923 = -t905 * t735 + t910 * t741;
t805 = Ifges(6,4) * t823 + Ifges(6,2) * t822 + Ifges(6,6) * t839;
t806 = Ifges(6,1) * t823 + Ifges(6,4) * t822 + Ifges(6,5) * t839;
t774 = mrSges(6,1) * t788 - mrSges(6,2) * t789 + Ifges(6,5) * t798 + Ifges(6,6) * t797 + Ifges(6,3) * t813 + t823 * t805 - t822 * t806;
t804 = Ifges(6,5) * t823 + Ifges(6,6) * t822 + Ifges(6,3) * t839;
t777 = -mrSges(6,1) * t790 + mrSges(6,3) * t789 + Ifges(6,4) * t798 + Ifges(6,2) * t797 + Ifges(6,6) * t813 - t823 * t804 + t839 * t806;
t778 = mrSges(6,2) * t790 - mrSges(6,3) * t788 + Ifges(6,1) * t798 + Ifges(6,4) * t797 + Ifges(6,5) * t813 + t822 * t804 - t839 * t805;
t824 = Ifges(5,5) * t856 + Ifges(5,6) * t855 + Ifges(5,3) * t884;
t826 = Ifges(5,1) * t856 + Ifges(5,4) * t855 + Ifges(5,5) * t884;
t752 = -mrSges(5,1) * t809 + mrSges(5,3) * t796 + Ifges(5,4) * t820 + Ifges(5,2) * t819 + Ifges(5,6) * t883 - pkin(4) * t766 - t897 * t774 - t856 * t824 + t884 * t826 + (pkin(11) * t772 + t777 * t906 + t778 * t901) * t899;
t825 = Ifges(5,4) * t856 + Ifges(5,2) * t855 + Ifges(5,6) * t884;
t753 = mrSges(5,2) * t809 - mrSges(5,3) * t795 + Ifges(5,1) * t820 + Ifges(5,4) * t819 + Ifges(5,5) * t883 - t901 * t777 + t906 * t778 + t855 * t824 - t884 * t825 + (-t766 * t897 - t767 * t899) * pkin(11);
t848 = Ifges(4,5) * t871 + Ifges(4,6) * t870 + Ifges(4,3) * t888;
t850 = Ifges(4,1) * t871 + Ifges(4,4) * t870 + Ifges(4,5) * t888;
t730 = -mrSges(4,1) * t836 + mrSges(4,3) * t812 + Ifges(4,4) * t845 + Ifges(4,2) * t844 + Ifges(4,6) * t885 - pkin(3) * t915 + pkin(10) * t921 + t907 * t752 + t902 * t753 - t871 * t848 + t888 * t850;
t849 = Ifges(4,4) * t871 + Ifges(4,2) * t870 + Ifges(4,6) * t888;
t733 = mrSges(4,2) * t836 - mrSges(4,3) * t811 + Ifges(4,1) * t845 + Ifges(4,4) * t844 + Ifges(4,5) * t885 - pkin(10) * t760 - t902 * t752 + t907 * t753 + t870 * t848 - t888 * t849;
t861 = Ifges(3,3) * t890 + (Ifges(3,5) * t904 + Ifges(3,6) * t909) * t928;
t863 = Ifges(3,5) * t890 + (Ifges(3,1) * t904 + Ifges(3,4) * t909) * t928;
t727 = -mrSges(3,1) * t864 + mrSges(3,3) * t854 + Ifges(3,4) * t878 + Ifges(3,2) * t879 + Ifges(3,6) * t889 - pkin(2) * t913 + pkin(9) * t922 + t908 * t730 + t903 * t733 - t861 * t925 + t890 * t863;
t862 = Ifges(3,6) * t890 + (Ifges(3,4) * t904 + Ifges(3,2) * t909) * t928;
t729 = mrSges(3,2) * t864 - mrSges(3,3) * t853 + Ifges(3,1) * t878 + Ifges(3,4) * t879 + Ifges(3,5) * t889 - pkin(9) * t751 - t903 * t730 + t908 * t733 + t861 * t924 - t890 * t862;
t914 = mrSges(5,1) * t795 - mrSges(5,2) * t796 + Ifges(5,5) * t820 + Ifges(5,6) * t819 + Ifges(5,3) * t883 + pkin(4) * t767 + t772 * t938 + t899 * t774 + t777 * t934 + t778 * t935 + t856 * t825 - t855 * t826;
t912 = mrSges(4,1) * t811 - mrSges(4,2) * t812 + Ifges(4,5) * t845 + Ifges(4,6) * t844 + Ifges(4,3) * t885 + pkin(3) * t760 + t871 * t849 - t870 * t850 + t914;
t732 = (t862 * t904 - t863 * t909) * t928 + pkin(2) * t751 + Ifges(3,3) * t889 + Ifges(3,5) * t878 + Ifges(3,6) * t879 + mrSges(3,1) * t853 - mrSges(3,2) * t854 + t912;
t916 = mrSges(2,1) * t886 - mrSges(2,2) * t887 + Ifges(2,3) * qJDD(1) + pkin(1) * t738 + t727 * t932 + t729 * t933 + t900 * t732 + t743 * t939;
t725 = -mrSges(2,2) * g(3) - mrSges(2,3) * t886 + Ifges(2,5) * qJDD(1) - t911 * Ifges(2,6) - t904 * t727 + t909 * t729 + (-t737 * t898 - t738 * t900) * pkin(8);
t724 = mrSges(2,1) * g(3) + mrSges(2,3) * t887 + t911 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t737 - t898 * t732 + (pkin(8) * t743 + t727 * t909 + t729 * t904) * t900;
t1 = [-m(1) * g(1) + t923; -m(1) * g(2) + t929; (-m(1) - m(2)) * g(3) + t737; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t929 - t905 * t724 + t910 * t725; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t923 + t910 * t724 + t905 * t725; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t916; t916; t732; t912; t914; t774;];
tauJB = t1;
