% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:25:01
% EndTime: 2019-05-05 10:25:27
% DurationCPUTime: 25.54s
% Computational Cost: add. (439727->344), mult. (915244->442), div. (0->0), fcn. (688119->14), ass. (0->148)
t885 = sin(pkin(12));
t887 = cos(pkin(12));
t869 = g(1) * t885 - g(2) * t887;
t870 = -g(1) * t887 - g(2) * t885;
t884 = -g(3) + qJDD(1);
t886 = sin(pkin(6));
t888 = cos(pkin(6));
t893 = sin(qJ(2));
t898 = cos(qJ(2));
t841 = -t893 * t870 + (t869 * t888 + t884 * t886) * t898;
t923 = t888 * t893;
t924 = t886 * t893;
t842 = t869 * t923 + t898 * t870 + t884 * t924;
t899 = qJD(2) ^ 2;
t836 = -pkin(2) * t899 + qJDD(2) * pkin(8) + t842;
t853 = -t869 * t886 + t884 * t888;
t892 = sin(qJ(3));
t897 = cos(qJ(3));
t822 = -t892 * t836 + t897 * t853;
t919 = qJD(2) * qJD(3);
t918 = t897 * t919;
t867 = qJDD(2) * t892 + t918;
t809 = (-t867 + t918) * pkin(9) + (t892 * t897 * t899 + qJDD(3)) * pkin(3) + t822;
t823 = t897 * t836 + t892 * t853;
t868 = qJDD(2) * t897 - t892 * t919;
t921 = qJD(2) * t892;
t874 = qJD(3) * pkin(3) - pkin(9) * t921;
t883 = t897 ^ 2;
t811 = -pkin(3) * t883 * t899 + pkin(9) * t868 - qJD(3) * t874 + t823;
t891 = sin(qJ(4));
t896 = cos(qJ(4));
t787 = t896 * t809 - t891 * t811;
t859 = (-t891 * t892 + t896 * t897) * qJD(2);
t832 = qJD(4) * t859 + t867 * t896 + t868 * t891;
t860 = (t891 * t897 + t892 * t896) * qJD(2);
t881 = qJDD(3) + qJDD(4);
t882 = qJD(3) + qJD(4);
t783 = (t859 * t882 - t832) * pkin(10) + (t859 * t860 + t881) * pkin(4) + t787;
t788 = t891 * t809 + t896 * t811;
t831 = -qJD(4) * t860 - t867 * t891 + t868 * t896;
t852 = pkin(4) * t882 - pkin(10) * t860;
t855 = t859 ^ 2;
t785 = -pkin(4) * t855 + pkin(10) * t831 - t852 * t882 + t788;
t890 = sin(qJ(5));
t895 = cos(qJ(5));
t780 = t890 * t783 + t895 * t785;
t846 = t859 * t890 + t860 * t895;
t803 = -qJD(5) * t846 + t831 * t895 - t832 * t890;
t845 = t859 * t895 - t860 * t890;
t819 = -mrSges(6,1) * t845 + mrSges(6,2) * t846;
t879 = qJD(5) + t882;
t834 = mrSges(6,1) * t879 - mrSges(6,3) * t846;
t878 = qJDD(5) + t881;
t820 = -pkin(5) * t845 - pkin(11) * t846;
t877 = t879 ^ 2;
t777 = -pkin(5) * t877 + pkin(11) * t878 + t820 * t845 + t780;
t906 = -qJDD(2) * pkin(2) - t841;
t821 = -t868 * pkin(3) + t874 * t921 + (-pkin(9) * t883 - pkin(8)) * t899 + t906;
t793 = -t831 * pkin(4) - t855 * pkin(10) + t860 * t852 + t821;
t804 = qJD(5) * t845 + t831 * t890 + t832 * t895;
t781 = (-t845 * t879 - t804) * pkin(11) + (t846 * t879 - t803) * pkin(5) + t793;
t889 = sin(qJ(6));
t894 = cos(qJ(6));
t774 = -t777 * t889 + t781 * t894;
t825 = -t846 * t889 + t879 * t894;
t791 = qJD(6) * t825 + t804 * t894 + t878 * t889;
t801 = qJDD(6) - t803;
t826 = t846 * t894 + t879 * t889;
t810 = -mrSges(7,1) * t825 + mrSges(7,2) * t826;
t840 = qJD(6) - t845;
t812 = -mrSges(7,2) * t840 + mrSges(7,3) * t825;
t770 = m(7) * t774 + mrSges(7,1) * t801 - mrSges(7,3) * t791 - t810 * t826 + t812 * t840;
t775 = t777 * t894 + t781 * t889;
t790 = -qJD(6) * t826 - t804 * t889 + t878 * t894;
t813 = mrSges(7,1) * t840 - mrSges(7,3) * t826;
t771 = m(7) * t775 - mrSges(7,2) * t801 + mrSges(7,3) * t790 + t810 * t825 - t813 * t840;
t913 = -t770 * t889 + t894 * t771;
t757 = m(6) * t780 - mrSges(6,2) * t878 + mrSges(6,3) * t803 + t819 * t845 - t834 * t879 + t913;
t779 = t783 * t895 - t785 * t890;
t833 = -mrSges(6,2) * t879 + mrSges(6,3) * t845;
t776 = -pkin(5) * t878 - pkin(11) * t877 + t820 * t846 - t779;
t907 = -m(7) * t776 + t790 * mrSges(7,1) - mrSges(7,2) * t791 + t825 * t812 - t813 * t826;
t766 = m(6) * t779 + mrSges(6,1) * t878 - mrSges(6,3) * t804 - t819 * t846 + t833 * t879 + t907;
t751 = t890 * t757 + t895 * t766;
t847 = -mrSges(5,1) * t859 + mrSges(5,2) * t860;
t850 = -mrSges(5,2) * t882 + mrSges(5,3) * t859;
t748 = m(5) * t787 + mrSges(5,1) * t881 - mrSges(5,3) * t832 - t847 * t860 + t850 * t882 + t751;
t851 = mrSges(5,1) * t882 - mrSges(5,3) * t860;
t914 = t895 * t757 - t766 * t890;
t749 = m(5) * t788 - mrSges(5,2) * t881 + mrSges(5,3) * t831 + t847 * t859 - t851 * t882 + t914;
t742 = t896 * t748 + t891 * t749;
t857 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t892 + Ifges(4,2) * t897) * qJD(2);
t858 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t892 + Ifges(4,4) * t897) * qJD(2);
t838 = Ifges(5,4) * t860 + Ifges(5,2) * t859 + Ifges(5,6) * t882;
t839 = Ifges(5,1) * t860 + Ifges(5,4) * t859 + Ifges(5,5) * t882;
t794 = Ifges(7,5) * t826 + Ifges(7,6) * t825 + Ifges(7,3) * t840;
t796 = Ifges(7,1) * t826 + Ifges(7,4) * t825 + Ifges(7,5) * t840;
t763 = -mrSges(7,1) * t776 + mrSges(7,3) * t775 + Ifges(7,4) * t791 + Ifges(7,2) * t790 + Ifges(7,6) * t801 - t794 * t826 + t796 * t840;
t795 = Ifges(7,4) * t826 + Ifges(7,2) * t825 + Ifges(7,6) * t840;
t764 = mrSges(7,2) * t776 - mrSges(7,3) * t774 + Ifges(7,1) * t791 + Ifges(7,4) * t790 + Ifges(7,5) * t801 + t794 * t825 - t795 * t840;
t815 = Ifges(6,4) * t846 + Ifges(6,2) * t845 + Ifges(6,6) * t879;
t816 = Ifges(6,1) * t846 + Ifges(6,4) * t845 + Ifges(6,5) * t879;
t905 = -mrSges(6,1) * t779 + mrSges(6,2) * t780 - Ifges(6,5) * t804 - Ifges(6,6) * t803 - Ifges(6,3) * t878 - pkin(5) * t907 - pkin(11) * t913 - t894 * t763 - t889 * t764 - t846 * t815 + t845 * t816;
t902 = -mrSges(5,1) * t787 + mrSges(5,2) * t788 - Ifges(5,5) * t832 - Ifges(5,6) * t831 - Ifges(5,3) * t881 - pkin(4) * t751 - t860 * t838 + t859 * t839 + t905;
t926 = mrSges(4,1) * t822 - mrSges(4,2) * t823 + Ifges(4,5) * t867 + Ifges(4,6) * t868 + Ifges(4,3) * qJDD(3) + pkin(3) * t742 + (t857 * t892 - t858 * t897) * qJD(2) - t902;
t835 = -t899 * pkin(8) + t906;
t871 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t921;
t920 = qJD(2) * t897;
t872 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t920;
t759 = t894 * t770 + t889 * t771;
t909 = m(6) * t793 - t803 * mrSges(6,1) + t804 * mrSges(6,2) - t845 * t833 + t846 * t834 + t759;
t904 = m(5) * t821 - t831 * mrSges(5,1) + mrSges(5,2) * t832 - t859 * t850 + t851 * t860 + t909;
t901 = -m(4) * t835 + t868 * mrSges(4,1) - mrSges(4,2) * t867 - t871 * t921 + t872 * t920 - t904;
t754 = m(3) * t841 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t899 + t901;
t925 = t754 * t898;
t866 = (-mrSges(4,1) * t897 + mrSges(4,2) * t892) * qJD(2);
t740 = m(4) * t822 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t867 + qJD(3) * t872 - t866 * t921 + t742;
t915 = -t748 * t891 + t896 * t749;
t741 = m(4) * t823 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t868 - qJD(3) * t871 + t866 * t920 + t915;
t916 = -t740 * t892 + t897 * t741;
t732 = m(3) * t842 - mrSges(3,1) * t899 - qJDD(2) * mrSges(3,2) + t916;
t735 = t897 * t740 + t892 * t741;
t734 = m(3) * t853 + t735;
t723 = t732 * t923 - t734 * t886 + t888 * t925;
t721 = m(2) * t869 + t723;
t727 = t898 * t732 - t754 * t893;
t726 = m(2) * t870 + t727;
t922 = t887 * t721 + t885 * t726;
t722 = t732 * t924 + t888 * t734 + t886 * t925;
t917 = -t721 * t885 + t887 * t726;
t912 = m(2) * t884 + t722;
t814 = Ifges(6,5) * t846 + Ifges(6,6) * t845 + Ifges(6,3) * t879;
t743 = mrSges(6,2) * t793 - mrSges(6,3) * t779 + Ifges(6,1) * t804 + Ifges(6,4) * t803 + Ifges(6,5) * t878 - pkin(11) * t759 - t763 * t889 + t764 * t894 + t814 * t845 - t815 * t879;
t903 = mrSges(7,1) * t774 - mrSges(7,2) * t775 + Ifges(7,5) * t791 + Ifges(7,6) * t790 + Ifges(7,3) * t801 + t795 * t826 - t796 * t825;
t744 = -mrSges(6,1) * t793 + mrSges(6,3) * t780 + Ifges(6,4) * t804 + Ifges(6,2) * t803 + Ifges(6,6) * t878 - pkin(5) * t759 - t814 * t846 + t816 * t879 - t903;
t837 = Ifges(5,5) * t860 + Ifges(5,6) * t859 + Ifges(5,3) * t882;
t728 = -mrSges(5,1) * t821 + mrSges(5,3) * t788 + Ifges(5,4) * t832 + Ifges(5,2) * t831 + Ifges(5,6) * t881 - pkin(4) * t909 + pkin(10) * t914 + t890 * t743 + t895 * t744 - t860 * t837 + t882 * t839;
t736 = mrSges(5,2) * t821 - mrSges(5,3) * t787 + Ifges(5,1) * t832 + Ifges(5,4) * t831 + Ifges(5,5) * t881 - pkin(10) * t751 + t743 * t895 - t744 * t890 + t837 * t859 - t838 * t882;
t856 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t892 + Ifges(4,6) * t897) * qJD(2);
t717 = -mrSges(4,1) * t835 + mrSges(4,3) * t823 + Ifges(4,4) * t867 + Ifges(4,2) * t868 + Ifges(4,6) * qJDD(3) - pkin(3) * t904 + pkin(9) * t915 + qJD(3) * t858 + t896 * t728 + t891 * t736 - t856 * t921;
t718 = mrSges(4,2) * t835 - mrSges(4,3) * t822 + Ifges(4,1) * t867 + Ifges(4,4) * t868 + Ifges(4,5) * qJDD(3) - pkin(9) * t742 - qJD(3) * t857 - t728 * t891 + t736 * t896 + t856 * t920;
t716 = mrSges(3,2) * t853 - mrSges(3,3) * t841 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t899 - pkin(8) * t735 - t717 * t892 + t718 * t897;
t719 = -mrSges(3,1) * t853 + mrSges(3,3) * t842 + t899 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t735 - t926;
t908 = pkin(7) * t727 + t716 * t893 + t719 * t898;
t715 = mrSges(3,1) * t841 - mrSges(3,2) * t842 + Ifges(3,3) * qJDD(2) + pkin(2) * t901 + pkin(8) * t916 + t897 * t717 + t892 * t718;
t714 = mrSges(2,2) * t884 - mrSges(2,3) * t869 + t898 * t716 - t893 * t719 + (-t722 * t886 - t723 * t888) * pkin(7);
t713 = -mrSges(2,1) * t884 + mrSges(2,3) * t870 - pkin(1) * t722 - t886 * t715 + t888 * t908;
t1 = [-m(1) * g(1) + t917; -m(1) * g(2) + t922; -m(1) * g(3) + t912; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t922 - t885 * t713 + t887 * t714; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t917 + t887 * t713 + t885 * t714; -mrSges(1,1) * g(2) + mrSges(2,1) * t869 + mrSges(1,2) * g(1) - mrSges(2,2) * t870 + pkin(1) * t723 + t888 * t715 + t886 * t908; t912; t715; t926; -t902; -t905; t903;];
tauJB  = t1;
