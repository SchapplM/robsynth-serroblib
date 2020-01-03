% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:18:48
% EndTime: 2019-12-31 19:19:21
% DurationCPUTime: 33.51s
% Computational Cost: add. (472211->331), mult. (1471615->462), div. (0->0), fcn. (1228642->14), ass. (0->158)
t853 = sin(pkin(11));
t855 = sin(pkin(5));
t856 = cos(pkin(11));
t858 = cos(pkin(5));
t861 = sin(qJ(3));
t857 = cos(pkin(6));
t865 = cos(qJ(3));
t897 = t857 * t865;
t854 = sin(pkin(6));
t902 = t854 * t865;
t871 = (-t853 * t861 + t856 * t897) * t855 + t858 * t902;
t826 = t871 * qJD(1);
t898 = t857 * t861;
t903 = t854 * t861;
t873 = t858 * t903 + (t853 * t865 + t856 * t898) * t855;
t827 = t873 * qJD(1);
t815 = -t827 * qJD(3) + t871 * qJDD(1);
t900 = t855 * t857;
t838 = (t854 * t858 + t856 * t900) * qJD(1) * pkin(8);
t862 = sin(qJ(1));
t866 = cos(qJ(1));
t850 = -g(1) * t866 - g(2) * t862;
t867 = qJD(1) ^ 2;
t906 = qJ(2) * t855;
t842 = -pkin(1) * t867 + qJDD(1) * t906 + t850;
t909 = pkin(8) * t853;
t883 = -pkin(2) * t856 - t854 * t909;
t895 = qJD(1) * t855;
t907 = pkin(8) * qJDD(1);
t878 = qJD(1) * t883 * t895 + t857 * t907;
t849 = t862 * g(1) - g(2) * t866;
t841 = qJDD(1) * pkin(1) + t867 * t906 + t849;
t891 = qJD(2) * t895;
t899 = t856 * t858;
t901 = t855 * t856;
t884 = -g(3) * t901 + t841 * t899 - 0.2e1 * t853 * t891;
t797 = (pkin(2) * qJDD(1) + qJD(1) * t838) * t858 + (-t878 * t855 - t842) * t853 + t884;
t843 = (pkin(2) * t858 - t900 * t909) * qJD(1);
t904 = t853 * t858;
t892 = t841 * t904 + (t842 + 0.2e1 * t891) * t856;
t798 = (-qJD(1) * t843 + t854 * t907) * t858 + (-g(3) * t853 + t878 * t856) * t855 + t892;
t890 = -g(3) * t858 + qJDD(2);
t806 = (-t841 + t883 * qJDD(1) + (-t838 * t856 + t843 * t853) * qJD(1)) * t855 + t890;
t773 = -t861 * t798 + (t797 * t857 + t806 * t854) * t865;
t908 = Ifges(3,3) * t858;
t905 = t853 * t855;
t774 = t797 * t898 + t865 * t798 + t806 * t903;
t813 = -mrSges(4,1) * t826 + mrSges(4,2) * t827;
t879 = -t854 * t901 + t857 * t858;
t839 = t879 * qJD(1) + qJD(3);
t823 = mrSges(4,1) * t839 - mrSges(4,3) * t827;
t836 = t879 * qJDD(1) + qJDD(3);
t814 = -pkin(3) * t826 - pkin(9) * t827;
t835 = t839 ^ 2;
t770 = -pkin(3) * t835 + pkin(9) * t836 + t814 * t826 + t774;
t782 = -t854 * t797 + t857 * t806;
t816 = t826 * qJD(3) + t873 * qJDD(1);
t772 = (-t826 * t839 - t816) * pkin(9) + (t827 * t839 - t815) * pkin(3) + t782;
t860 = sin(qJ(4));
t864 = cos(qJ(4));
t767 = t864 * t770 + t860 * t772;
t820 = -t827 * t860 + t839 * t864;
t821 = t827 * t864 + t839 * t860;
t800 = -pkin(4) * t820 - pkin(10) * t821;
t812 = qJDD(4) - t815;
t825 = qJD(4) - t826;
t824 = t825 ^ 2;
t764 = -pkin(4) * t824 + pkin(10) * t812 + t800 * t820 + t767;
t769 = -t836 * pkin(3) - t835 * pkin(9) + t827 * t814 - t773;
t792 = -qJD(4) * t821 - t816 * t860 + t836 * t864;
t793 = qJD(4) * t820 + t816 * t864 + t836 * t860;
t765 = (-t820 * t825 - t793) * pkin(10) + (t821 * t825 - t792) * pkin(4) + t769;
t859 = sin(qJ(5));
t863 = cos(qJ(5));
t761 = -t764 * t859 + t765 * t863;
t804 = -t821 * t859 + t825 * t863;
t777 = qJD(5) * t804 + t793 * t863 + t812 * t859;
t805 = t821 * t863 + t825 * t859;
t783 = -mrSges(6,1) * t804 + mrSges(6,2) * t805;
t817 = qJD(5) - t820;
t784 = -mrSges(6,2) * t817 + mrSges(6,3) * t804;
t790 = qJDD(5) - t792;
t758 = m(6) * t761 + mrSges(6,1) * t790 - mrSges(6,3) * t777 - t783 * t805 + t784 * t817;
t762 = t764 * t863 + t765 * t859;
t776 = -qJD(5) * t805 - t793 * t859 + t812 * t863;
t785 = mrSges(6,1) * t817 - mrSges(6,3) * t805;
t759 = m(6) * t762 - mrSges(6,2) * t790 + mrSges(6,3) * t776 + t783 * t804 - t785 * t817;
t752 = -t758 * t859 + t863 * t759;
t799 = -mrSges(5,1) * t820 + mrSges(5,2) * t821;
t808 = mrSges(5,1) * t825 - mrSges(5,3) * t821;
t750 = m(5) * t767 - mrSges(5,2) * t812 + mrSges(5,3) * t792 + t799 * t820 - t808 * t825 + t752;
t766 = -t770 * t860 + t772 * t864;
t763 = -pkin(4) * t812 - pkin(10) * t824 + t800 * t821 - t766;
t760 = -m(6) * t763 + t776 * mrSges(6,1) - mrSges(6,2) * t777 + t804 * t784 - t785 * t805;
t807 = -mrSges(5,2) * t825 + mrSges(5,3) * t820;
t756 = m(5) * t766 + mrSges(5,1) * t812 - mrSges(5,3) * t793 - t799 * t821 + t807 * t825 + t760;
t888 = t864 * t750 - t860 * t756;
t741 = m(4) * t774 - mrSges(4,2) * t836 + mrSges(4,3) * t815 + t813 * t826 - t823 * t839 + t888;
t744 = t860 * t750 + t864 * t756;
t822 = -mrSges(4,2) * t839 + mrSges(4,3) * t826;
t743 = m(4) * t782 - mrSges(4,1) * t815 + mrSges(4,2) * t816 - t822 * t826 + t827 * t823 + t744;
t751 = t758 * t863 + t759 * t859;
t870 = -m(5) * t769 + t792 * mrSges(5,1) - mrSges(5,2) * t793 + t820 * t807 - t808 * t821 - t751;
t747 = m(4) * t773 + mrSges(4,1) * t836 - mrSges(4,3) * t816 - t813 * t827 + t822 * t839 + t870;
t730 = t741 * t898 - t743 * t854 + t747 * t897;
t818 = -t842 * t853 + t884;
t887 = -mrSges(3,1) * t856 + mrSges(3,2) * t853;
t840 = t887 * t895;
t881 = -mrSges(3,2) * t858 + mrSges(3,3) * t901;
t845 = t881 * qJD(1);
t882 = mrSges(3,1) * t858 - mrSges(3,3) * t905;
t726 = m(3) * t818 + t882 * qJDD(1) + (-t840 * t905 + t845 * t858) * qJD(1) + t730;
t729 = t741 * t903 + t857 * t743 + t747 * t902;
t828 = -t841 * t855 + t890;
t844 = t882 * qJD(1);
t728 = m(3) * t828 + (t887 * qJDD(1) + (t844 * t853 - t845 * t856) * qJD(1)) * t855 + t729;
t735 = t865 * t741 - t747 * t861;
t819 = -g(3) * t905 + t892;
t734 = m(3) * t819 + t881 * qJDD(1) + (t840 * t901 - t844 * t858) * qJD(1) + t735;
t715 = t726 * t899 - t728 * t855 + t734 * t904;
t712 = m(2) * t849 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t867 + t715;
t720 = -t726 * t853 + t856 * t734;
t718 = m(2) * t850 - mrSges(2,1) * t867 - qJDD(1) * mrSges(2,2) + t720;
t896 = t866 * t712 + t862 * t718;
t714 = t726 * t901 + t858 * t728 + t734 * t905;
t889 = -t712 * t862 + t866 * t718;
t886 = Ifges(3,5) * t853 + Ifges(3,6) * t856;
t778 = Ifges(6,5) * t805 + Ifges(6,6) * t804 + Ifges(6,3) * t817;
t780 = Ifges(6,1) * t805 + Ifges(6,4) * t804 + Ifges(6,5) * t817;
t753 = -mrSges(6,1) * t763 + mrSges(6,3) * t762 + Ifges(6,4) * t777 + Ifges(6,2) * t776 + Ifges(6,6) * t790 - t778 * t805 + t780 * t817;
t779 = Ifges(6,4) * t805 + Ifges(6,2) * t804 + Ifges(6,6) * t817;
t754 = mrSges(6,2) * t763 - mrSges(6,3) * t761 + Ifges(6,1) * t777 + Ifges(6,4) * t776 + Ifges(6,5) * t790 + t778 * t804 - t779 * t817;
t786 = Ifges(5,5) * t821 + Ifges(5,6) * t820 + Ifges(5,3) * t825;
t787 = Ifges(5,4) * t821 + Ifges(5,2) * t820 + Ifges(5,6) * t825;
t736 = mrSges(5,2) * t769 - mrSges(5,3) * t766 + Ifges(5,1) * t793 + Ifges(5,4) * t792 + Ifges(5,5) * t812 - pkin(10) * t751 - t753 * t859 + t754 * t863 + t786 * t820 - t787 * t825;
t788 = Ifges(5,1) * t821 + Ifges(5,4) * t820 + Ifges(5,5) * t825;
t869 = mrSges(6,1) * t761 - mrSges(6,2) * t762 + Ifges(6,5) * t777 + Ifges(6,6) * t776 + Ifges(6,3) * t790 + t779 * t805 - t780 * t804;
t737 = -mrSges(5,1) * t769 + mrSges(5,3) * t767 + Ifges(5,4) * t793 + Ifges(5,2) * t792 + Ifges(5,6) * t812 - pkin(4) * t751 - t786 * t821 + t788 * t825 - t869;
t809 = Ifges(4,5) * t827 + Ifges(4,6) * t826 + Ifges(4,3) * t839;
t810 = Ifges(4,4) * t827 + Ifges(4,2) * t826 + Ifges(4,6) * t839;
t722 = mrSges(4,2) * t782 - mrSges(4,3) * t773 + Ifges(4,1) * t816 + Ifges(4,4) * t815 + Ifges(4,5) * t836 - pkin(9) * t744 + t736 * t864 - t737 * t860 + t809 * t826 - t810 * t839;
t811 = Ifges(4,1) * t827 + Ifges(4,4) * t826 + Ifges(4,5) * t839;
t868 = mrSges(5,1) * t766 - mrSges(5,2) * t767 + Ifges(5,5) * t793 + Ifges(5,6) * t792 + Ifges(5,3) * t812 + pkin(4) * t760 + pkin(10) * t752 + t863 * t753 + t859 * t754 + t821 * t787 - t820 * t788;
t723 = -mrSges(4,1) * t782 + mrSges(4,3) * t774 + Ifges(4,4) * t816 + Ifges(4,2) * t815 + Ifges(4,6) * t836 - pkin(3) * t744 - t827 * t809 + t839 * t811 - t868;
t877 = pkin(8) * t735 + t722 * t861 + t723 * t865;
t876 = Ifges(3,5) * t858 + (Ifges(3,1) * t853 + Ifges(3,4) * t856) * t855;
t875 = Ifges(3,6) * t858 + (Ifges(3,4) * t853 + Ifges(3,2) * t856) * t855;
t721 = mrSges(4,1) * t773 - mrSges(4,2) * t774 + Ifges(4,5) * t816 + Ifges(4,6) * t815 + Ifges(4,3) * t836 + pkin(3) * t870 + pkin(9) * t888 + t860 * t736 + t864 * t737 + t827 * t810 - t826 * t811;
t832 = t875 * qJD(1);
t833 = t876 * qJD(1);
t706 = qJDD(1) * t908 + mrSges(3,1) * t818 - mrSges(3,2) * t819 + pkin(2) * t730 + t721 * t857 + t877 * t854 + (t886 * qJDD(1) + (t832 * t853 - t833 * t856) * qJD(1)) * t855;
t831 = (t886 * t855 + t908) * qJD(1);
t708 = -mrSges(3,1) * t828 + mrSges(3,3) * t819 - pkin(2) * t729 - t721 * t854 + (-t831 * t905 + t833 * t858) * qJD(1) + t877 * t857 + t875 * qJDD(1);
t710 = mrSges(3,2) * t828 - mrSges(3,3) * t818 + t722 * t865 - t723 * t861 + (t831 * t901 - t832 * t858) * qJD(1) + (-t729 * t854 - t730 * t857) * pkin(8) + t876 * qJDD(1);
t874 = mrSges(2,1) * t849 - mrSges(2,2) * t850 + Ifges(2,3) * qJDD(1) + pkin(1) * t715 + t858 * t706 + t708 * t901 + t710 * t905 + t720 * t906;
t704 = -mrSges(2,2) * g(3) - mrSges(2,3) * t849 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t867 - t708 * t853 + t710 * t856 + (-t714 * t855 - t715 * t858) * qJ(2);
t703 = mrSges(2,1) * g(3) + mrSges(2,3) * t850 + Ifges(2,5) * t867 + Ifges(2,6) * qJDD(1) - pkin(1) * t714 - t706 * t855 + (qJ(2) * t720 + t708 * t856 + t710 * t853) * t858;
t1 = [-m(1) * g(1) + t889; -m(1) * g(2) + t896; (-m(1) - m(2)) * g(3) + t714; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t896 - t862 * t703 + t866 * t704; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t889 + t866 * t703 + t862 * t704; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t874; t874; t728; t721; t868; t869;];
tauJB = t1;
