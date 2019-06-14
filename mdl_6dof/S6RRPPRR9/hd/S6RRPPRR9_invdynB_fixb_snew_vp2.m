% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 11:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:32:14
% EndTime: 2019-05-06 11:32:27
% DurationCPUTime: 7.88s
% Computational Cost: add. (85222->366), mult. (199025->446), div. (0->0), fcn. (132937->10), ass. (0->157)
t914 = -2 * qJD(3);
t913 = -2 * qJD(4);
t912 = Ifges(3,1) + Ifges(5,3) + Ifges(4,2);
t887 = Ifges(3,4) + Ifges(4,6) - Ifges(5,6);
t886 = Ifges(3,5) - Ifges(4,4) + Ifges(5,5);
t911 = -Ifges(3,2) - Ifges(5,2) - Ifges(4,3);
t885 = Ifges(3,6) - Ifges(4,5) - Ifges(5,4);
t910 = Ifges(3,3) + Ifges(4,1) + Ifges(5,1);
t853 = sin(qJ(1));
t857 = cos(qJ(1));
t835 = t853 * g(1) - g(2) * t857;
t848 = sin(pkin(6));
t858 = qJD(1) ^ 2;
t812 = pkin(8) * t848 * t858 + qJDD(1) * pkin(1) + t835;
t836 = -g(1) * t857 - g(2) * t853;
t888 = qJDD(1) * t848;
t813 = -pkin(1) * t858 + pkin(8) * t888 + t836;
t852 = sin(qJ(2));
t849 = cos(pkin(6));
t856 = cos(qJ(2));
t897 = t849 * t856;
t899 = t848 * t856;
t764 = -g(3) * t899 + t812 * t897 - t852 * t813;
t892 = qJD(1) * t848;
t814 = (-pkin(2) * t856 - qJ(3) * t852) * t892;
t842 = qJD(1) * t849 + qJD(2);
t840 = t842 ^ 2;
t841 = qJDD(1) * t849 + qJDD(2);
t878 = t852 * t892;
t750 = -t841 * pkin(2) - t840 * qJ(3) + t814 * t878 + qJDD(3) - t764;
t891 = qJD(1) * t856;
t819 = (qJD(2) * t891 + qJDD(1) * t852) * t848;
t877 = t848 * t891;
t873 = t842 * t877;
t901 = t848 ^ 2 * t858;
t882 = t856 * t901;
t737 = t750 - (t852 * t882 + t841) * qJ(4) - (-t819 + t873) * pkin(3) + t842 * t913;
t898 = t849 * t852;
t895 = t812 * t898 + t856 * t813;
t909 = pkin(2) * t840 - t841 * qJ(3) - t814 * t877 + t842 * t914 - t895;
t908 = g(3) * t849;
t907 = mrSges(3,1) - mrSges(4,2);
t906 = -mrSges(3,3) - mrSges(5,1);
t905 = qJ(3) * t842;
t807 = pkin(3) * t878 - qJ(4) * t842;
t904 = t807 * t852;
t810 = mrSges(5,1) * t877 + mrSges(5,2) * t842;
t903 = t810 * t856;
t902 = t812 * t848;
t900 = t848 * t852;
t781 = -t902 - t908;
t805 = mrSges(3,1) * t842 - mrSges(3,3) * t878;
t806 = -mrSges(3,2) * t842 + mrSges(3,3) * t877;
t820 = -qJD(2) * t878 + t856 * t888;
t869 = t878 * t914 - t908 + (t842 * t878 - t820) * pkin(2);
t745 = -t902 + (-t819 - t873) * qJ(3) + t869;
t809 = -mrSges(4,1) * t877 - mrSges(4,3) * t842;
t866 = -t820 * qJ(4) + t877 * t913 + t869;
t847 = t856 ^ 2;
t883 = t847 * t901;
t736 = -pkin(3) * t883 - qJ(3) * t819 + (-t812 + (-t856 * t905 - t904) * qJD(1)) * t848 + t866;
t818 = pkin(4) * t877 - pkin(9) * t842;
t846 = t852 ^ 2;
t730 = (-pkin(3) * t847 - pkin(4) * t846) * t901 + (pkin(9) - qJ(3)) * t819 + (-t812 + (-t904 + (-t818 - t905) * t856) * qJD(1)) * t848 + t866;
t860 = -qJ(4) * t883 + t842 * t807 + qJDD(4) - t909;
t734 = -pkin(9) * t841 + (pkin(3) + pkin(4)) * t820 + (pkin(9) * t882 + (pkin(4) * qJD(1) * t842 - g(3)) * t848) * t852 + t860;
t851 = sin(qJ(5));
t855 = cos(qJ(5));
t727 = t855 * t730 + t851 * t734;
t796 = t842 * t855 + t851 * t878;
t762 = -qJD(5) * t796 + t819 * t855 - t841 * t851;
t795 = -t842 * t851 + t855 * t878;
t766 = -mrSges(6,1) * t795 + mrSges(6,2) * t796;
t827 = qJD(5) + t877;
t771 = mrSges(6,1) * t827 - mrSges(6,3) * t796;
t803 = qJDD(5) + t820;
t767 = -pkin(5) * t795 - pkin(10) * t796;
t824 = t827 ^ 2;
t725 = -pkin(5) * t824 + pkin(10) * t803 + t767 * t795 + t727;
t733 = -pkin(9) * t846 * t901 - pkin(4) * t819 + t842 * t818 - t737;
t763 = qJD(5) * t795 + t819 * t851 + t841 * t855;
t728 = t733 + (t796 * t827 - t762) * pkin(5) + (-t795 * t827 - t763) * pkin(10);
t850 = sin(qJ(6));
t854 = cos(qJ(6));
t722 = -t725 * t850 + t728 * t854;
t768 = -t796 * t850 + t827 * t854;
t742 = qJD(6) * t768 + t763 * t854 + t803 * t850;
t769 = t796 * t854 + t827 * t850;
t751 = -mrSges(7,1) * t768 + mrSges(7,2) * t769;
t794 = qJD(6) - t795;
t752 = -mrSges(7,2) * t794 + mrSges(7,3) * t768;
t759 = qJDD(6) - t762;
t720 = m(7) * t722 + mrSges(7,1) * t759 - mrSges(7,3) * t742 - t751 * t769 + t752 * t794;
t723 = t725 * t854 + t728 * t850;
t741 = -qJD(6) * t769 - t763 * t850 + t803 * t854;
t753 = mrSges(7,1) * t794 - mrSges(7,3) * t769;
t721 = m(7) * t723 - mrSges(7,2) * t759 + mrSges(7,3) * t741 + t751 * t768 - t753 * t794;
t874 = -t720 * t850 + t854 * t721;
t711 = m(6) * t727 - mrSges(6,2) * t803 + mrSges(6,3) * t762 + t766 * t795 - t771 * t827 + t874;
t726 = -t730 * t851 + t734 * t855;
t770 = -mrSges(6,2) * t827 + mrSges(6,3) * t795;
t724 = -pkin(5) * t803 - pkin(10) * t824 + t767 * t796 - t726;
t864 = -m(7) * t724 + t741 * mrSges(7,1) - mrSges(7,2) * t742 + t768 * t752 - t753 * t769;
t716 = m(6) * t726 + mrSges(6,1) * t803 - mrSges(6,3) * t763 - t766 * t796 + t770 * t827 + t864;
t875 = t855 * t711 - t716 * t851;
t872 = m(5) * t736 - t820 * mrSges(5,3) + t875;
t865 = m(4) * t745 - t819 * mrSges(4,3) + t809 * t877 + t872;
t808 = mrSges(5,1) * t878 - mrSges(5,3) * t842;
t811 = mrSges(4,1) * t878 + mrSges(4,2) * t842;
t894 = -t808 - t811;
t698 = m(3) * t781 - t907 * t820 + (mrSges(3,2) - mrSges(5,2)) * t819 + ((-t806 - t810) * t856 + (t805 + t894) * t852) * t892 + t865;
t884 = g(3) * t900;
t765 = -t884 + t895;
t816 = (-mrSges(3,1) * t856 + mrSges(3,2) * t852) * t892;
t744 = t884 + t909;
t815 = (mrSges(4,2) * t856 - mrSges(4,3) * t852) * t892;
t704 = t851 * t711 + t855 * t716;
t739 = pkin(3) * t820 + t860 - t884;
t817 = (-mrSges(5,2) * t852 - mrSges(5,3) * t856) * t892;
t871 = -m(5) * t739 - t841 * mrSges(5,2) - t842 * t808 - t817 * t877 - t704;
t862 = -m(4) * t744 + t841 * mrSges(4,3) + t842 * t811 + t815 * t877 - t871;
t702 = (mrSges(4,1) - t906) * t820 + m(3) * t765 - mrSges(3,2) * t841 - t805 * t842 + t816 * t877 + t862;
t712 = t854 * t720 + t850 * t721;
t868 = -m(6) * t733 + t762 * mrSges(6,1) - t763 * mrSges(6,2) + t795 * t770 - t796 * t771 - t712;
t861 = m(5) * t737 - t841 * mrSges(5,3) - t842 * t810 + t868;
t859 = -m(4) * t750 - t819 * mrSges(4,1) - t861;
t893 = -t815 - t817;
t707 = m(3) * t764 + t906 * t819 + t907 * t841 + (-t816 + t893) * t878 + (t806 - t809) * t842 + t859;
t690 = -t698 * t848 + t702 * t898 + t707 * t897;
t688 = m(2) * t835 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t858 + t690;
t694 = t856 * t702 - t707 * t852;
t693 = m(2) * t836 - mrSges(2,1) * t858 - qJDD(1) * mrSges(2,2) + t694;
t896 = t857 * t688 + t853 * t693;
t689 = t849 * t698 + t702 * t900 + t707 * t899;
t881 = (t886 * t852 + t885 * t856) * t892 + t910 * t842;
t880 = (-t887 * t852 + t911 * t856) * t892 - t885 * t842;
t879 = (-t912 * t852 - t887 * t856) * t892 - t886 * t842;
t876 = -t688 * t853 + t857 * t693;
t746 = Ifges(7,5) * t769 + Ifges(7,6) * t768 + Ifges(7,3) * t794;
t748 = Ifges(7,1) * t769 + Ifges(7,4) * t768 + Ifges(7,5) * t794;
t713 = -mrSges(7,1) * t724 + mrSges(7,3) * t723 + Ifges(7,4) * t742 + Ifges(7,2) * t741 + Ifges(7,6) * t759 - t746 * t769 + t748 * t794;
t747 = Ifges(7,4) * t769 + Ifges(7,2) * t768 + Ifges(7,6) * t794;
t714 = mrSges(7,2) * t724 - mrSges(7,3) * t722 + Ifges(7,1) * t742 + Ifges(7,4) * t741 + Ifges(7,5) * t759 + t746 * t768 - t747 * t794;
t754 = Ifges(6,5) * t796 + Ifges(6,6) * t795 + Ifges(6,3) * t827;
t755 = Ifges(6,4) * t796 + Ifges(6,2) * t795 + Ifges(6,6) * t827;
t695 = mrSges(6,2) * t733 - mrSges(6,3) * t726 + Ifges(6,1) * t763 + Ifges(6,4) * t762 + Ifges(6,5) * t803 - pkin(10) * t712 - t713 * t850 + t714 * t854 + t754 * t795 - t755 * t827;
t756 = Ifges(6,1) * t796 + Ifges(6,4) * t795 + Ifges(6,5) * t827;
t696 = -mrSges(6,1) * t733 - mrSges(7,1) * t722 + mrSges(7,2) * t723 + mrSges(6,3) * t727 + Ifges(6,4) * t763 - Ifges(7,5) * t742 + Ifges(6,2) * t762 + Ifges(6,6) * t803 - Ifges(7,6) * t741 - Ifges(7,3) * t759 - pkin(5) * t712 - t747 * t769 + t748 * t768 - t754 * t796 + t756 * t827;
t703 = mrSges(4,2) * t820 - mrSges(5,2) * t819 + (t894 * t852 - t903) * t892 + t865;
t708 = t819 * mrSges(5,1) + t817 * t878 + t861;
t685 = mrSges(3,2) * t781 - mrSges(3,3) * t764 + pkin(4) * t868 - mrSges(4,3) * t745 + mrSges(4,1) * t750 - mrSges(5,2) * t736 + mrSges(5,1) * t737 + pkin(3) * t708 - qJ(3) * t703 + pkin(9) * t875 + t851 * t695 + t855 * t696 + t880 * t842 + t886 * t841 + t887 * t820 + t912 * t819 + t881 * t877;
t686 = (pkin(3) * mrSges(5,1) - t911) * t820 - t795 * t756 + t796 * t755 - mrSges(3,1) * t781 + Ifges(6,6) * t762 + Ifges(6,5) * t763 + mrSges(3,3) * t765 - mrSges(4,1) * t744 + mrSges(4,2) * t745 - mrSges(5,3) * t736 + mrSges(5,1) * t739 + mrSges(6,1) * t726 - mrSges(6,2) * t727 + pkin(4) * t704 - pkin(2) * t703 + t850 * t714 + pkin(5) * t864 + t854 * t713 + Ifges(6,3) * t803 + (qJ(4) * t903 + (qJ(4) * t808 - t881) * t852) * t892 + t885 * t841 + (qJ(4) * mrSges(5,2) + t887) * t819 - t879 * t842 - pkin(3) * t871 - qJ(4) * t872 + pkin(10) * t874;
t867 = pkin(8) * t694 + t685 * t852 + t686 * t856;
t684 = mrSges(3,1) * t764 - mrSges(3,2) * t765 - mrSges(4,3) * t744 + mrSges(4,2) * t750 - mrSges(5,3) * t737 + mrSges(5,2) * t739 - qJ(4) * t708 - pkin(9) * t704 - t851 * t696 + pkin(2) * (-t842 * t809 + t859) + qJ(3) * t862 + t855 * t695 + (-mrSges(4,2) * pkin(2) + t910) * t841 + (-mrSges(5,1) * pkin(2) + t886) * t819 + (qJ(3) * (mrSges(4,1) + mrSges(5,1)) + t885) * t820 + (t879 * t856 + (pkin(2) * t893 - t880) * t852) * t892;
t683 = -mrSges(2,2) * g(3) - mrSges(2,3) * t835 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t858 + t685 * t856 - t686 * t852 + (-t689 * t848 - t690 * t849) * pkin(8);
t682 = mrSges(2,1) * g(3) + mrSges(2,3) * t836 + Ifges(2,5) * t858 + Ifges(2,6) * qJDD(1) - pkin(1) * t689 - t684 * t848 + t867 * t849;
t1 = [-m(1) * g(1) + t876; -m(1) * g(2) + t896; (-m(1) - m(2)) * g(3) + t689; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t896 - t853 * t682 + t857 * t683; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t876 + t857 * t682 + t853 * t683; -mrSges(1,1) * g(2) + mrSges(2,1) * t835 + mrSges(1,2) * g(1) - mrSges(2,2) * t836 + Ifges(2,3) * qJDD(1) + pkin(1) * t690 + t684 * t849 + t867 * t848;];
tauB  = t1;
