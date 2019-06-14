% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:27:19
% EndTime: 2019-05-06 04:27:31
% DurationCPUTime: 11.04s
% Computational Cost: add. (189408->340), mult. (375657->411), div. (0->0), fcn. (254079->10), ass. (0->139)
t852 = sin(qJ(1));
t857 = cos(qJ(1));
t833 = -t857 * g(1) - t852 * g(2);
t887 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t833;
t886 = -pkin(1) - pkin(7);
t885 = mrSges(2,1) - mrSges(3,2);
t884 = -Ifges(3,4) + Ifges(2,5);
t883 = (Ifges(3,5) - Ifges(2,6));
t832 = t852 * g(1) - t857 * g(2);
t859 = qJD(1) ^ 2;
t870 = -t859 * qJ(2) + qJDD(2) - t832;
t801 = qJDD(1) * t886 + t870;
t851 = sin(qJ(3));
t856 = cos(qJ(3));
t794 = -g(3) * t856 + t851 * t801;
t824 = (mrSges(4,1) * t851 + mrSges(4,2) * t856) * qJD(1);
t880 = qJD(1) * qJD(3);
t836 = t856 * t880;
t826 = -t851 * qJDD(1) - t836;
t881 = qJD(1) * t856;
t831 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t881;
t838 = t851 * qJD(1);
t800 = t859 * t886 - t887;
t878 = t851 * t880;
t827 = qJDD(1) * t856 - t878;
t773 = (-t827 + t878) * pkin(8) + (-t826 + t836) * pkin(3) + t800;
t825 = (pkin(3) * t851 - pkin(8) * t856) * qJD(1);
t858 = qJD(3) ^ 2;
t776 = -pkin(3) * t858 + qJDD(3) * pkin(8) - t825 * t838 + t794;
t850 = sin(qJ(4));
t855 = cos(qJ(4));
t754 = t855 * t773 - t850 * t776;
t822 = qJD(3) * t855 - t850 * t881;
t787 = qJD(4) * t822 + qJDD(3) * t850 + t827 * t855;
t821 = qJDD(4) - t826;
t823 = qJD(3) * t850 + t855 * t881;
t835 = t838 + qJD(4);
t744 = (t822 * t835 - t787) * pkin(9) + (t822 * t823 + t821) * pkin(4) + t754;
t755 = t850 * t773 + t855 * t776;
t786 = -qJD(4) * t823 + qJDD(3) * t855 - t827 * t850;
t799 = pkin(4) * t835 - pkin(9) * t823;
t820 = t822 ^ 2;
t746 = -pkin(4) * t820 + pkin(9) * t786 - t799 * t835 + t755;
t849 = sin(qJ(5));
t854 = cos(qJ(5));
t731 = t854 * t744 - t849 * t746;
t789 = t822 * t854 - t823 * t849;
t761 = qJD(5) * t789 + t786 * t849 + t787 * t854;
t790 = t822 * t849 + t823 * t854;
t814 = qJDD(5) + t821;
t834 = qJD(5) + t835;
t728 = (t789 * t834 - t761) * pkin(10) + (t789 * t790 + t814) * pkin(5) + t731;
t732 = t849 * t744 + t854 * t746;
t760 = -qJD(5) * t790 + t786 * t854 - t787 * t849;
t779 = pkin(5) * t834 - pkin(10) * t790;
t788 = t789 ^ 2;
t729 = -pkin(5) * t788 + pkin(10) * t760 - t779 * t834 + t732;
t848 = sin(qJ(6));
t853 = cos(qJ(6));
t726 = t728 * t853 - t729 * t848;
t768 = t789 * t853 - t790 * t848;
t740 = qJD(6) * t768 + t760 * t848 + t761 * t853;
t769 = t789 * t848 + t790 * t853;
t752 = -mrSges(7,1) * t768 + mrSges(7,2) * t769;
t829 = qJD(6) + t834;
t762 = -mrSges(7,2) * t829 + mrSges(7,3) * t768;
t810 = qJDD(6) + t814;
t723 = m(7) * t726 + mrSges(7,1) * t810 - mrSges(7,3) * t740 - t752 * t769 + t762 * t829;
t727 = t728 * t848 + t729 * t853;
t739 = -qJD(6) * t769 + t760 * t853 - t761 * t848;
t763 = mrSges(7,1) * t829 - mrSges(7,3) * t769;
t724 = m(7) * t727 - mrSges(7,2) * t810 + mrSges(7,3) * t739 + t752 * t768 - t763 * t829;
t715 = t853 * t723 + t848 * t724;
t770 = -mrSges(6,1) * t789 + mrSges(6,2) * t790;
t777 = -mrSges(6,2) * t834 + mrSges(6,3) * t789;
t712 = m(6) * t731 + mrSges(6,1) * t814 - mrSges(6,3) * t761 - t770 * t790 + t777 * t834 + t715;
t778 = mrSges(6,1) * t834 - mrSges(6,3) * t790;
t873 = -t723 * t848 + t853 * t724;
t713 = m(6) * t732 - mrSges(6,2) * t814 + mrSges(6,3) * t760 + t770 * t789 - t778 * t834 + t873;
t708 = t854 * t712 + t849 * t713;
t792 = -mrSges(5,1) * t822 + mrSges(5,2) * t823;
t795 = -mrSges(5,2) * t835 + mrSges(5,3) * t822;
t706 = m(5) * t754 + mrSges(5,1) * t821 - mrSges(5,3) * t787 - t792 * t823 + t795 * t835 + t708;
t796 = mrSges(5,1) * t835 - mrSges(5,3) * t823;
t874 = -t712 * t849 + t854 * t713;
t707 = m(5) * t755 - mrSges(5,2) * t821 + mrSges(5,3) * t786 + t792 * t822 - t796 * t835 + t874;
t875 = -t706 * t850 + t855 * t707;
t698 = m(4) * t794 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t826 - qJD(3) * t831 - t824 * t838 + t875;
t793 = g(3) * t851 + t801 * t856;
t830 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t838;
t775 = -qJDD(3) * pkin(3) - pkin(8) * t858 + t825 * t881 - t793;
t753 = -pkin(4) * t786 - pkin(9) * t820 + t823 * t799 + t775;
t734 = -pkin(5) * t760 - pkin(10) * t788 + t779 * t790 + t753;
t872 = m(7) * t734 - t739 * mrSges(7,1) + t740 * mrSges(7,2) - t768 * t762 + t769 * t763;
t864 = m(6) * t753 - t760 * mrSges(6,1) + t761 * mrSges(6,2) - t789 * t777 + t790 * t778 + t872;
t861 = -m(5) * t775 + t786 * mrSges(5,1) - t787 * mrSges(5,2) + t822 * t795 - t823 * t796 - t864;
t718 = m(4) * t793 + qJDD(3) * mrSges(4,1) - t827 * mrSges(4,3) + qJD(3) * t830 - t824 * t881 + t861;
t692 = t851 * t698 + t856 * t718;
t807 = -qJDD(1) * pkin(1) + t870;
t869 = -m(3) * t807 + (t859 * mrSges(3,3)) - t692;
t688 = m(2) * t832 - (t859 * mrSges(2,2)) + qJDD(1) * t885 + t869;
t804 = t859 * pkin(1) + t887;
t700 = t855 * t706 + t850 * t707;
t868 = -m(4) * t800 + mrSges(4,1) * t826 - t827 * mrSges(4,2) - t830 * t838 - t831 * t881 - t700;
t865 = -m(3) * t804 + (t859 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t868;
t695 = m(2) * t833 - (mrSges(2,1) * t859) - qJDD(1) * mrSges(2,2) + t865;
t882 = t857 * t688 + t852 * t695;
t877 = -t688 * t852 + t857 * t695;
t876 = t856 * t698 - t851 * t718;
t748 = Ifges(7,4) * t769 + Ifges(7,2) * t768 + Ifges(7,6) * t829;
t749 = Ifges(7,1) * t769 + Ifges(7,4) * t768 + Ifges(7,5) * t829;
t867 = -mrSges(7,1) * t726 + mrSges(7,2) * t727 - Ifges(7,5) * t740 - Ifges(7,6) * t739 - Ifges(7,3) * t810 - t769 * t748 + t768 * t749;
t747 = Ifges(7,5) * t769 + Ifges(7,6) * t768 + Ifges(7,3) * t829;
t716 = -mrSges(7,1) * t734 + mrSges(7,3) * t727 + Ifges(7,4) * t740 + Ifges(7,2) * t739 + Ifges(7,6) * t810 - t747 * t769 + t749 * t829;
t717 = mrSges(7,2) * t734 - mrSges(7,3) * t726 + Ifges(7,1) * t740 + Ifges(7,4) * t739 + Ifges(7,5) * t810 + t747 * t768 - t748 * t829;
t764 = Ifges(6,5) * t790 + Ifges(6,6) * t789 + Ifges(6,3) * t834;
t766 = Ifges(6,1) * t790 + Ifges(6,4) * t789 + Ifges(6,5) * t834;
t701 = -mrSges(6,1) * t753 + mrSges(6,3) * t732 + Ifges(6,4) * t761 + Ifges(6,2) * t760 + Ifges(6,6) * t814 - pkin(5) * t872 + pkin(10) * t873 + t853 * t716 + t848 * t717 - t790 * t764 + t834 * t766;
t765 = Ifges(6,4) * t790 + Ifges(6,2) * t789 + Ifges(6,6) * t834;
t702 = mrSges(6,2) * t753 - mrSges(6,3) * t731 + Ifges(6,1) * t761 + Ifges(6,4) * t760 + Ifges(6,5) * t814 - pkin(10) * t715 - t716 * t848 + t717 * t853 + t764 * t789 - t765 * t834;
t780 = Ifges(5,5) * t823 + Ifges(5,6) * t822 + Ifges(5,3) * t835;
t782 = Ifges(5,1) * t823 + Ifges(5,4) * t822 + Ifges(5,5) * t835;
t684 = -mrSges(5,1) * t775 + mrSges(5,3) * t755 + Ifges(5,4) * t787 + Ifges(5,2) * t786 + Ifges(5,6) * t821 - pkin(4) * t864 + pkin(9) * t874 + t854 * t701 + t849 * t702 - t823 * t780 + t835 * t782;
t781 = Ifges(5,4) * t823 + Ifges(5,2) * t822 + Ifges(5,6) * t835;
t686 = mrSges(5,2) * t775 - mrSges(5,3) * t754 + Ifges(5,1) * t787 + Ifges(5,4) * t786 + Ifges(5,5) * t821 - pkin(9) * t708 - t701 * t849 + t702 * t854 + t780 * t822 - t781 * t835;
t812 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t856 - Ifges(4,2) * t851) * qJD(1);
t813 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t856 - Ifges(4,4) * t851) * qJD(1);
t866 = mrSges(4,1) * t793 - mrSges(4,2) * t794 + Ifges(4,5) * t827 + Ifges(4,6) * t826 + Ifges(4,3) * qJDD(3) + pkin(3) * t861 + pkin(8) * t875 + t855 * t684 + t850 * t686 + t812 * t881 + t813 * t838;
t811 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t856 - Ifges(4,6) * t851) * qJD(1);
t681 = mrSges(4,2) * t800 - mrSges(4,3) * t793 + Ifges(4,1) * t827 + Ifges(4,4) * t826 + Ifges(4,5) * qJDD(3) - pkin(8) * t700 - qJD(3) * t812 - t684 * t850 + t686 * t855 - t811 * t838;
t862 = -mrSges(6,1) * t731 + mrSges(6,2) * t732 - Ifges(6,5) * t761 - Ifges(6,6) * t760 - Ifges(6,3) * t814 - pkin(5) * t715 - t790 * t765 + t789 * t766 + t867;
t860 = mrSges(5,1) * t754 - mrSges(5,2) * t755 + Ifges(5,5) * t787 + Ifges(5,6) * t786 + Ifges(5,3) * t821 + pkin(4) * t708 + t823 * t781 - t822 * t782 - t862;
t682 = -mrSges(4,1) * t800 + mrSges(4,3) * t794 + Ifges(4,4) * t827 + Ifges(4,2) * t826 + Ifges(4,6) * qJDD(3) - pkin(3) * t700 + qJD(3) * t813 - t811 * t881 - t860;
t690 = qJDD(1) * mrSges(3,2) - t869;
t863 = mrSges(2,1) * t832 - mrSges(2,2) * t833 + mrSges(3,2) * t807 - mrSges(3,3) * t804 - pkin(1) * t690 - pkin(7) * t692 + qJ(2) * t865 + t856 * t681 - t682 * t851 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t691 = -m(3) * g(3) + t876;
t679 = (t883 * t859) - mrSges(2,3) * t832 + mrSges(3,1) * t807 + t884 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t866 - qJ(2) * t691 + pkin(2) * t692;
t678 = -mrSges(3,1) * t804 + mrSges(2,3) * t833 - pkin(1) * t691 - pkin(2) * t868 - pkin(7) * t876 + g(3) * t885 - qJDD(1) * t883 - t851 * t681 - t856 * t682 + t859 * t884;
t1 = [-m(1) * g(1) + t877; -m(1) * g(2) + t882; (-m(1) - m(2) - m(3)) * g(3) + t876; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t882 - t852 * t678 + t857 * t679; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t877 + t857 * t678 + t852 * t679; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t863; t863; t690; t866; t860; -t862; -t867;];
tauJB  = t1;
