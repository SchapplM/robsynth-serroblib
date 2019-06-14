% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-05 17:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:08:00
% EndTime: 2019-05-05 17:08:09
% DurationCPUTime: 8.34s
% Computational Cost: add. (124231->337), mult. (281235->416), div. (0->0), fcn. (186856->10), ass. (0->133)
t883 = -2 * qJD(4);
t850 = sin(qJ(1));
t853 = cos(qJ(1));
t829 = t850 * g(1) - t853 * g(2);
t855 = qJD(1) ^ 2;
t863 = -t855 * qJ(2) + qJDD(2) - t829;
t882 = -pkin(1) - pkin(7);
t801 = qJDD(1) * t882 + t863;
t849 = sin(qJ(3));
t852 = cos(qJ(3));
t788 = t849 * g(3) + t852 * t801;
t874 = qJD(1) * qJD(3);
t872 = t849 * t874;
t824 = qJDD(1) * t852 - t872;
t760 = (-t824 - t872) * qJ(4) + (-t849 * t852 * t855 + qJDD(3)) * pkin(3) + t788;
t789 = -g(3) * t852 + t849 * t801;
t823 = -qJDD(1) * t849 - t852 * t874;
t876 = qJD(1) * t852;
t827 = qJD(3) * pkin(3) - qJ(4) * t876;
t841 = t849 ^ 2;
t761 = -pkin(3) * t841 * t855 + qJ(4) * t823 - qJD(3) * t827 + t789;
t845 = sin(pkin(9));
t847 = cos(pkin(9));
t877 = qJD(1) * t849;
t811 = -t845 * t877 + t847 * t876;
t744 = t760 * t847 - t845 * t761 + t811 * t883;
t830 = -t853 * g(1) - t850 * g(2);
t864 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t830;
t881 = mrSges(2,1) - mrSges(3,2);
t880 = Ifges(2,5) - Ifges(3,4);
t879 = -Ifges(2,6) + Ifges(3,5);
t810 = (t845 * t852 + t847 * t849) * qJD(1);
t745 = t845 * t760 + t847 * t761 + t810 * t883;
t780 = mrSges(5,1) * t810 + mrSges(5,2) * t811;
t785 = -t847 * t823 + t824 * t845;
t800 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t811;
t779 = pkin(4) * t810 - qJ(5) * t811;
t854 = qJD(3) ^ 2;
t736 = -pkin(4) * t854 + qJDD(3) * qJ(5) - t779 * t810 + t745;
t765 = -t823 * pkin(3) + qJDD(4) + t827 * t876 + (-qJ(4) * t841 + t882) * t855 + t864;
t786 = t823 * t845 + t824 * t847;
t739 = (qJD(3) * t810 - t786) * qJ(5) + (qJD(3) * t811 + t785) * pkin(4) + t765;
t844 = sin(pkin(10));
t846 = cos(pkin(10));
t794 = qJD(3) * t844 + t811 * t846;
t731 = -0.2e1 * qJD(5) * t794 - t844 * t736 + t846 * t739;
t774 = qJDD(3) * t844 + t786 * t846;
t793 = qJD(3) * t846 - t811 * t844;
t729 = (t793 * t810 - t774) * pkin(8) + (t793 * t794 + t785) * pkin(5) + t731;
t732 = 0.2e1 * qJD(5) * t793 + t846 * t736 + t844 * t739;
t771 = pkin(5) * t810 - pkin(8) * t794;
t773 = qJDD(3) * t846 - t786 * t844;
t792 = t793 ^ 2;
t730 = -pkin(5) * t792 + pkin(8) * t773 - t771 * t810 + t732;
t848 = sin(qJ(6));
t851 = cos(qJ(6));
t727 = t729 * t851 - t730 * t848;
t763 = t793 * t851 - t794 * t848;
t742 = t763 * qJD(6) + t773 * t848 + t774 * t851;
t764 = t793 * t848 + t794 * t851;
t750 = -mrSges(7,1) * t763 + mrSges(7,2) * t764;
t808 = qJD(6) + t810;
t751 = -mrSges(7,2) * t808 + t763 * mrSges(7,3);
t784 = qJDD(6) + t785;
t723 = m(7) * t727 + mrSges(7,1) * t784 - t742 * mrSges(7,3) - t764 * t750 + t751 * t808;
t728 = t729 * t848 + t730 * t851;
t741 = -t764 * qJD(6) + t773 * t851 - t774 * t848;
t752 = mrSges(7,1) * t808 - t764 * mrSges(7,3);
t724 = m(7) * t728 - mrSges(7,2) * t784 + t741 * mrSges(7,3) + t763 * t750 - t752 * t808;
t715 = t851 * t723 + t848 * t724;
t767 = -mrSges(6,1) * t793 + mrSges(6,2) * t794;
t769 = -mrSges(6,2) * t810 + mrSges(6,3) * t793;
t713 = m(6) * t731 + mrSges(6,1) * t785 - mrSges(6,3) * t774 - t767 * t794 + t769 * t810 + t715;
t770 = mrSges(6,1) * t810 - mrSges(6,3) * t794;
t867 = -t723 * t848 + t851 * t724;
t714 = m(6) * t732 - mrSges(6,2) * t785 + mrSges(6,3) * t773 + t767 * t793 - t770 * t810 + t867;
t868 = -t713 * t844 + t846 * t714;
t706 = m(5) * t745 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t785 - qJD(3) * t800 - t780 * t810 + t868;
t735 = -qJDD(3) * pkin(4) - qJ(5) * t854 + t811 * t779 + qJDD(5) - t744;
t733 = -pkin(5) * t773 - pkin(8) * t792 + t771 * t794 + t735;
t861 = m(7) * t733 - t741 * mrSges(7,1) + t742 * mrSges(7,2) - t763 * t751 + t764 * t752;
t726 = m(6) * t735 - t773 * mrSges(6,1) + mrSges(6,2) * t774 - t793 * t769 + t770 * t794 + t861;
t799 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t810;
t719 = m(5) * t744 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t786 + qJD(3) * t799 - t780 * t811 - t726;
t695 = t845 * t706 + t847 * t719;
t822 = (mrSges(4,1) * t849 + mrSges(4,2) * t852) * qJD(1);
t826 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t877;
t692 = m(4) * t788 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t824 + qJD(3) * t826 - t822 * t876 + t695;
t828 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t876;
t869 = t847 * t706 - t719 * t845;
t693 = m(4) * t789 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t823 - qJD(3) * t828 - t822 * t877 + t869;
t689 = t852 * t692 + t849 * t693;
t809 = -qJDD(1) * pkin(1) + t863;
t862 = -m(3) * t809 + t855 * mrSges(3,3) - t689;
t685 = m(2) * t829 - t855 * mrSges(2,2) + qJDD(1) * t881 + t862;
t804 = t855 * pkin(1) - t864;
t709 = t846 * t713 + t844 * t714;
t707 = m(5) * t765 + mrSges(5,1) * t785 + t786 * mrSges(5,2) + t799 * t810 + t811 * t800 + t709;
t798 = t855 * t882 + t864;
t860 = -m(4) * t798 + mrSges(4,1) * t823 - t824 * mrSges(4,2) - t826 * t877 - t828 * t876 - t707;
t857 = -m(3) * t804 + t855 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t860;
t702 = m(2) * t830 - mrSges(2,1) * t855 - qJDD(1) * mrSges(2,2) + t857;
t878 = t853 * t685 + t850 * t702;
t871 = -t685 * t850 + t853 * t702;
t870 = -t849 * t692 + t852 * t693;
t746 = Ifges(7,5) * t764 + Ifges(7,6) * t763 + Ifges(7,3) * t808;
t748 = Ifges(7,1) * t764 + Ifges(7,4) * t763 + Ifges(7,5) * t808;
t716 = -mrSges(7,1) * t733 + mrSges(7,3) * t728 + Ifges(7,4) * t742 + Ifges(7,2) * t741 + Ifges(7,6) * t784 - t764 * t746 + t748 * t808;
t747 = Ifges(7,4) * t764 + Ifges(7,2) * t763 + Ifges(7,6) * t808;
t717 = mrSges(7,2) * t733 - mrSges(7,3) * t727 + Ifges(7,1) * t742 + Ifges(7,4) * t741 + Ifges(7,5) * t784 + t763 * t746 - t747 * t808;
t753 = Ifges(6,5) * t794 + Ifges(6,6) * t793 + Ifges(6,3) * t810;
t755 = Ifges(6,1) * t794 + Ifges(6,4) * t793 + Ifges(6,5) * t810;
t697 = -mrSges(6,1) * t735 + mrSges(6,3) * t732 + Ifges(6,4) * t774 + Ifges(6,2) * t773 + Ifges(6,6) * t785 - pkin(5) * t861 + pkin(8) * t867 + t851 * t716 + t848 * t717 - t794 * t753 + t810 * t755;
t754 = Ifges(6,4) * t794 + Ifges(6,2) * t793 + Ifges(6,6) * t810;
t699 = mrSges(6,2) * t735 - mrSges(6,3) * t731 + Ifges(6,1) * t774 + Ifges(6,4) * t773 + Ifges(6,5) * t785 - pkin(8) * t715 - t716 * t848 + t717 * t851 + t753 * t793 - t754 * t810;
t775 = Ifges(5,5) * t811 - Ifges(5,6) * t810 + Ifges(5,3) * qJD(3);
t776 = Ifges(5,4) * t811 - Ifges(5,2) * t810 + Ifges(5,6) * qJD(3);
t683 = mrSges(5,2) * t765 - mrSges(5,3) * t744 + Ifges(5,1) * t786 - Ifges(5,4) * t785 + Ifges(5,5) * qJDD(3) - qJ(5) * t709 - qJD(3) * t776 - t697 * t844 + t699 * t846 - t775 * t810;
t777 = Ifges(5,1) * t811 - Ifges(5,4) * t810 + Ifges(5,5) * qJD(3);
t858 = mrSges(7,1) * t727 - mrSges(7,2) * t728 + Ifges(7,5) * t742 + Ifges(7,6) * t741 + Ifges(7,3) * t784 + t764 * t747 - t763 * t748;
t690 = -pkin(4) * t709 - pkin(5) * t715 - t858 + Ifges(5,6) * qJDD(3) + (-Ifges(5,2) - Ifges(6,3)) * t785 - mrSges(6,1) * t731 + mrSges(6,2) * t732 + mrSges(5,3) * t745 - mrSges(5,1) * t765 - Ifges(6,6) * t773 - Ifges(6,5) * t774 + qJD(3) * t777 + Ifges(5,4) * t786 + t793 * t755 - t794 * t754 - t811 * t775;
t812 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t852 - Ifges(4,6) * t849) * qJD(1);
t814 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t852 - Ifges(4,4) * t849) * qJD(1);
t680 = -mrSges(4,1) * t798 + mrSges(4,3) * t789 + Ifges(4,4) * t824 + Ifges(4,2) * t823 + Ifges(4,6) * qJDD(3) - pkin(3) * t707 + qJ(4) * t869 + qJD(3) * t814 + t845 * t683 + t847 * t690 - t812 * t876;
t813 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t852 - Ifges(4,2) * t849) * qJD(1);
t682 = mrSges(4,2) * t798 - mrSges(4,3) * t788 + Ifges(4,1) * t824 + Ifges(4,4) * t823 + Ifges(4,5) * qJDD(3) - qJ(4) * t695 - qJD(3) * t813 + t683 * t847 - t690 * t845 - t812 * t877;
t687 = qJDD(1) * mrSges(3,2) - t862;
t859 = mrSges(2,1) * t829 - mrSges(2,2) * t830 + mrSges(3,2) * t809 - mrSges(3,3) * t804 - pkin(1) * t687 - pkin(7) * t689 + qJ(2) * t857 - t680 * t849 + t852 * t682 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t856 = mrSges(5,1) * t744 - mrSges(4,2) * t789 - mrSges(5,2) * t745 + Ifges(5,5) * t786 - Ifges(5,6) * t785 + pkin(3) * t695 - pkin(4) * t726 + qJ(5) * t868 + t846 * t697 + t844 * t699 + t811 * t776 + t810 * t777 + mrSges(4,1) * t788 + t814 * t877 + t813 * t876 + Ifges(4,6) * t823 + Ifges(4,5) * t824 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t688 = -m(3) * g(3) + t870;
t679 = t856 + pkin(2) * t689 - qJ(2) * t688 + t879 * t855 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + mrSges(3,1) * t809 + t880 * qJDD(1) - mrSges(2,3) * t829;
t678 = -mrSges(3,1) * t804 + mrSges(2,3) * t830 - pkin(1) * t688 - pkin(2) * t860 - pkin(7) * t870 + g(3) * t881 - qJDD(1) * t879 - t852 * t680 - t849 * t682 + t855 * t880;
t1 = [-m(1) * g(1) + t871; -m(1) * g(2) + t878; (-m(1) - m(2) - m(3)) * g(3) + t870; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t878 - t850 * t678 + t853 * t679; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t871 + t853 * t678 + t850 * t679; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t859; t859; t687; t856; t707; t726; t858;];
tauJB  = t1;
