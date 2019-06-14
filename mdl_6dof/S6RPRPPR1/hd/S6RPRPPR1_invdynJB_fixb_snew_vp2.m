% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-05-05 16:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:24:34
% EndTime: 2019-05-05 16:24:48
% DurationCPUTime: 13.93s
% Computational Cost: add. (211087->342), mult. (471760->433), div. (0->0), fcn. (317058->12), ass. (0->136)
t876 = -2 * qJD(4);
t847 = sin(qJ(1));
t850 = cos(qJ(1));
t829 = t847 * g(1) - g(2) * t850;
t820 = qJDD(1) * pkin(1) + t829;
t830 = -g(1) * t850 - g(2) * t847;
t852 = qJD(1) ^ 2;
t822 = -pkin(1) * t852 + t830;
t842 = sin(pkin(9));
t844 = cos(pkin(9));
t795 = t842 * t820 + t844 * t822;
t787 = -pkin(2) * t852 + qJDD(1) * pkin(7) + t795;
t839 = -g(3) + qJDD(2);
t846 = sin(qJ(3));
t849 = cos(qJ(3));
t775 = -t846 * t787 + t849 * t839;
t868 = qJD(1) * qJD(3);
t867 = t849 * t868;
t823 = qJDD(1) * t846 + t867;
t762 = (-t823 + t867) * qJ(4) + (t846 * t849 * t852 + qJDD(3)) * pkin(3) + t775;
t776 = t849 * t787 + t846 * t839;
t824 = qJDD(1) * t849 - t846 * t868;
t871 = qJD(1) * t846;
t826 = qJD(3) * pkin(3) - qJ(4) * t871;
t838 = t849 ^ 2;
t765 = -pkin(3) * t838 * t852 + qJ(4) * t824 - qJD(3) * t826 + t776;
t841 = sin(pkin(10));
t873 = cos(pkin(10));
t809 = (t841 * t849 + t846 * t873) * qJD(1);
t746 = t762 * t873 - t841 * t765 + t809 * t876;
t870 = qJD(1) * t849;
t808 = t841 * t871 - t873 * t870;
t747 = t841 * t762 + t873 * t765 + t808 * t876;
t789 = pkin(4) * t808 - qJ(5) * t809;
t851 = qJD(3) ^ 2;
t745 = -pkin(4) * t851 + qJDD(3) * qJ(5) - t789 * t808 + t747;
t794 = t844 * t820 - t842 * t822;
t858 = -qJDD(1) * pkin(2) - t794;
t766 = -t824 * pkin(3) + qJDD(4) + t826 * t871 + (-qJ(4) * t838 - pkin(7)) * t852 + t858;
t796 = t823 * t841 - t873 * t824;
t797 = t823 * t873 + t841 * t824;
t750 = (qJD(3) * t808 - t797) * qJ(5) + (qJD(3) * t809 + t796) * pkin(4) + t766;
t840 = sin(pkin(11));
t843 = cos(pkin(11));
t802 = qJD(3) * t840 + t809 * t843;
t740 = -0.2e1 * qJD(5) * t802 - t840 * t745 + t843 * t750;
t782 = qJDD(3) * t840 + t797 * t843;
t801 = qJD(3) * t843 - t809 * t840;
t738 = (t801 * t808 - t782) * pkin(8) + (t801 * t802 + t796) * pkin(5) + t740;
t741 = 0.2e1 * qJD(5) * t801 + t843 * t745 + t840 * t750;
t778 = pkin(5) * t808 - pkin(8) * t802;
t781 = qJDD(3) * t843 - t797 * t840;
t800 = t801 ^ 2;
t739 = -pkin(5) * t800 + pkin(8) * t781 - t778 * t808 + t741;
t845 = sin(qJ(6));
t848 = cos(qJ(6));
t736 = t738 * t848 - t739 * t845;
t771 = t801 * t848 - t802 * t845;
t753 = qJD(6) * t771 + t781 * t845 + t782 * t848;
t772 = t801 * t845 + t802 * t848;
t758 = -mrSges(7,1) * t771 + mrSges(7,2) * t772;
t807 = qJD(6) + t808;
t763 = -mrSges(7,2) * t807 + mrSges(7,3) * t771;
t793 = qJDD(6) + t796;
t733 = m(7) * t736 + mrSges(7,1) * t793 - mrSges(7,3) * t753 - t758 * t772 + t763 * t807;
t737 = t738 * t845 + t739 * t848;
t752 = -qJD(6) * t772 + t781 * t848 - t782 * t845;
t764 = mrSges(7,1) * t807 - mrSges(7,3) * t772;
t734 = m(7) * t737 - mrSges(7,2) * t793 + mrSges(7,3) * t752 + t758 * t771 - t764 * t807;
t725 = t848 * t733 + t845 * t734;
t773 = -mrSges(6,1) * t801 + mrSges(6,2) * t802;
t861 = -mrSges(6,2) * t808 + mrSges(6,3) * t801;
t723 = m(6) * t740 + t796 * mrSges(6,1) - t782 * mrSges(6,3) - t802 * t773 + t808 * t861 + t725;
t777 = mrSges(6,1) * t808 - mrSges(6,3) * t802;
t862 = -t733 * t845 + t848 * t734;
t724 = m(6) * t741 - mrSges(6,2) * t796 + mrSges(6,3) * t781 + t773 * t801 - t777 * t808 + t862;
t719 = -t723 * t840 + t843 * t724;
t790 = mrSges(5,1) * t808 + mrSges(5,2) * t809;
t804 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t809;
t716 = m(5) * t747 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t796 - qJD(3) * t804 - t790 * t808 + t719;
t744 = -qJDD(3) * pkin(4) - t851 * qJ(5) + t809 * t789 + qJDD(5) - t746;
t742 = -t781 * pkin(5) - t800 * pkin(8) + t802 * t778 + t744;
t857 = m(7) * t742 - t752 * mrSges(7,1) + mrSges(7,2) * t753 - t771 * t763 + t764 * t772;
t735 = m(6) * t744 - t781 * mrSges(6,1) + mrSges(6,2) * t782 + t777 * t802 - t801 * t861 + t857;
t803 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t808;
t729 = m(5) * t746 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t797 + qJD(3) * t803 - t790 * t809 - t735;
t708 = t841 * t716 + t873 * t729;
t754 = Ifges(7,5) * t772 + Ifges(7,6) * t771 + Ifges(7,3) * t807;
t756 = Ifges(7,1) * t772 + Ifges(7,4) * t771 + Ifges(7,5) * t807;
t726 = -mrSges(7,1) * t742 + mrSges(7,3) * t737 + Ifges(7,4) * t753 + Ifges(7,2) * t752 + Ifges(7,6) * t793 - t754 * t772 + t756 * t807;
t755 = Ifges(7,4) * t772 + Ifges(7,2) * t771 + Ifges(7,6) * t807;
t727 = mrSges(7,2) * t742 - mrSges(7,3) * t736 + Ifges(7,1) * t753 + Ifges(7,4) * t752 + Ifges(7,5) * t793 + t754 * t771 - t755 * t807;
t767 = Ifges(6,5) * t802 + Ifges(6,6) * t801 + Ifges(6,3) * t808;
t769 = Ifges(6,1) * t802 + Ifges(6,4) * t801 + Ifges(6,5) * t808;
t709 = -mrSges(6,1) * t744 + mrSges(6,3) * t741 + Ifges(6,4) * t782 + Ifges(6,2) * t781 + Ifges(6,6) * t796 - pkin(5) * t857 + pkin(8) * t862 + t848 * t726 + t845 * t727 - t802 * t767 + t808 * t769;
t768 = Ifges(6,4) * t802 + Ifges(6,2) * t801 + Ifges(6,6) * t808;
t710 = mrSges(6,2) * t744 - mrSges(6,3) * t740 + Ifges(6,1) * t782 + Ifges(6,4) * t781 + Ifges(6,5) * t796 - pkin(8) * t725 - t726 * t845 + t727 * t848 + t767 * t801 - t768 * t808;
t784 = Ifges(5,4) * t809 - Ifges(5,2) * t808 + Ifges(5,6) * qJD(3);
t785 = Ifges(5,1) * t809 - Ifges(5,4) * t808 + Ifges(5,5) * qJD(3);
t814 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t846 + Ifges(4,2) * t849) * qJD(1);
t815 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t846 + Ifges(4,4) * t849) * qJD(1);
t875 = (t814 * t846 - t815 * t849) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t775 + mrSges(5,1) * t746 - mrSges(4,2) * t776 - mrSges(5,2) * t747 + Ifges(4,5) * t823 + Ifges(5,5) * t797 + Ifges(4,6) * t824 - Ifges(5,6) * t796 + pkin(3) * t708 - pkin(4) * t735 + qJ(5) * t719 + t843 * t709 + t840 * t710 + t809 * t784 + t808 * t785;
t821 = (-mrSges(4,1) * t849 + mrSges(4,2) * t846) * qJD(1);
t828 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t870;
t706 = m(4) * t775 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t823 + qJD(3) * t828 - t821 * t871 + t708;
t827 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t871;
t863 = t873 * t716 - t729 * t841;
t707 = m(4) * t776 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t824 - qJD(3) * t827 + t821 * t870 + t863;
t864 = -t706 * t846 + t849 * t707;
t698 = m(3) * t795 - mrSges(3,1) * t852 - qJDD(1) * mrSges(3,2) + t864;
t718 = t843 * t723 + t840 * t724;
t717 = m(5) * t766 + t796 * mrSges(5,1) + mrSges(5,2) * t797 + t808 * t803 + t804 * t809 + t718;
t786 = -t852 * pkin(7) + t858;
t854 = -m(4) * t786 + t824 * mrSges(4,1) - mrSges(4,2) * t823 - t827 * t871 + t828 * t870 - t717;
t712 = m(3) * t794 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t852 + t854;
t694 = t842 * t698 + t844 * t712;
t691 = m(2) * t829 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t852 + t694;
t865 = t844 * t698 - t712 * t842;
t692 = m(2) * t830 - mrSges(2,1) * t852 - qJDD(1) * mrSges(2,2) + t865;
t872 = t850 * t691 + t847 * t692;
t701 = t849 * t706 + t846 * t707;
t699 = m(3) * t839 + t701;
t866 = -t691 * t847 + t850 * t692;
t783 = Ifges(5,5) * t809 - Ifges(5,6) * t808 + Ifges(5,3) * qJD(3);
t695 = mrSges(5,2) * t766 - mrSges(5,3) * t746 + Ifges(5,1) * t797 - Ifges(5,4) * t796 + Ifges(5,5) * qJDD(3) - qJ(5) * t718 - qJD(3) * t784 - t709 * t840 + t710 * t843 - t783 * t808;
t855 = mrSges(7,1) * t736 - mrSges(7,2) * t737 + Ifges(7,5) * t753 + Ifges(7,6) * t752 + Ifges(7,3) * t793 + t772 * t755 - t771 * t756;
t702 = -t809 * t783 + t801 * t769 - t802 * t768 + Ifges(5,4) * t797 - Ifges(6,6) * t781 - Ifges(6,5) * t782 + qJD(3) * t785 - mrSges(5,1) * t766 + mrSges(5,3) * t747 - mrSges(6,1) * t740 + mrSges(6,2) * t741 - pkin(5) * t725 - pkin(4) * t718 + Ifges(5,6) * qJDD(3) + (-Ifges(5,2) - Ifges(6,3)) * t796 - t855;
t813 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t846 + Ifges(4,6) * t849) * qJD(1);
t684 = -mrSges(4,1) * t786 + mrSges(4,3) * t776 + Ifges(4,4) * t823 + Ifges(4,2) * t824 + Ifges(4,6) * qJDD(3) - pkin(3) * t717 + qJ(4) * t863 + qJD(3) * t815 + t841 * t695 + t702 * t873 - t813 * t871;
t687 = mrSges(4,2) * t786 - mrSges(4,3) * t775 + Ifges(4,1) * t823 + Ifges(4,4) * t824 + Ifges(4,5) * qJDD(3) - qJ(4) * t708 - qJD(3) * t814 + t695 * t873 - t841 * t702 + t813 * t870;
t856 = mrSges(2,1) * t829 + mrSges(3,1) * t794 - mrSges(2,2) * t830 - mrSges(3,2) * t795 + pkin(1) * t694 + pkin(2) * t854 + pkin(7) * t864 + t849 * t684 + t846 * t687 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t685 = -mrSges(3,1) * t839 + mrSges(3,3) * t795 + t852 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t701 - t875;
t682 = mrSges(3,2) * t839 - mrSges(3,3) * t794 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t852 - pkin(7) * t701 - t684 * t846 + t687 * t849;
t681 = -mrSges(2,2) * g(3) - mrSges(2,3) * t829 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t852 - qJ(2) * t694 + t682 * t844 - t685 * t842;
t680 = mrSges(2,1) * g(3) + mrSges(2,3) * t830 + t852 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t699 + qJ(2) * t865 + t842 * t682 + t844 * t685;
t1 = [-m(1) * g(1) + t866; -m(1) * g(2) + t872; (-m(1) - m(2)) * g(3) + t699; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t872 - t847 * t680 + t850 * t681; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t866 + t850 * t680 + t847 * t681; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t856; t856; t699; t875; t717; t735; t855;];
tauJB  = t1;
