% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPRRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:16:15
% EndTime: 2019-05-04 20:16:24
% DurationCPUTime: 8.46s
% Computational Cost: add. (140062->280), mult. (255551->351), div. (0->0), fcn. (184446->14), ass. (0->138)
t856 = Ifges(5,1) + Ifges(6,2);
t845 = Ifges(5,4) + Ifges(6,6);
t844 = Ifges(5,5) - Ifges(6,4);
t855 = Ifges(5,2) + Ifges(6,3);
t843 = Ifges(5,6) - Ifges(6,5);
t854 = Ifges(5,3) + Ifges(6,1);
t793 = sin(pkin(11));
t797 = cos(pkin(11));
t778 = -t797 * g(1) - t793 * g(2);
t792 = sin(pkin(12));
t796 = cos(pkin(12));
t777 = t793 * g(1) - t797 * g(2);
t791 = -g(3) + qJDD(1);
t795 = sin(pkin(6));
t799 = cos(pkin(6));
t821 = t777 * t799 + t791 * t795;
t738 = -t792 * t778 + t821 * t796;
t739 = t796 * t778 + t821 * t792;
t753 = -t795 * t777 + t799 * t791 + qJDD(2);
t805 = cos(qJ(3));
t798 = cos(pkin(7));
t802 = sin(qJ(3));
t837 = t798 * t802;
t794 = sin(pkin(7));
t838 = t794 * t802;
t731 = t738 * t837 + t805 * t739 + t753 * t838;
t807 = qJD(3) ^ 2;
t729 = -t807 * pkin(3) + qJDD(3) * pkin(9) + t731;
t801 = sin(qJ(4));
t726 = t801 * t729;
t733 = -t794 * t738 + t798 * t753;
t804 = cos(qJ(4));
t836 = t804 * t733;
t723 = -t726 + t836;
t772 = (mrSges(6,2) * t804 - mrSges(6,3) * t801) * qJD(3);
t773 = (-mrSges(5,1) * t804 + mrSges(5,2) * t801) * qJD(3);
t828 = qJD(3) * qJD(4);
t826 = t804 * t828;
t774 = t801 * qJDD(3) + t826;
t830 = qJD(3) * t804;
t780 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t830;
t781 = -mrSges(6,1) * t830 - qJD(4) * mrSges(6,3);
t771 = (-pkin(4) * t804 - qJ(5) * t801) * qJD(3);
t806 = qJD(4) ^ 2;
t829 = t801 * qJD(3);
t820 = -t806 * qJ(5) + t771 * t829 + qJDD(5) + t726;
t847 = pkin(10) * t807;
t848 = -pkin(4) - pkin(10);
t719 = t774 * pkin(5) + t848 * qJDD(4) + (-pkin(5) * t828 - t801 * t847 - t733) * t804 + t820;
t827 = t801 * t828;
t775 = t804 * qJDD(3) - t827;
t783 = pkin(5) * t829 - qJD(4) * pkin(10);
t790 = t804 ^ 2;
t730 = -t802 * t739 + (t738 * t798 + t753 * t794) * t805;
t813 = -qJDD(3) * pkin(3) - t730;
t849 = -2 * qJD(5);
t810 = pkin(4) * t827 + t829 * t849 + (-t774 - t826) * qJ(5) + t813;
t722 = -t783 * t829 + (-pkin(5) * t790 - pkin(9)) * t807 + t848 * t775 + t810;
t800 = sin(qJ(6));
t803 = cos(qJ(6));
t715 = t803 * t719 - t800 * t722;
t769 = -t800 * qJD(4) - t803 * t830;
t748 = t769 * qJD(6) + t803 * qJDD(4) - t800 * t775;
t770 = t803 * qJD(4) - t800 * t830;
t749 = -t769 * mrSges(7,1) + t770 * mrSges(7,2);
t785 = qJD(6) + t829;
t751 = -t785 * mrSges(7,2) + t769 * mrSges(7,3);
t768 = qJDD(6) + t774;
t712 = m(7) * t715 + t768 * mrSges(7,1) - t748 * mrSges(7,3) - t770 * t749 + t785 * t751;
t716 = t800 * t719 + t803 * t722;
t747 = -t770 * qJD(6) - t800 * qJDD(4) - t803 * t775;
t752 = t785 * mrSges(7,1) - t770 * mrSges(7,3);
t713 = m(7) * t716 - t768 * mrSges(7,2) + t747 * mrSges(7,3) + t769 * t749 - t785 * t752;
t703 = t803 * t712 + t800 * t713;
t721 = -qJDD(4) * pkin(4) + t820 - t836;
t815 = -m(6) * t721 - t774 * mrSges(6,1) - t703;
t700 = m(5) * t723 - t774 * mrSges(5,3) + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t780 - t781) * qJD(4) + (-t772 - t773) * t829 + t815;
t724 = t804 * t729 + t801 * t733;
t779 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t829;
t812 = -t806 * pkin(4) + qJDD(4) * qJ(5) + t771 * t830 + t724;
t720 = qJD(4) * t849 - t812;
t782 = mrSges(6,1) * t829 + qJD(4) * mrSges(6,2);
t718 = -t790 * t847 + t775 * pkin(5) + ((2 * qJD(5)) + t783) * qJD(4) + t812;
t816 = -m(7) * t718 + t747 * mrSges(7,1) - t748 * mrSges(7,2) + t769 * t751 - t770 * t752;
t811 = -m(6) * t720 + qJDD(4) * mrSges(6,3) + qJD(4) * t782 + t772 * t830 - t816;
t706 = t773 * t830 + m(5) * t724 - qJDD(4) * mrSges(5,2) - qJD(4) * t779 + (mrSges(5,3) + mrSges(6,1)) * t775 + t811;
t824 = -t801 * t700 + t804 * t706;
t692 = m(4) * t731 - t807 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t824;
t695 = t804 * t700 + t801 * t706;
t694 = m(4) * t733 + t695;
t846 = t807 * pkin(9);
t728 = t813 - t846;
t725 = -t775 * pkin(4) + t810 - t846;
t834 = -t800 * t712 + t803 * t713;
t819 = -m(6) * t725 - t775 * mrSges(6,2) + t782 * t829 - t834;
t809 = -m(5) * t728 + t780 * t830 + t775 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t774 + (-t779 * t801 - t781 * t804) * qJD(3) + t819;
t698 = m(4) * t730 + qJDD(3) * mrSges(4,1) - t807 * mrSges(4,2) + t809;
t839 = t698 * t805;
t681 = t692 * t837 - t794 * t694 + t798 * t839;
t677 = m(3) * t738 + t681;
t688 = t805 * t692 - t802 * t698;
t687 = m(3) * t739 + t688;
t853 = t677 * t796 + t687 * t792;
t701 = qJDD(4) * mrSges(6,2) + qJD(4) * t781 + t772 * t829 - t815;
t740 = Ifges(7,5) * t770 + Ifges(7,6) * t769 + Ifges(7,3) * t785;
t742 = Ifges(7,1) * t770 + Ifges(7,4) * t769 + Ifges(7,5) * t785;
t707 = -mrSges(7,1) * t718 + mrSges(7,3) * t716 + Ifges(7,4) * t748 + Ifges(7,2) * t747 + Ifges(7,6) * t768 - t770 * t740 + t785 * t742;
t741 = Ifges(7,4) * t770 + Ifges(7,2) * t769 + Ifges(7,6) * t785;
t708 = mrSges(7,2) * t718 - mrSges(7,3) * t715 + Ifges(7,1) * t748 + Ifges(7,4) * t747 + Ifges(7,5) * t768 + t769 * t740 - t785 * t741;
t831 = t844 * qJD(4) + (t801 * t856 + t845 * t804) * qJD(3);
t832 = t843 * qJD(4) + (t801 * t845 + t804 * t855) * qJD(3);
t852 = (t832 * t801 - t831 * t804) * qJD(3) + t854 * qJDD(4) + t844 * t774 + t843 * t775 + mrSges(5,1) * t723 - mrSges(5,2) * t724 + mrSges(6,2) * t721 - mrSges(6,3) * t720 - pkin(4) * t701 - pkin(10) * t703 + qJ(5) * (t775 * mrSges(6,1) + t811) - t800 * t707 + t803 * t708;
t680 = t692 * t838 + t798 * t694 + t794 * t839;
t679 = m(3) * t753 + t680;
t667 = -t795 * t679 + t799 * t853;
t665 = m(2) * t777 + t667;
t673 = -t792 * t677 + t796 * t687;
t672 = m(2) * t778 + t673;
t835 = t797 * t665 + t793 * t672;
t833 = t854 * qJD(4) + (t801 * t844 + t804 * t843) * qJD(3);
t666 = t799 * t679 + t795 * t853;
t825 = -t793 * t665 + t797 * t672;
t823 = m(2) * t791 + t666;
t702 = -t774 * mrSges(6,3) + t781 * t830 - t819;
t682 = -mrSges(5,1) * t728 - mrSges(6,1) * t720 + mrSges(6,2) * t725 + mrSges(5,3) * t724 - pkin(4) * t702 - pkin(5) * t816 - pkin(10) * t834 + t831 * qJD(4) + t843 * qJDD(4) - t803 * t707 - t800 * t708 + t845 * t774 + t775 * t855 - t833 * t829;
t814 = mrSges(7,1) * t715 - mrSges(7,2) * t716 + Ifges(7,5) * t748 + Ifges(7,6) * t747 + Ifges(7,3) * t768 + t770 * t741 - t769 * t742;
t683 = mrSges(6,1) * t721 + mrSges(5,2) * t728 - mrSges(5,3) * t723 - mrSges(6,3) * t725 + pkin(5) * t703 - qJ(5) * t702 - t832 * qJD(4) + t844 * qJDD(4) + t774 * t856 + t845 * t775 + t833 * t830 + t814;
t669 = mrSges(4,2) * t733 - mrSges(4,3) * t730 + Ifges(4,5) * qJDD(3) - t807 * Ifges(4,6) - pkin(9) * t695 - t801 * t682 + t804 * t683;
t674 = -mrSges(4,1) * t733 + mrSges(4,3) * t731 + t807 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t695 - t852;
t818 = pkin(8) * t688 + t669 * t802 + t674 * t805;
t668 = mrSges(4,1) * t730 - mrSges(4,2) * t731 + Ifges(4,3) * qJDD(3) + pkin(3) * t809 + pkin(9) * t824 + t804 * t682 + t801 * t683;
t662 = -mrSges(3,1) * t753 + mrSges(3,3) * t739 - pkin(2) * t680 - t794 * t668 + t818 * t798;
t663 = mrSges(3,2) * t753 - mrSges(3,3) * t738 + t805 * t669 - t802 * t674 + (-t680 * t794 - t681 * t798) * pkin(8);
t817 = qJ(2) * t673 + t662 * t796 + t663 * t792;
t661 = mrSges(3,1) * t738 - mrSges(3,2) * t739 + pkin(2) * t681 + t798 * t668 + t818 * t794;
t660 = mrSges(2,2) * t791 - mrSges(2,3) * t777 - t792 * t662 + t796 * t663 + (-t666 * t795 - t667 * t799) * qJ(2);
t659 = -mrSges(2,1) * t791 + mrSges(2,3) * t778 - pkin(1) * t666 - t795 * t661 + t817 * t799;
t1 = [-m(1) * g(1) + t825; -m(1) * g(2) + t835; -m(1) * g(3) + t823; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t835 - t793 * t659 + t797 * t660; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t825 + t797 * t659 + t793 * t660; -mrSges(1,1) * g(2) + mrSges(2,1) * t777 + mrSges(1,2) * g(1) - mrSges(2,2) * t778 + pkin(1) * t667 + t799 * t661 + t817 * t795; t823; t679; t668; t852; t701; t814;];
tauJB  = t1;
