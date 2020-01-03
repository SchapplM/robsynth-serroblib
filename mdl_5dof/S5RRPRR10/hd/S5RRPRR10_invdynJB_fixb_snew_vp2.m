% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:24:50
% EndTime: 2019-12-31 20:25:04
% DurationCPUTime: 14.00s
% Computational Cost: add. (196979->325), mult. (516522->428), div. (0->0), fcn. (396194->12), ass. (0->137)
t822 = -2 * qJD(3);
t786 = sin(pkin(10));
t788 = cos(pkin(10));
t792 = sin(qJ(2));
t796 = cos(qJ(2));
t787 = sin(pkin(5));
t814 = qJD(1) * t787;
t765 = (t786 * t792 - t788 * t796) * t814;
t812 = qJD(1) * qJD(2);
t774 = (qJDD(1) * t792 + t796 * t812) * t787;
t789 = cos(pkin(5));
t780 = t789 * qJDD(1) + qJDD(2);
t781 = t789 * qJD(1) + qJD(2);
t793 = sin(qJ(1));
t797 = cos(qJ(1));
t777 = t793 * g(1) - t797 * g(2);
t798 = qJD(1) ^ 2;
t821 = pkin(7) * t787;
t771 = qJDD(1) * pkin(1) + t798 * t821 + t777;
t778 = -t797 * g(1) - t793 * g(2);
t772 = -t798 * pkin(1) + qJDD(1) * t821 + t778;
t816 = t789 * t796;
t804 = t771 * t816 - t792 * t772;
t820 = t787 ^ 2 * t798;
t716 = t780 * pkin(2) - t774 * qJ(3) + (pkin(2) * t792 * t820 + (qJ(3) * qJD(1) * t781 - g(3)) * t787) * t796 + t804;
t817 = t789 * t792;
t819 = t787 * t792;
t742 = -g(3) * t819 + t771 * t817 + t796 * t772;
t810 = t792 * t814;
t768 = t781 * pkin(2) - qJ(3) * t810;
t775 = (qJDD(1) * t796 - t792 * t812) * t787;
t811 = t796 ^ 2 * t820;
t719 = -pkin(2) * t811 + t775 * qJ(3) - t781 * t768 + t742;
t766 = (t786 * t796 + t788 * t792) * t814;
t701 = t788 * t716 - t786 * t719 + t766 * t822;
t818 = t787 * t796;
t702 = t786 * t716 + t788 * t719 + t765 * t822;
t743 = t765 * mrSges(4,1) + t766 * mrSges(4,2);
t747 = -t786 * t774 + t788 * t775;
t753 = t781 * mrSges(4,1) - t766 * mrSges(4,3);
t744 = t765 * pkin(3) - t766 * pkin(8);
t779 = t781 ^ 2;
t700 = -t779 * pkin(3) + t780 * pkin(8) - t765 * t744 + t702;
t757 = -t789 * g(3) - t787 * t771;
t729 = -t775 * pkin(2) - qJ(3) * t811 + t768 * t810 + qJDD(3) + t757;
t748 = t788 * t774 + t786 * t775;
t704 = (t765 * t781 - t748) * pkin(8) + (t766 * t781 - t747) * pkin(3) + t729;
t791 = sin(qJ(4));
t795 = cos(qJ(4));
t697 = t795 * t700 + t791 * t704;
t750 = -t791 * t766 + t795 * t781;
t751 = t795 * t766 + t791 * t781;
t731 = -t750 * pkin(4) - t751 * pkin(9);
t746 = qJDD(4) - t747;
t764 = qJD(4) + t765;
t763 = t764 ^ 2;
t694 = -t763 * pkin(4) + t746 * pkin(9) + t750 * t731 + t697;
t699 = -t780 * pkin(3) - t779 * pkin(8) + t766 * t744 - t701;
t726 = -t751 * qJD(4) - t791 * t748 + t795 * t780;
t727 = t750 * qJD(4) + t795 * t748 + t791 * t780;
t695 = (-t750 * t764 - t727) * pkin(9) + (t751 * t764 - t726) * pkin(4) + t699;
t790 = sin(qJ(5));
t794 = cos(qJ(5));
t691 = -t790 * t694 + t794 * t695;
t733 = -t790 * t751 + t794 * t764;
t707 = t733 * qJD(5) + t794 * t727 + t790 * t746;
t734 = t794 * t751 + t790 * t764;
t712 = -t733 * mrSges(6,1) + t734 * mrSges(6,2);
t749 = qJD(5) - t750;
t717 = -t749 * mrSges(6,2) + t733 * mrSges(6,3);
t725 = qJDD(5) - t726;
t688 = m(6) * t691 + t725 * mrSges(6,1) - t707 * mrSges(6,3) - t734 * t712 + t749 * t717;
t692 = t794 * t694 + t790 * t695;
t706 = -t734 * qJD(5) - t790 * t727 + t794 * t746;
t718 = t749 * mrSges(6,1) - t734 * mrSges(6,3);
t689 = m(6) * t692 - t725 * mrSges(6,2) + t706 * mrSges(6,3) + t733 * t712 - t749 * t718;
t682 = -t790 * t688 + t794 * t689;
t730 = -t750 * mrSges(5,1) + t751 * mrSges(5,2);
t736 = t764 * mrSges(5,1) - t751 * mrSges(5,3);
t680 = m(5) * t697 - t746 * mrSges(5,2) + t726 * mrSges(5,3) + t750 * t730 - t764 * t736 + t682;
t696 = -t791 * t700 + t795 * t704;
t693 = -t746 * pkin(4) - t763 * pkin(9) + t751 * t731 - t696;
t690 = -m(6) * t693 + t706 * mrSges(6,1) - t707 * mrSges(6,2) + t733 * t717 - t734 * t718;
t735 = -t764 * mrSges(5,2) + t750 * mrSges(5,3);
t686 = m(5) * t696 + t746 * mrSges(5,1) - t727 * mrSges(5,3) - t751 * t730 + t764 * t735 + t690;
t806 = t795 * t680 - t791 * t686;
t671 = m(4) * t702 - t780 * mrSges(4,2) + t747 * mrSges(4,3) - t765 * t743 - t781 * t753 + t806;
t752 = -t781 * mrSges(4,2) - t765 * mrSges(4,3);
t681 = t794 * t688 + t790 * t689;
t801 = -m(5) * t699 + t726 * mrSges(5,1) - t727 * mrSges(5,2) + t750 * t735 - t751 * t736 - t681;
t677 = m(4) * t701 + t780 * mrSges(4,1) - t748 * mrSges(4,3) - t766 * t743 + t781 * t752 + t801;
t666 = t786 * t671 + t788 * t677;
t741 = -g(3) * t818 + t804;
t809 = t796 * t814;
t770 = -t781 * mrSges(3,2) + mrSges(3,3) * t809;
t773 = (-mrSges(3,1) * t796 + mrSges(3,2) * t792) * t814;
t664 = m(3) * t741 + t780 * mrSges(3,1) - t774 * mrSges(3,3) + t781 * t770 - t773 * t810 + t666;
t769 = t781 * mrSges(3,1) - mrSges(3,3) * t810;
t807 = t788 * t671 - t786 * t677;
t665 = m(3) * t742 - t780 * mrSges(3,2) + t775 * mrSges(3,3) - t781 * t769 + t773 * t809 + t807;
t675 = t791 * t680 + t795 * t686;
t674 = m(4) * t729 - t747 * mrSges(4,1) + t748 * mrSges(4,2) + t765 * t752 + t766 * t753 + t675;
t673 = m(3) * t757 - t775 * mrSges(3,1) + t774 * mrSges(3,2) + (t769 * t792 - t770 * t796) * t814 + t674;
t651 = t664 * t816 + t665 * t817 - t787 * t673;
t648 = m(2) * t777 + qJDD(1) * mrSges(2,1) - t798 * mrSges(2,2) + t651;
t657 = -t792 * t664 + t796 * t665;
t655 = m(2) * t778 - t798 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t657;
t815 = t797 * t648 + t793 * t655;
t650 = t664 * t818 + t665 * t819 + t789 * t673;
t808 = -t793 * t648 + t797 * t655;
t708 = Ifges(6,5) * t734 + Ifges(6,6) * t733 + Ifges(6,3) * t749;
t710 = Ifges(6,1) * t734 + Ifges(6,4) * t733 + Ifges(6,5) * t749;
t683 = -mrSges(6,1) * t693 + mrSges(6,3) * t692 + Ifges(6,4) * t707 + Ifges(6,2) * t706 + Ifges(6,6) * t725 - t734 * t708 + t749 * t710;
t709 = Ifges(6,4) * t734 + Ifges(6,2) * t733 + Ifges(6,6) * t749;
t684 = mrSges(6,2) * t693 - mrSges(6,3) * t691 + Ifges(6,1) * t707 + Ifges(6,4) * t706 + Ifges(6,5) * t725 + t733 * t708 - t749 * t709;
t720 = Ifges(5,5) * t751 + Ifges(5,6) * t750 + Ifges(5,3) * t764;
t721 = Ifges(5,4) * t751 + Ifges(5,2) * t750 + Ifges(5,6) * t764;
t667 = mrSges(5,2) * t699 - mrSges(5,3) * t696 + Ifges(5,1) * t727 + Ifges(5,4) * t726 + Ifges(5,5) * t746 - pkin(9) * t681 - t790 * t683 + t794 * t684 + t750 * t720 - t764 * t721;
t722 = Ifges(5,1) * t751 + Ifges(5,4) * t750 + Ifges(5,5) * t764;
t800 = mrSges(6,1) * t691 - mrSges(6,2) * t692 + Ifges(6,5) * t707 + Ifges(6,6) * t706 + Ifges(6,3) * t725 + t734 * t709 - t733 * t710;
t668 = -mrSges(5,1) * t699 + mrSges(5,3) * t697 + Ifges(5,4) * t727 + Ifges(5,2) * t726 + Ifges(5,6) * t746 - pkin(4) * t681 - t751 * t720 + t764 * t722 - t800;
t737 = Ifges(4,5) * t766 - Ifges(4,6) * t765 + Ifges(4,3) * t781;
t738 = Ifges(4,4) * t766 - Ifges(4,2) * t765 + Ifges(4,6) * t781;
t652 = mrSges(4,2) * t729 - mrSges(4,3) * t701 + Ifges(4,1) * t748 + Ifges(4,4) * t747 + Ifges(4,5) * t780 - pkin(8) * t675 + t795 * t667 - t791 * t668 - t765 * t737 - t781 * t738;
t739 = Ifges(4,1) * t766 - Ifges(4,4) * t765 + Ifges(4,5) * t781;
t799 = mrSges(5,1) * t696 - mrSges(5,2) * t697 + Ifges(5,5) * t727 + Ifges(5,6) * t726 + Ifges(5,3) * t746 + pkin(4) * t690 + pkin(9) * t682 + t794 * t683 + t790 * t684 + t751 * t721 - t750 * t722;
t658 = -mrSges(4,1) * t729 + mrSges(4,3) * t702 + Ifges(4,4) * t748 + Ifges(4,2) * t747 + Ifges(4,6) * t780 - pkin(3) * t675 - t766 * t737 + t781 * t739 - t799;
t754 = Ifges(3,3) * t781 + (Ifges(3,5) * t792 + Ifges(3,6) * t796) * t814;
t756 = Ifges(3,5) * t781 + (Ifges(3,1) * t792 + Ifges(3,4) * t796) * t814;
t642 = -mrSges(3,1) * t757 + mrSges(3,3) * t742 + Ifges(3,4) * t774 + Ifges(3,2) * t775 + Ifges(3,6) * t780 - pkin(2) * t674 + qJ(3) * t807 + t786 * t652 + t788 * t658 - t754 * t810 + t781 * t756;
t755 = Ifges(3,6) * t781 + (Ifges(3,4) * t792 + Ifges(3,2) * t796) * t814;
t644 = mrSges(3,2) * t757 - mrSges(3,3) * t741 + Ifges(3,1) * t774 + Ifges(3,4) * t775 + Ifges(3,5) * t780 - qJ(3) * t666 + t788 * t652 - t786 * t658 + t754 * t809 - t781 * t755;
t646 = Ifges(3,5) * t774 + Ifges(3,6) * t775 + mrSges(3,1) * t741 - mrSges(3,2) * t742 + Ifges(4,5) * t748 + Ifges(4,6) * t747 + t766 * t738 + t765 * t739 + mrSges(4,1) * t701 - mrSges(4,2) * t702 + t791 * t667 + t795 * t668 + pkin(3) * t801 + pkin(8) * t806 + pkin(2) * t666 + (Ifges(3,3) + Ifges(4,3)) * t780 + (t755 * t792 - t756 * t796) * t814;
t802 = mrSges(2,1) * t777 - mrSges(2,2) * t778 + Ifges(2,3) * qJDD(1) + pkin(1) * t651 + t642 * t818 + t644 * t819 + t789 * t646 + t657 * t821;
t640 = -mrSges(2,2) * g(3) - mrSges(2,3) * t777 + Ifges(2,5) * qJDD(1) - t798 * Ifges(2,6) - t792 * t642 + t796 * t644 + (-t650 * t787 - t651 * t789) * pkin(7);
t639 = mrSges(2,1) * g(3) + mrSges(2,3) * t778 + t798 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t650 - t787 * t646 + (pkin(7) * t657 + t642 * t796 + t644 * t792) * t789;
t1 = [-m(1) * g(1) + t808; -m(1) * g(2) + t815; (-m(1) - m(2)) * g(3) + t650; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t815 - t793 * t639 + t797 * t640; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t808 + t797 * t639 + t793 * t640; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t802; t802; t646; t674; t799; t800;];
tauJB = t1;
