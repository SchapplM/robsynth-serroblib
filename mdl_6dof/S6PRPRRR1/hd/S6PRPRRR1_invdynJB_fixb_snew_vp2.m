% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:05:49
% EndTime: 2019-05-05 00:06:01
% DurationCPUTime: 11.62s
% Computational Cost: add. (209097->300), mult. (393492->387), div. (0->0), fcn. (280826->14), ass. (0->134)
t799 = sin(qJ(5));
t800 = sin(qJ(4));
t803 = cos(qJ(5));
t804 = cos(qJ(4));
t768 = (t799 * t800 - t803 * t804) * qJD(2);
t793 = sin(pkin(11));
t796 = cos(pkin(11));
t777 = g(1) * t793 - g(2) * t796;
t778 = -g(1) * t796 - g(2) * t793;
t791 = -g(3) + qJDD(1);
t801 = sin(qJ(2));
t797 = cos(pkin(6));
t805 = cos(qJ(2));
t828 = t797 * t805;
t794 = sin(pkin(6));
t830 = t794 * t805;
t748 = t777 * t828 - t778 * t801 + t791 * t830;
t743 = qJDD(2) * pkin(2) + t748;
t829 = t797 * t801;
t831 = t794 * t801;
t749 = t777 * t829 + t805 * t778 + t791 * t831;
t806 = qJD(2) ^ 2;
t744 = -pkin(2) * t806 + t749;
t792 = sin(pkin(12));
t795 = cos(pkin(12));
t722 = t792 * t743 + t795 * t744;
t720 = -pkin(3) * t806 + qJDD(2) * pkin(8) + t722;
t761 = -t777 * t794 + t797 * t791;
t760 = qJDD(3) + t761;
t716 = -t800 * t720 + t804 * t760;
t824 = qJD(2) * qJD(4);
t823 = t804 * t824;
t775 = qJDD(2) * t800 + t823;
t713 = (-t775 + t823) * pkin(9) + (t800 * t804 * t806 + qJDD(4)) * pkin(4) + t716;
t717 = t804 * t720 + t800 * t760;
t776 = qJDD(2) * t804 - t800 * t824;
t826 = qJD(2) * t800;
t783 = qJD(4) * pkin(4) - pkin(9) * t826;
t790 = t804 ^ 2;
t714 = -pkin(4) * t790 * t806 + pkin(9) * t776 - qJD(4) * t783 + t717;
t709 = t799 * t713 + t803 * t714;
t769 = (t799 * t804 + t800 * t803) * qJD(2);
t736 = -qJD(5) * t769 - t775 * t799 + t776 * t803;
t751 = mrSges(6,1) * t768 + mrSges(6,2) * t769;
t789 = qJD(4) + qJD(5);
t759 = mrSges(6,1) * t789 - mrSges(6,3) * t769;
t788 = qJDD(4) + qJDD(5);
t752 = pkin(5) * t768 - pkin(10) * t769;
t787 = t789 ^ 2;
t706 = -pkin(5) * t787 + pkin(10) * t788 - t752 * t768 + t709;
t721 = t795 * t743 - t792 * t744;
t814 = -qJDD(2) * pkin(3) - t721;
t715 = -t776 * pkin(4) + t783 * t826 + (-pkin(9) * t790 - pkin(8)) * t806 + t814;
t737 = -qJD(5) * t768 + t775 * t803 + t776 * t799;
t710 = (t768 * t789 - t737) * pkin(10) + (t769 * t789 - t736) * pkin(5) + t715;
t798 = sin(qJ(6));
t802 = cos(qJ(6));
t703 = -t706 * t798 + t710 * t802;
t753 = -t769 * t798 + t789 * t802;
t725 = qJD(6) * t753 + t737 * t802 + t788 * t798;
t754 = t769 * t802 + t789 * t798;
t730 = -mrSges(7,1) * t753 + mrSges(7,2) * t754;
t735 = qJDD(6) - t736;
t762 = qJD(6) + t768;
t738 = -mrSges(7,2) * t762 + mrSges(7,3) * t753;
t699 = m(7) * t703 + mrSges(7,1) * t735 - mrSges(7,3) * t725 - t730 * t754 + t738 * t762;
t704 = t706 * t802 + t710 * t798;
t724 = -qJD(6) * t754 - t737 * t798 + t788 * t802;
t739 = mrSges(7,1) * t762 - mrSges(7,3) * t754;
t700 = m(7) * t704 - mrSges(7,2) * t735 + mrSges(7,3) * t724 + t730 * t753 - t739 * t762;
t818 = -t699 * t798 + t802 * t700;
t686 = m(6) * t709 - mrSges(6,2) * t788 + mrSges(6,3) * t736 - t751 * t768 - t759 * t789 + t818;
t708 = t713 * t803 - t714 * t799;
t758 = -mrSges(6,2) * t789 - mrSges(6,3) * t768;
t705 = -pkin(5) * t788 - pkin(10) * t787 + t752 * t769 - t708;
t812 = -m(7) * t705 + t724 * mrSges(7,1) - mrSges(7,2) * t725 + t753 * t738 - t739 * t754;
t695 = m(6) * t708 + mrSges(6,1) * t788 - mrSges(6,3) * t737 - t751 * t769 + t758 * t789 + t812;
t681 = t799 * t686 + t803 * t695;
t766 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t800 + Ifges(5,2) * t804) * qJD(2);
t767 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t800 + Ifges(5,4) * t804) * qJD(2);
t726 = Ifges(7,5) * t754 + Ifges(7,6) * t753 + Ifges(7,3) * t762;
t728 = Ifges(7,1) * t754 + Ifges(7,4) * t753 + Ifges(7,5) * t762;
t692 = -mrSges(7,1) * t705 + mrSges(7,3) * t704 + Ifges(7,4) * t725 + Ifges(7,2) * t724 + Ifges(7,6) * t735 - t726 * t754 + t728 * t762;
t727 = Ifges(7,4) * t754 + Ifges(7,2) * t753 + Ifges(7,6) * t762;
t693 = mrSges(7,2) * t705 - mrSges(7,3) * t703 + Ifges(7,1) * t725 + Ifges(7,4) * t724 + Ifges(7,5) * t735 + t726 * t753 - t727 * t762;
t746 = Ifges(6,4) * t769 - Ifges(6,2) * t768 + Ifges(6,6) * t789;
t747 = Ifges(6,1) * t769 - Ifges(6,4) * t768 + Ifges(6,5) * t789;
t810 = -mrSges(6,1) * t708 + mrSges(6,2) * t709 - Ifges(6,5) * t737 - Ifges(6,6) * t736 - Ifges(6,3) * t788 - pkin(5) * t812 - pkin(10) * t818 - t802 * t692 - t798 * t693 - t769 * t746 - t768 * t747;
t832 = mrSges(5,1) * t716 - mrSges(5,2) * t717 + Ifges(5,5) * t775 + Ifges(5,6) * t776 + Ifges(5,3) * qJDD(4) + pkin(4) * t681 + (t766 * t800 - t767 * t804) * qJD(2) - t810;
t774 = (-mrSges(5,1) * t804 + mrSges(5,2) * t800) * qJD(2);
t825 = qJD(2) * t804;
t780 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t825;
t679 = m(5) * t716 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t775 + qJD(4) * t780 - t774 * t826 + t681;
t779 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t826;
t819 = t803 * t686 - t695 * t799;
t680 = m(5) * t717 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t776 - qJD(4) * t779 + t774 * t825 + t819;
t820 = -t679 * t800 + t804 * t680;
t669 = m(4) * t722 - mrSges(4,1) * t806 - qJDD(2) * mrSges(4,2) + t820;
t719 = -t806 * pkin(8) + t814;
t688 = t802 * t699 + t798 * t700;
t811 = m(6) * t715 - t736 * mrSges(6,1) + mrSges(6,2) * t737 + t768 * t758 + t759 * t769 + t688;
t808 = -m(5) * t719 + t776 * mrSges(5,1) - mrSges(5,2) * t775 - t779 * t826 + t780 * t825 - t811;
t683 = m(4) * t721 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t806 + t808;
t666 = t792 * t669 + t795 * t683;
t664 = m(3) * t748 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t806 + t666;
t821 = t795 * t669 - t683 * t792;
t665 = m(3) * t749 - mrSges(3,1) * t806 - qJDD(2) * mrSges(3,2) + t821;
t673 = t804 * t679 + t800 * t680;
t672 = m(4) * t760 + t673;
t671 = m(3) * t761 + t672;
t651 = t664 * t828 + t665 * t829 - t671 * t794;
t649 = m(2) * t777 + t651;
t655 = -t664 * t801 + t805 * t665;
t654 = m(2) * t778 + t655;
t827 = t796 * t649 + t793 * t654;
t650 = t664 * t830 + t665 * t831 + t797 * t671;
t822 = -t649 * t793 + t796 * t654;
t817 = m(2) * t791 + t650;
t745 = Ifges(6,5) * t769 - Ifges(6,6) * t768 + Ifges(6,3) * t789;
t674 = mrSges(6,2) * t715 - mrSges(6,3) * t708 + Ifges(6,1) * t737 + Ifges(6,4) * t736 + Ifges(6,5) * t788 - pkin(10) * t688 - t692 * t798 + t693 * t802 - t745 * t768 - t746 * t789;
t809 = mrSges(7,1) * t703 - mrSges(7,2) * t704 + Ifges(7,5) * t725 + Ifges(7,6) * t724 + Ifges(7,3) * t735 + t727 * t754 - t728 * t753;
t675 = -mrSges(6,1) * t715 + mrSges(6,3) * t709 + Ifges(6,4) * t737 + Ifges(6,2) * t736 + Ifges(6,6) * t788 - pkin(5) * t688 - t745 * t769 + t747 * t789 - t809;
t765 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t800 + Ifges(5,6) * t804) * qJD(2);
t657 = -mrSges(5,1) * t719 + mrSges(5,3) * t717 + Ifges(5,4) * t775 + Ifges(5,2) * t776 + Ifges(5,6) * qJDD(4) - pkin(4) * t811 + pkin(9) * t819 + qJD(4) * t767 + t799 * t674 + t803 * t675 - t765 * t826;
t658 = mrSges(5,2) * t719 - mrSges(5,3) * t716 + Ifges(5,1) * t775 + Ifges(5,4) * t776 + Ifges(5,5) * qJDD(4) - pkin(9) * t681 - qJD(4) * t766 + t674 * t803 - t675 * t799 + t765 * t825;
t647 = mrSges(4,2) * t760 - mrSges(4,3) * t721 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t806 - pkin(8) * t673 - t657 * t800 + t658 * t804;
t656 = -mrSges(4,1) * t760 + mrSges(4,3) * t722 + t806 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t673 - t832;
t644 = -mrSges(3,1) * t761 + mrSges(3,3) * t749 + t806 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t672 + qJ(3) * t821 + t792 * t647 + t795 * t656;
t645 = mrSges(3,2) * t761 - mrSges(3,3) * t748 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t806 - qJ(3) * t666 + t647 * t795 - t656 * t792;
t813 = pkin(7) * t655 + t644 * t805 + t645 * t801;
t646 = mrSges(3,1) * t748 - mrSges(3,2) * t749 + mrSges(4,1) * t721 - mrSges(4,2) * t722 + t800 * t658 + t804 * t657 + pkin(3) * t808 + pkin(8) * t820 + pkin(2) * t666 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t643 = mrSges(2,2) * t791 - mrSges(2,3) * t777 - t801 * t644 + t805 * t645 + (-t650 * t794 - t651 * t797) * pkin(7);
t642 = -mrSges(2,1) * t791 + mrSges(2,3) * t778 - pkin(1) * t650 - t794 * t646 + t797 * t813;
t1 = [-m(1) * g(1) + t822; -m(1) * g(2) + t827; -m(1) * g(3) + t817; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t827 - t793 * t642 + t796 * t643; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t822 + t796 * t642 + t793 * t643; -mrSges(1,1) * g(2) + mrSges(2,1) * t777 + mrSges(1,2) * g(1) - mrSges(2,2) * t778 + pkin(1) * t651 + t797 * t646 + t794 * t813; t817; t646; t672; t832; -t810; t809;];
tauJB  = t1;
