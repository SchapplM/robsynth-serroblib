% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 01:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:31:17
% EndTime: 2019-05-05 01:31:26
% DurationCPUTime: 6.71s
% Computational Cost: add. (114754->293), mult. (215171->365), div. (0->0), fcn. (143388->12), ass. (0->131)
t782 = sin(pkin(11));
t784 = cos(pkin(11));
t764 = g(1) * t782 - g(2) * t784;
t765 = -g(1) * t784 - g(2) * t782;
t779 = -g(3) + qJDD(1);
t793 = cos(qJ(2));
t785 = cos(pkin(6));
t789 = sin(qJ(2));
t818 = t785 * t789;
t783 = sin(pkin(6));
t819 = t783 * t789;
t722 = t764 * t818 + t793 * t765 + t779 * t819;
t825 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t722;
t721 = -t789 * t765 + (t764 * t785 + t779 * t783) * t793;
t824 = -pkin(2) - pkin(8);
t823 = mrSges(3,1) - mrSges(4,2);
t822 = (-Ifges(4,4) + Ifges(3,5));
t821 = Ifges(4,5) - Ifges(3,6);
t795 = qJD(2) ^ 2;
t800 = -t795 * qJ(3) + qJDD(3) - t721;
t716 = qJDD(2) * t824 + t800;
t739 = -t764 * t783 + t779 * t785;
t788 = sin(qJ(4));
t792 = cos(qJ(4));
t710 = t788 * t716 + t792 * t739;
t760 = (mrSges(5,1) * t788 + mrSges(5,2) * t792) * qJD(2);
t815 = qJD(2) * qJD(4);
t771 = t792 * t815;
t762 = -t788 * qJDD(2) - t771;
t816 = qJD(2) * t792;
t767 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t816;
t773 = t788 * qJD(2);
t761 = (pkin(4) * t788 - pkin(9) * t792) * qJD(2);
t794 = qJD(4) ^ 2;
t701 = -pkin(4) * t794 + qJDD(4) * pkin(9) - t761 * t773 + t710;
t715 = t795 * t824 - t825;
t813 = t788 * t815;
t763 = qJDD(2) * t792 - t813;
t704 = (-t763 + t813) * pkin(9) + (-t762 + t771) * pkin(4) + t715;
t787 = sin(qJ(5));
t791 = cos(qJ(5));
t690 = -t787 * t701 + t791 * t704;
t758 = qJD(4) * t791 - t787 * t816;
t729 = qJD(5) * t758 + qJDD(4) * t787 + t763 * t791;
t755 = qJDD(5) - t762;
t759 = qJD(4) * t787 + t791 * t816;
t770 = t773 + qJD(5);
t688 = (t758 * t770 - t729) * pkin(10) + (t758 * t759 + t755) * pkin(5) + t690;
t691 = t791 * t701 + t787 * t704;
t728 = -qJD(5) * t759 + qJDD(4) * t791 - t763 * t787;
t737 = pkin(5) * t770 - pkin(10) * t759;
t754 = t758 ^ 2;
t689 = -pkin(5) * t754 + pkin(10) * t728 - t737 * t770 + t691;
t786 = sin(qJ(6));
t790 = cos(qJ(6));
t686 = t688 * t790 - t689 * t786;
t730 = t758 * t790 - t759 * t786;
t698 = qJD(6) * t730 + t728 * t786 + t729 * t790;
t731 = t758 * t786 + t759 * t790;
t712 = -mrSges(7,1) * t730 + mrSges(7,2) * t731;
t769 = qJD(6) + t770;
t719 = -mrSges(7,2) * t769 + mrSges(7,3) * t730;
t748 = qJDD(6) + t755;
t682 = m(7) * t686 + mrSges(7,1) * t748 - t698 * mrSges(7,3) - t712 * t731 + t719 * t769;
t687 = t688 * t786 + t689 * t790;
t697 = -qJD(6) * t731 + t728 * t790 - t729 * t786;
t720 = mrSges(7,1) * t769 - mrSges(7,3) * t731;
t683 = m(7) * t687 - mrSges(7,2) * t748 + t697 * mrSges(7,3) + t712 * t730 - t720 * t769;
t675 = t790 * t682 + t786 * t683;
t732 = -mrSges(6,1) * t758 + mrSges(6,2) * t759;
t735 = -mrSges(6,2) * t770 + mrSges(6,3) * t758;
t673 = m(6) * t690 + mrSges(6,1) * t755 - mrSges(6,3) * t729 - t732 * t759 + t735 * t770 + t675;
t736 = mrSges(6,1) * t770 - mrSges(6,3) * t759;
t809 = -t682 * t786 + t790 * t683;
t674 = m(6) * t691 - mrSges(6,2) * t755 + mrSges(6,3) * t728 + t732 * t758 - t736 * t770 + t809;
t810 = -t673 * t787 + t791 * t674;
t667 = m(5) * t710 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t762 - qJD(4) * t767 - t760 * t773 + t810;
t709 = t716 * t792 - t788 * t739;
t766 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t773;
t700 = -qJDD(4) * pkin(4) - pkin(9) * t794 + t761 * t816 - t709;
t692 = -pkin(5) * t728 - pkin(10) * t754 + t737 * t759 + t700;
t802 = m(7) * t692 - t697 * mrSges(7,1) + t698 * mrSges(7,2) - t730 * t719 + t720 * t731;
t797 = -m(6) * t700 + t728 * mrSges(6,1) - mrSges(6,2) * t729 + t758 * t735 - t736 * t759 - t802;
t678 = m(5) * t709 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t763 + qJD(4) * t766 - t760 * t816 + t797;
t659 = t788 * t667 + t792 * t678;
t718 = -qJDD(2) * pkin(2) + t800;
t804 = -m(4) * t718 + (t795 * mrSges(4,3)) - t659;
t652 = m(3) * t721 - (t795 * mrSges(3,2)) + qJDD(2) * t823 + t804;
t820 = t652 * t793;
t811 = t792 * t667 - t678 * t788;
t657 = m(4) * t739 + t811;
t655 = m(3) * t739 + t657;
t717 = t795 * pkin(2) + t825;
t669 = t791 * t673 + t787 * t674;
t803 = -m(5) * t715 + mrSges(5,1) * t762 - t763 * mrSges(5,2) - t766 * t773 - t767 * t816 - t669;
t798 = -m(4) * t717 + (t795 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t803;
t665 = m(3) * t722 - (mrSges(3,1) * t795) - qJDD(2) * mrSges(3,2) + t798;
t643 = -t655 * t783 + t665 * t818 + t785 * t820;
t641 = m(2) * t764 + t643;
t648 = -t652 * t789 + t793 * t665;
t647 = m(2) * t765 + t648;
t817 = t784 * t641 + t782 * t647;
t642 = t785 * t655 + t665 * t819 + t783 * t820;
t812 = -t641 * t782 + t784 * t647;
t807 = m(2) * t779 + t642;
t706 = Ifges(7,5) * t731 + Ifges(7,6) * t730 + Ifges(7,3) * t769;
t708 = Ifges(7,1) * t731 + Ifges(7,4) * t730 + Ifges(7,5) * t769;
t676 = -mrSges(7,1) * t692 + mrSges(7,3) * t687 + Ifges(7,4) * t698 + Ifges(7,2) * t697 + Ifges(7,6) * t748 - t706 * t731 + t708 * t769;
t707 = Ifges(7,4) * t731 + Ifges(7,2) * t730 + Ifges(7,6) * t769;
t677 = mrSges(7,2) * t692 - mrSges(7,3) * t686 + Ifges(7,1) * t698 + Ifges(7,4) * t697 + Ifges(7,5) * t748 + t706 * t730 - t707 * t769;
t723 = Ifges(6,5) * t759 + Ifges(6,6) * t758 + Ifges(6,3) * t770;
t725 = Ifges(6,1) * t759 + Ifges(6,4) * t758 + Ifges(6,5) * t770;
t658 = -mrSges(6,1) * t700 + mrSges(6,3) * t691 + Ifges(6,4) * t729 + Ifges(6,2) * t728 + Ifges(6,6) * t755 - pkin(5) * t802 + pkin(10) * t809 + t790 * t676 + t786 * t677 - t759 * t723 + t770 * t725;
t724 = Ifges(6,4) * t759 + Ifges(6,2) * t758 + Ifges(6,6) * t770;
t661 = mrSges(6,2) * t700 - mrSges(6,3) * t690 + Ifges(6,1) * t729 + Ifges(6,4) * t728 + Ifges(6,5) * t755 - pkin(10) * t675 - t676 * t786 + t677 * t790 + t723 * t758 - t724 * t770;
t745 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t792 - Ifges(5,6) * t788) * qJD(2);
t746 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t792 - Ifges(5,2) * t788) * qJD(2);
t644 = mrSges(5,2) * t715 - mrSges(5,3) * t709 + Ifges(5,1) * t763 + Ifges(5,4) * t762 + Ifges(5,5) * qJDD(4) - pkin(9) * t669 - qJD(4) * t746 - t658 * t787 + t661 * t791 - t745 * t773;
t747 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t792 - Ifges(5,4) * t788) * qJD(2);
t801 = -mrSges(7,1) * t686 + mrSges(7,2) * t687 - Ifges(7,5) * t698 - Ifges(7,6) * t697 - Ifges(7,3) * t748 - t731 * t707 + t730 * t708;
t796 = mrSges(6,1) * t690 - mrSges(6,2) * t691 + Ifges(6,5) * t729 + Ifges(6,6) * t728 + Ifges(6,3) * t755 + pkin(5) * t675 + t759 * t724 - t758 * t725 - t801;
t649 = -mrSges(5,1) * t715 + mrSges(5,3) * t710 + Ifges(5,4) * t763 + Ifges(5,2) * t762 + Ifges(5,6) * qJDD(4) - pkin(4) * t669 + qJD(4) * t747 - t745 * t816 - t796;
t638 = -mrSges(4,1) * t717 + mrSges(3,3) * t722 - pkin(2) * t657 - pkin(3) * t803 - pkin(8) * t811 - qJDD(2) * t821 - t788 * t644 - t792 * t649 - t739 * t823 + (t795 * t822);
t799 = mrSges(5,1) * t709 - mrSges(5,2) * t710 + Ifges(5,5) * t763 + Ifges(5,6) * t762 + Ifges(5,3) * qJDD(4) + pkin(4) * t797 + pkin(9) * t810 + t791 * t658 + t787 * t661 + t746 * t816 + t747 * t773;
t639 = t822 * qJDD(2) + t799 - mrSges(3,3) * t721 + mrSges(4,1) * t718 + pkin(3) * t659 - qJ(3) * t657 + t821 * t795 + (mrSges(3,2) - mrSges(4,3)) * t739;
t805 = pkin(7) * t648 + t638 * t793 + t639 * t789;
t653 = qJDD(2) * mrSges(4,2) - t804;
t637 = mrSges(3,1) * t721 - mrSges(3,2) * t722 + mrSges(4,2) * t718 - mrSges(4,3) * t717 + t792 * t644 - t788 * t649 - pkin(8) * t659 - pkin(2) * t653 + qJ(3) * t798 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t636 = mrSges(2,2) * t779 - mrSges(2,3) * t764 - t789 * t638 + t793 * t639 + (-t642 * t783 - t643 * t785) * pkin(7);
t635 = -mrSges(2,1) * t779 + mrSges(2,3) * t765 - pkin(1) * t642 - t783 * t637 + t785 * t805;
t1 = [-m(1) * g(1) + t812; -m(1) * g(2) + t817; -m(1) * g(3) + t807; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t817 - t782 * t635 + t784 * t636; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t812 + t784 * t635 + t782 * t636; -mrSges(1,1) * g(2) + mrSges(2,1) * t764 + mrSges(1,2) * g(1) - mrSges(2,2) * t765 + pkin(1) * t643 + t785 * t637 + t783 * t805; t807; t637; t653; t799; t796; -t801;];
tauJB  = t1;
