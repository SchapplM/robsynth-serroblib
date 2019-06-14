% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 02:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_invdynB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:52:31
% EndTime: 2019-05-05 01:54:00
% DurationCPUTime: 92.53s
% Computational Cost: add. (1530962->357), mult. (4323227->502), div. (0->0), fcn. (3689140->18), ass. (0->177)
t751 = sin(pkin(14));
t754 = sin(pkin(7));
t756 = cos(pkin(14));
t759 = cos(pkin(7));
t763 = sin(qJ(4));
t758 = cos(pkin(8));
t767 = cos(qJ(4));
t801 = t758 * t767;
t753 = sin(pkin(8));
t808 = t753 * t767;
t771 = t754 * (-t751 * t763 + t756 * t801) + t759 * t808;
t723 = t771 * qJD(2);
t802 = t758 * t763;
t809 = t753 * t763;
t773 = t759 * t809 + (t751 * t767 + t756 * t802) * t754;
t724 = t773 * qJD(2);
t707 = -t724 * qJD(4) + qJDD(2) * t771;
t752 = sin(pkin(13));
t757 = cos(pkin(13));
t746 = t752 * g(1) - t757 * g(2);
t747 = -t757 * g(1) - t752 * g(2);
t750 = -g(3) + qJDD(1);
t768 = cos(qJ(2));
t760 = cos(pkin(6));
t764 = sin(qJ(2));
t800 = t760 * t764;
t755 = sin(pkin(6));
t805 = t755 * t764;
t720 = t746 * t800 + t768 * t747 + t750 * t805;
t769 = qJD(2) ^ 2;
t812 = qJ(3) * t754;
t718 = -t769 * pkin(2) + qJDD(2) * t812 + t720;
t814 = pkin(10) * t751;
t785 = -pkin(3) * t756 - t753 * t814;
t797 = qJD(2) * t754;
t733 = t785 * t797;
t806 = t754 * t758;
t780 = pkin(10) * (t753 * t759 + t756 * t806);
t734 = qJD(2) * t780;
t799 = t760 * t768;
t804 = t755 * t768;
t719 = t746 * t799 - t764 * t747 + t750 * t804;
t717 = qJDD(2) * pkin(2) + t769 * t812 + t719;
t736 = -t755 * t746 + t760 * t750;
t792 = qJD(3) * t797;
t803 = t756 * t759;
t807 = t754 * t756;
t793 = t717 * t803 + t736 * t807 - 0.2e1 * t751 * t792;
t676 = (pkin(3) * qJDD(2) + qJD(2) * t734) * t759 + (-t718 + (-pkin(10) * qJDD(2) * t758 - qJD(2) * t733) * t754) * t751 + t793;
t810 = t751 * t759;
t811 = t751 * t754;
t692 = t717 * t810 + t736 * t811 + (t718 + 0.2e1 * t792) * t756;
t738 = (pkin(3) * t759 - t806 * t814) * qJD(2);
t677 = (t733 * t807 - t738 * t759) * qJD(2) + qJDD(2) * t780 + t692;
t795 = t759 * t736 + qJDD(3);
t690 = (-t717 + t785 * qJDD(2) + (-t734 * t756 + t738 * t751) * qJD(2)) * t754 + t795;
t662 = -t763 * t677 + (t676 * t758 + t690 * t753) * t767;
t813 = Ifges(4,3) * t759;
t663 = t676 * t802 + t767 * t677 + t690 * t809;
t705 = -t723 * mrSges(5,1) + t724 * mrSges(5,2);
t781 = -t753 * t807 + t758 * t759;
t735 = qJD(2) * t781 + qJD(4);
t713 = t735 * mrSges(5,1) - t724 * mrSges(5,3);
t732 = qJDD(2) * t781 + qJDD(4);
t706 = -t723 * pkin(4) - t724 * pkin(11);
t731 = t735 ^ 2;
t659 = -t731 * pkin(4) + t732 * pkin(11) + t723 * t706 + t663;
t664 = -t753 * t676 + t758 * t690;
t708 = t723 * qJD(4) + qJDD(2) * t773;
t661 = (-t723 * t735 - t708) * pkin(11) + (t724 * t735 - t707) * pkin(4) + t664;
t762 = sin(qJ(5));
t766 = cos(qJ(5));
t655 = t766 * t659 + t762 * t661;
t711 = t766 * t724 + t762 * t735;
t686 = -t711 * qJD(5) - t762 * t708 + t766 * t732;
t710 = -t762 * t724 + t766 * t735;
t693 = -t710 * mrSges(6,1) + t711 * mrSges(6,2);
t722 = qJD(5) - t723;
t699 = t722 * mrSges(6,1) - t711 * mrSges(6,3);
t704 = qJDD(5) - t707;
t694 = -t710 * pkin(5) - t711 * pkin(12);
t721 = t722 ^ 2;
t653 = -t721 * pkin(5) + t704 * pkin(12) + t710 * t694 + t655;
t658 = -t732 * pkin(4) - t731 * pkin(11) + t724 * t706 - t662;
t687 = t710 * qJD(5) + t766 * t708 + t762 * t732;
t656 = (-t710 * t722 - t687) * pkin(12) + (t711 * t722 - t686) * pkin(5) + t658;
t761 = sin(qJ(6));
t765 = cos(qJ(6));
t650 = -t761 * t653 + t765 * t656;
t696 = -t761 * t711 + t765 * t722;
t667 = t696 * qJD(6) + t765 * t687 + t761 * t704;
t697 = t765 * t711 + t761 * t722;
t675 = -t696 * mrSges(7,1) + t697 * mrSges(7,2);
t709 = qJD(6) - t710;
t678 = -t709 * mrSges(7,2) + t696 * mrSges(7,3);
t684 = qJDD(6) - t686;
t648 = m(7) * t650 + t684 * mrSges(7,1) - t667 * mrSges(7,3) - t697 * t675 + t709 * t678;
t651 = t765 * t653 + t761 * t656;
t666 = -t697 * qJD(6) - t761 * t687 + t765 * t704;
t679 = t709 * mrSges(7,1) - t697 * mrSges(7,3);
t649 = m(7) * t651 - t684 * mrSges(7,2) + t666 * mrSges(7,3) + t696 * t675 - t709 * t679;
t789 = -t761 * t648 + t765 * t649;
t641 = m(6) * t655 - t704 * mrSges(6,2) + t686 * mrSges(6,3) + t710 * t693 - t722 * t699 + t789;
t654 = -t762 * t659 + t766 * t661;
t698 = -t722 * mrSges(6,2) + t710 * mrSges(6,3);
t652 = -t704 * pkin(5) - t721 * pkin(12) + t711 * t694 - t654;
t776 = -m(7) * t652 + t666 * mrSges(7,1) - t667 * mrSges(7,2) + t696 * t678 - t697 * t679;
t646 = m(6) * t654 + t704 * mrSges(6,1) - t687 * mrSges(6,3) - t711 * t693 + t722 * t698 + t776;
t790 = t766 * t641 - t762 * t646;
t632 = m(5) * t663 - t732 * mrSges(5,2) + t707 * mrSges(5,3) + t723 * t705 - t735 * t713 + t790;
t635 = t762 * t641 + t766 * t646;
t712 = -t735 * mrSges(5,2) + t723 * mrSges(5,3);
t634 = m(5) * t664 - t707 * mrSges(5,1) + t708 * mrSges(5,2) - t723 * t712 + t724 * t713 + t635;
t642 = t765 * t648 + t761 * t649;
t770 = -m(6) * t658 + t686 * mrSges(6,1) - t687 * mrSges(6,2) + t710 * t698 - t711 * t699 - t642;
t638 = m(5) * t662 + t732 * mrSges(5,1) - t708 * mrSges(5,3) - t724 * t705 + t735 * t712 + t770;
t621 = t632 * t802 - t753 * t634 + t638 * t801;
t691 = -t751 * t718 + t793;
t788 = -mrSges(4,1) * t756 + mrSges(4,2) * t751;
t737 = t788 * t797;
t783 = -mrSges(4,2) * t759 + mrSges(4,3) * t807;
t742 = t783 * qJD(2);
t784 = mrSges(4,1) * t759 - mrSges(4,3) * t811;
t617 = m(4) * t691 + t784 * qJDD(2) + (-t737 * t811 + t742 * t759) * qJD(2) + t621;
t620 = t632 * t809 + t758 * t634 + t638 * t808;
t703 = -t754 * t717 + t795;
t741 = t784 * qJD(2);
t619 = m(4) * t703 + (t788 * qJDD(2) + (t741 * t751 - t742 * t756) * qJD(2)) * t754 + t620;
t626 = t767 * t632 - t763 * t638;
t625 = m(4) * t692 + t783 * qJDD(2) + (t737 * t807 - t741 * t759) * qJD(2) + t626;
t606 = t617 * t803 - t754 * t619 + t625 * t810;
t602 = m(3) * t719 + qJDD(2) * mrSges(3,1) - t769 * mrSges(3,2) + t606;
t605 = t617 * t807 + t759 * t619 + t625 * t811;
t604 = m(3) * t736 + t605;
t611 = -t751 * t617 + t756 * t625;
t610 = m(3) * t720 - t769 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t611;
t592 = t602 * t799 - t755 * t604 + t610 * t800;
t590 = m(2) * t746 + t592;
t599 = -t764 * t602 + t768 * t610;
t598 = m(2) * t747 + t599;
t798 = t757 * t590 + t752 * t598;
t591 = t602 * t804 + t760 * t604 + t610 * t805;
t791 = -t752 * t590 + t757 * t598;
t787 = Ifges(4,5) * t751 + Ifges(4,6) * t756;
t668 = Ifges(7,5) * t697 + Ifges(7,6) * t696 + Ifges(7,3) * t709;
t670 = Ifges(7,1) * t697 + Ifges(7,4) * t696 + Ifges(7,5) * t709;
t643 = -mrSges(7,1) * t652 + mrSges(7,3) * t651 + Ifges(7,4) * t667 + Ifges(7,2) * t666 + Ifges(7,6) * t684 - t697 * t668 + t709 * t670;
t669 = Ifges(7,4) * t697 + Ifges(7,2) * t696 + Ifges(7,6) * t709;
t644 = mrSges(7,2) * t652 - mrSges(7,3) * t650 + Ifges(7,1) * t667 + Ifges(7,4) * t666 + Ifges(7,5) * t684 + t696 * t668 - t709 * t669;
t680 = Ifges(6,5) * t711 + Ifges(6,6) * t710 + Ifges(6,3) * t722;
t681 = Ifges(6,4) * t711 + Ifges(6,2) * t710 + Ifges(6,6) * t722;
t627 = mrSges(6,2) * t658 - mrSges(6,3) * t654 + Ifges(6,1) * t687 + Ifges(6,4) * t686 + Ifges(6,5) * t704 - pkin(12) * t642 - t761 * t643 + t765 * t644 + t710 * t680 - t722 * t681;
t682 = Ifges(6,1) * t711 + Ifges(6,4) * t710 + Ifges(6,5) * t722;
t628 = -mrSges(6,1) * t658 - mrSges(7,1) * t650 + mrSges(7,2) * t651 + mrSges(6,3) * t655 + Ifges(6,4) * t687 - Ifges(7,5) * t667 + Ifges(6,2) * t686 + Ifges(6,6) * t704 - Ifges(7,6) * t666 - Ifges(7,3) * t684 - pkin(5) * t642 - t697 * t669 + t696 * t670 - t711 * t680 + t722 * t682;
t701 = Ifges(5,4) * t724 + Ifges(5,2) * t723 + Ifges(5,6) * t735;
t702 = Ifges(5,1) * t724 + Ifges(5,4) * t723 + Ifges(5,5) * t735;
t612 = mrSges(5,1) * t662 - mrSges(5,2) * t663 + Ifges(5,5) * t708 + Ifges(5,6) * t707 + Ifges(5,3) * t732 + pkin(4) * t770 + pkin(11) * t790 + t762 * t627 + t766 * t628 + t724 * t701 - t723 * t702;
t774 = Ifges(4,6) * t759 + (Ifges(4,4) * t751 + Ifges(4,2) * t756) * t754;
t728 = t774 * qJD(2);
t775 = Ifges(4,5) * t759 + (Ifges(4,1) * t751 + Ifges(4,4) * t756) * t754;
t729 = t775 * qJD(2);
t700 = Ifges(5,5) * t724 + Ifges(5,6) * t723 + Ifges(5,3) * t735;
t613 = mrSges(5,2) * t664 - mrSges(5,3) * t662 + Ifges(5,1) * t708 + Ifges(5,4) * t707 + Ifges(5,5) * t732 - pkin(11) * t635 + t766 * t627 - t762 * t628 + t723 * t700 - t735 * t701;
t614 = Ifges(5,4) * t708 + Ifges(5,2) * t707 + Ifges(5,6) * t732 - t724 * t700 + t735 * t702 - mrSges(5,1) * t664 + mrSges(5,3) * t663 - Ifges(6,5) * t687 - Ifges(6,6) * t686 - Ifges(6,3) * t704 - t711 * t681 + t710 * t682 - mrSges(6,1) * t654 + mrSges(6,2) * t655 - t761 * t644 - t765 * t643 - pkin(5) * t776 - pkin(12) * t789 - pkin(4) * t635;
t778 = pkin(10) * t626 + t613 * t763 + t614 * t767;
t593 = qJDD(2) * t813 + mrSges(4,1) * t691 - mrSges(4,2) * t692 + pkin(3) * t621 + t758 * t612 + t778 * t753 + (t787 * qJDD(2) + (t728 * t751 - t729 * t756) * qJD(2)) * t754;
t727 = (t754 * t787 + t813) * qJD(2);
t594 = -mrSges(4,1) * t703 + mrSges(4,3) * t692 - pkin(3) * t620 - t753 * t612 + (-t727 * t811 + t729 * t759) * qJD(2) + t778 * t758 + t774 * qJDD(2);
t595 = mrSges(4,2) * t703 - mrSges(4,3) * t691 + t767 * t613 - t763 * t614 + (t727 * t807 - t728 * t759) * qJD(2) + (-t620 * t753 - t621 * t758) * pkin(10) + t775 * qJDD(2);
t777 = qJ(3) * t611 + t594 * t756 + t595 * t751;
t587 = -mrSges(3,1) * t736 + mrSges(3,3) * t720 + t769 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t605 - t754 * t593 + t759 * t777;
t588 = mrSges(3,2) * t736 - mrSges(3,3) * t719 + Ifges(3,5) * qJDD(2) - t769 * Ifges(3,6) - t751 * t594 + t756 * t595 + (-t605 * t754 - t606 * t759) * qJ(3);
t779 = pkin(9) * t599 + t587 * t768 + t588 * t764;
t586 = mrSges(3,1) * t719 - mrSges(3,2) * t720 + Ifges(3,3) * qJDD(2) + pkin(2) * t606 + t759 * t593 + t754 * t777;
t585 = mrSges(2,2) * t750 - mrSges(2,3) * t746 - t764 * t587 + t768 * t588 + (-t591 * t755 - t592 * t760) * pkin(9);
t584 = -mrSges(2,1) * t750 + mrSges(2,3) * t747 - pkin(1) * t591 - t755 * t586 + t760 * t779;
t1 = [-m(1) * g(1) + t791; -m(1) * g(2) + t798; -m(1) * g(3) + m(2) * t750 + t591; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t798 - t752 * t584 + t757 * t585; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t791 + t757 * t584 + t752 * t585; -mrSges(1,1) * g(2) + mrSges(2,1) * t746 + mrSges(1,2) * g(1) - mrSges(2,2) * t747 + pkin(1) * t592 + t760 * t586 + t755 * t779;];
tauB  = t1;
