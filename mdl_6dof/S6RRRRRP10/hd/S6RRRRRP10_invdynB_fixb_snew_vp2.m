% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 06:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:20:43
% EndTime: 2019-05-08 06:21:29
% DurationCPUTime: 26.86s
% Computational Cost: add. (453047->375), mult. (961382->474), div. (0->0), fcn. (767874->12), ass. (0->155)
t801 = Ifges(6,1) + Ifges(7,1);
t794 = Ifges(6,4) - Ifges(7,5);
t793 = Ifges(7,4) + Ifges(6,5);
t800 = Ifges(6,2) + Ifges(7,3);
t792 = Ifges(6,6) - Ifges(7,6);
t799 = -Ifges(6,3) - Ifges(7,2);
t754 = sin(pkin(6));
t759 = sin(qJ(2));
t763 = cos(qJ(2));
t779 = qJD(1) * qJD(2);
t743 = (-qJDD(1) * t763 + t759 * t779) * t754;
t798 = cos(qJ(5));
t797 = pkin(8) * t754;
t755 = cos(pkin(6));
t796 = g(3) * t755;
t795 = -mrSges(6,3) - mrSges(7,2);
t791 = t754 * t759;
t790 = t754 * t763;
t789 = t755 * t759;
t788 = t755 * t763;
t760 = sin(qJ(1));
t764 = cos(qJ(1));
t747 = t760 * g(1) - g(2) * t764;
t765 = qJD(1) ^ 2;
t738 = qJDD(1) * pkin(1) + t765 * t797 + t747;
t748 = -g(1) * t764 - g(2) * t760;
t739 = -pkin(1) * t765 + qJDD(1) * t797 + t748;
t782 = t738 * t789 + t763 * t739;
t711 = -g(3) * t791 + t782;
t751 = qJD(1) * t755 + qJD(2);
t781 = qJD(1) * t754;
t777 = t759 * t781;
t736 = mrSges(3,1) * t751 - mrSges(3,3) * t777;
t740 = (-mrSges(3,1) * t763 + mrSges(3,2) * t759) * t781;
t750 = qJDD(1) * t755 + qJDD(2);
t741 = (-pkin(2) * t763 - pkin(9) * t759) * t781;
t749 = t751 ^ 2;
t780 = qJD(1) * t763;
t690 = -pkin(2) * t749 + pkin(9) * t750 + (-g(3) * t759 + t741 * t780) * t754 + t782;
t742 = (qJDD(1) * t759 + t763 * t779) * t754;
t691 = pkin(2) * t743 - pkin(9) * t742 - t796 + (-t738 + (pkin(2) * t759 - pkin(9) * t763) * t751 * qJD(1)) * t754;
t758 = sin(qJ(3));
t762 = cos(qJ(3));
t658 = t762 * t690 + t758 * t691;
t731 = t751 * t758 + t762 * t777;
t708 = -t731 * qJD(3) - t758 * t742 + t750 * t762;
t730 = t751 * t762 - t758 * t777;
t712 = -mrSges(4,1) * t730 + mrSges(4,2) * t731;
t776 = t754 * t780;
t746 = qJD(3) - t776;
t719 = mrSges(4,1) * t746 - mrSges(4,3) * t731;
t735 = qJDD(3) + t743;
t713 = -pkin(3) * t730 - pkin(10) * t731;
t744 = t746 ^ 2;
t653 = -pkin(3) * t744 + pkin(10) * t735 + t713 * t730 + t658;
t710 = -g(3) * t790 + t738 * t788 - t759 * t739;
t689 = -pkin(2) * t750 - pkin(9) * t749 + t741 * t777 - t710;
t709 = qJD(3) * t730 + t742 * t762 + t750 * t758;
t656 = (-t730 * t746 - t709) * pkin(10) + (t731 * t746 - t708) * pkin(3) + t689;
t757 = sin(qJ(4));
t761 = cos(qJ(4));
t643 = -t653 * t757 + t761 * t656;
t716 = -t731 * t757 + t746 * t761;
t676 = qJD(4) * t716 + t709 * t761 + t735 * t757;
t706 = qJDD(4) - t708;
t717 = t731 * t761 + t746 * t757;
t729 = qJD(4) - t730;
t640 = (t716 * t729 - t676) * pkin(11) + (t716 * t717 + t706) * pkin(4) + t643;
t644 = t761 * t653 + t757 * t656;
t675 = -qJD(4) * t717 - t709 * t757 + t735 * t761;
t698 = pkin(4) * t729 - pkin(11) * t717;
t715 = t716 ^ 2;
t642 = -pkin(4) * t715 + pkin(11) * t675 - t698 * t729 + t644;
t756 = sin(qJ(5));
t636 = t756 * t640 + t798 * t642;
t693 = t756 * t716 + t798 * t717;
t649 = qJD(5) * t693 - t798 * t675 + t676 * t756;
t727 = qJD(5) + t729;
t679 = mrSges(6,1) * t727 - mrSges(6,3) * t693;
t692 = -t798 * t716 + t717 * t756;
t701 = qJDD(5) + t706;
t668 = pkin(5) * t692 - qJ(6) * t693;
t725 = t727 ^ 2;
t633 = -pkin(5) * t725 + qJ(6) * t701 + 0.2e1 * qJD(6) * t727 - t668 * t692 + t636;
t680 = -mrSges(7,1) * t727 + mrSges(7,2) * t693;
t778 = m(7) * t633 + t701 * mrSges(7,3) + t727 * t680;
t669 = mrSges(7,1) * t692 - mrSges(7,3) * t693;
t783 = -mrSges(6,1) * t692 - mrSges(6,2) * t693 - t669;
t626 = m(6) * t636 - mrSges(6,2) * t701 + t795 * t649 - t679 * t727 + t783 * t692 + t778;
t635 = t798 * t640 - t756 * t642;
t650 = -t692 * qJD(5) + t756 * t675 + t798 * t676;
t678 = -mrSges(6,2) * t727 - mrSges(6,3) * t692;
t634 = -t701 * pkin(5) - t725 * qJ(6) + t693 * t668 + qJDD(6) - t635;
t677 = -mrSges(7,2) * t692 + mrSges(7,3) * t727;
t771 = -m(7) * t634 + t701 * mrSges(7,1) + t727 * t677;
t628 = m(6) * t635 + mrSges(6,1) * t701 + t795 * t650 + t678 * t727 + t783 * t693 + t771;
t623 = t756 * t626 + t798 * t628;
t694 = -mrSges(5,1) * t716 + mrSges(5,2) * t717;
t696 = -mrSges(5,2) * t729 + mrSges(5,3) * t716;
t619 = m(5) * t643 + mrSges(5,1) * t706 - mrSges(5,3) * t676 - t694 * t717 + t696 * t729 + t623;
t697 = mrSges(5,1) * t729 - mrSges(5,3) * t717;
t772 = t798 * t626 - t628 * t756;
t620 = m(5) * t644 - mrSges(5,2) * t706 + mrSges(5,3) * t675 + t694 * t716 - t697 * t729 + t772;
t773 = -t619 * t757 + t761 * t620;
t616 = m(4) * t658 - mrSges(4,2) * t735 + mrSges(4,3) * t708 + t712 * t730 - t719 * t746 + t773;
t657 = -t758 * t690 + t691 * t762;
t718 = -mrSges(4,2) * t746 + mrSges(4,3) * t730;
t652 = -pkin(3) * t735 - pkin(10) * t744 + t731 * t713 - t657;
t645 = -pkin(4) * t675 - pkin(11) * t715 + t717 * t698 + t652;
t638 = -0.2e1 * qJD(6) * t693 + (t692 * t727 - t650) * qJ(6) + (t693 * t727 + t649) * pkin(5) + t645;
t631 = m(7) * t638 + t649 * mrSges(7,1) - t650 * mrSges(7,3) + t692 * t677 - t693 * t680;
t768 = m(6) * t645 + t649 * mrSges(6,1) + mrSges(6,2) * t650 + t692 * t678 + t679 * t693 + t631;
t766 = -m(5) * t652 + t675 * mrSges(5,1) - mrSges(5,2) * t676 + t716 * t696 - t697 * t717 - t768;
t630 = m(4) * t657 + mrSges(4,1) * t735 - mrSges(4,3) * t709 - t712 * t731 + t718 * t746 + t766;
t774 = t762 * t616 - t630 * t758;
t607 = m(3) * t711 - mrSges(3,2) * t750 - mrSges(3,3) * t743 - t736 * t751 + t740 * t776 + t774;
t610 = t758 * t616 + t762 * t630;
t723 = -t738 * t754 - t796;
t737 = -mrSges(3,2) * t751 + mrSges(3,3) * t776;
t609 = m(3) * t723 + mrSges(3,1) * t743 + mrSges(3,2) * t742 + (t736 * t759 - t737 * t763) * t781 + t610;
t617 = t619 * t761 + t620 * t757;
t767 = -m(4) * t689 + t708 * mrSges(4,1) - mrSges(4,2) * t709 + t730 * t718 - t719 * t731 - t617;
t613 = m(3) * t710 + mrSges(3,1) * t750 - mrSges(3,3) * t742 + t737 * t751 - t740 * t777 + t767;
t596 = t607 * t789 - t609 * t754 + t613 * t788;
t594 = m(2) * t747 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t765 + t596;
t600 = t763 * t607 - t613 * t759;
t599 = m(2) * t748 - mrSges(2,1) * t765 - qJDD(1) * mrSges(2,2) + t600;
t787 = t764 * t594 + t760 * t599;
t786 = t692 * t800 - t693 * t794 - t727 * t792;
t785 = t692 * t792 - t693 * t793 + t727 * t799;
t784 = -t794 * t692 + t693 * t801 + t793 * t727;
t595 = t607 * t791 + t755 * t609 + t613 * t790;
t775 = -t594 * t760 + t764 * t599;
t621 = -mrSges(6,1) * t645 - mrSges(7,1) * t638 + mrSges(7,2) * t633 + mrSges(6,3) * t636 - pkin(5) * t631 - t649 * t800 + t794 * t650 + t785 * t693 + t792 * t701 + t784 * t727;
t622 = mrSges(6,2) * t645 + mrSges(7,2) * t634 - mrSges(6,3) * t635 - mrSges(7,3) * t638 - qJ(6) * t631 - t794 * t649 + t650 * t801 + t785 * t692 + t793 * t701 + t786 * t727;
t681 = Ifges(5,5) * t717 + Ifges(5,6) * t716 + Ifges(5,3) * t729;
t683 = Ifges(5,1) * t717 + Ifges(5,4) * t716 + Ifges(5,5) * t729;
t602 = -mrSges(5,1) * t652 + mrSges(5,3) * t644 + Ifges(5,4) * t676 + Ifges(5,2) * t675 + Ifges(5,6) * t706 - pkin(4) * t768 + pkin(11) * t772 + t798 * t621 + t756 * t622 - t717 * t681 + t729 * t683;
t682 = Ifges(5,4) * t717 + Ifges(5,2) * t716 + Ifges(5,6) * t729;
t603 = mrSges(5,2) * t652 - mrSges(5,3) * t643 + Ifges(5,1) * t676 + Ifges(5,4) * t675 + Ifges(5,5) * t706 - pkin(11) * t623 - t756 * t621 + t798 * t622 + t716 * t681 - t729 * t682;
t702 = Ifges(4,5) * t731 + Ifges(4,6) * t730 + Ifges(4,3) * t746;
t703 = Ifges(4,4) * t731 + Ifges(4,2) * t730 + Ifges(4,6) * t746;
t592 = mrSges(4,2) * t689 - mrSges(4,3) * t657 + Ifges(4,1) * t709 + Ifges(4,4) * t708 + Ifges(4,5) * t735 - pkin(10) * t617 - t602 * t757 + t603 * t761 + t702 * t730 - t703 * t746;
t704 = Ifges(4,1) * t731 + Ifges(4,4) * t730 + Ifges(4,5) * t746;
t601 = -pkin(5) * t771 + (qJ(6) * mrSges(7,2) + t792) * t649 + (pkin(5) * mrSges(7,2) - t793) * t650 - pkin(3) * t617 + (qJ(6) * t669 - t784) * t692 + (pkin(5) * t669 + t786) * t693 + t799 * t701 + Ifges(4,6) * t735 - t731 * t702 - t717 * t682 + t716 * t683 - Ifges(5,3) * t706 + Ifges(4,2) * t708 + Ifges(4,4) * t709 - mrSges(4,1) * t689 - Ifges(5,6) * t675 - Ifges(5,5) * t676 + mrSges(4,3) * t658 - mrSges(5,1) * t643 + mrSges(5,2) * t644 + mrSges(6,2) * t636 - mrSges(6,1) * t635 + mrSges(7,1) * t634 - mrSges(7,3) * t633 - pkin(4) * t623 + t746 * t704 - qJ(6) * t778;
t720 = Ifges(3,3) * t751 + (Ifges(3,5) * t759 + Ifges(3,6) * t763) * t781;
t721 = Ifges(3,6) * t751 + (Ifges(3,4) * t759 + Ifges(3,2) * t763) * t781;
t590 = mrSges(3,2) * t723 - mrSges(3,3) * t710 + Ifges(3,1) * t742 - Ifges(3,4) * t743 + Ifges(3,5) * t750 - pkin(9) * t610 + t592 * t762 - t601 * t758 + t720 * t776 - t721 * t751;
t722 = Ifges(3,5) * t751 + (Ifges(3,1) * t759 + Ifges(3,4) * t763) * t781;
t591 = Ifges(3,4) * t742 - Ifges(3,2) * t743 + Ifges(3,6) * t750 - t720 * t777 + t751 * t722 - mrSges(3,1) * t723 + mrSges(3,3) * t711 - Ifges(4,5) * t709 - Ifges(4,6) * t708 - Ifges(4,3) * t735 - t731 * t703 + t730 * t704 - mrSges(4,1) * t657 + mrSges(4,2) * t658 - t757 * t603 - t761 * t602 - pkin(3) * t766 - pkin(10) * t773 - pkin(2) * t610;
t769 = pkin(8) * t600 + t590 * t759 + t591 * t763;
t589 = Ifges(3,5) * t742 - Ifges(3,6) * t743 + Ifges(3,3) * t750 + mrSges(3,1) * t710 - mrSges(3,2) * t711 + t758 * t592 + t762 * t601 + pkin(2) * t767 + pkin(9) * t774 + (t721 * t759 - t722 * t763) * t781;
t588 = -mrSges(2,2) * g(3) - mrSges(2,3) * t747 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t765 + t590 * t763 - t591 * t759 + (-t595 * t754 - t596 * t755) * pkin(8);
t587 = mrSges(2,1) * g(3) + mrSges(2,3) * t748 + Ifges(2,5) * t765 + Ifges(2,6) * qJDD(1) - pkin(1) * t595 - t589 * t754 + t769 * t755;
t1 = [-m(1) * g(1) + t775; -m(1) * g(2) + t787; (-m(1) - m(2)) * g(3) + t595; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t787 - t760 * t587 + t764 * t588; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t775 + t764 * t587 + t760 * t588; -mrSges(1,1) * g(2) + mrSges(2,1) * t747 + mrSges(1,2) * g(1) - mrSges(2,2) * t748 + Ifges(2,3) * qJDD(1) + pkin(1) * t596 + t589 * t755 + t769 * t754;];
tauB  = t1;
