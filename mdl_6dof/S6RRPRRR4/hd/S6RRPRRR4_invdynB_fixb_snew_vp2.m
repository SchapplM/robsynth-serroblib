% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 20:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:41:52
% EndTime: 2019-05-06 20:42:56
% DurationCPUTime: 45.50s
% Computational Cost: add. (705984->397), mult. (1842582->516), div. (0->0), fcn. (1482014->14), ass. (0->162)
t772 = sin(pkin(6));
t778 = sin(qJ(2));
t783 = cos(qJ(2));
t801 = qJD(1) * qJD(2);
t760 = (qJDD(1) * t778 + t783 * t801) * t772;
t774 = cos(pkin(6));
t766 = t774 * qJDD(1) + qJDD(2);
t767 = t774 * qJD(1) + qJD(2);
t779 = sin(qJ(1));
t784 = cos(qJ(1));
t763 = t779 * g(1) - t784 * g(2);
t785 = qJD(1) ^ 2;
t810 = pkin(8) * t772;
t757 = qJDD(1) * pkin(1) + t785 * t810 + t763;
t764 = -t784 * g(1) - t779 * g(2);
t758 = -t785 * pkin(1) + qJDD(1) * t810 + t764;
t805 = t774 * t783;
t791 = t757 * t805 - t778 * t758;
t809 = t772 ^ 2 * t785;
t698 = t766 * pkin(2) - t760 * qJ(3) + (pkin(2) * t778 * t809 + (qJ(3) * qJD(1) * t767 - g(3)) * t772) * t783 + t791;
t806 = t774 * t778;
t808 = t772 * t778;
t725 = -g(3) * t808 + t757 * t806 + t783 * t758;
t803 = qJD(1) * t772;
t799 = t778 * t803;
t754 = t767 * pkin(2) - qJ(3) * t799;
t761 = (qJDD(1) * t783 - t778 * t801) * t772;
t800 = t783 ^ 2 * t809;
t701 = -pkin(2) * t800 + t761 * qJ(3) - t767 * t754 + t725;
t771 = sin(pkin(12));
t773 = cos(pkin(12));
t751 = (t771 * t783 + t773 * t778) * t803;
t679 = -0.2e1 * qJD(3) * t751 + t773 * t698 - t771 * t701;
t807 = t772 * t783;
t798 = t783 * t803;
t750 = -t771 * t799 + t773 * t798;
t680 = 0.2e1 * qJD(3) * t750 + t771 * t698 + t773 * t701;
t726 = -t750 * mrSges(4,1) + t751 * mrSges(4,2);
t731 = -t771 * t760 + t773 * t761;
t737 = t767 * mrSges(4,1) - t751 * mrSges(4,3);
t727 = -t750 * pkin(3) - t751 * pkin(9);
t765 = t767 ^ 2;
t670 = -t765 * pkin(3) + t766 * pkin(9) + t750 * t727 + t680;
t741 = -t774 * g(3) - t772 * t757;
t711 = -t761 * pkin(2) - qJ(3) * t800 + t754 * t799 + qJDD(3) + t741;
t732 = t773 * t760 + t771 * t761;
t683 = (-t750 * t767 - t732) * pkin(9) + (t751 * t767 - t731) * pkin(3) + t711;
t777 = sin(qJ(4));
t782 = cos(qJ(4));
t662 = -t777 * t670 + t782 * t683;
t734 = -t777 * t751 + t782 * t767;
t709 = t734 * qJD(4) + t782 * t732 + t777 * t766;
t730 = qJDD(4) - t731;
t735 = t782 * t751 + t777 * t767;
t749 = qJD(4) - t750;
t659 = (t734 * t749 - t709) * pkin(10) + (t734 * t735 + t730) * pkin(4) + t662;
t663 = t782 * t670 + t777 * t683;
t708 = -t735 * qJD(4) - t777 * t732 + t782 * t766;
t719 = t749 * pkin(4) - t735 * pkin(10);
t733 = t734 ^ 2;
t661 = -t733 * pkin(4) + t708 * pkin(10) - t749 * t719 + t663;
t776 = sin(qJ(5));
t781 = cos(qJ(5));
t656 = t776 * t659 + t781 * t661;
t714 = t776 * t734 + t781 * t735;
t677 = -t714 * qJD(5) + t781 * t708 - t776 * t709;
t713 = t781 * t734 - t776 * t735;
t691 = -t713 * mrSges(6,1) + t714 * mrSges(6,2);
t744 = qJD(5) + t749;
t700 = t744 * mrSges(6,1) - t714 * mrSges(6,3);
t728 = qJDD(5) + t730;
t692 = -t713 * pkin(5) - t714 * pkin(11);
t743 = t744 ^ 2;
t654 = -t743 * pkin(5) + t728 * pkin(11) + t713 * t692 + t656;
t669 = -t766 * pkin(3) - t765 * pkin(9) + t751 * t727 - t679;
t664 = -t708 * pkin(4) - t733 * pkin(10) + t735 * t719 + t669;
t678 = t713 * qJD(5) + t776 * t708 + t781 * t709;
t657 = (-t713 * t744 - t678) * pkin(11) + (t714 * t744 - t677) * pkin(5) + t664;
t775 = sin(qJ(6));
t780 = cos(qJ(6));
t651 = -t775 * t654 + t780 * t657;
t694 = -t775 * t714 + t780 * t744;
t667 = t694 * qJD(6) + t780 * t678 + t775 * t728;
t676 = qJDD(6) - t677;
t695 = t780 * t714 + t775 * t744;
t684 = -t694 * mrSges(7,1) + t695 * mrSges(7,2);
t712 = qJD(6) - t713;
t685 = -t712 * mrSges(7,2) + t694 * mrSges(7,3);
t649 = m(7) * t651 + t676 * mrSges(7,1) - t667 * mrSges(7,3) - t695 * t684 + t712 * t685;
t652 = t780 * t654 + t775 * t657;
t666 = -t695 * qJD(6) - t775 * t678 + t780 * t728;
t686 = t712 * mrSges(7,1) - t695 * mrSges(7,3);
t650 = m(7) * t652 - t676 * mrSges(7,2) + t666 * mrSges(7,3) + t694 * t684 - t712 * t686;
t793 = -t775 * t649 + t780 * t650;
t640 = m(6) * t656 - t728 * mrSges(6,2) + t677 * mrSges(6,3) + t713 * t691 - t744 * t700 + t793;
t655 = t781 * t659 - t776 * t661;
t699 = -t744 * mrSges(6,2) + t713 * mrSges(6,3);
t653 = -t728 * pkin(5) - t743 * pkin(11) + t714 * t692 - t655;
t789 = -m(7) * t653 + t666 * mrSges(7,1) - t667 * mrSges(7,2) + t694 * t685 - t695 * t686;
t645 = m(6) * t655 + t728 * mrSges(6,1) - t678 * mrSges(6,3) - t714 * t691 + t744 * t699 + t789;
t635 = t776 * t640 + t781 * t645;
t715 = -t734 * mrSges(5,1) + t735 * mrSges(5,2);
t717 = -t749 * mrSges(5,2) + t734 * mrSges(5,3);
t633 = m(5) * t662 + t730 * mrSges(5,1) - t709 * mrSges(5,3) - t735 * t715 + t749 * t717 + t635;
t718 = t749 * mrSges(5,1) - t735 * mrSges(5,3);
t794 = t781 * t640 - t776 * t645;
t634 = m(5) * t663 - t730 * mrSges(5,2) + t708 * mrSges(5,3) + t734 * t715 - t749 * t718 + t794;
t795 = -t777 * t633 + t782 * t634;
t624 = m(4) * t680 - t766 * mrSges(4,2) + t731 * mrSges(4,3) + t750 * t726 - t767 * t737 + t795;
t736 = -t767 * mrSges(4,2) + t750 * mrSges(4,3);
t641 = t780 * t649 + t775 * t650;
t787 = m(6) * t664 - t677 * mrSges(6,1) + t678 * mrSges(6,2) - t713 * t699 + t714 * t700 + t641;
t786 = -m(5) * t669 + t708 * mrSges(5,1) - t709 * mrSges(5,2) + t734 * t717 - t735 * t718 - t787;
t637 = m(4) * t679 + t766 * mrSges(4,1) - t732 * mrSges(4,3) - t751 * t726 + t767 * t736 + t786;
t621 = t771 * t624 + t773 * t637;
t724 = -g(3) * t807 + t791;
t756 = -t767 * mrSges(3,2) + mrSges(3,3) * t798;
t759 = (-mrSges(3,1) * t783 + mrSges(3,2) * t778) * t803;
t619 = m(3) * t724 + t766 * mrSges(3,1) - t760 * mrSges(3,3) + t767 * t756 - t759 * t799 + t621;
t755 = t767 * mrSges(3,1) - mrSges(3,3) * t799;
t796 = t773 * t624 - t771 * t637;
t620 = m(3) * t725 - t766 * mrSges(3,2) + t761 * mrSges(3,3) - t767 * t755 + t759 * t798 + t796;
t627 = t782 * t633 + t777 * t634;
t788 = m(4) * t711 - t731 * mrSges(4,1) + t732 * mrSges(4,2) - t750 * t736 + t751 * t737 + t627;
t626 = m(3) * t741 - t761 * mrSges(3,1) + t760 * mrSges(3,2) + (t755 * t778 - t756 * t783) * t803 + t788;
t606 = t619 * t805 + t620 * t806 - t772 * t626;
t604 = m(2) * t763 + qJDD(1) * mrSges(2,1) - t785 * mrSges(2,2) + t606;
t610 = -t778 * t619 + t783 * t620;
t609 = m(2) * t764 - t785 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t610;
t804 = t784 * t604 + t779 * t609;
t605 = t619 * t807 + t620 * t808 + t774 * t626;
t797 = -t779 * t604 + t784 * t609;
t671 = Ifges(7,5) * t695 + Ifges(7,6) * t694 + Ifges(7,3) * t712;
t673 = Ifges(7,1) * t695 + Ifges(7,4) * t694 + Ifges(7,5) * t712;
t642 = -mrSges(7,1) * t653 + mrSges(7,3) * t652 + Ifges(7,4) * t667 + Ifges(7,2) * t666 + Ifges(7,6) * t676 - t695 * t671 + t712 * t673;
t672 = Ifges(7,4) * t695 + Ifges(7,2) * t694 + Ifges(7,6) * t712;
t643 = mrSges(7,2) * t653 - mrSges(7,3) * t651 + Ifges(7,1) * t667 + Ifges(7,4) * t666 + Ifges(7,5) * t676 + t694 * t671 - t712 * t672;
t687 = Ifges(6,5) * t714 + Ifges(6,6) * t713 + Ifges(6,3) * t744;
t688 = Ifges(6,4) * t714 + Ifges(6,2) * t713 + Ifges(6,6) * t744;
t628 = mrSges(6,2) * t664 - mrSges(6,3) * t655 + Ifges(6,1) * t678 + Ifges(6,4) * t677 + Ifges(6,5) * t728 - pkin(11) * t641 - t775 * t642 + t780 * t643 + t713 * t687 - t744 * t688;
t689 = Ifges(6,1) * t714 + Ifges(6,4) * t713 + Ifges(6,5) * t744;
t629 = -mrSges(6,1) * t664 - mrSges(7,1) * t651 + mrSges(7,2) * t652 + mrSges(6,3) * t656 + Ifges(6,4) * t678 - Ifges(7,5) * t667 + Ifges(6,2) * t677 + Ifges(6,6) * t728 - Ifges(7,6) * t666 - Ifges(7,3) * t676 - pkin(5) * t641 - t695 * t672 + t694 * t673 - t714 * t687 + t744 * t689;
t702 = Ifges(5,5) * t735 + Ifges(5,6) * t734 + Ifges(5,3) * t749;
t704 = Ifges(5,1) * t735 + Ifges(5,4) * t734 + Ifges(5,5) * t749;
t612 = -mrSges(5,1) * t669 + mrSges(5,3) * t663 + Ifges(5,4) * t709 + Ifges(5,2) * t708 + Ifges(5,6) * t730 - pkin(4) * t787 + pkin(10) * t794 + t776 * t628 + t781 * t629 - t735 * t702 + t749 * t704;
t703 = Ifges(5,4) * t735 + Ifges(5,2) * t734 + Ifges(5,6) * t749;
t613 = mrSges(5,2) * t669 - mrSges(5,3) * t662 + Ifges(5,1) * t709 + Ifges(5,4) * t708 + Ifges(5,5) * t730 - pkin(10) * t635 + t781 * t628 - t776 * t629 + t734 * t702 - t749 * t703;
t720 = Ifges(4,5) * t751 + Ifges(4,6) * t750 + Ifges(4,3) * t767;
t721 = Ifges(4,4) * t751 + Ifges(4,2) * t750 + Ifges(4,6) * t767;
t602 = mrSges(4,2) * t711 - mrSges(4,3) * t679 + Ifges(4,1) * t732 + Ifges(4,4) * t731 + Ifges(4,5) * t766 - pkin(9) * t627 - t777 * t612 + t782 * t613 + t750 * t720 - t767 * t721;
t722 = Ifges(4,1) * t751 + Ifges(4,4) * t750 + Ifges(4,5) * t767;
t611 = -pkin(5) * t789 - pkin(11) * t793 - t780 * t642 - t775 * t643 + Ifges(4,6) * t766 + t767 * t722 - t751 * t720 + t734 * t704 - t735 * t703 + Ifges(4,2) * t731 + Ifges(4,4) * t732 - Ifges(6,3) * t728 - Ifges(5,3) * t730 - t714 * t688 + t713 * t689 - mrSges(4,1) * t711 - Ifges(5,6) * t708 - Ifges(5,5) * t709 + mrSges(4,3) * t680 - Ifges(6,6) * t677 - Ifges(6,5) * t678 - mrSges(5,1) * t662 + mrSges(5,2) * t663 + mrSges(6,2) * t656 - mrSges(6,1) * t655 - pkin(4) * t635 - pkin(3) * t627;
t738 = Ifges(3,3) * t767 + (Ifges(3,5) * t778 + Ifges(3,6) * t783) * t803;
t740 = Ifges(3,5) * t767 + (Ifges(3,1) * t778 + Ifges(3,4) * t783) * t803;
t599 = -mrSges(3,1) * t741 + mrSges(3,3) * t725 + Ifges(3,4) * t760 + Ifges(3,2) * t761 + Ifges(3,6) * t766 - pkin(2) * t788 + qJ(3) * t796 + t771 * t602 + t773 * t611 - t738 * t799 + t767 * t740;
t739 = Ifges(3,6) * t767 + (Ifges(3,4) * t778 + Ifges(3,2) * t783) * t803;
t600 = mrSges(3,2) * t741 - mrSges(3,3) * t724 + Ifges(3,1) * t760 + Ifges(3,4) * t761 + Ifges(3,5) * t766 - qJ(3) * t621 + t773 * t602 - t771 * t611 + t738 * t798 - t767 * t739;
t790 = pkin(8) * t610 + t599 * t783 + t600 * t778;
t601 = Ifges(3,5) * t760 + Ifges(3,6) * t761 + mrSges(3,1) * t724 - mrSges(3,2) * t725 + Ifges(4,5) * t732 + Ifges(4,6) * t731 + t751 * t721 - t750 * t722 + mrSges(4,1) * t679 - mrSges(4,2) * t680 + t777 * t613 + t782 * t612 + pkin(3) * t786 + pkin(9) * t795 + pkin(2) * t621 + (Ifges(3,3) + Ifges(4,3)) * t766 + (t739 * t778 - t740 * t783) * t803;
t598 = -mrSges(2,2) * g(3) - mrSges(2,3) * t763 + Ifges(2,5) * qJDD(1) - t785 * Ifges(2,6) - t778 * t599 + t783 * t600 + (-t605 * t772 - t606 * t774) * pkin(8);
t597 = mrSges(2,1) * g(3) + mrSges(2,3) * t764 + t785 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t605 - t772 * t601 + t790 * t774;
t1 = [-m(1) * g(1) + t797; -m(1) * g(2) + t804; (-m(1) - m(2)) * g(3) + t605; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t804 - t779 * t597 + t784 * t598; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t797 + t784 * t597 + t779 * t598; -mrSges(1,1) * g(2) + mrSges(2,1) * t763 + mrSges(1,2) * g(1) - mrSges(2,2) * t764 + Ifges(2,3) * qJDD(1) + pkin(1) * t606 + t774 * t601 + t790 * t772;];
tauB  = t1;
