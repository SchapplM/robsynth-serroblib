% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 12:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 12:21:01
% EndTime: 2019-05-08 12:22:35
% DurationCPUTime: 67.78s
% Computational Cost: add. (1150183->399), mult. (2450788->516), div. (0->0), fcn. (1995776->14), ass. (0->165)
t762 = sin(pkin(6));
t768 = sin(qJ(2));
t774 = cos(qJ(2));
t790 = qJD(1) * qJD(2);
t751 = (-qJDD(1) * t774 + t768 * t790) * t762;
t800 = pkin(8) * t762;
t763 = cos(pkin(6));
t799 = t763 * g(3);
t798 = t762 * t768;
t797 = t762 * t774;
t796 = t763 * t768;
t795 = t763 * t774;
t769 = sin(qJ(1));
t775 = cos(qJ(1));
t755 = t769 * g(1) - t775 * g(2);
t776 = qJD(1) ^ 2;
t746 = qJDD(1) * pkin(1) + t776 * t800 + t755;
t756 = -t775 * g(1) - t769 * g(2);
t747 = -t776 * pkin(1) + qJDD(1) * t800 + t756;
t793 = t746 * t796 + t774 * t747;
t720 = -g(3) * t798 + t793;
t759 = t763 * qJD(1) + qJD(2);
t792 = qJD(1) * t762;
t789 = t768 * t792;
t744 = t759 * mrSges(3,1) - mrSges(3,3) * t789;
t748 = (-mrSges(3,1) * t774 + mrSges(3,2) * t768) * t792;
t758 = t763 * qJDD(1) + qJDD(2);
t749 = (-pkin(2) * t774 - pkin(9) * t768) * t792;
t757 = t759 ^ 2;
t791 = qJD(1) * t774;
t699 = -t757 * pkin(2) + t758 * pkin(9) + (-g(3) * t768 + t749 * t791) * t762 + t793;
t750 = (qJDD(1) * t768 + t774 * t790) * t762;
t700 = t751 * pkin(2) - t750 * pkin(9) - t799 + (-t746 + (pkin(2) * t768 - pkin(9) * t774) * t759 * qJD(1)) * t762;
t767 = sin(qJ(3));
t773 = cos(qJ(3));
t676 = t773 * t699 + t767 * t700;
t739 = t767 * t759 + t773 * t789;
t717 = -t739 * qJD(3) - t767 * t750 + t773 * t758;
t738 = t773 * t759 - t767 * t789;
t721 = -t738 * mrSges(4,1) + t739 * mrSges(4,2);
t788 = t762 * t791;
t754 = qJD(3) - t788;
t727 = t754 * mrSges(4,1) - t739 * mrSges(4,3);
t743 = qJDD(3) + t751;
t722 = -t738 * pkin(3) - t739 * pkin(10);
t752 = t754 ^ 2;
t669 = -t752 * pkin(3) + t743 * pkin(10) + t738 * t722 + t676;
t719 = -g(3) * t797 + t746 * t795 - t768 * t747;
t698 = -t758 * pkin(2) - t757 * pkin(9) + t749 * t789 - t719;
t718 = t738 * qJD(3) + t773 * t750 + t767 * t758;
t672 = (-t738 * t754 - t718) * pkin(10) + (t739 * t754 - t717) * pkin(3) + t698;
t766 = sin(qJ(4));
t772 = cos(qJ(4));
t655 = -t766 * t669 + t772 * t672;
t724 = -t766 * t739 + t772 * t754;
t686 = t724 * qJD(4) + t772 * t718 + t766 * t743;
t715 = qJDD(4) - t717;
t725 = t772 * t739 + t766 * t754;
t737 = qJD(4) - t738;
t648 = (t724 * t737 - t686) * pkin(11) + (t724 * t725 + t715) * pkin(4) + t655;
t656 = t772 * t669 + t766 * t672;
t685 = -t725 * qJD(4) - t766 * t718 + t772 * t743;
t708 = t737 * pkin(4) - t725 * pkin(11);
t723 = t724 ^ 2;
t654 = -t723 * pkin(4) + t685 * pkin(11) - t737 * t708 + t656;
t765 = sin(qJ(5));
t771 = cos(qJ(5));
t642 = t771 * t648 - t765 * t654;
t702 = t771 * t724 - t765 * t725;
t666 = t702 * qJD(5) + t765 * t685 + t771 * t686;
t703 = t765 * t724 + t771 * t725;
t710 = qJDD(5) + t715;
t735 = qJD(5) + t737;
t640 = (t702 * t735 - t666) * pkin(12) + (t702 * t703 + t710) * pkin(5) + t642;
t643 = t765 * t648 + t771 * t654;
t665 = -t703 * qJD(5) + t771 * t685 - t765 * t686;
t689 = t735 * pkin(5) - t703 * pkin(12);
t701 = t702 ^ 2;
t641 = -t701 * pkin(5) + t665 * pkin(12) - t735 * t689 + t643;
t764 = sin(qJ(6));
t770 = cos(qJ(6));
t638 = t770 * t640 - t764 * t641;
t681 = t770 * t702 - t764 * t703;
t652 = t681 * qJD(6) + t764 * t665 + t770 * t666;
t682 = t764 * t702 + t770 * t703;
t663 = -t681 * mrSges(7,1) + t682 * mrSges(7,2);
t732 = qJD(6) + t735;
t673 = -t732 * mrSges(7,2) + t681 * mrSges(7,3);
t709 = qJDD(6) + t710;
t636 = m(7) * t638 + t709 * mrSges(7,1) - t652 * mrSges(7,3) - t682 * t663 + t732 * t673;
t639 = t764 * t640 + t770 * t641;
t651 = -t682 * qJD(6) + t770 * t665 - t764 * t666;
t674 = t732 * mrSges(7,1) - t682 * mrSges(7,3);
t637 = m(7) * t639 - t709 * mrSges(7,2) + t651 * mrSges(7,3) + t681 * t663 - t732 * t674;
t628 = t770 * t636 + t764 * t637;
t683 = -t702 * mrSges(6,1) + t703 * mrSges(6,2);
t687 = -t735 * mrSges(6,2) + t702 * mrSges(6,3);
t626 = m(6) * t642 + t710 * mrSges(6,1) - t666 * mrSges(6,3) - t703 * t683 + t735 * t687 + t628;
t688 = t735 * mrSges(6,1) - t703 * mrSges(6,3);
t783 = -t764 * t636 + t770 * t637;
t627 = m(6) * t643 - t710 * mrSges(6,2) + t665 * mrSges(6,3) + t702 * t683 - t735 * t688 + t783;
t622 = t771 * t626 + t765 * t627;
t704 = -t724 * mrSges(5,1) + t725 * mrSges(5,2);
t706 = -t737 * mrSges(5,2) + t724 * mrSges(5,3);
t620 = m(5) * t655 + t715 * mrSges(5,1) - t686 * mrSges(5,3) - t725 * t704 + t737 * t706 + t622;
t707 = t737 * mrSges(5,1) - t725 * mrSges(5,3);
t784 = -t765 * t626 + t771 * t627;
t621 = m(5) * t656 - t715 * mrSges(5,2) + t685 * mrSges(5,3) + t724 * t704 - t737 * t707 + t784;
t785 = -t766 * t620 + t772 * t621;
t615 = m(4) * t676 - t743 * mrSges(4,2) + t717 * mrSges(4,3) + t738 * t721 - t754 * t727 + t785;
t675 = -t767 * t699 + t773 * t700;
t726 = -t754 * mrSges(4,2) + t738 * mrSges(4,3);
t668 = -t743 * pkin(3) - t752 * pkin(10) + t739 * t722 - t675;
t657 = -t685 * pkin(4) - t723 * pkin(11) + t725 * t708 + t668;
t645 = -t665 * pkin(5) - t701 * pkin(12) + t703 * t689 + t657;
t782 = m(7) * t645 - t651 * mrSges(7,1) + t652 * mrSges(7,2) - t681 * t673 + t682 * t674;
t779 = m(6) * t657 - t665 * mrSges(6,1) + t666 * mrSges(6,2) - t702 * t687 + t703 * t688 + t782;
t777 = -m(5) * t668 + t685 * mrSges(5,1) - t686 * mrSges(5,2) + t724 * t706 - t725 * t707 - t779;
t632 = m(4) * t675 + t743 * mrSges(4,1) - t718 * mrSges(4,3) - t739 * t721 + t754 * t726 + t777;
t786 = t773 * t615 - t767 * t632;
t606 = m(3) * t720 - t758 * mrSges(3,2) - t751 * mrSges(3,3) - t759 * t744 + t748 * t788 + t786;
t609 = t767 * t615 + t773 * t632;
t731 = -t762 * t746 - t799;
t745 = -t759 * mrSges(3,2) + mrSges(3,3) * t788;
t608 = m(3) * t731 + t751 * mrSges(3,1) + t750 * mrSges(3,2) + (t744 * t768 - t745 * t774) * t792 + t609;
t616 = t772 * t620 + t766 * t621;
t778 = -m(4) * t698 + t717 * mrSges(4,1) - t718 * mrSges(4,2) + t738 * t726 - t739 * t727 - t616;
t612 = m(3) * t719 + t758 * mrSges(3,1) - t750 * mrSges(3,3) + t759 * t745 - t748 * t789 + t778;
t595 = t606 * t796 - t762 * t608 + t612 * t795;
t593 = m(2) * t755 + qJDD(1) * mrSges(2,1) - t776 * mrSges(2,2) + t595;
t599 = t774 * t606 - t768 * t612;
t598 = m(2) * t756 - t776 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t599;
t794 = t775 * t593 + t769 * t598;
t594 = t606 * t798 + t763 * t608 + t612 * t797;
t787 = -t769 * t593 + t775 * t598;
t658 = Ifges(7,5) * t682 + Ifges(7,6) * t681 + Ifges(7,3) * t732;
t660 = Ifges(7,1) * t682 + Ifges(7,4) * t681 + Ifges(7,5) * t732;
t629 = -mrSges(7,1) * t645 + mrSges(7,3) * t639 + Ifges(7,4) * t652 + Ifges(7,2) * t651 + Ifges(7,6) * t709 - t682 * t658 + t732 * t660;
t659 = Ifges(7,4) * t682 + Ifges(7,2) * t681 + Ifges(7,6) * t732;
t630 = mrSges(7,2) * t645 - mrSges(7,3) * t638 + Ifges(7,1) * t652 + Ifges(7,4) * t651 + Ifges(7,5) * t709 + t681 * t658 - t732 * t659;
t677 = Ifges(6,5) * t703 + Ifges(6,6) * t702 + Ifges(6,3) * t735;
t679 = Ifges(6,1) * t703 + Ifges(6,4) * t702 + Ifges(6,5) * t735;
t617 = -mrSges(6,1) * t657 + mrSges(6,3) * t643 + Ifges(6,4) * t666 + Ifges(6,2) * t665 + Ifges(6,6) * t710 - pkin(5) * t782 + pkin(12) * t783 + t770 * t629 + t764 * t630 - t703 * t677 + t735 * t679;
t678 = Ifges(6,4) * t703 + Ifges(6,2) * t702 + Ifges(6,6) * t735;
t618 = mrSges(6,2) * t657 - mrSges(6,3) * t642 + Ifges(6,1) * t666 + Ifges(6,4) * t665 + Ifges(6,5) * t710 - pkin(12) * t628 - t764 * t629 + t770 * t630 + t702 * t677 - t735 * t678;
t690 = Ifges(5,5) * t725 + Ifges(5,6) * t724 + Ifges(5,3) * t737;
t692 = Ifges(5,1) * t725 + Ifges(5,4) * t724 + Ifges(5,5) * t737;
t601 = -mrSges(5,1) * t668 + mrSges(5,3) * t656 + Ifges(5,4) * t686 + Ifges(5,2) * t685 + Ifges(5,6) * t715 - pkin(4) * t779 + pkin(11) * t784 + t771 * t617 + t765 * t618 - t725 * t690 + t737 * t692;
t691 = Ifges(5,4) * t725 + Ifges(5,2) * t724 + Ifges(5,6) * t737;
t602 = mrSges(5,2) * t668 - mrSges(5,3) * t655 + Ifges(5,1) * t686 + Ifges(5,4) * t685 + Ifges(5,5) * t715 - pkin(11) * t622 - t765 * t617 + t771 * t618 + t724 * t690 - t737 * t691;
t711 = Ifges(4,5) * t739 + Ifges(4,6) * t738 + Ifges(4,3) * t754;
t712 = Ifges(4,4) * t739 + Ifges(4,2) * t738 + Ifges(4,6) * t754;
t591 = mrSges(4,2) * t698 - mrSges(4,3) * t675 + Ifges(4,1) * t718 + Ifges(4,4) * t717 + Ifges(4,5) * t743 - pkin(10) * t616 - t766 * t601 + t772 * t602 + t738 * t711 - t754 * t712;
t713 = Ifges(4,1) * t739 + Ifges(4,4) * t738 + Ifges(4,5) * t754;
t600 = t754 * t713 + Ifges(4,6) * t743 - t739 * t711 - t725 * t691 + Ifges(4,4) * t718 + t724 * t692 - Ifges(5,3) * t715 + Ifges(4,2) * t717 - Ifges(7,3) * t709 - Ifges(6,3) * t710 - t703 * t678 + t702 * t679 - mrSges(4,1) * t698 - Ifges(5,6) * t685 - Ifges(5,5) * t686 + t681 * t660 - t682 * t659 + mrSges(4,3) * t676 - Ifges(6,6) * t665 - Ifges(6,5) * t666 - mrSges(5,1) * t655 + mrSges(5,2) * t656 - Ifges(7,6) * t651 - Ifges(7,5) * t652 + mrSges(6,2) * t643 - mrSges(6,1) * t642 + mrSges(7,2) * t639 - mrSges(7,1) * t638 - pkin(3) * t616 - pkin(5) * t628 - pkin(4) * t622;
t728 = Ifges(3,3) * t759 + (Ifges(3,5) * t768 + Ifges(3,6) * t774) * t792;
t729 = Ifges(3,6) * t759 + (Ifges(3,4) * t768 + Ifges(3,2) * t774) * t792;
t589 = mrSges(3,2) * t731 - mrSges(3,3) * t719 + Ifges(3,1) * t750 - Ifges(3,4) * t751 + Ifges(3,5) * t758 - pkin(9) * t609 + t773 * t591 - t767 * t600 + t728 * t788 - t759 * t729;
t730 = Ifges(3,5) * t759 + (Ifges(3,1) * t768 + Ifges(3,4) * t774) * t792;
t590 = Ifges(3,4) * t750 - Ifges(3,2) * t751 + Ifges(3,6) * t758 - t728 * t789 + t759 * t730 - mrSges(3,1) * t731 + mrSges(3,3) * t720 - Ifges(4,5) * t718 - Ifges(4,6) * t717 - Ifges(4,3) * t743 - t739 * t712 + t738 * t713 - mrSges(4,1) * t675 + mrSges(4,2) * t676 - t766 * t602 - t772 * t601 - pkin(3) * t777 - pkin(10) * t785 - pkin(2) * t609;
t780 = pkin(8) * t599 + t589 * t768 + t590 * t774;
t588 = Ifges(3,5) * t750 - Ifges(3,6) * t751 + Ifges(3,3) * t758 + mrSges(3,1) * t719 - mrSges(3,2) * t720 + t767 * t591 + t773 * t600 + pkin(2) * t778 + pkin(9) * t786 + (t729 * t768 - t730 * t774) * t792;
t587 = -mrSges(2,2) * g(3) - mrSges(2,3) * t755 + Ifges(2,5) * qJDD(1) - t776 * Ifges(2,6) + t774 * t589 - t768 * t590 + (-t594 * t762 - t595 * t763) * pkin(8);
t586 = mrSges(2,1) * g(3) + mrSges(2,3) * t756 + t776 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t594 - t762 * t588 + t780 * t763;
t1 = [-m(1) * g(1) + t787; -m(1) * g(2) + t794; (-m(1) - m(2)) * g(3) + t594; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t794 - t769 * t586 + t775 * t587; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t787 + t775 * t586 + t769 * t587; -mrSges(1,1) * g(2) + mrSges(2,1) * t755 + mrSges(1,2) * g(1) - mrSges(2,2) * t756 + Ifges(2,3) * qJDD(1) + pkin(1) * t595 + t763 * t588 + t780 * t762;];
tauB  = t1;
