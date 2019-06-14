% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 16:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR14_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 16:34:56
% EndTime: 2019-05-07 16:35:24
% DurationCPUTime: 20.44s
% Computational Cost: add. (324602->380), mult. (694993->473), div. (0->0), fcn. (541964->12), ass. (0->157)
t823 = Ifges(4,1) + Ifges(5,2);
t817 = Ifges(4,5) - Ifges(5,4);
t822 = Ifges(4,2) + Ifges(5,3);
t816 = Ifges(4,6) - Ifges(5,5);
t815 = -Ifges(5,6) - Ifges(4,4);
t821 = Ifges(4,3) + Ifges(5,1);
t773 = sin(pkin(6));
t778 = sin(qJ(2));
t782 = cos(qJ(2));
t800 = qJD(1) * qJD(2);
t759 = (-qJDD(1) * t782 + t778 * t800) * t773;
t803 = qJD(1) * t773;
t757 = (-pkin(2) * t782 - pkin(9) * t778) * t803;
t774 = cos(pkin(6));
t770 = t774 * qJD(1) + qJD(2);
t768 = t770 ^ 2;
t769 = t774 * qJDD(1) + qJDD(2);
t802 = qJD(1) * t782;
t779 = sin(qJ(1));
t783 = cos(qJ(1));
t766 = t779 * g(1) - t783 * g(2);
t784 = qJD(1) ^ 2;
t819 = pkin(8) * t773;
t754 = qJDD(1) * pkin(1) + t784 * t819 + t766;
t767 = -t783 * g(1) - t779 * g(2);
t755 = -t784 * pkin(1) + qJDD(1) * t819 + t767;
t811 = t774 * t778;
t804 = t754 * t811 + t782 * t755;
t693 = -t768 * pkin(2) + t769 * pkin(9) + (-g(3) * t778 + t757 * t802) * t773 + t804;
t758 = (qJDD(1) * t778 + t782 * t800) * t773;
t818 = t774 * g(3);
t694 = t759 * pkin(2) - t758 * pkin(9) - t818 + (-t754 + (pkin(2) * t778 - pkin(9) * t782) * t770 * qJD(1)) * t773;
t777 = sin(qJ(3));
t820 = cos(qJ(3));
t674 = t820 * t693 + t777 * t694;
t799 = t778 * t803;
t745 = -t820 * t770 + t777 * t799;
t746 = t777 * t770 + t820 * t799;
t721 = t745 * pkin(3) - t746 * qJ(4);
t751 = qJDD(3) + t759;
t798 = t773 * t802;
t764 = -qJD(3) + t798;
t763 = t764 ^ 2;
t670 = t763 * pkin(3) - t751 * qJ(4) + 0.2e1 * qJD(4) * t764 + t745 * t721 - t674;
t814 = t745 * t764;
t813 = t773 * t778;
t812 = t773 * t782;
t810 = t774 * t782;
t720 = -g(3) * t813 + t804;
t752 = t770 * mrSges(3,1) - mrSges(3,3) * t799;
t756 = (-mrSges(3,1) * t782 + mrSges(3,2) * t778) * t803;
t673 = -t777 * t693 + t820 * t694;
t718 = -t745 * qJD(3) + t820 * t758 + t777 * t769;
t722 = t745 * mrSges(4,1) + t746 * mrSges(4,2);
t729 = t764 * mrSges(4,2) - t745 * mrSges(4,3);
t731 = t745 * mrSges(5,1) + t764 * mrSges(5,3);
t671 = -t751 * pkin(3) - t763 * qJ(4) + t746 * t721 + qJDD(4) - t673;
t662 = (t745 * t746 - t751) * pkin(10) + (t718 - t814) * pkin(4) + t671;
t717 = t746 * qJD(3) + t777 * t758 - t820 * t769;
t733 = t746 * pkin(4) + t764 * pkin(10);
t744 = t745 ^ 2;
t719 = -g(3) * t812 + t754 * t810 - t778 * t755;
t692 = -t769 * pkin(2) - t768 * pkin(9) + t757 * t799 - t719;
t786 = (-t718 - t814) * qJ(4) + t692 + (-t764 * pkin(3) - 0.2e1 * qJD(4)) * t746;
t666 = -t744 * pkin(4) - t746 * t733 + (pkin(3) + pkin(10)) * t717 + t786;
t776 = sin(qJ(5));
t781 = cos(qJ(5));
t656 = t781 * t662 - t776 * t666;
t727 = t781 * t745 + t776 * t764;
t682 = t727 * qJD(5) + t776 * t717 + t781 * t751;
t714 = qJDD(5) + t718;
t728 = t776 * t745 - t781 * t764;
t743 = qJD(5) + t746;
t654 = (t727 * t743 - t682) * pkin(11) + (t727 * t728 + t714) * pkin(5) + t656;
t657 = t776 * t662 + t781 * t666;
t681 = -t728 * qJD(5) + t781 * t717 - t776 * t751;
t703 = t743 * pkin(5) - t728 * pkin(11);
t726 = t727 ^ 2;
t655 = -t726 * pkin(5) + t681 * pkin(11) - t743 * t703 + t657;
t775 = sin(qJ(6));
t780 = cos(qJ(6));
t652 = t780 * t654 - t775 * t655;
t695 = t780 * t727 - t775 * t728;
t669 = t695 * qJD(6) + t775 * t681 + t780 * t682;
t696 = t775 * t727 + t780 * t728;
t679 = -t695 * mrSges(7,1) + t696 * mrSges(7,2);
t741 = qJD(6) + t743;
t683 = -t741 * mrSges(7,2) + t695 * mrSges(7,3);
t704 = qJDD(6) + t714;
t650 = m(7) * t652 + t704 * mrSges(7,1) - t669 * mrSges(7,3) - t696 * t679 + t741 * t683;
t653 = t775 * t654 + t780 * t655;
t668 = -t696 * qJD(6) + t780 * t681 - t775 * t682;
t684 = t741 * mrSges(7,1) - t696 * mrSges(7,3);
t651 = m(7) * t653 - t704 * mrSges(7,2) + t668 * mrSges(7,3) + t695 * t679 - t741 * t684;
t641 = t780 * t650 + t775 * t651;
t697 = -t727 * mrSges(6,1) + t728 * mrSges(6,2);
t701 = -t743 * mrSges(6,2) + t727 * mrSges(6,3);
t639 = m(6) * t656 + t714 * mrSges(6,1) - t682 * mrSges(6,3) - t728 * t697 + t743 * t701 + t641;
t702 = t743 * mrSges(6,1) - t728 * mrSges(6,3);
t795 = -t775 * t650 + t780 * t651;
t640 = m(6) * t657 - t714 * mrSges(6,2) + t681 * mrSges(6,3) + t727 * t697 - t743 * t702 + t795;
t636 = t781 * t639 + t776 * t640;
t723 = -t745 * mrSges(5,2) - t746 * mrSges(5,3);
t789 = -m(5) * t671 - t718 * mrSges(5,1) - t746 * t723 - t636;
t634 = m(4) * t673 - t718 * mrSges(4,3) - t746 * t722 + (-t729 + t731) * t764 + (mrSges(4,1) - mrSges(5,2)) * t751 + t789;
t730 = -t764 * mrSges(4,1) - t746 * mrSges(4,3);
t732 = t746 * mrSges(5,1) - t764 * mrSges(5,2);
t665 = -t717 * pkin(4) - t744 * pkin(10) - t764 * t733 - t670;
t659 = -t681 * pkin(5) - t726 * pkin(11) + t728 * t703 + t665;
t790 = m(7) * t659 - t668 * mrSges(7,1) + t669 * mrSges(7,2) - t695 * t683 + t696 * t684;
t788 = -m(6) * t665 + t681 * mrSges(6,1) - t682 * mrSges(6,2) + t727 * t701 - t728 * t702 - t790;
t785 = -m(5) * t670 + t751 * mrSges(5,3) - t764 * t732 - t788;
t646 = m(4) * t674 + (-mrSges(4,3) - mrSges(5,1)) * t717 + (-t722 - t723) * t745 + t785 + t764 * t730 - t751 * mrSges(4,2);
t796 = -t777 * t634 + t820 * t646;
t624 = m(3) * t720 - t769 * mrSges(3,2) - t759 * mrSges(3,3) - t770 * t752 + t756 * t798 + t796;
t627 = t820 * t634 + t777 * t646;
t738 = -t773 * t754 - t818;
t753 = -t770 * mrSges(3,2) + mrSges(3,3) * t798;
t626 = m(3) * t738 + t759 * mrSges(3,1) + t758 * mrSges(3,2) + (t752 * t778 - t753 * t782) * t803 + t627;
t672 = t717 * pkin(3) + t786;
t808 = -t776 * t639 + t781 * t640;
t794 = -m(5) * t672 + t717 * mrSges(5,2) + t745 * t731 - t808;
t787 = -m(4) * t692 - t717 * mrSges(4,1) - t745 * t729 + (-t730 + t732) * t746 + (-mrSges(4,2) + mrSges(5,3)) * t718 + t794;
t632 = m(3) * t719 + t769 * mrSges(3,1) - t758 * mrSges(3,3) + t770 * t753 - t756 * t799 + t787;
t615 = t624 * t811 - t773 * t626 + t632 * t810;
t613 = m(2) * t766 + qJDD(1) * mrSges(2,1) - t784 * mrSges(2,2) + t615;
t620 = t782 * t624 - t778 * t632;
t619 = m(2) * t767 - t784 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t620;
t809 = t783 * t613 + t779 * t619;
t807 = t745 * t816 - t746 * t817 + t764 * t821;
t806 = t745 * t822 + t746 * t815 + t764 * t816;
t805 = -t815 * t745 - t746 * t823 + t817 * t764;
t614 = t624 * t813 + t774 * t626 + t632 * t812;
t797 = -t779 * t613 + t783 * t619;
t675 = Ifges(7,5) * t696 + Ifges(7,6) * t695 + Ifges(7,3) * t741;
t677 = Ifges(7,1) * t696 + Ifges(7,4) * t695 + Ifges(7,5) * t741;
t642 = -mrSges(7,1) * t659 + mrSges(7,3) * t653 + Ifges(7,4) * t669 + Ifges(7,2) * t668 + Ifges(7,6) * t704 - t696 * t675 + t741 * t677;
t676 = Ifges(7,4) * t696 + Ifges(7,2) * t695 + Ifges(7,6) * t741;
t643 = mrSges(7,2) * t659 - mrSges(7,3) * t652 + Ifges(7,1) * t669 + Ifges(7,4) * t668 + Ifges(7,5) * t704 + t695 * t675 - t741 * t676;
t685 = Ifges(6,5) * t728 + Ifges(6,6) * t727 + Ifges(6,3) * t743;
t687 = Ifges(6,1) * t728 + Ifges(6,4) * t727 + Ifges(6,5) * t743;
t628 = -mrSges(6,1) * t665 + mrSges(6,3) * t657 + Ifges(6,4) * t682 + Ifges(6,2) * t681 + Ifges(6,6) * t714 - pkin(5) * t790 + pkin(11) * t795 + t780 * t642 + t775 * t643 - t728 * t685 + t743 * t687;
t686 = Ifges(6,4) * t728 + Ifges(6,2) * t727 + Ifges(6,6) * t743;
t629 = mrSges(6,2) * t665 - mrSges(6,3) * t656 + Ifges(6,1) * t682 + Ifges(6,4) * t681 + Ifges(6,5) * t714 - pkin(11) * t641 - t775 * t642 + t780 * t643 + t727 * t685 - t743 * t686;
t635 = -t718 * mrSges(5,3) - t746 * t732 - t794;
t611 = -mrSges(4,1) * t692 - mrSges(5,1) * t670 + mrSges(5,2) * t672 + mrSges(4,3) * t674 - pkin(3) * t635 - pkin(4) * t788 - pkin(10) * t808 - t781 * t628 - t776 * t629 - t717 * t822 - t815 * t718 + t807 * t746 + t816 * t751 + t805 * t764;
t616 = t823 * t718 + t728 * t686 - t727 * t687 + Ifges(6,3) * t714 + Ifges(7,3) * t704 + mrSges(4,2) * t692 - t695 * t677 + t696 * t676 + Ifges(6,6) * t681 + Ifges(6,5) * t682 + Ifges(7,5) * t669 + mrSges(5,1) * t671 - mrSges(5,3) * t672 - mrSges(4,3) * t673 + Ifges(7,6) * t668 - mrSges(6,2) * t657 + mrSges(6,1) * t656 - mrSges(7,2) * t653 + mrSges(7,1) * t652 + pkin(5) * t641 + pkin(4) * t636 - qJ(4) * t635 + t815 * t717 + t817 * t751 - t806 * t764 + t807 * t745;
t735 = Ifges(3,3) * t770 + (Ifges(3,5) * t778 + Ifges(3,6) * t782) * t803;
t736 = Ifges(3,6) * t770 + (Ifges(3,4) * t778 + Ifges(3,2) * t782) * t803;
t609 = mrSges(3,2) * t738 - mrSges(3,3) * t719 + Ifges(3,1) * t758 - Ifges(3,4) * t759 + Ifges(3,5) * t769 - pkin(9) * t627 - t777 * t611 + t820 * t616 + t735 * t798 - t770 * t736;
t737 = Ifges(3,5) * t770 + (Ifges(3,1) * t778 + Ifges(3,4) * t782) * t803;
t610 = -mrSges(3,1) * t738 + mrSges(3,3) * t720 + mrSges(4,2) * t674 + mrSges(5,3) * t670 - mrSges(5,2) * t671 - mrSges(4,1) * t673 + pkin(10) * t636 - qJ(4) * t785 - pkin(2) * t627 - t781 * t629 - pkin(3) * (t764 * t731 + t789) + t776 * t628 + t770 * t737 + Ifges(3,6) * t769 + Ifges(3,4) * t758 - Ifges(3,2) * t759 - t735 * t799 + (pkin(3) * mrSges(5,2) - t821) * t751 + t806 * t746 + (qJ(4) * t723 + t805) * t745 - t817 * t718 + (qJ(4) * mrSges(5,1) + t816) * t717;
t791 = pkin(8) * t620 + t609 * t778 + t610 * t782;
t608 = Ifges(3,5) * t758 - Ifges(3,6) * t759 + Ifges(3,3) * t769 + mrSges(3,1) * t719 - mrSges(3,2) * t720 + t777 * t616 + t820 * t611 + pkin(2) * t787 + pkin(9) * t796 + (t736 * t778 - t737 * t782) * t803;
t607 = -mrSges(2,2) * g(3) - mrSges(2,3) * t766 + Ifges(2,5) * qJDD(1) - t784 * Ifges(2,6) + t782 * t609 - t778 * t610 + (-t614 * t773 - t615 * t774) * pkin(8);
t606 = mrSges(2,1) * g(3) + mrSges(2,3) * t767 + t784 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t614 - t773 * t608 + t791 * t774;
t1 = [-m(1) * g(1) + t797; -m(1) * g(2) + t809; (-m(1) - m(2)) * g(3) + t614; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t809 - t779 * t606 + t783 * t607; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t797 + t783 * t606 + t779 * t607; -mrSges(1,1) * g(2) + mrSges(2,1) * t766 + mrSges(1,2) * g(1) - mrSges(2,2) * t767 + Ifges(2,3) * qJDD(1) + pkin(1) * t615 + t774 * t608 + t791 * t773;];
tauB  = t1;
