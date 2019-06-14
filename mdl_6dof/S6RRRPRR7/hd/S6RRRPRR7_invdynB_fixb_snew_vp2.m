% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 11:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:31:50
% EndTime: 2019-05-07 11:33:27
% DurationCPUTime: 62.35s
% Computational Cost: add. (1006555->398), mult. (2231409->516), div. (0->0), fcn. (1825033->14), ass. (0->162)
t765 = cos(pkin(6));
t798 = g(3) * t765;
t763 = sin(pkin(6));
t769 = sin(qJ(2));
t797 = t763 * t769;
t774 = cos(qJ(2));
t796 = t763 * t774;
t795 = t765 * t769;
t794 = t765 * t774;
t770 = sin(qJ(1));
t775 = cos(qJ(1));
t754 = t770 * g(1) - g(2) * t775;
t776 = qJD(1) ^ 2;
t745 = pkin(8) * t763 * t776 + qJDD(1) * pkin(1) + t754;
t755 = -g(1) * t775 - g(2) * t770;
t789 = qJDD(1) * t763;
t746 = -pkin(1) * t776 + pkin(8) * t789 + t755;
t792 = t745 * t795 + t774 * t746;
t720 = -g(3) * t797 + t792;
t759 = qJD(1) * t765 + qJD(2);
t791 = qJD(1) * t763;
t788 = t769 * t791;
t743 = mrSges(3,1) * t759 - mrSges(3,3) * t788;
t747 = (-mrSges(3,1) * t774 + mrSges(3,2) * t769) * t791;
t750 = -qJD(2) * t788 + t774 * t789;
t758 = qJDD(1) * t765 + qJDD(2);
t748 = (-pkin(2) * t774 - pkin(9) * t769) * t791;
t757 = t759 ^ 2;
t790 = qJD(1) * t774;
t706 = -pkin(2) * t757 + pkin(9) * t758 + (-g(3) * t769 + t748 * t790) * t763 + t792;
t749 = (qJD(2) * t790 + qJDD(1) * t769) * t763;
t707 = -pkin(2) * t750 - pkin(9) * t749 - t798 + (-t745 + (pkin(2) * t769 - pkin(9) * t774) * t759 * qJD(1)) * t763;
t768 = sin(qJ(3));
t773 = cos(qJ(3));
t684 = -t706 * t768 + t773 * t707;
t737 = t759 * t773 - t768 * t788;
t718 = qJD(3) * t737 + t749 * t773 + t758 * t768;
t738 = t759 * t768 + t773 * t788;
t742 = qJDD(3) - t750;
t787 = t763 * t790;
t753 = qJD(3) - t787;
t670 = (t737 * t753 - t718) * qJ(4) + (t737 * t738 + t742) * pkin(3) + t684;
t685 = t773 * t706 + t768 * t707;
t717 = -qJD(3) * t738 - t749 * t768 + t758 * t773;
t728 = pkin(3) * t753 - qJ(4) * t738;
t736 = t737 ^ 2;
t672 = -pkin(3) * t736 + qJ(4) * t717 - t728 * t753 + t685;
t762 = sin(pkin(12));
t764 = cos(pkin(12));
t725 = t737 * t762 + t738 * t764;
t652 = -0.2e1 * qJD(4) * t725 + t764 * t670 - t672 * t762;
t692 = t717 * t762 + t718 * t764;
t724 = t737 * t764 - t738 * t762;
t649 = (t724 * t753 - t692) * pkin(10) + (t724 * t725 + t742) * pkin(4) + t652;
t653 = 0.2e1 * qJD(4) * t724 + t762 * t670 + t764 * t672;
t691 = t717 * t764 - t718 * t762;
t710 = pkin(4) * t753 - pkin(10) * t725;
t723 = t724 ^ 2;
t651 = -pkin(4) * t723 + pkin(10) * t691 - t710 * t753 + t653;
t767 = sin(qJ(5));
t772 = cos(qJ(5));
t646 = t767 * t649 + t772 * t651;
t700 = t724 * t767 + t725 * t772;
t668 = -qJD(5) * t700 + t691 * t772 - t692 * t767;
t699 = t724 * t772 - t725 * t767;
t682 = -mrSges(6,1) * t699 + mrSges(6,2) * t700;
t752 = qJD(5) + t753;
t689 = mrSges(6,1) * t752 - mrSges(6,3) * t700;
t741 = qJDD(5) + t742;
t683 = -pkin(5) * t699 - pkin(11) * t700;
t751 = t752 ^ 2;
t644 = -pkin(5) * t751 + pkin(11) * t741 + t683 * t699 + t646;
t719 = -g(3) * t796 + t745 * t794 - t769 * t746;
t705 = -pkin(2) * t758 - pkin(9) * t757 + t748 * t788 - t719;
t676 = -pkin(3) * t717 - qJ(4) * t736 + t738 * t728 + qJDD(4) + t705;
t655 = -pkin(4) * t691 - pkin(10) * t723 + t725 * t710 + t676;
t669 = qJD(5) * t699 + t691 * t767 + t692 * t772;
t647 = t655 + (-t699 * t752 - t669) * pkin(11) + (t700 * t752 - t668) * pkin(5);
t766 = sin(qJ(6));
t771 = cos(qJ(6));
t641 = -t644 * t766 + t647 * t771;
t686 = -t700 * t766 + t752 * t771;
t658 = qJD(6) * t686 + t669 * t771 + t741 * t766;
t667 = qJDD(6) - t668;
t687 = t700 * t771 + t752 * t766;
t673 = -mrSges(7,1) * t686 + mrSges(7,2) * t687;
t698 = qJD(6) - t699;
t674 = -mrSges(7,2) * t698 + mrSges(7,3) * t686;
t639 = m(7) * t641 + mrSges(7,1) * t667 - mrSges(7,3) * t658 - t673 * t687 + t674 * t698;
t642 = t644 * t771 + t647 * t766;
t657 = -qJD(6) * t687 - t669 * t766 + t741 * t771;
t675 = mrSges(7,1) * t698 - mrSges(7,3) * t687;
t640 = m(7) * t642 - mrSges(7,2) * t667 + mrSges(7,3) * t657 + t673 * t686 - t675 * t698;
t782 = -t639 * t766 + t771 * t640;
t627 = m(6) * t646 - mrSges(6,2) * t741 + mrSges(6,3) * t668 + t682 * t699 - t689 * t752 + t782;
t645 = t649 * t772 - t651 * t767;
t688 = -mrSges(6,2) * t752 + mrSges(6,3) * t699;
t643 = -pkin(5) * t741 - pkin(11) * t751 + t683 * t700 - t645;
t779 = -m(7) * t643 + t657 * mrSges(7,1) - mrSges(7,2) * t658 + t686 * t674 - t675 * t687;
t635 = m(6) * t645 + mrSges(6,1) * t741 - mrSges(6,3) * t669 - t682 * t700 + t688 * t752 + t779;
t624 = t767 * t627 + t772 * t635;
t701 = -mrSges(5,1) * t724 + mrSges(5,2) * t725;
t708 = -mrSges(5,2) * t753 + mrSges(5,3) * t724;
t622 = m(5) * t652 + mrSges(5,1) * t742 - mrSges(5,3) * t692 - t701 * t725 + t708 * t753 + t624;
t709 = mrSges(5,1) * t753 - mrSges(5,3) * t725;
t783 = t772 * t627 - t635 * t767;
t623 = m(5) * t653 - mrSges(5,2) * t742 + mrSges(5,3) * t691 + t701 * t724 - t709 * t753 + t783;
t616 = t764 * t622 + t762 * t623;
t726 = -mrSges(4,1) * t737 + mrSges(4,2) * t738;
t727 = -mrSges(4,2) * t753 + mrSges(4,3) * t737;
t614 = m(4) * t684 + mrSges(4,1) * t742 - mrSges(4,3) * t718 - t726 * t738 + t727 * t753 + t616;
t729 = mrSges(4,1) * t753 - mrSges(4,3) * t738;
t784 = -t622 * t762 + t764 * t623;
t615 = m(4) * t685 - mrSges(4,2) * t742 + mrSges(4,3) * t717 + t726 * t737 - t729 * t753 + t784;
t785 = -t614 * t768 + t773 * t615;
t606 = m(3) * t720 - mrSges(3,2) * t758 + mrSges(3,3) * t750 - t743 * t759 + t747 * t787 + t785;
t609 = t773 * t614 + t768 * t615;
t733 = -t745 * t763 - t798;
t744 = -mrSges(3,2) * t759 + mrSges(3,3) * t787;
t608 = m(3) * t733 - mrSges(3,1) * t750 + mrSges(3,2) * t749 + (t743 * t769 - t744 * t774) * t791 + t609;
t631 = t771 * t639 + t766 * t640;
t781 = m(6) * t655 - t668 * mrSges(6,1) + t669 * mrSges(6,2) - t699 * t688 + t700 * t689 + t631;
t778 = m(5) * t676 - t691 * mrSges(5,1) + t692 * mrSges(5,2) - t724 * t708 + t725 * t709 + t781;
t777 = -m(4) * t705 + t717 * mrSges(4,1) - t718 * mrSges(4,2) + t737 * t727 - t738 * t729 - t778;
t630 = m(3) * t719 + t758 * mrSges(3,1) - t749 * mrSges(3,3) + t759 * t744 - t747 * t788 + t777;
t597 = t606 * t795 - t608 * t763 + t630 * t794;
t595 = m(2) * t754 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t776 + t597;
t601 = t774 * t606 - t630 * t769;
t600 = m(2) * t755 - mrSges(2,1) * t776 - qJDD(1) * mrSges(2,2) + t601;
t793 = t775 * t595 + t770 * t600;
t596 = t606 * t797 + t765 * t608 + t630 * t796;
t786 = -t595 * t770 + t775 * t600;
t659 = Ifges(7,5) * t687 + Ifges(7,6) * t686 + Ifges(7,3) * t698;
t661 = Ifges(7,1) * t687 + Ifges(7,4) * t686 + Ifges(7,5) * t698;
t632 = -mrSges(7,1) * t643 + mrSges(7,3) * t642 + Ifges(7,4) * t658 + Ifges(7,2) * t657 + Ifges(7,6) * t667 - t659 * t687 + t661 * t698;
t660 = Ifges(7,4) * t687 + Ifges(7,2) * t686 + Ifges(7,6) * t698;
t633 = mrSges(7,2) * t643 - mrSges(7,3) * t641 + Ifges(7,1) * t658 + Ifges(7,4) * t657 + Ifges(7,5) * t667 + t659 * t686 - t660 * t698;
t677 = Ifges(6,5) * t700 + Ifges(6,6) * t699 + Ifges(6,3) * t752;
t678 = Ifges(6,4) * t700 + Ifges(6,2) * t699 + Ifges(6,6) * t752;
t617 = mrSges(6,2) * t655 - mrSges(6,3) * t645 + Ifges(6,1) * t669 + Ifges(6,4) * t668 + Ifges(6,5) * t741 - pkin(11) * t631 - t632 * t766 + t633 * t771 + t677 * t699 - t678 * t752;
t679 = Ifges(6,1) * t700 + Ifges(6,4) * t699 + Ifges(6,5) * t752;
t618 = -mrSges(6,1) * t655 - mrSges(7,1) * t641 + mrSges(7,2) * t642 + mrSges(6,3) * t646 + Ifges(6,4) * t669 - Ifges(7,5) * t658 + Ifges(6,2) * t668 + Ifges(6,6) * t741 - Ifges(7,6) * t657 - Ifges(7,3) * t667 - pkin(5) * t631 - t660 * t687 + t661 * t686 - t677 * t700 + t679 * t752;
t693 = Ifges(5,5) * t725 + Ifges(5,6) * t724 + Ifges(5,3) * t753;
t695 = Ifges(5,1) * t725 + Ifges(5,4) * t724 + Ifges(5,5) * t753;
t602 = -mrSges(5,1) * t676 + mrSges(5,3) * t653 + Ifges(5,4) * t692 + Ifges(5,2) * t691 + Ifges(5,6) * t742 - pkin(4) * t781 + pkin(10) * t783 + t767 * t617 + t772 * t618 - t725 * t693 + t753 * t695;
t694 = Ifges(5,4) * t725 + Ifges(5,2) * t724 + Ifges(5,6) * t753;
t610 = mrSges(5,2) * t676 - mrSges(5,3) * t652 + Ifges(5,1) * t692 + Ifges(5,4) * t691 + Ifges(5,5) * t742 - pkin(10) * t624 + t617 * t772 - t618 * t767 + t693 * t724 - t694 * t753;
t711 = Ifges(4,5) * t738 + Ifges(4,6) * t737 + Ifges(4,3) * t753;
t713 = Ifges(4,1) * t738 + Ifges(4,4) * t737 + Ifges(4,5) * t753;
t591 = -mrSges(4,1) * t705 + mrSges(4,3) * t685 + Ifges(4,4) * t718 + Ifges(4,2) * t717 + Ifges(4,6) * t742 - pkin(3) * t778 + qJ(4) * t784 + t764 * t602 + t762 * t610 - t738 * t711 + t753 * t713;
t712 = Ifges(4,4) * t738 + Ifges(4,2) * t737 + Ifges(4,6) * t753;
t592 = mrSges(4,2) * t705 - mrSges(4,3) * t684 + Ifges(4,1) * t718 + Ifges(4,4) * t717 + Ifges(4,5) * t742 - qJ(4) * t616 - t602 * t762 + t610 * t764 + t711 * t737 - t712 * t753;
t730 = Ifges(3,3) * t759 + (Ifges(3,5) * t769 + Ifges(3,6) * t774) * t791;
t731 = Ifges(3,6) * t759 + (Ifges(3,4) * t769 + Ifges(3,2) * t774) * t791;
t590 = mrSges(3,2) * t733 - mrSges(3,3) * t719 + Ifges(3,1) * t749 + Ifges(3,4) * t750 + Ifges(3,5) * t758 - pkin(9) * t609 - t591 * t768 + t592 * t773 + t730 * t787 - t731 * t759;
t732 = Ifges(3,5) * t759 + (Ifges(3,1) * t769 + Ifges(3,4) * t774) * t791;
t593 = -t730 * t788 - pkin(5) * t779 - pkin(11) * t782 + (-Ifges(4,3) - Ifges(5,3)) * t742 - pkin(2) * t609 - pkin(3) * t616 - t771 * t632 - t766 * t633 + Ifges(3,6) * t758 + t759 * t732 + Ifges(3,2) * t750 - Ifges(6,3) * t741 + Ifges(3,4) * t749 - mrSges(3,1) * t733 + t737 * t713 - t738 * t712 + t724 * t695 - t725 * t694 - Ifges(4,6) * t717 - Ifges(4,5) * t718 + mrSges(3,3) * t720 + t699 * t679 - t700 * t678 - Ifges(5,6) * t691 - Ifges(5,5) * t692 - mrSges(4,1) * t684 + mrSges(4,2) * t685 - Ifges(6,6) * t668 - Ifges(6,5) * t669 - mrSges(5,1) * t652 + mrSges(5,2) * t653 + mrSges(6,2) * t646 - mrSges(6,1) * t645 - pkin(4) * t624;
t780 = pkin(8) * t601 + t590 * t769 + t593 * t774;
t589 = Ifges(3,5) * t749 + Ifges(3,6) * t750 + Ifges(3,3) * t758 + mrSges(3,1) * t719 - mrSges(3,2) * t720 + t768 * t592 + t773 * t591 + pkin(2) * t777 + pkin(9) * t785 + (t731 * t769 - t732 * t774) * t791;
t588 = -mrSges(2,2) * g(3) - mrSges(2,3) * t754 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t776 + t590 * t774 - t593 * t769 + (-t596 * t763 - t597 * t765) * pkin(8);
t587 = mrSges(2,1) * g(3) + mrSges(2,3) * t755 + Ifges(2,5) * t776 + Ifges(2,6) * qJDD(1) - pkin(1) * t596 - t589 * t763 + t765 * t780;
t1 = [-m(1) * g(1) + t786; -m(1) * g(2) + t793; (-m(1) - m(2)) * g(3) + t596; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t793 - t770 * t587 + t775 * t588; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t786 + t775 * t587 + t770 * t588; -mrSges(1,1) * g(2) + mrSges(2,1) * t754 + mrSges(1,2) * g(1) - mrSges(2,2) * t755 + Ifges(2,3) * qJDD(1) + pkin(1) * t597 + t589 * t765 + t763 * t780;];
tauB  = t1;
