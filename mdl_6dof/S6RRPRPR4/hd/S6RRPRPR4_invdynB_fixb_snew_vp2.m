% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:33:01
% EndTime: 2019-05-06 13:34:10
% DurationCPUTime: 46.81s
% Computational Cost: add. (671737->397), mult. (1784526->517), div. (0->0), fcn. (1425174->14), ass. (0->163)
t803 = -2 * qJD(3);
t761 = sin(pkin(11));
t764 = cos(pkin(11));
t768 = sin(qJ(2));
t772 = cos(qJ(2));
t762 = sin(pkin(6));
t794 = qJD(1) * t762;
t741 = (t761 * t768 - t764 * t772) * t794;
t792 = qJD(1) * qJD(2);
t750 = (qJDD(1) * t768 + t772 * t792) * t762;
t765 = cos(pkin(6));
t755 = qJDD(1) * t765 + qJDD(2);
t756 = qJD(1) * t765 + qJD(2);
t769 = sin(qJ(1));
t773 = cos(qJ(1));
t752 = t769 * g(1) - g(2) * t773;
t774 = qJD(1) ^ 2;
t801 = pkin(8) * t762;
t747 = qJDD(1) * pkin(1) + t774 * t801 + t752;
t753 = -g(1) * t773 - g(2) * t769;
t748 = -pkin(1) * t774 + qJDD(1) * t801 + t753;
t796 = t765 * t772;
t782 = t747 * t796 - t768 * t748;
t800 = t762 ^ 2 * t774;
t690 = t755 * pkin(2) - t750 * qJ(3) + (pkin(2) * t768 * t800 + (qJ(3) * qJD(1) * t756 - g(3)) * t762) * t772 + t782;
t797 = t765 * t768;
t799 = t762 * t768;
t718 = -g(3) * t799 + t747 * t797 + t772 * t748;
t790 = t768 * t794;
t744 = pkin(2) * t756 - qJ(3) * t790;
t751 = (qJDD(1) * t772 - t768 * t792) * t762;
t791 = t772 ^ 2 * t800;
t693 = -pkin(2) * t791 + qJ(3) * t751 - t744 * t756 + t718;
t742 = (t761 * t772 + t764 * t768) * t794;
t667 = t764 * t690 - t761 * t693 + t742 * t803;
t802 = 2 * qJD(5);
t798 = t762 * t772;
t668 = t761 * t690 + t764 * t693 + t741 * t803;
t719 = mrSges(4,1) * t741 + mrSges(4,2) * t742;
t723 = -t750 * t761 + t751 * t764;
t729 = mrSges(4,1) * t756 - mrSges(4,3) * t742;
t720 = pkin(3) * t741 - pkin(9) * t742;
t754 = t756 ^ 2;
t662 = -pkin(3) * t754 + pkin(9) * t755 - t720 * t741 + t668;
t733 = -t765 * g(3) - t762 * t747;
t704 = -t751 * pkin(2) - qJ(3) * t791 + t744 * t790 + qJDD(3) + t733;
t724 = t750 * t764 + t751 * t761;
t671 = (t741 * t756 - t724) * pkin(9) + (t742 * t756 - t723) * pkin(3) + t704;
t767 = sin(qJ(4));
t771 = cos(qJ(4));
t654 = -t767 * t662 + t771 * t671;
t726 = -t742 * t767 + t756 * t771;
t701 = qJD(4) * t726 + t724 * t771 + t755 * t767;
t722 = qJDD(4) - t723;
t727 = t742 * t771 + t756 * t767;
t740 = qJD(4) + t741;
t651 = (t726 * t740 - t701) * qJ(5) + (t726 * t727 + t722) * pkin(4) + t654;
t655 = t771 * t662 + t767 * t671;
t700 = -qJD(4) * t727 - t724 * t767 + t755 * t771;
t711 = pkin(4) * t740 - qJ(5) * t727;
t725 = t726 ^ 2;
t653 = -pkin(4) * t725 + qJ(5) * t700 - t711 * t740 + t655;
t760 = sin(pkin(12));
t763 = cos(pkin(12));
t706 = t726 * t763 - t727 * t760;
t648 = t760 * t651 + t763 * t653 + t706 * t802;
t675 = t700 * t763 - t701 * t760;
t707 = t726 * t760 + t727 * t763;
t683 = -mrSges(6,1) * t706 + mrSges(6,2) * t707;
t692 = mrSges(6,1) * t740 - mrSges(6,3) * t707;
t684 = -pkin(5) * t706 - pkin(10) * t707;
t739 = t740 ^ 2;
t646 = -pkin(5) * t739 + pkin(10) * t722 + t684 * t706 + t648;
t661 = -t755 * pkin(3) - t754 * pkin(9) + t742 * t720 - t667;
t656 = -t700 * pkin(4) - t725 * qJ(5) + t727 * t711 + qJDD(5) + t661;
t676 = t700 * t760 + t701 * t763;
t649 = (-t706 * t740 - t676) * pkin(10) + (t707 * t740 - t675) * pkin(5) + t656;
t766 = sin(qJ(6));
t770 = cos(qJ(6));
t643 = -t646 * t766 + t649 * t770;
t686 = -t707 * t766 + t740 * t770;
t659 = qJD(6) * t686 + t676 * t770 + t722 * t766;
t687 = t707 * t770 + t740 * t766;
t672 = -mrSges(7,1) * t686 + mrSges(7,2) * t687;
t674 = qJDD(6) - t675;
t705 = qJD(6) - t706;
t677 = -mrSges(7,2) * t705 + mrSges(7,3) * t686;
t641 = m(7) * t643 + mrSges(7,1) * t674 - mrSges(7,3) * t659 - t672 * t687 + t677 * t705;
t644 = t646 * t770 + t649 * t766;
t658 = -qJD(6) * t687 - t676 * t766 + t722 * t770;
t678 = mrSges(7,1) * t705 - mrSges(7,3) * t687;
t642 = m(7) * t644 - mrSges(7,2) * t674 + mrSges(7,3) * t658 + t672 * t686 - t678 * t705;
t784 = -t641 * t766 + t770 * t642;
t632 = m(6) * t648 - mrSges(6,2) * t722 + mrSges(6,3) * t675 + t683 * t706 - t692 * t740 + t784;
t781 = -t763 * t651 + t760 * t653;
t647 = -0.2e1 * qJD(5) * t707 - t781;
t691 = -mrSges(6,2) * t740 + mrSges(6,3) * t706;
t645 = -t722 * pkin(5) - t739 * pkin(10) + (t802 + t684) * t707 + t781;
t778 = -m(7) * t645 + t658 * mrSges(7,1) - mrSges(7,2) * t659 + t686 * t677 - t678 * t687;
t637 = m(6) * t647 + mrSges(6,1) * t722 - mrSges(6,3) * t676 - t683 * t707 + t691 * t740 + t778;
t627 = t760 * t632 + t763 * t637;
t708 = -mrSges(5,1) * t726 + mrSges(5,2) * t727;
t710 = -mrSges(5,2) * t740 + mrSges(5,3) * t726;
t625 = m(5) * t654 + mrSges(5,1) * t722 - mrSges(5,3) * t701 - t708 * t727 + t710 * t740 + t627;
t712 = mrSges(5,1) * t740 - mrSges(5,3) * t727;
t785 = t763 * t632 - t637 * t760;
t626 = m(5) * t655 - mrSges(5,2) * t722 + mrSges(5,3) * t700 + t708 * t726 - t712 * t740 + t785;
t786 = -t625 * t767 + t771 * t626;
t616 = m(4) * t668 - mrSges(4,2) * t755 + mrSges(4,3) * t723 - t719 * t741 - t729 * t756 + t786;
t728 = -mrSges(4,2) * t756 - mrSges(4,3) * t741;
t633 = t770 * t641 + t766 * t642;
t776 = m(6) * t656 - t675 * mrSges(6,1) + mrSges(6,2) * t676 - t706 * t691 + t692 * t707 + t633;
t775 = -m(5) * t661 + t700 * mrSges(5,1) - mrSges(5,2) * t701 + t726 * t710 - t712 * t727 - t776;
t629 = m(4) * t667 + mrSges(4,1) * t755 - mrSges(4,3) * t724 - t719 * t742 + t728 * t756 + t775;
t613 = t761 * t616 + t764 * t629;
t717 = -g(3) * t798 + t782;
t789 = t772 * t794;
t746 = -mrSges(3,2) * t756 + mrSges(3,3) * t789;
t749 = (-mrSges(3,1) * t772 + mrSges(3,2) * t768) * t794;
t611 = m(3) * t717 + mrSges(3,1) * t755 - mrSges(3,3) * t750 + t746 * t756 - t749 * t790 + t613;
t745 = mrSges(3,1) * t756 - mrSges(3,3) * t790;
t787 = t764 * t616 - t629 * t761;
t612 = m(3) * t718 - mrSges(3,2) * t755 + mrSges(3,3) * t751 - t745 * t756 + t749 * t789 + t787;
t619 = t771 * t625 + t767 * t626;
t777 = m(4) * t704 - t723 * mrSges(4,1) + t724 * mrSges(4,2) + t741 * t728 + t742 * t729 + t619;
t618 = m(3) * t733 - t751 * mrSges(3,1) + t750 * mrSges(3,2) + (t745 * t768 - t746 * t772) * t794 + t777;
t598 = t611 * t796 + t612 * t797 - t618 * t762;
t596 = m(2) * t752 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t774 + t598;
t602 = -t611 * t768 + t772 * t612;
t601 = m(2) * t753 - mrSges(2,1) * t774 - qJDD(1) * mrSges(2,2) + t602;
t795 = t773 * t596 + t769 * t601;
t597 = t611 * t798 + t612 * t799 + t765 * t618;
t788 = -t596 * t769 + t773 * t601;
t663 = Ifges(7,5) * t687 + Ifges(7,6) * t686 + Ifges(7,3) * t705;
t665 = Ifges(7,1) * t687 + Ifges(7,4) * t686 + Ifges(7,5) * t705;
t634 = -mrSges(7,1) * t645 + mrSges(7,3) * t644 + Ifges(7,4) * t659 + Ifges(7,2) * t658 + Ifges(7,6) * t674 - t663 * t687 + t665 * t705;
t664 = Ifges(7,4) * t687 + Ifges(7,2) * t686 + Ifges(7,6) * t705;
t635 = mrSges(7,2) * t645 - mrSges(7,3) * t643 + Ifges(7,1) * t659 + Ifges(7,4) * t658 + Ifges(7,5) * t674 + t663 * t686 - t664 * t705;
t679 = Ifges(6,5) * t707 + Ifges(6,6) * t706 + Ifges(6,3) * t740;
t680 = Ifges(6,4) * t707 + Ifges(6,2) * t706 + Ifges(6,6) * t740;
t620 = mrSges(6,2) * t656 - mrSges(6,3) * t647 + Ifges(6,1) * t676 + Ifges(6,4) * t675 + Ifges(6,5) * t722 - pkin(10) * t633 - t634 * t766 + t635 * t770 + t679 * t706 - t680 * t740;
t681 = Ifges(6,1) * t707 + Ifges(6,4) * t706 + Ifges(6,5) * t740;
t621 = -mrSges(6,1) * t656 - mrSges(7,1) * t643 + mrSges(7,2) * t644 + mrSges(6,3) * t648 + Ifges(6,4) * t676 - Ifges(7,5) * t659 + Ifges(6,2) * t675 + Ifges(6,6) * t722 - Ifges(7,6) * t658 - Ifges(7,3) * t674 - pkin(5) * t633 - t664 * t687 + t665 * t686 - t679 * t707 + t681 * t740;
t694 = Ifges(5,5) * t727 + Ifges(5,6) * t726 + Ifges(5,3) * t740;
t696 = Ifges(5,1) * t727 + Ifges(5,4) * t726 + Ifges(5,5) * t740;
t604 = -mrSges(5,1) * t661 + mrSges(5,3) * t655 + Ifges(5,4) * t701 + Ifges(5,2) * t700 + Ifges(5,6) * t722 - pkin(4) * t776 + qJ(5) * t785 + t760 * t620 + t763 * t621 - t727 * t694 + t740 * t696;
t695 = Ifges(5,4) * t727 + Ifges(5,2) * t726 + Ifges(5,6) * t740;
t605 = mrSges(5,2) * t661 - mrSges(5,3) * t654 + Ifges(5,1) * t701 + Ifges(5,4) * t700 + Ifges(5,5) * t722 - qJ(5) * t627 + t620 * t763 - t621 * t760 + t694 * t726 - t695 * t740;
t713 = Ifges(4,5) * t742 - Ifges(4,6) * t741 + Ifges(4,3) * t756;
t714 = Ifges(4,4) * t742 - Ifges(4,2) * t741 + Ifges(4,6) * t756;
t594 = mrSges(4,2) * t704 - mrSges(4,3) * t667 + Ifges(4,1) * t724 + Ifges(4,4) * t723 + Ifges(4,5) * t755 - pkin(9) * t619 - t604 * t767 + t605 * t771 - t713 * t741 - t714 * t756;
t715 = Ifges(4,1) * t742 - Ifges(4,4) * t741 + Ifges(4,5) * t756;
t603 = -pkin(10) * t784 - pkin(3) * t619 - pkin(5) * t778 + (-Ifges(5,3) - Ifges(6,3)) * t722 - pkin(4) * t627 - t770 * t634 - t766 * t635 + t756 * t715 + Ifges(4,6) * t755 - t742 * t713 + t726 * t696 - t727 * t695 + Ifges(4,2) * t723 + Ifges(4,4) * t724 - t707 * t680 - mrSges(4,1) * t704 + t706 * t681 - Ifges(5,6) * t700 - Ifges(5,5) * t701 - Ifges(6,6) * t675 - Ifges(6,5) * t676 + mrSges(4,3) * t668 + mrSges(5,2) * t655 - mrSges(5,1) * t654 + mrSges(6,2) * t648 - mrSges(6,1) * t647;
t730 = Ifges(3,3) * t756 + (Ifges(3,5) * t768 + Ifges(3,6) * t772) * t794;
t732 = Ifges(3,5) * t756 + (Ifges(3,1) * t768 + Ifges(3,4) * t772) * t794;
t591 = -mrSges(3,1) * t733 + mrSges(3,3) * t718 + Ifges(3,4) * t750 + Ifges(3,2) * t751 + Ifges(3,6) * t755 - pkin(2) * t777 + qJ(3) * t787 + t761 * t594 + t764 * t603 - t730 * t790 + t756 * t732;
t731 = Ifges(3,6) * t756 + (Ifges(3,4) * t768 + Ifges(3,2) * t772) * t794;
t592 = mrSges(3,2) * t733 - mrSges(3,3) * t717 + Ifges(3,1) * t750 + Ifges(3,4) * t751 + Ifges(3,5) * t755 - qJ(3) * t613 + t594 * t764 - t603 * t761 + t730 * t789 - t731 * t756;
t779 = pkin(8) * t602 + t591 * t772 + t592 * t768;
t593 = Ifges(3,5) * t750 + Ifges(3,6) * t751 + mrSges(3,1) * t717 - mrSges(3,2) * t718 + Ifges(4,5) * t724 + Ifges(4,6) * t723 + t742 * t714 + t741 * t715 + mrSges(4,1) * t667 - mrSges(4,2) * t668 + t767 * t605 + t771 * t604 + pkin(3) * t775 + pkin(9) * t786 + pkin(2) * t613 + (Ifges(3,3) + Ifges(4,3)) * t755 + (t731 * t768 - t732 * t772) * t794;
t590 = -mrSges(2,2) * g(3) - mrSges(2,3) * t752 + Ifges(2,5) * qJDD(1) - t774 * Ifges(2,6) - t768 * t591 + t772 * t592 + (-t597 * t762 - t598 * t765) * pkin(8);
t589 = mrSges(2,1) * g(3) + mrSges(2,3) * t753 + t774 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t597 - t762 * t593 + t765 * t779;
t1 = [-m(1) * g(1) + t788; -m(1) * g(2) + t795; (-m(1) - m(2)) * g(3) + t597; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t795 - t769 * t589 + t773 * t590; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t788 + t773 * t589 + t769 * t590; -mrSges(1,1) * g(2) + mrSges(2,1) * t752 + mrSges(1,2) * g(1) - mrSges(2,2) * t753 + Ifges(2,3) * qJDD(1) + pkin(1) * t598 + t765 * t593 + t762 * t779;];
tauB  = t1;
