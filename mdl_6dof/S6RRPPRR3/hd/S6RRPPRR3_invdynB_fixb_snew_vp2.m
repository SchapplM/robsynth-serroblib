% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 10:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:01:31
% EndTime: 2019-05-06 10:02:17
% DurationCPUTime: 43.79s
% Computational Cost: add. (646296->396), mult. (1745992->516), div. (0->0), fcn. (1391951->14), ass. (0->161)
t801 = -2 * qJD(3);
t763 = sin(pkin(6));
t768 = sin(qJ(2));
t772 = cos(qJ(2));
t790 = qJD(1) * qJD(2);
t750 = (qJDD(1) * t768 + t772 * t790) * t763;
t765 = cos(pkin(6));
t756 = qJDD(1) * t765 + qJDD(2);
t757 = qJD(1) * t765 + qJD(2);
t769 = sin(qJ(1));
t773 = cos(qJ(1));
t753 = t769 * g(1) - g(2) * t773;
t774 = qJD(1) ^ 2;
t800 = pkin(8) * t763;
t747 = qJDD(1) * pkin(1) + t774 * t800 + t753;
t754 = -g(1) * t773 - g(2) * t769;
t748 = -pkin(1) * t774 + qJDD(1) * t800 + t754;
t794 = t765 * t772;
t781 = t747 * t794 - t748 * t768;
t798 = t763 ^ 2 * t774;
t688 = pkin(2) * t756 - qJ(3) * t750 + (pkin(2) * t768 * t798 + (qJ(3) * qJD(1) * t757 - g(3)) * t763) * t772 + t781;
t795 = t765 * t768;
t797 = t763 * t768;
t715 = -g(3) * t797 + t747 * t795 + t772 * t748;
t792 = qJD(1) * t763;
t788 = t768 * t792;
t744 = pkin(2) * t757 - qJ(3) * t788;
t751 = (qJDD(1) * t772 - t768 * t790) * t763;
t789 = t772 ^ 2 * t798;
t691 = -pkin(2) * t789 + qJ(3) * t751 - t744 * t757 + t715;
t762 = sin(pkin(11));
t799 = cos(pkin(11));
t741 = (t762 * t772 + t768 * t799) * t792;
t665 = t688 * t799 - t762 * t691 + t741 * t801;
t796 = t763 * t772;
t787 = t772 * t792;
t740 = t762 * t788 - t799 * t787;
t666 = t762 * t688 + t799 * t691 + t740 * t801;
t717 = mrSges(4,1) * t740 + mrSges(4,2) * t741;
t720 = t750 * t762 - t799 * t751;
t728 = mrSges(4,1) * t757 - mrSges(4,3) * t741;
t716 = pkin(3) * t740 - qJ(4) * t741;
t755 = t757 ^ 2;
t660 = -pkin(3) * t755 + qJ(4) * t756 - t716 * t740 + t666;
t732 = -g(3) * t765 - t747 * t763;
t699 = -pkin(2) * t751 - qJ(3) * t789 + t744 * t788 + qJDD(3) + t732;
t721 = t750 * t799 + t762 * t751;
t669 = (t740 * t757 - t721) * qJ(4) + (t741 * t757 + t720) * pkin(3) + t699;
t761 = sin(pkin(12));
t764 = cos(pkin(12));
t726 = t741 * t764 + t757 * t761;
t652 = -0.2e1 * qJD(4) * t726 - t660 * t761 + t764 * t669;
t709 = t721 * t764 + t756 * t761;
t725 = -t741 * t761 + t757 * t764;
t649 = (t725 * t740 - t709) * pkin(9) + (t725 * t726 + t720) * pkin(4) + t652;
t653 = 0.2e1 * qJD(4) * t725 + t764 * t660 + t761 * t669;
t706 = pkin(4) * t740 - pkin(9) * t726;
t708 = -t721 * t761 + t756 * t764;
t724 = t725 ^ 2;
t651 = -pkin(4) * t724 + pkin(9) * t708 - t706 * t740 + t653;
t767 = sin(qJ(5));
t771 = cos(qJ(5));
t646 = t767 * t649 + t771 * t651;
t701 = t725 * t767 + t726 * t771;
t673 = -qJD(5) * t701 + t708 * t771 - t709 * t767;
t700 = t725 * t771 - t726 * t767;
t681 = -mrSges(6,1) * t700 + mrSges(6,2) * t701;
t739 = qJD(5) + t740;
t690 = mrSges(6,1) * t739 - mrSges(6,3) * t701;
t719 = qJDD(5) + t720;
t682 = -pkin(5) * t700 - pkin(10) * t701;
t738 = t739 ^ 2;
t644 = -pkin(5) * t738 + pkin(10) * t719 + t682 * t700 + t646;
t659 = -t756 * pkin(3) - t755 * qJ(4) + t741 * t716 + qJDD(4) - t665;
t654 = -t708 * pkin(4) - t724 * pkin(9) + t726 * t706 + t659;
t674 = qJD(5) * t700 + t708 * t767 + t709 * t771;
t647 = (-t700 * t739 - t674) * pkin(10) + (t701 * t739 - t673) * pkin(5) + t654;
t766 = sin(qJ(6));
t770 = cos(qJ(6));
t641 = -t644 * t766 + t647 * t770;
t684 = -t701 * t766 + t739 * t770;
t657 = qJD(6) * t684 + t674 * t770 + t719 * t766;
t685 = t701 * t770 + t739 * t766;
t670 = -mrSges(7,1) * t684 + mrSges(7,2) * t685;
t672 = qJDD(6) - t673;
t698 = qJD(6) - t700;
t675 = -mrSges(7,2) * t698 + mrSges(7,3) * t684;
t639 = m(7) * t641 + mrSges(7,1) * t672 - mrSges(7,3) * t657 - t670 * t685 + t675 * t698;
t642 = t644 * t770 + t647 * t766;
t656 = -qJD(6) * t685 - t674 * t766 + t719 * t770;
t676 = mrSges(7,1) * t698 - mrSges(7,3) * t685;
t640 = m(7) * t642 - mrSges(7,2) * t672 + mrSges(7,3) * t656 + t670 * t684 - t676 * t698;
t782 = -t639 * t766 + t770 * t640;
t630 = m(6) * t646 - mrSges(6,2) * t719 + mrSges(6,3) * t673 + t681 * t700 - t690 * t739 + t782;
t645 = t649 * t771 - t651 * t767;
t689 = -mrSges(6,2) * t739 + mrSges(6,3) * t700;
t643 = -pkin(5) * t719 - pkin(10) * t738 + t682 * t701 - t645;
t778 = -m(7) * t643 + t656 * mrSges(7,1) - mrSges(7,2) * t657 + t684 * t675 - t676 * t685;
t635 = m(6) * t645 + mrSges(6,1) * t719 - mrSges(6,3) * t674 - t681 * t701 + t689 * t739 + t778;
t625 = t767 * t630 + t771 * t635;
t702 = -mrSges(5,1) * t725 + mrSges(5,2) * t726;
t704 = -mrSges(5,2) * t740 + mrSges(5,3) * t725;
t623 = m(5) * t652 + mrSges(5,1) * t720 - mrSges(5,3) * t709 - t702 * t726 + t704 * t740 + t625;
t705 = mrSges(5,1) * t740 - mrSges(5,3) * t726;
t783 = t771 * t630 - t635 * t767;
t624 = m(5) * t653 - mrSges(5,2) * t720 + mrSges(5,3) * t708 + t702 * t725 - t705 * t740 + t783;
t784 = -t623 * t761 + t764 * t624;
t614 = m(4) * t666 - mrSges(4,2) * t756 - mrSges(4,3) * t720 - t717 * t740 - t728 * t757 + t784;
t727 = -mrSges(4,2) * t757 - mrSges(4,3) * t740;
t631 = t770 * t639 + t766 * t640;
t776 = m(6) * t654 - t673 * mrSges(6,1) + mrSges(6,2) * t674 - t700 * t689 + t690 * t701 + t631;
t775 = -m(5) * t659 + t708 * mrSges(5,1) - mrSges(5,2) * t709 + t725 * t704 - t705 * t726 - t776;
t627 = m(4) * t665 + mrSges(4,1) * t756 - mrSges(4,3) * t721 - t717 * t741 + t727 * t757 + t775;
t611 = t762 * t614 + t799 * t627;
t714 = -g(3) * t796 + t781;
t746 = -mrSges(3,2) * t757 + mrSges(3,3) * t787;
t749 = (-mrSges(3,1) * t772 + mrSges(3,2) * t768) * t792;
t609 = m(3) * t714 + mrSges(3,1) * t756 - mrSges(3,3) * t750 + t746 * t757 - t749 * t788 + t611;
t745 = mrSges(3,1) * t757 - mrSges(3,3) * t788;
t785 = t799 * t614 - t627 * t762;
t610 = m(3) * t715 - mrSges(3,2) * t756 + mrSges(3,3) * t751 - t745 * t757 + t749 * t787 + t785;
t617 = t764 * t623 + t761 * t624;
t777 = m(4) * t699 + mrSges(4,1) * t720 + t721 * mrSges(4,2) + t727 * t740 + t741 * t728 + t617;
t616 = m(3) * t732 - mrSges(3,1) * t751 + mrSges(3,2) * t750 + (t745 * t768 - t746 * t772) * t792 + t777;
t596 = t609 * t794 + t610 * t795 - t616 * t763;
t594 = m(2) * t753 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t774 + t596;
t600 = -t609 * t768 + t772 * t610;
t599 = m(2) * t754 - mrSges(2,1) * t774 - qJDD(1) * mrSges(2,2) + t600;
t793 = t773 * t594 + t769 * t599;
t595 = t609 * t796 + t610 * t797 + t765 * t616;
t786 = -t594 * t769 + t773 * t599;
t661 = Ifges(7,5) * t685 + Ifges(7,6) * t684 + Ifges(7,3) * t698;
t663 = Ifges(7,1) * t685 + Ifges(7,4) * t684 + Ifges(7,5) * t698;
t632 = -mrSges(7,1) * t643 + mrSges(7,3) * t642 + Ifges(7,4) * t657 + Ifges(7,2) * t656 + Ifges(7,6) * t672 - t661 * t685 + t663 * t698;
t662 = Ifges(7,4) * t685 + Ifges(7,2) * t684 + Ifges(7,6) * t698;
t633 = mrSges(7,2) * t643 - mrSges(7,3) * t641 + Ifges(7,1) * t657 + Ifges(7,4) * t656 + Ifges(7,5) * t672 + t661 * t684 - t662 * t698;
t677 = Ifges(6,5) * t701 + Ifges(6,6) * t700 + Ifges(6,3) * t739;
t678 = Ifges(6,4) * t701 + Ifges(6,2) * t700 + Ifges(6,6) * t739;
t618 = mrSges(6,2) * t654 - mrSges(6,3) * t645 + Ifges(6,1) * t674 + Ifges(6,4) * t673 + Ifges(6,5) * t719 - pkin(10) * t631 - t632 * t766 + t633 * t770 + t677 * t700 - t678 * t739;
t679 = Ifges(6,1) * t701 + Ifges(6,4) * t700 + Ifges(6,5) * t739;
t619 = -mrSges(6,1) * t654 - mrSges(7,1) * t641 + mrSges(7,2) * t642 + mrSges(6,3) * t646 + Ifges(6,4) * t674 - Ifges(7,5) * t657 + Ifges(6,2) * t673 + Ifges(6,6) * t719 - Ifges(7,6) * t656 - Ifges(7,3) * t672 - pkin(5) * t631 - t662 * t685 + t663 * t684 - t677 * t701 + t679 * t739;
t692 = Ifges(5,5) * t726 + Ifges(5,6) * t725 + Ifges(5,3) * t740;
t694 = Ifges(5,1) * t726 + Ifges(5,4) * t725 + Ifges(5,5) * t740;
t602 = -mrSges(5,1) * t659 + mrSges(5,3) * t653 + Ifges(5,4) * t709 + Ifges(5,2) * t708 + Ifges(5,6) * t720 - pkin(4) * t776 + pkin(9) * t783 + t767 * t618 + t771 * t619 - t726 * t692 + t740 * t694;
t693 = Ifges(5,4) * t726 + Ifges(5,2) * t725 + Ifges(5,6) * t740;
t603 = mrSges(5,2) * t659 - mrSges(5,3) * t652 + Ifges(5,1) * t709 + Ifges(5,4) * t708 + Ifges(5,5) * t720 - pkin(9) * t625 + t618 * t771 - t619 * t767 + t692 * t725 - t693 * t740;
t710 = Ifges(4,5) * t741 - Ifges(4,6) * t740 + Ifges(4,3) * t757;
t711 = Ifges(4,4) * t741 - Ifges(4,2) * t740 + Ifges(4,6) * t757;
t592 = mrSges(4,2) * t699 - mrSges(4,3) * t665 + Ifges(4,1) * t721 - Ifges(4,4) * t720 + Ifges(4,5) * t756 - qJ(4) * t617 - t602 * t761 + t603 * t764 - t710 * t740 - t711 * t757;
t712 = Ifges(4,1) * t741 - Ifges(4,4) * t740 + Ifges(4,5) * t757;
t601 = -pkin(4) * t625 + (-Ifges(5,3) - Ifges(4,2)) * t720 - pkin(10) * t782 - pkin(3) * t617 - t770 * t632 - t766 * t633 + Ifges(4,6) * t756 + t757 * t712 - t741 * t710 + t725 * t694 - t726 * t693 - Ifges(6,3) * t719 + Ifges(4,4) * t721 - Ifges(5,6) * t708 - Ifges(5,5) * t709 - t701 * t678 - mrSges(4,1) * t699 + t700 * t679 - Ifges(6,6) * t673 - Ifges(6,5) * t674 + mrSges(4,3) * t666 - mrSges(5,1) * t652 + mrSges(5,2) * t653 - mrSges(6,1) * t645 + mrSges(6,2) * t646 - pkin(5) * t778;
t729 = Ifges(3,3) * t757 + (Ifges(3,5) * t768 + Ifges(3,6) * t772) * t792;
t731 = Ifges(3,5) * t757 + (Ifges(3,1) * t768 + Ifges(3,4) * t772) * t792;
t589 = -mrSges(3,1) * t732 + mrSges(3,3) * t715 + Ifges(3,4) * t750 + Ifges(3,2) * t751 + Ifges(3,6) * t756 - pkin(2) * t777 + qJ(3) * t785 + t762 * t592 + t601 * t799 - t729 * t788 + t757 * t731;
t730 = Ifges(3,6) * t757 + (Ifges(3,4) * t768 + Ifges(3,2) * t772) * t792;
t590 = mrSges(3,2) * t732 - mrSges(3,3) * t714 + Ifges(3,1) * t750 + Ifges(3,4) * t751 + Ifges(3,5) * t756 - qJ(3) * t611 + t592 * t799 - t762 * t601 + t729 * t787 - t757 * t730;
t779 = pkin(8) * t600 + t589 * t772 + t590 * t768;
t591 = Ifges(3,5) * t750 + Ifges(3,6) * t751 + mrSges(3,1) * t714 - mrSges(3,2) * t715 + Ifges(4,5) * t721 - Ifges(4,6) * t720 + t741 * t711 + t740 * t712 + mrSges(4,1) * t665 - mrSges(4,2) * t666 + t761 * t603 + t764 * t602 + pkin(3) * t775 + qJ(4) * t784 + pkin(2) * t611 + (Ifges(3,3) + Ifges(4,3)) * t756 + (t730 * t768 - t731 * t772) * t792;
t588 = -mrSges(2,2) * g(3) - mrSges(2,3) * t753 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t774 - t589 * t768 + t590 * t772 + (-t595 * t763 - t596 * t765) * pkin(8);
t587 = mrSges(2,1) * g(3) + mrSges(2,3) * t754 + Ifges(2,5) * t774 + Ifges(2,6) * qJDD(1) - pkin(1) * t595 - t591 * t763 + t765 * t779;
t1 = [-m(1) * g(1) + t786; -m(1) * g(2) + t793; (-m(1) - m(2)) * g(3) + t595; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t793 - t769 * t587 + t773 * t588; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t786 + t773 * t587 + t769 * t588; -mrSges(1,1) * g(2) + mrSges(2,1) * t753 + mrSges(1,2) * g(1) - mrSges(2,2) * t754 + Ifges(2,3) * qJDD(1) + pkin(1) * t596 + t591 * t765 + t763 * t779;];
tauB  = t1;
