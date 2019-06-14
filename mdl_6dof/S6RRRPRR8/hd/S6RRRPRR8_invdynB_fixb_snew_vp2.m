% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR8
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
% Datum: 2019-05-07 12:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:16:55
% EndTime: 2019-05-07 12:18:10
% DurationCPUTime: 55.10s
% Computational Cost: add. (899065->398), mult. (1980273->517), div. (0->0), fcn. (1611086->14), ass. (0->163)
t764 = sin(pkin(6));
t770 = sin(qJ(2));
t775 = cos(qJ(2));
t792 = qJD(1) * qJD(2);
t753 = (-qJDD(1) * t775 + t770 * t792) * t764;
t795 = qJD(1) * t764;
t751 = (-pkin(2) * t775 - pkin(9) * t770) * t795;
t766 = cos(pkin(6));
t760 = qJD(1) * t766 + qJD(2);
t758 = t760 ^ 2;
t759 = qJDD(1) * t766 + qJDD(2);
t794 = qJD(1) * t775;
t771 = sin(qJ(1));
t776 = cos(qJ(1));
t756 = t771 * g(1) - g(2) * t776;
t777 = qJD(1) ^ 2;
t803 = pkin(8) * t764;
t748 = qJDD(1) * pkin(1) + t777 * t803 + t756;
t757 = -g(1) * t776 - g(2) * t771;
t749 = -pkin(1) * t777 + qJDD(1) * t803 + t757;
t799 = t766 * t770;
t796 = t748 * t799 + t775 * t749;
t705 = -t758 * pkin(2) + t759 * pkin(9) + (-g(3) * t770 + t751 * t794) * t764 + t796;
t752 = (qJDD(1) * t770 + t775 * t792) * t764;
t802 = t766 * g(3);
t706 = t753 * pkin(2) - t752 * pkin(9) - t802 + (-t748 + (pkin(2) * t770 - pkin(9) * t775) * t760 * qJD(1)) * t764;
t769 = sin(qJ(3));
t774 = cos(qJ(3));
t677 = -t769 * t705 + t774 * t706;
t791 = t770 * t795;
t741 = t760 * t774 - t769 * t791;
t720 = qJD(3) * t741 + t752 * t774 + t759 * t769;
t742 = t760 * t769 + t774 * t791;
t745 = qJDD(3) + t753;
t790 = t764 * t794;
t755 = qJD(3) - t790;
t662 = (t741 * t755 - t720) * qJ(4) + (t741 * t742 + t745) * pkin(3) + t677;
t678 = t774 * t705 + t769 * t706;
t719 = -qJD(3) * t742 - t752 * t769 + t759 * t774;
t731 = pkin(3) * t755 - qJ(4) * t742;
t740 = t741 ^ 2;
t669 = -pkin(3) * t740 + qJ(4) * t719 - t731 * t755 + t678;
t763 = sin(pkin(12));
t765 = cos(pkin(12));
t728 = t741 * t763 + t742 * t765;
t653 = -0.2e1 * qJD(4) * t728 + t662 * t765 - t763 * t669;
t801 = t764 * t770;
t800 = t764 * t775;
t798 = t766 * t775;
t723 = -g(3) * t801 + t796;
t746 = mrSges(3,1) * t760 - mrSges(3,3) * t791;
t750 = (-mrSges(3,1) * t775 + mrSges(3,2) * t770) * t795;
t727 = t741 * t765 - t763 * t742;
t654 = 0.2e1 * qJD(4) * t727 + t763 * t662 + t765 * t669;
t693 = t719 * t765 - t763 * t720;
t699 = -mrSges(5,1) * t727 + mrSges(5,2) * t728;
t711 = mrSges(5,1) * t755 - mrSges(5,3) * t728;
t700 = -pkin(4) * t727 - pkin(10) * t728;
t754 = t755 ^ 2;
t652 = -pkin(4) * t754 + pkin(10) * t745 + t700 * t727 + t654;
t722 = -g(3) * t800 + t748 * t798 - t770 * t749;
t704 = -t759 * pkin(2) - t758 * pkin(9) + t751 * t791 - t722;
t671 = -t719 * pkin(3) - t740 * qJ(4) + t742 * t731 + qJDD(4) + t704;
t694 = t719 * t763 + t720 * t765;
t657 = (-t727 * t755 - t694) * pkin(10) + (t728 * t755 - t693) * pkin(4) + t671;
t768 = sin(qJ(5));
t773 = cos(qJ(5));
t647 = -t768 * t652 + t773 * t657;
t708 = -t728 * t768 + t755 * t773;
t674 = qJD(5) * t708 + t694 * t773 + t745 * t768;
t692 = qJDD(5) - t693;
t709 = t728 * t773 + t755 * t768;
t726 = qJD(5) - t727;
t645 = (t708 * t726 - t674) * pkin(11) + (t708 * t709 + t692) * pkin(5) + t647;
t648 = t773 * t652 + t768 * t657;
t673 = -qJD(5) * t709 - t694 * t768 + t745 * t773;
t689 = pkin(5) * t726 - pkin(11) * t709;
t707 = t708 ^ 2;
t646 = -pkin(5) * t707 + pkin(11) * t673 - t689 * t726 + t648;
t767 = sin(qJ(6));
t772 = cos(qJ(6));
t643 = t645 * t772 - t646 * t767;
t684 = t708 * t772 - t709 * t767;
t660 = qJD(6) * t684 + t673 * t767 + t674 * t772;
t685 = t708 * t767 + t709 * t772;
t670 = -mrSges(7,1) * t684 + mrSges(7,2) * t685;
t721 = qJD(6) + t726;
t675 = -mrSges(7,2) * t721 + mrSges(7,3) * t684;
t690 = qJDD(6) + t692;
t641 = m(7) * t643 + mrSges(7,1) * t690 - mrSges(7,3) * t660 - t670 * t685 + t675 * t721;
t644 = t645 * t767 + t646 * t772;
t659 = -qJD(6) * t685 + t673 * t772 - t674 * t767;
t676 = mrSges(7,1) * t721 - mrSges(7,3) * t685;
t642 = m(7) * t644 - mrSges(7,2) * t690 + mrSges(7,3) * t659 + t670 * t684 - t676 * t721;
t633 = t772 * t641 + t767 * t642;
t686 = -mrSges(6,1) * t708 + mrSges(6,2) * t709;
t687 = -mrSges(6,2) * t726 + mrSges(6,3) * t708;
t631 = m(6) * t647 + mrSges(6,1) * t692 - mrSges(6,3) * t674 - t686 * t709 + t687 * t726 + t633;
t688 = mrSges(6,1) * t726 - mrSges(6,3) * t709;
t785 = -t641 * t767 + t772 * t642;
t632 = m(6) * t648 - mrSges(6,2) * t692 + mrSges(6,3) * t673 + t686 * t708 - t688 * t726 + t785;
t786 = -t631 * t768 + t773 * t632;
t626 = m(5) * t654 - mrSges(5,2) * t745 + mrSges(5,3) * t693 + t699 * t727 - t711 * t755 + t786;
t710 = -mrSges(5,2) * t755 + mrSges(5,3) * t727;
t651 = -pkin(4) * t745 - pkin(10) * t754 + t728 * t700 - t653;
t649 = -pkin(5) * t673 - pkin(11) * t707 + t689 * t709 + t651;
t781 = m(7) * t649 - t659 * mrSges(7,1) + mrSges(7,2) * t660 - t684 * t675 + t676 * t685;
t779 = -m(6) * t651 + t673 * mrSges(6,1) - mrSges(6,2) * t674 + t708 * t687 - t688 * t709 - t781;
t637 = m(5) * t653 + mrSges(5,1) * t745 - mrSges(5,3) * t694 - t699 * t728 + t710 * t755 + t779;
t618 = t763 * t626 + t765 * t637;
t729 = -mrSges(4,1) * t741 + mrSges(4,2) * t742;
t730 = -mrSges(4,2) * t755 + mrSges(4,3) * t741;
t616 = m(4) * t677 + mrSges(4,1) * t745 - mrSges(4,3) * t720 - t729 * t742 + t730 * t755 + t618;
t732 = mrSges(4,1) * t755 - mrSges(4,3) * t742;
t787 = t765 * t626 - t637 * t763;
t617 = m(4) * t678 - mrSges(4,2) * t745 + mrSges(4,3) * t719 + t729 * t741 - t732 * t755 + t787;
t788 = -t616 * t769 + t774 * t617;
t608 = m(3) * t723 - mrSges(3,2) * t759 - mrSges(3,3) * t753 - t746 * t760 + t750 * t790 + t788;
t611 = t774 * t616 + t769 * t617;
t736 = -t764 * t748 - t802;
t747 = -mrSges(3,2) * t760 + mrSges(3,3) * t790;
t610 = m(3) * t736 + t753 * mrSges(3,1) + t752 * mrSges(3,2) + (t746 * t770 - t747 * t775) * t795 + t611;
t627 = t773 * t631 + t768 * t632;
t780 = m(5) * t671 - t693 * mrSges(5,1) + mrSges(5,2) * t694 - t727 * t710 + t711 * t728 + t627;
t778 = -m(4) * t704 + t719 * mrSges(4,1) - mrSges(4,2) * t720 + t741 * t730 - t732 * t742 - t780;
t623 = m(3) * t722 + mrSges(3,1) * t759 - mrSges(3,3) * t752 + t747 * t760 - t750 * t791 + t778;
t599 = t608 * t799 - t610 * t764 + t623 * t798;
t597 = m(2) * t756 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t777 + t599;
t603 = t775 * t608 - t623 * t770;
t602 = m(2) * t757 - mrSges(2,1) * t777 - qJDD(1) * mrSges(2,2) + t603;
t797 = t776 * t597 + t771 * t602;
t598 = t608 * t801 + t766 * t610 + t623 * t800;
t789 = -t597 * t771 + t776 * t602;
t663 = Ifges(7,5) * t685 + Ifges(7,6) * t684 + Ifges(7,3) * t721;
t665 = Ifges(7,1) * t685 + Ifges(7,4) * t684 + Ifges(7,5) * t721;
t634 = -mrSges(7,1) * t649 + mrSges(7,3) * t644 + Ifges(7,4) * t660 + Ifges(7,2) * t659 + Ifges(7,6) * t690 - t663 * t685 + t665 * t721;
t664 = Ifges(7,4) * t685 + Ifges(7,2) * t684 + Ifges(7,6) * t721;
t635 = mrSges(7,2) * t649 - mrSges(7,3) * t643 + Ifges(7,1) * t660 + Ifges(7,4) * t659 + Ifges(7,5) * t690 + t663 * t684 - t664 * t721;
t679 = Ifges(6,5) * t709 + Ifges(6,6) * t708 + Ifges(6,3) * t726;
t681 = Ifges(6,1) * t709 + Ifges(6,4) * t708 + Ifges(6,5) * t726;
t619 = -mrSges(6,1) * t651 + mrSges(6,3) * t648 + Ifges(6,4) * t674 + Ifges(6,2) * t673 + Ifges(6,6) * t692 - pkin(5) * t781 + pkin(11) * t785 + t772 * t634 + t767 * t635 - t709 * t679 + t726 * t681;
t680 = Ifges(6,4) * t709 + Ifges(6,2) * t708 + Ifges(6,6) * t726;
t620 = mrSges(6,2) * t651 - mrSges(6,3) * t647 + Ifges(6,1) * t674 + Ifges(6,4) * t673 + Ifges(6,5) * t692 - pkin(11) * t633 - t634 * t767 + t635 * t772 + t679 * t708 - t680 * t726;
t695 = Ifges(5,5) * t728 + Ifges(5,6) * t727 + Ifges(5,3) * t755;
t696 = Ifges(5,4) * t728 + Ifges(5,2) * t727 + Ifges(5,6) * t755;
t604 = mrSges(5,2) * t671 - mrSges(5,3) * t653 + Ifges(5,1) * t694 + Ifges(5,4) * t693 + Ifges(5,5) * t745 - pkin(10) * t627 - t619 * t768 + t620 * t773 + t695 * t727 - t696 * t755;
t697 = Ifges(5,1) * t728 + Ifges(5,4) * t727 + Ifges(5,5) * t755;
t612 = Ifges(5,4) * t694 + Ifges(5,2) * t693 + Ifges(5,6) * t745 - t728 * t695 + t755 * t697 - mrSges(5,1) * t671 + mrSges(5,3) * t654 - Ifges(6,5) * t674 - Ifges(6,6) * t673 - Ifges(6,3) * t692 - t709 * t680 + t708 * t681 - mrSges(6,1) * t647 + mrSges(6,2) * t648 - Ifges(7,5) * t660 - Ifges(7,6) * t659 - Ifges(7,3) * t690 - t685 * t664 + t684 * t665 - mrSges(7,1) * t643 + mrSges(7,2) * t644 - pkin(5) * t633 - pkin(4) * t627;
t713 = Ifges(4,5) * t742 + Ifges(4,6) * t741 + Ifges(4,3) * t755;
t715 = Ifges(4,1) * t742 + Ifges(4,4) * t741 + Ifges(4,5) * t755;
t593 = -mrSges(4,1) * t704 + mrSges(4,3) * t678 + Ifges(4,4) * t720 + Ifges(4,2) * t719 + Ifges(4,6) * t745 - pkin(3) * t780 + qJ(4) * t787 + t763 * t604 + t765 * t612 - t742 * t713 + t755 * t715;
t714 = Ifges(4,4) * t742 + Ifges(4,2) * t741 + Ifges(4,6) * t755;
t595 = mrSges(4,2) * t704 - mrSges(4,3) * t677 + Ifges(4,1) * t720 + Ifges(4,4) * t719 + Ifges(4,5) * t745 - qJ(4) * t618 + t604 * t765 - t612 * t763 + t713 * t741 - t714 * t755;
t733 = Ifges(3,3) * t760 + (Ifges(3,5) * t770 + Ifges(3,6) * t775) * t795;
t734 = Ifges(3,6) * t760 + (Ifges(3,4) * t770 + Ifges(3,2) * t775) * t795;
t592 = mrSges(3,2) * t736 - mrSges(3,3) * t722 + Ifges(3,1) * t752 - Ifges(3,4) * t753 + Ifges(3,5) * t759 - pkin(9) * t611 - t593 * t769 + t595 * t774 + t733 * t790 - t734 * t760;
t735 = Ifges(3,5) * t760 + (Ifges(3,1) * t770 + Ifges(3,4) * t775) * t795;
t594 = -t733 * t791 - pkin(10) * t786 - pkin(4) * t779 - pkin(2) * t611 + (-Ifges(4,3) - Ifges(5,3)) * t745 - pkin(3) * t618 - t768 * t620 - t773 * t619 + t760 * t735 + Ifges(3,6) * t759 + Ifges(3,4) * t752 - Ifges(3,2) * t753 - mrSges(3,1) * t736 + t741 * t715 - t742 * t714 + t727 * t697 - t728 * t696 - Ifges(4,6) * t719 - Ifges(4,5) * t720 + mrSges(3,3) * t723 - Ifges(5,6) * t693 - Ifges(5,5) * t694 - mrSges(4,1) * t677 + mrSges(4,2) * t678 - mrSges(5,1) * t653 + mrSges(5,2) * t654;
t782 = pkin(8) * t603 + t592 * t770 + t594 * t775;
t591 = Ifges(3,5) * t752 - Ifges(3,6) * t753 + Ifges(3,3) * t759 + mrSges(3,1) * t722 - mrSges(3,2) * t723 + t769 * t595 + t774 * t593 + pkin(2) * t778 + pkin(9) * t788 + (t734 * t770 - t735 * t775) * t795;
t590 = -mrSges(2,2) * g(3) - mrSges(2,3) * t756 + Ifges(2,5) * qJDD(1) - t777 * Ifges(2,6) + t775 * t592 - t770 * t594 + (-t598 * t764 - t599 * t766) * pkin(8);
t589 = mrSges(2,1) * g(3) + mrSges(2,3) * t757 + t777 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t598 - t764 * t591 + t766 * t782;
t1 = [-m(1) * g(1) + t789; -m(1) * g(2) + t797; (-m(1) - m(2)) * g(3) + t598; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t797 - t771 * t589 + t776 * t590; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t789 + t776 * t589 + t771 * t590; -mrSges(1,1) * g(2) + mrSges(2,1) * t756 + mrSges(1,2) * g(1) - mrSges(2,2) * t757 + Ifges(2,3) * qJDD(1) + pkin(1) * t599 + t766 * t591 + t764 * t782;];
tauB  = t1;
