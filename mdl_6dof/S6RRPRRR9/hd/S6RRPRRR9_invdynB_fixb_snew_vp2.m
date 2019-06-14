% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR9
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
% Datum: 2019-05-06 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:03:22
% EndTime: 2019-05-06 23:04:23
% DurationCPUTime: 59.55s
% Computational Cost: add. (963306->398), mult. (2193199->516), div. (0->0), fcn. (1802379->14), ass. (0->162)
t763 = cos(pkin(6));
t796 = t763 * g(3);
t761 = sin(pkin(6));
t767 = sin(qJ(2));
t795 = t761 * t767;
t772 = cos(qJ(2));
t794 = t761 * t772;
t793 = t763 * t767;
t792 = t763 * t772;
t768 = sin(qJ(1));
t773 = cos(qJ(1));
t752 = t768 * g(1) - g(2) * t773;
t774 = qJD(1) ^ 2;
t743 = pkin(8) * t761 * t774 + qJDD(1) * pkin(1) + t752;
t753 = -g(1) * t773 - g(2) * t768;
t787 = qJDD(1) * t761;
t744 = -pkin(1) * t774 + pkin(8) * t787 + t753;
t790 = t743 * t793 + t772 * t744;
t715 = -g(3) * t795 + t790;
t757 = qJD(1) * t763 + qJD(2);
t789 = qJD(1) * t761;
t786 = t767 * t789;
t741 = mrSges(3,1) * t757 - mrSges(3,3) * t786;
t746 = (-mrSges(3,1) * t772 + mrSges(3,2) * t767) * t789;
t748 = -qJD(2) * t786 + t772 * t787;
t756 = qJDD(1) * t763 + qJDD(2);
t745 = (-pkin(2) * t772 - qJ(3) * t767) * t789;
t755 = t757 ^ 2;
t788 = qJD(1) * t772;
t704 = -t755 * pkin(2) + t756 * qJ(3) + (-g(3) * t767 + t745 * t788) * t761 + t790;
t747 = (qJD(2) * t788 + qJDD(1) * t767) * t761;
t705 = -t748 * pkin(2) - t796 - t747 * qJ(3) + (-t743 + (pkin(2) * t767 - qJ(3) * t772) * t757 * qJD(1)) * t761;
t760 = sin(pkin(12));
t762 = cos(pkin(12));
t736 = t757 * t760 + t762 * t786;
t679 = -0.2e1 * qJD(3) * t736 - t760 * t704 + t762 * t705;
t724 = t747 * t762 + t756 * t760;
t735 = t757 * t762 - t760 * t786;
t785 = t761 * t788;
t668 = (-t735 * t785 - t724) * pkin(9) + (t735 * t736 - t748) * pkin(3) + t679;
t680 = 0.2e1 * qJD(3) * t735 + t762 * t704 + t760 * t705;
t723 = -t747 * t760 + t756 * t762;
t725 = -pkin(3) * t785 - pkin(9) * t736;
t734 = t735 ^ 2;
t670 = -pkin(3) * t734 + pkin(9) * t723 + t725 * t785 + t680;
t766 = sin(qJ(4));
t771 = cos(qJ(4));
t650 = t771 * t668 - t766 * t670;
t717 = t735 * t771 - t736 * t766;
t688 = qJD(4) * t717 + t723 * t766 + t724 * t771;
t718 = t735 * t766 + t736 * t771;
t740 = qJDD(4) - t748;
t751 = qJD(4) - t785;
t647 = (t717 * t751 - t688) * pkin(10) + (t717 * t718 + t740) * pkin(4) + t650;
t651 = t766 * t668 + t771 * t670;
t687 = -qJD(4) * t718 + t723 * t771 - t724 * t766;
t708 = pkin(4) * t751 - pkin(10) * t718;
t716 = t717 ^ 2;
t649 = -pkin(4) * t716 + pkin(10) * t687 - t708 * t751 + t651;
t765 = sin(qJ(5));
t770 = cos(qJ(5));
t644 = t765 * t647 + t770 * t649;
t698 = t717 * t765 + t718 * t770;
t664 = -qJD(5) * t698 + t687 * t770 - t688 * t765;
t697 = t717 * t770 - t718 * t765;
t681 = -mrSges(6,1) * t697 + mrSges(6,2) * t698;
t750 = qJD(5) + t751;
t690 = mrSges(6,1) * t750 - mrSges(6,3) * t698;
t739 = qJDD(5) + t740;
t682 = -pkin(5) * t697 - pkin(11) * t698;
t749 = t750 ^ 2;
t642 = -pkin(5) * t749 + pkin(11) * t739 + t682 * t697 + t644;
t714 = -g(3) * t794 + t743 * t792 - t767 * t744;
t703 = -t756 * pkin(2) - t755 * qJ(3) + t745 * t786 + qJDD(3) - t714;
t683 = -t723 * pkin(3) - t734 * pkin(9) + t736 * t725 + t703;
t656 = -t687 * pkin(4) - t716 * pkin(10) + t718 * t708 + t683;
t665 = qJD(5) * t697 + t687 * t765 + t688 * t770;
t645 = t656 + (-t697 * t750 - t665) * pkin(11) + (t698 * t750 - t664) * pkin(5);
t764 = sin(qJ(6));
t769 = cos(qJ(6));
t639 = -t642 * t764 + t645 * t769;
t684 = -t698 * t764 + t750 * t769;
t654 = qJD(6) * t684 + t665 * t769 + t739 * t764;
t659 = qJDD(6) - t664;
t685 = t698 * t769 + t750 * t764;
t671 = -mrSges(7,1) * t684 + mrSges(7,2) * t685;
t696 = qJD(6) - t697;
t672 = -mrSges(7,2) * t696 + mrSges(7,3) * t684;
t637 = m(7) * t639 + mrSges(7,1) * t659 - mrSges(7,3) * t654 - t671 * t685 + t672 * t696;
t640 = t642 * t769 + t645 * t764;
t653 = -qJD(6) * t685 - t665 * t764 + t739 * t769;
t673 = mrSges(7,1) * t696 - mrSges(7,3) * t685;
t638 = m(7) * t640 - mrSges(7,2) * t659 + mrSges(7,3) * t653 + t671 * t684 - t673 * t696;
t780 = -t637 * t764 + t769 * t638;
t628 = m(6) * t644 - mrSges(6,2) * t739 + mrSges(6,3) * t664 + t681 * t697 - t690 * t750 + t780;
t643 = t647 * t770 - t649 * t765;
t689 = -mrSges(6,2) * t750 + mrSges(6,3) * t697;
t641 = -pkin(5) * t739 - pkin(11) * t749 + t682 * t698 - t643;
t777 = -m(7) * t641 + t653 * mrSges(7,1) - mrSges(7,2) * t654 + t684 * t672 - t673 * t685;
t633 = m(6) * t643 + mrSges(6,1) * t739 - mrSges(6,3) * t665 - t681 * t698 + t689 * t750 + t777;
t622 = t765 * t628 + t770 * t633;
t699 = -mrSges(5,1) * t717 + mrSges(5,2) * t718;
t706 = -mrSges(5,2) * t751 + mrSges(5,3) * t717;
t620 = m(5) * t650 + mrSges(5,1) * t740 - mrSges(5,3) * t688 - t699 * t718 + t706 * t751 + t622;
t707 = mrSges(5,1) * t751 - mrSges(5,3) * t718;
t781 = t770 * t628 - t633 * t765;
t621 = m(5) * t651 - mrSges(5,2) * t740 + mrSges(5,3) * t687 + t699 * t717 - t707 * t751 + t781;
t614 = t771 * t620 + t766 * t621;
t719 = -mrSges(4,1) * t735 + mrSges(4,2) * t736;
t721 = mrSges(4,2) * t785 + mrSges(4,3) * t735;
t612 = m(4) * t679 - mrSges(4,1) * t748 - mrSges(4,3) * t724 - t719 * t736 - t721 * t785 + t614;
t722 = -mrSges(4,1) * t785 - mrSges(4,3) * t736;
t782 = -t620 * t766 + t771 * t621;
t613 = m(4) * t680 + mrSges(4,2) * t748 + mrSges(4,3) * t723 + t719 * t735 + t722 * t785 + t782;
t783 = -t612 * t760 + t762 * t613;
t604 = m(3) * t715 - mrSges(3,2) * t756 + mrSges(3,3) * t748 - t741 * t757 + t746 * t785 + t783;
t607 = t762 * t612 + t760 * t613;
t729 = -t761 * t743 - t796;
t742 = -mrSges(3,2) * t757 + mrSges(3,3) * t785;
t606 = m(3) * t729 - t748 * mrSges(3,1) + t747 * mrSges(3,2) + (t741 * t767 - t742 * t772) * t789 + t607;
t629 = t769 * t637 + t764 * t638;
t779 = m(6) * t656 - t664 * mrSges(6,1) + t665 * mrSges(6,2) - t697 * t689 + t698 * t690 + t629;
t776 = m(5) * t683 - t687 * mrSges(5,1) + t688 * mrSges(5,2) - t717 * t706 + t718 * t707 + t779;
t775 = -m(4) * t703 + t723 * mrSges(4,1) - t724 * mrSges(4,2) + t735 * t721 - t736 * t722 - t776;
t625 = m(3) * t714 + t756 * mrSges(3,1) - t747 * mrSges(3,3) + t757 * t742 - t746 * t786 + t775;
t595 = t604 * t793 - t606 * t761 + t625 * t792;
t593 = m(2) * t752 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t774 + t595;
t599 = t772 * t604 - t625 * t767;
t598 = m(2) * t753 - mrSges(2,1) * t774 - qJDD(1) * mrSges(2,2) + t599;
t791 = t773 * t593 + t768 * t598;
t594 = t604 * t795 + t763 * t606 + t625 * t794;
t784 = -t593 * t768 + t773 * t598;
t660 = Ifges(7,5) * t685 + Ifges(7,6) * t684 + Ifges(7,3) * t696;
t662 = Ifges(7,1) * t685 + Ifges(7,4) * t684 + Ifges(7,5) * t696;
t630 = -mrSges(7,1) * t641 + mrSges(7,3) * t640 + Ifges(7,4) * t654 + Ifges(7,2) * t653 + Ifges(7,6) * t659 - t660 * t685 + t662 * t696;
t661 = Ifges(7,4) * t685 + Ifges(7,2) * t684 + Ifges(7,6) * t696;
t631 = mrSges(7,2) * t641 - mrSges(7,3) * t639 + Ifges(7,1) * t654 + Ifges(7,4) * t653 + Ifges(7,5) * t659 + t660 * t684 - t661 * t696;
t674 = Ifges(6,5) * t698 + Ifges(6,6) * t697 + Ifges(6,3) * t750;
t675 = Ifges(6,4) * t698 + Ifges(6,2) * t697 + Ifges(6,6) * t750;
t615 = mrSges(6,2) * t656 - mrSges(6,3) * t643 + Ifges(6,1) * t665 + Ifges(6,4) * t664 + Ifges(6,5) * t739 - pkin(11) * t629 - t630 * t764 + t631 * t769 + t674 * t697 - t675 * t750;
t676 = Ifges(6,1) * t698 + Ifges(6,4) * t697 + Ifges(6,5) * t750;
t616 = -mrSges(6,1) * t656 - mrSges(7,1) * t639 + mrSges(7,2) * t640 + mrSges(6,3) * t644 + Ifges(6,4) * t665 - Ifges(7,5) * t654 + Ifges(6,2) * t664 + Ifges(6,6) * t739 - Ifges(7,6) * t653 - Ifges(7,3) * t659 - pkin(5) * t629 - t661 * t685 + t662 * t684 - t674 * t698 + t676 * t750;
t691 = Ifges(5,5) * t718 + Ifges(5,6) * t717 + Ifges(5,3) * t751;
t693 = Ifges(5,1) * t718 + Ifges(5,4) * t717 + Ifges(5,5) * t751;
t600 = -mrSges(5,1) * t683 + mrSges(5,3) * t651 + Ifges(5,4) * t688 + Ifges(5,2) * t687 + Ifges(5,6) * t740 - pkin(4) * t779 + pkin(10) * t781 + t765 * t615 + t770 * t616 - t718 * t691 + t751 * t693;
t692 = Ifges(5,4) * t718 + Ifges(5,2) * t717 + Ifges(5,6) * t751;
t608 = mrSges(5,2) * t683 - mrSges(5,3) * t650 + Ifges(5,1) * t688 + Ifges(5,4) * t687 + Ifges(5,5) * t740 - pkin(10) * t622 + t615 * t770 - t616 * t765 + t691 * t717 - t692 * t751;
t709 = Ifges(4,5) * t736 + Ifges(4,6) * t735 - Ifges(4,3) * t785;
t711 = Ifges(4,1) * t736 + Ifges(4,4) * t735 - Ifges(4,5) * t785;
t589 = -mrSges(4,1) * t703 + mrSges(4,3) * t680 + Ifges(4,4) * t724 + Ifges(4,2) * t723 - Ifges(4,6) * t748 - pkin(3) * t776 + pkin(9) * t782 + t771 * t600 + t766 * t608 - t736 * t709 - t711 * t785;
t710 = Ifges(4,4) * t736 + Ifges(4,2) * t735 - Ifges(4,6) * t785;
t590 = mrSges(4,2) * t703 - mrSges(4,3) * t679 + Ifges(4,1) * t724 + Ifges(4,4) * t723 - Ifges(4,5) * t748 - pkin(9) * t614 - t600 * t766 + t608 * t771 + t709 * t735 + t710 * t785;
t726 = Ifges(3,3) * t757 + (Ifges(3,5) * t767 + Ifges(3,6) * t772) * t789;
t727 = Ifges(3,6) * t757 + (Ifges(3,4) * t767 + Ifges(3,2) * t772) * t789;
t588 = mrSges(3,2) * t729 - mrSges(3,3) * t714 + Ifges(3,1) * t747 + Ifges(3,4) * t748 + Ifges(3,5) * t756 - qJ(3) * t607 - t589 * t760 + t590 * t762 + t726 * t785 - t727 * t757;
t728 = Ifges(3,5) * t757 + (Ifges(3,1) * t767 + Ifges(3,4) * t772) * t789;
t591 = t735 * t711 - t736 * t710 - Ifges(6,3) * t739 - Ifges(5,3) * t740 - mrSges(6,1) * t643 + Ifges(3,4) * t747 + mrSges(5,2) * t651 + mrSges(3,3) * t715 + t717 * t693 - t718 * t692 - mrSges(4,1) * t679 + mrSges(4,2) * t680 - t726 * t786 - Ifges(5,6) * t687 - Ifges(5,5) * t688 - pkin(3) * t614 + t697 * t676 - t698 * t675 - Ifges(4,6) * t723 - Ifges(4,5) * t724 - mrSges(3,1) * t729 + mrSges(6,2) * t644 - Ifges(6,6) * t664 - Ifges(6,5) * t665 - pkin(11) * t780 - pkin(4) * t622 - mrSges(5,1) * t650 - pkin(5) * t777 + Ifges(3,6) * t756 + t757 * t728 - t764 * t631 - t769 * t630 - pkin(2) * t607 + (Ifges(3,2) + Ifges(4,3)) * t748;
t778 = pkin(8) * t599 + t588 * t767 + t591 * t772;
t587 = Ifges(3,5) * t747 + Ifges(3,6) * t748 + Ifges(3,3) * t756 + mrSges(3,1) * t714 - mrSges(3,2) * t715 + t760 * t590 + t762 * t589 + pkin(2) * t775 + qJ(3) * t783 + (t727 * t767 - t728 * t772) * t789;
t586 = -mrSges(2,2) * g(3) - mrSges(2,3) * t752 + Ifges(2,5) * qJDD(1) - t774 * Ifges(2,6) + t772 * t588 - t767 * t591 + (-t594 * t761 - t595 * t763) * pkin(8);
t585 = mrSges(2,1) * g(3) + mrSges(2,3) * t753 + t774 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t594 - t761 * t587 + t763 * t778;
t1 = [-m(1) * g(1) + t784; -m(1) * g(2) + t791; (-m(1) - m(2)) * g(3) + t594; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t791 - t768 * t585 + t773 * t586; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t784 + t773 * t585 + t768 * t586; -mrSges(1,1) * g(2) + mrSges(2,1) * t752 + mrSges(1,2) * g(1) - mrSges(2,2) * t753 + Ifges(2,3) * qJDD(1) + pkin(1) * t595 + t763 * t587 + t761 * t778;];
tauB  = t1;
