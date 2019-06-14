% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 22:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:15:15
% EndTime: 2019-05-07 22:16:34
% DurationCPUTime: 56.64s
% Computational Cost: add. (959062->398), mult. (2052889->516), div. (0->0), fcn. (1677758->14), ass. (0->162)
t796 = cos(qJ(4));
t763 = cos(pkin(6));
t795 = t763 * g(3);
t761 = sin(pkin(6));
t767 = sin(qJ(2));
t794 = t761 * t767;
t771 = cos(qJ(2));
t793 = t761 * t771;
t792 = t763 * t767;
t791 = t763 * t771;
t768 = sin(qJ(1));
t772 = cos(qJ(1));
t752 = t768 * g(1) - g(2) * t772;
t773 = qJD(1) ^ 2;
t743 = pkin(8) * t761 * t773 + qJDD(1) * pkin(1) + t752;
t753 = -g(1) * t772 - g(2) * t768;
t786 = qJDD(1) * t761;
t744 = -pkin(1) * t773 + pkin(8) * t786 + t753;
t789 = t743 * t792 + t771 * t744;
t719 = -g(3) * t794 + t789;
t757 = qJD(1) * t763 + qJD(2);
t788 = qJD(1) * t761;
t785 = t767 * t788;
t741 = mrSges(3,1) * t757 - mrSges(3,3) * t785;
t745 = (-mrSges(3,1) * t771 + mrSges(3,2) * t767) * t788;
t748 = -qJD(2) * t785 + t771 * t786;
t756 = qJDD(1) * t763 + qJDD(2);
t746 = (-pkin(2) * t771 - pkin(9) * t767) * t788;
t755 = t757 ^ 2;
t787 = qJD(1) * t771;
t700 = -t755 * pkin(2) + t756 * pkin(9) + (-g(3) * t767 + t746 * t787) * t761 + t789;
t747 = (qJD(2) * t787 + qJDD(1) * t767) * t761;
t701 = -t748 * pkin(2) - t747 * pkin(9) - t795 + (-t743 + (pkin(2) * t767 - pkin(9) * t771) * t757 * qJD(1)) * t761;
t766 = sin(qJ(3));
t770 = cos(qJ(3));
t671 = -t766 * t700 + t770 * t701;
t735 = t757 * t770 - t766 * t785;
t717 = qJD(3) * t735 + t747 * t770 + t756 * t766;
t736 = t757 * t766 + t770 * t785;
t740 = qJDD(3) - t748;
t784 = t761 * t787;
t751 = qJD(3) - t784;
t658 = (t735 * t751 - t717) * pkin(10) + (t735 * t736 + t740) * pkin(3) + t671;
t672 = t770 * t700 + t766 * t701;
t716 = -qJD(3) * t736 - t747 * t766 + t756 * t770;
t726 = pkin(3) * t751 - pkin(10) * t736;
t734 = t735 ^ 2;
t665 = -pkin(3) * t734 + pkin(10) * t716 - t726 * t751 + t672;
t765 = sin(qJ(4));
t653 = t765 * t658 + t665 * t796;
t722 = t765 * t735 + t736 * t796;
t681 = qJD(4) * t722 - t716 * t796 + t717 * t765;
t721 = -t735 * t796 + t736 * t765;
t695 = mrSges(5,1) * t721 + mrSges(5,2) * t722;
t750 = qJD(4) + t751;
t708 = mrSges(5,1) * t750 - mrSges(5,3) * t722;
t739 = qJDD(4) + t740;
t694 = pkin(4) * t721 - qJ(5) * t722;
t749 = t750 ^ 2;
t648 = -pkin(4) * t749 + qJ(5) * t739 - t694 * t721 + t653;
t718 = -g(3) * t793 + t743 * t791 - t767 * t744;
t699 = -t756 * pkin(2) - t755 * pkin(9) + t746 * t785 - t718;
t667 = -t716 * pkin(3) - t734 * pkin(10) + t736 * t726 + t699;
t682 = -t721 * qJD(4) + t765 * t716 + t717 * t796;
t651 = (t721 * t750 - t682) * qJ(5) + (t722 * t750 + t681) * pkin(4) + t667;
t760 = sin(pkin(12));
t762 = cos(pkin(12));
t706 = t722 * t762 + t750 * t760;
t643 = -0.2e1 * qJD(5) * t706 - t760 * t648 + t762 * t651;
t674 = t682 * t762 + t739 * t760;
t705 = -t722 * t760 + t750 * t762;
t641 = (t705 * t721 - t674) * pkin(11) + (t705 * t706 + t681) * pkin(5) + t643;
t644 = 0.2e1 * qJD(5) * t705 + t762 * t648 + t760 * t651;
t673 = -t682 * t760 + t739 * t762;
t689 = pkin(5) * t721 - pkin(11) * t706;
t704 = t705 ^ 2;
t642 = -pkin(5) * t704 + pkin(11) * t673 - t689 * t721 + t644;
t764 = sin(qJ(6));
t769 = cos(qJ(6));
t639 = t641 * t769 - t642 * t764;
t684 = t705 * t769 - t706 * t764;
t656 = qJD(6) * t684 + t673 * t764 + t674 * t769;
t685 = t705 * t764 + t706 * t769;
t666 = -mrSges(7,1) * t684 + mrSges(7,2) * t685;
t720 = qJD(6) + t721;
t669 = -mrSges(7,2) * t720 + mrSges(7,3) * t684;
t680 = qJDD(6) + t681;
t637 = m(7) * t639 + mrSges(7,1) * t680 - mrSges(7,3) * t656 - t666 * t685 + t669 * t720;
t640 = t641 * t764 + t642 * t769;
t655 = -qJD(6) * t685 + t673 * t769 - t674 * t764;
t670 = mrSges(7,1) * t720 - mrSges(7,3) * t685;
t638 = m(7) * t640 - mrSges(7,2) * t680 + mrSges(7,3) * t655 + t666 * t684 - t670 * t720;
t629 = t769 * t637 + t764 * t638;
t686 = -mrSges(6,1) * t705 + mrSges(6,2) * t706;
t687 = -mrSges(6,2) * t721 + mrSges(6,3) * t705;
t627 = m(6) * t643 + mrSges(6,1) * t681 - mrSges(6,3) * t674 - t686 * t706 + t687 * t721 + t629;
t688 = mrSges(6,1) * t721 - mrSges(6,3) * t706;
t779 = -t637 * t764 + t769 * t638;
t628 = m(6) * t644 - mrSges(6,2) * t681 + mrSges(6,3) * t673 + t686 * t705 - t688 * t721 + t779;
t780 = -t627 * t760 + t762 * t628;
t622 = m(5) * t653 - mrSges(5,2) * t739 - mrSges(5,3) * t681 - t695 * t721 - t708 * t750 + t780;
t652 = t658 * t796 - t765 * t665;
t707 = -mrSges(5,2) * t750 - mrSges(5,3) * t721;
t647 = -t739 * pkin(4) - t749 * qJ(5) + t722 * t694 + qJDD(5) - t652;
t645 = -t673 * pkin(5) - t704 * pkin(11) + t706 * t689 + t647;
t777 = m(7) * t645 - t655 * mrSges(7,1) + mrSges(7,2) * t656 - t684 * t669 + t670 * t685;
t775 = -m(6) * t647 + t673 * mrSges(6,1) - mrSges(6,2) * t674 + t705 * t687 - t688 * t706 - t777;
t633 = m(5) * t652 + mrSges(5,1) * t739 - mrSges(5,3) * t682 - t695 * t722 + t707 * t750 + t775;
t614 = t765 * t622 + t796 * t633;
t723 = -mrSges(4,1) * t735 + mrSges(4,2) * t736;
t724 = -mrSges(4,2) * t751 + mrSges(4,3) * t735;
t612 = m(4) * t671 + mrSges(4,1) * t740 - mrSges(4,3) * t717 - t723 * t736 + t724 * t751 + t614;
t725 = mrSges(4,1) * t751 - mrSges(4,3) * t736;
t781 = t622 * t796 - t633 * t765;
t613 = m(4) * t672 - mrSges(4,2) * t740 + mrSges(4,3) * t716 + t723 * t735 - t725 * t751 + t781;
t782 = -t612 * t766 + t770 * t613;
t604 = m(3) * t719 - mrSges(3,2) * t756 + mrSges(3,3) * t748 - t741 * t757 + t745 * t784 + t782;
t607 = t770 * t612 + t766 * t613;
t730 = -t761 * t743 - t795;
t742 = -mrSges(3,2) * t757 + mrSges(3,3) * t784;
t606 = m(3) * t730 - t748 * mrSges(3,1) + t747 * mrSges(3,2) + (t741 * t767 - t742 * t771) * t788 + t607;
t623 = t762 * t627 + t760 * t628;
t776 = m(5) * t667 + t681 * mrSges(5,1) + mrSges(5,2) * t682 + t721 * t707 + t708 * t722 + t623;
t774 = -m(4) * t699 + t716 * mrSges(4,1) - mrSges(4,2) * t717 + t735 * t724 - t725 * t736 - t776;
t619 = m(3) * t718 + mrSges(3,1) * t756 - mrSges(3,3) * t747 + t742 * t757 - t745 * t785 + t774;
t595 = t604 * t792 - t606 * t761 + t619 * t791;
t593 = m(2) * t752 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t773 + t595;
t599 = t771 * t604 - t619 * t767;
t598 = m(2) * t753 - mrSges(2,1) * t773 - qJDD(1) * mrSges(2,2) + t599;
t790 = t772 * t593 + t768 * t598;
t594 = t604 * t794 + t763 * t606 + t619 * t793;
t783 = -t593 * t768 + t772 * t598;
t659 = Ifges(7,5) * t685 + Ifges(7,6) * t684 + Ifges(7,3) * t720;
t661 = Ifges(7,1) * t685 + Ifges(7,4) * t684 + Ifges(7,5) * t720;
t630 = -mrSges(7,1) * t645 + mrSges(7,3) * t640 + Ifges(7,4) * t656 + Ifges(7,2) * t655 + Ifges(7,6) * t680 - t659 * t685 + t661 * t720;
t660 = Ifges(7,4) * t685 + Ifges(7,2) * t684 + Ifges(7,6) * t720;
t631 = mrSges(7,2) * t645 - mrSges(7,3) * t639 + Ifges(7,1) * t656 + Ifges(7,4) * t655 + Ifges(7,5) * t680 + t659 * t684 - t660 * t720;
t675 = Ifges(6,5) * t706 + Ifges(6,6) * t705 + Ifges(6,3) * t721;
t677 = Ifges(6,1) * t706 + Ifges(6,4) * t705 + Ifges(6,5) * t721;
t615 = -mrSges(6,1) * t647 + mrSges(6,3) * t644 + Ifges(6,4) * t674 + Ifges(6,2) * t673 + Ifges(6,6) * t681 - pkin(5) * t777 + pkin(11) * t779 + t769 * t630 + t764 * t631 - t706 * t675 + t721 * t677;
t676 = Ifges(6,4) * t706 + Ifges(6,2) * t705 + Ifges(6,6) * t721;
t616 = mrSges(6,2) * t647 - mrSges(6,3) * t643 + Ifges(6,1) * t674 + Ifges(6,4) * t673 + Ifges(6,5) * t681 - pkin(11) * t629 - t630 * t764 + t631 * t769 + t675 * t705 - t676 * t721;
t690 = Ifges(5,5) * t722 - Ifges(5,6) * t721 + Ifges(5,3) * t750;
t691 = Ifges(5,4) * t722 - Ifges(5,2) * t721 + Ifges(5,6) * t750;
t600 = mrSges(5,2) * t667 - mrSges(5,3) * t652 + Ifges(5,1) * t682 - Ifges(5,4) * t681 + Ifges(5,5) * t739 - qJ(5) * t623 - t615 * t760 + t616 * t762 - t690 * t721 - t691 * t750;
t692 = Ifges(5,1) * t722 - Ifges(5,4) * t721 + Ifges(5,5) * t750;
t608 = Ifges(5,4) * t682 + Ifges(5,6) * t739 - t722 * t690 + t750 * t692 - mrSges(5,1) * t667 + mrSges(5,3) * t653 - Ifges(6,5) * t674 - Ifges(6,6) * t673 - t706 * t676 + t705 * t677 - mrSges(6,1) * t643 + mrSges(6,2) * t644 - Ifges(7,5) * t656 - Ifges(7,6) * t655 - Ifges(7,3) * t680 - t685 * t660 + t684 * t661 - mrSges(7,1) * t639 + mrSges(7,2) * t640 - pkin(5) * t629 - pkin(4) * t623 + (-Ifges(5,2) - Ifges(6,3)) * t681;
t710 = Ifges(4,5) * t736 + Ifges(4,6) * t735 + Ifges(4,3) * t751;
t712 = Ifges(4,1) * t736 + Ifges(4,4) * t735 + Ifges(4,5) * t751;
t589 = -mrSges(4,1) * t699 + mrSges(4,3) * t672 + Ifges(4,4) * t717 + Ifges(4,2) * t716 + Ifges(4,6) * t740 - pkin(3) * t776 + pkin(10) * t781 + t765 * t600 + t608 * t796 - t736 * t710 + t751 * t712;
t711 = Ifges(4,4) * t736 + Ifges(4,2) * t735 + Ifges(4,6) * t751;
t591 = mrSges(4,2) * t699 - mrSges(4,3) * t671 + Ifges(4,1) * t717 + Ifges(4,4) * t716 + Ifges(4,5) * t740 - pkin(10) * t614 + t600 * t796 - t765 * t608 + t735 * t710 - t751 * t711;
t727 = Ifges(3,3) * t757 + (Ifges(3,5) * t767 + Ifges(3,6) * t771) * t788;
t728 = Ifges(3,6) * t757 + (Ifges(3,4) * t767 + Ifges(3,2) * t771) * t788;
t588 = mrSges(3,2) * t730 - mrSges(3,3) * t718 + Ifges(3,1) * t747 + Ifges(3,4) * t748 + Ifges(3,5) * t756 - pkin(9) * t607 - t589 * t766 + t591 * t770 + t727 * t784 - t728 * t757;
t729 = Ifges(3,5) * t757 + (Ifges(3,1) * t767 + Ifges(3,4) * t771) * t788;
t590 = -pkin(3) * t614 - t727 * t785 - pkin(4) * t775 - qJ(5) * t780 - pkin(2) * t607 + t757 * t729 - t760 * t616 - t762 * t615 + Ifges(3,6) * t756 - Ifges(4,3) * t740 + Ifges(3,4) * t747 + Ifges(3,2) * t748 + t735 * t712 - t736 * t711 - Ifges(5,3) * t739 - t722 * t691 - mrSges(3,1) * t730 - t721 * t692 - Ifges(4,6) * t716 - Ifges(4,5) * t717 + mrSges(3,3) * t719 + Ifges(5,6) * t681 - Ifges(5,5) * t682 - mrSges(4,1) * t671 + mrSges(4,2) * t672 + mrSges(5,2) * t653 - mrSges(5,1) * t652;
t778 = pkin(8) * t599 + t588 * t767 + t590 * t771;
t587 = Ifges(3,5) * t747 + Ifges(3,6) * t748 + Ifges(3,3) * t756 + mrSges(3,1) * t718 - mrSges(3,2) * t719 + t766 * t591 + t770 * t589 + pkin(2) * t774 + pkin(9) * t782 + (t728 * t767 - t729 * t771) * t788;
t586 = -mrSges(2,2) * g(3) - mrSges(2,3) * t752 + Ifges(2,5) * qJDD(1) - t773 * Ifges(2,6) + t771 * t588 - t767 * t590 + (-t594 * t761 - t595 * t763) * pkin(8);
t585 = mrSges(2,1) * g(3) + mrSges(2,3) * t753 + t773 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t594 - t761 * t587 + t763 * t778;
t1 = [-m(1) * g(1) + t783; -m(1) * g(2) + t790; (-m(1) - m(2)) * g(3) + t594; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t790 - t768 * t585 + t772 * t586; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t783 + t772 * t585 + t768 * t586; -mrSges(1,1) * g(2) + mrSges(2,1) * t752 + mrSges(1,2) * g(1) - mrSges(2,2) * t753 + Ifges(2,3) * qJDD(1) + pkin(1) * t595 + t763 * t587 + t761 * t778;];
tauB  = t1;
