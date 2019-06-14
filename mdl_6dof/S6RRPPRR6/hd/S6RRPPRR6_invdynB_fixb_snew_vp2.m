% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 10:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:53:05
% EndTime: 2019-05-06 10:53:16
% DurationCPUTime: 8.30s
% Computational Cost: add. (113346->363), mult. (258837->447), div. (0->0), fcn. (171019->10), ass. (0->140)
t789 = Ifges(3,1) + Ifges(4,1);
t784 = Ifges(3,4) - Ifges(4,5);
t783 = Ifges(3,5) + Ifges(4,4);
t788 = Ifges(3,2) + Ifges(4,3);
t782 = Ifges(3,6) - Ifges(4,6);
t787 = Ifges(3,3) + Ifges(4,2);
t786 = 2 * qJD(3);
t785 = mrSges(3,3) + mrSges(4,2);
t756 = cos(qJ(2));
t759 = qJD(1) ^ 2;
t781 = t756 ^ 2 * t759;
t753 = sin(qJ(1));
t757 = cos(qJ(1));
t731 = -t757 * g(1) - t753 * g(2);
t709 = -t759 * pkin(1) + qJDD(1) * pkin(7) + t731;
t752 = sin(qJ(2));
t688 = -t752 * g(3) + t756 * t709;
t720 = (-mrSges(3,1) * t756 + mrSges(3,2) * t752) * qJD(1);
t774 = qJD(1) * qJD(2);
t773 = t752 * t774;
t722 = t756 * qJDD(1) - t773;
t776 = qJD(1) * t752;
t726 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t776;
t718 = (-pkin(2) * t756 - qJ(3) * t752) * qJD(1);
t758 = qJD(2) ^ 2;
t775 = qJD(1) * t756;
t667 = -t758 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t786 + t718 * t775 + t688;
t719 = (-mrSges(4,1) * t756 - mrSges(4,3) * t752) * qJD(1);
t727 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t776;
t725 = -qJD(2) * pkin(3) - qJ(4) * t776;
t663 = -pkin(3) * t781 - t722 * qJ(4) + qJD(2) * t725 + t667;
t687 = -t756 * g(3) - t752 * t709;
t670 = -qJDD(2) * pkin(2) - t758 * qJ(3) + t718 * t776 + qJDD(3) - t687;
t772 = t756 * t774;
t721 = t752 * qJDD(1) + t772;
t664 = (-t721 + t772) * qJ(4) + (-t752 * t756 * t759 - qJDD(2)) * pkin(3) + t670;
t747 = sin(pkin(10));
t748 = cos(pkin(10));
t701 = (-t747 * t756 + t748 * t752) * qJD(1);
t632 = -0.2e1 * qJD(4) * t701 - t747 * t663 + t748 * t664;
t686 = t748 * t721 - t747 * t722;
t700 = (-t747 * t752 - t748 * t756) * qJD(1);
t629 = (-qJD(2) * t700 - t686) * pkin(8) + (t700 * t701 - qJDD(2)) * pkin(4) + t632;
t633 = 0.2e1 * qJD(4) * t700 + t748 * t663 + t747 * t664;
t685 = -t747 * t721 - t748 * t722;
t691 = -qJD(2) * pkin(4) - t701 * pkin(8);
t699 = t700 ^ 2;
t631 = -t699 * pkin(4) + t685 * pkin(8) + qJD(2) * t691 + t633;
t751 = sin(qJ(5));
t755 = cos(qJ(5));
t626 = t751 * t629 + t755 * t631;
t678 = t751 * t700 + t755 * t701;
t646 = -t678 * qJD(5) + t755 * t685 - t751 * t686;
t677 = t755 * t700 - t751 * t701;
t661 = -t677 * mrSges(6,1) + t678 * mrSges(6,2);
t741 = -qJD(2) + qJD(5);
t672 = t741 * mrSges(6,1) - t678 * mrSges(6,3);
t740 = -qJDD(2) + qJDD(5);
t662 = -t677 * pkin(5) - t678 * pkin(9);
t739 = t741 ^ 2;
t624 = -t739 * pkin(5) + t740 * pkin(9) + t677 * t662 + t626;
t730 = t753 * g(1) - t757 * g(2);
t708 = -qJDD(1) * pkin(1) - t759 * pkin(7) - t730;
t765 = -t722 * pkin(2) + t708 + (-t721 - t772) * qJ(3);
t650 = -pkin(2) * t773 + t722 * pkin(3) - qJ(4) * t781 + qJDD(4) - t765 + (t725 + t786) * t776;
t635 = -t685 * pkin(4) - t699 * pkin(8) + t701 * t691 + t650;
t647 = t677 * qJD(5) + t751 * t685 + t755 * t686;
t627 = (-t677 * t741 - t647) * pkin(9) + t635 + (t678 * t741 - t646) * pkin(5);
t750 = sin(qJ(6));
t754 = cos(qJ(6));
t621 = -t750 * t624 + t754 * t627;
t668 = -t750 * t678 + t754 * t741;
t638 = t668 * qJD(6) + t754 * t647 + t750 * t740;
t645 = qJDD(6) - t646;
t669 = t754 * t678 + t750 * t741;
t648 = -t668 * mrSges(7,1) + t669 * mrSges(7,2);
t673 = qJD(6) - t677;
t651 = -t673 * mrSges(7,2) + t668 * mrSges(7,3);
t619 = m(7) * t621 + t645 * mrSges(7,1) - t638 * mrSges(7,3) - t669 * t648 + t673 * t651;
t622 = t754 * t624 + t750 * t627;
t637 = -t669 * qJD(6) - t750 * t647 + t754 * t740;
t652 = t673 * mrSges(7,1) - t669 * mrSges(7,3);
t620 = m(7) * t622 - t645 * mrSges(7,2) + t637 * mrSges(7,3) + t668 * t648 - t673 * t652;
t767 = -t750 * t619 + t754 * t620;
t610 = m(6) * t626 - t740 * mrSges(6,2) + t646 * mrSges(6,3) + t677 * t661 - t741 * t672 + t767;
t625 = t755 * t629 - t751 * t631;
t671 = -t741 * mrSges(6,2) + t677 * mrSges(6,3);
t623 = -t740 * pkin(5) - t739 * pkin(9) + t678 * t662 - t625;
t763 = -m(7) * t623 + t637 * mrSges(7,1) - t638 * mrSges(7,2) + t668 * t651 - t669 * t652;
t615 = m(6) * t625 + t740 * mrSges(6,1) - t647 * mrSges(6,3) - t678 * t661 + t741 * t671 + t763;
t604 = t751 * t610 + t755 * t615;
t682 = -t700 * mrSges(5,1) + t701 * mrSges(5,2);
t689 = qJD(2) * mrSges(5,2) + t700 * mrSges(5,3);
t602 = m(5) * t632 - qJDD(2) * mrSges(5,1) - t686 * mrSges(5,3) - qJD(2) * t689 - t701 * t682 + t604;
t690 = -qJD(2) * mrSges(5,1) - t701 * mrSges(5,3);
t768 = t755 * t610 - t751 * t615;
t603 = m(5) * t633 + qJDD(2) * mrSges(5,2) + t685 * mrSges(5,3) + qJD(2) * t690 + t700 * t682 + t768;
t769 = -t747 * t602 + t748 * t603;
t764 = m(4) * t667 + qJDD(2) * mrSges(4,3) + qJD(2) * t727 + t719 * t775 + t769;
t596 = m(3) * t688 - qJDD(2) * mrSges(3,2) - qJD(2) * t726 + t720 * t775 + t785 * t722 + t764;
t728 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t775;
t598 = t748 * t602 + t747 * t603;
t729 = mrSges(4,2) * t775 + qJD(2) * mrSges(4,3);
t762 = -m(4) * t670 + qJDD(2) * mrSges(4,1) + qJD(2) * t729 - t598;
t597 = m(3) * t687 + qJDD(2) * mrSges(3,1) + qJD(2) * t728 - t785 * t721 + (-t719 - t720) * t776 + t762;
t770 = t756 * t596 - t752 * t597;
t590 = m(2) * t731 - t759 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t770;
t665 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t776 + t765;
t611 = t754 * t619 + t750 * t620;
t766 = m(6) * t635 - t646 * mrSges(6,1) + t647 * mrSges(6,2) - t677 * t671 + t678 * t672 + t611;
t761 = -m(5) * t650 + t685 * mrSges(5,1) - t686 * mrSges(5,2) + t700 * t689 - t701 * t690 - t766;
t607 = m(4) * t665 - t722 * mrSges(4,1) - t721 * mrSges(4,3) - t727 * t776 - t729 * t775 + t761;
t760 = -m(3) * t708 + t722 * mrSges(3,1) - t721 * mrSges(3,2) - t726 * t776 + t728 * t775 - t607;
t606 = m(2) * t730 + qJDD(1) * mrSges(2,1) - t759 * mrSges(2,2) + t760;
t780 = t753 * t590 + t757 * t606;
t591 = t752 * t596 + t756 * t597;
t779 = t787 * qJD(2) + (t783 * t752 + t782 * t756) * qJD(1);
t778 = -t782 * qJD(2) + (-t784 * t752 - t788 * t756) * qJD(1);
t777 = t783 * qJD(2) + (t789 * t752 + t784 * t756) * qJD(1);
t771 = t757 * t590 - t753 * t606;
t676 = Ifges(5,1) * t701 + Ifges(5,4) * t700 - Ifges(5,5) * qJD(2);
t675 = Ifges(5,4) * t701 + Ifges(5,2) * t700 - Ifges(5,6) * qJD(2);
t674 = Ifges(5,5) * t701 + Ifges(5,6) * t700 - Ifges(5,3) * qJD(2);
t655 = Ifges(6,1) * t678 + Ifges(6,4) * t677 + Ifges(6,5) * t741;
t654 = Ifges(6,4) * t678 + Ifges(6,2) * t677 + Ifges(6,6) * t741;
t653 = Ifges(6,5) * t678 + Ifges(6,6) * t677 + Ifges(6,3) * t741;
t641 = Ifges(7,1) * t669 + Ifges(7,4) * t668 + Ifges(7,5) * t673;
t640 = Ifges(7,4) * t669 + Ifges(7,2) * t668 + Ifges(7,6) * t673;
t639 = Ifges(7,5) * t669 + Ifges(7,6) * t668 + Ifges(7,3) * t673;
t613 = mrSges(7,2) * t623 - mrSges(7,3) * t621 + Ifges(7,1) * t638 + Ifges(7,4) * t637 + Ifges(7,5) * t645 + t668 * t639 - t673 * t640;
t612 = -mrSges(7,1) * t623 + mrSges(7,3) * t622 + Ifges(7,4) * t638 + Ifges(7,2) * t637 + Ifges(7,6) * t645 - t669 * t639 + t673 * t641;
t600 = -mrSges(6,1) * t635 - mrSges(7,1) * t621 + mrSges(7,2) * t622 + mrSges(6,3) * t626 + Ifges(6,4) * t647 - Ifges(7,5) * t638 + Ifges(6,2) * t646 + Ifges(6,6) * t740 - Ifges(7,6) * t637 - Ifges(7,3) * t645 - pkin(5) * t611 - t669 * t640 + t668 * t641 - t678 * t653 + t741 * t655;
t599 = mrSges(6,2) * t635 - mrSges(6,3) * t625 + Ifges(6,1) * t647 + Ifges(6,4) * t646 + Ifges(6,5) * t740 - pkin(9) * t611 - t750 * t612 + t754 * t613 + t677 * t653 - t741 * t654;
t592 = mrSges(5,2) * t650 - mrSges(5,3) * t632 + Ifges(5,1) * t686 + Ifges(5,4) * t685 - Ifges(5,5) * qJDD(2) - pkin(8) * t604 + qJD(2) * t675 + t755 * t599 - t751 * t600 + t700 * t674;
t587 = -mrSges(5,1) * t650 + mrSges(5,3) * t633 + Ifges(5,4) * t686 + Ifges(5,2) * t685 - Ifges(5,6) * qJDD(2) - pkin(4) * t766 + pkin(8) * t768 - qJD(2) * t676 + t751 * t599 + t755 * t600 - t701 * t674;
t586 = mrSges(3,2) * t708 + mrSges(4,2) * t670 - mrSges(3,3) * t687 - mrSges(4,3) * t665 - qJ(3) * t607 - qJ(4) * t598 + t778 * qJD(2) + t783 * qJDD(2) - t747 * t587 + t748 * t592 + t789 * t721 + t784 * t722 + t779 * t775;
t585 = -mrSges(3,1) * t708 - mrSges(4,1) * t665 + mrSges(4,2) * t667 + mrSges(3,3) * t688 - pkin(2) * t607 - pkin(3) * t761 - qJ(4) * t769 + t777 * qJD(2) + t782 * qJDD(2) - t748 * t587 - t747 * t592 + t784 * t721 + t788 * t722 - t779 * t776;
t584 = -qJ(3) * t764 - pkin(2) * t762 + pkin(5) * t763 + (t777 * t756 + (pkin(2) * t719 + t778) * t752) * qJD(1) + (-qJ(3) * mrSges(4,2) - t782) * t722 + (pkin(2) * mrSges(4,2) - t783) * t721 - mrSges(4,3) * t667 + t750 * t613 + t754 * t612 + mrSges(6,1) * t625 + pkin(9) * t767 - t677 * t655 + t678 * t654 - mrSges(6,2) * t626 + Ifges(6,6) * t646 + Ifges(6,5) * t647 + mrSges(5,1) * t632 - mrSges(5,2) * t633 - t700 * t676 + t701 * t675 + t759 * Ifges(2,5) + mrSges(2,3) * t731 + Ifges(6,3) * t740 + mrSges(2,1) * g(3) + pkin(4) * t604 + (-Ifges(5,3) - t787) * qJDD(2) + pkin(3) * t598 + Ifges(5,6) * t685 + Ifges(5,5) * t686 - mrSges(3,1) * t687 + mrSges(3,2) * t688 + mrSges(4,1) * t670 + Ifges(2,6) * qJDD(1) - pkin(1) * t591;
t583 = -mrSges(2,2) * g(3) - mrSges(2,3) * t730 + Ifges(2,5) * qJDD(1) - t759 * Ifges(2,6) - pkin(7) * t591 - t752 * t585 + t756 * t586;
t1 = [-m(1) * g(1) + t771; -m(1) * g(2) + t780; (-m(1) - m(2)) * g(3) + t591; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t780 + t757 * t583 - t753 * t584; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t771 + t753 * t583 + t757 * t584; -mrSges(1,1) * g(2) + mrSges(2,1) * t730 + mrSges(1,2) * g(1) - mrSges(2,2) * t731 + Ifges(2,3) * qJDD(1) + pkin(1) * t760 + pkin(7) * t770 + t756 * t585 + t752 * t586;];
tauB  = t1;
