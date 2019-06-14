% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:41:54
% EndTime: 2019-05-07 18:42:24
% DurationCPUTime: 25.77s
% Computational Cost: add. (440505->375), mult. (945818->475), div. (0->0), fcn. (752114->12), ass. (0->155)
t781 = Ifges(6,1) + Ifges(7,1);
t774 = Ifges(6,4) - Ifges(7,5);
t773 = Ifges(6,5) + Ifges(7,4);
t780 = Ifges(6,2) + Ifges(7,3);
t779 = -Ifges(7,2) - Ifges(6,3);
t772 = Ifges(6,6) - Ifges(7,6);
t733 = sin(pkin(6));
t737 = sin(qJ(2));
t741 = cos(qJ(2));
t758 = qJD(1) * qJD(2);
t722 = (-qJDD(1) * t741 + t737 * t758) * t733;
t778 = -2 * qJD(5);
t777 = pkin(8) * t733;
t734 = cos(pkin(6));
t776 = t734 * g(3);
t775 = -mrSges(6,3) - mrSges(7,2);
t771 = cos(pkin(11));
t770 = t733 * t737;
t769 = t733 * t741;
t768 = t734 * t737;
t767 = t734 * t741;
t738 = sin(qJ(1));
t742 = cos(qJ(1));
t725 = t738 * g(1) - t742 * g(2);
t743 = qJD(1) ^ 2;
t717 = qJDD(1) * pkin(1) + t743 * t777 + t725;
t726 = -t742 * g(1) - t738 * g(2);
t718 = -t743 * pkin(1) + qJDD(1) * t777 + t726;
t761 = t717 * t768 + t741 * t718;
t693 = -g(3) * t770 + t761;
t729 = t734 * qJD(1) + qJD(2);
t760 = qJD(1) * t733;
t756 = t737 * t760;
t715 = t729 * mrSges(3,1) - mrSges(3,3) * t756;
t719 = (-mrSges(3,1) * t741 + mrSges(3,2) * t737) * t760;
t728 = t734 * qJDD(1) + qJDD(2);
t720 = (-pkin(2) * t741 - pkin(9) * t737) * t760;
t727 = t729 ^ 2;
t759 = qJD(1) * t741;
t672 = -t727 * pkin(2) + t728 * pkin(9) + (-g(3) * t737 + t720 * t759) * t733 + t761;
t721 = (qJDD(1) * t737 + t741 * t758) * t733;
t673 = t722 * pkin(2) - t721 * pkin(9) - t776 + (-t717 + (pkin(2) * t737 - pkin(9) * t741) * t729 * qJD(1)) * t733;
t736 = sin(qJ(3));
t740 = cos(qJ(3));
t640 = t740 * t672 + t736 * t673;
t711 = t736 * t729 + t740 * t756;
t690 = -t711 * qJD(3) - t736 * t721 + t740 * t728;
t710 = t740 * t729 - t736 * t756;
t694 = -t710 * mrSges(4,1) + t711 * mrSges(4,2);
t755 = t733 * t759;
t724 = qJD(3) - t755;
t701 = t724 * mrSges(4,1) - t711 * mrSges(4,3);
t714 = qJDD(3) + t722;
t695 = -t710 * pkin(3) - t711 * pkin(10);
t723 = t724 ^ 2;
t630 = -t723 * pkin(3) + t714 * pkin(10) + t710 * t695 + t640;
t692 = -g(3) * t769 + t717 * t767 - t737 * t718;
t671 = -t728 * pkin(2) - t727 * pkin(9) + t720 * t756 - t692;
t691 = t710 * qJD(3) + t740 * t721 + t736 * t728;
t633 = (-t710 * t724 - t691) * pkin(10) + (t711 * t724 - t690) * pkin(3) + t671;
t735 = sin(qJ(4));
t739 = cos(qJ(4));
t625 = -t735 * t630 + t739 * t633;
t698 = -t735 * t711 + t739 * t724;
t658 = t698 * qJD(4) + t739 * t691 + t735 * t714;
t688 = qJDD(4) - t690;
t699 = t739 * t711 + t735 * t724;
t709 = qJD(4) - t710;
t622 = (t698 * t709 - t658) * qJ(5) + (t698 * t699 + t688) * pkin(4) + t625;
t626 = t739 * t630 + t735 * t633;
t657 = -t699 * qJD(4) - t735 * t691 + t739 * t714;
t680 = pkin(4) * t709 - qJ(5) * t699;
t697 = t698 ^ 2;
t624 = -pkin(4) * t697 + qJ(5) * t657 - t680 * t709 + t626;
t732 = sin(pkin(11));
t675 = -t771 * t698 + t732 * t699;
t618 = t732 * t622 + t771 * t624 + t675 * t778;
t637 = -t771 * t657 + t732 * t658;
t676 = t732 * t698 + t771 * t699;
t661 = mrSges(6,1) * t709 - mrSges(6,3) * t676;
t650 = pkin(5) * t675 - qJ(6) * t676;
t708 = t709 ^ 2;
t615 = -pkin(5) * t708 + qJ(6) * t688 + 0.2e1 * qJD(6) * t709 - t650 * t675 + t618;
t662 = -mrSges(7,1) * t709 + mrSges(7,2) * t676;
t757 = m(7) * t615 + t688 * mrSges(7,3) + t709 * t662;
t651 = mrSges(7,1) * t675 - mrSges(7,3) * t676;
t762 = -mrSges(6,1) * t675 - mrSges(6,2) * t676 - t651;
t608 = m(6) * t618 - mrSges(6,2) * t688 + t775 * t637 - t661 * t709 + t762 * t675 + t757;
t749 = t771 * t622 - t732 * t624;
t617 = t676 * t778 + t749;
t638 = t732 * t657 + t771 * t658;
t660 = -mrSges(6,2) * t709 - mrSges(6,3) * t675;
t616 = -t688 * pkin(5) - t708 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t650) * t676 - t749;
t659 = -mrSges(7,2) * t675 + mrSges(7,3) * t709;
t750 = -m(7) * t616 + t688 * mrSges(7,1) + t709 * t659;
t610 = m(6) * t617 + mrSges(6,1) * t688 + t775 * t638 + t660 * t709 + t762 * t676 + t750;
t603 = t732 * t608 + t771 * t610;
t677 = -mrSges(5,1) * t698 + mrSges(5,2) * t699;
t679 = -mrSges(5,2) * t709 + mrSges(5,3) * t698;
t601 = m(5) * t625 + mrSges(5,1) * t688 - mrSges(5,3) * t658 - t677 * t699 + t679 * t709 + t603;
t681 = mrSges(5,1) * t709 - mrSges(5,3) * t699;
t751 = t771 * t608 - t732 * t610;
t602 = m(5) * t626 - t688 * mrSges(5,2) + t657 * mrSges(5,3) + t698 * t677 - t709 * t681 + t751;
t752 = -t735 * t601 + t739 * t602;
t598 = m(4) * t640 - t714 * mrSges(4,2) + t690 * mrSges(4,3) + t710 * t694 - t724 * t701 + t752;
t639 = -t736 * t672 + t740 * t673;
t700 = -t724 * mrSges(4,2) + t710 * mrSges(4,3);
t629 = -t714 * pkin(3) - t723 * pkin(10) + t711 * t695 - t639;
t627 = -t657 * pkin(4) - t697 * qJ(5) + t699 * t680 + qJDD(5) + t629;
t620 = -0.2e1 * qJD(6) * t676 + (t675 * t709 - t638) * qJ(6) + (t676 * t709 + t637) * pkin(5) + t627;
t613 = m(7) * t620 + t637 * mrSges(7,1) - t638 * mrSges(7,3) + t675 * t659 - t676 * t662;
t746 = m(6) * t627 + t637 * mrSges(6,1) + t638 * mrSges(6,2) + t675 * t660 + t676 * t661 + t613;
t744 = -m(5) * t629 + t657 * mrSges(5,1) - t658 * mrSges(5,2) + t698 * t679 - t699 * t681 - t746;
t612 = m(4) * t639 + t714 * mrSges(4,1) - t691 * mrSges(4,3) - t711 * t694 + t724 * t700 + t744;
t753 = t740 * t598 - t736 * t612;
t589 = m(3) * t693 - t728 * mrSges(3,2) - t722 * mrSges(3,3) - t729 * t715 + t719 * t755 + t753;
t592 = t736 * t598 + t740 * t612;
t705 = -t733 * t717 - t776;
t716 = -t729 * mrSges(3,2) + mrSges(3,3) * t755;
t591 = m(3) * t705 + t722 * mrSges(3,1) + t721 * mrSges(3,2) + (t715 * t737 - t716 * t741) * t760 + t592;
t599 = t739 * t601 + t735 * t602;
t745 = -m(4) * t671 + t690 * mrSges(4,1) - t691 * mrSges(4,2) + t710 * t700 - t711 * t701 - t599;
t595 = m(3) * t692 + t728 * mrSges(3,1) - t721 * mrSges(3,3) + t729 * t716 - t719 * t756 + t745;
t578 = t589 * t768 - t733 * t591 + t595 * t767;
t576 = m(2) * t725 + qJDD(1) * mrSges(2,1) - t743 * mrSges(2,2) + t578;
t582 = t741 * t589 - t737 * t595;
t581 = m(2) * t726 - t743 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t582;
t766 = t742 * t576 + t738 * t581;
t765 = t675 * t780 - t676 * t774 - t709 * t772;
t764 = t675 * t772 - t676 * t773 + t709 * t779;
t763 = -t774 * t675 + t676 * t781 + t773 * t709;
t577 = t589 * t770 + t734 * t591 + t595 * t769;
t754 = -t738 * t576 + t742 * t581;
t604 = -mrSges(6,1) * t627 - mrSges(7,1) * t620 + mrSges(7,2) * t615 + mrSges(6,3) * t618 - pkin(5) * t613 - t637 * t780 + t774 * t638 + t764 * t676 + t772 * t688 + t763 * t709;
t605 = mrSges(6,2) * t627 + mrSges(7,2) * t616 - mrSges(6,3) * t617 - mrSges(7,3) * t620 - qJ(6) * t613 - t774 * t637 + t638 * t781 + t764 * t675 + t773 * t688 + t765 * t709;
t663 = Ifges(5,5) * t699 + Ifges(5,6) * t698 + Ifges(5,3) * t709;
t665 = Ifges(5,1) * t699 + Ifges(5,4) * t698 + Ifges(5,5) * t709;
t584 = -mrSges(5,1) * t629 + mrSges(5,3) * t626 + Ifges(5,4) * t658 + Ifges(5,2) * t657 + Ifges(5,6) * t688 - pkin(4) * t746 + qJ(5) * t751 + t771 * t604 + t732 * t605 - t699 * t663 + t709 * t665;
t664 = Ifges(5,4) * t699 + Ifges(5,2) * t698 + Ifges(5,6) * t709;
t585 = mrSges(5,2) * t629 - mrSges(5,3) * t625 + Ifges(5,1) * t658 + Ifges(5,4) * t657 + Ifges(5,5) * t688 - qJ(5) * t603 - t732 * t604 + t771 * t605 + t698 * t663 - t709 * t664;
t684 = Ifges(4,5) * t711 + Ifges(4,6) * t710 + Ifges(4,3) * t724;
t685 = Ifges(4,4) * t711 + Ifges(4,2) * t710 + Ifges(4,6) * t724;
t574 = mrSges(4,2) * t671 - mrSges(4,3) * t639 + Ifges(4,1) * t691 + Ifges(4,4) * t690 + Ifges(4,5) * t714 - pkin(10) * t599 - t735 * t584 + t739 * t585 + t710 * t684 - t724 * t685;
t686 = Ifges(4,1) * t711 + Ifges(4,4) * t710 + Ifges(4,5) * t724;
t583 = -qJ(6) * t757 - pkin(5) * t750 + (-Ifges(5,3) + t779) * t688 + (qJ(6) * mrSges(7,2) + t772) * t637 + (pkin(5) * mrSges(7,2) - t773) * t638 + (qJ(6) * t651 - t763) * t675 + (pkin(5) * t651 + t765) * t676 + mrSges(6,2) * t618 - mrSges(6,1) * t617 + mrSges(7,1) * t616 - mrSges(7,3) * t615 - pkin(4) * t603 - pkin(3) * t599 + mrSges(5,2) * t626 - mrSges(5,1) * t625 + t724 * t686 - t711 * t684 + Ifges(4,6) * t714 + t698 * t665 - t699 * t664 + Ifges(4,2) * t690 + Ifges(4,4) * t691 - mrSges(4,1) * t671 - Ifges(5,6) * t657 - Ifges(5,5) * t658 + mrSges(4,3) * t640;
t702 = Ifges(3,3) * t729 + (Ifges(3,5) * t737 + Ifges(3,6) * t741) * t760;
t703 = Ifges(3,6) * t729 + (Ifges(3,4) * t737 + Ifges(3,2) * t741) * t760;
t572 = mrSges(3,2) * t705 - mrSges(3,3) * t692 + Ifges(3,1) * t721 - Ifges(3,4) * t722 + Ifges(3,5) * t728 - pkin(9) * t592 + t740 * t574 - t736 * t583 + t702 * t755 - t729 * t703;
t704 = Ifges(3,5) * t729 + (Ifges(3,1) * t737 + Ifges(3,4) * t741) * t760;
t573 = Ifges(3,4) * t721 - Ifges(3,2) * t722 + Ifges(3,6) * t728 - t702 * t756 + t729 * t704 - mrSges(3,1) * t705 + mrSges(3,3) * t693 - Ifges(4,5) * t691 - Ifges(4,6) * t690 - Ifges(4,3) * t714 - t711 * t685 + t710 * t686 - mrSges(4,1) * t639 + mrSges(4,2) * t640 - t735 * t585 - t739 * t584 - pkin(3) * t744 - pkin(10) * t752 - pkin(2) * t592;
t747 = pkin(8) * t582 + t572 * t737 + t573 * t741;
t571 = Ifges(3,5) * t721 - Ifges(3,6) * t722 + Ifges(3,3) * t728 + mrSges(3,1) * t692 - mrSges(3,2) * t693 + t736 * t574 + t740 * t583 + pkin(2) * t745 + pkin(9) * t753 + (t703 * t737 - t704 * t741) * t760;
t570 = -mrSges(2,2) * g(3) - mrSges(2,3) * t725 + Ifges(2,5) * qJDD(1) - t743 * Ifges(2,6) + t741 * t572 - t737 * t573 + (-t577 * t733 - t578 * t734) * pkin(8);
t569 = mrSges(2,1) * g(3) + mrSges(2,3) * t726 + t743 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t577 - t733 * t571 + t734 * t747;
t1 = [-m(1) * g(1) + t754; -m(1) * g(2) + t766; (-m(1) - m(2)) * g(3) + t577; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t766 - t738 * t569 + t742 * t570; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t754 + t742 * t569 + t738 * t570; -mrSges(1,1) * g(2) + mrSges(2,1) * t725 + mrSges(1,2) * g(1) - mrSges(2,2) * t726 + Ifges(2,3) * qJDD(1) + pkin(1) * t578 + t734 * t571 + t747 * t733;];
tauB  = t1;
