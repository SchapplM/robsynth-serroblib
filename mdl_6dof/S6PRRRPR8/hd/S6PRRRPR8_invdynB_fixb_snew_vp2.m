% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 09:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:06:20
% EndTime: 2019-05-05 09:06:37
% DurationCPUTime: 16.62s
% Computational Cost: add. (280093->335), mult. (579240->424), div. (0->0), fcn. (450665->14), ass. (0->150)
t754 = Ifges(5,1) + Ifges(6,2);
t748 = Ifges(5,4) + Ifges(6,6);
t747 = Ifges(5,5) - Ifges(6,4);
t753 = -Ifges(5,2) - Ifges(6,3);
t746 = Ifges(5,6) - Ifges(6,5);
t752 = Ifges(5,3) + Ifges(6,1);
t703 = sin(pkin(7));
t710 = sin(qJ(3));
t713 = cos(qJ(3));
t731 = qJD(2) * qJD(3);
t685 = (-qJDD(2) * t713 + t710 * t731) * t703;
t702 = sin(pkin(12));
t705 = cos(pkin(12));
t694 = t702 * g(1) - t705 * g(2);
t695 = -t705 * g(1) - t702 * g(2);
t701 = -g(3) + qJDD(1);
t711 = sin(qJ(2));
t707 = cos(pkin(6));
t714 = cos(qJ(2));
t738 = t707 * t714;
t704 = sin(pkin(6));
t741 = t704 * t714;
t656 = t694 * t738 - t711 * t695 + t701 * t741;
t715 = qJD(2) ^ 2;
t749 = pkin(9) * t703;
t650 = qJDD(2) * pkin(2) + t715 * t749 + t656;
t739 = t707 * t711;
t742 = t704 * t711;
t657 = t694 * t739 + t714 * t695 + t701 * t742;
t651 = -t715 * pkin(2) + qJDD(2) * t749 + t657;
t677 = -t704 * t694 + t707 * t701;
t706 = cos(pkin(7));
t618 = -t710 * t651 + (t650 * t706 + t677 * t703) * t713;
t751 = -2 * qJD(5);
t750 = cos(qJ(4));
t700 = t706 * qJD(2) + qJD(3);
t732 = qJD(2) * t703;
t729 = t713 * t732;
t681 = -t700 * mrSges(4,2) + mrSges(4,3) * t729;
t682 = (-mrSges(4,1) * t713 + mrSges(4,2) * t710) * t732;
t684 = (qJDD(2) * t710 + t713 * t731) * t703;
t699 = t706 * qJDD(2) + qJDD(3);
t683 = (-pkin(3) * t713 - pkin(10) * t710) * t732;
t698 = t700 ^ 2;
t730 = t710 * t732;
t614 = -t699 * pkin(3) - t698 * pkin(10) + t683 * t730 - t618;
t709 = sin(qJ(4));
t676 = t709 * t700 + t750 * t730;
t645 = t676 * qJD(4) + t709 * t684 - t750 * t699;
t675 = -t750 * t700 + t709 * t730;
t646 = -t675 * qJD(4) + t750 * t684 + t709 * t699;
t692 = -qJD(4) + t729;
t660 = t692 * mrSges(5,2) - t675 * mrSges(5,3);
t661 = -t692 * mrSges(5,1) - t676 * mrSges(5,3);
t663 = t676 * mrSges(6,1) - t692 * mrSges(6,2);
t744 = t675 * t692;
t716 = (-t646 - t744) * qJ(5) + t614 + (-pkin(4) * t692 + t751) * t676;
t609 = t645 * pkin(4) + t716;
t662 = t675 * mrSges(6,1) + t692 * mrSges(6,3);
t740 = t706 * t710;
t743 = t703 * t710;
t619 = t650 * t740 + t713 * t651 + t677 * t743;
t615 = -t698 * pkin(3) + t699 * pkin(10) + t683 * t729 + t619;
t671 = t706 * t677;
t617 = t685 * pkin(3) - t684 * pkin(10) + t671 + (-t650 + (pkin(3) * t710 - pkin(10) * t713) * t700 * qJD(2)) * t703;
t610 = -t709 * t615 + t750 * t617;
t652 = t675 * pkin(4) - t676 * qJ(5);
t679 = qJDD(4) + t685;
t691 = t692 ^ 2;
t608 = -t679 * pkin(4) - t691 * qJ(5) + t676 * t652 + qJDD(5) - t610;
t603 = (t675 * t676 - t679) * pkin(11) + (t646 - t744) * pkin(5) + t608;
t664 = t676 * pkin(5) + t692 * pkin(11);
t674 = t675 ^ 2;
t606 = -t674 * pkin(5) - t676 * t664 + (pkin(4) + pkin(11)) * t645 + t716;
t708 = sin(qJ(6));
t712 = cos(qJ(6));
t601 = t712 * t603 - t708 * t606;
t658 = t712 * t675 + t708 * t692;
t622 = t658 * qJD(6) + t708 * t645 + t712 * t679;
t659 = t708 * t675 - t712 * t692;
t627 = -t658 * mrSges(7,1) + t659 * mrSges(7,2);
t673 = qJD(6) + t676;
t631 = -t673 * mrSges(7,2) + t658 * mrSges(7,3);
t642 = qJDD(6) + t646;
t599 = m(7) * t601 + t642 * mrSges(7,1) - t622 * mrSges(7,3) - t659 * t627 + t673 * t631;
t602 = t708 * t603 + t712 * t606;
t621 = -t659 * qJD(6) + t712 * t645 - t708 * t679;
t632 = t673 * mrSges(7,1) - t659 * mrSges(7,3);
t600 = m(7) * t602 - t642 * mrSges(7,2) + t621 * mrSges(7,3) + t658 * t627 - t673 * t632;
t736 = -t708 * t599 + t712 * t600;
t725 = -m(6) * t609 + t645 * mrSges(6,2) + t675 * t662 - t736;
t717 = -m(5) * t614 - t645 * mrSges(5,1) - t675 * t660 + (-t661 + t663) * t676 + (-mrSges(5,2) + mrSges(6,3)) * t646 + t725;
t587 = m(4) * t618 + t699 * mrSges(4,1) - t684 * mrSges(4,3) + t700 * t681 - t682 * t730 + t717;
t745 = t587 * t713;
t680 = t700 * mrSges(4,1) - mrSges(4,3) * t730;
t653 = t675 * mrSges(5,1) + t676 * mrSges(5,2);
t591 = t712 * t599 + t708 * t600;
t654 = -t675 * mrSges(6,2) - t676 * mrSges(6,3);
t720 = -m(6) * t608 - t646 * mrSges(6,1) - t676 * t654 - t591;
t589 = m(5) * t610 - t646 * mrSges(5,3) - t676 * t653 + (-t660 + t662) * t692 + (mrSges(5,1) - mrSges(6,2)) * t679 + t720;
t611 = t750 * t615 + t709 * t617;
t719 = -t691 * pkin(4) + t679 * qJ(5) - t675 * t652 + t611;
t607 = 0.2e1 * qJD(5) * t692 - t719;
t605 = -t645 * pkin(5) - t674 * pkin(11) + (t751 - t664) * t692 + t719;
t721 = -m(7) * t605 + t621 * mrSges(7,1) - t622 * mrSges(7,2) + t658 * t631 - t659 * t632;
t718 = -m(6) * t607 + t679 * mrSges(6,3) - t692 * t663 - t721;
t596 = m(5) * t611 - t679 * mrSges(5,2) + t692 * t661 + (-t653 - t654) * t675 + (-mrSges(5,3) - mrSges(6,1)) * t645 + t718;
t727 = -t709 * t589 + t750 * t596;
t581 = m(4) * t619 - t699 * mrSges(4,2) - t685 * mrSges(4,3) - t700 * t680 + t682 * t729 + t727;
t584 = t750 * t589 + t709 * t596;
t628 = -t703 * t650 + t671;
t583 = m(4) * t628 + t685 * mrSges(4,1) + t684 * mrSges(4,2) + (t680 * t710 - t681 * t713) * t732 + t584;
t570 = t581 * t740 - t703 * t583 + t706 * t745;
t566 = m(3) * t656 + qJDD(2) * mrSges(3,1) - t715 * mrSges(3,2) + t570;
t569 = t581 * t743 + t706 * t583 + t703 * t745;
t568 = m(3) * t677 + t569;
t577 = t713 * t581 - t710 * t587;
t576 = m(3) * t657 - t715 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t577;
t556 = t566 * t738 - t704 * t568 + t576 * t739;
t554 = m(2) * t694 + t556;
t562 = -t711 * t566 + t714 * t576;
t561 = m(2) * t695 + t562;
t737 = t705 * t554 + t702 * t561;
t735 = t746 * t675 - t747 * t676 + t752 * t692;
t734 = t753 * t675 + t748 * t676 - t746 * t692;
t733 = t748 * t675 - t754 * t676 + t747 * t692;
t555 = t566 * t741 + t707 * t568 + t576 * t742;
t728 = -t702 * t554 + t705 * t561;
t590 = -t646 * mrSges(6,3) - t676 * t663 - t725;
t623 = Ifges(7,5) * t659 + Ifges(7,6) * t658 + Ifges(7,3) * t673;
t625 = Ifges(7,1) * t659 + Ifges(7,4) * t658 + Ifges(7,5) * t673;
t592 = -mrSges(7,1) * t605 + mrSges(7,3) * t602 + Ifges(7,4) * t622 + Ifges(7,2) * t621 + Ifges(7,6) * t642 - t659 * t623 + t673 * t625;
t624 = Ifges(7,4) * t659 + Ifges(7,2) * t658 + Ifges(7,6) * t673;
t593 = mrSges(7,2) * t605 - mrSges(7,3) * t601 + Ifges(7,1) * t622 + Ifges(7,4) * t621 + Ifges(7,5) * t642 + t658 * t623 - t673 * t624;
t571 = -mrSges(5,1) * t614 - mrSges(6,1) * t607 + mrSges(6,2) * t609 + mrSges(5,3) * t611 - pkin(4) * t590 - pkin(5) * t721 - pkin(11) * t736 - t712 * t592 - t708 * t593 + t753 * t645 + t748 * t646 + t735 * t676 + t746 * t679 + t733 * t692;
t572 = mrSges(6,1) * t608 + mrSges(7,1) * t601 + mrSges(5,2) * t614 - mrSges(7,2) * t602 - mrSges(5,3) * t610 - mrSges(6,3) * t609 + Ifges(7,5) * t622 + Ifges(7,6) * t621 + Ifges(7,3) * t642 + pkin(5) * t591 - qJ(5) * t590 + t659 * t624 - t658 * t625 + t734 * t692 + t747 * t679 + t735 * t675 + t754 * t646 - t748 * t645;
t668 = Ifges(4,6) * t700 + (Ifges(4,4) * t710 + Ifges(4,2) * t713) * t732;
t669 = Ifges(4,5) * t700 + (Ifges(4,1) * t710 + Ifges(4,4) * t713) * t732;
t557 = Ifges(4,5) * t684 - Ifges(4,6) * t685 + Ifges(4,3) * t699 + mrSges(4,1) * t618 - mrSges(4,2) * t619 + t709 * t572 + t750 * t571 + pkin(3) * t717 + pkin(10) * t727 + (t668 * t710 - t669 * t713) * t732;
t667 = Ifges(4,3) * t700 + (Ifges(4,5) * t710 + Ifges(4,6) * t713) * t732;
t558 = mrSges(4,2) * t628 - mrSges(4,3) * t618 + Ifges(4,1) * t684 - Ifges(4,4) * t685 + Ifges(4,5) * t699 - pkin(10) * t584 - t709 * t571 + t750 * t572 + t667 * t729 - t700 * t668;
t563 = -t667 * t730 - pkin(4) * (t692 * t662 + t720) - t712 * t593 + t708 * t592 + Ifges(4,6) * t699 + t700 * t669 + Ifges(4,4) * t684 - Ifges(4,2) * t685 - mrSges(4,1) * t628 + mrSges(4,3) * t619 + mrSges(6,3) * t607 - mrSges(6,2) * t608 - mrSges(5,1) * t610 + mrSges(5,2) * t611 + pkin(11) * t591 - pkin(3) * t584 - qJ(5) * t718 + (pkin(4) * mrSges(6,2) - t752) * t679 - t734 * t676 + (qJ(5) * t654 + t733) * t675 - t747 * t646 + (qJ(5) * mrSges(6,1) + t746) * t645;
t722 = pkin(9) * t577 + t558 * t710 + t563 * t713;
t551 = -mrSges(3,1) * t677 + mrSges(3,3) * t657 + t715 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t569 - t703 * t557 + t722 * t706;
t552 = mrSges(3,2) * t677 - mrSges(3,3) * t656 + Ifges(3,5) * qJDD(2) - t715 * Ifges(3,6) + t713 * t558 - t710 * t563 + (-t569 * t703 - t570 * t706) * pkin(9);
t723 = pkin(8) * t562 + t551 * t714 + t552 * t711;
t550 = mrSges(3,1) * t656 - mrSges(3,2) * t657 + Ifges(3,3) * qJDD(2) + pkin(2) * t570 + t706 * t557 + t722 * t703;
t549 = mrSges(2,2) * t701 - mrSges(2,3) * t694 - t711 * t551 + t714 * t552 + (-t555 * t704 - t556 * t707) * pkin(8);
t548 = -mrSges(2,1) * t701 + mrSges(2,3) * t695 - pkin(1) * t555 - t704 * t550 + t723 * t707;
t1 = [-m(1) * g(1) + t728; -m(1) * g(2) + t737; -m(1) * g(3) + m(2) * t701 + t555; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t737 - t702 * t548 + t705 * t549; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t728 + t705 * t548 + t702 * t549; -mrSges(1,1) * g(2) + mrSges(2,1) * t694 + mrSges(1,2) * g(1) - mrSges(2,2) * t695 + pkin(1) * t556 + t707 * t550 + t723 * t704;];
tauB  = t1;
