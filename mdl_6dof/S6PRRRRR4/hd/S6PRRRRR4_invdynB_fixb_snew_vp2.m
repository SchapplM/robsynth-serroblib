% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRR4
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:25:40
% EndTime: 2019-05-05 11:26:26
% DurationCPUTime: 44.84s
% Computational Cost: add. (806756->353), mult. (1665953->467), div. (0->0), fcn. (1356893->16), ass. (0->156)
t703 = sin(pkin(13));
t706 = cos(pkin(13));
t694 = g(1) * t703 - g(2) * t706;
t695 = -g(1) * t706 - g(2) * t703;
t702 = -g(3) + qJDD(1);
t713 = sin(qJ(2));
t708 = cos(pkin(6));
t718 = cos(qJ(2));
t736 = t708 * t718;
t705 = sin(pkin(6));
t739 = t705 * t718;
t665 = t694 * t736 - t695 * t713 + t702 * t739;
t704 = sin(pkin(7));
t719 = qJD(2) ^ 2;
t662 = pkin(9) * t704 * t719 + qJDD(2) * pkin(2) + t665;
t737 = t708 * t713;
t740 = t705 * t713;
t666 = t694 * t737 + t718 * t695 + t702 * t740;
t732 = qJDD(2) * t704;
t663 = -pkin(2) * t719 + pkin(9) * t732 + t666;
t679 = -t694 * t705 + t702 * t708;
t707 = cos(pkin(7));
t712 = sin(qJ(3));
t717 = cos(qJ(3));
t633 = -t712 * t663 + (t662 * t707 + t679 * t704) * t717;
t701 = qJD(2) * t707 + qJD(3);
t733 = qJD(2) * t717;
t730 = t704 * t733;
t683 = -mrSges(4,2) * t701 + mrSges(4,3) * t730;
t734 = qJD(2) * t704;
t684 = (-mrSges(4,1) * t717 + mrSges(4,2) * t712) * t734;
t686 = (qJD(3) * t733 + qJDD(2) * t712) * t704;
t700 = qJDD(2) * t707 + qJDD(3);
t685 = (-pkin(3) * t717 - pkin(10) * t712) * t734;
t699 = t701 ^ 2;
t731 = t712 * t734;
t623 = -t700 * pkin(3) - t699 * pkin(10) + t685 * t731 - t633;
t711 = sin(qJ(4));
t716 = cos(qJ(4));
t678 = t701 * t711 + t716 * t731;
t654 = -qJD(4) * t678 - t686 * t711 + t700 * t716;
t677 = t701 * t716 - t711 * t731;
t655 = qJD(4) * t677 + t686 * t716 + t700 * t711;
t693 = qJD(4) - t730;
t667 = -mrSges(5,2) * t693 + mrSges(5,3) * t677;
t668 = mrSges(5,1) * t693 - mrSges(5,3) * t678;
t738 = t707 * t712;
t741 = t704 * t712;
t634 = t662 * t738 + t717 * t663 + t679 * t741;
t624 = -pkin(3) * t699 + pkin(10) * t700 + t685 * t730 + t634;
t675 = t707 * t679;
t687 = -qJD(3) * t731 + t717 * t732;
t631 = -t687 * pkin(3) - t686 * pkin(10) + t675 + (-t662 + (pkin(3) * t712 - pkin(10) * t717) * t701 * qJD(2)) * t704;
t612 = -t711 * t624 + t716 * t631;
t681 = qJDD(4) - t687;
t609 = (t677 * t693 - t655) * pkin(11) + (t677 * t678 + t681) * pkin(4) + t612;
t613 = t716 * t624 + t711 * t631;
t669 = pkin(4) * t693 - pkin(11) * t678;
t676 = t677 ^ 2;
t611 = -pkin(4) * t676 + pkin(11) * t654 - t669 * t693 + t613;
t710 = sin(qJ(5));
t715 = cos(qJ(5));
t606 = t710 * t609 + t715 * t611;
t660 = t677 * t715 - t678 * t710;
t661 = t677 * t710 + t678 * t715;
t642 = -pkin(5) * t660 - pkin(12) * t661;
t680 = qJDD(5) + t681;
t692 = qJD(5) + t693;
t691 = t692 ^ 2;
t604 = -pkin(5) * t691 + pkin(12) * t680 + t642 * t660 + t606;
t614 = -t654 * pkin(4) - t676 * pkin(11) + t678 * t669 + t623;
t627 = -qJD(5) * t661 + t654 * t715 - t655 * t710;
t628 = qJD(5) * t660 + t654 * t710 + t655 * t715;
t607 = (-t660 * t692 - t628) * pkin(12) + (t661 * t692 - t627) * pkin(5) + t614;
t709 = sin(qJ(6));
t714 = cos(qJ(6));
t601 = -t604 * t709 + t607 * t714;
t644 = -t661 * t709 + t692 * t714;
t617 = qJD(6) * t644 + t628 * t714 + t680 * t709;
t626 = qJDD(6) - t627;
t645 = t661 * t714 + t692 * t709;
t632 = -mrSges(7,1) * t644 + mrSges(7,2) * t645;
t657 = qJD(6) - t660;
t635 = -mrSges(7,2) * t657 + mrSges(7,3) * t644;
t599 = m(7) * t601 + mrSges(7,1) * t626 - mrSges(7,3) * t617 - t632 * t645 + t635 * t657;
t602 = t604 * t714 + t607 * t709;
t616 = -qJD(6) * t645 - t628 * t709 + t680 * t714;
t636 = mrSges(7,1) * t657 - mrSges(7,3) * t645;
t600 = m(7) * t602 - mrSges(7,2) * t626 + mrSges(7,3) * t616 + t632 * t644 - t636 * t657;
t591 = t714 * t599 + t709 * t600;
t646 = -mrSges(6,2) * t692 + mrSges(6,3) * t660;
t647 = mrSges(6,1) * t692 - mrSges(6,3) * t661;
t721 = m(6) * t614 - t627 * mrSges(6,1) + mrSges(6,2) * t628 - t660 * t646 + t647 * t661 + t591;
t720 = -m(5) * t623 + t654 * mrSges(5,1) - mrSges(5,2) * t655 + t677 * t667 - t668 * t678 - t721;
t587 = m(4) * t633 + mrSges(4,1) * t700 - mrSges(4,3) * t686 + t683 * t701 - t684 * t731 + t720;
t742 = t587 * t717;
t682 = mrSges(4,1) * t701 - mrSges(4,3) * t731;
t641 = -mrSges(6,1) * t660 + mrSges(6,2) * t661;
t726 = -t599 * t709 + t714 * t600;
t590 = m(6) * t606 - mrSges(6,2) * t680 + mrSges(6,3) * t627 + t641 * t660 - t647 * t692 + t726;
t605 = t609 * t715 - t611 * t710;
t603 = -pkin(5) * t680 - pkin(12) * t691 + t642 * t661 - t605;
t722 = -m(7) * t603 + t616 * mrSges(7,1) - mrSges(7,2) * t617 + t644 * t635 - t636 * t645;
t595 = m(6) * t605 + mrSges(6,1) * t680 - mrSges(6,3) * t628 - t641 * t661 + t646 * t692 + t722;
t584 = t710 * t590 + t715 * t595;
t664 = -mrSges(5,1) * t677 + mrSges(5,2) * t678;
t582 = m(5) * t612 + mrSges(5,1) * t681 - mrSges(5,3) * t655 - t664 * t678 + t667 * t693 + t584;
t727 = t715 * t590 - t595 * t710;
t583 = m(5) * t613 - mrSges(5,2) * t681 + mrSges(5,3) * t654 + t664 * t677 - t668 * t693 + t727;
t728 = -t582 * t711 + t716 * t583;
t573 = m(4) * t634 - mrSges(4,2) * t700 + mrSges(4,3) * t687 - t682 * t701 + t684 * t730 + t728;
t576 = t716 * t582 + t711 * t583;
t643 = -t704 * t662 + t675;
t575 = m(4) * t643 - t687 * mrSges(4,1) + t686 * mrSges(4,2) + (t682 * t712 - t683 * t717) * t734 + t576;
t562 = t573 * t738 - t575 * t704 + t707 * t742;
t558 = m(3) * t665 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t719 + t562;
t561 = t573 * t741 + t707 * t575 + t704 * t742;
t560 = m(3) * t679 + t561;
t569 = t717 * t573 - t587 * t712;
t568 = m(3) * t666 - mrSges(3,1) * t719 - qJDD(2) * mrSges(3,2) + t569;
t548 = t558 * t736 - t560 * t705 + t568 * t737;
t546 = m(2) * t694 + t548;
t554 = -t558 * t713 + t718 * t568;
t553 = m(2) * t695 + t554;
t735 = t706 * t546 + t703 * t553;
t547 = t558 * t739 + t708 * t560 + t568 * t740;
t729 = -t546 * t703 + t706 * t553;
t618 = Ifges(7,5) * t645 + Ifges(7,6) * t644 + Ifges(7,3) * t657;
t620 = Ifges(7,1) * t645 + Ifges(7,4) * t644 + Ifges(7,5) * t657;
t592 = -mrSges(7,1) * t603 + mrSges(7,3) * t602 + Ifges(7,4) * t617 + Ifges(7,2) * t616 + Ifges(7,6) * t626 - t618 * t645 + t620 * t657;
t619 = Ifges(7,4) * t645 + Ifges(7,2) * t644 + Ifges(7,6) * t657;
t593 = mrSges(7,2) * t603 - mrSges(7,3) * t601 + Ifges(7,1) * t617 + Ifges(7,4) * t616 + Ifges(7,5) * t626 + t618 * t644 - t619 * t657;
t637 = Ifges(6,5) * t661 + Ifges(6,6) * t660 + Ifges(6,3) * t692;
t638 = Ifges(6,4) * t661 + Ifges(6,2) * t660 + Ifges(6,6) * t692;
t577 = mrSges(6,2) * t614 - mrSges(6,3) * t605 + Ifges(6,1) * t628 + Ifges(6,4) * t627 + Ifges(6,5) * t680 - pkin(12) * t591 - t592 * t709 + t593 * t714 + t637 * t660 - t638 * t692;
t639 = Ifges(6,1) * t661 + Ifges(6,4) * t660 + Ifges(6,5) * t692;
t578 = -mrSges(6,1) * t614 - mrSges(7,1) * t601 + mrSges(7,2) * t602 + mrSges(6,3) * t606 + Ifges(6,4) * t628 - Ifges(7,5) * t617 + Ifges(6,2) * t627 + Ifges(6,6) * t680 - Ifges(7,6) * t616 - Ifges(7,3) * t626 - pkin(5) * t591 - t619 * t645 + t620 * t644 - t637 * t661 + t639 * t692;
t648 = Ifges(5,5) * t678 + Ifges(5,6) * t677 + Ifges(5,3) * t693;
t650 = Ifges(5,1) * t678 + Ifges(5,4) * t677 + Ifges(5,5) * t693;
t563 = -mrSges(5,1) * t623 + mrSges(5,3) * t613 + Ifges(5,4) * t655 + Ifges(5,2) * t654 + Ifges(5,6) * t681 - pkin(4) * t721 + pkin(11) * t727 + t710 * t577 + t715 * t578 - t678 * t648 + t693 * t650;
t649 = Ifges(5,4) * t678 + Ifges(5,2) * t677 + Ifges(5,6) * t693;
t564 = mrSges(5,2) * t623 - mrSges(5,3) * t612 + Ifges(5,1) * t655 + Ifges(5,4) * t654 + Ifges(5,5) * t681 - pkin(11) * t584 + t577 * t715 - t578 * t710 + t648 * t677 - t649 * t693;
t672 = Ifges(4,6) * t701 + (Ifges(4,4) * t712 + Ifges(4,2) * t717) * t734;
t673 = Ifges(4,5) * t701 + (Ifges(4,1) * t712 + Ifges(4,4) * t717) * t734;
t549 = Ifges(4,5) * t686 + Ifges(4,6) * t687 + Ifges(4,3) * t700 + mrSges(4,1) * t633 - mrSges(4,2) * t634 + t711 * t564 + t716 * t563 + pkin(3) * t720 + pkin(10) * t728 + (t672 * t712 - t673 * t717) * t734;
t671 = Ifges(4,3) * t701 + (Ifges(4,5) * t712 + Ifges(4,6) * t717) * t734;
t550 = mrSges(4,2) * t643 - mrSges(4,3) * t633 + Ifges(4,1) * t686 + Ifges(4,4) * t687 + Ifges(4,5) * t700 - pkin(10) * t576 - t563 * t711 + t564 * t716 + t671 * t730 - t672 * t701;
t555 = -pkin(5) * t722 - t714 * t592 - t709 * t593 + Ifges(4,6) * t700 + t701 * t673 - Ifges(5,3) * t681 + Ifges(4,4) * t686 + Ifges(4,2) * t687 + t677 * t650 - t678 * t649 - Ifges(6,3) * t680 - pkin(12) * t726 - t671 * t731 + t660 * t639 - t661 * t638 - Ifges(5,6) * t654 - Ifges(5,5) * t655 - mrSges(4,1) * t643 + mrSges(4,3) * t634 - Ifges(6,6) * t627 - Ifges(6,5) * t628 + mrSges(5,2) * t613 + mrSges(6,2) * t606 - mrSges(5,1) * t612 - mrSges(6,1) * t605 - pkin(4) * t584 - pkin(3) * t576;
t723 = pkin(9) * t569 + t550 * t712 + t555 * t717;
t543 = -mrSges(3,1) * t679 + mrSges(3,3) * t666 + t719 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t561 - t704 * t549 + t707 * t723;
t544 = mrSges(3,2) * t679 - mrSges(3,3) * t665 + Ifges(3,5) * qJDD(2) - t719 * Ifges(3,6) + t717 * t550 - t712 * t555 + (-t561 * t704 - t562 * t707) * pkin(9);
t724 = pkin(8) * t554 + t543 * t718 + t544 * t713;
t542 = mrSges(3,1) * t665 - mrSges(3,2) * t666 + Ifges(3,3) * qJDD(2) + pkin(2) * t562 + t707 * t549 + t704 * t723;
t541 = mrSges(2,2) * t702 - mrSges(2,3) * t694 - t713 * t543 + t718 * t544 + (-t547 * t705 - t548 * t708) * pkin(8);
t540 = -mrSges(2,1) * t702 + mrSges(2,3) * t695 - pkin(1) * t547 - t705 * t542 + t708 * t724;
t1 = [-m(1) * g(1) + t729; -m(1) * g(2) + t735; -m(1) * g(3) + m(2) * t702 + t547; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t735 - t703 * t540 + t706 * t541; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t729 + t706 * t540 + t703 * t541; -mrSges(1,1) * g(2) + mrSges(2,1) * t694 + mrSges(1,2) * g(1) - mrSges(2,2) * t695 + pkin(1) * t548 + t708 * t542 + t705 * t724;];
tauB  = t1;
