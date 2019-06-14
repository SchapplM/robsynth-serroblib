% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:28:23
% EndTime: 2019-05-05 09:28:35
% DurationCPUTime: 10.30s
% Computational Cost: add. (168171->320), mult. (328668->400), div. (0->0), fcn. (234407->12), ass. (0->134)
t713 = Ifges(6,1) + Ifges(7,1);
t708 = Ifges(6,4) + Ifges(7,4);
t707 = Ifges(6,5) + Ifges(7,5);
t712 = Ifges(6,2) + Ifges(7,2);
t711 = Ifges(6,6) + Ifges(7,6);
t710 = Ifges(6,3) + Ifges(7,3);
t672 = sin(qJ(4));
t673 = sin(qJ(3));
t676 = cos(qJ(4));
t677 = cos(qJ(3));
t646 = (t672 * t673 - t676 * t677) * qJD(2);
t667 = sin(pkin(11));
t669 = cos(pkin(11));
t656 = g(1) * t667 - g(2) * t669;
t657 = -g(1) * t669 - g(2) * t667;
t666 = -g(3) + qJDD(1);
t668 = sin(pkin(6));
t670 = cos(pkin(6));
t674 = sin(qJ(2));
t678 = cos(qJ(2));
t628 = -t674 * t657 + (t656 * t670 + t666 * t668) * t678;
t709 = -mrSges(6,2) - mrSges(7,2);
t679 = qJD(2) ^ 2;
t683 = -qJDD(2) * pkin(2) - t628;
t623 = -t679 * pkin(8) + t683;
t695 = qJD(2) * qJD(3);
t692 = t677 * t695;
t654 = qJDD(2) * t673 + t692;
t655 = qJDD(2) * t677 - t673 * t695;
t697 = qJD(2) * t673;
t658 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t697;
t696 = qJD(2) * t677;
t659 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t696;
t703 = t670 * t674;
t704 = t668 * t674;
t629 = t656 * t703 + t678 * t657 + t666 * t704;
t624 = -pkin(2) * t679 + qJDD(2) * pkin(8) + t629;
t640 = -t656 * t668 + t666 * t670;
t606 = -t673 * t624 + t677 * t640;
t589 = (-t654 + t692) * pkin(9) + (t673 * t677 * t679 + qJDD(3)) * pkin(3) + t606;
t607 = t677 * t624 + t673 * t640;
t661 = qJD(3) * pkin(3) - pkin(9) * t697;
t665 = t677 ^ 2;
t590 = -pkin(3) * t665 * t679 + pkin(9) * t655 - qJD(3) * t661 + t607;
t585 = t672 * t589 + t676 * t590;
t647 = (t672 * t677 + t673 * t676) * qJD(2);
t632 = pkin(4) * t646 - pkin(10) * t647;
t664 = qJD(3) + qJD(4);
t662 = t664 ^ 2;
t663 = qJDD(3) + qJDD(4);
t580 = -pkin(4) * t662 + pkin(10) * t663 - t632 * t646 + t585;
t597 = -t655 * pkin(3) + t661 * t697 + (-pkin(9) * t665 - pkin(8)) * t679 + t683;
t616 = -qJD(4) * t647 - t654 * t672 + t655 * t676;
t617 = -qJD(4) * t646 + t654 * t676 + t655 * t672;
t583 = (t646 * t664 - t617) * pkin(10) + (t647 * t664 - t616) * pkin(4) + t597;
t671 = sin(qJ(5));
t675 = cos(qJ(5));
t575 = -t671 * t580 + t675 * t583;
t634 = -t647 * t671 + t664 * t675;
t595 = qJD(5) * t634 + t617 * t675 + t663 * t671;
t635 = t647 * t675 + t664 * t671;
t609 = -mrSges(7,1) * t634 + mrSges(7,2) * t635;
t610 = -mrSges(6,1) * t634 + mrSges(6,2) * t635;
t614 = qJDD(5) - t616;
t641 = qJD(5) + t646;
t619 = -mrSges(6,2) * t641 + mrSges(6,3) * t634;
t572 = -0.2e1 * qJD(6) * t635 + (t634 * t641 - t595) * qJ(6) + (t634 * t635 + t614) * pkin(5) + t575;
t618 = -mrSges(7,2) * t641 + mrSges(7,3) * t634;
t694 = m(7) * t572 + t614 * mrSges(7,1) + t641 * t618;
t564 = m(6) * t575 + t614 * mrSges(6,1) + t641 * t619 + (-t609 - t610) * t635 + (-mrSges(6,3) - mrSges(7,3)) * t595 + t694;
t576 = t675 * t580 + t671 * t583;
t594 = -qJD(5) * t635 - t617 * t671 + t663 * t675;
t620 = pkin(5) * t641 - qJ(6) * t635;
t633 = t634 ^ 2;
t574 = -pkin(5) * t633 + qJ(6) * t594 + 0.2e1 * qJD(6) * t634 - t620 * t641 + t576;
t693 = m(7) * t574 + t594 * mrSges(7,3) + t634 * t609;
t621 = mrSges(7,1) * t641 - mrSges(7,3) * t635;
t698 = -mrSges(6,1) * t641 + mrSges(6,3) * t635 - t621;
t567 = m(6) * t576 + t594 * mrSges(6,3) + t634 * t610 + t709 * t614 + t698 * t641 + t693;
t562 = t675 * t564 + t671 * t567;
t638 = -mrSges(5,2) * t664 - mrSges(5,3) * t646;
t639 = mrSges(5,1) * t664 - mrSges(5,3) * t647;
t682 = m(5) * t597 - t616 * mrSges(5,1) + mrSges(5,2) * t617 + t646 * t638 + t639 * t647 + t562;
t680 = -m(4) * t623 + t655 * mrSges(4,1) - mrSges(4,2) * t654 - t658 * t697 + t659 * t696 - t682;
t557 = m(3) * t628 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t679 + t680;
t705 = t557 * t678;
t631 = mrSges(5,1) * t646 + mrSges(5,2) * t647;
t688 = -t564 * t671 + t675 * t567;
t560 = m(5) * t585 - mrSges(5,2) * t663 + mrSges(5,3) * t616 - t631 * t646 - t639 * t664 + t688;
t584 = t589 * t676 - t672 * t590;
t579 = -pkin(4) * t663 - pkin(10) * t662 + t647 * t632 - t584;
t577 = -pkin(5) * t594 - qJ(6) * t633 + t620 * t635 + qJDD(6) + t579;
t687 = m(7) * t577 - t594 * mrSges(7,1) - t634 * t618;
t681 = -m(6) * t579 + t594 * mrSges(6,1) + t709 * t595 + t634 * t619 + t698 * t635 - t687;
t569 = m(5) * t584 + t663 * mrSges(5,1) - t617 * mrSges(5,3) - t647 * t631 + t664 * t638 + t681;
t553 = t672 * t560 + t676 * t569;
t653 = (-mrSges(4,1) * t677 + mrSges(4,2) * t673) * qJD(2);
t551 = m(4) * t606 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t654 + qJD(3) * t659 - t653 * t697 + t553;
t689 = t676 * t560 - t569 * t672;
t552 = m(4) * t607 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t655 - qJD(3) * t658 + t653 * t696 + t689;
t690 = -t551 * t673 + t677 * t552;
t543 = m(3) * t629 - mrSges(3,1) * t679 - qJDD(2) * mrSges(3,2) + t690;
t546 = t677 * t551 + t673 * t552;
t545 = m(3) * t640 + t546;
t533 = t543 * t703 - t545 * t668 + t670 * t705;
t531 = m(2) * t656 + t533;
t538 = t678 * t543 - t557 * t674;
t537 = m(2) * t657 + t538;
t702 = t669 * t531 + t667 * t537;
t701 = t634 * t711 + t635 * t707 + t641 * t710;
t700 = -t634 * t712 - t635 * t708 - t641 * t711;
t699 = t708 * t634 + t635 * t713 + t707 * t641;
t532 = t543 * t704 + t670 * t545 + t668 * t705;
t691 = -t531 * t667 + t669 * t537;
t554 = -mrSges(6,1) * t579 + mrSges(6,3) * t576 - mrSges(7,1) * t577 + mrSges(7,3) * t574 - pkin(5) * t687 + qJ(6) * t693 + (-qJ(6) * t621 + t699) * t641 + (-pkin(5) * t621 - t701) * t635 + (-mrSges(7,2) * qJ(6) + t711) * t614 + (-mrSges(7,2) * pkin(5) + t708) * t595 + t712 * t594;
t570 = -t595 * mrSges(7,3) - t635 * t609 + t694;
t561 = mrSges(6,2) * t579 + mrSges(7,2) * t577 - mrSges(6,3) * t575 - mrSges(7,3) * t572 - qJ(6) * t570 + t708 * t594 + t595 * t713 + t707 * t614 + t701 * t634 + t700 * t641;
t625 = Ifges(5,5) * t647 - Ifges(5,6) * t646 + Ifges(5,3) * t664;
t626 = Ifges(5,4) * t647 - Ifges(5,2) * t646 + Ifges(5,6) * t664;
t539 = mrSges(5,2) * t597 - mrSges(5,3) * t584 + Ifges(5,1) * t617 + Ifges(5,4) * t616 + Ifges(5,5) * t663 - pkin(10) * t562 - t554 * t671 + t561 * t675 - t625 * t646 - t626 * t664;
t627 = Ifges(5,1) * t647 - Ifges(5,4) * t646 + Ifges(5,5) * t664;
t547 = -mrSges(5,1) * t597 - mrSges(6,1) * t575 - mrSges(7,1) * t572 + mrSges(6,2) * t576 + mrSges(7,2) * t574 + mrSges(5,3) * t585 + Ifges(5,4) * t617 + Ifges(5,2) * t616 + Ifges(5,6) * t663 - pkin(4) * t562 - pkin(5) * t570 - t647 * t625 + t664 * t627 + t700 * t635 + t699 * t634 - t710 * t614 - t707 * t595 - t711 * t594;
t643 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t673 + Ifges(4,6) * t677) * qJD(2);
t645 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t673 + Ifges(4,4) * t677) * qJD(2);
t529 = -mrSges(4,1) * t623 + mrSges(4,3) * t607 + Ifges(4,4) * t654 + Ifges(4,2) * t655 + Ifges(4,6) * qJDD(3) - pkin(3) * t682 + pkin(9) * t689 + qJD(3) * t645 + t672 * t539 + t676 * t547 - t643 * t697;
t644 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t673 + Ifges(4,2) * t677) * qJD(2);
t534 = mrSges(4,2) * t623 - mrSges(4,3) * t606 + Ifges(4,1) * t654 + Ifges(4,4) * t655 + Ifges(4,5) * qJDD(3) - pkin(9) * t553 - qJD(3) * t644 + t539 * t676 - t547 * t672 + t643 * t696;
t527 = mrSges(3,2) * t640 - mrSges(3,3) * t628 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t679 - pkin(8) * t546 - t529 * t673 + t534 * t677;
t528 = -pkin(2) * t546 + mrSges(3,3) * t629 - mrSges(3,1) * t640 + Ifges(3,6) * qJDD(2) - pkin(3) * t553 - Ifges(4,5) * t654 - Ifges(4,6) * t655 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t606 + mrSges(4,2) * t607 - t671 * t561 - t675 * t554 - pkin(4) * t681 - pkin(10) * t688 - Ifges(5,5) * t617 - Ifges(5,6) * t616 - Ifges(5,3) * t663 - mrSges(5,1) * t584 + mrSges(5,2) * t585 - t647 * t626 - t646 * t627 + t679 * Ifges(3,5) + (-t644 * t673 + t645 * t677) * qJD(2);
t684 = pkin(7) * t538 + t527 * t674 + t528 * t678;
t526 = mrSges(3,1) * t628 - mrSges(3,2) * t629 + Ifges(3,3) * qJDD(2) + pkin(2) * t680 + pkin(8) * t690 + t677 * t529 + t673 * t534;
t525 = mrSges(2,2) * t666 - mrSges(2,3) * t656 + t678 * t527 - t674 * t528 + (-t532 * t668 - t533 * t670) * pkin(7);
t524 = -mrSges(2,1) * t666 + mrSges(2,3) * t657 - pkin(1) * t532 - t668 * t526 + t684 * t670;
t1 = [-m(1) * g(1) + t691; -m(1) * g(2) + t702; -m(1) * g(3) + m(2) * t666 + t532; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t702 - t667 * t524 + t669 * t525; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t691 + t669 * t524 + t667 * t525; -mrSges(1,1) * g(2) + mrSges(2,1) * t656 + mrSges(1,2) * g(1) - mrSges(2,2) * t657 + pkin(1) * t533 + t670 * t526 + t684 * t668;];
tauB  = t1;
