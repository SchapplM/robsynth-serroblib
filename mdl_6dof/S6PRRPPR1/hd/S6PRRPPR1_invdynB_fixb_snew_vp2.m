% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-05-05 02:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:31:06
% EndTime: 2019-05-05 02:31:26
% DurationCPUTime: 18.90s
% Computational Cost: add. (301217->340), mult. (670352->442), div. (0->0), fcn. (483635->14), ass. (0->140)
t708 = -2 * qJD(4);
t672 = sin(pkin(10));
t675 = cos(pkin(10));
t661 = g(1) * t672 - g(2) * t675;
t662 = -g(1) * t675 - g(2) * t672;
t669 = -g(3) + qJDD(1);
t682 = cos(qJ(2));
t676 = cos(pkin(6));
t679 = sin(qJ(2));
t704 = t676 * t679;
t673 = sin(pkin(6));
t705 = t673 * t679;
t625 = t661 * t704 + t682 * t662 + t669 * t705;
t684 = qJD(2) ^ 2;
t616 = -pkin(2) * t684 + qJDD(2) * pkin(8) + t625;
t642 = -t661 * t673 + t669 * t676;
t678 = sin(qJ(3));
t681 = cos(qJ(3));
t602 = -t678 * t616 + t681 * t642;
t699 = qJD(2) * qJD(3);
t698 = t681 * t699;
t659 = qJDD(2) * t678 + t698;
t597 = (-t659 + t698) * qJ(4) + (t678 * t681 * t684 + qJDD(3)) * pkin(3) + t602;
t603 = t681 * t616 + t678 * t642;
t660 = qJDD(2) * t681 - t678 * t699;
t702 = qJD(2) * t678;
t663 = qJD(3) * pkin(3) - qJ(4) * t702;
t668 = t681 ^ 2;
t598 = -pkin(3) * t668 * t684 + qJ(4) * t660 - qJD(3) * t663 + t603;
t671 = sin(pkin(11));
t707 = cos(pkin(11));
t647 = (t671 * t681 + t678 * t707) * qJD(2);
t581 = t597 * t707 - t671 * t598 + t647 * t708;
t624 = -t679 * t662 + (t661 * t676 + t669 * t673) * t682;
t688 = -qJDD(2) * pkin(2) - t624;
t615 = -t684 * pkin(8) + t688;
t664 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t702;
t701 = qJD(2) * t681;
t665 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t701;
t646 = t671 * t702 - t707 * t701;
t582 = t671 * t597 + t707 * t598 + t646 * t708;
t627 = pkin(4) * t646 - qJ(5) * t647;
t683 = qJD(3) ^ 2;
t580 = -pkin(4) * t683 + qJDD(3) * qJ(5) - t627 * t646 + t582;
t599 = -t660 * pkin(3) + qJDD(4) + t663 * t702 + (-qJ(4) * t668 - pkin(8)) * t684 + t688;
t631 = t659 * t671 - t707 * t660;
t632 = t659 * t707 + t671 * t660;
t585 = (qJD(3) * t646 - t632) * qJ(5) + (qJD(3) * t647 + t631) * pkin(4) + t599;
t670 = sin(pkin(12));
t674 = cos(pkin(12));
t637 = qJD(3) * t670 + t647 * t674;
t575 = -0.2e1 * qJD(5) * t637 - t670 * t580 + t674 * t585;
t620 = qJDD(3) * t670 + t632 * t674;
t636 = qJD(3) * t674 - t647 * t670;
t573 = (t636 * t646 - t620) * pkin(9) + (t636 * t637 + t631) * pkin(5) + t575;
t576 = 0.2e1 * qJD(5) * t636 + t674 * t580 + t670 * t585;
t617 = pkin(5) * t646 - pkin(9) * t637;
t619 = qJDD(3) * t674 - t632 * t670;
t635 = t636 ^ 2;
t574 = -pkin(5) * t635 + pkin(9) * t619 - t617 * t646 + t576;
t677 = sin(qJ(6));
t680 = cos(qJ(6));
t571 = t573 * t680 - t574 * t677;
t608 = t636 * t680 - t637 * t677;
t588 = qJD(6) * t608 + t619 * t677 + t620 * t680;
t609 = t636 * t677 + t637 * t680;
t593 = -mrSges(7,1) * t608 + mrSges(7,2) * t609;
t645 = qJD(6) + t646;
t600 = -mrSges(7,2) * t645 + mrSges(7,3) * t608;
t630 = qJDD(6) + t631;
t569 = m(7) * t571 + mrSges(7,1) * t630 - mrSges(7,3) * t588 - t593 * t609 + t600 * t645;
t572 = t573 * t677 + t574 * t680;
t587 = -qJD(6) * t609 + t619 * t680 - t620 * t677;
t601 = mrSges(7,1) * t645 - mrSges(7,3) * t609;
t570 = m(7) * t572 - mrSges(7,2) * t630 + mrSges(7,3) * t587 + t593 * t608 - t601 * t645;
t561 = t680 * t569 + t677 * t570;
t610 = -mrSges(6,1) * t636 + mrSges(6,2) * t637;
t613 = -mrSges(6,2) * t646 + mrSges(6,3) * t636;
t559 = m(6) * t575 + mrSges(6,1) * t631 - mrSges(6,3) * t620 - t610 * t637 + t613 * t646 + t561;
t614 = mrSges(6,1) * t646 - mrSges(6,3) * t637;
t693 = -t569 * t677 + t680 * t570;
t560 = m(6) * t576 - mrSges(6,2) * t631 + mrSges(6,3) * t619 + t610 * t636 - t614 * t646 + t693;
t555 = t674 * t559 + t670 * t560;
t640 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t646;
t641 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t647;
t687 = m(5) * t599 + t631 * mrSges(5,1) + mrSges(5,2) * t632 + t646 * t640 + t641 * t647 + t555;
t685 = -m(4) * t615 + t660 * mrSges(4,1) - mrSges(4,2) * t659 - t664 * t702 + t665 * t701 - t687;
t551 = m(3) * t624 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t684 + t685;
t706 = t551 * t682;
t628 = mrSges(5,1) * t646 + mrSges(5,2) * t647;
t694 = -t559 * t670 + t674 * t560;
t554 = m(5) * t582 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t631 - qJD(3) * t641 - t628 * t646 + t694;
t579 = -qJDD(3) * pkin(4) - t683 * qJ(5) + t647 * t627 + qJDD(5) - t581;
t577 = -t619 * pkin(5) - t635 * pkin(9) + t637 * t617 + t579;
t689 = m(7) * t577 - t587 * mrSges(7,1) + mrSges(7,2) * t588 - t608 * t600 + t601 * t609;
t686 = -m(6) * t579 + t619 * mrSges(6,1) - mrSges(6,2) * t620 + t636 * t613 - t614 * t637 - t689;
t565 = m(5) * t581 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t632 + qJD(3) * t640 - t628 * t647 + t686;
t546 = t671 * t554 + t707 * t565;
t658 = (-mrSges(4,1) * t681 + mrSges(4,2) * t678) * qJD(2);
t544 = m(4) * t602 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t659 + qJD(3) * t665 - t658 * t702 + t546;
t695 = t707 * t554 - t565 * t671;
t545 = m(4) * t603 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t660 - qJD(3) * t664 + t658 * t701 + t695;
t696 = -t544 * t678 + t681 * t545;
t536 = m(3) * t625 - mrSges(3,1) * t684 - qJDD(2) * mrSges(3,2) + t696;
t539 = t681 * t544 + t678 * t545;
t538 = m(3) * t642 + t539;
t527 = t536 * t704 - t538 * t673 + t676 * t706;
t525 = m(2) * t661 + t527;
t531 = t682 * t536 - t551 * t679;
t530 = m(2) * t662 + t531;
t703 = t675 * t525 + t672 * t530;
t526 = t536 * t705 + t676 * t538 + t673 * t706;
t697 = -t525 * t672 + t675 * t530;
t589 = Ifges(7,5) * t609 + Ifges(7,6) * t608 + Ifges(7,3) * t645;
t591 = Ifges(7,1) * t609 + Ifges(7,4) * t608 + Ifges(7,5) * t645;
t562 = -mrSges(7,1) * t577 + mrSges(7,3) * t572 + Ifges(7,4) * t588 + Ifges(7,2) * t587 + Ifges(7,6) * t630 - t589 * t609 + t591 * t645;
t590 = Ifges(7,4) * t609 + Ifges(7,2) * t608 + Ifges(7,6) * t645;
t563 = mrSges(7,2) * t577 - mrSges(7,3) * t571 + Ifges(7,1) * t588 + Ifges(7,4) * t587 + Ifges(7,5) * t630 + t589 * t608 - t590 * t645;
t604 = Ifges(6,5) * t637 + Ifges(6,6) * t636 + Ifges(6,3) * t646;
t606 = Ifges(6,1) * t637 + Ifges(6,4) * t636 + Ifges(6,5) * t646;
t547 = -mrSges(6,1) * t579 + mrSges(6,3) * t576 + Ifges(6,4) * t620 + Ifges(6,2) * t619 + Ifges(6,6) * t631 - pkin(5) * t689 + pkin(9) * t693 + t680 * t562 + t677 * t563 - t637 * t604 + t646 * t606;
t605 = Ifges(6,4) * t637 + Ifges(6,2) * t636 + Ifges(6,6) * t646;
t548 = mrSges(6,2) * t579 - mrSges(6,3) * t575 + Ifges(6,1) * t620 + Ifges(6,4) * t619 + Ifges(6,5) * t631 - pkin(9) * t561 - t562 * t677 + t563 * t680 + t604 * t636 - t605 * t646;
t621 = Ifges(5,5) * t647 - Ifges(5,6) * t646 + Ifges(5,3) * qJD(3);
t622 = Ifges(5,4) * t647 - Ifges(5,2) * t646 + Ifges(5,6) * qJD(3);
t532 = mrSges(5,2) * t599 - mrSges(5,3) * t581 + Ifges(5,1) * t632 - Ifges(5,4) * t631 + Ifges(5,5) * qJDD(3) - qJ(5) * t555 - qJD(3) * t622 - t547 * t670 + t548 * t674 - t621 * t646;
t623 = Ifges(5,1) * t647 - Ifges(5,4) * t646 + Ifges(5,5) * qJD(3);
t540 = Ifges(5,4) * t632 + Ifges(5,6) * qJDD(3) - t647 * t621 + qJD(3) * t623 - mrSges(5,1) * t599 + mrSges(5,3) * t582 - Ifges(6,5) * t620 - Ifges(6,6) * t619 - t637 * t605 + t636 * t606 - mrSges(6,1) * t575 + mrSges(6,2) * t576 - Ifges(7,5) * t588 - Ifges(7,6) * t587 - Ifges(7,3) * t630 - t609 * t590 + t608 * t591 - mrSges(7,1) * t571 + mrSges(7,2) * t572 - pkin(5) * t561 - pkin(4) * t555 + (-Ifges(5,2) - Ifges(6,3)) * t631;
t649 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t678 + Ifges(4,6) * t681) * qJD(2);
t651 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t678 + Ifges(4,4) * t681) * qJD(2);
t521 = -mrSges(4,1) * t615 + mrSges(4,3) * t603 + Ifges(4,4) * t659 + Ifges(4,2) * t660 + Ifges(4,6) * qJDD(3) - pkin(3) * t687 + qJ(4) * t695 + qJD(3) * t651 + t671 * t532 + t540 * t707 - t649 * t702;
t650 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t678 + Ifges(4,2) * t681) * qJD(2);
t523 = mrSges(4,2) * t615 - mrSges(4,3) * t602 + Ifges(4,1) * t659 + Ifges(4,4) * t660 + Ifges(4,5) * qJDD(3) - qJ(4) * t546 - qJD(3) * t650 + t532 * t707 - t671 * t540 + t649 * t701;
t520 = mrSges(3,2) * t642 - mrSges(3,3) * t624 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t684 - pkin(8) * t539 - t521 * t678 + t523 * t681;
t522 = Ifges(3,6) * qJDD(2) - pkin(2) * t539 + mrSges(3,3) * t625 - mrSges(3,1) * t642 - pkin(3) * t546 - Ifges(4,5) * t659 - Ifges(4,6) * t660 - mrSges(4,1) * t602 + mrSges(4,2) * t603 - t670 * t548 - t674 * t547 - pkin(4) * t686 - qJ(5) * t694 - Ifges(5,5) * t632 + Ifges(5,6) * t631 - mrSges(5,1) * t581 + mrSges(5,2) * t582 - t647 * t622 - t646 * t623 + t684 * Ifges(3,5) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t650 * t678 + t651 * t681) * qJD(2);
t690 = pkin(7) * t531 + t520 * t679 + t522 * t682;
t519 = mrSges(3,1) * t624 - mrSges(3,2) * t625 + Ifges(3,3) * qJDD(2) + pkin(2) * t685 + pkin(8) * t696 + t681 * t521 + t678 * t523;
t518 = mrSges(2,2) * t669 - mrSges(2,3) * t661 + t682 * t520 - t679 * t522 + (-t526 * t673 - t527 * t676) * pkin(7);
t517 = -mrSges(2,1) * t669 + mrSges(2,3) * t662 - pkin(1) * t526 - t673 * t519 + t676 * t690;
t1 = [-m(1) * g(1) + t697; -m(1) * g(2) + t703; -m(1) * g(3) + m(2) * t669 + t526; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t703 - t672 * t517 + t675 * t518; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t697 + t675 * t517 + t672 * t518; -mrSges(1,1) * g(2) + mrSges(2,1) * t661 + mrSges(1,2) * g(1) - mrSges(2,2) * t662 + pkin(1) * t527 + t676 * t519 + t673 * t690;];
tauB  = t1;
