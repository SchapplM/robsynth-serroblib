% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:13:42
% EndTime: 2019-05-05 07:14:03
% DurationCPUTime: 20.10s
% Computational Cost: add. (355376->341), mult. (715659->442), div. (0->0), fcn. (522227->14), ass. (0->141)
t673 = sin(pkin(11));
t676 = cos(pkin(11));
t660 = g(1) * t673 - g(2) * t676;
t661 = -g(1) * t676 - g(2) * t673;
t671 = -g(3) + qJDD(1);
t674 = sin(pkin(6));
t677 = cos(pkin(6));
t681 = sin(qJ(2));
t684 = cos(qJ(2));
t629 = -t681 * t661 + (t660 * t677 + t671 * t674) * t684;
t706 = cos(qJ(4));
t685 = qJD(2) ^ 2;
t689 = -qJDD(2) * pkin(2) - t629;
t624 = -t685 * pkin(8) + t689;
t680 = sin(qJ(3));
t683 = cos(qJ(3));
t699 = qJD(2) * qJD(3);
t698 = t683 * t699;
t658 = qJDD(2) * t680 + t698;
t659 = qJDD(2) * t683 - t680 * t699;
t701 = qJD(2) * t680;
t662 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t701;
t700 = qJD(2) * t683;
t663 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t700;
t703 = t677 * t681;
t704 = t674 * t681;
t630 = t660 * t703 + t684 * t661 + t671 * t704;
t625 = -pkin(2) * t685 + qJDD(2) * pkin(8) + t630;
t643 = -t660 * t674 + t671 * t677;
t607 = -t680 * t625 + t683 * t643;
t598 = (-t658 + t698) * pkin(9) + (t680 * t683 * t685 + qJDD(3)) * pkin(3) + t607;
t608 = t683 * t625 + t680 * t643;
t665 = qJD(3) * pkin(3) - pkin(9) * t701;
t670 = t683 ^ 2;
t599 = -pkin(3) * t670 * t685 + pkin(9) * t659 - qJD(3) * t665 + t608;
t679 = sin(qJ(4));
t586 = t679 * t598 + t599 * t706;
t649 = t679 * t701 - t706 * t700;
t650 = (t679 * t683 + t680 * t706) * qJD(2);
t632 = pkin(4) * t649 - qJ(5) * t650;
t669 = qJD(3) + qJD(4);
t667 = t669 ^ 2;
t668 = qJDD(3) + qJDD(4);
t581 = -pkin(4) * t667 + qJ(5) * t668 - t632 * t649 + t586;
t602 = -t659 * pkin(3) + t665 * t701 + (-pkin(9) * t670 - pkin(8)) * t685 + t689;
t619 = qJD(4) * t650 + t658 * t679 - t659 * t706;
t620 = -t649 * qJD(4) + t658 * t706 + t679 * t659;
t584 = (t649 * t669 - t620) * qJ(5) + (t650 * t669 + t619) * pkin(4) + t602;
t672 = sin(pkin(12));
t675 = cos(pkin(12));
t638 = t650 * t675 + t669 * t672;
t576 = -0.2e1 * qJD(5) * t638 - t672 * t581 + t675 * t584;
t613 = t620 * t675 + t668 * t672;
t637 = -t650 * t672 + t669 * t675;
t574 = (t637 * t649 - t613) * pkin(10) + (t637 * t638 + t619) * pkin(5) + t576;
t577 = 0.2e1 * qJD(5) * t637 + t675 * t581 + t672 * t584;
t612 = -t620 * t672 + t668 * t675;
t623 = pkin(5) * t649 - pkin(10) * t638;
t636 = t637 ^ 2;
t575 = -pkin(5) * t636 + pkin(10) * t612 - t623 * t649 + t577;
t678 = sin(qJ(6));
t682 = cos(qJ(6));
t572 = t574 * t682 - t575 * t678;
t610 = t637 * t682 - t638 * t678;
t589 = qJD(6) * t610 + t612 * t678 + t613 * t682;
t611 = t637 * t678 + t638 * t682;
t594 = -mrSges(7,1) * t610 + mrSges(7,2) * t611;
t644 = qJD(6) + t649;
t600 = -mrSges(7,2) * t644 + mrSges(7,3) * t610;
t617 = qJDD(6) + t619;
t570 = m(7) * t572 + mrSges(7,1) * t617 - mrSges(7,3) * t589 - t594 * t611 + t600 * t644;
t573 = t574 * t678 + t575 * t682;
t588 = -qJD(6) * t611 + t612 * t682 - t613 * t678;
t601 = mrSges(7,1) * t644 - mrSges(7,3) * t611;
t571 = m(7) * t573 - mrSges(7,2) * t617 + mrSges(7,3) * t588 + t594 * t610 - t601 * t644;
t562 = t682 * t570 + t678 * t571;
t614 = -mrSges(6,1) * t637 + mrSges(6,2) * t638;
t621 = -mrSges(6,2) * t649 + mrSges(6,3) * t637;
t560 = m(6) * t576 + mrSges(6,1) * t619 - mrSges(6,3) * t613 - t614 * t638 + t621 * t649 + t562;
t622 = mrSges(6,1) * t649 - mrSges(6,3) * t638;
t693 = -t570 * t678 + t682 * t571;
t561 = m(6) * t577 - mrSges(6,2) * t619 + mrSges(6,3) * t612 + t614 * t637 - t622 * t649 + t693;
t556 = t675 * t560 + t672 * t561;
t641 = -mrSges(5,2) * t669 - mrSges(5,3) * t649;
t642 = mrSges(5,1) * t669 - mrSges(5,3) * t650;
t688 = m(5) * t602 + t619 * mrSges(5,1) + mrSges(5,2) * t620 + t649 * t641 + t642 * t650 + t556;
t686 = -m(4) * t624 + t659 * mrSges(4,1) - mrSges(4,2) * t658 - t662 * t701 + t663 * t700 - t688;
t552 = m(3) * t629 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t685 + t686;
t705 = t552 * t684;
t633 = mrSges(5,1) * t649 + mrSges(5,2) * t650;
t694 = -t560 * t672 + t675 * t561;
t555 = m(5) * t586 - mrSges(5,2) * t668 - mrSges(5,3) * t619 - t633 * t649 - t642 * t669 + t694;
t585 = t598 * t706 - t679 * t599;
t580 = -t668 * pkin(4) - t667 * qJ(5) + t650 * t632 + qJDD(5) - t585;
t578 = -t612 * pkin(5) - t636 * pkin(10) + t638 * t623 + t580;
t690 = m(7) * t578 - t588 * mrSges(7,1) + mrSges(7,2) * t589 - t610 * t600 + t601 * t611;
t687 = -m(6) * t580 + t612 * mrSges(6,1) - mrSges(6,2) * t613 + t637 * t621 - t622 * t638 - t690;
t566 = m(5) * t585 + mrSges(5,1) * t668 - mrSges(5,3) * t620 - t633 * t650 + t641 * t669 + t687;
t547 = t679 * t555 + t566 * t706;
t657 = (-mrSges(4,1) * t683 + mrSges(4,2) * t680) * qJD(2);
t545 = m(4) * t607 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t658 + qJD(3) * t663 - t657 * t701 + t547;
t695 = t555 * t706 - t566 * t679;
t546 = m(4) * t608 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t659 - qJD(3) * t662 + t657 * t700 + t695;
t696 = -t545 * t680 + t683 * t546;
t537 = m(3) * t630 - mrSges(3,1) * t685 - qJDD(2) * mrSges(3,2) + t696;
t540 = t683 * t545 + t680 * t546;
t539 = m(3) * t643 + t540;
t528 = t537 * t703 - t539 * t674 + t677 * t705;
t526 = m(2) * t660 + t528;
t532 = t684 * t537 - t552 * t681;
t531 = m(2) * t661 + t532;
t702 = t676 * t526 + t673 * t531;
t527 = t537 * t704 + t677 * t539 + t674 * t705;
t697 = -t526 * t673 + t676 * t531;
t590 = Ifges(7,5) * t611 + Ifges(7,6) * t610 + Ifges(7,3) * t644;
t592 = Ifges(7,1) * t611 + Ifges(7,4) * t610 + Ifges(7,5) * t644;
t563 = -mrSges(7,1) * t578 + mrSges(7,3) * t573 + Ifges(7,4) * t589 + Ifges(7,2) * t588 + Ifges(7,6) * t617 - t590 * t611 + t592 * t644;
t591 = Ifges(7,4) * t611 + Ifges(7,2) * t610 + Ifges(7,6) * t644;
t564 = mrSges(7,2) * t578 - mrSges(7,3) * t572 + Ifges(7,1) * t589 + Ifges(7,4) * t588 + Ifges(7,5) * t617 + t590 * t610 - t591 * t644;
t603 = Ifges(6,5) * t638 + Ifges(6,6) * t637 + Ifges(6,3) * t649;
t605 = Ifges(6,1) * t638 + Ifges(6,4) * t637 + Ifges(6,5) * t649;
t548 = -mrSges(6,1) * t580 + mrSges(6,3) * t577 + Ifges(6,4) * t613 + Ifges(6,2) * t612 + Ifges(6,6) * t619 - pkin(5) * t690 + pkin(10) * t693 + t682 * t563 + t678 * t564 - t638 * t603 + t649 * t605;
t604 = Ifges(6,4) * t638 + Ifges(6,2) * t637 + Ifges(6,6) * t649;
t549 = mrSges(6,2) * t580 - mrSges(6,3) * t576 + Ifges(6,1) * t613 + Ifges(6,4) * t612 + Ifges(6,5) * t619 - pkin(10) * t562 - t563 * t678 + t564 * t682 + t603 * t637 - t604 * t649;
t626 = Ifges(5,5) * t650 - Ifges(5,6) * t649 + Ifges(5,3) * t669;
t627 = Ifges(5,4) * t650 - Ifges(5,2) * t649 + Ifges(5,6) * t669;
t533 = mrSges(5,2) * t602 - mrSges(5,3) * t585 + Ifges(5,1) * t620 - Ifges(5,4) * t619 + Ifges(5,5) * t668 - qJ(5) * t556 - t548 * t672 + t549 * t675 - t626 * t649 - t627 * t669;
t628 = Ifges(5,1) * t650 - Ifges(5,4) * t649 + Ifges(5,5) * t669;
t541 = Ifges(5,4) * t620 + Ifges(5,6) * t668 - t650 * t626 + t669 * t628 - mrSges(5,1) * t602 + mrSges(5,3) * t586 - Ifges(6,5) * t613 - Ifges(6,6) * t612 - t638 * t604 + t637 * t605 - mrSges(6,1) * t576 + mrSges(6,2) * t577 - Ifges(7,5) * t589 - Ifges(7,6) * t588 - Ifges(7,3) * t617 - t611 * t591 + t610 * t592 - mrSges(7,1) * t572 + mrSges(7,2) * t573 - pkin(5) * t562 - pkin(4) * t556 + (-Ifges(5,2) - Ifges(6,3)) * t619;
t646 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t680 + Ifges(4,6) * t683) * qJD(2);
t648 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t680 + Ifges(4,4) * t683) * qJD(2);
t522 = -mrSges(4,1) * t624 + mrSges(4,3) * t608 + Ifges(4,4) * t658 + Ifges(4,2) * t659 + Ifges(4,6) * qJDD(3) - pkin(3) * t688 + pkin(9) * t695 + qJD(3) * t648 + t679 * t533 + t541 * t706 - t646 * t701;
t647 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t680 + Ifges(4,2) * t683) * qJD(2);
t524 = mrSges(4,2) * t624 - mrSges(4,3) * t607 + Ifges(4,1) * t658 + Ifges(4,4) * t659 + Ifges(4,5) * qJDD(3) - pkin(9) * t547 - qJD(3) * t647 + t533 * t706 - t679 * t541 + t646 * t700;
t521 = mrSges(3,2) * t643 - mrSges(3,3) * t629 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t685 - pkin(8) * t540 - t522 * t680 + t524 * t683;
t523 = t685 * Ifges(3,5) - t650 * t627 - t649 * t628 + Ifges(3,6) * qJDD(2) - pkin(2) * t540 + mrSges(3,3) * t630 - mrSges(3,1) * t643 - pkin(3) * t547 - Ifges(4,5) * t658 - Ifges(4,6) * t659 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t607 + mrSges(4,2) * t608 - qJ(5) * t694 - mrSges(5,1) * t585 + mrSges(5,2) * t586 - t672 * t549 - t675 * t548 - pkin(4) * t687 - Ifges(5,5) * t620 + Ifges(5,6) * t619 - Ifges(5,3) * t668 + (-t647 * t680 + t648 * t683) * qJD(2);
t691 = pkin(7) * t532 + t521 * t681 + t523 * t684;
t520 = mrSges(3,1) * t629 - mrSges(3,2) * t630 + Ifges(3,3) * qJDD(2) + pkin(2) * t686 + pkin(8) * t696 + t683 * t522 + t680 * t524;
t519 = mrSges(2,2) * t671 - mrSges(2,3) * t660 + t684 * t521 - t681 * t523 + (-t527 * t674 - t528 * t677) * pkin(7);
t518 = -mrSges(2,1) * t671 + mrSges(2,3) * t661 - pkin(1) * t527 - t674 * t520 + t677 * t691;
t1 = [-m(1) * g(1) + t697; -m(1) * g(2) + t702; -m(1) * g(3) + m(2) * t671 + t527; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t702 - t673 * t518 + t676 * t519; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t697 + t676 * t518 + t673 * t519; -mrSges(1,1) * g(2) + mrSges(2,1) * t660 + mrSges(1,2) * g(1) - mrSges(2,2) * t661 + pkin(1) * t528 + t677 * t520 + t674 * t691;];
tauB  = t1;
