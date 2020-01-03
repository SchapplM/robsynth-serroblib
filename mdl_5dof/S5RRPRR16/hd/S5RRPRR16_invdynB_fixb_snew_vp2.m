% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR16_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:45:09
% EndTime: 2019-12-31 20:45:19
% DurationCPUTime: 5.79s
% Computational Cost: add. (61828->309), mult. (139973->388), div. (0->0), fcn. (97738->10), ass. (0->133)
t681 = -2 * qJD(3);
t680 = Ifges(3,1) + Ifges(4,2);
t673 = Ifges(3,4) + Ifges(4,6);
t672 = Ifges(3,5) - Ifges(4,4);
t679 = Ifges(3,2) + Ifges(4,3);
t671 = Ifges(3,6) - Ifges(4,5);
t678 = Ifges(3,3) + Ifges(4,1);
t632 = cos(pkin(5));
t626 = qJD(1) * t632 + qJD(2);
t635 = sin(qJ(2));
t631 = sin(pkin(5));
t660 = qJD(1) * t631;
t654 = t635 * t660;
t677 = (pkin(2) * t626 + t681) * t654;
t636 = sin(qJ(1));
t640 = cos(qJ(1));
t621 = t636 * g(1) - g(2) * t640;
t641 = qJD(1) ^ 2;
t606 = pkin(7) * t631 * t641 + qJDD(1) * pkin(1) + t621;
t622 = -g(1) * t640 - g(2) * t636;
t657 = qJDD(1) * t631;
t607 = -pkin(1) * t641 + pkin(7) * t657 + t622;
t639 = cos(qJ(2));
t667 = t632 * t635;
t669 = t631 * t635;
t573 = -g(3) * t669 + t606 * t667 + t639 * t607;
t608 = (-pkin(2) * t639 - qJ(3) * t635) * t660;
t624 = t626 ^ 2;
t625 = qJDD(1) * t632 + qJDD(2);
t659 = qJD(1) * t639;
t653 = t631 * t659;
t554 = pkin(2) * t624 - t625 * qJ(3) - t608 * t653 + t626 * t681 - t573;
t676 = -pkin(2) - pkin(8);
t675 = g(3) * t632;
t674 = mrSges(3,1) - mrSges(4,2);
t670 = t631 ^ 2 * t641;
t668 = t631 * t639;
t666 = t632 * t639;
t586 = -t606 * t631 - t675;
t602 = mrSges(3,1) * t626 - mrSges(3,3) * t654;
t603 = -mrSges(3,2) * t626 + mrSges(3,3) * t653;
t605 = mrSges(4,1) * t654 + mrSges(4,2) * t626;
t612 = (qJD(2) * t659 + qJDD(1) * t635) * t631;
t613 = -qJD(2) * t654 + t639 * t657;
t555 = -pkin(2) * t613 + (-t626 * t653 - t612) * qJ(3) + t586 + t677;
t604 = -mrSges(4,1) * t653 - mrSges(4,3) * t626;
t611 = pkin(3) * t654 - pkin(8) * t626;
t655 = t639 ^ 2 * t670;
t547 = -pkin(3) * t655 - t675 - qJ(3) * t612 + t676 * t613 + (-t606 + (-qJ(3) * t626 * t639 - t611 * t635) * qJD(1)) * t631 + t677;
t661 = g(3) * t668 + t635 * t607;
t648 = -qJ(3) * t624 + t608 * t654 + qJDD(3) + t661;
t549 = pkin(3) * t612 + t676 * t625 + (-pkin(3) * t626 * t660 - pkin(8) * t635 * t670 - t606 * t632) * t639 + t648;
t634 = sin(qJ(4));
t638 = cos(qJ(4));
t543 = t638 * t547 + t634 * t549;
t596 = t626 * t638 - t634 * t653;
t570 = -qJD(4) * t596 - t613 * t638 - t625 * t634;
t595 = -t626 * t634 - t638 * t653;
t574 = -mrSges(5,1) * t595 + mrSges(5,2) * t596;
t617 = qJD(4) + t654;
t579 = mrSges(5,1) * t617 - mrSges(5,3) * t596;
t601 = qJDD(4) + t612;
t575 = -pkin(4) * t595 - pkin(9) * t596;
t615 = t617 ^ 2;
t540 = -pkin(4) * t615 + pkin(9) * t601 + t575 * t595 + t543;
t546 = pkin(3) * t613 - pkin(8) * t655 + t626 * t611 - t554;
t571 = qJD(4) * t595 - t613 * t634 + t625 * t638;
t541 = (-t595 * t617 - t571) * pkin(9) + (t596 * t617 - t570) * pkin(4) + t546;
t633 = sin(qJ(5));
t637 = cos(qJ(5));
t537 = -t540 * t633 + t541 * t637;
t576 = -t596 * t633 + t617 * t637;
t552 = qJD(5) * t576 + t571 * t637 + t601 * t633;
t577 = t596 * t637 + t617 * t633;
t561 = -mrSges(6,1) * t576 + mrSges(6,2) * t577;
t594 = qJD(5) - t595;
t562 = -mrSges(6,2) * t594 + mrSges(6,3) * t576;
t568 = qJDD(5) - t570;
t535 = m(6) * t537 + mrSges(6,1) * t568 - mrSges(6,3) * t552 - t561 * t577 + t562 * t594;
t538 = t540 * t637 + t541 * t633;
t551 = -qJD(5) * t577 - t571 * t633 + t601 * t637;
t563 = mrSges(6,1) * t594 - mrSges(6,3) * t577;
t536 = m(6) * t538 - mrSges(6,2) * t568 + mrSges(6,3) * t551 + t561 * t576 - t563 * t594;
t650 = -t535 * t633 + t637 * t536;
t527 = m(5) * t543 - mrSges(5,2) * t601 + mrSges(5,3) * t570 + t574 * t595 - t579 * t617 + t650;
t542 = -t547 * t634 + t549 * t638;
t578 = -mrSges(5,2) * t617 + mrSges(5,3) * t595;
t539 = -pkin(4) * t601 - pkin(9) * t615 + t575 * t596 - t542;
t645 = -m(6) * t539 + t551 * mrSges(6,1) - mrSges(6,2) * t552 + t576 * t562 - t563 * t577;
t531 = m(5) * t542 + mrSges(5,1) * t601 - mrSges(5,3) * t571 - t574 * t596 + t578 * t617 + t645;
t651 = t638 * t527 - t531 * t634;
t649 = m(4) * t555 - t612 * mrSges(4,3) + t604 * t653 + t651;
t518 = m(3) * t586 + mrSges(3,2) * t612 - t674 * t613 + (-t603 * t639 + (t602 - t605) * t635) * t660 + t649;
t656 = t606 * t666;
t572 = t656 - t661;
t609 = (mrSges(4,2) * t639 - mrSges(4,3) * t635) * t660;
t610 = (-mrSges(3,1) * t639 + mrSges(3,2) * t635) * t660;
t521 = t527 * t634 + t531 * t638;
t560 = -pkin(2) * t625 + t648 - t656;
t646 = -m(4) * t560 - t612 * mrSges(4,1) - t521;
t519 = m(3) * t572 - mrSges(3,3) * t612 + (t603 - t604) * t626 + t674 * t625 + (-t609 - t610) * t654 + t646;
t528 = t637 * t535 + t633 * t536;
t643 = -m(5) * t546 + mrSges(5,1) * t570 - t571 * mrSges(5,2) + t595 * t578 - t596 * t579 - t528;
t642 = -m(4) * t554 + t625 * mrSges(4,3) + t626 * t605 + t609 * t653 - t643;
t525 = t642 + (mrSges(3,3) + mrSges(4,1)) * t613 + m(3) * t573 - mrSges(3,2) * t625 - t602 * t626 + t610 * t653;
t508 = -t518 * t631 + t519 * t666 + t525 * t667;
t506 = m(2) * t621 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t641 + t508;
t512 = -t519 * t635 + t639 * t525;
t511 = m(2) * t622 - mrSges(2,1) * t641 - qJDD(1) * mrSges(2,2) + t512;
t665 = t640 * t506 + t636 * t511;
t664 = (t672 * t635 + t671 * t639) * t660 + t678 * t626;
t663 = (-t673 * t635 - t679 * t639) * t660 - t671 * t626;
t662 = (t680 * t635 + t673 * t639) * t660 + t672 * t626;
t507 = t632 * t518 + t519 * t668 + t525 * t669;
t652 = -t506 * t636 + t640 * t511;
t556 = Ifges(6,5) * t577 + Ifges(6,6) * t576 + Ifges(6,3) * t594;
t558 = Ifges(6,1) * t577 + Ifges(6,4) * t576 + Ifges(6,5) * t594;
t529 = -mrSges(6,1) * t539 + mrSges(6,3) * t538 + Ifges(6,4) * t552 + Ifges(6,2) * t551 + Ifges(6,6) * t568 - t556 * t577 + t558 * t594;
t557 = Ifges(6,4) * t577 + Ifges(6,2) * t576 + Ifges(6,6) * t594;
t530 = mrSges(6,2) * t539 - mrSges(6,3) * t537 + Ifges(6,1) * t552 + Ifges(6,4) * t551 + Ifges(6,5) * t568 + t556 * t576 - t557 * t594;
t564 = Ifges(5,5) * t596 + Ifges(5,6) * t595 + Ifges(5,3) * t617;
t565 = Ifges(5,4) * t596 + Ifges(5,2) * t595 + Ifges(5,6) * t617;
t513 = mrSges(5,2) * t546 - mrSges(5,3) * t542 + Ifges(5,1) * t571 + Ifges(5,4) * t570 + Ifges(5,5) * t601 - pkin(9) * t528 - t529 * t633 + t530 * t637 + t564 * t595 - t565 * t617;
t566 = Ifges(5,1) * t596 + Ifges(5,4) * t595 + Ifges(5,5) * t617;
t514 = -mrSges(5,1) * t546 - mrSges(6,1) * t537 + mrSges(6,2) * t538 + mrSges(5,3) * t543 + Ifges(5,4) * t571 - Ifges(6,5) * t552 + Ifges(5,2) * t570 + Ifges(5,6) * t601 - Ifges(6,6) * t551 - Ifges(6,3) * t568 - pkin(4) * t528 - t557 * t577 + t558 * t576 - t564 * t596 + t566 * t617;
t520 = mrSges(4,2) * t613 - t605 * t654 + t649;
t503 = -mrSges(3,1) * t586 - mrSges(4,1) * t554 + mrSges(4,2) * t555 + mrSges(3,3) * t573 - pkin(2) * t520 - pkin(3) * t643 - pkin(8) * t651 - t634 * t513 - t638 * t514 + t673 * t612 + t679 * t613 + t671 * t625 + t662 * t626 - t664 * t654;
t504 = pkin(4) * t645 + mrSges(3,2) * t586 + pkin(9) * t650 + t633 * t530 + t637 * t529 + pkin(3) * t521 - qJ(3) * t520 - mrSges(4,3) * t555 + mrSges(4,1) * t560 + Ifges(5,6) * t570 + Ifges(5,5) * t571 - mrSges(3,3) * t572 + Ifges(5,3) * t601 + mrSges(5,1) * t542 - mrSges(5,2) * t543 - t595 * t566 + t596 * t565 + t663 * t626 + t672 * t625 + t673 * t613 + t680 * t612 + t664 * t653;
t647 = pkin(7) * t512 + t503 * t639 + t504 * t635;
t502 = mrSges(3,1) * t572 - mrSges(3,2) * t573 + mrSges(4,2) * t560 - mrSges(4,3) * t554 + t638 * t513 - t634 * t514 - pkin(8) * t521 + pkin(2) * (-t604 * t626 + t646) + qJ(3) * t642 + (-pkin(2) * mrSges(4,2) + t678) * t625 + (qJ(3) * mrSges(4,1) + t671) * t613 + t672 * t612 + (-t662 * t639 + (-pkin(2) * t609 - t663) * t635) * t660;
t501 = -mrSges(2,2) * g(3) - mrSges(2,3) * t621 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t641 - t503 * t635 + t504 * t639 + (-t507 * t631 - t508 * t632) * pkin(7);
t500 = mrSges(2,1) * g(3) + mrSges(2,3) * t622 + Ifges(2,5) * t641 + Ifges(2,6) * qJDD(1) - pkin(1) * t507 - t502 * t631 + t647 * t632;
t1 = [-m(1) * g(1) + t652; -m(1) * g(2) + t665; (-m(1) - m(2)) * g(3) + t507; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t665 - t636 * t500 + t640 * t501; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t652 + t640 * t500 + t636 * t501; -mrSges(1,1) * g(2) + mrSges(2,1) * t621 + mrSges(1,2) * g(1) - mrSges(2,2) * t622 + Ifges(2,3) * qJDD(1) + pkin(1) * t508 + t502 * t632 + t647 * t631;];
tauB = t1;
