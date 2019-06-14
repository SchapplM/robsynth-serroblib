% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:08:01
% EndTime: 2019-05-06 01:08:12
% DurationCPUTime: 8.08s
% Computational Cost: add. (97659->320), mult. (190563->391), div. (0->0), fcn. (125735->10), ass. (0->126)
t681 = Ifges(6,1) + Ifges(7,1);
t675 = Ifges(6,4) - Ifges(7,5);
t680 = -Ifges(6,5) - Ifges(7,4);
t679 = Ifges(6,2) + Ifges(7,3);
t673 = Ifges(6,6) - Ifges(7,6);
t678 = -Ifges(6,3) - Ifges(7,2);
t644 = sin(qJ(4));
t645 = sin(qJ(3));
t647 = cos(qJ(4));
t648 = cos(qJ(3));
t615 = (t644 * t645 - t647 * t648) * qJD(1);
t677 = cos(qJ(5));
t676 = -mrSges(6,3) - mrSges(7,2);
t646 = sin(qJ(1));
t649 = cos(qJ(1));
t628 = t646 * g(1) - g(2) * t649;
t620 = qJDD(1) * pkin(1) + t628;
t629 = -g(1) * t649 - g(2) * t646;
t650 = qJD(1) ^ 2;
t622 = -pkin(1) * t650 + t629;
t641 = sin(pkin(10));
t642 = cos(pkin(10));
t602 = t641 * t620 + t642 * t622;
t598 = -pkin(2) * t650 + qJDD(1) * pkin(7) + t602;
t640 = -g(3) + qJDD(2);
t584 = -t645 * t598 + t648 * t640;
t665 = qJD(1) * qJD(3);
t662 = t648 * t665;
t623 = qJDD(1) * t645 + t662;
t563 = (-t623 + t662) * pkin(8) + (t645 * t648 * t650 + qJDD(3)) * pkin(3) + t584;
t585 = t648 * t598 + t645 * t640;
t624 = qJDD(1) * t648 - t645 * t665;
t667 = qJD(1) * t645;
t627 = qJD(3) * pkin(3) - pkin(8) * t667;
t639 = t648 ^ 2;
t564 = -pkin(3) * t639 * t650 + pkin(8) * t624 - qJD(3) * t627 + t585;
t555 = t644 * t563 + t647 * t564;
t616 = (t644 * t648 + t645 * t647) * qJD(1);
t586 = -qJD(4) * t616 - t623 * t644 + t624 * t647;
t599 = mrSges(5,1) * t615 + mrSges(5,2) * t616;
t638 = qJD(3) + qJD(4);
t606 = mrSges(5,1) * t638 - mrSges(5,3) * t616;
t637 = qJDD(3) + qJDD(4);
t600 = pkin(4) * t615 - pkin(9) * t616;
t636 = t638 ^ 2;
t551 = -pkin(4) * t636 + pkin(9) * t637 - t600 * t615 + t555;
t601 = t642 * t620 - t641 * t622;
t654 = -qJDD(1) * pkin(2) - t601;
t573 = -t624 * pkin(3) + t627 * t667 + (-pkin(8) * t639 - pkin(7)) * t650 + t654;
t587 = -qJD(4) * t615 + t623 * t647 + t624 * t644;
t553 = (t615 * t638 - t587) * pkin(9) + (t616 * t638 - t586) * pkin(4) + t573;
t643 = sin(qJ(5));
t548 = t677 * t551 + t643 * t553;
t604 = t677 * t616 + t643 * t638;
t558 = t604 * qJD(5) + t643 * t587 - t677 * t637;
t583 = qJDD(5) - t586;
t608 = qJD(5) + t615;
t590 = mrSges(6,1) * t608 - mrSges(6,3) * t604;
t603 = t643 * t616 - t677 * t638;
t576 = pkin(5) * t603 - qJ(6) * t604;
t607 = t608 ^ 2;
t544 = -pkin(5) * t607 + qJ(6) * t583 + 0.2e1 * qJD(6) * t608 - t576 * t603 + t548;
t591 = -mrSges(7,1) * t608 + mrSges(7,2) * t604;
t664 = m(7) * t544 + t583 * mrSges(7,3) + t608 * t591;
t577 = mrSges(7,1) * t603 - mrSges(7,3) * t604;
t668 = -mrSges(6,1) * t603 - mrSges(6,2) * t604 - t577;
t539 = m(6) * t548 - t583 * mrSges(6,2) + t676 * t558 - t608 * t590 + t668 * t603 + t664;
t547 = -t643 * t551 + t677 * t553;
t559 = -t603 * qJD(5) + t677 * t587 + t643 * t637;
t589 = -mrSges(6,2) * t608 - mrSges(6,3) * t603;
t545 = -t583 * pkin(5) - t607 * qJ(6) + t604 * t576 + qJDD(6) - t547;
t588 = -mrSges(7,2) * t603 + mrSges(7,3) * t608;
t656 = -m(7) * t545 + t583 * mrSges(7,1) + t608 * t588;
t541 = m(6) * t547 + t583 * mrSges(6,1) + t676 * t559 + t608 * t589 + t668 * t604 + t656;
t657 = t677 * t539 - t541 * t643;
t531 = m(5) * t555 - mrSges(5,2) * t637 + mrSges(5,3) * t586 - t599 * t615 - t606 * t638 + t657;
t554 = t647 * t563 - t644 * t564;
t605 = -mrSges(5,2) * t638 - mrSges(5,3) * t615;
t550 = -t637 * pkin(4) - t636 * pkin(9) + t616 * t600 - t554;
t546 = -0.2e1 * qJD(6) * t604 + (t603 * t608 - t559) * qJ(6) + (t604 * t608 + t558) * pkin(5) + t550;
t542 = m(7) * t546 + mrSges(7,1) * t558 - t559 * mrSges(7,3) + t588 * t603 - t604 * t591;
t652 = -m(6) * t550 - t558 * mrSges(6,1) - mrSges(6,2) * t559 - t603 * t589 - t590 * t604 - t542;
t536 = m(5) * t554 + mrSges(5,1) * t637 - mrSges(5,3) * t587 - t599 * t616 + t605 * t638 + t652;
t526 = t644 * t531 + t647 * t536;
t621 = (-mrSges(4,1) * t648 + mrSges(4,2) * t645) * qJD(1);
t666 = qJD(1) * t648;
t626 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t666;
t524 = m(4) * t584 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t623 + qJD(3) * t626 - t621 * t667 + t526;
t625 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t667;
t658 = t647 * t531 - t536 * t644;
t525 = m(4) * t585 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t624 - qJD(3) * t625 + t621 * t666 + t658;
t659 = -t524 * t645 + t648 * t525;
t517 = m(3) * t602 - mrSges(3,1) * t650 - qJDD(1) * mrSges(3,2) + t659;
t597 = -t650 * pkin(7) + t654;
t534 = t643 * t539 + t677 * t541;
t653 = m(5) * t573 - t586 * mrSges(5,1) + mrSges(5,2) * t587 + t615 * t605 + t606 * t616 + t534;
t651 = -m(4) * t597 + t624 * mrSges(4,1) - mrSges(4,2) * t623 - t625 * t667 + t626 * t666 - t653;
t528 = m(3) * t601 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t650 + t651;
t514 = t641 * t517 + t642 * t528;
t512 = m(2) * t628 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t650 + t514;
t660 = t642 * t517 - t528 * t641;
t513 = m(2) * t629 - mrSges(2,1) * t650 - qJDD(1) * mrSges(2,2) + t660;
t672 = t649 * t512 + t646 * t513;
t518 = t648 * t524 + t645 * t525;
t671 = t603 * t679 - t604 * t675 - t608 * t673;
t670 = t603 * t673 + t604 * t680 + t608 * t678;
t669 = -t675 * t603 + t604 * t681 - t680 * t608;
t663 = m(3) * t640 + t518;
t661 = -t512 * t646 + t649 * t513;
t614 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t645 + Ifges(4,4) * t648) * qJD(1);
t613 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t645 + Ifges(4,2) * t648) * qJD(1);
t612 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t645 + Ifges(4,6) * t648) * qJD(1);
t595 = Ifges(5,1) * t616 - Ifges(5,4) * t615 + Ifges(5,5) * t638;
t594 = Ifges(5,4) * t616 - Ifges(5,2) * t615 + Ifges(5,6) * t638;
t593 = Ifges(5,5) * t616 - Ifges(5,6) * t615 + Ifges(5,3) * t638;
t533 = mrSges(6,2) * t550 + mrSges(7,2) * t545 - mrSges(6,3) * t547 - mrSges(7,3) * t546 - qJ(6) * t542 - t675 * t558 + t559 * t681 - t583 * t680 + t670 * t603 + t671 * t608;
t532 = -mrSges(6,1) * t550 - mrSges(7,1) * t546 + mrSges(7,2) * t544 + mrSges(6,3) * t548 - pkin(5) * t542 - t558 * t679 + t675 * t559 + t673 * t583 + t670 * t604 + t669 * t608;
t520 = Ifges(5,4) * t587 + Ifges(5,2) * t586 + Ifges(5,6) * t637 - t616 * t593 + t638 * t595 - mrSges(5,1) * t573 + mrSges(5,3) * t555 - mrSges(6,1) * t547 + mrSges(6,2) * t548 + mrSges(7,1) * t545 - mrSges(7,3) * t544 - pkin(5) * t656 - qJ(6) * t664 - pkin(4) * t534 + (pkin(5) * t577 + t671) * t604 + (qJ(6) * t577 - t669) * t603 + t678 * t583 + (mrSges(7,2) * pkin(5) + t680) * t559 + (mrSges(7,2) * qJ(6) + t673) * t558;
t519 = mrSges(5,2) * t573 - mrSges(5,3) * t554 + Ifges(5,1) * t587 + Ifges(5,4) * t586 + Ifges(5,5) * t637 - pkin(9) * t534 - t643 * t532 + t677 * t533 - t615 * t593 - t638 * t594;
t508 = mrSges(4,2) * t597 - mrSges(4,3) * t584 + Ifges(4,1) * t623 + Ifges(4,4) * t624 + Ifges(4,5) * qJDD(3) - pkin(8) * t526 - qJD(3) * t613 + t519 * t647 - t520 * t644 + t612 * t666;
t507 = -mrSges(4,1) * t597 + mrSges(4,3) * t585 + Ifges(4,4) * t623 + Ifges(4,2) * t624 + Ifges(4,6) * qJDD(3) - pkin(3) * t653 + pkin(8) * t658 + qJD(3) * t614 + t644 * t519 + t647 * t520 - t612 * t667;
t506 = -pkin(2) * t518 - mrSges(3,1) * t640 + mrSges(3,3) * t602 - pkin(3) * t526 - Ifges(4,5) * t623 - Ifges(4,6) * t624 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t584 + mrSges(4,2) * t585 - t677 * t532 - pkin(4) * t652 - pkin(9) * t657 - Ifges(5,5) * t587 - Ifges(5,6) * t586 - Ifges(5,3) * t637 - mrSges(5,1) * t554 + mrSges(5,2) * t555 - t643 * t533 - t616 * t594 - t615 * t595 + t650 * Ifges(3,5) + Ifges(3,6) * qJDD(1) + (-t613 * t645 + t614 * t648) * qJD(1);
t505 = mrSges(3,2) * t640 - mrSges(3,3) * t601 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t650 - pkin(7) * t518 - t507 * t645 + t508 * t648;
t504 = -mrSges(2,2) * g(3) - mrSges(2,3) * t628 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t650 - qJ(2) * t514 + t505 * t642 - t506 * t641;
t503 = mrSges(2,1) * g(3) + mrSges(2,3) * t629 + t650 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t663 + qJ(2) * t660 + t641 * t505 + t642 * t506;
t1 = [-m(1) * g(1) + t661; -m(1) * g(2) + t672; (-m(1) - m(2)) * g(3) + t663; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t672 - t646 * t503 + t649 * t504; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t661 + t649 * t503 + t646 * t504; pkin(1) * t514 + mrSges(2,1) * t628 - mrSges(2,2) * t629 + t645 * t508 + t648 * t507 + pkin(2) * t651 + pkin(7) * t659 + mrSges(3,1) * t601 - mrSges(3,2) * t602 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
