% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:01:10
% EndTime: 2019-05-05 23:01:18
% DurationCPUTime: 7.85s
% Computational Cost: add. (117753->339), mult. (253785->418), div. (0->0), fcn. (176710->10), ass. (0->132)
t646 = sin(qJ(1));
t650 = cos(qJ(1));
t629 = -t650 * g(1) - t646 * g(2);
t659 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t629;
t677 = 2 * qJD(5);
t676 = -pkin(1) - pkin(7);
t675 = mrSges(2,1) - mrSges(3,2);
t674 = Ifges(2,5) - Ifges(3,4);
t673 = (-Ifges(2,6) + Ifges(3,5));
t628 = t646 * g(1) - g(2) * t650;
t651 = qJD(1) ^ 2;
t658 = -t651 * qJ(2) + qJDD(2) - t628;
t609 = qJDD(1) * t676 + t658;
t645 = sin(qJ(3));
t649 = cos(qJ(3));
t600 = g(3) * t645 + t609 * t649;
t669 = qJD(1) * qJD(3);
t667 = t645 * t669;
t624 = qJDD(1) * t649 - t667;
t579 = (-t624 - t667) * pkin(8) + (-t645 * t649 * t651 + qJDD(3)) * pkin(3) + t600;
t601 = -g(3) * t649 + t609 * t645;
t623 = -qJDD(1) * t645 - t649 * t669;
t670 = qJD(1) * t649;
t627 = qJD(3) * pkin(3) - pkin(8) * t670;
t638 = t645 ^ 2;
t580 = -pkin(3) * t638 * t651 + pkin(8) * t623 - qJD(3) * t627 + t601;
t644 = sin(qJ(4));
t648 = cos(qJ(4));
t561 = t579 * t648 - t644 * t580;
t616 = (-t644 * t649 - t645 * t648) * qJD(1);
t587 = qJD(4) * t616 + t623 * t644 + t624 * t648;
t617 = (-t644 * t645 + t648 * t649) * qJD(1);
t634 = qJDD(3) + qJDD(4);
t635 = qJD(3) + qJD(4);
t549 = (t616 * t635 - t587) * qJ(5) + (t616 * t617 + t634) * pkin(4) + t561;
t562 = t579 * t644 + t580 * t648;
t586 = -qJD(4) * t617 + t623 * t648 - t624 * t644;
t606 = pkin(4) * t635 - qJ(5) * t617;
t612 = t616 ^ 2;
t551 = -pkin(4) * t612 + qJ(5) * t586 - t606 * t635 + t562;
t641 = sin(pkin(10));
t642 = cos(pkin(10));
t597 = t616 * t642 - t617 * t641;
t546 = t549 * t641 + t551 * t642 + t597 * t677;
t565 = t586 * t642 - t587 * t641;
t598 = t616 * t641 + t617 * t642;
t574 = -mrSges(6,1) * t597 + mrSges(6,2) * t598;
t589 = mrSges(6,1) * t635 - mrSges(6,3) * t598;
t575 = -pkin(5) * t597 - pkin(9) * t598;
t633 = t635 ^ 2;
t544 = -pkin(5) * t633 + pkin(9) * t634 + t575 * t597 + t546;
t582 = -t623 * pkin(3) + t627 * t670 + (-pkin(8) * t638 + t676) * t651 + t659;
t556 = -t586 * pkin(4) - t612 * qJ(5) + t606 * t617 + qJDD(5) + t582;
t566 = t586 * t641 + t587 * t642;
t547 = (-t597 * t635 - t566) * pkin(9) + (t598 * t635 - t565) * pkin(5) + t556;
t643 = sin(qJ(6));
t647 = cos(qJ(6));
t541 = -t544 * t643 + t547 * t647;
t584 = -t598 * t643 + t635 * t647;
t554 = qJD(6) * t584 + t566 * t647 + t634 * t643;
t564 = qJDD(6) - t565;
t585 = t598 * t647 + t635 * t643;
t567 = -mrSges(7,1) * t584 + mrSges(7,2) * t585;
t591 = qJD(6) - t597;
t568 = -mrSges(7,2) * t591 + mrSges(7,3) * t584;
t539 = m(7) * t541 + mrSges(7,1) * t564 - mrSges(7,3) * t554 - t567 * t585 + t568 * t591;
t542 = t544 * t647 + t547 * t643;
t553 = -qJD(6) * t585 - t566 * t643 + t634 * t647;
t569 = mrSges(7,1) * t591 - mrSges(7,3) * t585;
t540 = m(7) * t542 - mrSges(7,2) * t564 + mrSges(7,3) * t553 + t567 * t584 - t569 * t591;
t662 = -t539 * t643 + t540 * t647;
t530 = m(6) * t546 - mrSges(6,2) * t634 + mrSges(6,3) * t565 + t574 * t597 - t589 * t635 + t662;
t661 = -t642 * t549 + t641 * t551;
t545 = -0.2e1 * qJD(5) * t598 - t661;
t588 = -mrSges(6,2) * t635 + mrSges(6,3) * t597;
t543 = -t634 * pkin(5) - t633 * pkin(9) + (t677 + t575) * t598 + t661;
t656 = -m(7) * t543 + mrSges(7,1) * t553 - mrSges(7,2) * t554 + t568 * t584 - t569 * t585;
t535 = m(6) * t545 + mrSges(6,1) * t634 - mrSges(6,3) * t566 - t574 * t598 + t588 * t635 + t656;
t524 = t530 * t641 + t535 * t642;
t599 = -mrSges(5,1) * t616 + mrSges(5,2) * t617;
t605 = -mrSges(5,2) * t635 + mrSges(5,3) * t616;
t522 = m(5) * t561 + mrSges(5,1) * t634 - mrSges(5,3) * t587 - t599 * t617 + t605 * t635 + t524;
t607 = mrSges(5,1) * t635 - mrSges(5,3) * t617;
t663 = t530 * t642 - t535 * t641;
t523 = m(5) * t562 - mrSges(5,2) * t634 + mrSges(5,3) * t586 + t599 * t616 - t607 * t635 + t663;
t516 = t522 * t648 + t523 * t644;
t622 = (mrSges(4,1) * t645 + mrSges(4,2) * t649) * qJD(1);
t671 = qJD(1) * t645;
t625 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t671;
t514 = m(4) * t600 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t624 + qJD(3) * t625 - t622 * t670 + t516;
t626 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t670;
t664 = -t522 * t644 + t523 * t648;
t515 = m(4) * t601 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t623 - qJD(3) * t626 - t622 * t671 + t664;
t511 = t649 * t514 + t645 * t515;
t611 = -qJDD(1) * pkin(1) + t658;
t657 = -m(3) * t611 + (mrSges(3,3) * t651) - t511;
t509 = m(2) * t628 - (t651 * mrSges(2,2)) + t675 * qJDD(1) + t657;
t610 = t651 * pkin(1) - t659;
t608 = t651 * t676 + t659;
t531 = t539 * t647 + t540 * t643;
t655 = m(6) * t556 - mrSges(6,1) * t565 + mrSges(6,2) * t566 - t588 * t597 + t589 * t598 + t531;
t654 = m(5) * t582 - mrSges(5,1) * t586 + mrSges(5,2) * t587 - t605 * t616 + t607 * t617 + t655;
t653 = -m(4) * t608 + mrSges(4,1) * t623 - mrSges(4,2) * t624 - t625 * t671 - t626 * t670 - t654;
t652 = -m(3) * t610 + (mrSges(3,2) * t651) + qJDD(1) * mrSges(3,3) - t653;
t527 = m(2) * t629 - (mrSges(2,1) * t651) - qJDD(1) * mrSges(2,2) + t652;
t672 = t509 * t650 + t527 * t646;
t666 = -t509 * t646 + t527 * t650;
t665 = -t645 * t514 + t515 * t649;
t615 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t649 - Ifges(4,4) * t645) * qJD(1);
t614 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t649 - Ifges(4,2) * t645) * qJD(1);
t613 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t649 - Ifges(4,6) * t645) * qJD(1);
t594 = Ifges(5,1) * t617 + Ifges(5,4) * t616 + Ifges(5,5) * t635;
t593 = Ifges(5,4) * t617 + Ifges(5,2) * t616 + Ifges(5,6) * t635;
t592 = Ifges(5,5) * t617 + Ifges(5,6) * t616 + Ifges(5,3) * t635;
t572 = Ifges(6,1) * t598 + Ifges(6,4) * t597 + Ifges(6,5) * t635;
t571 = Ifges(6,4) * t598 + Ifges(6,2) * t597 + Ifges(6,6) * t635;
t570 = Ifges(6,5) * t598 + Ifges(6,6) * t597 + Ifges(6,3) * t635;
t559 = Ifges(7,1) * t585 + Ifges(7,4) * t584 + Ifges(7,5) * t591;
t558 = Ifges(7,4) * t585 + Ifges(7,2) * t584 + Ifges(7,6) * t591;
t557 = Ifges(7,5) * t585 + Ifges(7,6) * t584 + Ifges(7,3) * t591;
t533 = mrSges(7,2) * t543 - mrSges(7,3) * t541 + Ifges(7,1) * t554 + Ifges(7,4) * t553 + Ifges(7,5) * t564 + t557 * t584 - t558 * t591;
t532 = -mrSges(7,1) * t543 + mrSges(7,3) * t542 + Ifges(7,4) * t554 + Ifges(7,2) * t553 + Ifges(7,6) * t564 - t557 * t585 + t559 * t591;
t518 = -mrSges(6,1) * t556 - mrSges(7,1) * t541 + mrSges(7,2) * t542 + mrSges(6,3) * t546 + Ifges(6,4) * t566 - Ifges(7,5) * t554 + Ifges(6,2) * t565 + Ifges(6,6) * t634 - Ifges(7,6) * t553 - Ifges(7,3) * t564 - pkin(5) * t531 - t558 * t585 + t559 * t584 - t570 * t598 + t572 * t635;
t517 = mrSges(6,2) * t556 - mrSges(6,3) * t545 + Ifges(6,1) * t566 + Ifges(6,4) * t565 + Ifges(6,5) * t634 - pkin(9) * t531 - t532 * t643 + t533 * t647 + t570 * t597 - t571 * t635;
t512 = mrSges(5,2) * t582 - mrSges(5,3) * t561 + Ifges(5,1) * t587 + Ifges(5,4) * t586 + Ifges(5,5) * t634 - qJ(5) * t524 + t517 * t642 - t518 * t641 + t592 * t616 - t593 * t635;
t510 = -m(3) * g(3) + t665;
t507 = -mrSges(5,1) * t582 + mrSges(5,3) * t562 + Ifges(5,4) * t587 + Ifges(5,2) * t586 + Ifges(5,6) * t634 - pkin(4) * t655 + qJ(5) * t663 + t641 * t517 + t642 * t518 - t617 * t592 + t635 * t594;
t506 = mrSges(4,2) * t608 - mrSges(4,3) * t600 + Ifges(4,1) * t624 + Ifges(4,4) * t623 + Ifges(4,5) * qJDD(3) - pkin(8) * t516 - qJD(3) * t614 - t507 * t644 + t512 * t648 - t613 * t671;
t505 = -mrSges(4,1) * t608 + mrSges(4,3) * t601 + Ifges(4,4) * t624 + Ifges(4,2) * t623 + Ifges(4,6) * qJDD(3) - pkin(3) * t654 + pkin(8) * t664 + qJD(3) * t615 + t648 * t507 + t644 * t512 - t613 * t670;
t504 = pkin(5) * t656 + (t673 * t651) + t674 * qJDD(1) + pkin(4) * t524 + pkin(3) * t516 + Ifges(4,3) * qJDD(3) + pkin(2) * t511 - qJ(2) * t510 + pkin(9) * t662 + (Ifges(6,3) + Ifges(5,3)) * t634 + (t614 * t649 + t615 * t645) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t643 * t533 + t647 * t532 - mrSges(2,3) * t628 + Ifges(4,5) * t624 - t616 * t594 + t617 * t593 + Ifges(4,6) * t623 + mrSges(3,1) * t611 + t598 * t571 + mrSges(4,1) * t600 - mrSges(4,2) * t601 - t597 * t572 + Ifges(5,6) * t586 + Ifges(5,5) * t587 + Ifges(6,5) * t566 + Ifges(6,6) * t565 + mrSges(6,1) * t545 - mrSges(6,2) * t546 + mrSges(5,1) * t561 - mrSges(5,2) * t562;
t503 = -mrSges(3,1) * t610 + mrSges(2,3) * t629 - pkin(1) * t510 - pkin(2) * t653 - pkin(7) * t665 + g(3) * t675 - qJDD(1) * t673 - t649 * t505 - t645 * t506 + t651 * t674;
t1 = [-m(1) * g(1) + t666; -m(1) * g(2) + t672; (-m(1) - m(2) - m(3)) * g(3) + t665; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t672 - t503 * t646 + t504 * t650; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t666 + t650 * t503 + t646 * t504; pkin(1) * t657 + qJ(2) * t652 - t645 * t505 - pkin(7) * t511 + mrSges(2,1) * t628 - mrSges(2,2) * t629 + t649 * t506 + mrSges(3,2) * t611 - mrSges(3,3) * t610 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
