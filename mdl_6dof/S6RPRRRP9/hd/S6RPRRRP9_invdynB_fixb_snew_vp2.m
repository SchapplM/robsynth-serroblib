% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:51:15
% EndTime: 2019-05-06 01:51:25
% DurationCPUTime: 5.25s
% Computational Cost: add. (58524->316), mult. (114116->370), div. (0->0), fcn. (73408->8), ass. (0->123)
t687 = Ifges(6,1) + Ifges(7,1);
t680 = Ifges(6,4) + Ifges(7,4);
t678 = Ifges(6,5) + Ifges(7,5);
t686 = Ifges(6,2) + Ifges(7,2);
t677 = -Ifges(6,6) - Ifges(7,6);
t685 = -Ifges(6,3) - Ifges(7,3);
t648 = sin(qJ(1));
t652 = cos(qJ(1));
t633 = -g(1) * t652 - g(2) * t648;
t684 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t633;
t683 = -pkin(1) - pkin(7);
t682 = mrSges(2,1) - mrSges(3,2);
t681 = -Ifges(3,4) + Ifges(2,5);
t679 = (Ifges(3,5) - Ifges(2,6));
t632 = g(1) * t648 - t652 * g(2);
t654 = qJD(1) ^ 2;
t660 = -qJ(2) * t654 + qJDD(2) - t632;
t610 = qJDD(1) * t683 + t660;
t647 = sin(qJ(3));
t651 = cos(qJ(3));
t603 = -g(3) * t651 + t647 * t610;
t626 = (mrSges(4,1) * t647 + mrSges(4,2) * t651) * qJD(1);
t671 = qJD(1) * qJD(3);
t636 = t651 * t671;
t628 = -t647 * qJDD(1) - t636;
t672 = qJD(1) * t651;
t631 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t672;
t638 = t647 * qJD(1);
t609 = t654 * t683 - t684;
t667 = t647 * t671;
t629 = qJDD(1) * t651 - t667;
t580 = (-t629 + t667) * pkin(8) + (-t628 + t636) * pkin(3) + t609;
t627 = (pkin(3) * t647 - pkin(8) * t651) * qJD(1);
t653 = qJD(3) ^ 2;
t584 = -pkin(3) * t653 + qJDD(3) * pkin(8) - t627 * t638 + t603;
t646 = sin(qJ(4));
t650 = cos(qJ(4));
t558 = t650 * t580 - t584 * t646;
t624 = qJD(3) * t650 - t646 * t672;
t597 = qJD(4) * t624 + qJDD(3) * t646 + t629 * t650;
t623 = qJDD(4) - t628;
t625 = qJD(3) * t646 + t650 * t672;
t635 = t638 + qJD(4);
t554 = (t624 * t635 - t597) * pkin(9) + (t624 * t625 + t623) * pkin(4) + t558;
t559 = t646 * t580 + t650 * t584;
t596 = -qJD(4) * t625 + qJDD(3) * t650 - t629 * t646;
t608 = pkin(4) * t635 - pkin(9) * t625;
t622 = t624 ^ 2;
t556 = -pkin(4) * t622 + pkin(9) * t596 - t608 * t635 + t559;
t645 = sin(qJ(5));
t649 = cos(qJ(5));
t548 = t649 * t554 - t556 * t645;
t599 = t624 * t649 - t625 * t645;
t565 = qJD(5) * t599 + t596 * t645 + t597 * t649;
t600 = t624 * t645 + t625 * t649;
t576 = -mrSges(7,1) * t599 + mrSges(7,2) * t600;
t577 = -mrSges(6,1) * t599 + mrSges(6,2) * t600;
t634 = qJD(5) + t635;
t586 = -mrSges(6,2) * t634 + mrSges(6,3) * t599;
t618 = qJDD(5) + t623;
t545 = -0.2e1 * qJD(6) * t600 + (t599 * t634 - t565) * qJ(6) + (t599 * t600 + t618) * pkin(5) + t548;
t585 = -mrSges(7,2) * t634 + mrSges(7,3) * t599;
t669 = m(7) * t545 + t618 * mrSges(7,1) + t634 * t585;
t538 = m(6) * t548 + mrSges(6,1) * t618 + t586 * t634 + (-t576 - t577) * t600 + (-mrSges(6,3) - mrSges(7,3)) * t565 + t669;
t549 = t645 * t554 + t649 * t556;
t564 = -qJD(5) * t600 + t596 * t649 - t597 * t645;
t588 = mrSges(7,1) * t634 - mrSges(7,3) * t600;
t589 = mrSges(6,1) * t634 - mrSges(6,3) * t600;
t587 = pkin(5) * t634 - qJ(6) * t600;
t598 = t599 ^ 2;
t547 = -pkin(5) * t598 + qJ(6) * t564 + 0.2e1 * qJD(6) * t599 - t587 * t634 + t549;
t668 = m(7) * t547 + t564 * mrSges(7,3) + t599 * t576;
t541 = m(6) * t549 + mrSges(6,3) * t564 + t577 * t599 + (-t588 - t589) * t634 + (-mrSges(6,2) - mrSges(7,2)) * t618 + t668;
t536 = t649 * t538 + t645 * t541;
t601 = -mrSges(5,1) * t624 + mrSges(5,2) * t625;
t604 = -mrSges(5,2) * t635 + mrSges(5,3) * t624;
t533 = m(5) * t558 + mrSges(5,1) * t623 - mrSges(5,3) * t597 - t601 * t625 + t604 * t635 + t536;
t605 = mrSges(5,1) * t635 - mrSges(5,3) * t625;
t663 = -t538 * t645 + t649 * t541;
t534 = m(5) * t559 - mrSges(5,2) * t623 + mrSges(5,3) * t596 + t601 * t624 - t605 * t635 + t663;
t664 = -t533 * t646 + t650 * t534;
t527 = m(4) * t603 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t628 - qJD(3) * t631 - t626 * t638 + t664;
t602 = g(3) * t647 + t610 * t651;
t630 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t638;
t583 = -qJDD(3) * pkin(3) - pkin(8) * t653 + t627 * t672 - t602;
t557 = -pkin(4) * t596 - pkin(9) * t622 + t625 * t608 + t583;
t551 = -pkin(5) * t564 - qJ(6) * t598 + t587 * t600 + qJDD(6) + t557;
t662 = m(7) * t551 - t564 * mrSges(7,1) + t565 * mrSges(7,2) - t599 * t585 + t600 * t588;
t656 = m(6) * t557 - t564 * mrSges(6,1) + t565 * mrSges(6,2) - t599 * t586 + t600 * t589 + t662;
t655 = -m(5) * t583 + t596 * mrSges(5,1) - t597 * mrSges(5,2) + t624 * t604 - t625 * t605 - t656;
t542 = m(4) * t602 + qJDD(3) * mrSges(4,1) - t629 * mrSges(4,3) + qJD(3) * t630 - t626 * t672 + t655;
t522 = t527 * t647 + t542 * t651;
t612 = -qJDD(1) * pkin(1) + t660;
t659 = -m(3) * t612 + (t654 * mrSges(3,3)) - t522;
t520 = m(2) * t632 - (mrSges(2,2) * t654) + qJDD(1) * t682 + t659;
t611 = pkin(1) * t654 + t684;
t528 = t650 * t533 + t646 * t534;
t658 = -m(4) * t609 + mrSges(4,1) * t628 - t629 * mrSges(4,2) - t630 * t638 - t631 * t672 - t528;
t657 = -m(3) * t611 + (t654 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t658;
t525 = m(2) * t633 - (mrSges(2,1) * t654) - qJDD(1) * mrSges(2,2) + t657;
t676 = t652 * t520 + t648 * t525;
t675 = t599 * t677 - t600 * t678 + t634 * t685;
t674 = -t599 * t686 - t600 * t680 + t634 * t677;
t673 = t680 * t599 + t600 * t687 + t678 * t634;
t666 = -t520 * t648 + t652 * t525;
t665 = t651 * t527 - t542 * t647;
t617 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t651 - Ifges(4,4) * t647) * qJD(1);
t616 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t651 - Ifges(4,2) * t647) * qJD(1);
t615 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t651 - Ifges(4,6) * t647) * qJD(1);
t592 = Ifges(5,1) * t625 + Ifges(5,4) * t624 + Ifges(5,5) * t635;
t591 = Ifges(5,4) * t625 + Ifges(5,2) * t624 + Ifges(5,6) * t635;
t590 = Ifges(5,5) * t625 + Ifges(5,6) * t624 + Ifges(5,3) * t635;
t543 = -mrSges(7,3) * t565 - t576 * t600 + t669;
t535 = mrSges(6,2) * t557 + mrSges(7,2) * t551 - mrSges(6,3) * t548 - mrSges(7,3) * t545 - qJ(6) * t543 + t680 * t564 + t565 * t687 - t675 * t599 + t678 * t618 + t674 * t634;
t529 = -mrSges(6,1) * t557 + mrSges(6,3) * t549 - mrSges(7,1) * t551 + mrSges(7,3) * t547 - pkin(5) * t662 + qJ(6) * t668 + (-qJ(6) * t588 + t673) * t634 + (-qJ(6) * mrSges(7,2) - t677) * t618 + t675 * t600 + t680 * t565 + t686 * t564;
t521 = -m(3) * g(3) + t665;
t518 = mrSges(5,2) * t583 - mrSges(5,3) * t558 + Ifges(5,1) * t597 + Ifges(5,4) * t596 + Ifges(5,5) * t623 - pkin(9) * t536 - t529 * t645 + t535 * t649 + t590 * t624 - t591 * t635;
t517 = -mrSges(5,1) * t583 + mrSges(5,3) * t559 + Ifges(5,4) * t597 + Ifges(5,2) * t596 + Ifges(5,6) * t623 - pkin(4) * t656 + pkin(9) * t663 + t649 * t529 + t645 * t535 - t625 * t590 + t635 * t592;
t516 = t685 * t618 + Ifges(4,6) * qJDD(3) + t677 * t564 - t678 * t565 + t673 * t599 + t674 * t600 - t615 * t672 - t625 * t591 + Ifges(4,2) * t628 + Ifges(4,4) * t629 - Ifges(5,3) * t623 + t624 * t592 - mrSges(4,1) * t609 + qJD(3) * t617 + mrSges(4,3) * t603 - Ifges(5,6) * t596 - Ifges(5,5) * t597 - mrSges(5,1) * t558 + mrSges(5,2) * t559 - mrSges(6,1) * t548 + mrSges(6,2) * t549 - mrSges(7,1) * t545 + mrSges(7,2) * t547 - pkin(5) * t543 - pkin(4) * t536 - pkin(3) * t528;
t515 = mrSges(4,2) * t609 - mrSges(4,3) * t602 + Ifges(4,1) * t629 + Ifges(4,4) * t628 + Ifges(4,5) * qJDD(3) - pkin(8) * t528 - qJD(3) * t616 - t517 * t646 + t518 * t650 - t615 * t638;
t514 = -qJ(2) * t521 - mrSges(2,3) * t632 + pkin(2) * t522 + mrSges(3,1) * t612 + t650 * t517 + pkin(3) * t655 + pkin(8) * t664 + t646 * t518 + Ifges(4,5) * t629 + Ifges(4,6) * t628 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t602 - mrSges(4,2) * t603 + (t679 * t654) + t681 * qJDD(1) + (t651 * t616 + t647 * t617) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t513 = -mrSges(3,1) * t611 + mrSges(2,3) * t633 - pkin(1) * t521 - pkin(2) * t658 - pkin(7) * t665 + g(3) * t682 - qJDD(1) * t679 - t647 * t515 - t651 * t516 + t654 * t681;
t1 = [-m(1) * g(1) + t666; -m(1) * g(2) + t676; (-m(1) - m(2) - m(3)) * g(3) + t665; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t676 - t648 * t513 + t652 * t514; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t666 + t652 * t513 + t648 * t514; pkin(1) * t659 + qJ(2) * t657 + t651 * t515 - t647 * t516 - pkin(7) * t522 + mrSges(2,1) * t632 - mrSges(2,2) * t633 + mrSges(3,2) * t612 - mrSges(3,3) * t611 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
