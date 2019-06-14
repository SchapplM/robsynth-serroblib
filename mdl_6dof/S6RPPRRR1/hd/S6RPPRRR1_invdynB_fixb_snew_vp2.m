% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:10:21
% EndTime: 2019-05-05 15:10:32
% DurationCPUTime: 11.57s
% Computational Cost: add. (189212->322), mult. (429744->401), div. (0->0), fcn. (317092->12), ass. (0->136)
t660 = qJD(1) ^ 2;
t650 = cos(pkin(11));
t687 = pkin(3) * t650;
t648 = sin(pkin(11));
t686 = mrSges(4,2) * t648;
t645 = t650 ^ 2;
t685 = t645 * t660;
t655 = sin(qJ(1));
t659 = cos(qJ(1));
t631 = t655 * g(1) - t659 * g(2);
t629 = qJDD(1) * pkin(1) + t631;
t632 = -t659 * g(1) - t655 * g(2);
t630 = -t660 * pkin(1) + t632;
t649 = sin(pkin(10));
t651 = cos(pkin(10));
t617 = t649 * t629 + t651 * t630;
t609 = -t660 * pkin(2) + qJDD(1) * qJ(3) + t617;
t647 = -g(3) + qJDD(2);
t680 = qJD(1) * qJD(3);
t683 = t650 * t647 - 0.2e1 * t648 * t680;
t592 = (-pkin(7) * qJDD(1) + t660 * t687 - t609) * t648 + t683;
t596 = t648 * t647 + (t609 + 0.2e1 * t680) * t650;
t679 = qJDD(1) * t650;
t593 = -pkin(3) * t685 + pkin(7) * t679 + t596;
t654 = sin(qJ(4));
t658 = cos(qJ(4));
t571 = t658 * t592 - t654 * t593;
t667 = t648 * t658 + t650 * t654;
t666 = -t648 * t654 + t650 * t658;
t622 = t666 * qJD(1);
t681 = t622 * qJD(4);
t615 = t667 * qJDD(1) + t681;
t623 = t667 * qJD(1);
t560 = (-t615 + t681) * pkin(8) + (t622 * t623 + qJDD(4)) * pkin(4) + t571;
t572 = t654 * t592 + t658 * t593;
t614 = -t623 * qJD(4) + t666 * qJDD(1);
t620 = qJD(4) * pkin(4) - t623 * pkin(8);
t621 = t622 ^ 2;
t562 = -t621 * pkin(4) + t614 * pkin(8) - qJD(4) * t620 + t572;
t653 = sin(qJ(5));
t657 = cos(qJ(5));
t558 = t653 * t560 + t657 * t562;
t608 = t653 * t622 + t657 * t623;
t577 = -t608 * qJD(5) + t657 * t614 - t653 * t615;
t607 = t657 * t622 - t653 * t623;
t587 = -t607 * mrSges(6,1) + t608 * mrSges(6,2);
t646 = qJD(4) + qJD(5);
t600 = t646 * mrSges(6,1) - t608 * mrSges(6,3);
t643 = qJDD(4) + qJDD(5);
t588 = -t607 * pkin(5) - t608 * pkin(9);
t642 = t646 ^ 2;
t555 = -t642 * pkin(5) + t643 * pkin(9) + t607 * t588 + t558;
t644 = t648 ^ 2;
t616 = t651 * t629 - t649 * t630;
t668 = qJDD(3) - t616;
t594 = (-pkin(2) - t687) * qJDD(1) + (-qJ(3) + (-t644 - t645) * pkin(7)) * t660 + t668;
t567 = -t614 * pkin(4) - t621 * pkin(8) + t623 * t620 + t594;
t578 = t607 * qJD(5) + t653 * t614 + t657 * t615;
t556 = (-t607 * t646 - t578) * pkin(9) + (t608 * t646 - t577) * pkin(5) + t567;
t652 = sin(qJ(6));
t656 = cos(qJ(6));
t552 = -t652 * t555 + t656 * t556;
t597 = -t652 * t608 + t656 * t646;
t565 = t597 * qJD(6) + t656 * t578 + t652 * t643;
t576 = qJDD(6) - t577;
t598 = t656 * t608 + t652 * t646;
t579 = -t597 * mrSges(7,1) + t598 * mrSges(7,2);
t602 = qJD(6) - t607;
t580 = -t602 * mrSges(7,2) + t597 * mrSges(7,3);
t550 = m(7) * t552 + t576 * mrSges(7,1) - t565 * mrSges(7,3) - t598 * t579 + t602 * t580;
t553 = t656 * t555 + t652 * t556;
t564 = -t598 * qJD(6) - t652 * t578 + t656 * t643;
t581 = t602 * mrSges(7,1) - t598 * mrSges(7,3);
t551 = m(7) * t553 - t576 * mrSges(7,2) + t564 * mrSges(7,3) + t597 * t579 - t602 * t581;
t672 = -t652 * t550 + t656 * t551;
t541 = m(6) * t558 - t643 * mrSges(6,2) + t577 * mrSges(6,3) + t607 * t587 - t646 * t600 + t672;
t557 = t657 * t560 - t653 * t562;
t599 = -t646 * mrSges(6,2) + t607 * mrSges(6,3);
t554 = -t643 * pkin(5) - t642 * pkin(9) + t608 * t588 - t557;
t663 = -m(7) * t554 + t564 * mrSges(7,1) - t565 * mrSges(7,2) + t597 * t580 - t598 * t581;
t546 = m(6) * t557 + t643 * mrSges(6,1) - t578 * mrSges(6,3) - t608 * t587 + t646 * t599 + t663;
t536 = t653 * t541 + t657 * t546;
t612 = -t622 * mrSges(5,1) + t623 * mrSges(5,2);
t618 = -qJD(4) * mrSges(5,2) + t622 * mrSges(5,3);
t534 = m(5) * t571 + qJDD(4) * mrSges(5,1) - t615 * mrSges(5,3) + qJD(4) * t618 - t623 * t612 + t536;
t619 = qJD(4) * mrSges(5,1) - t623 * mrSges(5,3);
t673 = t657 * t541 - t653 * t546;
t535 = m(5) * t572 - qJDD(4) * mrSges(5,2) + t614 * mrSges(5,3) - qJD(4) * t619 + t622 * t612 + t673;
t528 = t658 * t534 + t654 * t535;
t595 = -t648 * t609 + t683;
t665 = mrSges(4,3) * qJDD(1) + t660 * (-mrSges(4,1) * t650 + t686);
t526 = m(4) * t595 - t665 * t648 + t528;
t674 = -t654 * t534 + t658 * t535;
t527 = m(4) * t596 + t665 * t650 + t674;
t675 = -t648 * t526 + t650 * t527;
t520 = m(3) * t617 - t660 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t675;
t603 = -qJDD(1) * pkin(2) - t660 * qJ(3) + t668;
t542 = t656 * t550 + t652 * t551;
t664 = m(6) * t567 - t577 * mrSges(6,1) + t578 * mrSges(6,2) - t607 * t599 + t608 * t600 + t542;
t662 = m(5) * t594 - t614 * mrSges(5,1) + t615 * mrSges(5,2) - t622 * t618 + t623 * t619 + t664;
t661 = -m(4) * t603 + mrSges(4,1) * t679 - t662 + (t644 * t660 + t685) * mrSges(4,3);
t538 = t661 + m(3) * t616 - t660 * mrSges(3,2) + (mrSges(3,1) - t686) * qJDD(1);
t516 = t649 * t520 + t651 * t538;
t514 = m(2) * t631 + qJDD(1) * mrSges(2,1) - t660 * mrSges(2,2) + t516;
t676 = t651 * t520 - t649 * t538;
t515 = m(2) * t632 - t660 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t676;
t684 = t659 * t514 + t655 * t515;
t521 = t650 * t526 + t648 * t527;
t669 = Ifges(4,5) * t648 + Ifges(4,6) * t650;
t682 = t660 * t669;
t678 = m(3) * t647 + t521;
t677 = -t655 * t514 + t659 * t515;
t671 = Ifges(4,1) * t648 + Ifges(4,4) * t650;
t670 = Ifges(4,4) * t648 + Ifges(4,2) * t650;
t606 = Ifges(5,1) * t623 + Ifges(5,4) * t622 + Ifges(5,5) * qJD(4);
t605 = Ifges(5,4) * t623 + Ifges(5,2) * t622 + Ifges(5,6) * qJD(4);
t604 = Ifges(5,5) * t623 + Ifges(5,6) * t622 + Ifges(5,3) * qJD(4);
t584 = Ifges(6,1) * t608 + Ifges(6,4) * t607 + Ifges(6,5) * t646;
t583 = Ifges(6,4) * t608 + Ifges(6,2) * t607 + Ifges(6,6) * t646;
t582 = Ifges(6,5) * t608 + Ifges(6,6) * t607 + Ifges(6,3) * t646;
t570 = Ifges(7,1) * t598 + Ifges(7,4) * t597 + Ifges(7,5) * t602;
t569 = Ifges(7,4) * t598 + Ifges(7,2) * t597 + Ifges(7,6) * t602;
t568 = Ifges(7,5) * t598 + Ifges(7,6) * t597 + Ifges(7,3) * t602;
t544 = mrSges(7,2) * t554 - mrSges(7,3) * t552 + Ifges(7,1) * t565 + Ifges(7,4) * t564 + Ifges(7,5) * t576 + t597 * t568 - t602 * t569;
t543 = -mrSges(7,1) * t554 + mrSges(7,3) * t553 + Ifges(7,4) * t565 + Ifges(7,2) * t564 + Ifges(7,6) * t576 - t598 * t568 + t602 * t570;
t530 = -mrSges(6,1) * t567 - mrSges(7,1) * t552 + mrSges(7,2) * t553 + mrSges(6,3) * t558 + Ifges(6,4) * t578 - Ifges(7,5) * t565 + Ifges(6,2) * t577 + Ifges(6,6) * t643 - Ifges(7,6) * t564 - Ifges(7,3) * t576 - pkin(5) * t542 - t598 * t569 + t597 * t570 - t608 * t582 + t646 * t584;
t529 = mrSges(6,2) * t567 - mrSges(6,3) * t557 + Ifges(6,1) * t578 + Ifges(6,4) * t577 + Ifges(6,5) * t643 - pkin(9) * t542 - t652 * t543 + t656 * t544 + t607 * t582 - t646 * t583;
t522 = mrSges(5,2) * t594 - mrSges(5,3) * t571 + Ifges(5,1) * t615 + Ifges(5,4) * t614 + Ifges(5,5) * qJDD(4) - pkin(8) * t536 - qJD(4) * t605 + t657 * t529 - t653 * t530 + t622 * t604;
t517 = -mrSges(5,1) * t594 + mrSges(5,3) * t572 + Ifges(5,4) * t615 + Ifges(5,2) * t614 + Ifges(5,6) * qJDD(4) - pkin(4) * t664 + pkin(8) * t673 + qJD(4) * t606 + t653 * t529 + t657 * t530 - t623 * t604;
t510 = (Ifges(3,6) - t669) * qJDD(1) - pkin(9) * t672 + (-t648 * t670 + t650 * t671 + Ifges(3,5)) * t660 - pkin(5) * t663 - mrSges(4,1) * t595 + mrSges(4,2) * t596 - pkin(3) * t528 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t571 + mrSges(5,2) * t572 - Ifges(6,6) * t577 - pkin(2) * t521 - mrSges(6,1) * t557 + mrSges(6,2) * t558 - Ifges(5,6) * t614 - Ifges(5,5) * t615 + mrSges(3,3) * t617 - t656 * t543 + t607 * t584 - t608 * t583 - Ifges(6,3) * t643 - Ifges(6,5) * t578 + t622 * t606 - t623 * t605 - pkin(4) * t536 - mrSges(3,1) * t647 - t652 * t544;
t509 = mrSges(4,2) * t603 - mrSges(4,3) * t595 - pkin(7) * t528 + t671 * qJDD(1) - t654 * t517 + t658 * t522 + t650 * t682;
t508 = -mrSges(4,1) * t603 + mrSges(4,3) * t596 - pkin(3) * t662 + pkin(7) * t674 + t670 * qJDD(1) + t658 * t517 + t654 * t522 - t648 * t682;
t507 = mrSges(3,2) * t647 - mrSges(3,3) * t616 + Ifges(3,5) * qJDD(1) - t660 * Ifges(3,6) - qJ(3) * t521 - t648 * t508 + t650 * t509;
t506 = -mrSges(2,2) * g(3) - mrSges(2,3) * t631 + Ifges(2,5) * qJDD(1) - t660 * Ifges(2,6) - qJ(2) * t516 + t651 * t507 - t649 * t510;
t505 = mrSges(2,1) * g(3) + mrSges(2,3) * t632 + t660 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t678 + qJ(2) * t676 + t649 * t507 + t651 * t510;
t1 = [-m(1) * g(1) + t677; -m(1) * g(2) + t684; (-m(1) - m(2)) * g(3) + t678; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t684 - t655 * t505 + t659 * t506; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t677 + t659 * t505 + t655 * t506; pkin(1) * t516 + mrSges(2,1) * t631 - mrSges(2,2) * t632 + t650 * t508 + pkin(2) * t661 + qJ(3) * t675 + t648 * t509 + mrSges(3,1) * t616 - mrSges(3,2) * t617 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t686 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
