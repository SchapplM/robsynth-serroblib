% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 22:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:33:18
% EndTime: 2019-05-04 22:33:27
% DurationCPUTime: 5.54s
% Computational Cost: add. (68781->284), mult. (128264->347), div. (0->0), fcn. (80752->12), ass. (0->130)
t663 = Ifges(5,1) + Ifges(6,2);
t656 = Ifges(5,4) + Ifges(6,6);
t655 = Ifges(5,5) - Ifges(6,4);
t662 = Ifges(5,2) + Ifges(6,3);
t654 = Ifges(5,6) - Ifges(6,5);
t661 = Ifges(5,3) + Ifges(6,1);
t660 = -2 * qJD(5);
t659 = -pkin(4) - pkin(9);
t624 = qJD(2) ^ 2;
t658 = pkin(9) * t624;
t657 = t624 * pkin(8);
t613 = sin(pkin(6));
t619 = sin(qJ(2));
t653 = t613 * t619;
t622 = cos(qJ(2));
t652 = t613 * t622;
t616 = cos(pkin(6));
t651 = t616 * t619;
t650 = t616 * t622;
t612 = sin(pkin(10));
t615 = cos(pkin(10));
t594 = t612 * g(1) - t615 * g(2);
t610 = -g(3) + qJDD(1);
t568 = -t613 * t594 + t616 * t610;
t567 = qJDD(3) + t568;
t621 = cos(qJ(4));
t649 = t621 * t567;
t595 = -t615 * g(1) - t612 * g(2);
t553 = t594 * t650 - t619 * t595 + t610 * t652;
t551 = qJDD(2) * pkin(2) + t553;
t554 = t594 * t651 + t622 * t595 + t610 * t653;
t552 = -t624 * pkin(2) + t554;
t611 = sin(pkin(11));
t614 = cos(pkin(11));
t547 = t611 * t551 + t614 * t552;
t545 = -t624 * pkin(3) + qJDD(2) * pkin(8) + t547;
t618 = sin(qJ(4));
t542 = t618 * t545;
t540 = -t542 + t649;
t589 = (mrSges(6,2) * t621 - mrSges(6,3) * t618) * qJD(2);
t590 = (-mrSges(5,1) * t621 + mrSges(5,2) * t618) * qJD(2);
t641 = qJD(2) * qJD(4);
t638 = t621 * t641;
t591 = t618 * qJDD(2) + t638;
t643 = qJD(2) * t621;
t597 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t643;
t598 = -mrSges(6,1) * t643 - qJD(4) * mrSges(6,3);
t588 = (-pkin(4) * t621 - qJ(5) * t618) * qJD(2);
t623 = qJD(4) ^ 2;
t642 = t618 * qJD(2);
t634 = -t623 * qJ(5) + t588 * t642 + qJDD(5) + t542;
t535 = t591 * pkin(5) + t659 * qJDD(4) + (-pkin(5) * t641 - t618 * t658 - t567) * t621 + t634;
t639 = t618 * t641;
t592 = t621 * qJDD(2) - t639;
t602 = pkin(5) * t642 - qJD(4) * pkin(9);
t609 = t621 ^ 2;
t546 = t614 * t551 - t611 * t552;
t633 = -qJDD(2) * pkin(3) - t546;
t626 = pkin(4) * t639 + t642 * t660 + (-t591 - t638) * qJ(5) + t633;
t536 = -t602 * t642 + (-pkin(5) * t609 - pkin(8)) * t624 + t659 * t592 + t626;
t617 = sin(qJ(6));
t620 = cos(qJ(6));
t531 = t620 * t535 - t617 * t536;
t586 = -t617 * qJD(4) - t620 * t643;
t561 = t586 * qJD(6) + t620 * qJDD(4) - t617 * t592;
t587 = t620 * qJD(4) - t617 * t643;
t562 = -t586 * mrSges(7,1) + t587 * mrSges(7,2);
t604 = qJD(6) + t642;
t565 = -t604 * mrSges(7,2) + t586 * mrSges(7,3);
t584 = qJDD(6) + t591;
t529 = m(7) * t531 + t584 * mrSges(7,1) - t561 * mrSges(7,3) - t587 * t562 + t604 * t565;
t532 = t617 * t535 + t620 * t536;
t560 = -t587 * qJD(6) - t617 * qJDD(4) - t620 * t592;
t566 = t604 * mrSges(7,1) - t587 * mrSges(7,3);
t530 = m(7) * t532 - t584 * mrSges(7,2) + t560 * mrSges(7,3) + t586 * t562 - t604 * t566;
t521 = t620 * t529 + t617 * t530;
t538 = -qJDD(4) * pkin(4) + t634 - t649;
t629 = -m(6) * t538 - t591 * mrSges(6,1) - t521;
t519 = m(5) * t540 - t591 * mrSges(5,3) + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t597 - t598) * qJD(4) + (-t589 - t590) * t642 + t629;
t541 = t621 * t545 + t618 * t567;
t596 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t642;
t628 = -t623 * pkin(4) + qJDD(4) * qJ(5) + t588 * t643 + t541;
t537 = qJD(4) * t660 - t628;
t599 = mrSges(6,1) * t642 + qJD(4) * mrSges(6,2);
t534 = -t609 * t658 + t592 * pkin(5) + ((2 * qJD(5)) + t602) * qJD(4) + t628;
t630 = -m(7) * t534 + t560 * mrSges(7,1) - t561 * mrSges(7,2) + t586 * t565 - t587 * t566;
t627 = -m(6) * t537 + qJDD(4) * mrSges(6,3) + qJD(4) * t599 + t589 * t643 - t630;
t526 = t590 * t643 + m(5) * t541 - qJDD(4) * mrSges(5,2) - qJD(4) * t596 + (mrSges(5,3) + mrSges(6,1)) * t592 + t627;
t635 = -t618 * t519 + t621 * t526;
t512 = m(4) * t547 - t624 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t635;
t544 = t633 - t657;
t539 = -t592 * pkin(4) + t626 - t657;
t647 = -t617 * t529 + t620 * t530;
t632 = -m(6) * t539 - t592 * mrSges(6,2) + t599 * t642 - t647;
t625 = -m(5) * t544 + t597 * t643 + t592 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t591 + (-t596 * t618 - t598 * t621) * qJD(2) + t632;
t517 = m(4) * t546 + qJDD(2) * mrSges(4,1) - t624 * mrSges(4,2) + t625;
t509 = t611 * t512 + t614 * t517;
t507 = m(3) * t553 + qJDD(2) * mrSges(3,1) - t624 * mrSges(3,2) + t509;
t636 = t614 * t512 - t611 * t517;
t508 = m(3) * t554 - t624 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t636;
t515 = t621 * t519 + t618 * t526;
t640 = m(4) * t567 + t515;
t514 = m(3) * t568 + t640;
t494 = t507 * t650 + t508 * t651 - t613 * t514;
t492 = m(2) * t594 + t494;
t498 = -t619 * t507 + t622 * t508;
t497 = m(2) * t595 + t498;
t648 = t615 * t492 + t612 * t497;
t646 = t661 * qJD(4) + (t655 * t618 + t654 * t621) * qJD(2);
t645 = -t654 * qJD(4) + (-t656 * t618 - t662 * t621) * qJD(2);
t644 = t655 * qJD(4) + (t663 * t618 + t656 * t621) * qJD(2);
t493 = t507 * t652 + t508 * t653 + t616 * t514;
t637 = -t612 * t492 + t615 * t497;
t520 = -t591 * mrSges(6,3) + t598 * t643 - t632;
t555 = Ifges(7,5) * t587 + Ifges(7,6) * t586 + Ifges(7,3) * t604;
t557 = Ifges(7,1) * t587 + Ifges(7,4) * t586 + Ifges(7,5) * t604;
t522 = -mrSges(7,1) * t534 + mrSges(7,3) * t532 + Ifges(7,4) * t561 + Ifges(7,2) * t560 + Ifges(7,6) * t584 - t587 * t555 + t604 * t557;
t556 = Ifges(7,4) * t587 + Ifges(7,2) * t586 + Ifges(7,6) * t604;
t523 = mrSges(7,2) * t534 - mrSges(7,3) * t531 + Ifges(7,1) * t561 + Ifges(7,4) * t560 + Ifges(7,5) * t584 + t586 * t555 - t604 * t556;
t500 = -mrSges(5,1) * t544 - mrSges(6,1) * t537 + mrSges(6,2) * t539 + mrSges(5,3) * t541 - pkin(4) * t520 - pkin(5) * t630 - pkin(9) * t647 + t644 * qJD(4) + t654 * qJDD(4) - t620 * t522 - t617 * t523 + t656 * t591 + t662 * t592 - t646 * t642;
t501 = mrSges(6,1) * t538 + mrSges(7,1) * t531 + mrSges(5,2) * t544 - mrSges(7,2) * t532 - mrSges(5,3) * t540 - mrSges(6,3) * t539 + Ifges(7,5) * t561 + Ifges(7,6) * t560 + Ifges(7,3) * t584 + pkin(5) * t521 - qJ(5) * t520 + t587 * t556 - t586 * t557 + t656 * t592 + t663 * t591 + t655 * qJDD(4) + t645 * qJD(4) + t646 * t643;
t490 = mrSges(4,2) * t567 - mrSges(4,3) * t546 + Ifges(4,5) * qJDD(2) - t624 * Ifges(4,6) - pkin(8) * t515 - t618 * t500 + t621 * t501;
t499 = t624 * Ifges(4,5) - pkin(3) * t515 + mrSges(4,3) * t547 - mrSges(4,1) * t567 + Ifges(4,6) * qJDD(2) - pkin(4) * (-qJD(4) * t598 + t629) - qJ(5) * t627 - mrSges(6,2) * t538 + mrSges(6,3) * t537 - t620 * t523 + t617 * t522 + pkin(9) * t521 - mrSges(5,1) * t540 + mrSges(5,2) * t541 + (-qJ(5) * mrSges(6,1) - t654) * t592 - t655 * t591 + (pkin(4) * mrSges(6,2) - t661) * qJDD(4) + (t644 * t621 + (pkin(4) * t589 + t645) * t618) * qJD(2);
t487 = -mrSges(3,1) * t568 + mrSges(3,3) * t554 + t624 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t640 + qJ(3) * t636 + t611 * t490 + t614 * t499;
t488 = mrSges(3,2) * t568 - mrSges(3,3) * t553 + Ifges(3,5) * qJDD(2) - t624 * Ifges(3,6) - qJ(3) * t509 + t614 * t490 - t611 * t499;
t631 = pkin(7) * t498 + t487 * t622 + t488 * t619;
t489 = mrSges(3,1) * t553 - mrSges(3,2) * t554 + mrSges(4,1) * t546 - mrSges(4,2) * t547 + t618 * t501 + t621 * t500 + pkin(3) * t625 + pkin(8) * t635 + pkin(2) * t509 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t486 = mrSges(2,2) * t610 - mrSges(2,3) * t594 - t619 * t487 + t622 * t488 + (-t493 * t613 - t494 * t616) * pkin(7);
t485 = -mrSges(2,1) * t610 + mrSges(2,3) * t595 - pkin(1) * t493 - t613 * t489 + t616 * t631;
t1 = [-m(1) * g(1) + t637; -m(1) * g(2) + t648; -m(1) * g(3) + m(2) * t610 + t493; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t648 - t612 * t485 + t615 * t486; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t637 + t615 * t485 + t612 * t486; -mrSges(1,1) * g(2) + mrSges(2,1) * t594 + mrSges(1,2) * g(1) - mrSges(2,2) * t595 + pkin(1) * t494 + t616 * t489 + t613 * t631;];
tauB  = t1;
