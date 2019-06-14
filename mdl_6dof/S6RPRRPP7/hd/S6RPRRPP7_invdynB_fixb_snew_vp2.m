% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:46:31
% EndTime: 2019-05-05 21:46:36
% DurationCPUTime: 2.30s
% Computational Cost: add. (19569->293), mult. (36945->330), div. (0->0), fcn. (20710->6), ass. (0->112)
t686 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t665 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t664 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t685 = Ifges(5,2) + Ifges(6,3) + Ifges(7,2);
t663 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t684 = -Ifges(5,3) - Ifges(6,2) - Ifges(7,3);
t637 = sin(qJ(1));
t639 = cos(qJ(1));
t623 = -t639 * g(1) - t637 * g(2);
t683 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t623;
t635 = sin(qJ(4));
t638 = cos(qJ(3));
t669 = qJD(1) * t638;
t678 = cos(qJ(4));
t613 = -qJD(3) * t678 + t635 * t669;
t636 = sin(qJ(3));
t667 = qJD(1) * qJD(3);
t658 = t636 * t667;
t618 = t638 * qJDD(1) - t658;
t578 = -t613 * qJD(4) + t635 * qJDD(3) + t618 * t678;
t622 = t637 * g(1) - t639 * g(2);
t641 = qJD(1) ^ 2;
t648 = -t641 * qJ(2) + qJDD(2) - t622;
t679 = -pkin(1) - pkin(7);
t598 = qJDD(1) * t679 + t648;
t585 = t636 * g(3) + t638 * t598;
t616 = (pkin(3) * t636 - pkin(8) * t638) * qJD(1);
t640 = qJD(3) ^ 2;
t647 = qJDD(3) * pkin(3) + t640 * pkin(8) - t616 * t669 + t585;
t668 = t636 * qJD(1);
t625 = qJD(4) + t668;
t673 = t613 * t625;
t682 = (-t578 + t673) * qJ(5) - t647;
t614 = t635 * qJD(3) + t669 * t678;
t681 = -0.2e1 * t614;
t680 = 2 * qJD(5);
t677 = mrSges(2,1) - mrSges(3,2);
t676 = -mrSges(5,3) - mrSges(6,2);
t675 = -Ifges(3,4) + Ifges(2,5);
t674 = (Ifges(3,5) - Ifges(2,6));
t587 = t625 * mrSges(7,2) + t613 * mrSges(7,3);
t672 = t625 * t587;
t586 = -t638 * g(3) + t636 * t598;
t615 = (mrSges(4,1) * t636 + mrSges(4,2) * t638) * qJD(1);
t657 = t638 * t667;
t617 = -t636 * qJDD(1) - t657;
t620 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t669;
t597 = t641 * t679 - t683;
t549 = (-t618 + t658) * pkin(8) + (-t617 + t657) * pkin(3) + t597;
t553 = -t640 * pkin(3) + qJDD(3) * pkin(8) - t616 * t668 + t586;
t546 = t549 * t678 - t635 * t553;
t588 = -t625 * mrSges(5,2) - t613 * mrSges(5,3);
t612 = qJDD(4) - t617;
t581 = t613 * pkin(4) - t614 * qJ(5);
t624 = t625 ^ 2;
t544 = -t612 * pkin(4) - t624 * qJ(5) + t614 * t581 + qJDD(5) - t546;
t593 = -t613 * mrSges(6,2) + t625 * mrSges(6,3);
t537 = qJD(6) * t681 + (-t578 - t673) * qJ(6) + (t613 * t614 - t612) * pkin(5) + t544;
t583 = -t613 * mrSges(7,1) + t614 * mrSges(7,2);
t652 = -m(7) * t537 + t578 * mrSges(7,3) + t614 * t583;
t644 = -m(6) * t544 + t612 * mrSges(6,1) + t625 * t593 + t652;
t582 = t613 * mrSges(6,1) - t614 * mrSges(6,3);
t670 = -t613 * mrSges(5,1) - t614 * mrSges(5,2) - t582;
t533 = m(5) * t546 + (t587 + t588) * t625 + t670 * t614 + (mrSges(5,1) + mrSges(7,1)) * t612 + t676 * t578 + t644;
t547 = t635 * t549 + t678 * t553;
t577 = t614 * qJD(4) - qJDD(3) * t678 + t635 * t618;
t590 = -t625 * mrSges(7,1) - t614 * mrSges(7,3);
t591 = t625 * mrSges(5,1) - t614 * mrSges(5,3);
t543 = -t624 * pkin(4) + t612 * qJ(5) - t613 * t581 + t625 * t680 + t547;
t592 = -t625 * mrSges(6,1) + t614 * mrSges(6,2);
t589 = -t625 * pkin(5) - t614 * qJ(6);
t611 = t613 ^ 2;
t539 = -t611 * pkin(5) + t577 * qJ(6) + 0.2e1 * qJD(6) * t613 + t625 * t589 + t543;
t662 = m(7) * t539 + t577 * mrSges(7,3) + t613 * t583;
t649 = m(6) * t543 + t612 * mrSges(6,3) + t625 * t592 + t662;
t534 = m(5) * t547 + (t590 - t591) * t625 + t670 * t613 + (-mrSges(5,2) + mrSges(7,2)) * t612 + t676 * t577 + t649;
t654 = -t635 * t533 + t678 * t534;
t527 = m(4) * t586 - qJDD(3) * mrSges(4,2) + t617 * mrSges(4,3) - qJD(3) * t620 - t615 * t668 + t654;
t619 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t668;
t545 = qJD(5) * t681 + (t614 * t625 + t577) * pkin(4) + t682;
t541 = -t611 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t577 + (-pkin(4) * t625 + t589 + t680) * t614 - t682;
t651 = -m(7) * t541 + t577 * mrSges(7,1) - t578 * mrSges(7,2) + t613 * t587 - t614 * t590;
t535 = m(6) * t545 + t577 * mrSges(6,1) - t578 * mrSges(6,3) - t614 * t592 + t613 * t593 + t651;
t642 = m(5) * t647 - t577 * mrSges(5,1) - t578 * mrSges(5,2) - t613 * t588 - t614 * t591 - t535;
t529 = m(4) * t585 + qJDD(3) * mrSges(4,1) - t618 * mrSges(4,3) + qJD(3) * t619 - t615 * t669 + t642;
t521 = t636 * t527 + t638 * t529;
t600 = -qJDD(1) * pkin(1) + t648;
t646 = -m(3) * t600 + (t641 * mrSges(3,3)) - t521;
t519 = m(2) * t622 - (t641 * mrSges(2,2)) + qJDD(1) * t677 + t646;
t599 = t641 * pkin(1) + t683;
t528 = t678 * t533 + t635 * t534;
t645 = -m(4) * t597 + t617 * mrSges(4,1) - t618 * mrSges(4,2) - t619 * t668 - t620 * t669 - t528;
t643 = -m(3) * t599 + (t641 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t645;
t525 = m(2) * t623 - (t641 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t643;
t671 = t639 * t519 + t637 * t525;
t661 = t613 * t663 - t614 * t664 + t625 * t684;
t660 = t613 * t685 - t614 * t665 - t625 * t663;
t659 = t665 * t613 - t614 * t686 - t664 * t625;
t656 = -t637 * t519 + t639 * t525;
t655 = t638 * t527 - t636 * t529;
t604 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t638 - Ifges(4,4) * t636) * qJD(1);
t603 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t638 - Ifges(4,2) * t636) * qJD(1);
t602 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t638 - Ifges(4,6) * t636) * qJD(1);
t536 = -t612 * mrSges(7,1) - t652 - t672;
t522 = -mrSges(5,2) * t647 + mrSges(6,2) * t544 + mrSges(7,2) * t541 - mrSges(5,3) * t546 - mrSges(6,3) * t545 - mrSges(7,3) * t537 - qJ(5) * t535 - qJ(6) * t536 - t665 * t577 + t578 * t686 + t664 * t612 + t661 * t613 + t660 * t625;
t520 = -m(3) * g(3) + t655;
t517 = mrSges(5,1) * t647 + mrSges(5,3) * t547 - mrSges(6,1) * t545 + mrSges(6,2) * t543 + mrSges(7,1) * t541 - mrSges(7,3) * t539 - pkin(5) * t651 - qJ(6) * t662 - pkin(4) * t535 + (-qJ(6) * t590 - t659) * t625 + t661 * t614 + (-qJ(6) * mrSges(7,2) + t663) * t612 + t665 * t578 - t685 * t577;
t516 = -pkin(4) * (t644 + t672) - qJ(5) * (t625 * t590 + t649) + Ifges(4,6) * qJDD(3) + Ifges(4,2) * t617 + Ifges(4,4) * t618 + mrSges(4,3) * t586 - mrSges(4,1) * t597 + qJD(3) * t604 - mrSges(5,1) * t546 + mrSges(5,2) * t547 - mrSges(6,3) * t543 + mrSges(6,1) * t544 + pkin(5) * t536 + mrSges(7,1) * t537 - mrSges(7,2) * t539 - pkin(3) * t528 - t602 * t669 + (pkin(4) * t582 + t660) * t614 + (qJ(5) * t582 + t659) * t613 + (-pkin(4) * mrSges(7,1) - qJ(5) * mrSges(7,2) + t684) * t612 + (pkin(4) * mrSges(6,2) - t664) * t578 + (qJ(5) * mrSges(6,2) + t663) * t577;
t515 = mrSges(4,2) * t597 - mrSges(4,3) * t585 + Ifges(4,1) * t618 + Ifges(4,4) * t617 + Ifges(4,5) * qJDD(3) - pkin(8) * t528 - qJD(3) * t603 - t635 * t517 + t522 * t678 - t602 * t668;
t514 = -qJ(2) * t520 - mrSges(2,3) * t622 + pkin(2) * t521 + mrSges(3,1) * t600 + t678 * t517 + pkin(3) * t642 + pkin(8) * t654 + t635 * t522 + Ifges(4,5) * t618 + Ifges(4,6) * t617 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t585 - mrSges(4,2) * t586 + (t674 * t641) + t675 * qJDD(1) + (t638 * t603 + t636 * t604) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t513 = -mrSges(3,1) * t599 + mrSges(2,3) * t623 - pkin(1) * t520 - pkin(2) * t645 - pkin(7) * t655 + g(3) * t677 - qJDD(1) * t674 - t636 * t515 - t638 * t516 + t641 * t675;
t1 = [-m(1) * g(1) + t656; -m(1) * g(2) + t671; (-m(1) - m(2) - m(3)) * g(3) + t655; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t671 - t637 * t513 + t639 * t514; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t656 + t639 * t513 + t637 * t514; qJ(2) * t643 + pkin(1) * t646 - t636 * t516 - pkin(7) * t521 + mrSges(2,1) * t622 - mrSges(2,2) * t623 + t638 * t515 + mrSges(3,2) * t600 - mrSges(3,3) * t599 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
