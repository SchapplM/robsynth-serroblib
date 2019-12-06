% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:09:01
% EndTime: 2019-12-05 16:09:05
% DurationCPUTime: 2.31s
% Computational Cost: add. (20496->236), mult. (43824->292), div. (0->0), fcn. (26551->8), ass. (0->99)
t652 = Ifges(5,1) + Ifges(6,1);
t646 = Ifges(5,4) - Ifges(6,5);
t645 = Ifges(5,5) + Ifges(6,4);
t651 = -Ifges(5,2) - Ifges(6,3);
t650 = -Ifges(6,2) - Ifges(5,3);
t644 = Ifges(5,6) - Ifges(6,6);
t614 = sin(pkin(7));
t643 = cos(pkin(7));
t601 = -t643 * g(1) - t614 * g(2);
t612 = -g(3) + qJDD(1);
t616 = sin(qJ(2));
t618 = cos(qJ(2));
t582 = t618 * t601 + t616 * t612;
t620 = qJD(2) ^ 2;
t574 = -t620 * pkin(2) + qJDD(2) * pkin(6) + t582;
t600 = t614 * g(1) - t643 * g(2);
t615 = sin(qJ(3));
t617 = cos(qJ(3));
t553 = -t615 * t574 - t617 * t600;
t634 = qJD(2) * qJD(3);
t630 = t617 * t634;
t598 = t615 * qJDD(2) + t630;
t550 = (-t598 + t630) * qJ(4) + (t615 * t617 * t620 + qJDD(3)) * pkin(3) + t553;
t554 = t617 * t574 - t615 * t600;
t599 = t617 * qJDD(2) - t615 * t634;
t636 = qJD(2) * t615;
t602 = qJD(3) * pkin(3) - qJ(4) * t636;
t611 = t617 ^ 2;
t551 = -t611 * t620 * pkin(3) + t599 * qJ(4) - qJD(3) * t602 + t554;
t613 = sin(pkin(8));
t635 = qJD(2) * t617;
t642 = cos(pkin(8));
t584 = t613 * t636 - t642 * t635;
t648 = -2 * qJD(4);
t547 = t613 * t550 + t642 * t551 + t584 * t648;
t570 = t613 * t598 - t642 * t599;
t585 = (t613 * t617 + t642 * t615) * qJD(2);
t578 = qJD(3) * mrSges(5,1) - t585 * mrSges(5,3);
t564 = t584 * pkin(4) - t585 * qJ(5);
t619 = qJD(3) ^ 2;
t542 = -t619 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t584 * t564 + t547;
t579 = -qJD(3) * mrSges(6,1) + t585 * mrSges(6,2);
t632 = m(6) * t542 + qJDD(3) * mrSges(6,3) + qJD(3) * t579;
t565 = t584 * mrSges(6,1) - t585 * mrSges(6,3);
t637 = -t584 * mrSges(5,1) - t585 * mrSges(5,2) - t565;
t647 = -mrSges(5,3) - mrSges(6,2);
t534 = m(5) * t547 - qJDD(3) * mrSges(5,2) - qJD(3) * t578 + t647 * t570 + t637 * t584 + t632;
t623 = t642 * t550 - t613 * t551;
t546 = t585 * t648 + t623;
t571 = t642 * t598 + t613 * t599;
t577 = -qJD(3) * mrSges(5,2) - t584 * mrSges(5,3);
t543 = -qJDD(3) * pkin(4) - t619 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t564) * t585 - t623;
t580 = -t584 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t626 = -m(6) * t543 + qJDD(3) * mrSges(6,1) + qJD(3) * t580;
t535 = m(5) * t546 + qJDD(3) * mrSges(5,1) + qJD(3) * t577 + t647 * t571 + t637 * t585 + t626;
t528 = t613 * t534 + t642 * t535;
t539 = t571 * mrSges(6,2) + t585 * t565 - t626;
t587 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t615 + Ifges(4,2) * t617) * qJD(2);
t588 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t615 + Ifges(4,4) * t617) * qJD(2);
t638 = t645 * qJD(3) - t646 * t584 + t652 * t585;
t639 = t644 * qJD(3) + t651 * t584 + t646 * t585;
t649 = (t615 * t587 - t617 * t588) * qJD(2) + (Ifges(4,3) - t650) * qJDD(3) - t644 * t570 + t645 * t571 + t638 * t584 + t639 * t585 + mrSges(4,1) * t553 + mrSges(5,1) * t546 - mrSges(6,1) * t543 - mrSges(4,2) * t554 - mrSges(5,2) * t547 + mrSges(6,3) * t542 + Ifges(4,5) * t598 + Ifges(4,6) * t599 + pkin(3) * t528 - pkin(4) * t539 + qJ(5) * (-t570 * mrSges(6,2) - t584 * t565 + t632);
t597 = (-mrSges(4,1) * t617 + mrSges(4,2) * t615) * qJD(2);
t604 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t635;
t524 = m(4) * t553 + qJDD(3) * mrSges(4,1) - t598 * mrSges(4,3) + qJD(3) * t604 - t597 * t636 + t528;
t603 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t636;
t627 = t642 * t534 - t613 * t535;
t525 = m(4) * t554 - qJDD(3) * mrSges(4,2) + t599 * mrSges(4,3) - qJD(3) * t603 + t597 * t635 + t627;
t522 = -t615 * t524 + t617 * t525;
t518 = m(3) * t582 - t620 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t522;
t581 = -t616 * t601 + t618 * t612;
t624 = -qJDD(2) * pkin(2) - t581;
t552 = -t599 * pkin(3) + qJDD(4) + t602 * t636 + (-qJ(4) * t611 - pkin(6)) * t620 + t624;
t545 = -0.2e1 * qJD(5) * t585 + (qJD(3) * t584 - t571) * qJ(5) + (qJD(3) * t585 + t570) * pkin(4) + t552;
t540 = m(6) * t545 + t570 * mrSges(6,1) - t571 * mrSges(6,3) - t585 * t579 + t584 * t580;
t537 = m(5) * t552 + t570 * mrSges(5,1) + t571 * mrSges(5,2) + t584 * t577 + t585 * t578 + t540;
t573 = -t620 * pkin(6) + t624;
t536 = -m(4) * t573 + t599 * mrSges(4,1) - t598 * mrSges(4,2) - t603 * t636 + t604 * t635 - t537;
t533 = m(3) * t581 + qJDD(2) * mrSges(3,1) - t620 * mrSges(3,2) + t536;
t628 = t618 * t518 - t616 * t533;
t514 = m(2) * t601 + t628;
t521 = t617 * t524 + t615 * t525;
t520 = (m(2) + m(3)) * t600 - t521;
t641 = t614 * t514 + t643 * t520;
t515 = t616 * t518 + t618 * t533;
t640 = t650 * qJD(3) + t644 * t584 - t645 * t585;
t631 = m(2) * t612 + t515;
t629 = t643 * t514 - t614 * t520;
t526 = -mrSges(5,1) * t552 - mrSges(6,1) * t545 + mrSges(6,2) * t542 + mrSges(5,3) * t547 - pkin(4) * t540 + t638 * qJD(3) + t644 * qJDD(3) + t651 * t570 + t646 * t571 + t640 * t585;
t527 = mrSges(5,2) * t552 + mrSges(6,2) * t543 - mrSges(5,3) * t546 - mrSges(6,3) * t545 - qJ(5) * t540 - t639 * qJD(3) + t645 * qJDD(3) - t646 * t570 + t652 * t571 + t640 * t584;
t586 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t615 + Ifges(4,6) * t617) * qJD(2);
t510 = -mrSges(4,1) * t573 + mrSges(4,3) * t554 + Ifges(4,4) * t598 + Ifges(4,2) * t599 + Ifges(4,6) * qJDD(3) - pkin(3) * t537 + qJ(4) * t627 + qJD(3) * t588 + t642 * t526 + t613 * t527 - t586 * t636;
t511 = mrSges(4,2) * t573 - mrSges(4,3) * t553 + Ifges(4,1) * t598 + Ifges(4,4) * t599 + Ifges(4,5) * qJDD(3) - qJ(4) * t528 - qJD(3) * t587 - t613 * t526 + t642 * t527 + t586 * t635;
t622 = mrSges(3,1) * t581 - mrSges(3,2) * t582 + Ifges(3,3) * qJDD(2) + pkin(2) * t536 + pkin(6) * t522 + t617 * t510 + t615 * t511;
t509 = mrSges(3,1) * t600 + mrSges(3,3) * t582 + t620 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t521 - t649;
t508 = -mrSges(3,2) * t600 - mrSges(3,3) * t581 + Ifges(3,5) * qJDD(2) - t620 * Ifges(3,6) - pkin(6) * t521 - t615 * t510 + t617 * t511;
t507 = -mrSges(2,1) * t612 + mrSges(2,3) * t601 - pkin(1) * t515 - t622;
t506 = mrSges(2,2) * t612 - mrSges(2,3) * t600 - pkin(5) * t515 + t618 * t508 - t616 * t509;
t1 = [-m(1) * g(1) + t629; -m(1) * g(2) + t641; -m(1) * g(3) + t631; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t641 + t643 * t506 - t614 * t507; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t629 + t614 * t506 + t643 * t507; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t600 - mrSges(2,2) * t601 + t616 * t508 + t618 * t509 + pkin(1) * (m(3) * t600 - t521) + pkin(5) * t628; t631; t622; t649; t537; t539;];
tauJB = t1;
