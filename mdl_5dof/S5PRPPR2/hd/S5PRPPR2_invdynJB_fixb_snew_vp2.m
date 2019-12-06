% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:14
% EndTime: 2019-12-05 15:24:17
% DurationCPUTime: 2.19s
% Computational Cost: add. (22559->195), mult. (42594->248), div. (0->0), fcn. (27144->10), ass. (0->94)
t544 = qJD(2) ^ 2;
t536 = sin(pkin(7));
t539 = cos(pkin(7));
t522 = -g(1) * t539 - g(2) * t536;
t533 = -g(3) + qJDD(1);
t541 = sin(qJ(2));
t543 = cos(qJ(2));
t511 = -t522 * t541 + t543 * t533;
t507 = qJDD(2) * pkin(2) + t511;
t512 = t543 * t522 + t541 * t533;
t508 = -pkin(2) * t544 + t512;
t535 = sin(pkin(8));
t538 = cos(pkin(8));
t495 = t535 * t507 + t538 * t508;
t493 = -pkin(3) * t544 + qJDD(2) * qJ(4) + t495;
t534 = sin(pkin(9));
t521 = g(1) * t536 - t539 * g(2);
t520 = qJDD(3) - t521;
t537 = cos(pkin(9));
t562 = qJD(2) * qJD(4);
t564 = t537 * t520 - 0.2e1 * t534 * t562;
t569 = pkin(4) * t537;
t486 = (-pkin(6) * qJDD(2) + t544 * t569 - t493) * t534 + t564;
t489 = t534 * t520 + (t493 + 0.2e1 * t562) * t537;
t561 = qJDD(2) * t537;
t532 = t537 ^ 2;
t566 = t532 * t544;
t487 = -pkin(4) * t566 + pkin(6) * t561 + t489;
t540 = sin(qJ(5));
t542 = cos(qJ(5));
t485 = t486 * t540 + t487 * t542;
t494 = t507 * t538 - t535 * t508;
t552 = qJDD(4) - t494;
t570 = t534 ^ 2;
t490 = (-pkin(3) - t569) * qJDD(2) + (-qJ(4) + (-t532 - t570) * pkin(6)) * t544 + t552;
t550 = -t534 * t540 + t537 * t542;
t513 = t550 * qJD(2);
t551 = t534 * t542 + t537 * t540;
t514 = t551 * qJD(2);
t496 = Ifges(6,5) * t514 + Ifges(6,6) * t513 + Ifges(6,3) * qJD(5);
t498 = Ifges(6,1) * t514 + Ifges(6,4) * t513 + Ifges(6,5) * qJD(5);
t502 = -t514 * qJD(5) + t550 * qJDD(2);
t503 = t513 * qJD(5) + t551 * qJDD(2);
t474 = -mrSges(6,1) * t490 + mrSges(6,3) * t485 + Ifges(6,4) * t503 + Ifges(6,2) * t502 + Ifges(6,6) * qJDD(5) + qJD(5) * t498 - t496 * t514;
t484 = t486 * t542 - t487 * t540;
t497 = Ifges(6,4) * t514 + Ifges(6,2) * t513 + Ifges(6,6) * qJD(5);
t475 = mrSges(6,2) * t490 - mrSges(6,3) * t484 + Ifges(6,1) * t503 + Ifges(6,4) * t502 + Ifges(6,5) * qJDD(5) - qJD(5) * t497 + t496 * t513;
t492 = -qJDD(2) * pkin(3) - qJ(4) * t544 + t552;
t509 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t513;
t510 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t514;
t548 = m(6) * t490 - t502 * mrSges(6,1) + t503 * mrSges(6,2) - t513 * t509 + t514 * t510;
t554 = Ifges(5,4) * t534 + Ifges(5,2) * t537;
t500 = -mrSges(6,1) * t513 + mrSges(6,2) * t514;
t481 = m(6) * t484 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t503 + qJD(5) * t509 - t500 * t514;
t482 = m(6) * t485 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t502 - qJD(5) * t510 + t500 * t513;
t556 = -t481 * t540 + t542 * t482;
t553 = Ifges(5,5) * t534 + Ifges(5,6) * t537;
t563 = t544 * t553;
t457 = -mrSges(5,1) * t492 + mrSges(5,3) * t489 - pkin(4) * t548 + pkin(6) * t556 + t554 * qJDD(2) + t542 * t474 + t540 * t475 - t534 * t563;
t473 = t542 * t481 + t540 * t482;
t488 = -t493 * t534 + t564;
t567 = mrSges(5,2) * t534;
t549 = mrSges(5,3) * qJDD(2) + t544 * (-mrSges(5,1) * t537 + t567);
t471 = m(5) * t488 - t549 * t534 + t473;
t472 = m(5) * t489 + t549 * t537 + t556;
t467 = -t471 * t534 + t537 * t472;
t462 = m(4) * t495 - mrSges(4,1) * t544 - qJDD(2) * mrSges(4,2) + t467;
t547 = -m(5) * t492 + mrSges(5,1) * t561 - t548 + (t544 * t570 + t566) * mrSges(5,3);
t477 = m(4) * t494 - mrSges(4,2) * t544 + (mrSges(4,1) - t567) * qJDD(2) + t547;
t458 = t535 * t462 + t538 * t477;
t555 = Ifges(5,1) * t534 + Ifges(5,4) * t537;
t459 = mrSges(5,2) * t492 - mrSges(5,3) * t488 - pkin(6) * t473 + t555 * qJDD(2) - t474 * t540 + t475 * t542 + t537 * t563;
t483 = qJDD(2) * t567 - t547;
t571 = mrSges(3,1) * t511 + mrSges(4,1) * t494 - mrSges(3,2) * t512 - mrSges(4,2) * t495 + pkin(2) * t458 - pkin(3) * t483 + qJ(4) * t467 + t457 * t537 + t459 * t534 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t455 = m(3) * t511 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t544 + t458;
t557 = t538 * t462 - t535 * t477;
t456 = m(3) * t512 - mrSges(3,1) * t544 - qJDD(2) * mrSges(3,2) + t557;
t558 = -t455 * t541 + t543 * t456;
t449 = m(2) * t522 + t558;
t466 = t537 * t471 + t534 * t472;
t465 = m(4) * t520 + t466;
t464 = (m(2) + m(3)) * t521 - t465;
t565 = t536 * t449 + t539 * t464;
t450 = t543 * t455 + t541 * t456;
t560 = m(2) * t533 + t450;
t559 = t539 * t449 - t464 * t536;
t546 = mrSges(6,1) * t484 - mrSges(6,2) * t485 + Ifges(6,5) * t503 + Ifges(6,6) * t502 + Ifges(6,3) * qJDD(5) + t514 * t497 - t513 * t498;
t451 = -mrSges(4,1) * t520 - mrSges(5,1) * t488 + mrSges(5,2) * t489 + mrSges(4,3) * t495 - pkin(3) * t466 - pkin(4) * t473 + (Ifges(4,6) - t553) * qJDD(2) - t546 + (-t534 * t554 + t537 * t555 + Ifges(4,5)) * t544;
t446 = mrSges(4,2) * t520 - mrSges(4,3) * t494 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t544 - qJ(4) * t466 - t457 * t534 + t459 * t537;
t445 = -mrSges(3,2) * t521 - mrSges(3,3) * t511 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t544 - qJ(3) * t458 + t446 * t538 - t451 * t535;
t444 = -mrSges(2,1) * t533 + mrSges(2,3) * t522 - pkin(1) * t450 - t571;
t443 = mrSges(3,1) * t521 + mrSges(3,3) * t512 + t544 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t465 + qJ(3) * t557 + t535 * t446 + t538 * t451;
t442 = mrSges(2,2) * t533 - mrSges(2,3) * t521 - pkin(5) * t450 - t443 * t541 + t445 * t543;
t1 = [-m(1) * g(1) + t559; -m(1) * g(2) + t565; -m(1) * g(3) + t560; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t565 + t539 * t442 - t536 * t444; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t559 + t536 * t442 + t539 * t444; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t521 - mrSges(2,2) * t522 + t541 * t445 + t543 * t443 + pkin(1) * (m(3) * t521 - t465) + pkin(5) * t558; t560; t571; t465; t483; t546;];
tauJB = t1;
