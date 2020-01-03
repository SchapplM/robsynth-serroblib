% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:53
% EndTime: 2019-12-31 19:05:55
% DurationCPUTime: 1.92s
% Computational Cost: add. (27097->227), mult. (35035->276), div. (0->0), fcn. (16535->8), ass. (0->97)
t546 = qJD(1) ^ 2;
t541 = sin(qJ(1));
t545 = cos(qJ(1));
t520 = -t545 * g(1) - t541 * g(2);
t553 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t520;
t568 = -pkin(1) - pkin(2);
t499 = t568 * t546 + t553;
t519 = t541 * g(1) - t545 * g(2);
t552 = -t546 * qJ(2) + qJDD(2) - t519;
t502 = t568 * qJDD(1) + t552;
t540 = sin(qJ(3));
t544 = cos(qJ(3));
t488 = t544 * t499 + t540 * t502;
t530 = -qJD(1) + qJD(3);
t526 = t530 ^ 2;
t528 = -qJDD(1) + qJDD(3);
t485 = -(t526 * pkin(3)) + t528 * pkin(7) + t488;
t539 = sin(qJ(4));
t543 = cos(qJ(4));
t476 = t543 * g(3) - t539 * t485;
t561 = qJD(4) * t530;
t560 = t543 * t561;
t513 = t539 * t528 + t560;
t472 = (-t513 + t560) * pkin(8) + (t526 * t539 * t543 + qJDD(4)) * pkin(4) + t476;
t477 = t539 * g(3) + t543 * t485;
t514 = t543 * t528 - t539 * t561;
t564 = t530 * t539;
t517 = qJD(4) * pkin(4) - pkin(8) * t564;
t536 = t543 ^ 2;
t473 = -t536 * t526 * pkin(4) + t514 * pkin(8) - qJD(4) * t517 + t477;
t538 = sin(qJ(5));
t542 = cos(qJ(5));
t470 = t542 * t472 - t538 * t473;
t506 = (-t538 * t539 + t542 * t543) * t530;
t483 = t506 * qJD(5) + t542 * t513 + t538 * t514;
t507 = (t538 * t543 + t539 * t542) * t530;
t493 = -t506 * mrSges(6,1) + t507 * mrSges(6,2);
t529 = qJD(4) + qJD(5);
t494 = -t529 * mrSges(6,2) + t506 * mrSges(6,3);
t527 = qJDD(4) + qJDD(5);
t467 = m(6) * t470 + t527 * mrSges(6,1) - t483 * mrSges(6,3) - t507 * t493 + t529 * t494;
t471 = t538 * t472 + t542 * t473;
t482 = -t507 * qJD(5) - t538 * t513 + t542 * t514;
t495 = t529 * mrSges(6,1) - t507 * mrSges(6,3);
t468 = m(6) * t471 - t527 * mrSges(6,2) + t482 * mrSges(6,3) + t506 * t493 - t529 * t495;
t459 = t542 * t467 + t538 * t468;
t504 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t539 + Ifges(5,2) * t543) * t530;
t505 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t539 + Ifges(5,4) * t543) * t530;
t490 = Ifges(6,4) * t507 + Ifges(6,2) * t506 + Ifges(6,6) * t529;
t491 = Ifges(6,1) * t507 + Ifges(6,4) * t506 + Ifges(6,5) * t529;
t550 = -mrSges(6,1) * t470 + mrSges(6,2) * t471 - Ifges(6,5) * t483 - Ifges(6,6) * t482 - Ifges(6,3) * t527 - t507 * t490 + t506 * t491;
t570 = mrSges(5,1) * t476 - mrSges(5,2) * t477 + Ifges(5,5) * t513 + Ifges(5,6) * t514 + Ifges(5,3) * qJDD(4) + pkin(4) * t459 + (t539 * t504 - t543 * t505) * t530 - t550;
t569 = -m(3) - m(4);
t567 = -mrSges(2,1) - mrSges(3,1);
t566 = Ifges(3,4) + Ifges(2,5);
t565 = Ifges(2,6) - Ifges(3,6);
t563 = t530 * t543;
t508 = -t546 * pkin(1) + t553;
t512 = (-mrSges(5,1) * t543 + mrSges(5,2) * t539) * t530;
t516 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t563;
t457 = m(5) * t476 + qJDD(4) * mrSges(5,1) - t513 * mrSges(5,3) + qJD(4) * t516 - t512 * t564 + t459;
t515 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t564;
t557 = -t538 * t467 + t542 * t468;
t458 = m(5) * t477 - qJDD(4) * mrSges(5,2) + t514 * mrSges(5,3) - qJD(4) * t515 + t512 * t563 + t557;
t455 = -t539 * t457 + t543 * t458;
t452 = m(4) * t488 - (t526 * mrSges(4,1)) - t528 * mrSges(4,2) + t455;
t487 = -t540 * t499 + t544 * t502;
t554 = -t528 * pkin(3) - t487;
t484 = -t526 * pkin(7) + t554;
t474 = t517 * t564 - t514 * pkin(4) + (-pkin(8) * t536 - pkin(7)) * t526 + t554;
t551 = m(6) * t474 - t482 * mrSges(6,1) + t483 * mrSges(6,2) - t506 * t494 + t507 * t495;
t463 = -m(5) * t484 + t514 * mrSges(5,1) - t513 * mrSges(5,2) - t515 * t564 + t516 * t563 - t551;
t462 = m(4) * t487 + t528 * mrSges(4,1) - t526 * mrSges(4,2) + t463;
t558 = t544 * t452 - t540 * t462;
t555 = m(3) * t508 + qJDD(1) * mrSges(3,3) + t558;
t444 = m(2) * t520 - qJDD(1) * mrSges(2,2) + t567 * t546 + t555;
t449 = t540 * t452 + t544 * t462;
t511 = -qJDD(1) * pkin(1) + t552;
t448 = m(3) * t511 - qJDD(1) * mrSges(3,1) - t546 * mrSges(3,3) + t449;
t445 = m(2) * t519 + qJDD(1) * mrSges(2,1) - t546 * mrSges(2,2) - t448;
t562 = t541 * t444 + t545 * t445;
t559 = t545 * t444 - t541 * t445;
t454 = t543 * t457 + t539 * t458;
t489 = Ifges(6,5) * t507 + Ifges(6,6) * t506 + Ifges(6,3) * t529;
t460 = -mrSges(6,1) * t474 + mrSges(6,3) * t471 + Ifges(6,4) * t483 + Ifges(6,2) * t482 + Ifges(6,6) * t527 - t507 * t489 + t529 * t491;
t461 = mrSges(6,2) * t474 - mrSges(6,3) * t470 + Ifges(6,1) * t483 + Ifges(6,4) * t482 + Ifges(6,5) * t527 + t506 * t489 - t529 * t490;
t503 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t539 + Ifges(5,6) * t543) * t530;
t440 = -mrSges(5,1) * t484 + mrSges(5,3) * t477 + Ifges(5,4) * t513 + Ifges(5,2) * t514 + Ifges(5,6) * qJDD(4) - pkin(4) * t551 + pkin(8) * t557 + qJD(4) * t505 + t542 * t460 + t538 * t461 - t503 * t564;
t450 = mrSges(5,2) * t484 - mrSges(5,3) * t476 + Ifges(5,1) * t513 + Ifges(5,4) * t514 + Ifges(5,5) * qJDD(4) - pkin(8) * t459 - qJD(4) * t504 - t538 * t460 + t542 * t461 + t503 * t563;
t549 = mrSges(4,1) * t487 - mrSges(4,2) * t488 + Ifges(4,3) * t528 + pkin(3) * t463 + pkin(7) * t455 + t543 * t440 + t539 * t450;
t547 = -mrSges(3,1) * t511 - mrSges(2,2) * t520 - pkin(2) * t449 + qJ(2) * (-t546 * mrSges(3,1) + t555) - pkin(1) * t448 + mrSges(3,3) * t508 + mrSges(2,1) * t519 - t549 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);
t453 = t569 * g(3) - t454;
t439 = -mrSges(4,1) * g(3) + mrSges(4,3) * t488 + (t526 * Ifges(4,5)) + Ifges(4,6) * t528 - pkin(3) * t454 - t570;
t438 = mrSges(4,2) * g(3) - mrSges(4,3) * t487 + Ifges(4,5) * t528 - t526 * Ifges(4,6) - pkin(7) * t454 - t539 * t440 + t543 * t450;
t437 = mrSges(3,2) * t511 - mrSges(2,3) * t519 - pkin(6) * t449 - qJ(2) * t453 + t544 * t438 - t540 * t439 - t565 * t546 + t566 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t436 = mrSges(2,3) * t520 + mrSges(3,2) * t508 - t540 * t438 - t544 * t439 + pkin(2) * t454 - pkin(6) * t558 - pkin(1) * t453 + t566 * t546 + t565 * qJDD(1) + (pkin(2) * m(4) - t567) * g(3);
t1 = [-m(1) * g(1) + t559; -m(1) * g(2) + t562; (-m(1) - m(2) + t569) * g(3) - t454; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t562 - t541 * t436 + t545 * t437; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t559 + t545 * t436 + t541 * t437; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t547; t547; t448; t549; t570; -t550;];
tauJB = t1;
