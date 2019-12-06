% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR3
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:53
% EndTime: 2019-12-05 15:46:56
% DurationCPUTime: 2.34s
% Computational Cost: add. (26601->218), mult. (47684->279), div. (0->0), fcn. (29351->10), ass. (0->95)
t540 = sin(pkin(8));
t542 = cos(pkin(8));
t526 = -t542 * g(1) - t540 * g(2);
t538 = -g(3) + qJDD(1);
t545 = sin(qJ(2));
t548 = cos(qJ(2));
t509 = -t545 * t526 + t548 * t538;
t505 = qJDD(2) * pkin(2) + t509;
t510 = t548 * t526 + t545 * t538;
t549 = qJD(2) ^ 2;
t506 = -t549 * pkin(2) + t510;
t539 = sin(pkin(9));
t541 = cos(pkin(9));
t491 = t539 * t505 + t541 * t506;
t488 = -t549 * pkin(3) + qJDD(2) * pkin(6) + t491;
t525 = t540 * g(1) - t542 * g(2);
t524 = qJDD(3) - t525;
t544 = sin(qJ(4));
t547 = cos(qJ(4));
t484 = -t544 * t488 + t547 * t524;
t562 = qJD(2) * qJD(4);
t560 = t547 * t562;
t522 = t544 * qJDD(2) + t560;
t481 = (-t522 + t560) * pkin(7) + (t544 * t547 * t549 + qJDD(4)) * pkin(4) + t484;
t485 = t547 * t488 + t544 * t524;
t523 = t547 * qJDD(2) - t544 * t562;
t564 = qJD(2) * t544;
t529 = qJD(4) * pkin(4) - pkin(7) * t564;
t537 = t547 ^ 2;
t482 = -t537 * t549 * pkin(4) + t523 * pkin(7) - qJD(4) * t529 + t485;
t543 = sin(qJ(5));
t546 = cos(qJ(5));
t479 = t546 * t481 - t543 * t482;
t514 = (-t543 * t544 + t546 * t547) * qJD(2);
t496 = t514 * qJD(5) + t546 * t522 + t543 * t523;
t515 = (t543 * t547 + t544 * t546) * qJD(2);
t501 = -t514 * mrSges(6,1) + t515 * mrSges(6,2);
t536 = qJD(4) + qJD(5);
t507 = -t536 * mrSges(6,2) + t514 * mrSges(6,3);
t535 = qJDD(4) + qJDD(5);
t476 = m(6) * t479 + t535 * mrSges(6,1) - t496 * mrSges(6,3) - t515 * t501 + t536 * t507;
t480 = t543 * t481 + t546 * t482;
t495 = -t515 * qJD(5) - t543 * t522 + t546 * t523;
t508 = t536 * mrSges(6,1) - t515 * mrSges(6,3);
t477 = m(6) * t480 - t535 * mrSges(6,2) + t495 * mrSges(6,3) + t514 * t501 - t536 * t508;
t467 = t546 * t476 + t543 * t477;
t512 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t544 + Ifges(5,2) * t547) * qJD(2);
t513 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t544 + Ifges(5,4) * t547) * qJD(2);
t498 = Ifges(6,4) * t515 + Ifges(6,2) * t514 + Ifges(6,6) * t536;
t499 = Ifges(6,1) * t515 + Ifges(6,4) * t514 + Ifges(6,5) * t536;
t552 = -mrSges(6,1) * t479 + mrSges(6,2) * t480 - Ifges(6,5) * t496 - Ifges(6,6) * t495 - Ifges(6,3) * t535 - t515 * t498 + t514 * t499;
t568 = mrSges(5,1) * t484 - mrSges(5,2) * t485 + Ifges(5,5) * t522 + Ifges(5,6) * t523 + Ifges(5,3) * qJDD(4) + pkin(4) * t467 + (t544 * t512 - t547 * t513) * qJD(2) - t552;
t490 = t541 * t505 - t539 * t506;
t554 = -qJDD(2) * pkin(3) - t490;
t483 = t529 * t564 - t523 * pkin(4) + (-pkin(7) * t537 - pkin(6)) * t549 + t554;
t497 = Ifges(6,5) * t515 + Ifges(6,6) * t514 + Ifges(6,3) * t536;
t468 = -mrSges(6,1) * t483 + mrSges(6,3) * t480 + Ifges(6,4) * t496 + Ifges(6,2) * t495 + Ifges(6,6) * t535 - t515 * t497 + t536 * t499;
t469 = mrSges(6,2) * t483 - mrSges(6,3) * t479 + Ifges(6,1) * t496 + Ifges(6,4) * t495 + Ifges(6,5) * t535 + t514 * t497 - t536 * t498;
t487 = -t549 * pkin(6) + t554;
t511 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t544 + Ifges(5,6) * t547) * qJD(2);
t553 = m(6) * t483 - t495 * mrSges(6,1) + t496 * mrSges(6,2) - t514 * t507 + t515 * t508;
t556 = -t543 * t476 + t546 * t477;
t446 = -mrSges(5,1) * t487 + mrSges(5,3) * t485 + Ifges(5,4) * t522 + Ifges(5,2) * t523 + Ifges(5,6) * qJDD(4) - pkin(4) * t553 + pkin(7) * t556 + qJD(4) * t513 + t546 * t468 + t543 * t469 - t511 * t564;
t521 = (-mrSges(5,1) * t547 + mrSges(5,2) * t544) * qJD(2);
t563 = qJD(2) * t547;
t528 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t563;
t465 = m(5) * t484 + qJDD(4) * mrSges(5,1) - t522 * mrSges(5,3) + qJD(4) * t528 - t521 * t564 + t467;
t527 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t564;
t466 = m(5) * t485 - qJDD(4) * mrSges(5,2) + t523 * mrSges(5,3) - qJD(4) * t527 + t521 * t563 + t556;
t461 = -t544 * t465 + t547 * t466;
t456 = m(4) * t491 - t549 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t461;
t472 = -m(5) * t487 + t523 * mrSges(5,1) - t522 * mrSges(5,2) - t527 * t564 + t528 * t563 - t553;
t471 = m(4) * t490 + qJDD(2) * mrSges(4,1) - t549 * mrSges(4,2) + t472;
t452 = t539 * t456 + t541 * t471;
t453 = mrSges(5,2) * t487 - mrSges(5,3) * t484 + Ifges(5,1) * t522 + Ifges(5,4) * t523 + Ifges(5,5) * qJDD(4) - pkin(7) * t467 - qJD(4) * t512 - t543 * t468 + t546 * t469 + t511 * t563;
t567 = mrSges(3,1) * t509 + mrSges(4,1) * t490 - mrSges(3,2) * t510 - mrSges(4,2) * t491 + pkin(2) * t452 + pkin(3) * t472 + pkin(6) * t461 + t547 * t446 + t544 * t453 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t450 = m(3) * t509 + qJDD(2) * mrSges(3,1) - t549 * mrSges(3,2) + t452;
t557 = t541 * t456 - t539 * t471;
t451 = m(3) * t510 - t549 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t557;
t558 = -t545 * t450 + t548 * t451;
t443 = m(2) * t526 + t558;
t460 = t547 * t465 + t544 * t466;
t459 = m(4) * t524 + t460;
t458 = (m(2) + m(3)) * t525 - t459;
t565 = t540 * t443 + t542 * t458;
t444 = t548 * t450 + t545 * t451;
t561 = m(2) * t538 + t444;
t559 = t542 * t443 - t540 * t458;
t445 = -mrSges(4,1) * t524 + mrSges(4,3) * t491 + t549 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t460 - t568;
t440 = mrSges(4,2) * t524 - mrSges(4,3) * t490 + Ifges(4,5) * qJDD(2) - t549 * Ifges(4,6) - pkin(6) * t460 - t544 * t446 + t547 * t453;
t439 = -mrSges(3,2) * t525 - mrSges(3,3) * t509 + Ifges(3,5) * qJDD(2) - t549 * Ifges(3,6) - qJ(3) * t452 + t541 * t440 - t539 * t445;
t438 = -mrSges(2,1) * t538 + mrSges(2,3) * t526 - pkin(1) * t444 - t567;
t437 = mrSges(3,1) * t525 + mrSges(3,3) * t510 + t549 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t459 + qJ(3) * t557 + t539 * t440 + t541 * t445;
t436 = mrSges(2,2) * t538 - mrSges(2,3) * t525 - pkin(5) * t444 - t545 * t437 + t548 * t439;
t1 = [-m(1) * g(1) + t559; -m(1) * g(2) + t565; -m(1) * g(3) + t561; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t565 + t542 * t436 - t540 * t438; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t559 + t540 * t436 + t542 * t438; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t525 - mrSges(2,2) * t526 + t545 * t439 + t548 * t437 + pkin(1) * (m(3) * t525 - t459) + pkin(5) * t558; t561; t567; t459; t568; -t552;];
tauJB = t1;
