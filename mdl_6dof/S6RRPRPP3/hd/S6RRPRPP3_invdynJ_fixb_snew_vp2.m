% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 12:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:29:24
% EndTime: 2019-05-06 12:29:30
% DurationCPUTime: 3.07s
% Computational Cost: add. (22405->300), mult. (48097->353), div. (0->0), fcn. (32761->8), ass. (0->118)
t536 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t555 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t534 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t554 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t535 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t553 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t508 = sin(pkin(9));
t509 = cos(pkin(9));
t511 = sin(qJ(2));
t538 = t511 * qJD(1);
t493 = qJD(2) * t509 - t508 * t538;
t494 = qJD(2) * t508 + t509 * t538;
t510 = sin(qJ(4));
t544 = cos(qJ(4));
t466 = -t493 * t544 + t494 * t510;
t513 = cos(qJ(2));
t537 = qJD(1) * qJD(2);
t528 = t513 * t537;
t498 = qJDD(1) * t511 + t528;
t475 = qJDD(2) * t509 - t498 * t508;
t476 = qJDD(2) * t508 + t498 * t509;
t423 = -t466 * qJD(4) + t510 * t475 + t476 * t544;
t467 = t510 * t493 + t494 * t544;
t438 = -mrSges(7,2) * t467 + mrSges(7,3) * t466;
t516 = qJD(1) ^ 2;
t512 = sin(qJ(1));
t514 = cos(qJ(1));
t527 = g(1) * t512 - t514 * g(2);
t487 = -qJDD(1) * pkin(1) - pkin(7) * t516 - t527;
t505 = t511 * t537;
t499 = qJDD(1) * t513 - t505;
t444 = (-t498 - t528) * qJ(3) + (-t499 + t505) * pkin(2) + t487;
t523 = -g(1) * t514 - g(2) * t512;
t488 = -pkin(1) * t516 + qJDD(1) * pkin(7) + t523;
t471 = -g(3) * t511 + t513 * t488;
t496 = (-t513 * pkin(2) - t511 * qJ(3)) * qJD(1);
t515 = qJD(2) ^ 2;
t539 = qJD(1) * t513;
t450 = -pkin(2) * t515 + qJDD(2) * qJ(3) + t496 * t539 + t471;
t411 = -0.2e1 * qJD(3) * t494 + t509 * t444 - t450 * t508;
t407 = (-t493 * t539 - t476) * pkin(8) + (t493 * t494 - t499) * pkin(3) + t411;
t412 = 0.2e1 * qJD(3) * t493 + t508 * t444 + t509 * t450;
t477 = -pkin(3) * t539 - pkin(8) * t494;
t492 = t493 ^ 2;
t410 = -pkin(3) * t492 + pkin(8) * t475 + t477 * t539 + t412;
t404 = t407 * t544 - t510 * t410;
t439 = pkin(4) * t466 - qJ(5) * t467;
t495 = qJDD(4) - t499;
t503 = -qJD(4) + t539;
t502 = t503 ^ 2;
t401 = -t495 * pkin(4) - t502 * qJ(5) + t467 * t439 + qJDD(5) - t404;
t541 = t466 * t503;
t545 = 2 * qJD(6);
t395 = t503 * t545 + (t466 * t467 - t495) * qJ(6) + (t423 - t541) * pkin(5) + t401;
t456 = -mrSges(7,1) * t466 - mrSges(7,2) * t503;
t524 = -m(7) * t395 + t495 * mrSges(7,3) - t503 * t456;
t393 = mrSges(7,1) * t423 + t438 * t467 - t524;
t455 = mrSges(6,1) * t466 + mrSges(6,3) * t503;
t441 = -mrSges(6,2) * t466 - mrSges(6,3) * t467;
t548 = -m(6) * t401 - t423 * mrSges(6,1) - t467 * t441;
t390 = mrSges(6,2) * t495 - t455 * t503 + t393 - t548;
t422 = qJD(4) * t467 - t475 * t544 + t476 * t510;
t453 = pkin(5) * t467 + qJ(6) * t503;
t465 = t466 ^ 2;
t405 = t510 * t407 + t544 * t410;
t520 = -pkin(4) * t502 + qJ(5) * t495 - t439 * t466 + t405;
t546 = -2 * qJD(5);
t397 = -pkin(5) * t422 - qJ(6) * t465 + qJDD(6) + (t546 - t453) * t503 + t520;
t400 = 0.2e1 * qJD(5) * t503 - t520;
t457 = mrSges(6,1) * t467 - mrSges(6,2) * t503;
t454 = mrSges(7,1) * t467 + mrSges(7,3) * t503;
t532 = m(7) * t397 + t495 * mrSges(7,2) - t503 * t454;
t521 = -m(6) * t400 + t495 * mrSges(6,3) - t503 * t457 + t532;
t529 = -t536 * t466 + t467 * t554 - t535 * t503;
t530 = t555 * t466 + t536 * t467 - t534 * t503;
t540 = -t438 - t441;
t552 = -t534 * t422 + t529 * t466 + t553 * t495 + mrSges(5,1) * t404 - mrSges(5,2) * t405 + mrSges(6,2) * t401 + mrSges(7,2) * t397 - mrSges(6,3) * t400 - mrSges(7,3) * t395 - pkin(4) * t390 + qJ(5) * (t540 * t466 + (-mrSges(6,1) - mrSges(7,1)) * t422 + t521) - qJ(6) * t393 + t467 * t530 + t423 * t535;
t551 = pkin(1) * qJD(1);
t542 = -mrSges(7,1) - mrSges(5,3);
t440 = mrSges(5,1) * t466 + mrSges(5,2) * t467;
t451 = mrSges(5,2) * t503 - mrSges(5,3) * t466;
t385 = m(5) * t404 + (-t451 + t455) * t503 + (mrSges(5,1) - mrSges(6,2)) * t495 + (-t438 - t440) * t467 + t542 * t423 + t524 + t548;
t452 = -mrSges(5,1) * t503 - mrSges(5,3) * t467;
t388 = m(5) * t405 - mrSges(5,2) * t495 + t452 * t503 + (-t440 + t540) * t466 + (-mrSges(6,1) + t542) * t422 + t521;
t383 = t544 * t385 + t510 * t388;
t470 = -t513 * g(3) - t511 * t488;
t531 = t466 * t534 - t467 * t535 + t503 * t553;
t468 = -mrSges(4,1) * t493 + mrSges(4,2) * t494;
t473 = mrSges(4,2) * t539 + mrSges(4,3) * t493;
t381 = m(4) * t411 - mrSges(4,1) * t499 - mrSges(4,3) * t476 - t468 * t494 - t473 * t539 + t383;
t474 = -mrSges(4,1) * t539 - mrSges(4,3) * t494;
t525 = -t385 * t510 + t544 * t388;
t382 = m(4) * t412 + mrSges(4,2) * t499 + mrSges(4,3) * t475 + t468 * t493 + t474 * t539 + t525;
t526 = -t381 * t508 + t509 * t382;
t449 = -qJDD(2) * pkin(2) - qJ(3) * t515 + t496 * t538 + qJDD(3) - t470;
t413 = -pkin(3) * t475 - pkin(8) * t492 + t494 * t477 + t449;
t517 = (-t423 - t541) * qJ(5) + t413 + (-t503 * pkin(4) + t546) * t467;
t399 = -pkin(5) * t465 + t466 * t545 - t453 * t467 + (pkin(4) + qJ(6)) * t422 + t517;
t522 = m(7) * t399 - t423 * mrSges(7,2) + t422 * mrSges(7,3) - t467 * t454 + t466 * t456;
t377 = t381 * t509 + t382 * t508;
t403 = pkin(4) * t422 + t517;
t392 = m(6) * t403 - t422 * mrSges(6,2) - t423 * mrSges(6,3) - t466 * t455 - t467 * t457 + t522;
t519 = m(5) * t413 + t422 * mrSges(5,1) + t423 * mrSges(5,2) + t466 * t451 + t467 * t452 + t392;
t389 = m(4) * t449 - t475 * mrSges(4,1) + t476 * mrSges(4,2) - t493 * t473 + t494 * t474 + t519;
t501 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t539;
t500 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t538;
t497 = (-t513 * mrSges(3,1) + t511 * mrSges(3,2)) * qJD(1);
t486 = Ifges(3,5) * qJD(2) + (t511 * Ifges(3,1) + t513 * Ifges(3,4)) * qJD(1);
t485 = Ifges(3,6) * qJD(2) + (t511 * Ifges(3,4) + Ifges(3,2) * t513) * qJD(1);
t462 = Ifges(4,1) * t494 + Ifges(4,4) * t493 - Ifges(4,5) * t539;
t461 = Ifges(4,4) * t494 + Ifges(4,2) * t493 - Ifges(4,6) * t539;
t460 = Ifges(4,5) * t494 + Ifges(4,6) * t493 - Ifges(4,3) * t539;
t394 = -mrSges(7,1) * t422 - t438 * t466 + t532;
t379 = mrSges(6,1) * t401 + mrSges(7,1) * t395 + mrSges(5,2) * t413 - mrSges(7,2) * t399 - mrSges(5,3) * t404 - mrSges(6,3) * t403 + pkin(5) * t393 - qJ(5) * t392 - t536 * t422 + t423 * t554 + t531 * t466 + t535 * t495 + t530 * t503;
t378 = -mrSges(5,1) * t413 - mrSges(6,1) * t400 + mrSges(7,1) * t397 + mrSges(6,2) * t403 + mrSges(5,3) * t405 - mrSges(7,3) * t399 - pkin(4) * t392 + pkin(5) * t394 - qJ(6) * t522 + t555 * t422 + t536 * t423 + t531 * t467 + t534 * t495 - t529 * t503;
t376 = mrSges(4,2) * t449 - mrSges(4,3) * t411 + Ifges(4,1) * t476 + Ifges(4,4) * t475 - Ifges(4,5) * t499 - pkin(8) * t383 - t510 * t378 + t379 * t544 + t493 * t460 + t461 * t539;
t375 = -mrSges(4,1) * t449 + mrSges(4,3) * t412 + Ifges(4,4) * t476 + Ifges(4,2) * t475 - Ifges(4,6) * t499 - pkin(3) * t519 + pkin(8) * t525 + t378 * t544 + t510 * t379 - t494 * t460 - t462 * t539;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t527 - mrSges(2,2) * t523 + pkin(1) * (-m(3) * t487 + mrSges(3,1) * t499 - mrSges(3,2) * t498 - t377) + (mrSges(3,2) * t487 - mrSges(3,3) * t470 + Ifges(3,1) * t498 + Ifges(3,4) * t499 + Ifges(3,5) * qJDD(2) - qJ(3) * t377 - qJD(2) * t485 - t508 * t375 + t509 * t376 - t500 * t551 - pkin(7) * (m(3) * t470 + qJDD(2) * mrSges(3,1) - t498 * mrSges(3,3) + qJD(2) * t501 - t497 * t538 - t389)) * t511 + ((mrSges(3,3) * pkin(7) + Ifges(3,2) + Ifges(4,3)) * t499 + t501 * t551 + pkin(7) * (m(3) * t471 - qJDD(2) * mrSges(3,2) - qJD(2) * t500 + t497 * t539 + t526) + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t498 + t493 * t462 - t494 * t461 + qJD(2) * t486 - mrSges(3,1) * t487 - Ifges(4,6) * t475 - Ifges(4,5) * t476 + mrSges(3,3) * t471 - mrSges(4,1) * t411 + mrSges(4,2) * t412 - pkin(3) * t383 - pkin(2) * t377 - t552) * t513; Ifges(3,5) * t498 + Ifges(3,6) * t499 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t470 - mrSges(3,2) * t471 + t508 * t376 + t509 * t375 - pkin(2) * t389 + qJ(3) * t526 + (t511 * t485 - t513 * t486) * qJD(1); t389; t552; t390; t394;];
tauJ  = t1;
