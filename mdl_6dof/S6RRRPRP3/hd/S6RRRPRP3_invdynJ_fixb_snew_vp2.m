% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:37:50
% EndTime: 2019-05-07 07:37:59
% DurationCPUTime: 5.98s
% Computational Cost: add. (65338->326), mult. (134881->403), div. (0->0), fcn. (97006->10), ass. (0->128)
t556 = Ifges(6,1) + Ifges(7,1);
t548 = Ifges(6,4) - Ifges(7,5);
t547 = -Ifges(6,5) - Ifges(7,4);
t555 = Ifges(6,2) + Ifges(7,3);
t546 = Ifges(6,6) - Ifges(7,6);
t554 = -Ifges(6,3) - Ifges(7,2);
t519 = sin(qJ(3));
t520 = sin(qJ(2));
t522 = cos(qJ(2));
t552 = cos(qJ(3));
t498 = (t519 * t522 + t520 * t552) * qJD(1);
t514 = qJD(2) + qJD(3);
t516 = sin(pkin(10));
t517 = cos(pkin(10));
t485 = -t498 * t516 + t514 * t517;
t486 = t498 * t517 + t514 * t516;
t518 = sin(qJ(5));
t551 = cos(qJ(5));
t453 = -t485 * t551 + t518 * t486;
t538 = qJD(1) * t522;
t539 = qJD(1) * t520;
t497 = t519 * t539 - t552 * t538;
t537 = qJD(1) * qJD(2);
t504 = qJDD(1) * t520 + t522 * t537;
t505 = qJDD(1) * t522 - t520 * t537;
t471 = -t497 * qJD(3) + t504 * t552 + t519 * t505;
t513 = qJDD(2) + qJDD(3);
t455 = -t471 * t516 + t513 * t517;
t456 = t471 * t517 + t513 * t516;
t420 = -t453 * qJD(5) + t518 * t455 + t456 * t551;
t454 = t518 * t485 + t486 * t551;
t435 = mrSges(7,1) * t453 - mrSges(7,3) * t454;
t470 = qJD(3) * t498 + t504 * t519 - t552 * t505;
t508 = qJD(2) * pkin(2) - pkin(8) * t539;
t515 = t522 ^ 2;
t524 = qJD(1) ^ 2;
t521 = sin(qJ(1));
t523 = cos(qJ(1));
t535 = t521 * g(1) - t523 * g(2);
t529 = -qJDD(1) * pkin(1) - t535;
t472 = -t505 * pkin(2) + t508 * t539 + (-pkin(8) * t515 - pkin(7)) * t524 + t529;
t421 = (t497 * t514 - t471) * qJ(4) + (t498 * t514 + t470) * pkin(3) + t472;
t530 = -g(1) * t523 - g(2) * t521;
t500 = -pkin(1) * t524 + qJDD(1) * pkin(7) + t530;
t544 = t520 * t500;
t550 = pkin(2) * t524;
t462 = qJDD(2) * pkin(2) - t504 * pkin(8) - t544 + (pkin(8) * t537 + t520 * t550 - g(3)) * t522;
t488 = -g(3) * t520 + t522 * t500;
t463 = pkin(8) * t505 - qJD(2) * t508 - t515 * t550 + t488;
t439 = t519 * t462 + t552 * t463;
t479 = pkin(3) * t497 - qJ(4) * t498;
t512 = t514 ^ 2;
t424 = -pkin(3) * t512 + qJ(4) * t513 - t479 * t497 + t439;
t411 = -0.2e1 * qJD(4) * t486 + t517 * t421 - t516 * t424;
t408 = (t485 * t497 - t456) * pkin(9) + (t485 * t486 + t470) * pkin(4) + t411;
t412 = 0.2e1 * qJD(4) * t485 + t516 * t421 + t517 * t424;
t475 = pkin(4) * t497 - pkin(9) * t486;
t484 = t485 ^ 2;
t410 = -pkin(4) * t484 + pkin(9) * t455 - t475 * t497 + t412;
t403 = t408 * t551 - t518 * t410;
t434 = pkin(5) * t453 - qJ(6) * t454;
t469 = qJDD(5) + t470;
t493 = qJD(5) + t497;
t492 = t493 ^ 2;
t402 = -t469 * pkin(5) - t492 * qJ(6) + t454 * t434 + qJDD(6) - t403;
t442 = -mrSges(7,2) * t453 + mrSges(7,3) * t493;
t531 = -m(7) * t402 + t469 * mrSges(7,1) + t493 * t442;
t398 = t420 * mrSges(7,2) + t454 * t435 - t531;
t404 = t518 * t408 + t551 * t410;
t401 = -pkin(5) * t492 + qJ(6) * t469 + 0.2e1 * qJD(6) * t493 - t434 * t453 + t404;
t419 = t454 * qJD(5) - t455 * t551 + t518 * t456;
t445 = -mrSges(7,1) * t493 + mrSges(7,2) * t454;
t536 = m(7) * t401 + t469 * mrSges(7,3) + t493 * t445;
t541 = -t548 * t453 + t556 * t454 - t547 * t493;
t543 = t555 * t453 - t548 * t454 - t546 * t493;
t553 = t546 * t419 + t547 * t420 - t541 * t453 + t543 * t454 + t554 * t469 - mrSges(6,1) * t403 + mrSges(7,1) * t402 + mrSges(6,2) * t404 - mrSges(7,3) * t401 + pkin(5) * t398 - qJ(6) * (-t419 * mrSges(7,2) - t453 * t435 + t536);
t549 = -mrSges(6,3) - mrSges(7,2);
t480 = mrSges(4,1) * t497 + mrSges(4,2) * t498;
t490 = mrSges(4,1) * t514 - mrSges(4,3) * t498;
t444 = mrSges(6,1) * t493 - mrSges(6,3) * t454;
t540 = -mrSges(6,1) * t453 - mrSges(6,2) * t454 - t435;
t390 = m(6) * t404 - t469 * mrSges(6,2) + t549 * t419 - t493 * t444 + t540 * t453 + t536;
t443 = -mrSges(6,2) * t493 - mrSges(6,3) * t453;
t392 = m(6) * t403 + t469 * mrSges(6,1) + t549 * t420 + t493 * t443 + t540 * t454 + t531;
t385 = t518 * t390 + t551 * t392;
t458 = -mrSges(5,1) * t485 + mrSges(5,2) * t486;
t473 = -mrSges(5,2) * t497 + mrSges(5,3) * t485;
t383 = m(5) * t411 + mrSges(5,1) * t470 - mrSges(5,3) * t456 - t458 * t486 + t473 * t497 + t385;
t474 = mrSges(5,1) * t497 - mrSges(5,3) * t486;
t532 = t551 * t390 - t392 * t518;
t384 = m(5) * t412 - mrSges(5,2) * t470 + mrSges(5,3) * t455 + t458 * t485 - t474 * t497 + t532;
t533 = -t383 * t516 + t517 * t384;
t377 = m(4) * t439 - mrSges(4,2) * t513 - mrSges(4,3) * t470 - t480 * t497 - t490 * t514 + t533;
t438 = t552 * t462 - t519 * t463;
t423 = -t513 * pkin(3) - t512 * qJ(4) + t498 * t479 + qJDD(4) - t438;
t413 = -t455 * pkin(4) - t484 * pkin(9) + t486 * t475 + t423;
t406 = -0.2e1 * qJD(6) * t454 + (t453 * t493 - t420) * qJ(6) + (t454 * t493 + t419) * pkin(5) + t413;
t399 = m(7) * t406 + t419 * mrSges(7,1) - t420 * mrSges(7,3) + t453 * t442 - t454 * t445;
t526 = m(6) * t413 + t419 * mrSges(6,1) + mrSges(6,2) * t420 + t453 * t443 + t444 * t454 + t399;
t396 = m(5) * t423 - t455 * mrSges(5,1) + mrSges(5,2) * t456 - t485 * t473 + t474 * t486 + t526;
t489 = -mrSges(4,2) * t514 - mrSges(4,3) * t497;
t394 = m(4) * t438 + mrSges(4,1) * t513 - mrSges(4,3) * t471 - t480 * t498 + t489 * t514 - t396;
t374 = t519 * t377 + t552 * t394;
t379 = t517 * t383 + t516 * t384;
t542 = t546 * t453 + t547 * t454 + t554 * t493;
t534 = t552 * t377 - t394 * t519;
t386 = -mrSges(6,1) * t413 - mrSges(7,1) * t406 + mrSges(7,2) * t401 + mrSges(6,3) * t404 - pkin(5) * t399 - t555 * t419 + t548 * t420 + t542 * t454 + t546 * t469 + t541 * t493;
t387 = mrSges(6,2) * t413 + mrSges(7,2) * t402 - mrSges(6,3) * t403 - mrSges(7,3) * t406 - qJ(6) * t399 - t548 * t419 + t556 * t420 + t542 * t453 - t547 * t469 + t543 * t493;
t446 = Ifges(5,5) * t486 + Ifges(5,6) * t485 + Ifges(5,3) * t497;
t448 = Ifges(5,1) * t486 + Ifges(5,4) * t485 + Ifges(5,5) * t497;
t371 = -mrSges(5,1) * t423 + mrSges(5,3) * t412 + Ifges(5,4) * t456 + Ifges(5,2) * t455 + Ifges(5,6) * t470 - pkin(4) * t526 + pkin(9) * t532 + t386 * t551 + t518 * t387 - t486 * t446 + t497 * t448;
t447 = Ifges(5,4) * t486 + Ifges(5,2) * t485 + Ifges(5,6) * t497;
t373 = mrSges(5,2) * t423 - mrSges(5,3) * t411 + Ifges(5,1) * t456 + Ifges(5,4) * t455 + Ifges(5,5) * t470 - pkin(9) * t385 - t518 * t386 + t387 * t551 + t485 * t446 - t497 * t447;
t477 = Ifges(4,4) * t498 - Ifges(4,2) * t497 + Ifges(4,6) * t514;
t478 = Ifges(4,1) * t498 - Ifges(4,4) * t497 + Ifges(4,5) * t514;
t528 = mrSges(4,1) * t438 - mrSges(4,2) * t439 + Ifges(4,5) * t471 - Ifges(4,6) * t470 + Ifges(4,3) * t513 - pkin(3) * t396 + qJ(4) * t533 + t517 * t371 + t516 * t373 + t498 * t477 + t497 * t478;
t527 = m(4) * t472 + t470 * mrSges(4,1) + t471 * mrSges(4,2) + t497 * t489 + t498 * t490 + t379;
t507 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t538;
t506 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t539;
t503 = (-mrSges(3,1) * t522 + mrSges(3,2) * t520) * qJD(1);
t499 = -t524 * pkin(7) + t529;
t496 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t520 + Ifges(3,4) * t522) * qJD(1);
t495 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t520 + Ifges(3,2) * t522) * qJD(1);
t487 = -t522 * g(3) - t544;
t476 = Ifges(4,5) * t498 - Ifges(4,6) * t497 + Ifges(4,3) * t514;
t369 = (-Ifges(4,2) - Ifges(5,3)) * t470 + Ifges(4,6) * t513 + t514 * t478 - t498 * t476 + t485 * t448 - t486 * t447 + Ifges(4,4) * t471 - mrSges(4,1) * t472 - Ifges(5,5) * t456 - Ifges(5,6) * t455 + mrSges(4,3) * t439 - mrSges(5,1) * t411 + mrSges(5,2) * t412 - pkin(4) * t385 - pkin(3) * t379 + t553;
t368 = mrSges(4,2) * t472 - mrSges(4,3) * t438 + Ifges(4,1) * t471 - Ifges(4,4) * t470 + Ifges(4,5) * t513 - qJ(4) * t379 - t371 * t516 + t373 * t517 - t476 * t497 - t477 * t514;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t535 - mrSges(2,2) * t530 + t520 * (mrSges(3,2) * t499 - mrSges(3,3) * t487 + Ifges(3,1) * t504 + Ifges(3,4) * t505 + Ifges(3,5) * qJDD(2) - pkin(8) * t374 - qJD(2) * t495 + t368 * t552 - t519 * t369) + t522 * (-mrSges(3,1) * t499 + mrSges(3,3) * t488 + Ifges(3,4) * t504 + Ifges(3,2) * t505 + Ifges(3,6) * qJDD(2) - pkin(2) * t527 + pkin(8) * t534 + qJD(2) * t496 + t519 * t368 + t369 * t552) + pkin(1) * (-m(3) * t499 + t505 * mrSges(3,1) - t504 * mrSges(3,2) + (-t506 * t520 + t507 * t522) * qJD(1) - t527) + pkin(7) * (t522 * (m(3) * t488 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t505 - qJD(2) * t506 + t503 * t538 + t534) - t520 * (m(3) * t487 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t504 + qJD(2) * t507 - t503 * t539 + t374)); Ifges(3,3) * qJDD(2) + Ifges(3,5) * t504 + Ifges(3,6) * t505 + mrSges(3,1) * t487 - mrSges(3,2) * t488 + (t520 * t495 - t522 * t496) * qJD(1) + pkin(2) * t374 + t528; t528; t396; -t553; t398;];
tauJ  = t1;
