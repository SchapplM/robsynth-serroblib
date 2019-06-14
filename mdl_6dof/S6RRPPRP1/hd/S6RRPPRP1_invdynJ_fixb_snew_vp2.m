% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-05-06 09:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:02:58
% EndTime: 2019-05-06 09:03:05
% DurationCPUTime: 4.96s
% Computational Cost: add. (47254->325), mult. (111492->403), div. (0->0), fcn. (78832->10), ass. (0->126)
t555 = -2 * qJD(3);
t554 = Ifges(6,1) + Ifges(7,1);
t547 = Ifges(6,4) - Ifges(7,5);
t546 = -Ifges(6,5) - Ifges(7,4);
t553 = Ifges(6,2) + Ifges(7,3);
t545 = Ifges(6,6) - Ifges(7,6);
t552 = -Ifges(6,3) - Ifges(7,2);
t517 = sin(qJ(2));
t519 = cos(qJ(2));
t534 = qJD(1) * qJD(2);
t505 = qJDD(1) * t517 + t519 * t534;
t522 = qJD(1) ^ 2;
t518 = sin(qJ(1));
t520 = cos(qJ(1));
t527 = -g(1) * t520 - g(2) * t518;
t502 = -pkin(1) * t522 + qJDD(1) * pkin(7) + t527;
t542 = t517 * t502;
t549 = pkin(2) * t522;
t458 = qJDD(2) * pkin(2) - t505 * qJ(3) - t542 + (qJ(3) * t534 + t517 * t549 - g(3)) * t519;
t487 = -g(3) * t517 + t519 * t502;
t506 = qJDD(1) * t519 - t517 * t534;
t537 = qJD(1) * t517;
t507 = qJD(2) * pkin(2) - qJ(3) * t537;
t512 = t519 ^ 2;
t460 = qJ(3) * t506 - qJD(2) * t507 - t512 * t549 + t487;
t514 = sin(pkin(9));
t543 = cos(pkin(9));
t496 = (t514 * t519 + t543 * t517) * qJD(1);
t437 = t543 * t458 - t514 * t460 + t496 * t555;
t513 = sin(pkin(10));
t515 = cos(pkin(10));
t484 = qJD(2) * t515 - t496 * t513;
t485 = qJD(2) * t513 + t496 * t515;
t516 = sin(qJ(5));
t550 = cos(qJ(5));
t454 = -t550 * t484 + t516 * t485;
t479 = t543 * t505 + t514 * t506;
t468 = qJDD(2) * t515 - t479 * t513;
t469 = qJDD(2) * t513 + t479 * t515;
t427 = -t454 * qJD(5) + t516 * t468 + t550 * t469;
t455 = t516 * t484 + t550 * t485;
t440 = mrSges(7,1) * t454 - mrSges(7,3) * t455;
t536 = qJD(1) * t519;
t495 = t514 * t537 - t543 * t536;
t438 = t514 * t458 + t543 * t460 + t495 * t555;
t473 = pkin(3) * t495 - qJ(4) * t496;
t521 = qJD(2) ^ 2;
t419 = -pkin(3) * t521 + qJDD(2) * qJ(4) - t473 * t495 + t438;
t532 = t518 * g(1) - t520 * g(2);
t525 = -qJDD(1) * pkin(1) - t532;
t462 = -t506 * pkin(2) + qJDD(3) + t507 * t537 + (-qJ(3) * t512 - pkin(7)) * t522 + t525;
t478 = t505 * t514 - t543 * t506;
t422 = (qJD(2) * t495 - t479) * qJ(4) + (qJD(2) * t496 + t478) * pkin(3) + t462;
t414 = -0.2e1 * qJD(4) * t485 - t513 * t419 + t515 * t422;
t411 = (t484 * t495 - t469) * pkin(8) + (t484 * t485 + t478) * pkin(4) + t414;
t415 = 0.2e1 * qJD(4) * t484 + t515 * t419 + t513 * t422;
t465 = pkin(4) * t495 - pkin(8) * t485;
t483 = t484 ^ 2;
t413 = -pkin(4) * t483 + pkin(8) * t468 - t465 * t495 + t415;
t406 = t550 * t411 - t516 * t413;
t439 = pkin(5) * t454 - qJ(6) * t455;
t477 = qJDD(5) + t478;
t494 = qJD(5) + t495;
t493 = t494 ^ 2;
t405 = -t477 * pkin(5) - t493 * qJ(6) + t455 * t439 + qJDD(6) - t406;
t444 = -mrSges(7,2) * t454 + mrSges(7,3) * t494;
t528 = -m(7) * t405 + t477 * mrSges(7,1) + t494 * t444;
t401 = t427 * mrSges(7,2) + t455 * t440 - t528;
t407 = t516 * t411 + t550 * t413;
t404 = -pkin(5) * t493 + qJ(6) * t477 + 0.2e1 * qJD(6) * t494 - t439 * t454 + t407;
t426 = t455 * qJD(5) - t550 * t468 + t516 * t469;
t447 = -mrSges(7,1) * t494 + mrSges(7,2) * t455;
t533 = m(7) * t404 + t477 * mrSges(7,3) + t494 * t447;
t539 = -t547 * t454 + t554 * t455 - t546 * t494;
t540 = t553 * t454 - t547 * t455 - t545 * t494;
t551 = t545 * t426 + t546 * t427 - t539 * t454 + t540 * t455 + t552 * t477 - mrSges(6,1) * t406 + mrSges(7,1) * t405 + mrSges(6,2) * t407 - mrSges(7,3) * t404 + pkin(5) * t401 - qJ(6) * (-t426 * mrSges(7,2) - t454 * t440 + t533);
t548 = -mrSges(6,3) - mrSges(7,2);
t474 = mrSges(4,1) * t495 + mrSges(4,2) * t496;
t489 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t496;
t446 = mrSges(6,1) * t494 - mrSges(6,3) * t455;
t538 = -mrSges(6,1) * t454 - mrSges(6,2) * t455 - t440;
t394 = m(6) * t407 - t477 * mrSges(6,2) + t548 * t426 - t494 * t446 + t538 * t454 + t533;
t445 = -mrSges(6,2) * t494 - mrSges(6,3) * t454;
t396 = m(6) * t406 + t477 * mrSges(6,1) + t548 * t427 + t494 * t445 + t538 * t455 + t528;
t389 = t516 * t394 + t550 * t396;
t459 = -mrSges(5,1) * t484 + mrSges(5,2) * t485;
t463 = -mrSges(5,2) * t495 + mrSges(5,3) * t484;
t387 = m(5) * t414 + mrSges(5,1) * t478 - mrSges(5,3) * t469 - t459 * t485 + t463 * t495 + t389;
t464 = mrSges(5,1) * t495 - mrSges(5,3) * t485;
t529 = t550 * t394 - t396 * t516;
t388 = m(5) * t415 - mrSges(5,2) * t478 + mrSges(5,3) * t468 + t459 * t484 - t464 * t495 + t529;
t530 = -t387 * t513 + t515 * t388;
t381 = m(4) * t438 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t478 - qJD(2) * t489 - t474 * t495 + t530;
t418 = -qJDD(2) * pkin(3) - t521 * qJ(4) + t496 * t473 + qJDD(4) - t437;
t416 = -t468 * pkin(4) - t483 * pkin(8) + t485 * t465 + t418;
t409 = -0.2e1 * qJD(6) * t455 + (t454 * t494 - t427) * qJ(6) + (t455 * t494 + t426) * pkin(5) + t416;
t402 = m(7) * t409 + t426 * mrSges(7,1) - t427 * mrSges(7,3) + t454 * t444 - t455 * t447;
t524 = m(6) * t416 + t426 * mrSges(6,1) + mrSges(6,2) * t427 + t454 * t445 + t446 * t455 + t402;
t399 = m(5) * t418 - t468 * mrSges(5,1) + mrSges(5,2) * t469 - t484 * t463 + t464 * t485 + t524;
t488 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t495;
t398 = m(4) * t437 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t479 + qJD(2) * t488 - t474 * t496 - t399;
t378 = t514 * t381 + t543 * t398;
t383 = t515 * t387 + t513 * t388;
t541 = t545 * t454 + t546 * t455 + t552 * t494;
t531 = t543 * t381 - t398 * t514;
t382 = m(4) * t462 + t478 * mrSges(4,1) + t479 * mrSges(4,2) + t495 * t488 + t496 * t489 + t383;
t509 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t536;
t508 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t537;
t504 = (-mrSges(3,1) * t519 + mrSges(3,2) * t517) * qJD(1);
t501 = -t522 * pkin(7) + t525;
t499 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t517 + Ifges(3,4) * t519) * qJD(1);
t498 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t517 + Ifges(3,2) * t519) * qJD(1);
t486 = -t519 * g(3) - t542;
t472 = Ifges(4,1) * t496 - Ifges(4,4) * t495 + Ifges(4,5) * qJD(2);
t471 = Ifges(4,4) * t496 - Ifges(4,2) * t495 + Ifges(4,6) * qJD(2);
t470 = Ifges(4,5) * t496 - Ifges(4,6) * t495 + Ifges(4,3) * qJD(2);
t450 = Ifges(5,1) * t485 + Ifges(5,4) * t484 + Ifges(5,5) * t495;
t449 = Ifges(5,4) * t485 + Ifges(5,2) * t484 + Ifges(5,6) * t495;
t448 = Ifges(5,5) * t485 + Ifges(5,6) * t484 + Ifges(5,3) * t495;
t391 = mrSges(6,2) * t416 + mrSges(7,2) * t405 - mrSges(6,3) * t406 - mrSges(7,3) * t409 - qJ(6) * t402 - t547 * t426 + t554 * t427 + t541 * t454 - t546 * t477 + t540 * t494;
t390 = -mrSges(6,1) * t416 - mrSges(7,1) * t409 + mrSges(7,2) * t404 + mrSges(6,3) * t407 - pkin(5) * t402 - t553 * t426 + t547 * t427 + t541 * t455 + t545 * t477 + t539 * t494;
t377 = mrSges(5,2) * t418 - mrSges(5,3) * t414 + Ifges(5,1) * t469 + Ifges(5,4) * t468 + Ifges(5,5) * t478 - pkin(8) * t389 - t516 * t390 + t550 * t391 + t484 * t448 - t495 * t449;
t376 = -mrSges(5,1) * t418 + mrSges(5,3) * t415 + Ifges(5,4) * t469 + Ifges(5,2) * t468 + Ifges(5,6) * t478 - pkin(4) * t524 + pkin(8) * t529 + t550 * t390 + t516 * t391 - t485 * t448 + t495 * t450;
t375 = t551 + Ifges(4,6) * qJDD(2) - t496 * t470 + t484 * t450 - t485 * t449 + Ifges(4,4) * t479 - Ifges(5,6) * t468 - Ifges(5,5) * t469 + qJD(2) * t472 - mrSges(4,1) * t462 + mrSges(4,3) * t438 - mrSges(5,1) * t414 + mrSges(5,2) * t415 + (-Ifges(5,3) - Ifges(4,2)) * t478 - pkin(4) * t389 - pkin(3) * t383;
t374 = mrSges(4,2) * t462 - mrSges(4,3) * t437 + Ifges(4,1) * t479 - Ifges(4,4) * t478 + Ifges(4,5) * qJDD(2) - qJ(4) * t383 - qJD(2) * t471 - t376 * t513 + t377 * t515 - t470 * t495;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t532 - mrSges(2,2) * t527 + t517 * (mrSges(3,2) * t501 - mrSges(3,3) * t486 + Ifges(3,1) * t505 + Ifges(3,4) * t506 + Ifges(3,5) * qJDD(2) - qJ(3) * t378 - qJD(2) * t498 + t543 * t374 - t514 * t375) + t519 * (-mrSges(3,1) * t501 + mrSges(3,3) * t487 + Ifges(3,4) * t505 + Ifges(3,2) * t506 + Ifges(3,6) * qJDD(2) - pkin(2) * t382 + qJ(3) * t531 + qJD(2) * t499 + t514 * t374 + t543 * t375) + pkin(1) * (-m(3) * t501 + t506 * mrSges(3,1) - t505 * mrSges(3,2) + (-t508 * t517 + t509 * t519) * qJD(1) - t382) + pkin(7) * (t519 * (m(3) * t487 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t506 - qJD(2) * t508 + t504 * t536 + t531) - t517 * (m(3) * t486 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t505 + qJD(2) * t509 - t504 * t537 + t378)); Ifges(3,5) * t505 + Ifges(3,6) * t506 + mrSges(3,1) * t486 - mrSges(3,2) * t487 + Ifges(4,5) * t479 - Ifges(4,6) * t478 + t496 * t471 + t495 * t472 + mrSges(4,1) * t437 - mrSges(4,2) * t438 + t513 * t377 + t515 * t376 - pkin(3) * t399 + qJ(4) * t530 + pkin(2) * t378 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t517 * t498 - t519 * t499) * qJD(1); t382; t399; -t551; t401;];
tauJ  = t1;
