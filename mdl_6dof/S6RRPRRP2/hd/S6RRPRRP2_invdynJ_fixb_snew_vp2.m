% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:23:16
% EndTime: 2019-05-06 17:23:26
% DurationCPUTime: 5.68s
% Computational Cost: add. (55276->325), mult. (127985->404), div. (0->0), fcn. (94231->10), ass. (0->128)
t552 = Ifges(6,1) + Ifges(7,1);
t545 = Ifges(6,4) - Ifges(7,5);
t544 = -Ifges(6,5) - Ifges(7,4);
t551 = Ifges(6,2) + Ifges(7,3);
t543 = Ifges(6,6) - Ifges(7,6);
t550 = -Ifges(6,3) - Ifges(7,2);
t512 = sin(pkin(10));
t513 = cos(pkin(10));
t516 = sin(qJ(2));
t519 = cos(qJ(2));
t490 = (-t512 * t516 + t513 * t519) * qJD(1);
t491 = (t512 * t519 + t513 * t516) * qJD(1);
t515 = sin(qJ(4));
t518 = cos(qJ(4));
t475 = t490 * t518 - t491 * t515;
t534 = qJD(1) * qJD(2);
t499 = qJDD(1) * t516 + t519 * t534;
t500 = qJDD(1) * t519 - t516 * t534;
t479 = -t499 * t512 + t500 * t513;
t480 = t499 * t513 + t500 * t512;
t446 = qJD(4) * t475 + t479 * t515 + t480 * t518;
t476 = t490 * t515 + t491 * t518;
t510 = qJD(2) + qJD(4);
t514 = sin(qJ(5));
t548 = cos(qJ(5));
t466 = t514 * t476 - t510 * t548;
t509 = qJDD(2) + qJDD(4);
t422 = -t466 * qJD(5) + t446 * t548 + t514 * t509;
t467 = t476 * t548 + t514 * t510;
t448 = mrSges(7,1) * t466 - mrSges(7,3) * t467;
t521 = qJD(1) ^ 2;
t517 = sin(qJ(1));
t520 = cos(qJ(1));
t527 = -g(1) * t520 - g(2) * t517;
t496 = -pkin(1) * t521 + qJDD(1) * pkin(7) + t527;
t541 = t516 * t496;
t547 = pkin(2) * t521;
t463 = qJDD(2) * pkin(2) - t499 * qJ(3) - t541 + (qJ(3) * t534 + t516 * t547 - g(3)) * t519;
t482 = -g(3) * t516 + t519 * t496;
t536 = qJD(1) * t516;
t501 = qJD(2) * pkin(2) - qJ(3) * t536;
t511 = t519 ^ 2;
t464 = qJ(3) * t500 - qJD(2) * t501 - t511 * t547 + t482;
t435 = -0.2e1 * qJD(3) * t491 + t513 * t463 - t512 * t464;
t415 = (qJD(2) * t490 - t480) * pkin(8) + (t490 * t491 + qJDD(2)) * pkin(3) + t435;
t436 = 0.2e1 * qJD(3) * t490 + t512 * t463 + t513 * t464;
t485 = qJD(2) * pkin(3) - pkin(8) * t491;
t489 = t490 ^ 2;
t418 = -pkin(3) * t489 + pkin(8) * t479 - qJD(2) * t485 + t436;
t413 = t515 * t415 + t518 * t418;
t459 = -pkin(4) * t475 - pkin(9) * t476;
t508 = t510 ^ 2;
t408 = -pkin(4) * t508 + pkin(9) * t509 + t459 * t475 + t413;
t532 = t517 * g(1) - t520 * g(2);
t526 = -qJDD(1) * pkin(1) - t532;
t465 = -t500 * pkin(2) + qJDD(3) + t501 * t536 + (-qJ(3) * t511 - pkin(7)) * t521 + t526;
t432 = -t479 * pkin(3) - t489 * pkin(8) + t491 * t485 + t465;
t445 = -qJD(4) * t476 + t479 * t518 - t480 * t515;
t410 = (-t475 * t510 - t446) * pkin(9) + (t476 * t510 - t445) * pkin(4) + t432;
t404 = -t514 * t408 + t548 * t410;
t444 = qJDD(5) - t445;
t447 = pkin(5) * t466 - qJ(6) * t467;
t471 = qJD(5) - t475;
t470 = t471 ^ 2;
t402 = -t444 * pkin(5) - t470 * qJ(6) + t467 * t447 + qJDD(6) - t404;
t450 = -mrSges(7,2) * t466 + mrSges(7,3) * t471;
t528 = -m(7) * t402 + t444 * mrSges(7,1) + t471 * t450;
t398 = t422 * mrSges(7,2) + t467 * t448 - t528;
t405 = t548 * t408 + t514 * t410;
t401 = -pkin(5) * t470 + qJ(6) * t444 + 0.2e1 * qJD(6) * t471 - t447 * t466 + t405;
t421 = t467 * qJD(5) + t514 * t446 - t509 * t548;
t453 = -mrSges(7,1) * t471 + mrSges(7,2) * t467;
t533 = m(7) * t401 + t444 * mrSges(7,3) + t471 * t453;
t538 = t545 * t466 - t552 * t467 + t544 * t471;
t539 = t551 * t466 - t545 * t467 - t543 * t471;
t549 = -t543 * t421 - t544 * t422 - t550 * t444 - t538 * t466 - t539 * t467 + mrSges(6,1) * t404 - mrSges(7,1) * t402 - mrSges(6,2) * t405 + mrSges(7,3) * t401 - pkin(5) * t398 + qJ(6) * (-t421 * mrSges(7,2) - t466 * t448 + t533);
t546 = -mrSges(6,3) - mrSges(7,2);
t458 = -mrSges(5,1) * t475 + mrSges(5,2) * t476;
t469 = mrSges(5,1) * t510 - mrSges(5,3) * t476;
t452 = mrSges(6,1) * t471 - mrSges(6,3) * t467;
t537 = -mrSges(6,1) * t466 - mrSges(6,2) * t467 - t448;
t393 = m(6) * t405 - t444 * mrSges(6,2) + t546 * t421 - t471 * t452 + t537 * t466 + t533;
t451 = -mrSges(6,2) * t471 - mrSges(6,3) * t466;
t395 = m(6) * t404 + t444 * mrSges(6,1) + t546 * t422 + t471 * t451 + t537 * t467 + t528;
t529 = t548 * t393 - t395 * t514;
t381 = m(5) * t413 - mrSges(5,2) * t509 + mrSges(5,3) * t445 + t458 * t475 - t469 * t510 + t529;
t412 = t518 * t415 - t515 * t418;
t468 = -mrSges(5,2) * t510 + mrSges(5,3) * t475;
t407 = -t509 * pkin(4) - t508 * pkin(9) + t476 * t459 - t412;
t403 = -0.2e1 * qJD(6) * t467 + (t466 * t471 - t422) * qJ(6) + (t467 * t471 + t421) * pkin(5) + t407;
t399 = m(7) * t403 + mrSges(7,1) * t421 - t422 * mrSges(7,3) + t450 * t466 - t467 * t453;
t522 = -m(6) * t407 - t421 * mrSges(6,1) - mrSges(6,2) * t422 - t466 * t451 - t452 * t467 - t399;
t390 = m(5) * t412 + mrSges(5,1) * t509 - mrSges(5,3) * t446 - t458 * t476 + t468 * t510 + t522;
t378 = t515 * t381 + t518 * t390;
t478 = -mrSges(4,1) * t490 + mrSges(4,2) * t491;
t483 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t490;
t376 = m(4) * t435 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t480 + qJD(2) * t483 - t478 * t491 + t378;
t484 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t491;
t530 = t518 * t381 - t390 * t515;
t377 = m(4) * t436 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t479 - qJD(2) * t484 + t478 * t490 + t530;
t370 = t513 * t376 + t512 * t377;
t388 = t514 * t393 + t548 * t395;
t540 = t543 * t466 + t544 * t467 + t550 * t471;
t535 = qJD(1) * t519;
t531 = -t376 * t512 + t513 * t377;
t525 = -m(5) * t432 + t445 * mrSges(5,1) - t446 * mrSges(5,2) + t475 * t468 - t476 * t469 - t388;
t384 = -mrSges(6,1) * t407 - mrSges(7,1) * t403 + mrSges(7,2) * t401 + mrSges(6,3) * t405 - pkin(5) * t399 - t551 * t421 + t545 * t422 + t543 * t444 + t540 * t467 - t538 * t471;
t386 = mrSges(6,2) * t407 + mrSges(7,2) * t402 - mrSges(6,3) * t404 - mrSges(7,3) * t403 - qJ(6) * t399 - t545 * t421 + t552 * t422 - t544 * t444 + t540 * t466 + t539 * t471;
t455 = Ifges(5,4) * t476 + Ifges(5,2) * t475 + Ifges(5,6) * t510;
t456 = Ifges(5,1) * t476 + Ifges(5,4) * t475 + Ifges(5,5) * t510;
t524 = mrSges(5,1) * t412 - mrSges(5,2) * t413 + Ifges(5,5) * t446 + Ifges(5,6) * t445 + Ifges(5,3) * t509 + pkin(4) * t522 + pkin(9) * t529 + t548 * t384 + t514 * t386 + t476 * t455 - t475 * t456;
t382 = m(4) * t465 - t479 * mrSges(4,1) + t480 * mrSges(4,2) - t490 * t483 + t491 * t484 - t525;
t503 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t535;
t502 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t536;
t498 = (-mrSges(3,1) * t519 + mrSges(3,2) * t516) * qJD(1);
t495 = -t521 * pkin(7) + t526;
t494 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t516 + Ifges(3,4) * t519) * qJD(1);
t493 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t516 + Ifges(3,2) * t519) * qJD(1);
t481 = -t519 * g(3) - t541;
t474 = Ifges(4,1) * t491 + Ifges(4,4) * t490 + Ifges(4,5) * qJD(2);
t473 = Ifges(4,4) * t491 + Ifges(4,2) * t490 + Ifges(4,6) * qJD(2);
t472 = Ifges(4,5) * t491 + Ifges(4,6) * t490 + Ifges(4,3) * qJD(2);
t454 = Ifges(5,5) * t476 + Ifges(5,6) * t475 + Ifges(5,3) * t510;
t372 = -mrSges(5,1) * t432 + mrSges(5,3) * t413 + Ifges(5,4) * t446 + Ifges(5,2) * t445 + Ifges(5,6) * t509 - pkin(4) * t388 - t476 * t454 + t510 * t456 - t549;
t371 = mrSges(5,2) * t432 - mrSges(5,3) * t412 + Ifges(5,1) * t446 + Ifges(5,4) * t445 + Ifges(5,5) * t509 - pkin(9) * t388 - t514 * t384 + t386 * t548 + t475 * t454 - t510 * t455;
t369 = mrSges(4,2) * t465 - mrSges(4,3) * t435 + Ifges(4,1) * t480 + Ifges(4,4) * t479 + Ifges(4,5) * qJDD(2) - pkin(8) * t378 - qJD(2) * t473 + t371 * t518 - t372 * t515 + t472 * t490;
t368 = -mrSges(4,1) * t465 + mrSges(4,3) * t436 + Ifges(4,4) * t480 + Ifges(4,2) * t479 + Ifges(4,6) * qJDD(2) + pkin(3) * t525 + pkin(8) * t530 + qJD(2) * t474 + t515 * t371 + t518 * t372 - t491 * t472;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t532 - mrSges(2,2) * t527 + t516 * (mrSges(3,2) * t495 - mrSges(3,3) * t481 + Ifges(3,1) * t499 + Ifges(3,4) * t500 + Ifges(3,5) * qJDD(2) - qJ(3) * t370 - qJD(2) * t493 - t512 * t368 + t513 * t369) + t519 * (-mrSges(3,1) * t495 + mrSges(3,3) * t482 + Ifges(3,4) * t499 + Ifges(3,2) * t500 + Ifges(3,6) * qJDD(2) - pkin(2) * t382 + qJ(3) * t531 + qJD(2) * t494 + t513 * t368 + t512 * t369) + pkin(1) * (-m(3) * t495 - t499 * mrSges(3,2) + t500 * mrSges(3,1) - t382 + (-t502 * t516 + t503 * t519) * qJD(1)) + pkin(7) * (t519 * (m(3) * t482 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t500 - qJD(2) * t502 + t498 * t535 + t531) - t516 * (m(3) * t481 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t499 + qJD(2) * t503 - t498 * t536 + t370)); t491 * t473 + Ifges(3,5) * t499 + Ifges(3,6) * t500 - t490 * t474 + Ifges(4,6) * t479 + Ifges(4,5) * t480 + mrSges(3,1) * t481 - mrSges(3,2) * t482 + mrSges(4,1) * t435 - mrSges(4,2) * t436 + pkin(3) * t378 + t524 + (t516 * t493 - t519 * t494) * qJD(1) + pkin(2) * t370 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); t382; t524; t549; t398;];
tauJ  = t1;
