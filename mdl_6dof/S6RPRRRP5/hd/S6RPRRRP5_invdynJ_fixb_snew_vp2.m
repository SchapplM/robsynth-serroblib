% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:28:56
% EndTime: 2019-05-06 01:29:02
% DurationCPUTime: 5.00s
% Computational Cost: add. (48625->295), mult. (116686->363), div. (0->0), fcn. (91009->10), ass. (0->126)
t556 = Ifges(6,1) + Ifges(7,1);
t549 = Ifges(6,4) - Ifges(7,5);
t548 = -Ifges(6,5) - Ifges(7,4);
t555 = Ifges(6,2) + Ifges(7,3);
t547 = Ifges(6,6) - Ifges(7,6);
t554 = -Ifges(6,3) - Ifges(7,2);
t517 = qJD(1) ^ 2;
t508 = sin(pkin(10));
t509 = cos(pkin(10));
t512 = sin(qJ(3));
t515 = cos(qJ(3));
t524 = -t508 * t512 + t509 * t515;
t492 = t524 * qJD(1);
t525 = t508 * t515 + t509 * t512;
t493 = t525 * qJD(1);
t511 = sin(qJ(4));
t514 = cos(qJ(4));
t477 = t492 * t514 - t493 * t511;
t482 = -t493 * qJD(3) + qJDD(1) * t524;
t537 = t492 * qJD(3);
t483 = qJDD(1) * t525 + t537;
t447 = qJD(4) * t477 + t482 * t511 + t483 * t514;
t478 = t492 * t511 + t493 * t514;
t507 = qJD(3) + qJD(4);
t510 = sin(qJ(5));
t552 = cos(qJ(5));
t463 = t510 * t478 - t552 * t507;
t504 = qJDD(3) + qJDD(4);
t422 = -t463 * qJD(5) + t552 * t447 + t510 * t504;
t464 = t552 * t478 + t510 * t507;
t449 = mrSges(7,1) * t463 - mrSges(7,3) * t464;
t513 = sin(qJ(1));
t516 = cos(qJ(1));
t528 = -g(1) * t516 - g(2) * t513;
t494 = -pkin(1) * t517 + qJDD(1) * qJ(2) + t528;
t536 = qJD(1) * qJD(2);
t533 = -t509 * g(3) - 0.2e1 * t508 * t536;
t545 = pkin(7) * qJDD(1);
t551 = pkin(2) * t517;
t470 = (t509 * t551 - t494 - t545) * t508 + t533;
t485 = -g(3) * t508 + (t494 + 0.2e1 * t536) * t509;
t506 = t509 ^ 2;
t471 = -t506 * t551 + t509 * t545 + t485;
t451 = t515 * t470 - t512 * t471;
t418 = (-t483 + t537) * pkin(8) + (t492 * t493 + qJDD(3)) * pkin(3) + t451;
t452 = t512 * t470 + t515 * t471;
t488 = qJD(3) * pkin(3) - pkin(8) * t493;
t491 = t492 ^ 2;
t425 = -pkin(3) * t491 + pkin(8) * t482 - qJD(3) * t488 + t452;
t416 = t511 * t418 + t514 * t425;
t462 = -pkin(4) * t477 - pkin(9) * t478;
t503 = t507 ^ 2;
t411 = -pkin(4) * t503 + pkin(9) * t504 + t462 * t477 + t416;
t534 = t513 * g(1) - t516 * g(2);
t527 = qJDD(2) - t534;
t539 = -t508 ^ 2 - t506;
t481 = (-pkin(2) * t509 - pkin(1)) * qJDD(1) + (pkin(7) * t539 - qJ(2)) * t517 + t527;
t438 = -t482 * pkin(3) - t491 * pkin(8) + t493 * t488 + t481;
t446 = -qJD(4) * t478 + t482 * t514 - t483 * t511;
t413 = (-t477 * t507 - t447) * pkin(9) + (t478 * t507 - t446) * pkin(4) + t438;
t407 = -t510 * t411 + t552 * t413;
t445 = qJDD(5) - t446;
t448 = pkin(5) * t463 - qJ(6) * t464;
t473 = qJD(5) - t477;
t472 = t473 ^ 2;
t405 = -t445 * pkin(5) - t472 * qJ(6) + t464 * t448 + qJDD(6) - t407;
t453 = -mrSges(7,2) * t463 + mrSges(7,3) * t473;
t529 = -m(7) * t405 + t445 * mrSges(7,1) + t473 * t453;
t401 = t422 * mrSges(7,2) + t464 * t449 - t529;
t408 = t552 * t411 + t510 * t413;
t404 = -pkin(5) * t472 + qJ(6) * t445 + 0.2e1 * qJD(6) * t473 - t448 * t463 + t408;
t421 = t464 * qJD(5) + t510 * t447 - t552 * t504;
t456 = -mrSges(7,1) * t473 + mrSges(7,2) * t464;
t535 = m(7) * t404 + t445 * mrSges(7,3) + t473 * t456;
t541 = t549 * t463 - t556 * t464 + t548 * t473;
t542 = t555 * t463 - t549 * t464 - t547 * t473;
t553 = -t547 * t421 - t548 * t422 - t554 * t445 - t541 * t463 - t542 * t464 + mrSges(6,1) * t407 - mrSges(7,1) * t405 - mrSges(6,2) * t408 + mrSges(7,3) * t404 - pkin(5) * t401 + qJ(6) * (-t421 * mrSges(7,2) - t463 * t449 + t535);
t550 = -mrSges(6,3) - mrSges(7,2);
t461 = -mrSges(5,1) * t477 + mrSges(5,2) * t478;
t468 = mrSges(5,1) * t507 - mrSges(5,3) * t478;
t455 = mrSges(6,1) * t473 - mrSges(6,3) * t464;
t540 = -mrSges(6,1) * t463 - mrSges(6,2) * t464 - t449;
t396 = m(6) * t408 - t445 * mrSges(6,2) + t550 * t421 - t473 * t455 + t540 * t463 + t535;
t454 = -mrSges(6,2) * t473 - mrSges(6,3) * t463;
t398 = m(6) * t407 + t445 * mrSges(6,1) + t550 * t422 + t473 * t454 + t540 * t464 + t529;
t530 = t552 * t396 - t398 * t510;
t385 = m(5) * t416 - mrSges(5,2) * t504 + mrSges(5,3) * t446 + t461 * t477 - t468 * t507 + t530;
t415 = t514 * t418 - t511 * t425;
t467 = -mrSges(5,2) * t507 + mrSges(5,3) * t477;
t410 = -t504 * pkin(4) - t503 * pkin(9) + t478 * t462 - t415;
t406 = -0.2e1 * qJD(6) * t464 + (t463 * t473 - t422) * qJ(6) + (t464 * t473 + t421) * pkin(5) + t410;
t402 = m(7) * t406 + mrSges(7,1) * t421 - t422 * mrSges(7,3) + t453 * t463 - t464 * t456;
t519 = -m(6) * t410 - t421 * mrSges(6,1) - mrSges(6,2) * t422 - t463 * t454 - t455 * t464 - t402;
t393 = m(5) * t415 + mrSges(5,1) * t504 - mrSges(5,3) * t447 - t461 * t478 + t467 * t507 + t519;
t381 = t511 * t385 + t514 * t393;
t480 = -mrSges(4,1) * t492 + mrSges(4,2) * t493;
t486 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t492;
t379 = m(4) * t451 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t483 + qJD(3) * t486 - t480 * t493 + t381;
t487 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t493;
t531 = t514 * t385 - t393 * t511;
t380 = m(4) * t452 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t482 - qJD(3) * t487 + t480 * t492 + t531;
t544 = t515 * t379 + t512 * t380;
t391 = t510 * t396 + t552 * t398;
t543 = t547 * t463 + t548 * t464 + t554 * t473;
t532 = -t512 * t379 + t515 * t380;
t526 = -mrSges(3,1) * t509 + mrSges(3,2) * t508;
t523 = mrSges(3,3) * qJDD(1) + t517 * t526;
t522 = m(5) * t438 - t446 * mrSges(5,1) + t447 * mrSges(5,2) - t477 * t467 + t478 * t468 + t391;
t387 = -mrSges(6,1) * t410 - mrSges(7,1) * t406 + mrSges(7,2) * t404 + mrSges(6,3) * t408 - pkin(5) * t402 - t555 * t421 + t549 * t422 + t547 * t445 + t543 * t464 - t541 * t473;
t389 = mrSges(6,2) * t410 + mrSges(7,2) * t405 - mrSges(6,3) * t407 - mrSges(7,3) * t406 - qJ(6) * t402 - t549 * t421 + t556 * t422 - t548 * t445 + t543 * t463 + t542 * t473;
t458 = Ifges(5,4) * t478 + Ifges(5,2) * t477 + Ifges(5,6) * t507;
t459 = Ifges(5,1) * t478 + Ifges(5,4) * t477 + Ifges(5,5) * t507;
t521 = mrSges(5,1) * t415 - mrSges(5,2) * t416 + Ifges(5,5) * t447 + Ifges(5,6) * t446 + Ifges(5,3) * t504 + pkin(4) * t519 + pkin(9) * t530 + t552 * t387 + t510 * t389 + t478 * t458 - t459 * t477;
t518 = m(4) * t481 - t482 * mrSges(4,1) + t483 * mrSges(4,2) - t492 * t486 + t493 * t487 + t522;
t490 = -qJDD(1) * pkin(1) - t517 * qJ(2) + t527;
t484 = -t508 * t494 + t533;
t476 = Ifges(4,1) * t493 + Ifges(4,4) * t492 + Ifges(4,5) * qJD(3);
t475 = Ifges(4,4) * t493 + Ifges(4,2) * t492 + Ifges(4,6) * qJD(3);
t474 = Ifges(4,5) * t493 + Ifges(4,6) * t492 + Ifges(4,3) * qJD(3);
t457 = Ifges(5,5) * t478 + Ifges(5,6) * t477 + Ifges(5,3) * t507;
t382 = mrSges(3,3) * t517 * t539 + m(3) * t490 + qJDD(1) * t526 + t518;
t375 = -mrSges(5,1) * t438 + mrSges(5,3) * t416 + Ifges(5,4) * t447 + Ifges(5,2) * t446 + Ifges(5,6) * t504 - pkin(4) * t391 - t478 * t457 + t507 * t459 - t553;
t374 = mrSges(5,2) * t438 - mrSges(5,3) * t415 + Ifges(5,1) * t447 + Ifges(5,4) * t446 + Ifges(5,5) * t504 - pkin(9) * t391 - t510 * t387 + t552 * t389 + t477 * t457 - t507 * t458;
t373 = mrSges(4,2) * t481 - mrSges(4,3) * t451 + Ifges(4,1) * t483 + Ifges(4,4) * t482 + Ifges(4,5) * qJDD(3) - pkin(8) * t381 - qJD(3) * t475 + t374 * t514 - t375 * t511 + t474 * t492;
t372 = -mrSges(4,1) * t481 + mrSges(4,3) * t452 + Ifges(4,4) * t483 + Ifges(4,2) * t482 + Ifges(4,6) * qJDD(3) - pkin(3) * t522 + pkin(8) * t531 + qJD(3) * t476 + t511 * t374 + t514 * t375 - t493 * t474;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t534 - mrSges(2,2) * t528 + t508 * (mrSges(3,2) * t490 - mrSges(3,3) * t484 + t515 * t373 - t512 * t372 - pkin(7) * t544 + (Ifges(3,1) * t508 + Ifges(3,4) * t509) * qJDD(1)) + t509 * (-mrSges(3,1) * t490 + mrSges(3,3) * t485 + t512 * t373 + t515 * t372 - pkin(2) * t518 + pkin(7) * t532 + (Ifges(3,4) * t508 + Ifges(3,2) * t509) * qJDD(1)) - pkin(1) * t382 + qJ(2) * ((m(3) * t485 + t509 * t523 + t532) * t509 + (-m(3) * t484 + t508 * t523 - t544) * t508); t382; mrSges(4,1) * t451 - mrSges(4,2) * t452 + Ifges(4,5) * t483 + Ifges(4,6) * t482 + Ifges(4,3) * qJDD(3) + pkin(3) * t381 + t475 * t493 - t476 * t492 + t521; t521; t553; t401;];
tauJ  = t1;
