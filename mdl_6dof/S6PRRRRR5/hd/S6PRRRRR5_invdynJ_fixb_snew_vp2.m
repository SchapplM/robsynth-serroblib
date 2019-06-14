% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:57:31
% EndTime: 2019-05-05 11:57:38
% DurationCPUTime: 7.19s
% Computational Cost: add. (87650->302), mult. (180255->396), div. (0->0), fcn. (146072->16), ass. (0->135)
t494 = sin(pkin(13));
t497 = cos(pkin(13));
t486 = g(1) * t494 - g(2) * t497;
t493 = -g(3) + qJDD(1);
t496 = sin(pkin(6));
t499 = cos(pkin(6));
t532 = t486 * t499 + t493 * t496;
t495 = sin(pkin(7));
t503 = sin(qJ(3));
t508 = cos(qJ(3));
t523 = qJD(2) * qJD(3);
t479 = (-qJDD(2) * t508 + t503 * t523) * t495;
t487 = -g(1) * t497 - g(2) * t494;
t504 = sin(qJ(2));
t509 = cos(qJ(2));
t453 = -t487 * t504 + t509 * t532;
t510 = qJD(2) ^ 2;
t531 = pkin(9) * t495;
t449 = qJDD(2) * pkin(2) + t510 * t531 + t453;
t454 = t509 * t487 + t504 * t532;
t450 = -pkin(2) * t510 + qJDD(2) * t531 + t454;
t498 = cos(pkin(7));
t471 = -t486 * t496 + t493 * t499;
t529 = t471 * t495;
t417 = -t503 * t450 + (t449 * t498 + t529) * t508;
t492 = qJD(2) * t498 + qJD(3);
t524 = qJD(2) * t495;
t521 = t508 * t524;
t475 = -mrSges(4,2) * t492 + mrSges(4,3) * t521;
t476 = (-mrSges(4,1) * t508 + mrSges(4,2) * t503) * t524;
t478 = (qJDD(2) * t503 + t508 * t523) * t495;
t491 = qJDD(2) * t498 + qJDD(3);
t526 = t498 * t503;
t418 = t449 * t526 + t450 * t508 + t503 * t529;
t477 = (-pkin(3) * t508 - pkin(10) * t503) * t524;
t490 = t492 ^ 2;
t413 = -pkin(3) * t490 + pkin(10) * t491 + t477 * t521 + t418;
t465 = t498 * t471;
t416 = t479 * pkin(3) - t478 * pkin(10) + t465 + (-t449 + (pkin(3) * t503 - pkin(10) * t508) * t492 * qJD(2)) * t495;
t502 = sin(qJ(4));
t507 = cos(qJ(4));
t404 = t413 * t507 + t416 * t502;
t522 = t503 * t524;
t469 = t492 * t507 - t502 * t522;
t470 = t492 * t502 + t507 * t522;
t452 = -pkin(4) * t469 - pkin(11) * t470;
t473 = qJDD(4) + t479;
t485 = qJD(4) - t521;
t483 = t485 ^ 2;
t394 = -pkin(4) * t483 + pkin(11) * t473 + t452 * t469 + t404;
t412 = -t491 * pkin(3) - t490 * pkin(10) + t477 * t522 - t417;
t444 = -qJD(4) * t470 - t478 * t502 + t491 * t507;
t445 = qJD(4) * t469 + t478 * t507 + t491 * t502;
t397 = (-t469 * t485 - t445) * pkin(11) + (t470 * t485 - t444) * pkin(4) + t412;
t501 = sin(qJ(5));
t506 = cos(qJ(5));
t389 = -t501 * t394 + t397 * t506;
t456 = -t470 * t501 + t485 * t506;
t421 = qJD(5) * t456 + t445 * t506 + t473 * t501;
t442 = qJDD(5) - t444;
t457 = t470 * t506 + t485 * t501;
t468 = qJD(5) - t469;
t387 = (t456 * t468 - t421) * pkin(12) + (t456 * t457 + t442) * pkin(5) + t389;
t390 = t394 * t506 + t397 * t501;
t420 = -qJD(5) * t457 - t445 * t501 + t473 * t506;
t435 = pkin(5) * t468 - pkin(12) * t457;
t455 = t456 ^ 2;
t388 = -pkin(5) * t455 + pkin(12) * t420 - t435 * t468 + t390;
t500 = sin(qJ(6));
t505 = cos(qJ(6));
t385 = t387 * t505 - t388 * t500;
t428 = t456 * t505 - t457 * t500;
t402 = qJD(6) * t428 + t420 * t500 + t421 * t505;
t429 = t456 * t500 + t457 * t505;
t414 = -mrSges(7,1) * t428 + mrSges(7,2) * t429;
t466 = qJD(6) + t468;
t422 = -mrSges(7,2) * t466 + mrSges(7,3) * t428;
t437 = qJDD(6) + t442;
t381 = m(7) * t385 + mrSges(7,1) * t437 - mrSges(7,3) * t402 - t414 * t429 + t422 * t466;
t386 = t387 * t500 + t388 * t505;
t401 = -qJD(6) * t429 + t420 * t505 - t421 * t500;
t423 = mrSges(7,1) * t466 - mrSges(7,3) * t429;
t382 = m(7) * t386 - mrSges(7,2) * t437 + mrSges(7,3) * t401 + t414 * t428 - t423 * t466;
t373 = t381 * t505 + t382 * t500;
t430 = -mrSges(6,1) * t456 + mrSges(6,2) * t457;
t433 = -mrSges(6,2) * t468 + mrSges(6,3) * t456;
t371 = m(6) * t389 + mrSges(6,1) * t442 - mrSges(6,3) * t421 - t430 * t457 + t433 * t468 + t373;
t434 = mrSges(6,1) * t468 - mrSges(6,3) * t457;
t518 = -t381 * t500 + t382 * t505;
t372 = m(6) * t390 - mrSges(6,2) * t442 + mrSges(6,3) * t420 + t430 * t456 - t434 * t468 + t518;
t368 = t371 * t506 + t372 * t501;
t458 = -mrSges(5,2) * t485 + mrSges(5,3) * t469;
t459 = mrSges(5,1) * t485 - mrSges(5,3) * t470;
t513 = -m(5) * t412 + mrSges(5,1) * t444 - mrSges(5,2) * t445 + t458 * t469 - t459 * t470 - t368;
t364 = m(4) * t417 + mrSges(4,1) * t491 - mrSges(4,3) * t478 + t475 * t492 - t476 * t522 + t513;
t530 = t364 * t508;
t474 = mrSges(4,1) * t492 - mrSges(4,3) * t522;
t369 = -t371 * t501 + t372 * t506;
t451 = -mrSges(5,1) * t469 + mrSges(5,2) * t470;
t367 = m(5) * t404 - mrSges(5,2) * t473 + mrSges(5,3) * t444 + t451 * t469 - t459 * t485 + t369;
t403 = -t413 * t502 + t416 * t507;
t393 = -pkin(4) * t473 - pkin(11) * t483 + t452 * t470 - t403;
t391 = -pkin(5) * t420 - pkin(12) * t455 + t435 * t457 + t393;
t515 = m(7) * t391 - mrSges(7,1) * t401 + mrSges(7,2) * t402 - t422 * t428 + t423 * t429;
t383 = -m(6) * t393 + mrSges(6,1) * t420 - mrSges(6,2) * t421 + t433 * t456 - t434 * t457 - t515;
t377 = m(5) * t403 + mrSges(5,1) * t473 - mrSges(5,3) * t445 - t451 * t470 + t458 * t485 + t383;
t519 = t367 * t507 - t377 * t502;
t358 = m(4) * t418 - mrSges(4,2) * t491 - mrSges(4,3) * t479 - t474 * t492 + t476 * t521 + t519;
t525 = t358 * t526 + t498 * t530;
t360 = t367 * t502 + t377 * t507;
t520 = t358 * t508 - t364 * t503;
t407 = Ifges(7,4) * t429 + Ifges(7,2) * t428 + Ifges(7,6) * t466;
t408 = Ifges(7,1) * t429 + Ifges(7,4) * t428 + Ifges(7,5) * t466;
t514 = -mrSges(7,1) * t385 + mrSges(7,2) * t386 - Ifges(7,5) * t402 - Ifges(7,6) * t401 - Ifges(7,3) * t437 - t407 * t429 + t408 * t428;
t406 = Ifges(7,5) * t429 + Ifges(7,6) * t428 + Ifges(7,3) * t466;
t374 = -mrSges(7,1) * t391 + mrSges(7,3) * t386 + Ifges(7,4) * t402 + Ifges(7,2) * t401 + Ifges(7,6) * t437 - t406 * t429 + t408 * t466;
t375 = mrSges(7,2) * t391 - mrSges(7,3) * t385 + Ifges(7,1) * t402 + Ifges(7,4) * t401 + Ifges(7,5) * t437 + t406 * t428 - t407 * t466;
t424 = Ifges(6,5) * t457 + Ifges(6,6) * t456 + Ifges(6,3) * t468;
t426 = Ifges(6,1) * t457 + Ifges(6,4) * t456 + Ifges(6,5) * t468;
t361 = -mrSges(6,1) * t393 + mrSges(6,3) * t390 + Ifges(6,4) * t421 + Ifges(6,2) * t420 + Ifges(6,6) * t442 - pkin(5) * t515 + pkin(12) * t518 + t505 * t374 + t500 * t375 - t457 * t424 + t468 * t426;
t425 = Ifges(6,4) * t457 + Ifges(6,2) * t456 + Ifges(6,6) * t468;
t362 = mrSges(6,2) * t393 - mrSges(6,3) * t389 + Ifges(6,1) * t421 + Ifges(6,4) * t420 + Ifges(6,5) * t442 - pkin(12) * t373 - t374 * t500 + t375 * t505 + t424 * t456 - t425 * t468;
t439 = Ifges(5,4) * t470 + Ifges(5,2) * t469 + Ifges(5,6) * t485;
t440 = Ifges(5,1) * t470 + Ifges(5,4) * t469 + Ifges(5,5) * t485;
t512 = mrSges(5,1) * t403 - mrSges(5,2) * t404 + Ifges(5,5) * t445 + Ifges(5,6) * t444 + Ifges(5,3) * t473 + pkin(4) * t383 + pkin(11) * t369 + t506 * t361 + t501 * t362 + t470 * t439 - t469 * t440;
t511 = mrSges(6,1) * t389 - mrSges(6,2) * t390 + Ifges(6,5) * t421 + Ifges(6,6) * t420 + Ifges(6,3) * t442 + pkin(5) * t373 + t425 * t457 - t426 * t456 - t514;
t463 = Ifges(4,5) * t492 + (Ifges(4,1) * t503 + Ifges(4,4) * t508) * t524;
t462 = Ifges(4,6) * t492 + (Ifges(4,4) * t503 + Ifges(4,2) * t508) * t524;
t438 = Ifges(5,5) * t470 + Ifges(5,6) * t469 + Ifges(5,3) * t485;
t431 = -t495 * t449 + t465;
t359 = m(4) * t431 + t479 * mrSges(4,1) + t478 * mrSges(4,2) + (t474 * t503 - t475 * t508) * t524 + t360;
t355 = -mrSges(5,1) * t412 + mrSges(5,3) * t404 + Ifges(5,4) * t445 + Ifges(5,2) * t444 + Ifges(5,6) * t473 - pkin(4) * t368 - t438 * t470 + t440 * t485 - t511;
t354 = mrSges(5,2) * t412 - mrSges(5,3) * t403 + Ifges(5,1) * t445 + Ifges(5,4) * t444 + Ifges(5,5) * t473 - pkin(11) * t368 - t361 * t501 + t362 * t506 + t438 * t469 - t439 * t485;
t353 = Ifges(4,5) * t478 - Ifges(4,6) * t479 + Ifges(4,3) * t491 + mrSges(4,1) * t417 - mrSges(4,2) * t418 + t502 * t354 + t507 * t355 + pkin(3) * t513 + pkin(10) * t519 + (t462 * t503 - t463 * t508) * t524;
t1 = [m(2) * t493 + t499 * (m(3) * t471 + t498 * t359 + (t358 * t503 + t530) * t495) + (t504 * (m(3) * t454 - mrSges(3,1) * t510 - qJDD(2) * mrSges(3,2) + t520) + t509 * (m(3) * t453 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t510 - t359 * t495 + t525)) * t496; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t453 - mrSges(3,2) * t454 + t498 * t353 + pkin(2) * t525 + (t503 * (mrSges(4,2) * t431 - mrSges(4,3) * t417 + Ifges(4,1) * t478 - Ifges(4,4) * t479 + Ifges(4,5) * t491 - pkin(10) * t360 + t354 * t507 - t355 * t502 - t462 * t492) + t508 * (-mrSges(4,1) * t431 + mrSges(4,3) * t418 + Ifges(4,4) * t478 - Ifges(4,2) * t479 + Ifges(4,6) * t491 - pkin(3) * t360 + t492 * t463 - t512) - pkin(2) * t359 + pkin(9) * t520) * t495; t353; t512; t511; -t514;];
tauJ  = t1;
