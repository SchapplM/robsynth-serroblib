% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 08:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:43:59
% EndTime: 2019-05-05 08:44:06
% DurationCPUTime: 6.93s
% Computational Cost: add. (82184->300), mult. (172435->397), div. (0->0), fcn. (138973->16), ass. (0->132)
t486 = sin(pkin(12));
t490 = cos(pkin(12));
t477 = t486 * g(1) - t490 * g(2);
t484 = -g(3) + qJDD(1);
t488 = sin(pkin(6));
t492 = cos(pkin(6));
t523 = t477 * t492 + t484 * t488;
t487 = sin(pkin(7));
t495 = sin(qJ(3));
t498 = cos(qJ(3));
t513 = qJD(2) * qJD(3);
t469 = (-qJDD(2) * t498 + t495 * t513) * t487;
t478 = -t490 * g(1) - t486 * g(2);
t496 = sin(qJ(2));
t499 = cos(qJ(2));
t447 = -t496 * t478 + t523 * t499;
t500 = qJD(2) ^ 2;
t521 = t487 * pkin(9);
t440 = qJDD(2) * pkin(2) + t500 * t521 + t447;
t448 = t499 * t478 + t523 * t496;
t441 = -t500 * pkin(2) + qJDD(2) * t521 + t448;
t491 = cos(pkin(7));
t462 = -t488 * t477 + t492 * t484;
t519 = t462 * t487;
t411 = -t495 * t441 + t498 * (t440 * t491 + t519);
t522 = cos(qJ(4));
t483 = t491 * qJD(2) + qJD(3);
t514 = qJD(2) * t487;
t511 = t498 * t514;
t465 = -t483 * mrSges(4,2) + mrSges(4,3) * t511;
t466 = (-t498 * mrSges(4,1) + t495 * mrSges(4,2)) * t514;
t468 = (qJDD(2) * t495 + t498 * t513) * t487;
t482 = t491 * qJDD(2) + qJDD(3);
t516 = t491 * t495;
t412 = t440 * t516 + t498 * t441 + t495 * t519;
t467 = (-t498 * pkin(3) - t495 * pkin(10)) * t514;
t481 = t483 ^ 2;
t407 = -t481 * pkin(3) + t482 * pkin(10) + t467 * t511 + t412;
t458 = t491 * t462;
t410 = t469 * pkin(3) - t468 * pkin(10) + t458 + (-t440 + (pkin(3) * t495 - pkin(10) * t498) * t483 * qJD(2)) * t487;
t494 = sin(qJ(4));
t396 = t407 * t522 + t494 * t410;
t512 = t495 * t514;
t460 = -t483 * t522 + t494 * t512;
t461 = t494 * t483 + t512 * t522;
t442 = t460 * pkin(4) - t461 * qJ(5);
t463 = qJDD(4) + t469;
t476 = qJD(4) - t511;
t475 = t476 ^ 2;
t391 = -t475 * pkin(4) + t463 * qJ(5) - t460 * t442 + t396;
t406 = -t482 * pkin(3) - t481 * pkin(10) + t467 * t512 - t411;
t435 = t461 * qJD(4) + t494 * t468 - t482 * t522;
t436 = -t460 * qJD(4) + t468 * t522 + t494 * t482;
t394 = (t460 * t476 - t436) * qJ(5) + (t461 * t476 + t435) * pkin(4) + t406;
t485 = sin(pkin(13));
t489 = cos(pkin(13));
t450 = t489 * t461 + t485 * t476;
t386 = -0.2e1 * qJD(5) * t450 - t485 * t391 + t489 * t394;
t423 = t489 * t436 + t485 * t463;
t449 = -t485 * t461 + t489 * t476;
t384 = (t460 * t449 - t423) * pkin(11) + (t449 * t450 + t435) * pkin(5) + t386;
t387 = 0.2e1 * qJD(5) * t449 + t489 * t391 + t485 * t394;
t422 = -t485 * t436 + t489 * t463;
t428 = t460 * pkin(5) - t450 * pkin(11);
t446 = t449 ^ 2;
t385 = -t446 * pkin(5) + t422 * pkin(11) - t460 * t428 + t387;
t493 = sin(qJ(6));
t497 = cos(qJ(6));
t382 = t497 * t384 - t493 * t385;
t419 = t497 * t449 - t493 * t450;
t399 = t419 * qJD(6) + t493 * t422 + t497 * t423;
t420 = t493 * t449 + t497 * t450;
t408 = -t419 * mrSges(7,1) + t420 * mrSges(7,2);
t459 = qJD(6) + t460;
t413 = -t459 * mrSges(7,2) + t419 * mrSges(7,3);
t433 = qJDD(6) + t435;
t379 = m(7) * t382 + t433 * mrSges(7,1) - t399 * mrSges(7,3) - t420 * t408 + t459 * t413;
t383 = t493 * t384 + t497 * t385;
t398 = -t420 * qJD(6) + t497 * t422 - t493 * t423;
t414 = t459 * mrSges(7,1) - t420 * mrSges(7,3);
t380 = m(7) * t383 - t433 * mrSges(7,2) + t398 * mrSges(7,3) + t419 * t408 - t459 * t414;
t371 = t497 * t379 + t493 * t380;
t424 = -t449 * mrSges(6,1) + t450 * mrSges(6,2);
t507 = -t460 * mrSges(6,2) + t449 * mrSges(6,3);
t369 = m(6) * t386 + t435 * mrSges(6,1) - t423 * mrSges(6,3) - t450 * t424 + t460 * t507 + t371;
t427 = t460 * mrSges(6,1) - t450 * mrSges(6,3);
t508 = -t493 * t379 + t497 * t380;
t370 = m(6) * t387 - t435 * mrSges(6,2) + t422 * mrSges(6,3) + t449 * t424 - t460 * t427 + t508;
t366 = t489 * t369 + t485 * t370;
t451 = -t476 * mrSges(5,2) - t460 * mrSges(5,3);
t452 = t476 * mrSges(5,1) - t461 * mrSges(5,3);
t503 = -m(5) * t406 - t435 * mrSges(5,1) - t436 * mrSges(5,2) - t460 * t451 - t461 * t452 - t366;
t362 = m(4) * t411 + t482 * mrSges(4,1) - t468 * mrSges(4,3) + t483 * t465 - t466 * t512 + t503;
t520 = t362 * t498;
t464 = t483 * mrSges(4,1) - mrSges(4,3) * t512;
t367 = -t485 * t369 + t489 * t370;
t443 = t460 * mrSges(5,1) + t461 * mrSges(5,2);
t365 = m(5) * t396 - t463 * mrSges(5,2) - t435 * mrSges(5,3) - t460 * t443 - t476 * t452 + t367;
t395 = -t494 * t407 + t410 * t522;
t390 = -t463 * pkin(4) - t475 * qJ(5) + t461 * t442 + qJDD(5) - t395;
t388 = -t422 * pkin(5) - t446 * pkin(11) + t450 * t428 + t390;
t504 = m(7) * t388 - t398 * mrSges(7,1) + t399 * mrSges(7,2) - t419 * t413 + t420 * t414;
t381 = m(6) * t390 - t422 * mrSges(6,1) + t423 * mrSges(6,2) + t450 * t427 - t449 * t507 + t504;
t375 = m(5) * t395 + t463 * mrSges(5,1) - t436 * mrSges(5,3) - t461 * t443 + t476 * t451 - t381;
t509 = t365 * t522 - t494 * t375;
t356 = m(4) * t412 - t482 * mrSges(4,2) - t469 * mrSges(4,3) - t483 * t464 + t466 * t511 + t509;
t515 = t356 * t516 + t491 * t520;
t358 = t494 * t365 + t375 * t522;
t510 = t498 * t356 - t495 * t362;
t401 = Ifges(7,4) * t420 + Ifges(7,2) * t419 + Ifges(7,6) * t459;
t402 = Ifges(7,1) * t420 + Ifges(7,4) * t419 + Ifges(7,5) * t459;
t502 = mrSges(7,1) * t382 - mrSges(7,2) * t383 + Ifges(7,5) * t399 + Ifges(7,6) * t398 + Ifges(7,3) * t433 + t420 * t401 - t419 * t402;
t400 = Ifges(7,5) * t420 + Ifges(7,6) * t419 + Ifges(7,3) * t459;
t372 = -mrSges(7,1) * t388 + mrSges(7,3) * t383 + Ifges(7,4) * t399 + Ifges(7,2) * t398 + Ifges(7,6) * t433 - t420 * t400 + t459 * t402;
t373 = mrSges(7,2) * t388 - mrSges(7,3) * t382 + Ifges(7,1) * t399 + Ifges(7,4) * t398 + Ifges(7,5) * t433 + t419 * t400 - t459 * t401;
t415 = Ifges(6,5) * t450 + Ifges(6,6) * t449 + Ifges(6,3) * t460;
t417 = Ifges(6,1) * t450 + Ifges(6,4) * t449 + Ifges(6,5) * t460;
t359 = -mrSges(6,1) * t390 + mrSges(6,3) * t387 + Ifges(6,4) * t423 + Ifges(6,2) * t422 + Ifges(6,6) * t435 - pkin(5) * t504 + pkin(11) * t508 + t497 * t372 + t493 * t373 - t450 * t415 + t460 * t417;
t416 = Ifges(6,4) * t450 + Ifges(6,2) * t449 + Ifges(6,6) * t460;
t360 = mrSges(6,2) * t390 - mrSges(6,3) * t386 + Ifges(6,1) * t423 + Ifges(6,4) * t422 + Ifges(6,5) * t435 - pkin(11) * t371 - t493 * t372 + t497 * t373 + t449 * t415 - t460 * t416;
t430 = Ifges(5,4) * t461 - Ifges(5,2) * t460 + Ifges(5,6) * t476;
t431 = Ifges(5,1) * t461 - Ifges(5,4) * t460 + Ifges(5,5) * t476;
t501 = mrSges(5,1) * t395 - mrSges(5,2) * t396 + Ifges(5,5) * t436 - Ifges(5,6) * t435 + Ifges(5,3) * t463 - pkin(4) * t381 + qJ(5) * t367 + t489 * t359 + t485 * t360 + t461 * t430 + t460 * t431;
t456 = Ifges(4,5) * t483 + (t495 * Ifges(4,1) + t498 * Ifges(4,4)) * t514;
t455 = Ifges(4,6) * t483 + (t495 * Ifges(4,4) + t498 * Ifges(4,2)) * t514;
t429 = Ifges(5,5) * t461 - Ifges(5,6) * t460 + Ifges(5,3) * t476;
t425 = -t487 * t440 + t458;
t357 = m(4) * t425 + t469 * mrSges(4,1) + t468 * mrSges(4,2) + (t464 * t495 - t465 * t498) * t514 + t358;
t353 = (-Ifges(5,2) - Ifges(6,3)) * t435 - t502 + t476 * t431 - t461 * t429 + Ifges(5,6) * t463 + t449 * t417 - t450 * t416 + Ifges(5,4) * t436 - Ifges(6,6) * t422 - Ifges(6,5) * t423 - mrSges(5,1) * t406 + mrSges(5,3) * t396 - mrSges(6,1) * t386 + mrSges(6,2) * t387 - pkin(5) * t371 - pkin(4) * t366;
t352 = mrSges(5,2) * t406 - mrSges(5,3) * t395 + Ifges(5,1) * t436 - Ifges(5,4) * t435 + Ifges(5,5) * t463 - qJ(5) * t366 - t485 * t359 + t489 * t360 - t460 * t429 - t476 * t430;
t351 = Ifges(4,5) * t468 - Ifges(4,6) * t469 + Ifges(4,3) * t482 + mrSges(4,1) * t411 - mrSges(4,2) * t412 + t494 * t352 + t522 * t353 + pkin(3) * t503 + pkin(10) * t509 + (t495 * t455 - t498 * t456) * t514;
t1 = [m(2) * t484 + t492 * (m(3) * t462 + t491 * t357 + (t356 * t495 + t520) * t487) + (t496 * (m(3) * t448 - t500 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t510) + t499 * (m(3) * t447 + qJDD(2) * mrSges(3,1) - t500 * mrSges(3,2) - t487 * t357 + t515)) * t488; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t447 - mrSges(3,2) * t448 + t491 * t351 + pkin(2) * t515 + (t495 * (mrSges(4,2) * t425 - mrSges(4,3) * t411 + Ifges(4,1) * t468 - Ifges(4,4) * t469 + Ifges(4,5) * t482 - pkin(10) * t358 + t352 * t522 - t494 * t353 - t483 * t455) + t498 * (-mrSges(4,1) * t425 + mrSges(4,3) * t412 + Ifges(4,4) * t468 - Ifges(4,2) * t469 + Ifges(4,6) * t482 - pkin(3) * t358 + t483 * t456 - t501) - pkin(2) * t357 + pkin(9) * t510) * t487; t351; t501; t381; t502;];
tauJ  = t1;
