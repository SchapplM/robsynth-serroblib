% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:42:36
% EndTime: 2019-05-05 10:42:40
% DurationCPUTime: 4.05s
% Computational Cost: add. (46498->292), mult. (90910->374), div. (0->0), fcn. (66941->14), ass. (0->124)
t492 = sin(pkin(12));
t494 = cos(pkin(12));
t479 = g(1) * t492 - g(2) * t494;
t491 = -g(3) + qJDD(1);
t493 = sin(pkin(6));
t495 = cos(pkin(6));
t526 = t479 * t495 + t491 * t493;
t480 = -g(1) * t494 - g(2) * t492;
t500 = sin(qJ(2));
t505 = cos(qJ(2));
t448 = -t500 * t480 + t505 * t526;
t449 = t505 * t480 + t500 * t526;
t506 = qJD(2) ^ 2;
t444 = -pkin(2) * t506 + qJDD(2) * pkin(8) + t449;
t460 = -t479 * t493 + t491 * t495;
t499 = sin(qJ(3));
t504 = cos(qJ(3));
t424 = -t499 * t444 + t504 * t460;
t521 = qJD(2) * qJD(3);
t520 = t504 * t521;
t477 = qJDD(2) * t499 + t520;
t412 = (-t477 + t520) * pkin(9) + (t499 * t504 * t506 + qJDD(3)) * pkin(3) + t424;
t425 = t504 * t444 + t499 * t460;
t478 = qJDD(2) * t504 - t499 * t521;
t523 = qJD(2) * t499;
t484 = qJD(3) * pkin(3) - pkin(9) * t523;
t490 = t504 ^ 2;
t413 = -pkin(3) * t490 * t506 + pkin(9) * t478 - qJD(3) * t484 + t425;
t498 = sin(qJ(4));
t503 = cos(qJ(4));
t402 = t498 * t412 + t503 * t413;
t469 = (t498 * t504 + t499 * t503) * qJD(2);
t438 = -t469 * qJD(4) - t498 * t477 + t478 * t503;
t522 = qJD(2) * t504;
t468 = -t498 * t523 + t503 * t522;
t451 = -mrSges(5,1) * t468 + mrSges(5,2) * t469;
t489 = qJD(3) + qJD(4);
t459 = mrSges(5,1) * t489 - mrSges(5,3) * t469;
t488 = qJDD(3) + qJDD(4);
t452 = -pkin(4) * t468 - pkin(10) * t469;
t487 = t489 ^ 2;
t391 = -pkin(4) * t487 + pkin(10) * t488 + t452 * t468 + t402;
t512 = -qJDD(2) * pkin(2) - t448;
t419 = -t478 * pkin(3) + t484 * t523 + (-pkin(9) * t490 - pkin(8)) * t506 + t512;
t439 = qJD(4) * t468 + t477 * t503 + t478 * t498;
t399 = (-t468 * t489 - t439) * pkin(10) + (t469 * t489 - t438) * pkin(4) + t419;
t497 = sin(qJ(5));
t502 = cos(qJ(5));
t386 = -t497 * t391 + t502 * t399;
t454 = -t469 * t497 + t489 * t502;
t416 = qJD(5) * t454 + t439 * t502 + t488 * t497;
t436 = qJDD(5) - t438;
t455 = t469 * t502 + t489 * t497;
t463 = qJD(5) - t468;
t384 = (t454 * t463 - t416) * pkin(11) + (t454 * t455 + t436) * pkin(5) + t386;
t387 = t502 * t391 + t497 * t399;
t415 = -qJD(5) * t455 - t439 * t497 + t488 * t502;
t442 = pkin(5) * t463 - pkin(11) * t455;
t453 = t454 ^ 2;
t385 = -pkin(5) * t453 + pkin(11) * t415 - t442 * t463 + t387;
t496 = sin(qJ(6));
t501 = cos(qJ(6));
t382 = t384 * t501 - t385 * t496;
t426 = t454 * t501 - t455 * t496;
t396 = qJD(6) * t426 + t415 * t496 + t416 * t501;
t427 = t454 * t496 + t455 * t501;
t408 = -mrSges(7,1) * t426 + mrSges(7,2) * t427;
t461 = qJD(6) + t463;
t417 = -mrSges(7,2) * t461 + mrSges(7,3) * t426;
t431 = qJDD(6) + t436;
t378 = m(7) * t382 + mrSges(7,1) * t431 - mrSges(7,3) * t396 - t408 * t427 + t417 * t461;
t383 = t384 * t496 + t385 * t501;
t395 = -qJD(6) * t427 + t415 * t501 - t416 * t496;
t418 = mrSges(7,1) * t461 - mrSges(7,3) * t427;
t379 = m(7) * t383 - mrSges(7,2) * t431 + mrSges(7,3) * t395 + t408 * t426 - t418 * t461;
t370 = t501 * t378 + t496 * t379;
t428 = -mrSges(6,1) * t454 + mrSges(6,2) * t455;
t440 = -mrSges(6,2) * t463 + mrSges(6,3) * t454;
t368 = m(6) * t386 + mrSges(6,1) * t436 - mrSges(6,3) * t416 - t428 * t455 + t440 * t463 + t370;
t441 = mrSges(6,1) * t463 - mrSges(6,3) * t455;
t516 = -t378 * t496 + t501 * t379;
t369 = m(6) * t387 - mrSges(6,2) * t436 + mrSges(6,3) * t415 + t428 * t454 - t441 * t463 + t516;
t517 = -t368 * t497 + t502 * t369;
t362 = m(5) * t402 - mrSges(5,2) * t488 + mrSges(5,3) * t438 + t451 * t468 - t459 * t489 + t517;
t401 = t412 * t503 - t498 * t413;
t458 = -mrSges(5,2) * t489 + mrSges(5,3) * t468;
t390 = -pkin(4) * t488 - pkin(10) * t487 + t469 * t452 - t401;
t388 = -pkin(5) * t415 - pkin(11) * t453 + t442 * t455 + t390;
t514 = m(7) * t388 - t395 * mrSges(7,1) + mrSges(7,2) * t396 - t426 * t417 + t418 * t427;
t509 = -m(6) * t390 + t415 * mrSges(6,1) - mrSges(6,2) * t416 + t454 * t440 - t441 * t455 - t514;
t374 = m(5) * t401 + mrSges(5,1) * t488 - mrSges(5,3) * t439 - t451 * t469 + t458 * t489 + t509;
t355 = t498 * t362 + t503 * t374;
t364 = t502 * t368 + t497 * t369;
t476 = (-mrSges(4,1) * t504 + mrSges(4,2) * t499) * qJD(2);
t482 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t522;
t353 = m(4) * t424 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t477 + qJD(3) * t482 - t476 * t523 + t355;
t481 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t523;
t518 = t503 * t362 - t374 * t498;
t354 = m(4) * t425 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t478 - qJD(3) * t481 + t476 * t522 + t518;
t519 = -t353 * t499 + t504 * t354;
t405 = Ifges(7,4) * t427 + Ifges(7,2) * t426 + Ifges(7,6) * t461;
t406 = Ifges(7,1) * t427 + Ifges(7,4) * t426 + Ifges(7,5) * t461;
t513 = -mrSges(7,1) * t382 + mrSges(7,2) * t383 - Ifges(7,5) * t396 - Ifges(7,6) * t395 - Ifges(7,3) * t431 - t427 * t405 + t426 * t406;
t511 = m(5) * t419 - t438 * mrSges(5,1) + mrSges(5,2) * t439 - t468 * t458 + t459 * t469 + t364;
t404 = Ifges(7,5) * t427 + Ifges(7,6) * t426 + Ifges(7,3) * t461;
t371 = -mrSges(7,1) * t388 + mrSges(7,3) * t383 + Ifges(7,4) * t396 + Ifges(7,2) * t395 + Ifges(7,6) * t431 - t404 * t427 + t406 * t461;
t372 = mrSges(7,2) * t388 - mrSges(7,3) * t382 + Ifges(7,1) * t396 + Ifges(7,4) * t395 + Ifges(7,5) * t431 + t404 * t426 - t405 * t461;
t420 = Ifges(6,5) * t455 + Ifges(6,6) * t454 + Ifges(6,3) * t463;
t422 = Ifges(6,1) * t455 + Ifges(6,4) * t454 + Ifges(6,5) * t463;
t357 = -mrSges(6,1) * t390 + mrSges(6,3) * t387 + Ifges(6,4) * t416 + Ifges(6,2) * t415 + Ifges(6,6) * t436 - pkin(5) * t514 + pkin(11) * t516 + t501 * t371 + t496 * t372 - t455 * t420 + t463 * t422;
t421 = Ifges(6,4) * t455 + Ifges(6,2) * t454 + Ifges(6,6) * t463;
t359 = mrSges(6,2) * t390 - mrSges(6,3) * t386 + Ifges(6,1) * t416 + Ifges(6,4) * t415 + Ifges(6,5) * t436 - pkin(11) * t370 - t371 * t496 + t372 * t501 + t420 * t454 - t421 * t463;
t446 = Ifges(5,4) * t469 + Ifges(5,2) * t468 + Ifges(5,6) * t489;
t447 = Ifges(5,1) * t469 + Ifges(5,4) * t468 + Ifges(5,5) * t489;
t510 = mrSges(5,1) * t401 - mrSges(5,2) * t402 + Ifges(5,5) * t439 + Ifges(5,6) * t438 + Ifges(5,3) * t488 + pkin(4) * t509 + pkin(10) * t517 + t502 * t357 + t497 * t359 + t469 * t446 - t468 * t447;
t443 = -t506 * pkin(8) + t512;
t508 = -m(4) * t443 + t478 * mrSges(4,1) - mrSges(4,2) * t477 - t481 * t523 + t482 * t522 - t511;
t507 = mrSges(6,1) * t386 - mrSges(6,2) * t387 + Ifges(6,5) * t416 + Ifges(6,6) * t415 + Ifges(6,3) * t436 + pkin(5) * t370 + t455 * t421 - t454 * t422 - t513;
t467 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t499 + Ifges(4,4) * t504) * qJD(2);
t466 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t499 + Ifges(4,2) * t504) * qJD(2);
t445 = Ifges(5,5) * t469 + Ifges(5,6) * t468 + Ifges(5,3) * t489;
t351 = -mrSges(5,1) * t419 + mrSges(5,3) * t402 + Ifges(5,4) * t439 + Ifges(5,2) * t438 + Ifges(5,6) * t488 - pkin(4) * t364 - t469 * t445 + t489 * t447 - t507;
t350 = mrSges(5,2) * t419 - mrSges(5,3) * t401 + Ifges(5,1) * t439 + Ifges(5,4) * t438 + Ifges(5,5) * t488 - pkin(10) * t364 - t357 * t497 + t359 * t502 + t445 * t468 - t446 * t489;
t1 = [m(2) * t491 + t495 * (m(3) * t460 + t353 * t504 + t354 * t499) + (t500 * (m(3) * t449 - mrSges(3,1) * t506 - qJDD(2) * mrSges(3,2) + t519) + t505 * (m(3) * t448 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t506 + t508)) * t493; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t448 - mrSges(3,2) * t449 + t499 * (mrSges(4,2) * t443 - mrSges(4,3) * t424 + Ifges(4,1) * t477 + Ifges(4,4) * t478 + Ifges(4,5) * qJDD(3) - pkin(9) * t355 - qJD(3) * t466 + t350 * t503 - t351 * t498) + t504 * (-mrSges(4,1) * t443 + mrSges(4,3) * t425 + Ifges(4,4) * t477 + Ifges(4,2) * t478 + Ifges(4,6) * qJDD(3) - pkin(3) * t511 + pkin(9) * t518 + qJD(3) * t467 + t498 * t350 + t503 * t351) + pkin(2) * t508 + pkin(8) * t519; t510 + pkin(3) * t355 + (t466 * t499 - t467 * t504) * qJD(2) + Ifges(4,3) * qJDD(3) + Ifges(4,5) * t477 + Ifges(4,6) * t478 + mrSges(4,1) * t424 - mrSges(4,2) * t425; t510; t507; -t513;];
tauJ  = t1;
