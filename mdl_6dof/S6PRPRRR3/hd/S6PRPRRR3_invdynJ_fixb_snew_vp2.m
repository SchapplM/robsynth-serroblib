% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:40:43
% EndTime: 2019-05-05 00:40:46
% DurationCPUTime: 3.44s
% Computational Cost: add. (33812->263), mult. (76491->338), div. (0->0), fcn. (60140->14), ass. (0->122)
t483 = sin(pkin(11));
t486 = cos(pkin(11));
t467 = g(1) * t483 - g(2) * t486;
t481 = -g(3) + qJDD(1);
t484 = sin(pkin(6));
t487 = cos(pkin(6));
t523 = t467 * t487 + t481 * t484;
t496 = qJD(2) ^ 2;
t468 = -g(1) * t486 - g(2) * t483;
t491 = sin(qJ(2));
t495 = cos(qJ(2));
t445 = -t468 * t491 + t495 * t523;
t485 = cos(pkin(12));
t479 = t485 ^ 2;
t522 = 0.2e1 * t485;
t521 = pkin(3) * t485;
t482 = sin(pkin(12));
t520 = mrSges(4,2) * t482;
t518 = t479 * t496;
t446 = t495 * t468 + t491 * t523;
t438 = -pkin(2) * t496 + qJDD(2) * qJ(3) + t446;
t458 = -t467 * t484 + t481 * t487;
t513 = qJD(2) * qJD(3);
t515 = t458 * t485 - 0.2e1 * t482 * t513;
t418 = (-pkin(8) * qJDD(2) + t496 * t521 - t438) * t482 + t515;
t430 = t438 * t485 + t458 * t482 + t513 * t522;
t512 = qJDD(2) * t485;
t421 = -pkin(3) * t518 + pkin(8) * t512 + t430;
t490 = sin(qJ(4));
t494 = cos(qJ(4));
t395 = t418 * t494 - t490 * t421;
t506 = t482 * t494 + t485 * t490;
t505 = -t482 * t490 + t485 * t494;
t460 = t505 * qJD(2);
t514 = t460 * qJD(4);
t452 = qJDD(2) * t506 + t514;
t461 = t506 * qJD(2);
t392 = (-t452 + t514) * pkin(9) + (t460 * t461 + qJDD(4)) * pkin(4) + t395;
t396 = t418 * t490 + t421 * t494;
t451 = -t461 * qJD(4) + qJDD(2) * t505;
t457 = qJD(4) * pkin(4) - pkin(9) * t461;
t459 = t460 ^ 2;
t394 = -pkin(4) * t459 + pkin(9) * t451 - qJD(4) * t457 + t396;
t489 = sin(qJ(5));
t493 = cos(qJ(5));
t389 = t392 * t489 + t394 * t493;
t444 = t460 * t489 + t461 * t493;
t412 = -qJD(5) * t444 + t451 * t493 - t452 * t489;
t443 = t460 * t493 - t461 * t489;
t427 = -mrSges(6,1) * t443 + mrSges(6,2) * t444;
t480 = qJD(4) + qJD(5);
t437 = mrSges(6,1) * t480 - mrSges(6,3) * t444;
t477 = qJDD(4) + qJDD(5);
t428 = -pkin(5) * t443 - pkin(10) * t444;
t476 = t480 ^ 2;
t386 = -pkin(5) * t476 + pkin(10) * t477 + t428 * t443 + t389;
t478 = t482 ^ 2;
t501 = qJDD(3) - t445;
t431 = (-pkin(2) - t521) * qJDD(2) + (-qJ(3) + (-t478 - t479) * pkin(8)) * t496 + t501;
t401 = -t451 * pkin(4) - t459 * pkin(9) + t457 * t461 + t431;
t413 = qJD(5) * t443 + t451 * t489 + t452 * t493;
t390 = (-t443 * t480 - t413) * pkin(10) + (t444 * t480 - t412) * pkin(5) + t401;
t488 = sin(qJ(6));
t492 = cos(qJ(6));
t383 = -t386 * t488 + t390 * t492;
t432 = -t444 * t488 + t480 * t492;
t399 = qJD(6) * t432 + t413 * t492 + t477 * t488;
t411 = qJDD(6) - t412;
t433 = t444 * t492 + t480 * t488;
t414 = -mrSges(7,1) * t432 + mrSges(7,2) * t433;
t439 = qJD(6) - t443;
t419 = -mrSges(7,2) * t439 + mrSges(7,3) * t432;
t380 = m(7) * t383 + mrSges(7,1) * t411 - mrSges(7,3) * t399 - t414 * t433 + t419 * t439;
t384 = t386 * t492 + t390 * t488;
t398 = -qJD(6) * t433 - t413 * t488 + t477 * t492;
t420 = mrSges(7,1) * t439 - mrSges(7,3) * t433;
t381 = m(7) * t384 - mrSges(7,2) * t411 + mrSges(7,3) * t398 + t414 * t432 - t420 * t439;
t508 = -t380 * t488 + t381 * t492;
t367 = m(6) * t389 - mrSges(6,2) * t477 + mrSges(6,3) * t412 + t427 * t443 - t437 * t480 + t508;
t388 = t392 * t493 - t394 * t489;
t436 = -mrSges(6,2) * t480 + mrSges(6,3) * t443;
t385 = -pkin(5) * t477 - pkin(10) * t476 + t428 * t444 - t388;
t502 = -m(7) * t385 + mrSges(7,1) * t398 - mrSges(7,2) * t399 + t419 * t432 - t420 * t433;
t376 = m(6) * t388 + mrSges(6,1) * t477 - mrSges(6,3) * t413 - t427 * t444 + t436 * t480 + t502;
t364 = t367 * t489 + t376 * t493;
t449 = -mrSges(5,1) * t460 + mrSges(5,2) * t461;
t455 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t460;
t362 = m(5) * t395 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t452 + qJD(4) * t455 - t449 * t461 + t364;
t456 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t461;
t509 = t367 * t493 - t376 * t489;
t363 = m(5) * t396 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t451 - qJD(4) * t456 + t449 * t460 + t509;
t516 = t362 * t494 + t363 * t490;
t370 = t380 * t492 + t381 * t488;
t429 = -t438 * t482 + t515;
t504 = mrSges(4,3) * qJDD(2) + t496 * (-mrSges(4,1) * t485 + t520);
t355 = m(4) * t429 - t482 * t504 + t516;
t510 = -t490 * t362 + t363 * t494;
t356 = m(4) * t430 + t485 * t504 + t510;
t511 = -t355 * t482 + t356 * t485;
t503 = m(6) * t401 - mrSges(6,1) * t412 + mrSges(6,2) * t413 - t436 * t443 + t437 * t444 + t370;
t402 = Ifges(7,5) * t433 + Ifges(7,6) * t432 + Ifges(7,3) * t439;
t404 = Ifges(7,1) * t433 + Ifges(7,4) * t432 + Ifges(7,5) * t439;
t373 = -mrSges(7,1) * t385 + mrSges(7,3) * t384 + Ifges(7,4) * t399 + Ifges(7,2) * t398 + Ifges(7,6) * t411 - t402 * t433 + t404 * t439;
t403 = Ifges(7,4) * t433 + Ifges(7,2) * t432 + Ifges(7,6) * t439;
t374 = mrSges(7,2) * t385 - mrSges(7,3) * t383 + Ifges(7,1) * t399 + Ifges(7,4) * t398 + Ifges(7,5) * t411 + t402 * t432 - t403 * t439;
t423 = Ifges(6,4) * t444 + Ifges(6,2) * t443 + Ifges(6,6) * t480;
t424 = Ifges(6,1) * t444 + Ifges(6,4) * t443 + Ifges(6,5) * t480;
t500 = mrSges(6,1) * t388 - mrSges(6,2) * t389 + Ifges(6,5) * t413 + Ifges(6,6) * t412 + Ifges(6,3) * t477 + pkin(5) * t502 + pkin(10) * t508 + t373 * t492 + t374 * t488 + t423 * t444 - t424 * t443;
t499 = m(5) * t431 - mrSges(5,1) * t451 + t452 * mrSges(5,2) - t455 * t460 + t461 * t456 + t503;
t498 = mrSges(7,1) * t383 - mrSges(7,2) * t384 + Ifges(7,5) * t399 + Ifges(7,6) * t398 + Ifges(7,3) * t411 + t403 * t433 - t404 * t432;
t435 = -qJDD(2) * pkin(2) - t496 * qJ(3) + t501;
t497 = -m(4) * t435 + mrSges(4,1) * t512 - t499 + (t478 * t496 + t518) * mrSges(4,3);
t442 = Ifges(5,1) * t461 + Ifges(5,4) * t460 + Ifges(5,5) * qJD(4);
t441 = Ifges(5,4) * t461 + Ifges(5,2) * t460 + Ifges(5,6) * qJD(4);
t440 = Ifges(5,5) * t461 + Ifges(5,6) * t460 + Ifges(5,3) * qJD(4);
t422 = Ifges(6,5) * t444 + Ifges(6,6) * t443 + Ifges(6,3) * t480;
t368 = qJDD(2) * t520 - t497;
t358 = -mrSges(6,1) * t401 + mrSges(6,3) * t389 + Ifges(6,4) * t413 + Ifges(6,2) * t412 + Ifges(6,6) * t477 - pkin(5) * t370 - t422 * t444 + t424 * t480 - t498;
t357 = mrSges(6,2) * t401 - mrSges(6,3) * t388 + Ifges(6,1) * t413 + Ifges(6,4) * t412 + Ifges(6,5) * t477 - pkin(10) * t370 - t373 * t488 + t374 * t492 + t422 * t443 - t423 * t480;
t353 = mrSges(5,2) * t431 - mrSges(5,3) * t395 + Ifges(5,1) * t452 + Ifges(5,4) * t451 + Ifges(5,5) * qJDD(4) - pkin(9) * t364 - qJD(4) * t441 + t357 * t493 - t358 * t489 + t440 * t460;
t352 = -mrSges(5,1) * t431 + mrSges(5,3) * t396 + Ifges(5,4) * t452 + Ifges(5,2) * t451 + Ifges(5,6) * qJDD(4) - pkin(4) * t503 + pkin(9) * t509 + qJD(4) * t442 + t489 * t357 + t493 * t358 - t461 * t440;
t1 = [m(2) * t481 + t487 * (m(3) * t458 + t355 * t485 + t356 * t482) + (t491 * (m(3) * t446 - mrSges(3,1) * t496 - qJDD(2) * mrSges(3,2) + t511) + t495 * (-t496 * mrSges(3,2) + m(3) * t445 + t497 + (mrSges(3,1) - t520) * qJDD(2))) * t484; mrSges(3,1) * t445 - mrSges(3,2) * t446 + t482 * (mrSges(4,2) * t435 - mrSges(4,3) * t429 - pkin(8) * t516 - t490 * t352 + t494 * t353) + t485 * (-mrSges(4,1) * t435 + mrSges(4,3) * t430 - pkin(3) * t499 + pkin(8) * t510 + t494 * t352 + t490 * t353) - pkin(2) * t368 + qJ(3) * t511 + (Ifges(4,2) * t479 + Ifges(3,3) + (Ifges(4,1) * t482 + Ifges(4,4) * t522) * t482) * qJDD(2); t368; mrSges(5,1) * t395 - mrSges(5,2) * t396 + Ifges(5,5) * t452 + Ifges(5,6) * t451 + Ifges(5,3) * qJDD(4) + pkin(4) * t364 + t441 * t461 - t442 * t460 + t500; t500; t498;];
tauJ  = t1;
