% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:58:54
% EndTime: 2019-05-05 06:58:58
% DurationCPUTime: 3.97s
% Computational Cost: add. (41306->291), mult. (88300->376), div. (0->0), fcn. (65683->14), ass. (0->123)
t491 = sin(pkin(11));
t494 = cos(pkin(11));
t477 = g(1) * t491 - g(2) * t494;
t489 = -g(3) + qJDD(1);
t492 = sin(pkin(6));
t495 = cos(pkin(6));
t524 = t477 * t495 + t489 * t492;
t478 = -g(1) * t494 - g(2) * t491;
t499 = sin(qJ(2));
t503 = cos(qJ(2));
t449 = -t499 * t478 + t524 * t503;
t523 = 2 * qJD(5);
t450 = t503 * t478 + t524 * t499;
t504 = qJD(2) ^ 2;
t443 = -pkin(2) * t504 + qJDD(2) * pkin(8) + t450;
t461 = -t477 * t492 + t489 * t495;
t498 = sin(qJ(3));
t502 = cos(qJ(3));
t429 = -t443 * t498 + t502 * t461;
t518 = qJD(2) * qJD(3);
t517 = t502 * t518;
t475 = qJDD(2) * t498 + t517;
t414 = (-t475 + t517) * pkin(9) + (t498 * t502 * t504 + qJDD(3)) * pkin(3) + t429;
t430 = t502 * t443 + t498 * t461;
t476 = qJDD(2) * t502 - t498 * t518;
t520 = qJD(2) * t498;
t482 = qJD(3) * pkin(3) - pkin(9) * t520;
t488 = t502 ^ 2;
t417 = -pkin(3) * t488 * t504 + pkin(9) * t476 - qJD(3) * t482 + t430;
t497 = sin(qJ(4));
t501 = cos(qJ(4));
t394 = t501 * t414 - t417 * t497;
t467 = (-t498 * t497 + t502 * t501) * qJD(2);
t439 = qJD(4) * t467 + t475 * t501 + t476 * t497;
t468 = (t502 * t497 + t498 * t501) * qJD(2);
t486 = qJDD(3) + qJDD(4);
t487 = qJD(3) + qJD(4);
t390 = (t467 * t487 - t439) * qJ(5) + (t467 * t468 + t486) * pkin(4) + t394;
t395 = t497 * t414 + t501 * t417;
t438 = -qJD(4) * t468 - t475 * t497 + t476 * t501;
t459 = pkin(4) * t487 - qJ(5) * t468;
t463 = t467 ^ 2;
t392 = -pkin(4) * t463 + qJ(5) * t438 - t459 * t487 + t395;
t490 = sin(pkin(12));
t493 = cos(pkin(12));
t453 = t467 * t493 - t468 * t490;
t387 = t490 * t390 + t493 * t392 + t453 * t523;
t415 = t438 * t493 - t439 * t490;
t454 = t467 * t490 + t468 * t493;
t426 = -mrSges(6,1) * t453 + mrSges(6,2) * t454;
t441 = mrSges(6,1) * t487 - mrSges(6,3) * t454;
t427 = -pkin(5) * t453 - pkin(10) * t454;
t485 = t487 ^ 2;
t384 = -pkin(5) * t485 + pkin(10) * t486 + t427 * t453 + t387;
t509 = -qJDD(2) * pkin(2) - t449;
t428 = -pkin(3) * t476 + t482 * t520 + (-pkin(9) * t488 - pkin(8)) * t504 + t509;
t397 = -pkin(4) * t438 - qJ(5) * t463 + t468 * t459 + qJDD(5) + t428;
t416 = t438 * t490 + t439 * t493;
t388 = (-t453 * t487 - t416) * pkin(10) + (t454 * t487 - t415) * pkin(5) + t397;
t496 = sin(qJ(6));
t500 = cos(qJ(6));
t381 = -t384 * t496 + t388 * t500;
t435 = -t454 * t496 + t487 * t500;
t400 = qJD(6) * t435 + t416 * t500 + t486 * t496;
t413 = qJDD(6) - t415;
t436 = t454 * t500 + t487 * t496;
t418 = -mrSges(7,1) * t435 + mrSges(7,2) * t436;
t445 = qJD(6) - t453;
t419 = -mrSges(7,2) * t445 + mrSges(7,3) * t435;
t378 = m(7) * t381 + mrSges(7,1) * t413 - mrSges(7,3) * t400 - t418 * t436 + t419 * t445;
t382 = t384 * t500 + t388 * t496;
t399 = -qJD(6) * t436 - t416 * t496 + t486 * t500;
t420 = mrSges(7,1) * t445 - mrSges(7,3) * t436;
t379 = m(7) * t382 - mrSges(7,2) * t413 + mrSges(7,3) * t399 + t418 * t435 - t420 * t445;
t513 = -t378 * t496 + t500 * t379;
t365 = m(6) * t387 - mrSges(6,2) * t486 + mrSges(6,3) * t415 + t426 * t453 - t441 * t487 + t513;
t512 = -t390 * t493 + t392 * t490;
t386 = -0.2e1 * qJD(5) * t454 - t512;
t440 = -mrSges(6,2) * t487 + mrSges(6,3) * t453;
t383 = -pkin(5) * t486 - pkin(10) * t485 + (t523 + t427) * t454 + t512;
t510 = -m(7) * t383 + t399 * mrSges(7,1) - mrSges(7,2) * t400 + t435 * t419 - t420 * t436;
t374 = m(6) * t386 + mrSges(6,1) * t486 - mrSges(6,3) * t416 - t426 * t454 + t440 * t487 + t510;
t362 = t490 * t365 + t493 * t374;
t455 = -mrSges(5,1) * t467 + mrSges(5,2) * t468;
t458 = -mrSges(5,2) * t487 + mrSges(5,3) * t467;
t359 = m(5) * t394 + mrSges(5,1) * t486 - mrSges(5,3) * t439 - t455 * t468 + t458 * t487 + t362;
t460 = mrSges(5,1) * t487 - mrSges(5,3) * t468;
t514 = t493 * t365 - t374 * t490;
t360 = m(5) * t395 - mrSges(5,2) * t486 + mrSges(5,3) * t438 + t455 * t467 - t460 * t487 + t514;
t353 = t501 * t359 + t497 * t360;
t368 = t500 * t378 + t496 * t379;
t519 = qJD(2) * t502;
t474 = (-t502 * mrSges(4,1) + t498 * mrSges(4,2)) * qJD(2);
t480 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t519;
t351 = m(4) * t429 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t475 + qJD(3) * t480 - t474 * t520 + t353;
t479 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t520;
t515 = -t359 * t497 + t501 * t360;
t352 = m(4) * t430 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t476 - qJD(3) * t479 + t474 * t519 + t515;
t516 = -t351 * t498 + t502 * t352;
t366 = m(6) * t397 - t415 * mrSges(6,1) + t416 * mrSges(6,2) - t453 * t440 + t454 * t441 + t368;
t508 = m(5) * t428 - t438 * mrSges(5,1) + mrSges(5,2) * t439 - t467 * t458 + t460 * t468 + t366;
t402 = Ifges(7,4) * t436 + Ifges(7,2) * t435 + Ifges(7,6) * t445;
t403 = Ifges(7,1) * t436 + Ifges(7,4) * t435 + Ifges(7,5) * t445;
t507 = mrSges(7,1) * t381 - mrSges(7,2) * t382 + Ifges(7,5) * t400 + Ifges(7,6) * t399 + Ifges(7,3) * t413 + t402 * t436 - t403 * t435;
t401 = Ifges(7,5) * t436 + Ifges(7,6) * t435 + Ifges(7,3) * t445;
t371 = -mrSges(7,1) * t383 + mrSges(7,3) * t382 + Ifges(7,4) * t400 + Ifges(7,2) * t399 + Ifges(7,6) * t413 - t401 * t436 + t403 * t445;
t372 = mrSges(7,2) * t383 - mrSges(7,3) * t381 + Ifges(7,1) * t400 + Ifges(7,4) * t399 + Ifges(7,5) * t413 + t401 * t435 - t402 * t445;
t422 = Ifges(6,4) * t454 + Ifges(6,2) * t453 + Ifges(6,6) * t487;
t423 = Ifges(6,1) * t454 + Ifges(6,4) * t453 + Ifges(6,5) * t487;
t447 = Ifges(5,4) * t468 + Ifges(5,2) * t467 + Ifges(5,6) * t487;
t448 = Ifges(5,1) * t468 + Ifges(5,4) * t467 + Ifges(5,5) * t487;
t506 = mrSges(5,1) * t394 + mrSges(6,1) * t386 - mrSges(5,2) * t395 - mrSges(6,2) * t387 + pkin(4) * t362 + pkin(5) * t510 + pkin(10) * t513 + t500 * t371 + t496 * t372 + t454 * t422 - t453 * t423 - t467 * t448 + Ifges(6,6) * t415 + Ifges(6,5) * t416 + t468 * t447 + Ifges(5,6) * t438 + Ifges(5,5) * t439 + (Ifges(6,3) + Ifges(5,3)) * t486;
t442 = -pkin(8) * t504 + t509;
t505 = -m(4) * t442 + t476 * mrSges(4,1) - mrSges(4,2) * t475 - t479 * t520 + t480 * t519 - t508;
t466 = Ifges(4,5) * qJD(3) + (t498 * Ifges(4,1) + t502 * Ifges(4,4)) * qJD(2);
t465 = Ifges(4,6) * qJD(3) + (t498 * Ifges(4,4) + t502 * Ifges(4,2)) * qJD(2);
t446 = Ifges(5,5) * t468 + Ifges(5,6) * t467 + Ifges(5,3) * t487;
t421 = Ifges(6,5) * t454 + Ifges(6,6) * t453 + Ifges(6,3) * t487;
t355 = -mrSges(6,1) * t397 + mrSges(6,3) * t387 + Ifges(6,4) * t416 + Ifges(6,2) * t415 + Ifges(6,6) * t486 - pkin(5) * t368 - t421 * t454 + t423 * t487 - t507;
t354 = mrSges(6,2) * t397 - mrSges(6,3) * t386 + Ifges(6,1) * t416 + Ifges(6,4) * t415 + Ifges(6,5) * t486 - pkin(10) * t368 - t371 * t496 + t372 * t500 + t421 * t453 - t422 * t487;
t349 = mrSges(5,2) * t428 - mrSges(5,3) * t394 + Ifges(5,1) * t439 + Ifges(5,4) * t438 + Ifges(5,5) * t486 - qJ(5) * t362 + t354 * t493 - t355 * t490 + t446 * t467 - t447 * t487;
t348 = -mrSges(5,1) * t428 + mrSges(5,3) * t395 + Ifges(5,4) * t439 + Ifges(5,2) * t438 + Ifges(5,6) * t486 - pkin(4) * t366 + qJ(5) * t514 + t490 * t354 + t493 * t355 - t468 * t446 + t487 * t448;
t1 = [m(2) * t489 + t495 * (m(3) * t461 + t351 * t502 + t352 * t498) + (t499 * (m(3) * t450 - mrSges(3,1) * t504 - qJDD(2) * mrSges(3,2) + t516) + t503 * (m(3) * t449 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t504 + t505)) * t492; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t449 - mrSges(3,2) * t450 + t498 * (mrSges(4,2) * t442 - mrSges(4,3) * t429 + Ifges(4,1) * t475 + Ifges(4,4) * t476 + Ifges(4,5) * qJDD(3) - pkin(9) * t353 - qJD(3) * t465 - t497 * t348 + t349 * t501) + t502 * (-mrSges(4,1) * t442 + mrSges(4,3) * t430 + Ifges(4,4) * t475 + Ifges(4,2) * t476 + Ifges(4,6) * qJDD(3) - pkin(3) * t508 + pkin(9) * t515 + qJD(3) * t466 + t501 * t348 + t497 * t349) + pkin(2) * t505 + pkin(8) * t516; Ifges(4,3) * qJDD(3) + Ifges(4,5) * t475 + Ifges(4,6) * t476 + mrSges(4,1) * t429 - mrSges(4,2) * t430 + pkin(3) * t353 + (t498 * t465 - t502 * t466) * qJD(2) + t506; t506; t366; t507;];
tauJ  = t1;
