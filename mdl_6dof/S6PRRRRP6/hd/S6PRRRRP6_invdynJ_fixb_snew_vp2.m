% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:12:06
% EndTime: 2019-05-05 10:12:10
% DurationCPUTime: 3.64s
% Computational Cost: add. (37308->277), mult. (76364->354), div. (0->0), fcn. (60537->14), ass. (0->125)
t524 = Ifges(6,1) + Ifges(7,1);
t516 = Ifges(6,4) - Ifges(7,5);
t515 = -Ifges(6,5) - Ifges(7,4);
t523 = Ifges(6,2) + Ifges(7,3);
t514 = Ifges(6,6) - Ifges(7,6);
t522 = -Ifges(6,3) - Ifges(7,2);
t476 = sin(pkin(12));
t479 = cos(pkin(12));
t468 = g(1) * t476 - g(2) * t479;
t475 = -g(3) + qJDD(1);
t478 = sin(pkin(6));
t481 = cos(pkin(6));
t521 = t468 * t481 + t475 * t478;
t477 = sin(pkin(7));
t484 = sin(qJ(3));
t487 = cos(qJ(3));
t501 = qJD(2) * qJD(3);
t461 = (-qJDD(2) * t487 + t484 * t501) * t477;
t469 = -g(1) * t479 - g(2) * t476;
t485 = sin(qJ(2));
t488 = cos(qJ(2));
t437 = -t469 * t485 + t521 * t488;
t489 = qJD(2) ^ 2;
t518 = t477 * pkin(9);
t433 = qJDD(2) * pkin(2) + t489 * t518 + t437;
t438 = t488 * t469 + t521 * t485;
t434 = -pkin(2) * t489 + qJDD(2) * t518 + t438;
t480 = cos(pkin(7));
t454 = -t468 * t478 + t475 * t481;
t511 = t454 * t477;
t395 = -t484 * t434 + t487 * (t433 * t480 + t511);
t474 = qJD(2) * t480 + qJD(3);
t483 = sin(qJ(4));
t486 = cos(qJ(4));
t502 = qJD(2) * t477;
t499 = t484 * t502;
t452 = t474 * t486 - t483 * t499;
t460 = (qJDD(2) * t484 + t487 * t501) * t477;
t473 = qJDD(2) * t480 + qJDD(3);
t429 = qJD(4) * t452 + t460 * t486 + t473 * t483;
t453 = t474 * t483 + t486 * t499;
t498 = t487 * t502;
t467 = qJD(4) - t498;
t482 = sin(qJ(5));
t519 = cos(qJ(5));
t439 = t453 * t482 - t467 * t519;
t455 = qJDD(4) + t461;
t400 = -t439 * qJD(5) + t429 * t519 + t482 * t455;
t440 = t453 * t519 + t482 * t467;
t412 = mrSges(7,1) * t439 - mrSges(7,3) * t440;
t508 = t480 * t484;
t396 = t433 * t508 + t487 * t434 + t484 * t511;
t459 = (-t487 * pkin(3) - t484 * pkin(10)) * t502;
t472 = t474 ^ 2;
t392 = -pkin(3) * t472 + pkin(10) * t473 + t459 * t498 + t396;
t448 = t480 * t454;
t394 = pkin(3) * t461 - pkin(10) * t460 + t448 + (-t433 + (pkin(3) * t484 - pkin(10) * t487) * t474 * qJD(2)) * t477;
t388 = t486 * t392 + t483 * t394;
t436 = -pkin(4) * t452 - pkin(11) * t453;
t466 = t467 ^ 2;
t384 = -pkin(4) * t466 + pkin(11) * t455 + t436 * t452 + t388;
t391 = -pkin(3) * t473 - pkin(10) * t472 + t459 * t499 - t395;
t428 = -qJD(4) * t453 - t460 * t483 + t473 * t486;
t386 = (-t452 * t467 - t429) * pkin(11) + (t453 * t467 - t428) * pkin(4) + t391;
t380 = -t482 * t384 + t386 * t519;
t411 = pkin(5) * t439 - qJ(6) * t440;
t426 = qJDD(5) - t428;
t450 = qJD(5) - t452;
t449 = t450 ^ 2;
t378 = -t426 * pkin(5) - t449 * qJ(6) + t440 * t411 + qJDD(6) - t380;
t416 = -mrSges(7,2) * t439 + mrSges(7,3) * t450;
t495 = -m(7) * t378 + t426 * mrSges(7,1) + t450 * t416;
t374 = mrSges(7,2) * t400 + t412 * t440 - t495;
t381 = t384 * t519 + t482 * t386;
t377 = -pkin(5) * t449 + qJ(6) * t426 + 0.2e1 * qJD(6) * t450 - t411 * t439 + t381;
t399 = qJD(5) * t440 + t429 * t482 - t455 * t519;
t419 = -mrSges(7,1) * t450 + mrSges(7,2) * t440;
t500 = m(7) * t377 + t426 * mrSges(7,3) + t450 * t419;
t504 = t516 * t439 - t524 * t440 + t515 * t450;
t505 = t523 * t439 - t516 * t440 - t514 * t450;
t520 = -t399 * t514 - t400 * t515 - t522 * t426 - t439 * t504 - t440 * t505 + mrSges(6,1) * t380 - mrSges(7,1) * t378 - mrSges(6,2) * t381 + mrSges(7,3) * t377 - pkin(5) * t374 + qJ(6) * (-mrSges(7,2) * t399 - t412 * t439 + t500);
t517 = -mrSges(6,3) - mrSges(7,2);
t457 = -mrSges(4,2) * t474 + mrSges(4,3) * t498;
t458 = (-t487 * mrSges(4,1) + t484 * mrSges(4,2)) * t502;
t418 = mrSges(6,1) * t450 - mrSges(6,3) * t440;
t503 = -mrSges(6,1) * t439 - mrSges(6,2) * t440 - t412;
t370 = m(6) * t381 - mrSges(6,2) * t426 + t399 * t517 - t418 * t450 + t439 * t503 + t500;
t417 = -mrSges(6,2) * t450 - mrSges(6,3) * t439;
t371 = m(6) * t380 + mrSges(6,1) * t426 + t400 * t517 + t417 * t450 + t440 * t503 + t495;
t365 = t482 * t370 + t371 * t519;
t441 = -mrSges(5,2) * t467 + mrSges(5,3) * t452;
t442 = mrSges(5,1) * t467 - mrSges(5,3) * t453;
t491 = -m(5) * t391 + t428 * mrSges(5,1) - t429 * mrSges(5,2) + t452 * t441 - t453 * t442 - t365;
t359 = m(4) * t395 + t473 * mrSges(4,1) - t460 * mrSges(4,3) + t474 * t457 - t458 * t499 + t491;
t512 = t359 * t487;
t456 = mrSges(4,1) * t474 - mrSges(4,3) * t499;
t366 = t370 * t519 - t371 * t482;
t435 = -mrSges(5,1) * t452 + mrSges(5,2) * t453;
t362 = m(5) * t388 - mrSges(5,2) * t455 + mrSges(5,3) * t428 + t435 * t452 - t442 * t467 + t366;
t387 = -t483 * t392 + t394 * t486;
t383 = -pkin(4) * t455 - pkin(11) * t466 + t453 * t436 - t387;
t379 = -0.2e1 * qJD(6) * t440 + (t439 * t450 - t400) * qJ(6) + (t440 * t450 + t399) * pkin(5) + t383;
t375 = m(7) * t379 + mrSges(7,1) * t399 - t400 * mrSges(7,3) + t416 * t439 - t440 * t419;
t372 = -m(6) * t383 - t399 * mrSges(6,1) - mrSges(6,2) * t400 - t439 * t417 - t418 * t440 - t375;
t368 = m(5) * t387 + mrSges(5,1) * t455 - mrSges(5,3) * t429 - t435 * t453 + t441 * t467 + t372;
t496 = t486 * t362 - t368 * t483;
t355 = m(4) * t396 - mrSges(4,2) * t473 - mrSges(4,3) * t461 - t456 * t474 + t458 * t498 + t496;
t507 = t355 * t508 + t480 * t512;
t357 = t483 * t362 + t486 * t368;
t506 = t514 * t439 + t515 * t440 + t522 * t450;
t497 = t487 * t355 - t359 * t484;
t363 = -mrSges(6,1) * t383 - mrSges(7,1) * t379 + mrSges(7,2) * t377 + mrSges(6,3) * t381 - pkin(5) * t375 - t523 * t399 + t516 * t400 + t514 * t426 + t506 * t440 - t504 * t450;
t364 = mrSges(6,2) * t383 + mrSges(7,2) * t378 - mrSges(6,3) * t380 - mrSges(7,3) * t379 - qJ(6) * t375 - t516 * t399 + t524 * t400 - t515 * t426 + t506 * t439 + t505 * t450;
t423 = Ifges(5,4) * t453 + Ifges(5,2) * t452 + Ifges(5,6) * t467;
t424 = Ifges(5,1) * t453 + Ifges(5,4) * t452 + Ifges(5,5) * t467;
t490 = mrSges(5,1) * t387 - mrSges(5,2) * t388 + Ifges(5,5) * t429 + Ifges(5,6) * t428 + Ifges(5,3) * t455 + pkin(4) * t372 + pkin(11) * t366 + t363 * t519 + t482 * t364 + t453 * t423 - t452 * t424;
t446 = Ifges(4,5) * t474 + (t484 * Ifges(4,1) + t487 * Ifges(4,4)) * t502;
t445 = Ifges(4,6) * t474 + (t484 * Ifges(4,4) + t487 * Ifges(4,2)) * t502;
t422 = Ifges(5,5) * t453 + Ifges(5,6) * t452 + Ifges(5,3) * t467;
t414 = -t433 * t477 + t448;
t356 = m(4) * t414 + mrSges(4,1) * t461 + mrSges(4,2) * t460 + (t456 * t484 - t457 * t487) * t502 + t357;
t352 = -mrSges(5,1) * t391 + mrSges(5,3) * t388 + Ifges(5,4) * t429 + Ifges(5,2) * t428 + Ifges(5,6) * t455 - pkin(4) * t365 - t453 * t422 + t467 * t424 - t520;
t351 = mrSges(5,2) * t391 - mrSges(5,3) * t387 + Ifges(5,1) * t429 + Ifges(5,4) * t428 + Ifges(5,5) * t455 - pkin(11) * t365 - t482 * t363 + t364 * t519 + t452 * t422 - t467 * t423;
t350 = Ifges(4,5) * t460 - Ifges(4,6) * t461 + Ifges(4,3) * t473 + mrSges(4,1) * t395 - mrSges(4,2) * t396 + t483 * t351 + t486 * t352 + pkin(3) * t491 + pkin(10) * t496 + (t484 * t445 - t487 * t446) * t502;
t1 = [m(2) * t475 + t481 * (m(3) * t454 + t356 * t480 + (t355 * t484 + t512) * t477) + (t485 * (m(3) * t438 - mrSges(3,1) * t489 - qJDD(2) * mrSges(3,2) + t497) + t488 * (m(3) * t437 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t489 - t356 * t477 + t507)) * t478; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t437 - mrSges(3,2) * t438 + t480 * t350 + pkin(2) * t507 + (t484 * (mrSges(4,2) * t414 - mrSges(4,3) * t395 + Ifges(4,1) * t460 - Ifges(4,4) * t461 + Ifges(4,5) * t473 - pkin(10) * t357 + t351 * t486 - t352 * t483 - t445 * t474) + t487 * (-mrSges(4,1) * t414 + mrSges(4,3) * t396 + Ifges(4,4) * t460 - Ifges(4,2) * t461 + Ifges(4,6) * t473 - pkin(3) * t357 + t474 * t446 - t490) - pkin(2) * t356 + pkin(9) * t497) * t477; t350; t490; t520; t374;];
tauJ  = t1;
