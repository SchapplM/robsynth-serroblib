% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 21:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:55:55
% EndTime: 2019-05-05 21:55:59
% DurationCPUTime: 3.70s
% Computational Cost: add. (39583->293), mult. (84852->375), div. (0->0), fcn. (59253->12), ass. (0->121)
t505 = 2 * qJD(5);
t482 = sin(qJ(1));
t486 = cos(qJ(1));
t500 = t482 * g(1) - g(2) * t486;
t457 = qJDD(1) * pkin(1) + t500;
t487 = qJD(1) ^ 2;
t495 = -g(1) * t486 - g(2) * t482;
t459 = -pkin(1) * t487 + t495;
t476 = sin(pkin(10));
t478 = cos(pkin(10));
t441 = t476 * t457 + t478 * t459;
t438 = -pkin(2) * t487 + qJDD(1) * pkin(7) + t441;
t474 = -g(3) + qJDD(2);
t481 = sin(qJ(3));
t485 = cos(qJ(3));
t421 = -t481 * t438 + t485 * t474;
t502 = qJD(1) * qJD(3);
t501 = t485 * t502;
t460 = qJDD(1) * t481 + t501;
t412 = (-t460 + t501) * pkin(8) + (t481 * t485 * t487 + qJDD(3)) * pkin(3) + t421;
t422 = t485 * t438 + t481 * t474;
t461 = qJDD(1) * t485 - t481 * t502;
t504 = qJD(1) * t481;
t464 = qJD(3) * pkin(3) - pkin(8) * t504;
t473 = t485 ^ 2;
t413 = -pkin(3) * t473 * t487 + pkin(8) * t461 - qJD(3) * t464 + t422;
t480 = sin(qJ(4));
t484 = cos(qJ(4));
t385 = t484 * t412 - t480 * t413;
t452 = (-t480 * t481 + t484 * t485) * qJD(1);
t424 = qJD(4) * t452 + t460 * t484 + t461 * t480;
t453 = (t480 * t485 + t481 * t484) * qJD(1);
t471 = qJDD(3) + qJDD(4);
t472 = qJD(3) + qJD(4);
t376 = (t452 * t472 - t424) * qJ(5) + (t452 * t453 + t471) * pkin(4) + t385;
t386 = t480 * t412 + t484 * t413;
t423 = -qJD(4) * t453 - t460 * t480 + t461 * t484;
t443 = pkin(4) * t472 - qJ(5) * t453;
t445 = t452 ^ 2;
t378 = -pkin(4) * t445 + qJ(5) * t423 - t443 * t472 + t386;
t475 = sin(pkin(11));
t477 = cos(pkin(11));
t435 = t452 * t477 - t453 * t475;
t373 = t475 * t376 + t477 * t378 + t435 * t505;
t397 = t423 * t477 - t424 * t475;
t436 = t452 * t475 + t453 * t477;
t410 = -mrSges(6,1) * t435 + mrSges(6,2) * t436;
t426 = mrSges(6,1) * t472 - mrSges(6,3) * t436;
t411 = -pkin(5) * t435 - pkin(9) * t436;
t470 = t472 ^ 2;
t370 = -pkin(5) * t470 + pkin(9) * t471 + t411 * t435 + t373;
t440 = t478 * t457 - t476 * t459;
t493 = -qJDD(1) * pkin(2) - t440;
t414 = -t461 * pkin(3) + t464 * t504 + (-pkin(8) * t473 - pkin(7)) * t487 + t493;
t380 = -t423 * pkin(4) - t445 * qJ(5) + t453 * t443 + qJDD(5) + t414;
t398 = t423 * t475 + t424 * t477;
t374 = (-t435 * t472 - t398) * pkin(9) + (t436 * t472 - t397) * pkin(5) + t380;
t479 = sin(qJ(6));
t483 = cos(qJ(6));
t367 = -t370 * t479 + t374 * t483;
t419 = -t436 * t479 + t472 * t483;
t384 = qJD(6) * t419 + t398 * t483 + t471 * t479;
t396 = qJDD(6) - t397;
t420 = t436 * t483 + t472 * t479;
t399 = -mrSges(7,1) * t419 + mrSges(7,2) * t420;
t429 = qJD(6) - t435;
t400 = -mrSges(7,2) * t429 + mrSges(7,3) * t419;
t364 = m(7) * t367 + mrSges(7,1) * t396 - mrSges(7,3) * t384 - t399 * t420 + t400 * t429;
t368 = t370 * t483 + t374 * t479;
t383 = -qJD(6) * t420 - t398 * t479 + t471 * t483;
t401 = mrSges(7,1) * t429 - mrSges(7,3) * t420;
t365 = m(7) * t368 - mrSges(7,2) * t396 + mrSges(7,3) * t383 + t399 * t419 - t401 * t429;
t496 = -t364 * t479 + t483 * t365;
t351 = m(6) * t373 - mrSges(6,2) * t471 + mrSges(6,3) * t397 + t410 * t435 - t426 * t472 + t496;
t494 = -t477 * t376 + t475 * t378;
t372 = -0.2e1 * qJD(5) * t436 - t494;
t425 = -mrSges(6,2) * t472 + mrSges(6,3) * t435;
t369 = -t471 * pkin(5) - t470 * pkin(9) + (t505 + t411) * t436 + t494;
t492 = -m(7) * t369 + t383 * mrSges(7,1) - mrSges(7,2) * t384 + t419 * t400 - t401 * t420;
t360 = m(6) * t372 + mrSges(6,1) * t471 - mrSges(6,3) * t398 - t410 * t436 + t425 * t472 + t492;
t348 = t475 * t351 + t477 * t360;
t439 = -mrSges(5,1) * t452 + mrSges(5,2) * t453;
t442 = -mrSges(5,2) * t472 + mrSges(5,3) * t452;
t345 = m(5) * t385 + mrSges(5,1) * t471 - mrSges(5,3) * t424 - t439 * t453 + t442 * t472 + t348;
t444 = mrSges(5,1) * t472 - mrSges(5,3) * t453;
t497 = t477 * t351 - t360 * t475;
t346 = m(5) * t386 - mrSges(5,2) * t471 + mrSges(5,3) * t423 + t439 * t452 - t444 * t472 + t497;
t339 = t484 * t345 + t480 * t346;
t354 = t483 * t364 + t479 * t365;
t503 = qJD(1) * t485;
t458 = (-mrSges(4,1) * t485 + mrSges(4,2) * t481) * qJD(1);
t463 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t503;
t337 = m(4) * t421 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t460 + qJD(3) * t463 - t458 * t504 + t339;
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t504;
t498 = -t345 * t480 + t484 * t346;
t338 = m(4) * t422 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t461 - qJD(3) * t462 + t458 * t503 + t498;
t499 = -t337 * t481 + t485 * t338;
t352 = m(6) * t380 - t397 * mrSges(6,1) + t398 * mrSges(6,2) - t435 * t425 + t436 * t426 + t354;
t491 = m(5) * t414 - t423 * mrSges(5,1) + mrSges(5,2) * t424 - t452 * t442 + t444 * t453 + t352;
t388 = Ifges(7,4) * t420 + Ifges(7,2) * t419 + Ifges(7,6) * t429;
t389 = Ifges(7,1) * t420 + Ifges(7,4) * t419 + Ifges(7,5) * t429;
t490 = mrSges(7,1) * t367 - mrSges(7,2) * t368 + Ifges(7,5) * t384 + Ifges(7,6) * t383 + Ifges(7,3) * t396 + t388 * t420 - t389 * t419;
t387 = Ifges(7,5) * t420 + Ifges(7,6) * t419 + Ifges(7,3) * t429;
t357 = -mrSges(7,1) * t369 + mrSges(7,3) * t368 + Ifges(7,4) * t384 + Ifges(7,2) * t383 + Ifges(7,6) * t396 - t387 * t420 + t389 * t429;
t358 = mrSges(7,2) * t369 - mrSges(7,3) * t367 + Ifges(7,1) * t384 + Ifges(7,4) * t383 + Ifges(7,5) * t396 + t387 * t419 - t388 * t429;
t403 = Ifges(6,4) * t436 + Ifges(6,2) * t435 + Ifges(6,6) * t472;
t404 = Ifges(6,1) * t436 + Ifges(6,4) * t435 + Ifges(6,5) * t472;
t431 = Ifges(5,4) * t453 + Ifges(5,2) * t452 + Ifges(5,6) * t472;
t432 = Ifges(5,1) * t453 + Ifges(5,4) * t452 + Ifges(5,5) * t472;
t489 = mrSges(5,1) * t385 + mrSges(6,1) * t372 - mrSges(5,2) * t386 - mrSges(6,2) * t373 + pkin(4) * t348 + pkin(5) * t492 + pkin(9) * t496 + t483 * t357 + t479 * t358 + t436 * t403 - t435 * t404 - t452 * t432 + Ifges(6,6) * t397 + Ifges(6,5) * t398 + t453 * t431 + Ifges(5,6) * t423 + Ifges(5,5) * t424 + (Ifges(6,3) + Ifges(5,3)) * t471;
t437 = -t487 * pkin(7) + t493;
t488 = -m(4) * t437 + t461 * mrSges(4,1) - mrSges(4,2) * t460 - t462 * t504 + t463 * t503 - t491;
t451 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t481 + Ifges(4,4) * t485) * qJD(1);
t450 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t481 + Ifges(4,2) * t485) * qJD(1);
t430 = Ifges(5,5) * t453 + Ifges(5,6) * t452 + Ifges(5,3) * t472;
t402 = Ifges(6,5) * t436 + Ifges(6,6) * t435 + Ifges(6,3) * t472;
t341 = -mrSges(6,1) * t380 + mrSges(6,3) * t373 + Ifges(6,4) * t398 + Ifges(6,2) * t397 + Ifges(6,6) * t471 - pkin(5) * t354 - t402 * t436 + t404 * t472 - t490;
t340 = mrSges(6,2) * t380 - mrSges(6,3) * t372 + Ifges(6,1) * t398 + Ifges(6,4) * t397 + Ifges(6,5) * t471 - pkin(9) * t354 - t357 * t479 + t358 * t483 + t402 * t435 - t403 * t472;
t335 = mrSges(5,2) * t414 - mrSges(5,3) * t385 + Ifges(5,1) * t424 + Ifges(5,4) * t423 + Ifges(5,5) * t471 - qJ(5) * t348 + t340 * t477 - t341 * t475 + t430 * t452 - t431 * t472;
t334 = -mrSges(5,1) * t414 + mrSges(5,3) * t386 + Ifges(5,4) * t424 + Ifges(5,2) * t423 + Ifges(5,6) * t471 - pkin(4) * t352 + qJ(5) * t497 + t475 * t340 + t477 * t341 - t453 * t430 + t472 * t432;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t500 - mrSges(2,2) * t495 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t440 - mrSges(3,2) * t441 + t481 * (mrSges(4,2) * t437 - mrSges(4,3) * t421 + Ifges(4,1) * t460 + Ifges(4,4) * t461 + Ifges(4,5) * qJDD(3) - pkin(8) * t339 - qJD(3) * t450 - t480 * t334 + t484 * t335) + t485 * (-mrSges(4,1) * t437 + mrSges(4,3) * t422 + Ifges(4,4) * t460 + Ifges(4,2) * t461 + Ifges(4,6) * qJDD(3) - pkin(3) * t491 + pkin(8) * t498 + qJD(3) * t451 + t484 * t334 + t480 * t335) + pkin(2) * t488 + pkin(7) * t499 + pkin(1) * (t476 * (m(3) * t441 - mrSges(3,1) * t487 - qJDD(1) * mrSges(3,2) + t499) + t478 * (m(3) * t440 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t487 + t488)); m(3) * t474 + t337 * t485 + t338 * t481; Ifges(4,6) * t461 + Ifges(4,5) * t460 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t421 - mrSges(4,2) * t422 + pkin(3) * t339 + t489 + (t450 * t481 - t451 * t485) * qJD(1); t489; t352; t490;];
tauJ  = t1;
