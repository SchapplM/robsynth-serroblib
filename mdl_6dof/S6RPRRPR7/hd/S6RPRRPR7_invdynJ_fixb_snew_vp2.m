% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:01:06
% EndTime: 2019-05-05 23:01:10
% DurationCPUTime: 3.53s
% Computational Cost: add. (36194->290), mult. (78078->368), div. (0->0), fcn. (54431->10), ass. (0->116)
t476 = sin(qJ(1));
t480 = cos(qJ(1));
t491 = -g(1) * t480 - g(2) * t476;
t487 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t491;
t502 = 2 * qJD(5);
t501 = -pkin(1) - pkin(7);
t481 = qJD(1) ^ 2;
t495 = t476 * g(1) - t480 * g(2);
t486 = -qJ(2) * t481 + qJDD(2) - t495;
t445 = qJDD(1) * t501 + t486;
t475 = sin(qJ(3));
t479 = cos(qJ(3));
t437 = t475 * g(3) + t479 * t445;
t498 = qJD(1) * qJD(3);
t496 = t475 * t498;
t457 = qJDD(1) * t479 - t496;
t413 = (-t457 - t496) * pkin(8) + (-t475 * t479 * t481 + qJDD(3)) * pkin(3) + t437;
t438 = -g(3) * t479 + t475 * t445;
t456 = -qJDD(1) * t475 - t479 * t498;
t499 = qJD(1) * t479;
t460 = qJD(3) * pkin(3) - pkin(8) * t499;
t470 = t475 ^ 2;
t414 = -pkin(3) * t470 * t481 + pkin(8) * t456 - qJD(3) * t460 + t438;
t474 = sin(qJ(4));
t478 = cos(qJ(4));
t393 = t478 * t413 - t414 * t474;
t452 = (-t479 * t474 - t475 * t478) * qJD(1);
t424 = qJD(4) * t452 + t456 * t474 + t457 * t478;
t453 = (-t475 * t474 + t479 * t478) * qJD(1);
t467 = qJDD(3) + qJDD(4);
t468 = qJD(3) + qJD(4);
t379 = (t452 * t468 - t424) * qJ(5) + (t452 * t453 + t467) * pkin(4) + t393;
t394 = t474 * t413 + t478 * t414;
t423 = -qJD(4) * t453 + t456 * t478 - t457 * t474;
t442 = pkin(4) * t468 - qJ(5) * t453;
t448 = t452 ^ 2;
t381 = -pkin(4) * t448 + qJ(5) * t423 - t442 * t468 + t394;
t471 = sin(pkin(10));
t472 = cos(pkin(10));
t434 = t452 * t472 - t453 * t471;
t376 = t471 * t379 + t472 * t381 + t434 * t502;
t399 = t423 * t472 - t424 * t471;
t435 = t452 * t471 + t453 * t472;
t408 = -mrSges(6,1) * t434 + mrSges(6,2) * t435;
t426 = mrSges(6,1) * t468 - mrSges(6,3) * t435;
t409 = -pkin(5) * t434 - pkin(9) * t435;
t466 = t468 ^ 2;
t373 = -pkin(5) * t466 + pkin(9) * t467 + t409 * t434 + t376;
t417 = -pkin(3) * t456 + t460 * t499 + (-pkin(8) * t470 + t501) * t481 + t487;
t386 = -pkin(4) * t423 - qJ(5) * t448 + t453 * t442 + qJDD(5) + t417;
t400 = t423 * t471 + t424 * t472;
t377 = (-t434 * t468 - t400) * pkin(9) + (t435 * t468 - t399) * pkin(5) + t386;
t473 = sin(qJ(6));
t477 = cos(qJ(6));
t370 = -t373 * t473 + t377 * t477;
t421 = -t435 * t473 + t468 * t477;
t384 = qJD(6) * t421 + t400 * t477 + t467 * t473;
t398 = qJDD(6) - t399;
t422 = t435 * t477 + t468 * t473;
t401 = -mrSges(7,1) * t421 + mrSges(7,2) * t422;
t428 = qJD(6) - t434;
t402 = -mrSges(7,2) * t428 + mrSges(7,3) * t421;
t367 = m(7) * t370 + mrSges(7,1) * t398 - mrSges(7,3) * t384 - t401 * t422 + t402 * t428;
t371 = t373 * t477 + t377 * t473;
t383 = -qJD(6) * t422 - t400 * t473 + t467 * t477;
t403 = mrSges(7,1) * t428 - mrSges(7,3) * t422;
t368 = m(7) * t371 - mrSges(7,2) * t398 + mrSges(7,3) * t383 + t401 * t421 - t403 * t428;
t492 = -t367 * t473 + t477 * t368;
t354 = m(6) * t376 - mrSges(6,2) * t467 + mrSges(6,3) * t399 + t408 * t434 - t426 * t468 + t492;
t489 = -t379 * t472 + t381 * t471;
t375 = -0.2e1 * qJD(5) * t435 - t489;
t425 = -mrSges(6,2) * t468 + mrSges(6,3) * t434;
t372 = -pkin(5) * t467 - pkin(9) * t466 + (t502 + t409) * t435 + t489;
t485 = -m(7) * t372 + t383 * mrSges(7,1) - mrSges(7,2) * t384 + t421 * t402 - t403 * t422;
t363 = m(6) * t375 + mrSges(6,1) * t467 - mrSges(6,3) * t400 - t408 * t435 + t425 * t468 + t485;
t351 = t471 * t354 + t472 * t363;
t436 = -mrSges(5,1) * t452 + mrSges(5,2) * t453;
t441 = -mrSges(5,2) * t468 + mrSges(5,3) * t452;
t348 = m(5) * t393 + mrSges(5,1) * t467 - mrSges(5,3) * t424 - t436 * t453 + t441 * t468 + t351;
t443 = mrSges(5,1) * t468 - mrSges(5,3) * t453;
t493 = t472 * t354 - t363 * t471;
t349 = m(5) * t394 - mrSges(5,2) * t467 + mrSges(5,3) * t423 + t436 * t452 - t443 * t468 + t493;
t342 = t478 * t348 + t474 * t349;
t357 = t477 * t367 + t473 * t368;
t500 = qJD(1) * t475;
t494 = -t348 * t474 + t478 * t349;
t455 = (t475 * mrSges(4,1) + t479 * mrSges(4,2)) * qJD(1);
t458 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t500;
t459 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t499;
t490 = (m(4) * t437 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t457 + qJD(3) * t458 - t455 * t499 + t342) * t479 + t475 * (m(4) * t438 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t456 - qJD(3) * t459 - t455 * t500 + t494);
t355 = m(6) * t386 - t399 * mrSges(6,1) + t400 * mrSges(6,2) - t434 * t425 + t435 * t426 + t357;
t388 = Ifges(7,4) * t422 + Ifges(7,2) * t421 + Ifges(7,6) * t428;
t389 = Ifges(7,1) * t422 + Ifges(7,4) * t421 + Ifges(7,5) * t428;
t484 = mrSges(7,1) * t370 - mrSges(7,2) * t371 + Ifges(7,5) * t384 + Ifges(7,6) * t383 + Ifges(7,3) * t398 + t388 * t422 - t389 * t421;
t483 = m(5) * t417 - t423 * mrSges(5,1) + t424 * mrSges(5,2) - t452 * t441 + t453 * t443 + t355;
t387 = Ifges(7,5) * t422 + Ifges(7,6) * t421 + Ifges(7,3) * t428;
t360 = -mrSges(7,1) * t372 + mrSges(7,3) * t371 + Ifges(7,4) * t384 + Ifges(7,2) * t383 + Ifges(7,6) * t398 - t387 * t422 + t389 * t428;
t361 = mrSges(7,2) * t372 - mrSges(7,3) * t370 + Ifges(7,1) * t384 + Ifges(7,4) * t383 + Ifges(7,5) * t398 + t387 * t421 - t388 * t428;
t405 = Ifges(6,4) * t435 + Ifges(6,2) * t434 + Ifges(6,6) * t468;
t406 = Ifges(6,1) * t435 + Ifges(6,4) * t434 + Ifges(6,5) * t468;
t430 = Ifges(5,4) * t453 + Ifges(5,2) * t452 + Ifges(5,6) * t468;
t431 = Ifges(5,1) * t453 + Ifges(5,4) * t452 + Ifges(5,5) * t468;
t482 = mrSges(5,1) * t393 + mrSges(6,1) * t375 - mrSges(5,2) * t394 - mrSges(6,2) * t376 + pkin(4) * t351 + pkin(5) * t485 + pkin(9) * t492 + t477 * t360 + t473 * t361 + t435 * t405 - t434 * t406 - t452 * t431 + Ifges(6,6) * t399 + Ifges(6,5) * t400 + t453 * t430 + Ifges(5,6) * t423 + Ifges(5,5) * t424 + (Ifges(6,3) + Ifges(5,3)) * t467;
t451 = (Ifges(4,5) * qJD(3)) + (t479 * Ifges(4,1) - t475 * Ifges(4,4)) * qJD(1);
t450 = (Ifges(4,6) * qJD(3)) + (t479 * Ifges(4,4) - t475 * Ifges(4,2)) * qJD(1);
t447 = -qJDD(1) * pkin(1) + t486;
t446 = pkin(1) * t481 - t487;
t444 = t481 * t501 + t487;
t429 = Ifges(5,5) * t453 + Ifges(5,6) * t452 + Ifges(5,3) * t468;
t404 = Ifges(6,5) * t435 + Ifges(6,6) * t434 + Ifges(6,3) * t468;
t344 = -mrSges(6,1) * t386 + mrSges(6,3) * t376 + Ifges(6,4) * t400 + Ifges(6,2) * t399 + Ifges(6,6) * t467 - pkin(5) * t357 - t404 * t435 + t406 * t468 - t484;
t343 = mrSges(6,2) * t386 - mrSges(6,3) * t375 + Ifges(6,1) * t400 + Ifges(6,4) * t399 + Ifges(6,5) * t467 - pkin(9) * t357 - t360 * t473 + t361 * t477 + t404 * t434 - t405 * t468;
t339 = mrSges(5,2) * t417 - mrSges(5,3) * t393 + Ifges(5,1) * t424 + Ifges(5,4) * t423 + Ifges(5,5) * t467 - qJ(5) * t351 + t343 * t472 - t344 * t471 + t429 * t452 - t430 * t468;
t338 = m(3) * t447 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t481) + t490;
t337 = -mrSges(5,1) * t417 + mrSges(5,3) * t394 + Ifges(5,4) * t424 + Ifges(5,2) * t423 + Ifges(5,6) * t467 - pkin(4) * t355 + qJ(5) * t493 + t471 * t343 + t472 * t344 - t453 * t429 + t468 * t431;
t1 = [mrSges(2,1) * t495 - mrSges(2,2) * t491 + mrSges(3,2) * t447 - mrSges(3,3) * t446 + t479 * (mrSges(4,2) * t444 - mrSges(4,3) * t437 + Ifges(4,1) * t457 + Ifges(4,4) * t456 + Ifges(4,5) * qJDD(3) - pkin(8) * t342 - qJD(3) * t450 - t337 * t474 + t339 * t478) - t475 * (-mrSges(4,1) * t444 + mrSges(4,3) * t438 + Ifges(4,4) * t457 + Ifges(4,2) * t456 + Ifges(4,6) * qJDD(3) - pkin(3) * t483 + pkin(8) * t494 + qJD(3) * t451 + t478 * t337 + t474 * t339) - pkin(7) * t490 - pkin(1) * t338 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t446 + m(4) * t444 - t456 * mrSges(4,1) + t481 * mrSges(3,2) + t457 * mrSges(4,2) + t483 + qJDD(1) * mrSges(3,3) + (t458 * t475 + t459 * t479) * qJD(1)) * qJ(2); t338; t482 + Ifges(4,3) * qJDD(3) + Ifges(4,6) * t456 + Ifges(4,5) * t457 + mrSges(4,1) * t437 - mrSges(4,2) * t438 + (t479 * t450 + t475 * t451) * qJD(1) + pkin(3) * t342; t482; t355; t484;];
tauJ  = t1;
