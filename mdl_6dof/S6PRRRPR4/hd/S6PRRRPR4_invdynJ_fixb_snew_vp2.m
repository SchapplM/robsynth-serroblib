% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPR4
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
% Datum: 2019-05-05 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:42:31
% EndTime: 2019-05-05 07:42:36
% DurationCPUTime: 4.54s
% Computational Cost: add. (50417->290), mult. (101590->371), div. (0->0), fcn. (74206->14), ass. (0->121)
t471 = sin(pkin(11));
t474 = cos(pkin(11));
t461 = t471 * g(1) - t474 * g(2);
t469 = -g(3) + qJDD(1);
t472 = sin(pkin(6));
t475 = cos(pkin(6));
t504 = t461 * t475 + t469 * t472;
t462 = -t474 * g(1) - t471 * g(2);
t479 = sin(qJ(2));
t483 = cos(qJ(2));
t422 = t483 * t462 + t504 * t479;
t485 = qJD(2) ^ 2;
t417 = -t485 * pkin(2) + qJDD(2) * pkin(8) + t422;
t441 = -t472 * t461 + t475 * t469;
t478 = sin(qJ(3));
t482 = cos(qJ(3));
t412 = t482 * t417 + t478 * t441;
t458 = (-t482 * pkin(3) - t478 * pkin(9)) * qJD(2);
t484 = qJD(3) ^ 2;
t498 = t482 * qJD(2);
t396 = -t484 * pkin(3) + qJDD(3) * pkin(9) + t458 * t498 + t412;
t421 = -t479 * t462 + t483 * t504;
t416 = -qJDD(2) * pkin(2) - t485 * pkin(8) - t421;
t497 = qJD(2) * qJD(3);
t496 = t482 * t497;
t459 = t478 * qJDD(2) + t496;
t468 = t478 * t497;
t460 = t482 * qJDD(2) - t468;
t401 = (-t459 - t496) * pkin(9) + (-t460 + t468) * pkin(3) + t416;
t477 = sin(qJ(4));
t481 = cos(qJ(4));
t385 = -t477 * t396 + t481 * t401;
t499 = t478 * qJD(2);
t455 = t481 * qJD(3) - t477 * t499;
t432 = t455 * qJD(4) + t477 * qJDD(3) + t481 * t459;
t452 = qJDD(4) - t460;
t456 = t477 * qJD(3) + t481 * t499;
t467 = qJD(4) - t498;
t375 = (t455 * t467 - t432) * qJ(5) + (t455 * t456 + t452) * pkin(4) + t385;
t386 = t481 * t396 + t477 * t401;
t431 = -t456 * qJD(4) + t481 * qJDD(3) - t477 * t459;
t439 = t467 * pkin(4) - t456 * qJ(5);
t451 = t455 ^ 2;
t377 = -t451 * pkin(4) + t431 * qJ(5) - t467 * t439 + t386;
t470 = sin(pkin(12));
t473 = cos(pkin(12));
t435 = t470 * t455 + t473 * t456;
t369 = -0.2e1 * qJD(5) * t435 + t473 * t375 - t470 * t377;
t407 = t470 * t431 + t473 * t432;
t434 = t473 * t455 - t470 * t456;
t367 = (t434 * t467 - t407) * pkin(10) + (t434 * t435 + t452) * pkin(5) + t369;
t370 = 0.2e1 * qJD(5) * t434 + t470 * t375 + t473 * t377;
t406 = t473 * t431 - t470 * t432;
t420 = t467 * pkin(5) - t435 * pkin(10);
t433 = t434 ^ 2;
t368 = -t433 * pkin(5) + t406 * pkin(10) - t467 * t420 + t370;
t476 = sin(qJ(6));
t480 = cos(qJ(6));
t365 = t480 * t367 - t476 * t368;
t409 = t480 * t434 - t476 * t435;
t383 = t409 * qJD(6) + t476 * t406 + t480 * t407;
t410 = t476 * t434 + t480 * t435;
t393 = -t409 * mrSges(7,1) + t410 * mrSges(7,2);
t466 = qJD(6) + t467;
t399 = -t466 * mrSges(7,2) + t409 * mrSges(7,3);
t448 = qJDD(6) + t452;
t360 = m(7) * t365 + t448 * mrSges(7,1) - t383 * mrSges(7,3) - t410 * t393 + t466 * t399;
t366 = t476 * t367 + t480 * t368;
t382 = -t410 * qJD(6) + t480 * t406 - t476 * t407;
t400 = t466 * mrSges(7,1) - t410 * mrSges(7,3);
t361 = m(7) * t366 - t448 * mrSges(7,2) + t382 * mrSges(7,3) + t409 * t393 - t466 * t400;
t354 = t480 * t360 + t476 * t361;
t413 = -t434 * mrSges(6,1) + t435 * mrSges(6,2);
t418 = -t467 * mrSges(6,2) + t434 * mrSges(6,3);
t352 = m(6) * t369 + t452 * mrSges(6,1) - t407 * mrSges(6,3) - t435 * t413 + t467 * t418 + t354;
t419 = t467 * mrSges(6,1) - t435 * mrSges(6,3);
t492 = -t476 * t360 + t480 * t361;
t353 = m(6) * t370 - t452 * mrSges(6,2) + t406 * mrSges(6,3) + t434 * t413 - t467 * t419 + t492;
t348 = t473 * t352 + t470 * t353;
t404 = Ifges(6,4) * t435 + Ifges(6,2) * t434 + Ifges(6,6) * t467;
t405 = Ifges(6,1) * t435 + Ifges(6,4) * t434 + Ifges(6,5) * t467;
t424 = Ifges(5,4) * t456 + Ifges(5,2) * t455 + Ifges(5,6) * t467;
t425 = Ifges(5,1) * t456 + Ifges(5,4) * t455 + Ifges(5,5) * t467;
t389 = Ifges(7,4) * t410 + Ifges(7,2) * t409 + Ifges(7,6) * t466;
t390 = Ifges(7,1) * t410 + Ifges(7,4) * t409 + Ifges(7,5) * t466;
t489 = -mrSges(7,1) * t365 + mrSges(7,2) * t366 - Ifges(7,5) * t383 - Ifges(7,6) * t382 - Ifges(7,3) * t448 - t410 * t389 + t409 * t390;
t503 = mrSges(5,1) * t385 + mrSges(6,1) * t369 - mrSges(5,2) * t386 - mrSges(6,2) * t370 + Ifges(5,5) * t432 + Ifges(6,5) * t407 + Ifges(5,6) * t431 + Ifges(6,6) * t406 + pkin(4) * t348 + pkin(5) * t354 + t435 * t404 - t434 * t405 + t456 * t424 - t455 * t425 + (Ifges(5,3) + Ifges(6,3)) * t452 - t489;
t457 = (-t482 * mrSges(4,1) + t478 * mrSges(4,2)) * qJD(2);
t463 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t499;
t436 = -t455 * mrSges(5,1) + t456 * mrSges(5,2);
t438 = -t467 * mrSges(5,2) + t455 * mrSges(5,3);
t346 = m(5) * t385 + t452 * mrSges(5,1) - t432 * mrSges(5,3) - t456 * t436 + t467 * t438 + t348;
t440 = t467 * mrSges(5,1) - t456 * mrSges(5,3);
t493 = -t470 * t352 + t473 * t353;
t347 = m(5) * t386 - t452 * mrSges(5,2) + t431 * mrSges(5,3) + t455 * t436 - t467 * t440 + t493;
t494 = -t477 * t346 + t481 * t347;
t341 = m(4) * t412 - qJDD(3) * mrSges(4,2) + t460 * mrSges(4,3) - qJD(3) * t463 + t457 * t498 + t494;
t411 = -t478 * t417 + t482 * t441;
t464 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t498;
t395 = -qJDD(3) * pkin(3) - t484 * pkin(9) + t458 * t499 - t411;
t387 = -t431 * pkin(4) - t451 * qJ(5) + t456 * t439 + qJDD(5) + t395;
t372 = -t406 * pkin(5) - t433 * pkin(10) + t435 * t420 + t387;
t491 = m(7) * t372 - t382 * mrSges(7,1) + t383 * mrSges(7,2) - t409 * t399 + t410 * t400;
t363 = m(6) * t387 - t406 * mrSges(6,1) + t407 * mrSges(6,2) - t434 * t418 + t435 * t419 + t491;
t487 = -m(5) * t395 + t431 * mrSges(5,1) - t432 * mrSges(5,2) + t455 * t438 - t456 * t440 - t363;
t362 = m(4) * t411 + qJDD(3) * mrSges(4,1) - t459 * mrSges(4,3) + qJD(3) * t464 - t457 * t499 + t487;
t495 = t482 * t341 - t478 * t362;
t342 = t481 * t346 + t477 * t347;
t488 = -m(4) * t416 + t460 * mrSges(4,1) - t459 * mrSges(4,2) - t463 * t499 + t464 * t498 - t342;
t447 = Ifges(4,5) * qJD(3) + (t478 * Ifges(4,1) + t482 * Ifges(4,4)) * qJD(2);
t446 = Ifges(4,6) * qJD(3) + (t478 * Ifges(4,4) + t482 * Ifges(4,2)) * qJD(2);
t423 = Ifges(5,5) * t456 + Ifges(5,6) * t455 + Ifges(5,3) * t467;
t403 = Ifges(6,5) * t435 + Ifges(6,6) * t434 + Ifges(6,3) * t467;
t388 = Ifges(7,5) * t410 + Ifges(7,6) * t409 + Ifges(7,3) * t466;
t356 = mrSges(7,2) * t372 - mrSges(7,3) * t365 + Ifges(7,1) * t383 + Ifges(7,4) * t382 + Ifges(7,5) * t448 + t409 * t388 - t466 * t389;
t355 = -mrSges(7,1) * t372 + mrSges(7,3) * t366 + Ifges(7,4) * t383 + Ifges(7,2) * t382 + Ifges(7,6) * t448 - t410 * t388 + t466 * t390;
t344 = mrSges(6,2) * t387 - mrSges(6,3) * t369 + Ifges(6,1) * t407 + Ifges(6,4) * t406 + Ifges(6,5) * t452 - pkin(10) * t354 - t476 * t355 + t480 * t356 + t434 * t403 - t467 * t404;
t343 = -mrSges(6,1) * t387 + mrSges(6,3) * t370 + Ifges(6,4) * t407 + Ifges(6,2) * t406 + Ifges(6,6) * t452 - pkin(5) * t491 + pkin(10) * t492 + t480 * t355 + t476 * t356 - t435 * t403 + t467 * t405;
t339 = mrSges(5,2) * t395 - mrSges(5,3) * t385 + Ifges(5,1) * t432 + Ifges(5,4) * t431 + Ifges(5,5) * t452 - qJ(5) * t348 - t470 * t343 + t473 * t344 + t455 * t423 - t467 * t424;
t338 = -mrSges(5,1) * t395 + mrSges(5,3) * t386 + Ifges(5,4) * t432 + Ifges(5,2) * t431 + Ifges(5,6) * t452 - pkin(4) * t363 + qJ(5) * t493 + t473 * t343 + t470 * t344 - t456 * t423 + t467 * t425;
t1 = [m(2) * t469 + t475 * (m(3) * t441 + t478 * t341 + t482 * t362) + (t479 * (m(3) * t422 - t485 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t495) + t483 * (m(3) * t421 + qJDD(2) * mrSges(3,1) - t485 * mrSges(3,2) + t488)) * t472; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t421 - mrSges(3,2) * t422 + t478 * (mrSges(4,2) * t416 - mrSges(4,3) * t411 + Ifges(4,1) * t459 + Ifges(4,4) * t460 + Ifges(4,5) * qJDD(3) - pkin(9) * t342 - qJD(3) * t446 - t477 * t338 + t481 * t339) + t482 * (-mrSges(4,1) * t416 + mrSges(4,3) * t412 + Ifges(4,4) * t459 + Ifges(4,2) * t460 + Ifges(4,6) * qJDD(3) - pkin(3) * t342 + qJD(3) * t447 - t503) + pkin(2) * t488 + pkin(8) * t495; Ifges(4,5) * t459 + Ifges(4,6) * t460 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t411 - mrSges(4,2) * t412 + t477 * t339 + t481 * t338 + pkin(3) * t487 + pkin(9) * t494 + (t478 * t446 - t482 * t447) * qJD(2); t503; t363; -t489;];
tauJ  = t1;
