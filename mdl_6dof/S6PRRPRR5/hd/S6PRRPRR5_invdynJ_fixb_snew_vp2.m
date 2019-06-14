% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 05:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:23:00
% EndTime: 2019-05-05 05:23:05
% DurationCPUTime: 4.36s
% Computational Cost: add. (46278->290), mult. (96670->371), div. (0->0), fcn. (70818->14), ass. (0->121)
t471 = sin(pkin(11));
t474 = cos(pkin(11));
t461 = g(1) * t471 - g(2) * t474;
t469 = -g(3) + qJDD(1);
t472 = sin(pkin(6));
t475 = cos(pkin(6));
t502 = t461 * t475 + t469 * t472;
t462 = -g(1) * t474 - g(2) * t471;
t479 = sin(qJ(2));
t483 = cos(qJ(2));
t421 = -t462 * t479 + t483 * t502;
t422 = t483 * t462 + t479 * t502;
t485 = qJD(2) ^ 2;
t417 = -pkin(2) * t485 + qJDD(2) * pkin(8) + t422;
t439 = -t461 * t472 + t469 * t475;
t478 = sin(qJ(3));
t482 = cos(qJ(3));
t412 = t417 * t482 + t439 * t478;
t457 = (-pkin(3) * t482 - qJ(4) * t478) * qJD(2);
t484 = qJD(3) ^ 2;
t498 = t482 * qJD(2);
t396 = -pkin(3) * t484 + qJDD(3) * qJ(4) + t457 * t498 + t412;
t416 = -qJDD(2) * pkin(2) - t485 * pkin(8) - t421;
t497 = qJD(2) * qJD(3);
t496 = t482 * t497;
t459 = qJDD(2) * t478 + t496;
t468 = t478 * t497;
t460 = qJDD(2) * t482 - t468;
t402 = (-t459 - t496) * qJ(4) + (-t460 + t468) * pkin(3) + t416;
t470 = sin(pkin(12));
t473 = cos(pkin(12));
t499 = t478 * qJD(2);
t452 = qJD(3) * t470 + t473 * t499;
t385 = -0.2e1 * qJD(4) * t452 - t470 * t396 + t402 * t473;
t437 = qJDD(3) * t470 + t459 * t473;
t451 = qJD(3) * t473 - t470 * t499;
t375 = (-t451 * t498 - t437) * pkin(9) + (t451 * t452 - t460) * pkin(4) + t385;
t386 = 0.2e1 * qJD(4) * t451 + t396 * t473 + t402 * t470;
t436 = qJDD(3) * t473 - t459 * t470;
t438 = -pkin(4) * t498 - pkin(9) * t452;
t450 = t451 ^ 2;
t377 = -pkin(4) * t450 + pkin(9) * t436 + t438 * t498 + t386;
t477 = sin(qJ(5));
t481 = cos(qJ(5));
t369 = t375 * t481 - t477 * t377;
t429 = t451 * t481 - t452 * t477;
t404 = qJD(5) * t429 + t436 * t477 + t437 * t481;
t430 = t451 * t477 + t452 * t481;
t454 = qJDD(5) - t460;
t467 = qJD(5) - t498;
t367 = (t429 * t467 - t404) * pkin(10) + (t429 * t430 + t454) * pkin(5) + t369;
t370 = t375 * t477 + t377 * t481;
t403 = -qJD(5) * t430 + t436 * t481 - t437 * t477;
t420 = pkin(5) * t467 - pkin(10) * t430;
t428 = t429 ^ 2;
t368 = -pkin(5) * t428 + pkin(10) * t403 - t420 * t467 + t370;
t476 = sin(qJ(6));
t480 = cos(qJ(6));
t365 = t367 * t480 - t368 * t476;
t409 = t429 * t480 - t430 * t476;
t383 = qJD(6) * t409 + t403 * t476 + t404 * t480;
t410 = t429 * t476 + t430 * t480;
t393 = -mrSges(7,1) * t409 + mrSges(7,2) * t410;
t466 = qJD(6) + t467;
t399 = -mrSges(7,2) * t466 + mrSges(7,3) * t409;
t448 = qJDD(6) + t454;
t360 = m(7) * t365 + mrSges(7,1) * t448 - mrSges(7,3) * t383 - t393 * t410 + t399 * t466;
t366 = t367 * t476 + t368 * t480;
t382 = -qJD(6) * t410 + t403 * t480 - t404 * t476;
t400 = mrSges(7,1) * t466 - mrSges(7,3) * t410;
t361 = m(7) * t366 - mrSges(7,2) * t448 + mrSges(7,3) * t382 + t393 * t409 - t400 * t466;
t354 = t360 * t480 + t361 * t476;
t413 = -mrSges(6,1) * t429 + mrSges(6,2) * t430;
t418 = -mrSges(6,2) * t467 + mrSges(6,3) * t429;
t352 = m(6) * t369 + mrSges(6,1) * t454 - mrSges(6,3) * t404 - t413 * t430 + t418 * t467 + t354;
t419 = mrSges(6,1) * t467 - mrSges(6,3) * t430;
t492 = -t360 * t476 + t361 * t480;
t353 = m(6) * t370 - mrSges(6,2) * t454 + mrSges(6,3) * t403 + t413 * t429 - t419 * t467 + t492;
t348 = t352 * t481 + t353 * t477;
t458 = (-mrSges(4,1) * t482 + mrSges(4,2) * t478) * qJD(2);
t463 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t499;
t431 = -mrSges(5,1) * t451 + mrSges(5,2) * t452;
t434 = mrSges(5,2) * t498 + mrSges(5,3) * t451;
t346 = m(5) * t385 - mrSges(5,1) * t460 - mrSges(5,3) * t437 - t431 * t452 - t434 * t498 + t348;
t435 = -mrSges(5,1) * t498 - mrSges(5,3) * t452;
t493 = -t352 * t477 + t353 * t481;
t347 = m(5) * t386 + mrSges(5,2) * t460 + mrSges(5,3) * t436 + t431 * t451 + t435 * t498 + t493;
t494 = -t346 * t470 + t347 * t473;
t341 = m(4) * t412 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t460 - qJD(3) * t463 + t458 * t498 + t494;
t411 = -t417 * t478 + t439 * t482;
t395 = -qJDD(3) * pkin(3) - qJ(4) * t484 + t457 * t499 + qJDD(4) - t411;
t387 = -pkin(4) * t436 - pkin(9) * t450 + t438 * t452 + t395;
t372 = -pkin(5) * t403 - pkin(10) * t428 + t420 * t430 + t387;
t491 = m(7) * t372 - mrSges(7,1) * t382 + mrSges(7,2) * t383 - t399 * t409 + t400 * t410;
t488 = m(6) * t387 - mrSges(6,1) * t403 + mrSges(6,2) * t404 - t418 * t429 + t419 * t430 + t491;
t363 = m(5) * t395 - mrSges(5,1) * t436 + mrSges(5,2) * t437 - t434 * t451 + t435 * t452 + t488;
t464 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t498;
t362 = m(4) * t411 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t459 + qJD(3) * t464 - t458 * t499 - t363;
t495 = t341 * t482 - t362 * t478;
t342 = t346 * t473 + t347 * t470;
t389 = Ifges(7,4) * t410 + Ifges(7,2) * t409 + Ifges(7,6) * t466;
t390 = Ifges(7,1) * t410 + Ifges(7,4) * t409 + Ifges(7,5) * t466;
t489 = -mrSges(7,1) * t365 + mrSges(7,2) * t366 - Ifges(7,5) * t383 - Ifges(7,6) * t382 - Ifges(7,3) * t448 - t389 * t410 + t409 * t390;
t487 = -m(4) * t416 + mrSges(4,1) * t460 - mrSges(4,2) * t459 - t463 * t499 + t464 * t498 - t342;
t406 = Ifges(6,4) * t430 + Ifges(6,2) * t429 + Ifges(6,6) * t467;
t407 = Ifges(6,1) * t430 + Ifges(6,4) * t429 + Ifges(6,5) * t467;
t486 = mrSges(6,1) * t369 - mrSges(6,2) * t370 + Ifges(6,5) * t404 + Ifges(6,6) * t403 + Ifges(6,3) * t454 + pkin(5) * t354 + t430 * t406 - t429 * t407 - t489;
t447 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t478 + Ifges(4,4) * t482) * qJD(2);
t446 = Ifges(4,6) * qJD(3) + (t478 * Ifges(4,4) + Ifges(4,2) * t482) * qJD(2);
t425 = Ifges(5,1) * t452 + Ifges(5,4) * t451 - Ifges(5,5) * t498;
t424 = Ifges(5,4) * t452 + Ifges(5,2) * t451 - Ifges(5,6) * t498;
t423 = Ifges(5,5) * t452 + Ifges(5,6) * t451 - Ifges(5,3) * t498;
t405 = Ifges(6,5) * t430 + Ifges(6,6) * t429 + Ifges(6,3) * t467;
t388 = Ifges(7,5) * t410 + Ifges(7,6) * t409 + Ifges(7,3) * t466;
t356 = mrSges(7,2) * t372 - mrSges(7,3) * t365 + Ifges(7,1) * t383 + Ifges(7,4) * t382 + Ifges(7,5) * t448 + t388 * t409 - t389 * t466;
t355 = -mrSges(7,1) * t372 + mrSges(7,3) * t366 + Ifges(7,4) * t383 + Ifges(7,2) * t382 + Ifges(7,6) * t448 - t388 * t410 + t390 * t466;
t344 = mrSges(6,2) * t387 - mrSges(6,3) * t369 + Ifges(6,1) * t404 + Ifges(6,4) * t403 + Ifges(6,5) * t454 - pkin(10) * t354 - t355 * t476 + t356 * t480 + t405 * t429 - t406 * t467;
t343 = -mrSges(6,1) * t387 + mrSges(6,3) * t370 + Ifges(6,4) * t404 + Ifges(6,2) * t403 + Ifges(6,6) * t454 - pkin(5) * t491 + pkin(10) * t492 + t480 * t355 + t476 * t356 - t430 * t405 + t467 * t407;
t339 = mrSges(5,2) * t395 - mrSges(5,3) * t385 + Ifges(5,1) * t437 + Ifges(5,4) * t436 - Ifges(5,5) * t460 - pkin(9) * t348 - t343 * t477 + t344 * t481 + t423 * t451 + t424 * t498;
t338 = -mrSges(5,1) * t395 + mrSges(5,3) * t386 + Ifges(5,4) * t437 + Ifges(5,2) * t436 - Ifges(5,6) * t460 - pkin(4) * t488 + pkin(9) * t493 + t481 * t343 + t477 * t344 - t452 * t423 - t425 * t498;
t1 = [m(2) * t469 + t475 * (m(3) * t439 + t341 * t478 + t362 * t482) + (t479 * (m(3) * t422 - mrSges(3,1) * t485 - qJDD(2) * mrSges(3,2) + t495) + t483 * (m(3) * t421 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t485 + t487)) * t472; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t421 - mrSges(3,2) * t422 + t478 * (mrSges(4,2) * t416 - mrSges(4,3) * t411 + Ifges(4,1) * t459 + Ifges(4,4) * t460 + Ifges(4,5) * qJDD(3) - qJ(4) * t342 - qJD(3) * t446 - t338 * t470 + t339 * t473) + t482 * (-mrSges(4,1) * t416 - mrSges(5,1) * t385 + mrSges(5,2) * t386 + mrSges(4,3) * t412 + Ifges(4,4) * t459 - Ifges(5,5) * t437 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t436 - pkin(3) * t342 - pkin(4) * t348 + qJD(3) * t447 - t452 * t424 + t451 * t425 - t486 + (Ifges(4,2) + Ifges(5,3)) * t460) + pkin(2) * t487 + pkin(8) * t495; Ifges(4,5) * t459 + Ifges(4,6) * t460 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t411 - mrSges(4,2) * t412 + t470 * t339 + t473 * t338 - pkin(3) * t363 + qJ(4) * t494 + (t446 * t478 - t447 * t482) * qJD(2); t363; t486; -t489;];
tauJ  = t1;
