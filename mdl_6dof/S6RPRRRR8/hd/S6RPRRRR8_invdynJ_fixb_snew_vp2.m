% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:14:52
% EndTime: 2019-05-06 04:14:57
% DurationCPUTime: 3.66s
% Computational Cost: add. (40951->291), mult. (80134->366), div. (0->0), fcn. (55034->10), ass. (0->117)
t484 = sin(qJ(1));
t489 = cos(qJ(1));
t501 = -t489 * g(1) - t484 * g(2);
t498 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t501;
t511 = -pkin(1) - pkin(7);
t490 = qJD(1) ^ 2;
t505 = t484 * g(1) - t489 * g(2);
t497 = -t490 * qJ(2) + qJDD(2) - t505;
t450 = t511 * qJDD(1) + t497;
t483 = sin(qJ(3));
t488 = cos(qJ(3));
t440 = t483 * g(3) + t488 * t450;
t508 = qJD(1) * qJD(3);
t506 = t483 * t508;
t465 = qJDD(1) * t488 - t506;
t414 = (-t465 - t506) * pkin(8) + (-t483 * t488 * t490 + qJDD(3)) * pkin(3) + t440;
t441 = -g(3) * t488 + t483 * t450;
t464 = -qJDD(1) * t483 - t488 * t508;
t509 = qJD(1) * t488;
t468 = qJD(3) * pkin(3) - pkin(8) * t509;
t479 = t483 ^ 2;
t417 = -pkin(3) * t479 * t490 + pkin(8) * t464 - qJD(3) * t468 + t441;
t482 = sin(qJ(4));
t487 = cos(qJ(4));
t400 = t482 * t414 + t487 * t417;
t460 = (-t482 * t483 + t487 * t488) * qJD(1);
t429 = -t460 * qJD(4) + t464 * t487 - t482 * t465;
t510 = qJD(1) * t483;
t459 = -t482 * t509 - t487 * t510;
t438 = -mrSges(5,1) * t459 + mrSges(5,2) * t460;
t477 = qJD(3) + qJD(4);
t448 = mrSges(5,1) * t477 - mrSges(5,3) * t460;
t476 = qJDD(3) + qJDD(4);
t421 = -t464 * pkin(3) + t468 * t509 + (-pkin(8) * t479 + t511) * t490 + t498;
t430 = qJD(4) * t459 + t464 * t482 + t465 * t487;
t390 = (-t459 * t477 - t430) * pkin(9) + (t460 * t477 - t429) * pkin(4) + t421;
t439 = -pkin(4) * t459 - pkin(9) * t460;
t475 = t477 ^ 2;
t393 = -pkin(4) * t475 + pkin(9) * t476 + t439 * t459 + t400;
t481 = sin(qJ(5));
t486 = cos(qJ(5));
t379 = t486 * t390 - t481 * t393;
t443 = -t460 * t481 + t477 * t486;
t404 = qJD(5) * t443 + t430 * t486 + t476 * t481;
t428 = qJDD(5) - t429;
t444 = t460 * t486 + t477 * t481;
t455 = qJD(5) - t459;
t377 = (t443 * t455 - t404) * pkin(10) + (t443 * t444 + t428) * pkin(5) + t379;
t380 = t481 * t390 + t486 * t393;
t403 = -qJD(5) * t444 - t430 * t481 + t476 * t486;
t433 = pkin(5) * t455 - pkin(10) * t444;
t442 = t443 ^ 2;
t378 = -pkin(5) * t442 + pkin(10) * t403 - t433 * t455 + t380;
t480 = sin(qJ(6));
t485 = cos(qJ(6));
t375 = t377 * t485 - t378 * t480;
t415 = t443 * t485 - t444 * t480;
t386 = qJD(6) * t415 + t403 * t480 + t404 * t485;
t416 = t443 * t480 + t444 * t485;
t401 = -mrSges(7,1) * t415 + mrSges(7,2) * t416;
t452 = qJD(6) + t455;
t405 = -mrSges(7,2) * t452 + mrSges(7,3) * t415;
t423 = qJDD(6) + t428;
t371 = m(7) * t375 + mrSges(7,1) * t423 - mrSges(7,3) * t386 - t401 * t416 + t405 * t452;
t376 = t377 * t480 + t378 * t485;
t385 = -qJD(6) * t416 + t403 * t485 - t404 * t480;
t406 = mrSges(7,1) * t452 - mrSges(7,3) * t416;
t372 = m(7) * t376 - mrSges(7,2) * t423 + mrSges(7,3) * t385 + t401 * t415 - t406 * t452;
t363 = t485 * t371 + t480 * t372;
t418 = -mrSges(6,1) * t443 + mrSges(6,2) * t444;
t431 = -mrSges(6,2) * t455 + mrSges(6,3) * t443;
t361 = m(6) * t379 + mrSges(6,1) * t428 - mrSges(6,3) * t404 - t418 * t444 + t431 * t455 + t363;
t432 = mrSges(6,1) * t455 - mrSges(6,3) * t444;
t502 = -t371 * t480 + t485 * t372;
t362 = m(6) * t380 - mrSges(6,2) * t428 + mrSges(6,3) * t403 + t418 * t443 - t432 * t455 + t502;
t503 = -t361 * t481 + t486 * t362;
t355 = m(5) * t400 - mrSges(5,2) * t476 + mrSges(5,3) * t429 + t438 * t459 - t448 * t477 + t503;
t399 = t414 * t487 - t482 * t417;
t447 = -mrSges(5,2) * t477 + mrSges(5,3) * t459;
t392 = -pkin(4) * t476 - pkin(9) * t475 + t460 * t439 - t399;
t381 = -pkin(5) * t403 - pkin(10) * t442 + t433 * t444 + t392;
t496 = m(7) * t381 - t385 * mrSges(7,1) + mrSges(7,2) * t386 - t415 * t405 + t406 * t416;
t492 = -m(6) * t392 + t403 * mrSges(6,1) - mrSges(6,2) * t404 + t443 * t431 - t432 * t444 - t496;
t367 = m(5) * t399 + mrSges(5,1) * t476 - mrSges(5,3) * t430 - t438 * t460 + t447 * t477 + t492;
t350 = t482 * t355 + t487 * t367;
t357 = t486 * t361 + t481 * t362;
t504 = t487 * t355 - t367 * t482;
t463 = (mrSges(4,1) * t483 + mrSges(4,2) * t488) * qJD(1);
t466 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t510;
t467 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t509;
t500 = (m(4) * t440 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t465 + qJD(3) * t466 - t463 * t509 + t350) * t488 + (m(4) * t441 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t464 - qJD(3) * t467 - t463 * t510 + t504) * t483;
t495 = m(5) * t421 - t429 * mrSges(5,1) + t430 * mrSges(5,2) - t459 * t447 + t460 * t448 + t357;
t395 = Ifges(7,4) * t416 + Ifges(7,2) * t415 + Ifges(7,6) * t452;
t396 = Ifges(7,1) * t416 + Ifges(7,4) * t415 + Ifges(7,5) * t452;
t494 = -mrSges(7,1) * t375 + mrSges(7,2) * t376 - Ifges(7,5) * t386 - Ifges(7,6) * t385 - Ifges(7,3) * t423 - t416 * t395 + t415 * t396;
t394 = Ifges(7,5) * t416 + Ifges(7,6) * t415 + Ifges(7,3) * t452;
t364 = -mrSges(7,1) * t381 + mrSges(7,3) * t376 + Ifges(7,4) * t386 + Ifges(7,2) * t385 + Ifges(7,6) * t423 - t394 * t416 + t396 * t452;
t365 = mrSges(7,2) * t381 - mrSges(7,3) * t375 + Ifges(7,1) * t386 + Ifges(7,4) * t385 + Ifges(7,5) * t423 + t394 * t415 - t395 * t452;
t407 = Ifges(6,5) * t444 + Ifges(6,6) * t443 + Ifges(6,3) * t455;
t409 = Ifges(6,1) * t444 + Ifges(6,4) * t443 + Ifges(6,5) * t455;
t349 = -mrSges(6,1) * t392 + mrSges(6,3) * t380 + Ifges(6,4) * t404 + Ifges(6,2) * t403 + Ifges(6,6) * t428 - pkin(5) * t496 + pkin(10) * t502 + t485 * t364 + t480 * t365 - t444 * t407 + t455 * t409;
t408 = Ifges(6,4) * t444 + Ifges(6,2) * t443 + Ifges(6,6) * t455;
t352 = mrSges(6,2) * t392 - mrSges(6,3) * t379 + Ifges(6,1) * t404 + Ifges(6,4) * t403 + Ifges(6,5) * t428 - pkin(10) * t363 - t364 * t480 + t365 * t485 + t407 * t443 - t408 * t455;
t435 = Ifges(5,4) * t460 + Ifges(5,2) * t459 + Ifges(5,6) * t477;
t436 = Ifges(5,1) * t460 + Ifges(5,4) * t459 + Ifges(5,5) * t477;
t493 = mrSges(5,1) * t399 - mrSges(5,2) * t400 + Ifges(5,5) * t430 + Ifges(5,6) * t429 + Ifges(5,3) * t476 + pkin(4) * t492 + pkin(9) * t503 + t486 * t349 + t481 * t352 + t460 * t435 - t459 * t436;
t491 = mrSges(6,1) * t379 - mrSges(6,2) * t380 + Ifges(6,5) * t404 + Ifges(6,6) * t403 + Ifges(6,3) * t428 + pkin(5) * t363 + t444 * t408 - t443 * t409 - t494;
t458 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t488 - Ifges(4,4) * t483) * qJD(1);
t457 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t488 - Ifges(4,2) * t483) * qJD(1);
t454 = -qJDD(1) * pkin(1) + t497;
t451 = t490 * pkin(1) - t498;
t449 = t511 * t490 + t498;
t434 = Ifges(5,5) * t460 + Ifges(5,6) * t459 + Ifges(5,3) * t477;
t345 = -mrSges(5,1) * t421 + mrSges(5,3) * t400 + Ifges(5,4) * t430 + Ifges(5,2) * t429 + Ifges(5,6) * t476 - pkin(4) * t357 - t460 * t434 + t477 * t436 - t491;
t344 = m(3) * t454 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t490) + t500;
t343 = mrSges(5,2) * t421 - mrSges(5,3) * t399 + Ifges(5,1) * t430 + Ifges(5,4) * t429 + Ifges(5,5) * t476 - pkin(9) * t357 - t349 * t481 + t352 * t486 + t434 * t459 - t435 * t477;
t1 = [mrSges(2,1) * t505 - mrSges(2,2) * t501 + mrSges(3,2) * t454 - mrSges(3,3) * t451 + t488 * (mrSges(4,2) * t449 - mrSges(4,3) * t440 + Ifges(4,1) * t465 + Ifges(4,4) * t464 + Ifges(4,5) * qJDD(3) - pkin(8) * t350 - qJD(3) * t457 + t343 * t487 - t345 * t482) - t483 * (-mrSges(4,1) * t449 + mrSges(4,3) * t441 + Ifges(4,4) * t465 + Ifges(4,2) * t464 + Ifges(4,6) * qJDD(3) - pkin(3) * t495 + pkin(8) * t504 + qJD(3) * t458 + t482 * t343 + t487 * t345) - pkin(7) * t500 - pkin(1) * t344 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t451 + m(4) * t449 - t464 * mrSges(4,1) + t490 * mrSges(3,2) + t465 * mrSges(4,2) + t495 + qJDD(1) * mrSges(3,3) + (t466 * t483 + t467 * t488) * qJD(1)) * qJ(2); t344; t493 + Ifges(4,3) * qJDD(3) + (t457 * t488 + t458 * t483) * qJD(1) + Ifges(4,6) * t464 + Ifges(4,5) * t465 + mrSges(4,1) * t440 - mrSges(4,2) * t441 + pkin(3) * t350; t493; t491; -t494;];
tauJ  = t1;
