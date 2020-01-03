% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR12_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:48:27
% EndTime: 2019-12-31 22:48:36
% DurationCPUTime: 8.88s
% Computational Cost: add. (133801->306), mult. (332454->420), div. (0->0), fcn. (274768->14), ass. (0->141)
t450 = cos(pkin(5));
t445 = qJD(1) * t450 + qJD(2);
t447 = sin(pkin(6));
t449 = cos(pkin(6));
t448 = sin(pkin(5));
t459 = cos(qJ(2));
t478 = qJD(1) * t459;
t474 = t448 * t478;
t432 = (t445 * t447 + t449 * t474) * pkin(9);
t454 = sin(qJ(2));
t480 = qJD(1) * t448;
t492 = pkin(9) * t447;
t436 = (-pkin(2) * t459 - t454 * t492) * t480;
t477 = qJD(1) * qJD(2);
t442 = (qJDD(1) * t454 + t459 * t477) * t448;
t444 = qJDD(1) * t450 + qJDD(2);
t461 = qJD(1) ^ 2;
t455 = sin(qJ(1));
t460 = cos(qJ(1));
t470 = -g(1) * t460 - g(2) * t455;
t490 = t448 * pkin(8);
t440 = -pkin(1) * t461 + qJDD(1) * t490 + t470;
t473 = t455 * g(1) - g(2) * t460;
t439 = qJDD(1) * pkin(1) + t461 * t490 + t473;
t488 = t439 * t450;
t471 = -t454 * t440 + t459 * t488;
t479 = qJD(1) * t454;
t491 = pkin(9) * t449;
t393 = -t442 * t491 + t444 * pkin(2) + t445 * t432 + (-g(3) * t459 - t436 * t479) * t448 + t471;
t475 = t448 * t479;
t435 = pkin(2) * t445 - t475 * t491;
t443 = (qJDD(1) * t459 - t454 * t477) * t448;
t468 = t443 * t449 + t444 * t447;
t481 = t459 * t440 + t454 * t488;
t394 = -t445 * t435 + (-g(3) * t454 + t436 * t478) * t448 + t468 * pkin(9) + t481;
t489 = t450 * g(3);
t399 = -t442 * t492 - t443 * pkin(2) - t489 + (-t439 + (-t432 * t459 + t435 * t454) * qJD(1)) * t448;
t453 = sin(qJ(3));
t458 = cos(qJ(3));
t372 = -t453 * t394 + (t393 * t449 + t399 * t447) * t458;
t482 = t449 * t459;
t487 = t447 * t453;
t423 = t445 * t487 + (t453 * t482 + t454 * t458) * t480;
t409 = -t423 * qJD(3) - t453 * t442 + t458 * t468;
t486 = t447 * t458;
t422 = (-t453 * t454 + t458 * t482) * t480 + t445 * t486;
t485 = t448 * t454;
t484 = t448 * t459;
t483 = t449 * t453;
t373 = t393 * t483 + t458 * t394 + t399 * t487;
t412 = -pkin(3) * t422 - pkin(10) * t423;
t424 = -t443 * t447 + t444 * t449 + qJDD(3);
t433 = t445 * t449 - t447 * t474 + qJD(3);
t431 = t433 ^ 2;
t366 = -pkin(3) * t431 + pkin(10) * t424 + t412 * t422 + t373;
t378 = -t447 * t393 + t449 * t399;
t410 = t422 * qJD(3) + t458 * t442 + t453 * t468;
t368 = (-t422 * t433 - t410) * pkin(10) + (t423 * t433 - t409) * pkin(3) + t378;
t452 = sin(qJ(4));
t457 = cos(qJ(4));
t363 = t366 * t457 + t368 * t452;
t414 = -t423 * t452 + t433 * t457;
t415 = t423 * t457 + t433 * t452;
t396 = -pkin(4) * t414 - pkin(11) * t415;
t408 = qJDD(4) - t409;
t421 = qJD(4) - t422;
t420 = t421 ^ 2;
t360 = -pkin(4) * t420 + pkin(11) * t408 + t396 * t414 + t363;
t365 = -t424 * pkin(3) - t431 * pkin(10) + t423 * t412 - t372;
t382 = -qJD(4) * t415 - t410 * t452 + t424 * t457;
t383 = qJD(4) * t414 + t410 * t457 + t424 * t452;
t361 = (-t414 * t421 - t383) * pkin(11) + (t415 * t421 - t382) * pkin(4) + t365;
t451 = sin(qJ(5));
t456 = cos(qJ(5));
t357 = -t360 * t451 + t361 * t456;
t401 = -t415 * t451 + t421 * t456;
t371 = qJD(5) * t401 + t383 * t456 + t408 * t451;
t402 = t415 * t456 + t421 * t451;
t379 = -mrSges(6,1) * t401 + mrSges(6,2) * t402;
t381 = qJDD(5) - t382;
t413 = qJD(5) - t414;
t384 = -mrSges(6,2) * t413 + mrSges(6,3) * t401;
t354 = m(6) * t357 + mrSges(6,1) * t381 - mrSges(6,3) * t371 - t379 * t402 + t384 * t413;
t358 = t360 * t456 + t361 * t451;
t370 = -qJD(5) * t402 - t383 * t451 + t408 * t456;
t385 = mrSges(6,1) * t413 - mrSges(6,3) * t402;
t355 = m(6) * t358 - mrSges(6,2) * t381 + mrSges(6,3) * t370 + t379 * t401 - t385 * t413;
t348 = -t354 * t451 + t355 * t456;
t395 = -mrSges(5,1) * t414 + mrSges(5,2) * t415;
t404 = mrSges(5,1) * t421 - mrSges(5,3) * t415;
t346 = m(5) * t363 - mrSges(5,2) * t408 + mrSges(5,3) * t382 + t395 * t414 - t404 * t421 + t348;
t362 = -t366 * t452 + t368 * t457;
t359 = -pkin(4) * t408 - pkin(11) * t420 + t396 * t415 - t362;
t356 = -m(6) * t359 + t370 * mrSges(6,1) - mrSges(6,2) * t371 + t401 * t384 - t385 * t402;
t403 = -mrSges(5,2) * t421 + mrSges(5,3) * t414;
t352 = m(5) * t362 + mrSges(5,1) * t408 - mrSges(5,3) * t383 - t395 * t415 + t403 * t421 + t356;
t340 = t346 * t452 + t352 * t457;
t411 = -mrSges(4,1) * t422 + mrSges(4,2) * t423;
t417 = mrSges(4,1) * t433 - mrSges(4,3) * t423;
t472 = t346 * t457 - t352 * t452;
t337 = m(4) * t373 - mrSges(4,2) * t424 + mrSges(4,3) * t409 + t411 * t422 - t417 * t433 + t472;
t416 = -mrSges(4,2) * t433 + mrSges(4,3) * t422;
t339 = m(4) * t378 - mrSges(4,1) * t409 + mrSges(4,2) * t410 - t416 * t422 + t417 * t423 + t340;
t347 = t354 * t456 + t355 * t451;
t464 = -m(5) * t365 + t382 * mrSges(5,1) - mrSges(5,2) * t383 + t414 * t403 - t404 * t415 - t347;
t343 = m(4) * t372 + mrSges(4,1) * t424 - mrSges(4,3) * t410 - t411 * t423 + t416 * t433 + t464;
t328 = t337 * t487 + t339 * t449 + t343 * t486;
t331 = t337 * t458 - t343 * t453;
t329 = t343 * t449 * t458 + t337 * t483 - t339 * t447;
t374 = Ifges(6,5) * t402 + Ifges(6,6) * t401 + Ifges(6,3) * t413;
t376 = Ifges(6,1) * t402 + Ifges(6,4) * t401 + Ifges(6,5) * t413;
t349 = -mrSges(6,1) * t359 + mrSges(6,3) * t358 + Ifges(6,4) * t371 + Ifges(6,2) * t370 + Ifges(6,6) * t381 - t374 * t402 + t376 * t413;
t375 = Ifges(6,4) * t402 + Ifges(6,2) * t401 + Ifges(6,6) * t413;
t350 = mrSges(6,2) * t359 - mrSges(6,3) * t357 + Ifges(6,1) * t371 + Ifges(6,4) * t370 + Ifges(6,5) * t381 + t374 * t401 - t375 * t413;
t386 = Ifges(5,5) * t415 + Ifges(5,6) * t414 + Ifges(5,3) * t421;
t387 = Ifges(5,4) * t415 + Ifges(5,2) * t414 + Ifges(5,6) * t421;
t332 = mrSges(5,2) * t365 - mrSges(5,3) * t362 + Ifges(5,1) * t383 + Ifges(5,4) * t382 + Ifges(5,5) * t408 - pkin(11) * t347 - t349 * t451 + t350 * t456 + t386 * t414 - t387 * t421;
t388 = Ifges(5,1) * t415 + Ifges(5,4) * t414 + Ifges(5,5) * t421;
t463 = mrSges(6,1) * t357 - mrSges(6,2) * t358 + Ifges(6,5) * t371 + Ifges(6,6) * t370 + Ifges(6,3) * t381 + t375 * t402 - t376 * t401;
t333 = -mrSges(5,1) * t365 + mrSges(5,3) * t363 + Ifges(5,4) * t383 + Ifges(5,2) * t382 + Ifges(5,6) * t408 - pkin(4) * t347 - t386 * t415 + t388 * t421 - t463;
t405 = Ifges(4,5) * t423 + Ifges(4,6) * t422 + Ifges(4,3) * t433;
t406 = Ifges(4,4) * t423 + Ifges(4,2) * t422 + Ifges(4,6) * t433;
t325 = mrSges(4,2) * t378 - mrSges(4,3) * t372 + Ifges(4,1) * t410 + Ifges(4,4) * t409 + Ifges(4,5) * t424 - pkin(10) * t340 + t332 * t457 - t333 * t452 + t405 * t422 - t406 * t433;
t407 = Ifges(4,1) * t423 + Ifges(4,4) * t422 + Ifges(4,5) * t433;
t462 = mrSges(5,1) * t362 - mrSges(5,2) * t363 + Ifges(5,5) * t383 + Ifges(5,6) * t382 + Ifges(5,3) * t408 + pkin(4) * t356 + pkin(11) * t348 + t349 * t456 + t350 * t451 + t387 * t415 - t388 * t414;
t326 = -mrSges(4,1) * t378 + mrSges(4,3) * t373 + Ifges(4,4) * t410 + Ifges(4,2) * t409 + Ifges(4,6) * t424 - pkin(3) * t340 - t405 * t423 + t407 * t433 - t462;
t465 = pkin(9) * t331 + t325 * t453 + t326 * t458;
t441 = (-mrSges(3,1) * t459 + mrSges(3,2) * t454) * t480;
t438 = -mrSges(3,2) * t445 + mrSges(3,3) * t474;
t437 = mrSges(3,1) * t445 - mrSges(3,3) * t475;
t428 = -t448 * t439 - t489;
t427 = Ifges(3,5) * t445 + (Ifges(3,1) * t454 + Ifges(3,4) * t459) * t480;
t426 = Ifges(3,6) * t445 + (Ifges(3,4) * t454 + Ifges(3,2) * t459) * t480;
t425 = Ifges(3,3) * t445 + (Ifges(3,5) * t454 + Ifges(3,6) * t459) * t480;
t419 = -g(3) * t485 + t481;
t418 = -g(3) * t484 + t471;
t330 = m(3) * t419 - mrSges(3,2) * t444 + mrSges(3,3) * t443 - t437 * t445 + t441 * t474 + t331;
t327 = m(3) * t418 + mrSges(3,1) * t444 - mrSges(3,3) * t442 + t438 * t445 - t441 * t475 + t329;
t324 = mrSges(4,1) * t372 - mrSges(4,2) * t373 + Ifges(4,5) * t410 + Ifges(4,6) * t409 + Ifges(4,3) * t424 + pkin(3) * t464 + pkin(10) * t472 + t452 * t332 + t457 * t333 + t423 * t406 - t422 * t407;
t323 = mrSges(3,1) * t418 - mrSges(3,2) * t419 + Ifges(3,5) * t442 + Ifges(3,6) * t443 + Ifges(3,3) * t444 + pkin(2) * t329 + t449 * t324 + (t426 * t454 - t427 * t459) * t480 + t465 * t447;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t473 - mrSges(2,2) * t470 + (t425 * t474 + mrSges(3,2) * t428 - mrSges(3,3) * t418 + Ifges(3,1) * t442 + Ifges(3,4) * t443 + Ifges(3,5) * t444 + t458 * t325 - t453 * t326 - t445 * t426 + (-t328 * t447 - t329 * t449) * pkin(9)) * t485 + (-mrSges(3,1) * t428 + mrSges(3,3) * t419 + Ifges(3,4) * t442 + Ifges(3,2) * t443 + Ifges(3,6) * t444 - pkin(2) * t328 - t447 * t324 - t425 * t475 + t445 * t427 + t449 * t465) * t484 + t450 * t323 + pkin(1) * ((t327 * t459 + t330 * t454) * t450 + (-m(3) * t428 + t443 * mrSges(3,1) - t442 * mrSges(3,2) + (-t437 * t454 + t438 * t459) * t480 - t328) * t448) + (-t327 * t454 + t330 * t459) * t490; t323; t324; t462; t463;];
tauJ = t1;
