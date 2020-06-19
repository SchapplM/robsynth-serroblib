% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% m [7x1]
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
% Datum: 2020-06-19 10:42
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR10V2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 10:17:27
% EndTime: 2020-06-19 10:17:39
% DurationCPUTime: 10.22s
% Computational Cost: add. (51811->304), mult. (102278->388), div. (0->0), fcn. (78049->12), ass. (0->119)
t408 = sin(qJ(3));
t409 = sin(qJ(2));
t414 = cos(qJ(3));
t415 = cos(qJ(2));
t393 = (t408 * t409 - t414 * t415) * qJD(1);
t417 = qJD(1) ^ 2;
t410 = sin(qJ(1));
t416 = cos(qJ(1));
t427 = -g(1) * t416 - g(2) * t410;
t397 = -pkin(1) * t417 + t427;
t387 = -t415 * g(3) - t409 * t397;
t380 = (t409 * t415 * t417 + qJDD(2)) * pkin(2) + t387;
t388 = -t409 * g(3) + t415 * t397;
t381 = (-t415 ^ 2 * t417 - qJD(2) ^ 2) * pkin(2) + t388;
t359 = t408 * t380 + t414 * t381;
t394 = (t408 * t415 + t409 * t414) * qJD(1);
t377 = pkin(3) * t393 - pkin(5) * t394;
t404 = qJD(2) + qJD(3);
t402 = t404 ^ 2;
t403 = qJDD(2) + qJDD(3);
t343 = -pkin(3) * t402 + pkin(5) * t403 - t377 * t393 + t359;
t407 = sin(qJ(4));
t413 = cos(qJ(4));
t431 = qJD(1) * qJD(2);
t398 = qJDD(1) * t409 + t415 * t431;
t430 = t409 * t431;
t399 = qJDD(1) * t415 - t430;
t368 = -qJD(3) * t394 - t398 * t408 + t399 * t414;
t369 = -qJD(3) * t393 + t398 * t414 + t399 * t408;
t429 = t410 * g(1) - t416 * g(2);
t396 = -qJDD(1) * pkin(1) - t429;
t378 = (-t399 + t430) * pkin(2) + t396;
t419 = (t393 * t404 - t369) * pkin(5) + (t394 * t404 - t368) * pkin(3) + t378;
t327 = t413 * t343 + t407 * t419;
t358 = t380 * t414 - t381 * t408;
t342 = -pkin(3) * t403 - pkin(5) * t402 + t377 * t394 - t358;
t406 = sin(qJ(5));
t412 = cos(qJ(5));
t322 = t412 * t327 + t406 * t342;
t384 = t394 * t413 + t404 * t407;
t347 = -qJD(4) * t384 - t369 * t407 + t403 * t413;
t346 = qJDD(5) - t347;
t389 = qJD(4) + t393;
t363 = -t384 * t406 + t389 * t412;
t364 = t384 * t412 + t389 * t406;
t318 = (-t363 * t364 + t346) * pkin(6) + t322;
t326 = t407 * t343 - t413 * t419;
t383 = -t394 * t407 + t404 * t413;
t348 = qJD(4) * t383 + t369 * t413 + t403 * t407;
t367 = qJDD(4) - t368;
t333 = qJD(5) * t363 + t348 * t412 + t367 * t406;
t382 = qJD(5) - t383;
t319 = (-t363 * t382 - t333) * pkin(6) + t326;
t405 = sin(qJ(6));
t411 = cos(qJ(6));
t316 = -t318 * t405 + t319 * t411;
t349 = -t364 * t405 + t382 * t411;
t324 = qJD(6) * t349 + t333 * t411 + t346 * t405;
t332 = -qJD(5) * t364 - t348 * t406 + t367 * t412;
t331 = qJDD(6) - t332;
t350 = t364 * t411 + t382 * t405;
t335 = -mrSges(7,1) * t349 + mrSges(7,2) * t350;
t362 = qJD(6) - t363;
t336 = -mrSges(7,2) * t362 + mrSges(7,3) * t349;
t314 = m(7) * t316 + mrSges(7,1) * t331 - mrSges(7,3) * t324 - t335 * t350 + t336 * t362;
t317 = t318 * t411 + t319 * t405;
t323 = -qJD(6) * t350 - t333 * t405 + t346 * t411;
t337 = mrSges(7,1) * t362 - mrSges(7,3) * t350;
t315 = m(7) * t317 - mrSges(7,2) * t331 + mrSges(7,3) * t323 + t335 * t349 - t337 * t362;
t309 = -t314 * t405 + t411 * t315;
t344 = -mrSges(6,1) * t363 + mrSges(6,2) * t364;
t352 = mrSges(6,1) * t382 - mrSges(6,3) * t364;
t308 = m(6) * t322 - mrSges(6,2) * t346 + mrSges(6,3) * t332 + t344 * t363 - t352 * t382 + t309;
t321 = -t406 * t327 + t412 * t342;
t320 = (-t364 ^ 2 - t382 ^ 2) * pkin(6) - t321;
t351 = -mrSges(6,2) * t382 + mrSges(6,3) * t363;
t312 = m(6) * t321 - m(7) * t320 + mrSges(6,1) * t346 + mrSges(7,1) * t323 - mrSges(7,2) * t324 - mrSges(6,3) * t333 + t336 * t349 - t337 * t350 - t344 * t364 + t351 * t382;
t360 = -mrSges(5,1) * t383 + mrSges(5,2) * t384;
t371 = mrSges(5,1) * t389 - mrSges(5,3) * t384;
t303 = m(5) * t327 - mrSges(5,2) * t367 + mrSges(5,3) * t347 + t308 * t412 - t312 * t406 + t360 * t383 - t371 * t389;
t370 = -mrSges(5,2) * t389 + mrSges(5,3) * t383;
t426 = t411 * t314 + t405 * t315;
t306 = t367 * mrSges(5,1) + t332 * mrSges(6,1) - t333 * mrSges(6,2) - t348 * mrSges(5,3) + t363 * t351 - t364 * t352 - t384 * t360 + t389 * t370 + (-m(5) - m(6)) * t326 - t426;
t297 = t303 * t407 + t306 * t413;
t433 = qJD(1) * t409;
t432 = qJD(1) * t415;
t428 = t303 * t413 - t306 * t407;
t328 = Ifges(7,5) * t350 + Ifges(7,6) * t349 + Ifges(7,3) * t362;
t330 = Ifges(7,1) * t350 + Ifges(7,4) * t349 + Ifges(7,5) * t362;
t310 = -mrSges(7,1) * t320 + mrSges(7,3) * t317 + Ifges(7,4) * t324 + Ifges(7,2) * t323 + Ifges(7,6) * t331 - t328 * t350 + t330 * t362;
t329 = Ifges(7,4) * t350 + Ifges(7,2) * t349 + Ifges(7,6) * t362;
t311 = mrSges(7,2) * t320 - mrSges(7,3) * t316 + Ifges(7,1) * t324 + Ifges(7,4) * t323 + Ifges(7,5) * t331 + t328 * t349 - t329 * t362;
t338 = Ifges(6,5) * t364 + Ifges(6,6) * t363 + Ifges(6,3) * t382;
t339 = Ifges(6,4) * t364 + Ifges(6,2) * t363 + Ifges(6,6) * t382;
t300 = mrSges(6,2) * t326 - mrSges(6,3) * t321 + Ifges(6,1) * t333 + Ifges(6,4) * t332 + Ifges(6,5) * t346 - pkin(6) * t426 - t405 * t310 + t411 * t311 + t363 * t338 - t382 * t339;
t340 = Ifges(6,1) * t364 + Ifges(6,4) * t363 + Ifges(6,5) * t382;
t421 = mrSges(7,1) * t316 - mrSges(7,2) * t317 + Ifges(7,5) * t324 + Ifges(7,6) * t323 + Ifges(7,3) * t331 + t329 * t350 - t330 * t349;
t307 = -mrSges(6,1) * t326 + mrSges(6,3) * t322 + Ifges(6,4) * t333 + Ifges(6,2) * t332 + Ifges(6,6) * t346 - t338 * t364 + t340 * t382 - t421;
t353 = Ifges(5,5) * t384 + Ifges(5,6) * t383 + Ifges(5,3) * t389;
t354 = Ifges(5,4) * t384 + Ifges(5,2) * t383 + Ifges(5,6) * t389;
t295 = mrSges(5,2) * t342 + mrSges(5,3) * t326 + Ifges(5,1) * t348 + Ifges(5,4) * t347 + Ifges(5,5) * t367 + t300 * t412 - t307 * t406 + t353 * t383 - t354 * t389;
t355 = Ifges(5,1) * t384 + Ifges(5,4) * t383 + Ifges(5,5) * t389;
t418 = mrSges(6,1) * t321 - mrSges(6,2) * t322 + Ifges(6,5) * t333 + Ifges(6,6) * t332 + Ifges(6,3) * t346 + pkin(6) * t309 + t310 * t411 + t311 * t405 + t339 * t364 - t340 * t363;
t299 = -mrSges(5,1) * t342 + mrSges(5,3) * t327 + Ifges(5,4) * t348 + Ifges(5,2) * t347 + Ifges(5,6) * t367 - t353 * t384 + t355 * t389 - t418;
t373 = Ifges(4,4) * t394 - Ifges(4,2) * t393 + Ifges(4,6) * t404;
t374 = Ifges(4,1) * t394 - Ifges(4,4) * t393 + Ifges(4,5) * t404;
t422 = -m(5) * t342 + t347 * mrSges(5,1) - mrSges(5,2) * t348 - t308 * t406 - t312 * t412 + t383 * t370 - t371 * t384;
t424 = mrSges(4,1) * t358 - mrSges(4,2) * t359 + Ifges(4,5) * t369 + Ifges(4,6) * t368 + Ifges(4,3) * t403 + pkin(3) * t422 + pkin(5) * t428 + t295 * t407 + t299 * t413 + t394 * t373 + t393 * t374;
t385 = -mrSges(4,2) * t404 - mrSges(4,3) * t393;
t386 = mrSges(4,1) * t404 - mrSges(4,3) * t394;
t423 = m(4) * t378 - t368 * mrSges(4,1) + t369 * mrSges(4,2) + t393 * t385 + t394 * t386 + t297;
t420 = mrSges(5,1) * t326 + mrSges(5,2) * t327 - Ifges(5,5) * t348 - Ifges(5,6) * t347 - Ifges(5,3) * t367 - t300 * t406 - t307 * t412 - t354 * t384 + t355 * t383;
t392 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t409 + Ifges(3,4) * t415) * qJD(1);
t391 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t409 + Ifges(3,2) * t415) * qJD(1);
t375 = mrSges(4,1) * t393 + mrSges(4,2) * t394;
t372 = Ifges(4,5) * t394 - Ifges(4,6) * t393 + Ifges(4,3) * t404;
t293 = -mrSges(4,1) * t378 + mrSges(4,3) * t359 + Ifges(4,4) * t369 + Ifges(4,2) * t368 + Ifges(4,6) * t403 - pkin(3) * t297 - t372 * t394 + t374 * t404 + t420;
t292 = mrSges(4,2) * t378 - mrSges(4,3) * t358 + Ifges(4,1) * t369 + Ifges(4,4) * t368 + Ifges(4,5) * t403 - pkin(5) * t297 + t295 * t413 - t299 * t407 - t372 * t393 - t373 * t404;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t429 - mrSges(2,2) * t427 + t409 * (mrSges(3,2) * t396 - mrSges(3,3) * t387 + Ifges(3,1) * t398 + Ifges(3,4) * t399 + Ifges(3,5) * qJDD(2) - qJD(2) * t391 + t292 * t414 - t293 * t408) + t415 * (-mrSges(3,1) * t396 + mrSges(3,3) * t388 + Ifges(3,4) * t398 + Ifges(3,2) * t399 + Ifges(3,6) * qJDD(2) - pkin(2) * t423 + qJD(2) * t392 + t408 * t292 + t414 * t293) + pkin(1) * (-m(3) * t396 - t398 * mrSges(3,2) + t399 * mrSges(3,1) - (qJD(2) * mrSges(3,1) - mrSges(3,3) * t433) * t433 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t432) * t432 - t423); (t409 * t391 - t415 * t392) * qJD(1) + Ifges(3,3) * qJDD(2) + t424 + pkin(2) * (t408 * (m(4) * t359 - mrSges(4,2) * t403 + mrSges(4,3) * t368 - t375 * t393 - t386 * t404 + t428) + t414 * (m(4) * t358 + mrSges(4,1) * t403 - mrSges(4,3) * t369 - t375 * t394 + t385 * t404 + t422)) + Ifges(3,5) * t398 + Ifges(3,6) * t399 + mrSges(3,1) * t387 - mrSges(3,2) * t388; t424; -t420; t418; t421;];
tauJ = t1;
