% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:17:19
% EndTime: 2019-12-31 20:17:22
% DurationCPUTime: 2.58s
% Computational Cost: add. (25225->275), mult. (58519->356), div. (0->0), fcn. (41628->10), ass. (0->109)
t419 = qJD(1) ^ 2;
t434 = pkin(2) * t419;
t414 = sin(qJ(1));
t418 = cos(qJ(1));
t425 = -t418 * g(1) - t414 * g(2);
t395 = -t419 * pkin(1) + qJDD(1) * pkin(6) + t425;
t413 = sin(qJ(2));
t433 = t413 * t395;
t417 = cos(qJ(2));
t430 = qJD(1) * qJD(2);
t398 = t413 * qJDD(1) + t417 * t430;
t363 = qJDD(2) * pkin(2) - t398 * qJ(3) - t433 + (qJ(3) * t430 + t413 * t434 - g(3)) * t417;
t381 = -t413 * g(3) + t417 * t395;
t399 = t417 * qJDD(1) - t413 * t430;
t432 = qJD(1) * t413;
t400 = qJD(2) * pkin(2) - qJ(3) * t432;
t408 = t417 ^ 2;
t364 = t399 * qJ(3) - qJD(2) * t400 - t408 * t434 + t381;
t409 = sin(pkin(9));
t410 = cos(pkin(9));
t390 = (t417 * t409 + t413 * t410) * qJD(1);
t342 = -0.2e1 * qJD(3) * t390 + t410 * t363 - t409 * t364;
t379 = t410 * t398 + t409 * t399;
t389 = (-t413 * t409 + t417 * t410) * qJD(1);
t330 = (qJD(2) * t389 - t379) * pkin(7) + (t389 * t390 + qJDD(2)) * pkin(3) + t342;
t343 = 0.2e1 * qJD(3) * t389 + t409 * t363 + t410 * t364;
t378 = -t409 * t398 + t410 * t399;
t384 = qJD(2) * pkin(3) - t390 * pkin(7);
t388 = t389 ^ 2;
t332 = -t388 * pkin(3) + t378 * pkin(7) - qJD(2) * t384 + t343;
t412 = sin(qJ(4));
t416 = cos(qJ(4));
t328 = t412 * t330 + t416 * t332;
t375 = t412 * t389 + t416 * t390;
t349 = -t375 * qJD(4) + t416 * t378 - t412 * t379;
t374 = t416 * t389 - t412 * t390;
t358 = -t374 * mrSges(5,1) + t375 * mrSges(5,2);
t407 = qJD(2) + qJD(4);
t369 = t407 * mrSges(5,1) - t375 * mrSges(5,3);
t406 = qJDD(2) + qJDD(4);
t359 = -t374 * pkin(4) - t375 * pkin(8);
t405 = t407 ^ 2;
t324 = -t405 * pkin(4) + t406 * pkin(8) + t374 * t359 + t328;
t429 = t414 * g(1) - t418 * g(2);
t424 = -qJDD(1) * pkin(1) - t429;
t365 = -t399 * pkin(2) + qJDD(3) + t400 * t432 + (-qJ(3) * t408 - pkin(6)) * t419 + t424;
t341 = -t378 * pkin(3) - t388 * pkin(7) + t390 * t384 + t365;
t350 = t374 * qJD(4) + t412 * t378 + t416 * t379;
t325 = (-t374 * t407 - t350) * pkin(8) + (t375 * t407 - t349) * pkin(4) + t341;
t411 = sin(qJ(5));
t415 = cos(qJ(5));
t321 = -t411 * t324 + t415 * t325;
t366 = -t411 * t375 + t415 * t407;
t335 = t366 * qJD(5) + t415 * t350 + t411 * t406;
t348 = qJDD(5) - t349;
t367 = t415 * t375 + t411 * t407;
t351 = -t366 * mrSges(6,1) + t367 * mrSges(6,2);
t370 = qJD(5) - t374;
t352 = -t370 * mrSges(6,2) + t366 * mrSges(6,3);
t318 = m(6) * t321 + t348 * mrSges(6,1) - t335 * mrSges(6,3) - t367 * t351 + t370 * t352;
t322 = t415 * t324 + t411 * t325;
t334 = -t367 * qJD(5) - t411 * t350 + t415 * t406;
t353 = t370 * mrSges(6,1) - t367 * mrSges(6,3);
t319 = m(6) * t322 - t348 * mrSges(6,2) + t334 * mrSges(6,3) + t366 * t351 - t370 * t353;
t426 = -t411 * t318 + t415 * t319;
t305 = m(5) * t328 - t406 * mrSges(5,2) + t349 * mrSges(5,3) + t374 * t358 - t407 * t369 + t426;
t327 = t416 * t330 - t412 * t332;
t368 = -t407 * mrSges(5,2) + t374 * mrSges(5,3);
t323 = -t406 * pkin(4) - t405 * pkin(8) + t375 * t359 - t327;
t423 = -m(6) * t323 + t334 * mrSges(6,1) - t335 * mrSges(6,2) + t366 * t352 - t367 * t353;
t314 = m(5) * t327 + t406 * mrSges(5,1) - t350 * mrSges(5,3) - t375 * t358 + t407 * t368 + t423;
t302 = t412 * t305 + t416 * t314;
t377 = -t389 * mrSges(4,1) + t390 * mrSges(4,2);
t382 = -qJD(2) * mrSges(4,2) + t389 * mrSges(4,3);
t300 = m(4) * t342 + qJDD(2) * mrSges(4,1) - t379 * mrSges(4,3) + qJD(2) * t382 - t390 * t377 + t302;
t383 = qJD(2) * mrSges(4,1) - t390 * mrSges(4,3);
t427 = t416 * t305 - t412 * t314;
t301 = m(4) * t343 - qJDD(2) * mrSges(4,2) + t378 * mrSges(4,3) - qJD(2) * t383 + t389 * t377 + t427;
t294 = t410 * t300 + t409 * t301;
t308 = t415 * t318 + t411 * t319;
t431 = qJD(1) * t417;
t428 = -t409 * t300 + t410 * t301;
t422 = -m(5) * t341 + t349 * mrSges(5,1) - t350 * mrSges(5,2) + t374 * t368 - t375 * t369 - t308;
t336 = Ifges(6,5) * t367 + Ifges(6,6) * t366 + Ifges(6,3) * t370;
t338 = Ifges(6,1) * t367 + Ifges(6,4) * t366 + Ifges(6,5) * t370;
t311 = -mrSges(6,1) * t323 + mrSges(6,3) * t322 + Ifges(6,4) * t335 + Ifges(6,2) * t334 + Ifges(6,6) * t348 - t367 * t336 + t370 * t338;
t337 = Ifges(6,4) * t367 + Ifges(6,2) * t366 + Ifges(6,6) * t370;
t312 = mrSges(6,2) * t323 - mrSges(6,3) * t321 + Ifges(6,1) * t335 + Ifges(6,4) * t334 + Ifges(6,5) * t348 + t366 * t336 - t370 * t337;
t355 = Ifges(5,4) * t375 + Ifges(5,2) * t374 + Ifges(5,6) * t407;
t356 = Ifges(5,1) * t375 + Ifges(5,4) * t374 + Ifges(5,5) * t407;
t421 = mrSges(5,1) * t327 - mrSges(5,2) * t328 + Ifges(5,5) * t350 + Ifges(5,6) * t349 + Ifges(5,3) * t406 + pkin(4) * t423 + pkin(8) * t426 + t415 * t311 + t411 * t312 + t375 * t355 - t374 * t356;
t420 = mrSges(6,1) * t321 - mrSges(6,2) * t322 + Ifges(6,5) * t335 + Ifges(6,6) * t334 + Ifges(6,3) * t348 + t367 * t337 - t366 * t338;
t306 = m(4) * t365 - t378 * mrSges(4,1) + t379 * mrSges(4,2) - t389 * t382 + t390 * t383 - t422;
t402 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t431;
t401 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t432;
t397 = (-t417 * mrSges(3,1) + t413 * mrSges(3,2)) * qJD(1);
t394 = -t419 * pkin(6) + t424;
t393 = Ifges(3,5) * qJD(2) + (t413 * Ifges(3,1) + t417 * Ifges(3,4)) * qJD(1);
t392 = Ifges(3,6) * qJD(2) + (t413 * Ifges(3,4) + t417 * Ifges(3,2)) * qJD(1);
t380 = -t417 * g(3) - t433;
t373 = Ifges(4,1) * t390 + Ifges(4,4) * t389 + Ifges(4,5) * qJD(2);
t372 = Ifges(4,4) * t390 + Ifges(4,2) * t389 + Ifges(4,6) * qJD(2);
t371 = Ifges(4,5) * t390 + Ifges(4,6) * t389 + Ifges(4,3) * qJD(2);
t354 = Ifges(5,5) * t375 + Ifges(5,6) * t374 + Ifges(5,3) * t407;
t296 = -mrSges(5,1) * t341 + mrSges(5,3) * t328 + Ifges(5,4) * t350 + Ifges(5,2) * t349 + Ifges(5,6) * t406 - pkin(4) * t308 - t375 * t354 + t407 * t356 - t420;
t295 = mrSges(5,2) * t341 - mrSges(5,3) * t327 + Ifges(5,1) * t350 + Ifges(5,4) * t349 + Ifges(5,5) * t406 - pkin(8) * t308 - t411 * t311 + t415 * t312 + t374 * t354 - t407 * t355;
t293 = mrSges(4,2) * t365 - mrSges(4,3) * t342 + Ifges(4,1) * t379 + Ifges(4,4) * t378 + Ifges(4,5) * qJDD(2) - pkin(7) * t302 - qJD(2) * t372 + t416 * t295 - t412 * t296 + t389 * t371;
t292 = -mrSges(4,1) * t365 + mrSges(4,3) * t343 + Ifges(4,4) * t379 + Ifges(4,2) * t378 + Ifges(4,6) * qJDD(2) + pkin(3) * t422 + pkin(7) * t427 + qJD(2) * t373 + t412 * t295 + t416 * t296 - t390 * t371;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t429 - mrSges(2,2) * t425 + t413 * (mrSges(3,2) * t394 - mrSges(3,3) * t380 + Ifges(3,1) * t398 + Ifges(3,4) * t399 + Ifges(3,5) * qJDD(2) - qJ(3) * t294 - qJD(2) * t392 - t409 * t292 + t410 * t293) + t417 * (-mrSges(3,1) * t394 + mrSges(3,3) * t381 + Ifges(3,4) * t398 + Ifges(3,2) * t399 + Ifges(3,6) * qJDD(2) - pkin(2) * t306 + qJ(3) * t428 + qJD(2) * t393 + t410 * t292 + t409 * t293) + pkin(1) * (-t398 * mrSges(3,2) + t399 * mrSges(3,1) - m(3) * t394 + (-t401 * t413 + t402 * t417) * qJD(1) - t306) + pkin(6) * (t417 * (m(3) * t381 - qJDD(2) * mrSges(3,2) + t399 * mrSges(3,3) - qJD(2) * t401 + t397 * t431 + t428) - t413 * (m(3) * t380 + qJDD(2) * mrSges(3,1) - t398 * mrSges(3,3) + qJD(2) * t402 - t397 * t432 + t294)); t421 + Ifges(3,5) * t398 + Ifges(3,6) * t399 - t389 * t373 + t390 * t372 + Ifges(4,6) * t378 + Ifges(4,5) * t379 + mrSges(3,1) * t380 - mrSges(3,2) * t381 + mrSges(4,1) * t342 - mrSges(4,2) * t343 + pkin(3) * t302 + pkin(2) * t294 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t413 * t392 - t417 * t393) * qJD(1); t306; t421; t420;];
tauJ = t1;
