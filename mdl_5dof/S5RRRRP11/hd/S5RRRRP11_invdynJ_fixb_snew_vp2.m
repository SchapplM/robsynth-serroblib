% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:56
% EndTime: 2019-12-31 22:15:01
% DurationCPUTime: 2.75s
% Computational Cost: add. (24294->264), mult. (51880->339), div. (0->0), fcn. (39450->10), ass. (0->116)
t454 = Ifges(5,1) + Ifges(6,1);
t446 = Ifges(5,4) - Ifges(6,5);
t445 = -Ifges(5,5) - Ifges(6,4);
t453 = Ifges(5,2) + Ifges(6,3);
t444 = Ifges(5,6) - Ifges(6,6);
t452 = -Ifges(5,3) - Ifges(6,2);
t411 = sin(pkin(5));
t415 = sin(qJ(2));
t418 = cos(qJ(2));
t432 = qJD(1) * qJD(2);
t402 = (-qJDD(1) * t418 + t415 * t432) * t411;
t412 = cos(pkin(5));
t408 = t412 * qJD(1) + qJD(2);
t414 = sin(qJ(3));
t417 = cos(qJ(3));
t434 = qJD(1) * t411;
t430 = t415 * t434;
t390 = t417 * t408 - t414 * t430;
t401 = (qJDD(1) * t415 + t418 * t432) * t411;
t407 = t412 * qJDD(1) + qJDD(2);
t372 = t390 * qJD(3) + t417 * t401 + t414 * t407;
t391 = t414 * t408 + t417 * t430;
t433 = qJD(1) * t418;
t429 = t411 * t433;
t405 = qJD(3) - t429;
t413 = sin(qJ(4));
t450 = cos(qJ(4));
t377 = t413 * t391 - t450 * t405;
t394 = qJDD(3) + t402;
t338 = -t377 * qJD(4) + t450 * t372 + t413 * t394;
t378 = t450 * t391 + t413 * t405;
t356 = t377 * mrSges(6,1) - t378 * mrSges(6,3);
t400 = (-t418 * pkin(2) - t415 * pkin(8)) * t434;
t406 = t408 ^ 2;
t420 = qJD(1) ^ 2;
t416 = sin(qJ(1));
t419 = cos(qJ(1));
t425 = -t419 * g(1) - t416 * g(2);
t449 = t411 * pkin(7);
t398 = -t420 * pkin(1) + qJDD(1) * t449 + t425;
t428 = t416 * g(1) - t419 * g(2);
t397 = qJDD(1) * pkin(1) + t420 * t449 + t428;
t442 = t397 * t412;
t435 = t418 * t398 + t415 * t442;
t351 = -t406 * pkin(2) + t407 * pkin(8) + (-g(3) * t415 + t400 * t433) * t411 + t435;
t448 = t412 * g(3);
t352 = t402 * pkin(2) - t401 * pkin(8) - t448 + (-t397 + (pkin(2) * t415 - pkin(8) * t418) * t408 * qJD(1)) * t411;
t334 = t417 * t351 + t414 * t352;
t376 = -t390 * pkin(3) - t391 * pkin(9);
t404 = t405 ^ 2;
t330 = -t404 * pkin(3) + t394 * pkin(9) + t390 * t376 + t334;
t440 = t411 * t418;
t373 = -g(3) * t440 - t415 * t398 + t418 * t442;
t350 = -t407 * pkin(2) - t406 * pkin(8) + t400 * t430 - t373;
t371 = -t391 * qJD(3) - t414 * t401 + t417 * t407;
t332 = (-t390 * t405 - t372) * pkin(9) + (t391 * t405 - t371) * pkin(3) + t350;
t326 = -t413 * t330 + t450 * t332;
t355 = t377 * pkin(4) - t378 * qJ(5);
t369 = qJDD(4) - t371;
t388 = qJD(4) - t390;
t387 = t388 ^ 2;
t324 = -t369 * pkin(4) - t387 * qJ(5) + t378 * t355 + qJDD(5) - t326;
t359 = -t377 * mrSges(6,2) + t388 * mrSges(6,3);
t426 = -m(6) * t324 + t369 * mrSges(6,1) + t388 * t359;
t320 = t338 * mrSges(6,2) + t378 * t356 - t426;
t327 = t450 * t330 + t413 * t332;
t323 = -t387 * pkin(4) + t369 * qJ(5) + 0.2e1 * qJD(5) * t388 - t377 * t355 + t327;
t337 = t378 * qJD(4) + t413 * t372 - t450 * t394;
t362 = -t388 * mrSges(6,1) + t378 * mrSges(6,2);
t431 = m(6) * t323 + t369 * mrSges(6,3) + t388 * t362;
t437 = t446 * t377 - t454 * t378 + t445 * t388;
t438 = t453 * t377 - t446 * t378 - t444 * t388;
t451 = -t444 * t337 - t445 * t338 - t452 * t369 - t437 * t377 - t438 * t378 + mrSges(5,1) * t326 - mrSges(6,1) * t324 - mrSges(5,2) * t327 + mrSges(6,3) * t323 - pkin(4) * t320 + qJ(5) * (-t337 * mrSges(6,2) - t377 * t356 + t431);
t447 = -mrSges(5,3) - mrSges(6,2);
t441 = t411 * t415;
t361 = t388 * mrSges(5,1) - t378 * mrSges(5,3);
t436 = -t377 * mrSges(5,1) - t378 * mrSges(5,2) - t356;
t316 = m(5) * t327 - t369 * mrSges(5,2) + t447 * t337 - t388 * t361 + t436 * t377 + t431;
t360 = -t388 * mrSges(5,2) - t377 * mrSges(5,3);
t317 = m(5) * t326 + t369 * mrSges(5,1) + t447 * t338 + t388 * t360 + t436 * t378 + t426;
t312 = t450 * t316 - t413 * t317;
t375 = -t390 * mrSges(4,1) + t391 * mrSges(4,2);
t380 = t405 * mrSges(4,1) - t391 * mrSges(4,3);
t309 = m(4) * t334 - t394 * mrSges(4,2) + t371 * mrSges(4,3) + t390 * t375 - t405 * t380 + t312;
t333 = -t414 * t351 + t417 * t352;
t329 = -t394 * pkin(3) - t404 * pkin(9) + t391 * t376 - t333;
t325 = -0.2e1 * qJD(5) * t378 + (t377 * t388 - t338) * qJ(5) + (t378 * t388 + t337) * pkin(4) + t329;
t321 = m(6) * t325 + t337 * mrSges(6,1) - t338 * mrSges(6,3) + t377 * t359 - t378 * t362;
t318 = -m(5) * t329 - t337 * mrSges(5,1) - t338 * mrSges(5,2) - t377 * t360 - t378 * t361 - t321;
t379 = -t405 * mrSges(4,2) + t390 * mrSges(4,3);
t314 = m(4) * t333 + t394 * mrSges(4,1) - t372 * mrSges(4,3) - t391 * t375 + t405 * t379 + t318;
t304 = t414 * t309 + t417 * t314;
t439 = t444 * t377 + t445 * t378 + t452 * t388;
t427 = t417 * t309 - t414 * t314;
t311 = t413 * t316 + t450 * t317;
t422 = -m(4) * t350 + t371 * mrSges(4,1) - t372 * mrSges(4,2) + t390 * t379 - t391 * t380 - t311;
t306 = -mrSges(5,1) * t329 - mrSges(6,1) * t325 + mrSges(6,2) * t323 + mrSges(5,3) * t327 - pkin(4) * t321 - t453 * t337 + t446 * t338 + t444 * t369 + t439 * t378 - t437 * t388;
t310 = mrSges(5,2) * t329 + mrSges(6,2) * t324 - mrSges(5,3) * t326 - mrSges(6,3) * t325 - qJ(5) * t321 - t446 * t337 + t454 * t338 - t445 * t369 + t439 * t377 + t438 * t388;
t366 = Ifges(4,4) * t391 + Ifges(4,2) * t390 + Ifges(4,6) * t405;
t367 = Ifges(4,1) * t391 + Ifges(4,4) * t390 + Ifges(4,5) * t405;
t421 = mrSges(4,1) * t333 - mrSges(4,2) * t334 + Ifges(4,5) * t372 + Ifges(4,6) * t371 + Ifges(4,3) * t394 + pkin(3) * t318 + pkin(9) * t312 + t450 * t306 + t413 * t310 + t391 * t366 - t390 * t367;
t399 = (-t418 * mrSges(3,1) + t415 * mrSges(3,2)) * t434;
t396 = -t408 * mrSges(3,2) + mrSges(3,3) * t429;
t395 = t408 * mrSges(3,1) - mrSges(3,3) * t430;
t384 = -t411 * t397 - t448;
t383 = Ifges(3,5) * t408 + (t415 * Ifges(3,1) + t418 * Ifges(3,4)) * t434;
t382 = Ifges(3,6) * t408 + (t415 * Ifges(3,4) + t418 * Ifges(3,2)) * t434;
t381 = Ifges(3,3) * t408 + (t415 * Ifges(3,5) + t418 * Ifges(3,6)) * t434;
t374 = -g(3) * t441 + t435;
t365 = Ifges(4,5) * t391 + Ifges(4,6) * t390 + Ifges(4,3) * t405;
t305 = m(3) * t373 + t407 * mrSges(3,1) - t401 * mrSges(3,3) + t408 * t396 - t399 * t430 + t422;
t303 = m(3) * t374 - t407 * mrSges(3,2) - t402 * mrSges(3,3) - t408 * t395 + t399 * t429 + t427;
t302 = -mrSges(4,1) * t350 + mrSges(4,3) * t334 + Ifges(4,4) * t372 + Ifges(4,2) * t371 + Ifges(4,6) * t394 - pkin(3) * t311 - t391 * t365 + t405 * t367 - t451;
t301 = mrSges(4,2) * t350 - mrSges(4,3) * t333 + Ifges(4,1) * t372 + Ifges(4,4) * t371 + Ifges(4,5) * t394 - pkin(9) * t311 - t413 * t306 + t450 * t310 + t390 * t365 - t405 * t366;
t300 = Ifges(3,5) * t401 - Ifges(3,6) * t402 + Ifges(3,3) * t407 + mrSges(3,1) * t373 - mrSges(3,2) * t374 + t414 * t301 + t417 * t302 + pkin(2) * t422 + pkin(8) * t427 + (t415 * t382 - t418 * t383) * t434;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t428 - mrSges(2,2) * t425 + (mrSges(3,2) * t384 - mrSges(3,3) * t373 + Ifges(3,1) * t401 - Ifges(3,4) * t402 + Ifges(3,5) * t407 - pkin(8) * t304 + t417 * t301 - t414 * t302 + t381 * t429 - t408 * t382) * t441 + (-mrSges(3,1) * t384 + mrSges(3,3) * t374 + Ifges(3,4) * t401 - Ifges(3,2) * t402 + Ifges(3,6) * t407 - pkin(2) * t304 - t381 * t430 + t408 * t383 - t421) * t440 + t412 * t300 + pkin(1) * ((t415 * t303 + t418 * t305) * t412 + (-m(3) * t384 - t402 * mrSges(3,1) - t401 * mrSges(3,2) + (-t395 * t415 + t396 * t418) * t434 - t304) * t411) + (t418 * t303 - t415 * t305) * t449; t300; t421; t451; t320;];
tauJ = t1;
