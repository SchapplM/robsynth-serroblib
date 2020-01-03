% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:52
% EndTime: 2019-12-31 22:24:56
% DurationCPUTime: 2.99s
% Computational Cost: add. (32938->277), mult. (66211->354), div. (0->0), fcn. (46580->10), ass. (0->112)
t429 = qJD(1) ^ 2;
t446 = pkin(2) * t429;
t423 = sin(qJ(1));
t428 = cos(qJ(1));
t437 = -t428 * g(1) - t423 * g(2);
t403 = -t429 * pkin(1) + qJDD(1) * pkin(6) + t437;
t422 = sin(qJ(2));
t445 = t422 * t403;
t427 = cos(qJ(2));
t442 = qJD(1) * qJD(2);
t407 = t422 * qJDD(1) + t427 * t442;
t367 = qJDD(2) * pkin(2) - t407 * pkin(7) - t445 + (pkin(7) * t442 + t422 * t446 - g(3)) * t427;
t390 = -t422 * g(3) + t427 * t403;
t408 = t427 * qJDD(1) - t422 * t442;
t444 = qJD(1) * t422;
t411 = qJD(2) * pkin(2) - pkin(7) * t444;
t418 = t427 ^ 2;
t368 = t408 * pkin(7) - qJD(2) * t411 - t418 * t446 + t390;
t421 = sin(qJ(3));
t426 = cos(qJ(3));
t350 = t421 * t367 + t426 * t368;
t401 = (t427 * t421 + t422 * t426) * qJD(1);
t375 = -t401 * qJD(3) - t421 * t407 + t426 * t408;
t443 = qJD(1) * t427;
t400 = -t421 * t444 + t426 * t443;
t384 = -t400 * mrSges(4,1) + t401 * mrSges(4,2);
t417 = qJD(2) + qJD(3);
t392 = t417 * mrSges(4,1) - t401 * mrSges(4,3);
t416 = qJDD(2) + qJDD(3);
t376 = t400 * qJD(3) + t426 * t407 + t421 * t408;
t441 = t423 * g(1) - t428 * g(2);
t436 = -qJDD(1) * pkin(1) - t441;
t377 = -t408 * pkin(2) + t411 * t444 + (-pkin(7) * t418 - pkin(6)) * t429 + t436;
t339 = (-t400 * t417 - t376) * pkin(8) + (t401 * t417 - t375) * pkin(3) + t377;
t385 = -t400 * pkin(3) - t401 * pkin(8);
t415 = t417 ^ 2;
t342 = -t415 * pkin(3) + t416 * pkin(8) + t400 * t385 + t350;
t420 = sin(qJ(4));
t425 = cos(qJ(4));
t328 = t425 * t339 - t420 * t342;
t387 = -t420 * t401 + t425 * t417;
t353 = t387 * qJD(4) + t425 * t376 + t420 * t416;
t374 = qJDD(4) - t375;
t388 = t425 * t401 + t420 * t417;
t396 = qJD(4) - t400;
t326 = (t387 * t396 - t353) * pkin(9) + (t387 * t388 + t374) * pkin(4) + t328;
t329 = t420 * t339 + t425 * t342;
t352 = -t388 * qJD(4) - t420 * t376 + t425 * t416;
t380 = t396 * pkin(4) - t388 * pkin(9);
t386 = t387 ^ 2;
t327 = -t386 * pkin(4) + t352 * pkin(9) - t396 * t380 + t329;
t419 = sin(qJ(5));
t424 = cos(qJ(5));
t324 = t424 * t326 - t419 * t327;
t360 = t424 * t387 - t419 * t388;
t335 = t360 * qJD(5) + t419 * t352 + t424 * t353;
t361 = t419 * t387 + t424 * t388;
t347 = -t360 * mrSges(6,1) + t361 * mrSges(6,2);
t394 = qJD(5) + t396;
t354 = -t394 * mrSges(6,2) + t360 * mrSges(6,3);
t370 = qJDD(5) + t374;
t320 = m(6) * t324 + t370 * mrSges(6,1) - t335 * mrSges(6,3) - t361 * t347 + t394 * t354;
t325 = t419 * t326 + t424 * t327;
t334 = -t361 * qJD(5) + t424 * t352 - t419 * t353;
t355 = t394 * mrSges(6,1) - t361 * mrSges(6,3);
t321 = m(6) * t325 - t370 * mrSges(6,2) + t334 * mrSges(6,3) + t360 * t347 - t394 * t355;
t312 = t424 * t320 + t419 * t321;
t365 = -t387 * mrSges(5,1) + t388 * mrSges(5,2);
t378 = -t396 * mrSges(5,2) + t387 * mrSges(5,3);
t310 = m(5) * t328 + t374 * mrSges(5,1) - t353 * mrSges(5,3) - t388 * t365 + t396 * t378 + t312;
t379 = t396 * mrSges(5,1) - t388 * mrSges(5,3);
t438 = -t419 * t320 + t424 * t321;
t311 = m(5) * t329 - t374 * mrSges(5,2) + t352 * mrSges(5,3) + t387 * t365 - t396 * t379 + t438;
t439 = -t420 * t310 + t425 * t311;
t304 = m(4) * t350 - t416 * mrSges(4,2) + t375 * mrSges(4,3) + t400 * t384 - t417 * t392 + t439;
t349 = t426 * t367 - t421 * t368;
t391 = -t417 * mrSges(4,2) + t400 * mrSges(4,3);
t341 = -t416 * pkin(3) - t415 * pkin(8) + t401 * t385 - t349;
t330 = -t352 * pkin(4) - t386 * pkin(9) + t388 * t380 + t341;
t435 = m(6) * t330 - t334 * mrSges(6,1) + t335 * mrSges(6,2) - t360 * t354 + t361 * t355;
t431 = -m(5) * t341 + t352 * mrSges(5,1) - t353 * mrSges(5,2) + t387 * t378 - t388 * t379 - t435;
t316 = m(4) * t349 + t416 * mrSges(4,1) - t376 * mrSges(4,3) - t401 * t384 + t417 * t391 + t431;
t299 = t421 * t304 + t426 * t316;
t306 = t425 * t310 + t420 * t311;
t440 = t426 * t304 - t421 * t316;
t344 = Ifges(6,4) * t361 + Ifges(6,2) * t360 + Ifges(6,6) * t394;
t345 = Ifges(6,1) * t361 + Ifges(6,4) * t360 + Ifges(6,5) * t394;
t434 = -mrSges(6,1) * t324 + mrSges(6,2) * t325 - Ifges(6,5) * t335 - Ifges(6,6) * t334 - Ifges(6,3) * t370 - t361 * t344 + t360 * t345;
t343 = Ifges(6,5) * t361 + Ifges(6,6) * t360 + Ifges(6,3) * t394;
t313 = -mrSges(6,1) * t330 + mrSges(6,3) * t325 + Ifges(6,4) * t335 + Ifges(6,2) * t334 + Ifges(6,6) * t370 - t361 * t343 + t394 * t345;
t314 = mrSges(6,2) * t330 - mrSges(6,3) * t324 + Ifges(6,1) * t335 + Ifges(6,4) * t334 + Ifges(6,5) * t370 + t360 * t343 - t394 * t344;
t356 = Ifges(5,5) * t388 + Ifges(5,6) * t387 + Ifges(5,3) * t396;
t358 = Ifges(5,1) * t388 + Ifges(5,4) * t387 + Ifges(5,5) * t396;
t298 = -mrSges(5,1) * t341 + mrSges(5,3) * t329 + Ifges(5,4) * t353 + Ifges(5,2) * t352 + Ifges(5,6) * t374 - pkin(4) * t435 + pkin(9) * t438 + t424 * t313 + t419 * t314 - t388 * t356 + t396 * t358;
t357 = Ifges(5,4) * t388 + Ifges(5,2) * t387 + Ifges(5,6) * t396;
t301 = mrSges(5,2) * t341 - mrSges(5,3) * t328 + Ifges(5,1) * t353 + Ifges(5,4) * t352 + Ifges(5,5) * t374 - pkin(9) * t312 - t419 * t313 + t424 * t314 + t387 * t356 - t396 * t357;
t382 = Ifges(4,4) * t401 + Ifges(4,2) * t400 + Ifges(4,6) * t417;
t383 = Ifges(4,1) * t401 + Ifges(4,4) * t400 + Ifges(4,5) * t417;
t433 = mrSges(4,1) * t349 - mrSges(4,2) * t350 + Ifges(4,5) * t376 + Ifges(4,6) * t375 + Ifges(4,3) * t416 + pkin(3) * t431 + pkin(8) * t439 + t425 * t298 + t420 * t301 + t401 * t382 - t400 * t383;
t432 = m(4) * t377 - t375 * mrSges(4,1) + t376 * mrSges(4,2) - t400 * t391 + t401 * t392 + t306;
t430 = mrSges(5,1) * t328 - mrSges(5,2) * t329 + Ifges(5,5) * t353 + Ifges(5,6) * t352 + Ifges(5,3) * t374 + pkin(4) * t312 + t388 * t357 - t387 * t358 - t434;
t410 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t443;
t409 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t444;
t406 = (-t427 * mrSges(3,1) + t422 * mrSges(3,2)) * qJD(1);
t402 = -t429 * pkin(6) + t436;
t399 = Ifges(3,5) * qJD(2) + (t422 * Ifges(3,1) + t427 * Ifges(3,4)) * qJD(1);
t398 = Ifges(3,6) * qJD(2) + (t422 * Ifges(3,4) + t427 * Ifges(3,2)) * qJD(1);
t389 = -t427 * g(3) - t445;
t381 = Ifges(4,5) * t401 + Ifges(4,6) * t400 + Ifges(4,3) * t417;
t296 = -mrSges(4,1) * t377 + mrSges(4,3) * t350 + Ifges(4,4) * t376 + Ifges(4,2) * t375 + Ifges(4,6) * t416 - pkin(3) * t306 - t401 * t381 + t417 * t383 - t430;
t295 = mrSges(4,2) * t377 - mrSges(4,3) * t349 + Ifges(4,1) * t376 + Ifges(4,4) * t375 + Ifges(4,5) * t416 - pkin(8) * t306 - t420 * t298 + t425 * t301 + t400 * t381 - t417 * t382;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t441 - mrSges(2,2) * t437 + t422 * (mrSges(3,2) * t402 - mrSges(3,3) * t389 + Ifges(3,1) * t407 + Ifges(3,4) * t408 + Ifges(3,5) * qJDD(2) - pkin(7) * t299 - qJD(2) * t398 + t426 * t295 - t421 * t296) + t427 * (-mrSges(3,1) * t402 + mrSges(3,3) * t390 + Ifges(3,4) * t407 + Ifges(3,2) * t408 + Ifges(3,6) * qJDD(2) - pkin(2) * t432 + pkin(7) * t440 + qJD(2) * t399 + t421 * t295 + t426 * t296) + pkin(1) * (-m(3) * t402 + t408 * mrSges(3,1) - t407 * mrSges(3,2) + (-t409 * t422 + t410 * t427) * qJD(1) - t432) + pkin(6) * (t427 * (m(3) * t390 - qJDD(2) * mrSges(3,2) + t408 * mrSges(3,3) - qJD(2) * t409 + t406 * t443 + t440) - t422 * (m(3) * t389 + qJDD(2) * mrSges(3,1) - t407 * mrSges(3,3) + qJD(2) * t410 - t406 * t444 + t299)); t433 + Ifges(3,5) * t407 + Ifges(3,6) * t408 + mrSges(3,1) * t389 - mrSges(3,2) * t390 + Ifges(3,3) * qJDD(2) + pkin(2) * t299 + (t422 * t398 - t427 * t399) * qJD(1); t433; t430; -t434;];
tauJ = t1;
