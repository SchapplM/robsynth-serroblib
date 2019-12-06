% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR2
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:46
% EndTime: 2019-12-05 18:27:49
% DurationCPUTime: 3.32s
% Computational Cost: add. (33114->275), mult. (79033->356), div. (0->0), fcn. (56439->10), ass. (0->109)
t429 = qJD(1) ^ 2;
t444 = pkin(2) * t429;
t424 = sin(qJ(1));
t428 = cos(qJ(1));
t435 = -g(1) * t428 - g(2) * t424;
t403 = -pkin(1) * t429 + qJDD(1) * pkin(6) + t435;
t423 = sin(qJ(2));
t443 = t403 * t423;
t427 = cos(qJ(2));
t440 = qJD(1) * qJD(2);
t406 = qJDD(1) * t423 + t427 * t440;
t372 = qJDD(2) * pkin(2) - qJ(3) * t406 - t443 + (qJ(3) * t440 + t423 * t444 - g(3)) * t427;
t389 = -g(3) * t423 + t427 * t403;
t407 = qJDD(1) * t427 - t423 * t440;
t442 = qJD(1) * t423;
t408 = qJD(2) * pkin(2) - qJ(3) * t442;
t418 = t427 ^ 2;
t373 = qJ(3) * t407 - qJD(2) * t408 - t418 * t444 + t389;
t419 = sin(pkin(9));
t420 = cos(pkin(9));
t398 = (t419 * t427 + t420 * t423) * qJD(1);
t352 = -0.2e1 * qJD(3) * t398 + t420 * t372 - t373 * t419;
t387 = t406 * t420 + t407 * t419;
t397 = (-t419 * t423 + t420 * t427) * qJD(1);
t342 = (qJD(2) * t397 - t387) * pkin(7) + (t397 * t398 + qJDD(2)) * pkin(3) + t352;
t353 = 0.2e1 * qJD(3) * t397 + t419 * t372 + t420 * t373;
t386 = -t406 * t419 + t407 * t420;
t392 = qJD(2) * pkin(3) - pkin(7) * t398;
t396 = t397 ^ 2;
t344 = -pkin(3) * t396 + pkin(7) * t386 - qJD(2) * t392 + t353;
t422 = sin(qJ(4));
t426 = cos(qJ(4));
t330 = t426 * t342 - t344 * t422;
t382 = t397 * t426 - t398 * t422;
t359 = qJD(4) * t382 + t386 * t422 + t387 * t426;
t383 = t397 * t422 + t398 * t426;
t416 = qJDD(2) + qJDD(4);
t417 = qJD(2) + qJD(4);
t327 = (t382 * t417 - t359) * pkin(8) + (t382 * t383 + t416) * pkin(4) + t330;
t331 = t422 * t342 + t426 * t344;
t358 = -qJD(4) * t383 + t386 * t426 - t387 * t422;
t377 = pkin(4) * t417 - pkin(8) * t383;
t378 = t382 ^ 2;
t328 = -pkin(4) * t378 + pkin(8) * t358 - t377 * t417 + t331;
t421 = sin(qJ(5));
t425 = cos(qJ(5));
t325 = t327 * t425 - t328 * t421;
t366 = t382 * t425 - t383 * t421;
t338 = qJD(5) * t366 + t358 * t421 + t359 * t425;
t367 = t382 * t421 + t383 * t425;
t349 = -mrSges(6,1) * t366 + mrSges(6,2) * t367;
t414 = qJD(5) + t417;
t360 = -mrSges(6,2) * t414 + mrSges(6,3) * t366;
t413 = qJDD(5) + t416;
t322 = m(6) * t325 + mrSges(6,1) * t413 - mrSges(6,3) * t338 - t349 * t367 + t360 * t414;
t326 = t327 * t421 + t328 * t425;
t337 = -qJD(5) * t367 + t358 * t425 - t359 * t421;
t361 = mrSges(6,1) * t414 - mrSges(6,3) * t367;
t323 = m(6) * t326 - mrSges(6,2) * t413 + mrSges(6,3) * t337 + t349 * t366 - t361 * t414;
t315 = t425 * t322 + t421 * t323;
t368 = -mrSges(5,1) * t382 + mrSges(5,2) * t383;
t375 = -mrSges(5,2) * t417 + mrSges(5,3) * t382;
t312 = m(5) * t330 + mrSges(5,1) * t416 - mrSges(5,3) * t359 - t368 * t383 + t375 * t417 + t315;
t376 = mrSges(5,1) * t417 - mrSges(5,3) * t383;
t436 = -t322 * t421 + t425 * t323;
t313 = m(5) * t331 - mrSges(5,2) * t416 + mrSges(5,3) * t358 + t368 * t382 - t376 * t417 + t436;
t308 = t426 * t312 + t422 * t313;
t385 = -mrSges(4,1) * t397 + mrSges(4,2) * t398;
t390 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t397;
t306 = m(4) * t352 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t387 + qJD(2) * t390 - t385 * t398 + t308;
t391 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t398;
t437 = -t312 * t422 + t426 * t313;
t307 = m(4) * t353 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t386 - qJD(2) * t391 + t385 * t397 + t437;
t300 = t420 * t306 + t419 * t307;
t441 = t427 * qJD(1);
t439 = g(1) * t424 - t428 * g(2);
t438 = -t306 * t419 + t420 * t307;
t434 = -qJDD(1) * pkin(1) - t439;
t374 = -pkin(2) * t407 + qJDD(3) + t408 * t442 + (-qJ(3) * t418 - pkin(6)) * t429 + t434;
t351 = -pkin(3) * t386 - pkin(7) * t396 + t398 * t392 + t374;
t333 = -pkin(4) * t358 - pkin(8) * t378 + t377 * t383 + t351;
t433 = -m(6) * t333 + t337 * mrSges(6,1) - t338 * mrSges(6,2) + t366 * t360 - t367 * t361;
t346 = Ifges(6,4) * t367 + Ifges(6,2) * t366 + Ifges(6,6) * t414;
t347 = Ifges(6,1) * t367 + Ifges(6,4) * t366 + Ifges(6,5) * t414;
t432 = mrSges(6,1) * t325 - mrSges(6,2) * t326 + Ifges(6,5) * t338 + Ifges(6,6) * t337 + Ifges(6,3) * t413 + t367 * t346 - t366 * t347;
t431 = -m(5) * t351 + t358 * mrSges(5,1) - t359 * mrSges(5,2) + t382 * t375 - t383 * t376 + t433;
t363 = Ifges(5,4) * t383 + Ifges(5,2) * t382 + Ifges(5,6) * t417;
t364 = Ifges(5,1) * t383 + Ifges(5,4) * t382 + Ifges(5,5) * t417;
t430 = mrSges(5,1) * t330 - mrSges(5,2) * t331 + Ifges(5,5) * t359 + Ifges(5,6) * t358 + Ifges(5,3) * t416 + pkin(4) * t315 + t383 * t363 - t382 * t364 + t432;
t318 = m(4) * t374 - t386 * mrSges(4,1) + t387 * mrSges(4,2) - t397 * t390 + t398 * t391 - t431;
t410 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t441;
t409 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t442;
t405 = (-mrSges(3,1) * t427 + mrSges(3,2) * t423) * qJD(1);
t402 = -pkin(6) * t429 + t434;
t401 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t423 + Ifges(3,4) * t427) * qJD(1);
t400 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t423 + Ifges(3,2) * t427) * qJD(1);
t388 = -g(3) * t427 - t443;
t381 = Ifges(4,1) * t398 + Ifges(4,4) * t397 + Ifges(4,5) * qJD(2);
t380 = Ifges(4,4) * t398 + Ifges(4,2) * t397 + Ifges(4,6) * qJD(2);
t379 = Ifges(4,5) * t398 + Ifges(4,6) * t397 + Ifges(4,3) * qJD(2);
t362 = Ifges(5,5) * t383 + Ifges(5,6) * t382 + Ifges(5,3) * t417;
t345 = Ifges(6,5) * t367 + Ifges(6,6) * t366 + Ifges(6,3) * t414;
t317 = mrSges(6,2) * t333 - mrSges(6,3) * t325 + Ifges(6,1) * t338 + Ifges(6,4) * t337 + Ifges(6,5) * t413 + t345 * t366 - t346 * t414;
t316 = -mrSges(6,1) * t333 + mrSges(6,3) * t326 + Ifges(6,4) * t338 + Ifges(6,2) * t337 + Ifges(6,6) * t413 - t345 * t367 + t347 * t414;
t302 = mrSges(5,2) * t351 - mrSges(5,3) * t330 + Ifges(5,1) * t359 + Ifges(5,4) * t358 + Ifges(5,5) * t416 - pkin(8) * t315 - t316 * t421 + t317 * t425 + t362 * t382 - t363 * t417;
t301 = -mrSges(5,1) * t351 + mrSges(5,3) * t331 + Ifges(5,4) * t359 + Ifges(5,2) * t358 + Ifges(5,6) * t416 + pkin(4) * t433 + pkin(8) * t436 + t425 * t316 + t421 * t317 - t383 * t362 + t417 * t364;
t299 = mrSges(4,2) * t374 - mrSges(4,3) * t352 + Ifges(4,1) * t387 + Ifges(4,4) * t386 + Ifges(4,5) * qJDD(2) - pkin(7) * t308 - qJD(2) * t380 - t301 * t422 + t302 * t426 + t379 * t397;
t298 = -mrSges(4,1) * t374 + mrSges(4,3) * t353 + Ifges(4,4) * t387 + Ifges(4,2) * t386 + Ifges(4,6) * qJDD(2) + pkin(3) * t431 + pkin(7) * t437 + qJD(2) * t381 + t426 * t301 + t422 * t302 - t398 * t379;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t439 - mrSges(2,2) * t435 + t423 * (mrSges(3,2) * t402 - mrSges(3,3) * t388 + Ifges(3,1) * t406 + Ifges(3,4) * t407 + Ifges(3,5) * qJDD(2) - qJ(3) * t300 - qJD(2) * t400 - t419 * t298 + t420 * t299) + t427 * (-mrSges(3,1) * t402 + mrSges(3,3) * t389 + Ifges(3,4) * t406 + Ifges(3,2) * t407 + Ifges(3,6) * qJDD(2) - pkin(2) * t318 + qJ(3) * t438 + qJD(2) * t401 + t420 * t298 + t419 * t299) + pkin(1) * (-t406 * mrSges(3,2) + t407 * mrSges(3,1) - m(3) * t402 - t318 + (-t423 * t409 + t427 * t410) * qJD(1)) + pkin(6) * (t427 * (m(3) * t389 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t407 - qJD(2) * t409 + t405 * t441 + t438) - t423 * (m(3) * t388 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t406 + qJD(2) * t410 - t405 * t442 + t300)); Ifges(3,6) * t407 - t397 * t381 + t398 * t380 + Ifges(3,5) * t406 + Ifges(4,5) * t387 + mrSges(3,1) * t388 - mrSges(3,2) * t389 + Ifges(4,6) * t386 + mrSges(4,1) * t352 - mrSges(4,2) * t353 + pkin(3) * t308 + t430 + pkin(2) * t300 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t423 * t400 - t427 * t401) * qJD(1); t318; t430; t432;];
tauJ = t1;
