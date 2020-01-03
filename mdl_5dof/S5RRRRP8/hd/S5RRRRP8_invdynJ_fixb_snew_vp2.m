% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:00:01
% EndTime: 2019-12-31 22:00:06
% DurationCPUTime: 2.26s
% Computational Cost: add. (15234->252), mult. (30520->307), div. (0->0), fcn. (20495->8), ass. (0->100)
t411 = sin(qJ(3));
t415 = cos(qJ(3));
t412 = sin(qJ(2));
t434 = t412 * qJD(1);
t398 = t411 * qJD(2) + t415 * t434;
t416 = cos(qJ(2));
t432 = qJD(1) * qJD(2);
t429 = t416 * t432;
t401 = t412 * qJDD(1) + t429;
t372 = -t398 * qJD(3) + t415 * qJDD(2) - t411 * t401;
t397 = t415 * qJD(2) - t411 * t434;
t373 = t397 * qJD(3) + t411 * qJDD(2) + t415 * t401;
t410 = sin(qJ(4));
t414 = cos(qJ(4));
t375 = t414 * t397 - t410 * t398;
t339 = t375 * qJD(4) + t410 * t372 + t414 * t373;
t376 = t410 * t397 + t414 * t398;
t352 = -t375 * mrSges(6,1) + t376 * mrSges(6,2);
t419 = qJD(1) ^ 2;
t413 = sin(qJ(1));
t417 = cos(qJ(1));
t428 = t413 * g(1) - t417 * g(2);
t392 = -qJDD(1) * pkin(1) - t419 * pkin(6) - t428;
t407 = t412 * t432;
t402 = t416 * qJDD(1) - t407;
t356 = (-t401 - t429) * pkin(7) + (-t402 + t407) * pkin(2) + t392;
t424 = -t417 * g(1) - t413 * g(2);
t393 = -t419 * pkin(1) + qJDD(1) * pkin(6) + t424;
t381 = -t412 * g(3) + t416 * t393;
t400 = (-t416 * pkin(2) - t412 * pkin(7)) * qJD(1);
t418 = qJD(2) ^ 2;
t433 = t416 * qJD(1);
t360 = -t418 * pkin(2) + qJDD(2) * pkin(7) + t400 * t433 + t381;
t340 = t415 * t356 - t411 * t360;
t396 = qJDD(3) - t402;
t406 = qJD(3) - t433;
t324 = (t397 * t406 - t373) * pkin(8) + (t397 * t398 + t396) * pkin(3) + t340;
t341 = t411 * t356 + t415 * t360;
t382 = t406 * pkin(3) - t398 * pkin(8);
t395 = t397 ^ 2;
t326 = -t395 * pkin(3) + t372 * pkin(8) - t406 * t382 + t341;
t318 = t414 * t324 - t410 * t326;
t394 = qJDD(4) + t396;
t405 = qJD(4) + t406;
t314 = -0.2e1 * qJD(5) * t376 + (t375 * t405 - t339) * qJ(5) + (t375 * t376 + t394) * pkin(4) + t318;
t361 = -t405 * mrSges(6,2) + t375 * mrSges(6,3);
t431 = m(6) * t314 + t394 * mrSges(6,1) + t405 * t361;
t310 = -t339 * mrSges(6,3) - t376 * t352 + t431;
t319 = t410 * t324 + t414 * t326;
t338 = -t376 * qJD(4) + t414 * t372 - t410 * t373;
t363 = t405 * pkin(4) - t376 * qJ(5);
t374 = t375 ^ 2;
t316 = -t374 * pkin(4) + t338 * qJ(5) + 0.2e1 * qJD(5) * t375 - t405 * t363 + t319;
t437 = Ifges(5,4) + Ifges(6,4);
t442 = Ifges(5,5) + Ifges(6,5);
t443 = Ifges(5,1) + Ifges(6,1);
t438 = t437 * t375 + t443 * t376 + t442 * t405;
t441 = Ifges(5,6) + Ifges(6,6);
t445 = Ifges(5,2) + Ifges(6,2);
t439 = t445 * t375 + t437 * t376 + t441 * t405;
t440 = Ifges(5,3) + Ifges(6,3);
t446 = mrSges(5,1) * t318 + mrSges(6,1) * t314 - mrSges(5,2) * t319 - mrSges(6,2) * t316 + pkin(4) * t310 + t441 * t338 + t442 * t339 - t375 * t438 + t439 * t376 + t440 * t394;
t353 = -t375 * mrSges(5,1) + t376 * mrSges(5,2);
t362 = -t405 * mrSges(5,2) + t375 * mrSges(5,3);
t305 = m(5) * t318 + t394 * mrSges(5,1) + t405 * t362 + (-t352 - t353) * t376 + (-mrSges(5,3) - mrSges(6,3)) * t339 + t431;
t364 = t405 * mrSges(6,1) - t376 * mrSges(6,3);
t365 = t405 * mrSges(5,1) - t376 * mrSges(5,3);
t430 = m(6) * t316 + t338 * mrSges(6,3) + t375 * t352;
t308 = m(5) * t319 + t338 * mrSges(5,3) + t375 * t353 + (-t364 - t365) * t405 + (-mrSges(5,2) - mrSges(6,2)) * t394 + t430;
t303 = t414 * t305 + t410 * t308;
t367 = Ifges(4,4) * t398 + Ifges(4,2) * t397 + Ifges(4,6) * t406;
t368 = Ifges(4,1) * t398 + Ifges(4,4) * t397 + Ifges(4,5) * t406;
t444 = mrSges(4,1) * t340 - mrSges(4,2) * t341 + Ifges(4,5) * t373 + Ifges(4,6) * t372 + Ifges(4,3) * t396 + pkin(3) * t303 + t398 * t367 - t397 * t368 + t446;
t436 = -t441 * t375 - t442 * t376 - t440 * t405;
t380 = -t416 * g(3) - t412 * t393;
t377 = -t397 * mrSges(4,1) + t398 * mrSges(4,2);
t378 = -t406 * mrSges(4,2) + t397 * mrSges(4,3);
t300 = m(4) * t340 + t396 * mrSges(4,1) - t373 * mrSges(4,3) - t398 * t377 + t406 * t378 + t303;
t379 = t406 * mrSges(4,1) - t398 * mrSges(4,3);
t425 = -t410 * t305 + t414 * t308;
t301 = m(4) * t341 - t396 * mrSges(4,2) + t372 * mrSges(4,3) + t397 * t377 - t406 * t379 + t425;
t426 = -t411 * t300 + t415 * t301;
t359 = -qJDD(2) * pkin(2) - t418 * pkin(7) + t400 * t434 - t380;
t327 = -t372 * pkin(3) - t395 * pkin(8) + t398 * t382 + t359;
t321 = -t338 * pkin(4) - t374 * qJ(5) + t376 * t363 + qJDD(5) + t327;
t311 = m(6) * t321 - t338 * mrSges(6,1) + t339 * mrSges(6,2) - t375 * t361 + t376 * t364;
t297 = t415 * t300 + t411 * t301;
t423 = m(5) * t327 - t338 * mrSges(5,1) + t339 * mrSges(5,2) - t375 * t362 + t376 * t365 + t311;
t421 = -m(4) * t359 + t372 * mrSges(4,1) - t373 * mrSges(4,2) + t397 * t378 - t398 * t379 - t423;
t404 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t433;
t403 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t434;
t399 = (-t416 * mrSges(3,1) + t412 * mrSges(3,2)) * qJD(1);
t391 = Ifges(3,5) * qJD(2) + (t412 * Ifges(3,1) + t416 * Ifges(3,4)) * qJD(1);
t390 = Ifges(3,6) * qJD(2) + (t412 * Ifges(3,4) + t416 * Ifges(3,2)) * qJD(1);
t366 = Ifges(4,5) * t398 + Ifges(4,6) * t397 + Ifges(4,3) * t406;
t302 = mrSges(5,2) * t327 + mrSges(6,2) * t321 - mrSges(5,3) * t318 - mrSges(6,3) * t314 - qJ(5) * t310 + t437 * t338 + t339 * t443 - t436 * t375 + t442 * t394 - t439 * t405;
t298 = -mrSges(5,1) * t327 + mrSges(5,3) * t319 - mrSges(6,1) * t321 + mrSges(6,3) * t316 - pkin(4) * t311 + qJ(5) * t430 + (-qJ(5) * t364 + t438) * t405 + (-qJ(5) * mrSges(6,2) + t441) * t394 + t436 * t376 + t437 * t339 + t445 * t338;
t296 = mrSges(4,2) * t359 - mrSges(4,3) * t340 + Ifges(4,1) * t373 + Ifges(4,4) * t372 + Ifges(4,5) * t396 - pkin(8) * t303 - t410 * t298 + t414 * t302 + t397 * t366 - t406 * t367;
t295 = -mrSges(4,1) * t359 + mrSges(4,3) * t341 + Ifges(4,4) * t373 + Ifges(4,2) * t372 + Ifges(4,6) * t396 - pkin(3) * t423 + pkin(8) * t425 + t414 * t298 + t410 * t302 - t398 * t366 + t406 * t368;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t428 - mrSges(2,2) * t424 + t412 * (mrSges(3,2) * t392 - mrSges(3,3) * t380 + Ifges(3,1) * t401 + Ifges(3,4) * t402 + Ifges(3,5) * qJDD(2) - pkin(7) * t297 - qJD(2) * t390 - t411 * t295 + t415 * t296) + t416 * (-mrSges(3,1) * t392 + mrSges(3,3) * t381 + Ifges(3,4) * t401 + Ifges(3,2) * t402 + Ifges(3,6) * qJDD(2) - pkin(2) * t297 + qJD(2) * t391 - t444) + pkin(1) * (-m(3) * t392 + t402 * mrSges(3,1) - t401 * mrSges(3,2) + (-t403 * t412 + t404 * t416) * qJD(1) - t297) + pkin(6) * (t416 * (m(3) * t381 - qJDD(2) * mrSges(3,2) + t402 * mrSges(3,3) - qJD(2) * t403 + t399 * t433 + t426) - t412 * (m(3) * t380 + qJDD(2) * mrSges(3,1) - t401 * mrSges(3,3) + qJD(2) * t404 - t399 * t434 + t421)); Ifges(3,5) * t401 + Ifges(3,6) * t402 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t380 - mrSges(3,2) * t381 + t411 * t296 + t415 * t295 + pkin(2) * t421 + pkin(7) * t426 + (t412 * t390 - t416 * t391) * qJD(1); t444; t446; t311;];
tauJ = t1;
