% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR7
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:51
% EndTime: 2019-12-31 22:21:54
% DurationCPUTime: 2.94s
% Computational Cost: add. (30865->276), mult. (67044->355), div. (0->0), fcn. (48360->10), ass. (0->112)
t424 = qJD(1) ^ 2;
t441 = pkin(2) * t424;
t418 = sin(qJ(1));
t423 = cos(qJ(1));
t432 = -t423 * g(1) - t418 * g(2);
t397 = -t424 * pkin(1) + qJDD(1) * pkin(6) + t432;
t417 = sin(qJ(2));
t440 = t417 * t397;
t422 = cos(qJ(2));
t437 = qJD(1) * qJD(2);
t400 = t417 * qJDD(1) + t422 * t437;
t365 = qJDD(2) * pkin(2) - t400 * pkin(7) - t440 + (pkin(7) * t437 + t417 * t441 - g(3)) * t422;
t385 = -t417 * g(3) + t422 * t397;
t401 = t422 * qJDD(1) - t417 * t437;
t439 = qJD(1) * t417;
t404 = qJD(2) * pkin(2) - pkin(7) * t439;
t413 = t422 ^ 2;
t366 = t401 * pkin(7) - qJD(2) * t404 - t413 * t441 + t385;
t416 = sin(qJ(3));
t421 = cos(qJ(3));
t350 = t421 * t365 - t416 * t366;
t394 = (-t417 * t416 + t422 * t421) * qJD(1);
t372 = t394 * qJD(3) + t421 * t400 + t416 * t401;
t395 = (t422 * t416 + t417 * t421) * qJD(1);
t411 = qJDD(2) + qJDD(3);
t412 = qJD(2) + qJD(3);
t330 = (t394 * t412 - t372) * pkin(8) + (t394 * t395 + t411) * pkin(3) + t350;
t351 = t416 * t365 + t421 * t366;
t371 = -t395 * qJD(3) - t416 * t400 + t421 * t401;
t388 = t412 * pkin(3) - t395 * pkin(8);
t390 = t394 ^ 2;
t333 = -t390 * pkin(3) + t371 * pkin(8) - t412 * t388 + t351;
t415 = sin(qJ(4));
t420 = cos(qJ(4));
t328 = t415 * t330 + t420 * t333;
t382 = t415 * t394 + t420 * t395;
t346 = -t382 * qJD(4) + t420 * t371 - t415 * t372;
t381 = t420 * t394 - t415 * t395;
t359 = -t381 * mrSges(5,1) + t382 * mrSges(5,2);
t409 = qJD(4) + t412;
t375 = t409 * mrSges(5,1) - t382 * mrSges(5,3);
t408 = qJDD(4) + t411;
t360 = -t381 * pkin(4) - t382 * pkin(9);
t407 = t409 ^ 2;
t324 = -t407 * pkin(4) + t408 * pkin(9) + t381 * t360 + t328;
t436 = t418 * g(1) - t423 * g(2);
t431 = -qJDD(1) * pkin(1) - t436;
t373 = -t401 * pkin(2) + t404 * t439 + (-pkin(7) * t413 - pkin(6)) * t424 + t431;
t337 = -t371 * pkin(3) - t390 * pkin(8) + t395 * t388 + t373;
t347 = t381 * qJD(4) + t415 * t371 + t420 * t372;
t325 = (-t381 * t409 - t347) * pkin(9) + (t382 * t409 - t346) * pkin(4) + t337;
t414 = sin(qJ(5));
t419 = cos(qJ(5));
t321 = -t414 * t324 + t419 * t325;
t367 = -t414 * t382 + t419 * t409;
t335 = t367 * qJD(5) + t419 * t347 + t414 * t408;
t344 = qJDD(5) - t346;
t368 = t419 * t382 + t414 * t409;
t352 = -t367 * mrSges(6,1) + t368 * mrSges(6,2);
t379 = qJD(5) - t381;
t353 = -t379 * mrSges(6,2) + t367 * mrSges(6,3);
t318 = m(6) * t321 + t344 * mrSges(6,1) - t335 * mrSges(6,3) - t368 * t352 + t379 * t353;
t322 = t419 * t324 + t414 * t325;
t334 = -t368 * qJD(5) - t414 * t347 + t419 * t408;
t354 = t379 * mrSges(6,1) - t368 * mrSges(6,3);
t319 = m(6) * t322 - t344 * mrSges(6,2) + t334 * mrSges(6,3) + t367 * t352 - t379 * t354;
t433 = -t414 * t318 + t419 * t319;
t306 = m(5) * t328 - t408 * mrSges(5,2) + t346 * mrSges(5,3) + t381 * t359 - t409 * t375 + t433;
t327 = t420 * t330 - t415 * t333;
t374 = -t409 * mrSges(5,2) + t381 * mrSges(5,3);
t323 = -t408 * pkin(4) - t407 * pkin(9) + t382 * t360 - t327;
t430 = -m(6) * t323 + t334 * mrSges(6,1) - t335 * mrSges(6,2) + t367 * t353 - t368 * t354;
t314 = m(5) * t327 + t408 * mrSges(5,1) - t347 * mrSges(5,3) - t382 * t359 + t409 * t374 + t430;
t303 = t415 * t306 + t420 * t314;
t383 = -t394 * mrSges(4,1) + t395 * mrSges(4,2);
t386 = -t412 * mrSges(4,2) + t394 * mrSges(4,3);
t300 = m(4) * t350 + t411 * mrSges(4,1) - t372 * mrSges(4,3) - t395 * t383 + t412 * t386 + t303;
t387 = t412 * mrSges(4,1) - t395 * mrSges(4,3);
t434 = t420 * t306 - t415 * t314;
t301 = m(4) * t351 - t411 * mrSges(4,2) + t371 * mrSges(4,3) + t394 * t383 - t412 * t387 + t434;
t294 = t421 * t300 + t416 * t301;
t308 = t419 * t318 + t414 * t319;
t438 = qJD(1) * t422;
t435 = -t416 * t300 + t421 * t301;
t429 = -m(5) * t337 + t346 * mrSges(5,1) - t347 * mrSges(5,2) + t381 * t374 - t382 * t375 - t308;
t338 = Ifges(6,5) * t368 + Ifges(6,6) * t367 + Ifges(6,3) * t379;
t340 = Ifges(6,1) * t368 + Ifges(6,4) * t367 + Ifges(6,5) * t379;
t311 = -mrSges(6,1) * t323 + mrSges(6,3) * t322 + Ifges(6,4) * t335 + Ifges(6,2) * t334 + Ifges(6,6) * t344 - t368 * t338 + t379 * t340;
t339 = Ifges(6,4) * t368 + Ifges(6,2) * t367 + Ifges(6,6) * t379;
t312 = mrSges(6,2) * t323 - mrSges(6,3) * t321 + Ifges(6,1) * t335 + Ifges(6,4) * t334 + Ifges(6,5) * t344 + t367 * t338 - t379 * t339;
t356 = Ifges(5,4) * t382 + Ifges(5,2) * t381 + Ifges(5,6) * t409;
t357 = Ifges(5,1) * t382 + Ifges(5,4) * t381 + Ifges(5,5) * t409;
t428 = mrSges(5,1) * t327 - mrSges(5,2) * t328 + Ifges(5,5) * t347 + Ifges(5,6) * t346 + Ifges(5,3) * t408 + pkin(4) * t430 + pkin(9) * t433 + t419 * t311 + t414 * t312 + t382 * t356 - t381 * t357;
t427 = mrSges(6,1) * t321 - mrSges(6,2) * t322 + Ifges(6,5) * t335 + Ifges(6,6) * t334 + Ifges(6,3) * t344 + t368 * t339 - t367 * t340;
t426 = m(4) * t373 - t371 * mrSges(4,1) + t372 * mrSges(4,2) - t394 * t386 + t395 * t387 - t429;
t377 = Ifges(4,4) * t395 + Ifges(4,2) * t394 + Ifges(4,6) * t412;
t378 = Ifges(4,1) * t395 + Ifges(4,4) * t394 + Ifges(4,5) * t412;
t425 = mrSges(4,1) * t350 - mrSges(4,2) * t351 + Ifges(4,5) * t372 + Ifges(4,6) * t371 + Ifges(4,3) * t411 + pkin(3) * t303 + t395 * t377 - t394 * t378 + t428;
t403 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t438;
t402 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t439;
t399 = (-t422 * mrSges(3,1) + t417 * mrSges(3,2)) * qJD(1);
t396 = -t424 * pkin(6) + t431;
t393 = Ifges(3,5) * qJD(2) + (t417 * Ifges(3,1) + t422 * Ifges(3,4)) * qJD(1);
t392 = Ifges(3,6) * qJD(2) + (t417 * Ifges(3,4) + t422 * Ifges(3,2)) * qJD(1);
t384 = -t422 * g(3) - t440;
t376 = Ifges(4,5) * t395 + Ifges(4,6) * t394 + Ifges(4,3) * t412;
t355 = Ifges(5,5) * t382 + Ifges(5,6) * t381 + Ifges(5,3) * t409;
t296 = -mrSges(5,1) * t337 + mrSges(5,3) * t328 + Ifges(5,4) * t347 + Ifges(5,2) * t346 + Ifges(5,6) * t408 - pkin(4) * t308 - t382 * t355 + t409 * t357 - t427;
t295 = mrSges(5,2) * t337 - mrSges(5,3) * t327 + Ifges(5,1) * t347 + Ifges(5,4) * t346 + Ifges(5,5) * t408 - pkin(9) * t308 - t414 * t311 + t419 * t312 + t381 * t355 - t409 * t356;
t293 = mrSges(4,2) * t373 - mrSges(4,3) * t350 + Ifges(4,1) * t372 + Ifges(4,4) * t371 + Ifges(4,5) * t411 - pkin(8) * t303 + t420 * t295 - t415 * t296 + t394 * t376 - t412 * t377;
t292 = -mrSges(4,1) * t373 + mrSges(4,3) * t351 + Ifges(4,4) * t372 + Ifges(4,2) * t371 + Ifges(4,6) * t411 + pkin(3) * t429 + pkin(8) * t434 + t415 * t295 + t420 * t296 - t395 * t376 + t412 * t378;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t436 - mrSges(2,2) * t432 + t417 * (mrSges(3,2) * t396 - mrSges(3,3) * t384 + Ifges(3,1) * t400 + Ifges(3,4) * t401 + Ifges(3,5) * qJDD(2) - pkin(7) * t294 - qJD(2) * t392 - t416 * t292 + t421 * t293) + t422 * (-mrSges(3,1) * t396 + mrSges(3,3) * t385 + Ifges(3,4) * t400 + Ifges(3,2) * t401 + Ifges(3,6) * qJDD(2) - pkin(2) * t426 + pkin(7) * t435 + qJD(2) * t393 + t421 * t292 + t416 * t293) + pkin(1) * (t401 * mrSges(3,1) - m(3) * t396 - t400 * mrSges(3,2) + (-t402 * t417 + t403 * t422) * qJD(1) - t426) + pkin(6) * (t422 * (m(3) * t385 - qJDD(2) * mrSges(3,2) + t401 * mrSges(3,3) - qJD(2) * t402 + t399 * t438 + t435) - t417 * (m(3) * t384 + qJDD(2) * mrSges(3,1) - t400 * mrSges(3,3) + qJD(2) * t403 - t399 * t439 + t294)); Ifges(3,3) * qJDD(2) + Ifges(3,6) * t401 + Ifges(3,5) * t400 + mrSges(3,1) * t384 - mrSges(3,2) * t385 + pkin(2) * t294 + t425 + (t417 * t392 - t422 * t393) * qJD(1); t425; t428; t427;];
tauJ = t1;
