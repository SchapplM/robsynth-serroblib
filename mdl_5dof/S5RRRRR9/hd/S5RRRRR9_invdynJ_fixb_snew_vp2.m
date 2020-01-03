% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR9
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:28:18
% EndTime: 2019-12-31 22:28:26
% DurationCPUTime: 3.44s
% Computational Cost: add. (37997->276), mult. (77203->349), div. (0->0), fcn. (54238->10), ass. (0->111)
t420 = qJD(1) ^ 2;
t413 = sin(qJ(1));
t418 = cos(qJ(1));
t431 = t413 * g(1) - t418 * g(2);
t390 = -qJDD(1) * pkin(1) - t420 * pkin(6) - t431;
t412 = sin(qJ(2));
t417 = cos(qJ(2));
t433 = qJD(1) * qJD(2);
t432 = t417 * t433;
t399 = qJDD(1) * t412 + t432;
t406 = t412 * t433;
t400 = qJDD(1) * t417 - t406;
t357 = (-t399 - t432) * pkin(7) + (-t400 + t406) * pkin(2) + t390;
t427 = -g(1) * t418 - g(2) * t413;
t391 = -pkin(1) * t420 + qJDD(1) * pkin(6) + t427;
t379 = -g(3) * t412 + t417 * t391;
t398 = (-pkin(2) * t417 - pkin(7) * t412) * qJD(1);
t419 = qJD(2) ^ 2;
t434 = t417 * qJD(1);
t360 = -pkin(2) * t419 + qJDD(2) * pkin(7) + t398 * t434 + t379;
t411 = sin(qJ(3));
t416 = cos(qJ(3));
t344 = t416 * t357 - t411 * t360;
t435 = t412 * qJD(1);
t395 = qJD(2) * t416 - t411 * t435;
t371 = qJD(3) * t395 + qJDD(2) * t411 + t399 * t416;
t394 = qJDD(3) - t400;
t396 = qJD(2) * t411 + t416 * t435;
t405 = qJD(3) - t434;
t328 = (t395 * t405 - t371) * pkin(8) + (t395 * t396 + t394) * pkin(3) + t344;
t345 = t411 * t357 + t416 * t360;
t370 = -qJD(3) * t396 + qJDD(2) * t416 - t399 * t411;
t380 = pkin(3) * t405 - pkin(8) * t396;
t393 = t395 ^ 2;
t330 = -pkin(3) * t393 + pkin(8) * t370 - t380 * t405 + t345;
t410 = sin(qJ(4));
t415 = cos(qJ(4));
t315 = t415 * t328 - t410 * t330;
t373 = t395 * t415 - t396 * t410;
t343 = qJD(4) * t373 + t370 * t410 + t371 * t415;
t374 = t395 * t410 + t396 * t415;
t392 = qJDD(4) + t394;
t404 = qJD(4) + t405;
t312 = (t373 * t404 - t343) * pkin(9) + (t373 * t374 + t392) * pkin(4) + t315;
t316 = t410 * t328 + t415 * t330;
t342 = -qJD(4) * t374 + t370 * t415 - t371 * t410;
t363 = pkin(4) * t404 - pkin(9) * t374;
t372 = t373 ^ 2;
t313 = -pkin(4) * t372 + pkin(9) * t342 - t363 * t404 + t316;
t409 = sin(qJ(5));
t414 = cos(qJ(5));
t310 = t312 * t414 - t313 * t409;
t352 = t373 * t414 - t374 * t409;
t324 = qJD(5) * t352 + t342 * t409 + t343 * t414;
t353 = t373 * t409 + t374 * t414;
t336 = -mrSges(6,1) * t352 + mrSges(6,2) * t353;
t401 = qJD(5) + t404;
t346 = -mrSges(6,2) * t401 + mrSges(6,3) * t352;
t386 = qJDD(5) + t392;
t307 = m(6) * t310 + mrSges(6,1) * t386 - mrSges(6,3) * t324 - t336 * t353 + t346 * t401;
t311 = t312 * t409 + t313 * t414;
t323 = -qJD(5) * t353 + t342 * t414 - t343 * t409;
t347 = mrSges(6,1) * t401 - mrSges(6,3) * t353;
t308 = m(6) * t311 - mrSges(6,2) * t386 + mrSges(6,3) * t323 + t336 * t352 - t347 * t401;
t301 = t414 * t307 + t409 * t308;
t354 = -mrSges(5,1) * t373 + mrSges(5,2) * t374;
t361 = -mrSges(5,2) * t404 + mrSges(5,3) * t373;
t298 = m(5) * t315 + mrSges(5,1) * t392 - mrSges(5,3) * t343 - t354 * t374 + t361 * t404 + t301;
t362 = mrSges(5,1) * t404 - mrSges(5,3) * t374;
t428 = -t307 * t409 + t414 * t308;
t299 = m(5) * t316 - mrSges(5,2) * t392 + mrSges(5,3) * t342 + t354 * t373 - t362 * t404 + t428;
t294 = t298 * t415 + t299 * t410;
t378 = -t417 * g(3) - t412 * t391;
t375 = -mrSges(4,1) * t395 + mrSges(4,2) * t396;
t376 = -mrSges(4,2) * t405 + mrSges(4,3) * t395;
t292 = m(4) * t344 + mrSges(4,1) * t394 - mrSges(4,3) * t371 - t375 * t396 + t376 * t405 + t294;
t377 = mrSges(4,1) * t405 - mrSges(4,3) * t396;
t429 = -t298 * t410 + t299 * t415;
t293 = m(4) * t345 - mrSges(4,2) * t394 + mrSges(4,3) * t370 + t375 * t395 - t377 * t405 + t429;
t430 = -t292 * t411 + t293 * t416;
t359 = -qJDD(2) * pkin(2) - pkin(7) * t419 + t398 * t435 - t378;
t337 = -pkin(3) * t370 - pkin(8) * t393 + t396 * t380 + t359;
t318 = -pkin(4) * t342 - pkin(9) * t372 + t363 * t374 + t337;
t426 = m(6) * t318 - t323 * mrSges(6,1) + t324 * mrSges(6,2) - t352 * t346 + t353 * t347;
t288 = t416 * t292 + t411 * t293;
t332 = Ifges(6,4) * t353 + Ifges(6,2) * t352 + Ifges(6,6) * t401;
t333 = Ifges(6,1) * t353 + Ifges(6,4) * t352 + Ifges(6,5) * t401;
t425 = -mrSges(6,1) * t310 + mrSges(6,2) * t311 - Ifges(6,5) * t324 - Ifges(6,6) * t323 - Ifges(6,3) * t386 - t353 * t332 + t352 * t333;
t424 = m(5) * t337 - t342 * mrSges(5,1) + t343 * mrSges(5,2) - t373 * t361 + t374 * t362 + t426;
t349 = Ifges(5,4) * t374 + Ifges(5,2) * t373 + Ifges(5,6) * t404;
t350 = Ifges(5,1) * t374 + Ifges(5,4) * t373 + Ifges(5,5) * t404;
t423 = -mrSges(5,1) * t315 + mrSges(5,2) * t316 - Ifges(5,5) * t343 - Ifges(5,6) * t342 - Ifges(5,3) * t392 - pkin(4) * t301 - t374 * t349 + t373 * t350 + t425;
t422 = -m(4) * t359 + t370 * mrSges(4,1) - t371 * mrSges(4,2) + t395 * t376 - t396 * t377 - t424;
t365 = Ifges(4,4) * t396 + Ifges(4,2) * t395 + Ifges(4,6) * t405;
t366 = Ifges(4,1) * t396 + Ifges(4,4) * t395 + Ifges(4,5) * t405;
t421 = mrSges(4,1) * t344 - mrSges(4,2) * t345 + Ifges(4,5) * t371 + Ifges(4,6) * t370 + Ifges(4,3) * t394 + pkin(3) * t294 + t396 * t365 - t395 * t366 - t423;
t403 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t434;
t402 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t435;
t397 = (-mrSges(3,1) * t417 + mrSges(3,2) * t412) * qJD(1);
t389 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t412 + Ifges(3,4) * t417) * qJD(1);
t388 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t412 + Ifges(3,2) * t417) * qJD(1);
t364 = Ifges(4,5) * t396 + Ifges(4,6) * t395 + Ifges(4,3) * t405;
t348 = Ifges(5,5) * t374 + Ifges(5,6) * t373 + Ifges(5,3) * t404;
t331 = Ifges(6,5) * t353 + Ifges(6,6) * t352 + Ifges(6,3) * t401;
t303 = mrSges(6,2) * t318 - mrSges(6,3) * t310 + Ifges(6,1) * t324 + Ifges(6,4) * t323 + Ifges(6,5) * t386 + t331 * t352 - t332 * t401;
t302 = -mrSges(6,1) * t318 + mrSges(6,3) * t311 + Ifges(6,4) * t324 + Ifges(6,2) * t323 + Ifges(6,6) * t386 - t331 * t353 + t333 * t401;
t290 = mrSges(5,2) * t337 - mrSges(5,3) * t315 + Ifges(5,1) * t343 + Ifges(5,4) * t342 + Ifges(5,5) * t392 - pkin(9) * t301 - t302 * t409 + t303 * t414 + t348 * t373 - t349 * t404;
t289 = -mrSges(5,1) * t337 + mrSges(5,3) * t316 + Ifges(5,4) * t343 + Ifges(5,2) * t342 + Ifges(5,6) * t392 - pkin(4) * t426 + pkin(9) * t428 + t414 * t302 + t409 * t303 - t374 * t348 + t404 * t350;
t287 = mrSges(4,2) * t359 - mrSges(4,3) * t344 + Ifges(4,1) * t371 + Ifges(4,4) * t370 + Ifges(4,5) * t394 - pkin(8) * t294 - t289 * t410 + t290 * t415 + t364 * t395 - t365 * t405;
t286 = -mrSges(4,1) * t359 + mrSges(4,3) * t345 + Ifges(4,4) * t371 + Ifges(4,2) * t370 + Ifges(4,6) * t394 - pkin(3) * t424 + pkin(8) * t429 + t415 * t289 + t410 * t290 - t396 * t364 + t405 * t366;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t431 - mrSges(2,2) * t427 + t412 * (mrSges(3,2) * t390 - mrSges(3,3) * t378 + Ifges(3,1) * t399 + Ifges(3,4) * t400 + Ifges(3,5) * qJDD(2) - pkin(7) * t288 - qJD(2) * t388 - t411 * t286 + t416 * t287) + t417 * (-mrSges(3,1) * t390 + mrSges(3,3) * t379 + Ifges(3,4) * t399 + Ifges(3,2) * t400 + Ifges(3,6) * qJDD(2) - pkin(2) * t288 + qJD(2) * t389 - t421) + pkin(1) * (-m(3) * t390 + t400 * mrSges(3,1) - t399 * mrSges(3,2) + (-t402 * t412 + t403 * t417) * qJD(1) - t288) + pkin(6) * (t417 * (m(3) * t379 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t400 - qJD(2) * t402 + t397 * t434 + t430) - t412 * (m(3) * t378 + qJDD(2) * mrSges(3,1) - t399 * mrSges(3,3) + qJD(2) * t403 - t397 * t435 + t422)); Ifges(3,5) * t399 + Ifges(3,6) * t400 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t378 - mrSges(3,2) * t379 + t411 * t287 + t416 * t286 + pkin(2) * t422 + pkin(7) * t430 + (t412 * t388 - t417 * t389) * qJD(1); t421; -t423; -t425;];
tauJ = t1;
