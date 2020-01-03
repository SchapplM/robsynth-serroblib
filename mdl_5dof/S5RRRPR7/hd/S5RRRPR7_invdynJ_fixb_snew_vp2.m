% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:33
% EndTime: 2019-12-31 21:16:36
% DurationCPUTime: 2.81s
% Computational Cost: add. (29954->275), mult. (62360->355), div. (0->0), fcn. (43408->10), ass. (0->109)
t437 = cos(qJ(3));
t421 = qJD(1) ^ 2;
t436 = pkin(2) * t421;
t417 = sin(qJ(1));
t420 = cos(qJ(1));
t427 = -g(1) * t420 - g(2) * t417;
t396 = -pkin(1) * t421 + qJDD(1) * pkin(6) + t427;
t416 = sin(qJ(2));
t435 = t396 * t416;
t419 = cos(qJ(2));
t432 = qJD(1) * qJD(2);
t400 = qJDD(1) * t416 + t419 * t432;
t362 = qJDD(2) * pkin(2) - pkin(7) * t400 - t435 + (pkin(7) * t432 + t416 * t436 - g(3)) * t419;
t385 = -g(3) * t416 + t419 * t396;
t401 = qJDD(1) * t419 - t416 * t432;
t434 = qJD(1) * t416;
t404 = qJD(2) * pkin(2) - pkin(7) * t434;
t411 = t419 ^ 2;
t363 = pkin(7) * t401 - qJD(2) * t404 - t411 * t436 + t385;
t415 = sin(qJ(3));
t345 = t415 * t362 + t437 * t363;
t394 = (t415 * t419 + t437 * t416) * qJD(1);
t368 = qJD(3) * t394 + t400 * t415 - t437 * t401;
t433 = qJD(1) * t419;
t393 = t415 * t434 - t437 * t433;
t378 = mrSges(4,1) * t393 + mrSges(4,2) * t394;
t410 = qJD(2) + qJD(3);
t387 = mrSges(4,1) * t410 - mrSges(4,3) * t394;
t409 = qJDD(2) + qJDD(3);
t369 = -t393 * qJD(3) + t437 * t400 + t415 * t401;
t431 = g(1) * t417 - t420 * g(2);
t426 = -qJDD(1) * pkin(1) - t431;
t370 = -pkin(2) * t401 + t404 * t434 + (-pkin(7) * t411 - pkin(6)) * t421 + t426;
t334 = (t393 * t410 - t369) * qJ(4) + (t394 * t410 + t368) * pkin(3) + t370;
t377 = pkin(3) * t393 - qJ(4) * t394;
t408 = t410 ^ 2;
t337 = -pkin(3) * t408 + qJ(4) * t409 - t377 * t393 + t345;
t412 = sin(pkin(9));
t413 = cos(pkin(9));
t383 = t394 * t413 + t410 * t412;
t326 = -0.2e1 * qJD(4) * t383 + t413 * t334 - t337 * t412;
t356 = t369 * t413 + t409 * t412;
t382 = -t394 * t412 + t410 * t413;
t324 = (t382 * t393 - t356) * pkin(8) + (t382 * t383 + t368) * pkin(4) + t326;
t327 = 0.2e1 * qJD(4) * t382 + t412 * t334 + t413 * t337;
t355 = -t369 * t412 + t409 * t413;
t373 = pkin(4) * t393 - pkin(8) * t383;
t381 = t382 ^ 2;
t325 = -pkin(4) * t381 + pkin(8) * t355 - t373 * t393 + t327;
t414 = sin(qJ(5));
t418 = cos(qJ(5));
t322 = t324 * t418 - t325 * t414;
t353 = t382 * t418 - t383 * t414;
t333 = qJD(5) * t353 + t355 * t414 + t356 * t418;
t354 = t382 * t414 + t383 * t418;
t342 = -mrSges(6,1) * t353 + mrSges(6,2) * t354;
t389 = qJD(5) + t393;
t346 = -mrSges(6,2) * t389 + mrSges(6,3) * t353;
t367 = qJDD(5) + t368;
t318 = m(6) * t322 + mrSges(6,1) * t367 - mrSges(6,3) * t333 - t342 * t354 + t346 * t389;
t323 = t324 * t414 + t325 * t418;
t332 = -qJD(5) * t354 + t355 * t418 - t356 * t414;
t347 = mrSges(6,1) * t389 - mrSges(6,3) * t354;
t319 = m(6) * t323 - mrSges(6,2) * t367 + mrSges(6,3) * t332 + t342 * t353 - t347 * t389;
t310 = t418 * t318 + t414 * t319;
t358 = -mrSges(5,1) * t382 + mrSges(5,2) * t383;
t371 = -mrSges(5,2) * t393 + mrSges(5,3) * t382;
t308 = m(5) * t326 + mrSges(5,1) * t368 - mrSges(5,3) * t356 - t358 * t383 + t371 * t393 + t310;
t372 = mrSges(5,1) * t393 - mrSges(5,3) * t383;
t428 = -t318 * t414 + t418 * t319;
t309 = m(5) * t327 - mrSges(5,2) * t368 + mrSges(5,3) * t355 + t358 * t382 - t372 * t393 + t428;
t429 = -t308 * t412 + t413 * t309;
t302 = m(4) * t345 - mrSges(4,2) * t409 - mrSges(4,3) * t368 - t378 * t393 - t387 * t410 + t429;
t344 = t437 * t362 - t415 * t363;
t336 = -t409 * pkin(3) - t408 * qJ(4) + t394 * t377 + qJDD(4) - t344;
t328 = -t355 * pkin(4) - t381 * pkin(8) + t383 * t373 + t336;
t425 = m(6) * t328 - t332 * mrSges(6,1) + mrSges(6,2) * t333 - t353 * t346 + t347 * t354;
t321 = m(5) * t336 - t355 * mrSges(5,1) + mrSges(5,2) * t356 - t382 * t371 + t372 * t383 + t425;
t386 = -mrSges(4,2) * t410 - mrSges(4,3) * t393;
t314 = m(4) * t344 + mrSges(4,1) * t409 - mrSges(4,3) * t369 - t378 * t394 + t386 * t410 - t321;
t295 = t415 * t302 + t437 * t314;
t304 = t413 * t308 + t412 * t309;
t430 = t437 * t302 - t314 * t415;
t338 = Ifges(6,5) * t354 + Ifges(6,6) * t353 + Ifges(6,3) * t389;
t340 = Ifges(6,1) * t354 + Ifges(6,4) * t353 + Ifges(6,5) * t389;
t311 = -mrSges(6,1) * t328 + mrSges(6,3) * t323 + Ifges(6,4) * t333 + Ifges(6,2) * t332 + Ifges(6,6) * t367 - t338 * t354 + t340 * t389;
t339 = Ifges(6,4) * t354 + Ifges(6,2) * t353 + Ifges(6,6) * t389;
t312 = mrSges(6,2) * t328 - mrSges(6,3) * t322 + Ifges(6,1) * t333 + Ifges(6,4) * t332 + Ifges(6,5) * t367 + t338 * t353 - t339 * t389;
t348 = Ifges(5,5) * t383 + Ifges(5,6) * t382 + Ifges(5,3) * t393;
t350 = Ifges(5,1) * t383 + Ifges(5,4) * t382 + Ifges(5,5) * t393;
t297 = -mrSges(5,1) * t336 + mrSges(5,3) * t327 + Ifges(5,4) * t356 + Ifges(5,2) * t355 + Ifges(5,6) * t368 - pkin(4) * t425 + pkin(8) * t428 + t418 * t311 + t414 * t312 - t383 * t348 + t393 * t350;
t349 = Ifges(5,4) * t383 + Ifges(5,2) * t382 + Ifges(5,6) * t393;
t299 = mrSges(5,2) * t336 - mrSges(5,3) * t326 + Ifges(5,1) * t356 + Ifges(5,4) * t355 + Ifges(5,5) * t368 - pkin(8) * t310 - t311 * t414 + t312 * t418 + t348 * t382 - t349 * t393;
t375 = Ifges(4,4) * t394 - Ifges(4,2) * t393 + Ifges(4,6) * t410;
t376 = Ifges(4,1) * t394 - Ifges(4,4) * t393 + Ifges(4,5) * t410;
t424 = mrSges(4,1) * t344 - mrSges(4,2) * t345 + Ifges(4,5) * t369 - Ifges(4,6) * t368 + Ifges(4,3) * t409 - pkin(3) * t321 + qJ(4) * t429 + t413 * t297 + t412 * t299 + t394 * t375 + t376 * t393;
t423 = m(4) * t370 + mrSges(4,1) * t368 + mrSges(4,2) * t369 + t386 * t393 + t387 * t394 + t304;
t422 = mrSges(6,1) * t322 - mrSges(6,2) * t323 + Ifges(6,5) * t333 + Ifges(6,6) * t332 + Ifges(6,3) * t367 + t354 * t339 - t353 * t340;
t403 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t433;
t402 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t434;
t399 = (-t419 * mrSges(3,1) + t416 * mrSges(3,2)) * qJD(1);
t395 = -pkin(6) * t421 + t426;
t392 = Ifges(3,5) * qJD(2) + (t416 * Ifges(3,1) + t419 * Ifges(3,4)) * qJD(1);
t391 = Ifges(3,6) * qJD(2) + (t416 * Ifges(3,4) + t419 * Ifges(3,2)) * qJD(1);
t384 = -g(3) * t419 - t435;
t374 = Ifges(4,5) * t394 - Ifges(4,6) * t393 + Ifges(4,3) * t410;
t294 = -t422 + Ifges(4,6) * t409 + t410 * t376 - t394 * t374 + t382 * t350 - t383 * t349 + Ifges(4,4) * t369 - mrSges(4,1) * t370 - Ifges(5,6) * t355 - Ifges(5,5) * t356 + mrSges(4,3) * t345 - mrSges(5,1) * t326 + mrSges(5,2) * t327 - pkin(4) * t310 - pkin(3) * t304 + (-Ifges(4,2) - Ifges(5,3)) * t368;
t293 = mrSges(4,2) * t370 - mrSges(4,3) * t344 + Ifges(4,1) * t369 - Ifges(4,4) * t368 + Ifges(4,5) * t409 - qJ(4) * t304 - t297 * t412 + t299 * t413 - t374 * t393 - t375 * t410;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t431 - mrSges(2,2) * t427 + t416 * (mrSges(3,2) * t395 - mrSges(3,3) * t384 + Ifges(3,1) * t400 + Ifges(3,4) * t401 + Ifges(3,5) * qJDD(2) - pkin(7) * t295 - qJD(2) * t391 + t437 * t293 - t415 * t294) + t419 * (-mrSges(3,1) * t395 + mrSges(3,3) * t385 + Ifges(3,4) * t400 + Ifges(3,2) * t401 + Ifges(3,6) * qJDD(2) - pkin(2) * t423 + pkin(7) * t430 + qJD(2) * t392 + t415 * t293 + t437 * t294) + pkin(1) * (-m(3) * t395 + mrSges(3,1) * t401 - mrSges(3,2) * t400 + (-t402 * t416 + t403 * t419) * qJD(1) - t423) + pkin(6) * (t419 * (m(3) * t385 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t401 - qJD(2) * t402 + t399 * t433 + t430) - t416 * (m(3) * t384 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t400 + qJD(2) * t403 - t399 * t434 + t295)); Ifges(3,3) * qJDD(2) + t424 + mrSges(3,1) * t384 - mrSges(3,2) * t385 + Ifges(3,5) * t400 + Ifges(3,6) * t401 + pkin(2) * t295 + (t416 * t391 - t419 * t392) * qJD(1); t424; t321; t422;];
tauJ = t1;
