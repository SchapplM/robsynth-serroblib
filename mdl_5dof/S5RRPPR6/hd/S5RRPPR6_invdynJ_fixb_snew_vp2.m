% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:32:01
% EndTime: 2019-12-31 19:32:04
% DurationCPUTime: 2.68s
% Computational Cost: add. (22245->274), mult. (52389->355), div. (0->0), fcn. (35792->10), ass. (0->107)
t436 = -2 * qJD(3);
t413 = sin(qJ(2));
t416 = cos(qJ(2));
t429 = qJD(1) * qJD(2);
t401 = qJDD(1) * t413 + t416 * t429;
t419 = qJD(1) ^ 2;
t414 = sin(qJ(1));
t417 = cos(qJ(1));
t424 = -g(1) * t417 - g(2) * t414;
t398 = -pkin(1) * t419 + qJDD(1) * pkin(6) + t424;
t433 = t398 * t413;
t435 = pkin(2) * t419;
t359 = qJDD(2) * pkin(2) - qJ(3) * t401 - t433 + (qJ(3) * t429 + t413 * t435 - g(3)) * t416;
t384 = -g(3) * t413 + t416 * t398;
t402 = qJDD(1) * t416 - t413 * t429;
t432 = qJD(1) * t413;
t403 = qJD(2) * pkin(2) - qJ(3) * t432;
t408 = t416 ^ 2;
t361 = qJ(3) * t402 - qJD(2) * t403 - t408 * t435 + t384;
t410 = sin(pkin(8));
t434 = cos(pkin(8));
t392 = (t410 * t416 + t434 * t413) * qJD(1);
t345 = t434 * t359 - t410 * t361 + t392 * t436;
t431 = qJD(1) * t416;
t391 = t410 * t432 - t434 * t431;
t346 = t410 * t359 + t434 * t361 + t391 * t436;
t374 = mrSges(4,1) * t391 + mrSges(4,2) * t392;
t376 = t401 * t410 - t434 * t402;
t386 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t392;
t373 = pkin(3) * t391 - qJ(4) * t392;
t418 = qJD(2) ^ 2;
t334 = -pkin(3) * t418 + qJDD(2) * qJ(4) - t373 * t391 + t346;
t428 = g(1) * t414 - t417 * g(2);
t422 = -qJDD(1) * pkin(1) - t428;
t363 = -pkin(2) * t402 + qJDD(3) + t403 * t432 + (-qJ(3) * t408 - pkin(6)) * t419 + t422;
t377 = t434 * t401 + t410 * t402;
t337 = (qJD(2) * t391 - t377) * qJ(4) + (qJD(2) * t392 + t376) * pkin(3) + t363;
t409 = sin(pkin(9));
t411 = cos(pkin(9));
t382 = qJD(2) * t409 + t392 * t411;
t329 = -0.2e1 * qJD(4) * t382 - t334 * t409 + t411 * t337;
t369 = qJDD(2) * t409 + t377 * t411;
t381 = qJD(2) * t411 - t392 * t409;
t327 = (t381 * t391 - t369) * pkin(7) + (t381 * t382 + t376) * pkin(4) + t329;
t330 = 0.2e1 * qJD(4) * t381 + t411 * t334 + t409 * t337;
t366 = pkin(4) * t391 - pkin(7) * t382;
t368 = qJDD(2) * t411 - t377 * t409;
t380 = t381 ^ 2;
t328 = -pkin(4) * t380 + pkin(7) * t368 - t366 * t391 + t330;
t412 = sin(qJ(5));
t415 = cos(qJ(5));
t325 = t327 * t415 - t328 * t412;
t355 = t381 * t415 - t382 * t412;
t340 = qJD(5) * t355 + t368 * t412 + t369 * t415;
t356 = t381 * t412 + t382 * t415;
t347 = -mrSges(6,1) * t355 + mrSges(6,2) * t356;
t390 = qJD(5) + t391;
t348 = -mrSges(6,2) * t390 + mrSges(6,3) * t355;
t375 = qJDD(5) + t376;
t322 = m(6) * t325 + mrSges(6,1) * t375 - mrSges(6,3) * t340 - t347 * t356 + t348 * t390;
t326 = t327 * t412 + t328 * t415;
t339 = -qJD(5) * t356 + t368 * t415 - t369 * t412;
t349 = mrSges(6,1) * t390 - mrSges(6,3) * t356;
t323 = m(6) * t326 - mrSges(6,2) * t375 + mrSges(6,3) * t339 + t347 * t355 - t349 * t390;
t314 = t415 * t322 + t412 * t323;
t360 = -mrSges(5,1) * t381 + mrSges(5,2) * t382;
t364 = -mrSges(5,2) * t391 + mrSges(5,3) * t381;
t312 = m(5) * t329 + mrSges(5,1) * t376 - mrSges(5,3) * t369 - t360 * t382 + t364 * t391 + t314;
t365 = mrSges(5,1) * t391 - mrSges(5,3) * t382;
t425 = -t322 * t412 + t415 * t323;
t313 = m(5) * t330 - mrSges(5,2) * t376 + mrSges(5,3) * t368 + t360 * t381 - t365 * t391 + t425;
t426 = -t312 * t409 + t411 * t313;
t306 = m(4) * t346 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t376 - qJD(2) * t386 - t374 * t391 + t426;
t333 = -qJDD(2) * pkin(3) - t418 * qJ(4) + t392 * t373 + qJDD(4) - t345;
t331 = -t368 * pkin(4) - t380 * pkin(7) + t382 * t366 + t333;
t421 = m(6) * t331 - t339 * mrSges(6,1) + mrSges(6,2) * t340 - t355 * t348 + t349 * t356;
t324 = m(5) * t333 - t368 * mrSges(5,1) + mrSges(5,2) * t369 - t381 * t364 + t365 * t382 + t421;
t385 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t391;
t318 = m(4) * t345 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t377 + qJD(2) * t385 - t374 * t392 - t324;
t301 = t410 * t306 + t434 * t318;
t308 = t411 * t312 + t409 * t313;
t427 = t434 * t306 - t318 * t410;
t307 = m(4) * t363 + mrSges(4,1) * t376 + mrSges(4,2) * t377 + t385 * t391 + t386 * t392 + t308;
t342 = Ifges(6,4) * t356 + Ifges(6,2) * t355 + Ifges(6,6) * t390;
t343 = Ifges(6,1) * t356 + Ifges(6,4) * t355 + Ifges(6,5) * t390;
t420 = mrSges(6,1) * t325 - mrSges(6,2) * t326 + Ifges(6,5) * t340 + Ifges(6,6) * t339 + Ifges(6,3) * t375 + t356 * t342 - t355 * t343;
t405 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t431;
t404 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t432;
t400 = (-t416 * mrSges(3,1) + t413 * mrSges(3,2)) * qJD(1);
t397 = -pkin(6) * t419 + t422;
t395 = Ifges(3,5) * qJD(2) + (t413 * Ifges(3,1) + t416 * Ifges(3,4)) * qJD(1);
t394 = Ifges(3,6) * qJD(2) + (t413 * Ifges(3,4) + t416 * Ifges(3,2)) * qJD(1);
t383 = -g(3) * t416 - t433;
t372 = Ifges(4,1) * t392 - Ifges(4,4) * t391 + Ifges(4,5) * qJD(2);
t371 = Ifges(4,4) * t392 - Ifges(4,2) * t391 + Ifges(4,6) * qJD(2);
t370 = Ifges(4,5) * t392 - Ifges(4,6) * t391 + Ifges(4,3) * qJD(2);
t352 = Ifges(5,1) * t382 + Ifges(5,4) * t381 + Ifges(5,5) * t391;
t351 = Ifges(5,4) * t382 + Ifges(5,2) * t381 + Ifges(5,6) * t391;
t350 = Ifges(5,5) * t382 + Ifges(5,6) * t381 + Ifges(5,3) * t391;
t341 = Ifges(6,5) * t356 + Ifges(6,6) * t355 + Ifges(6,3) * t390;
t316 = mrSges(6,2) * t331 - mrSges(6,3) * t325 + Ifges(6,1) * t340 + Ifges(6,4) * t339 + Ifges(6,5) * t375 + t341 * t355 - t342 * t390;
t315 = -mrSges(6,1) * t331 + mrSges(6,3) * t326 + Ifges(6,4) * t340 + Ifges(6,2) * t339 + Ifges(6,6) * t375 - t341 * t356 + t343 * t390;
t303 = mrSges(5,2) * t333 - mrSges(5,3) * t329 + Ifges(5,1) * t369 + Ifges(5,4) * t368 + Ifges(5,5) * t376 - pkin(7) * t314 - t315 * t412 + t316 * t415 + t350 * t381 - t351 * t391;
t302 = -mrSges(5,1) * t333 + mrSges(5,3) * t330 + Ifges(5,4) * t369 + Ifges(5,2) * t368 + Ifges(5,6) * t376 - pkin(4) * t421 + pkin(7) * t425 + t415 * t315 + t412 * t316 - t382 * t350 + t391 * t352;
t300 = Ifges(4,6) * qJDD(2) - t392 * t370 + t381 * t352 - t382 * t351 - Ifges(5,5) * t369 + qJD(2) * t372 + Ifges(4,4) * t377 - mrSges(4,1) * t363 - Ifges(5,6) * t368 + mrSges(4,3) * t346 + mrSges(5,2) * t330 - mrSges(5,1) * t329 - pkin(4) * t314 - pkin(3) * t308 + (-Ifges(4,2) - Ifges(5,3)) * t376 - t420;
t299 = mrSges(4,2) * t363 - mrSges(4,3) * t345 + Ifges(4,1) * t377 - Ifges(4,4) * t376 + Ifges(4,5) * qJDD(2) - qJ(4) * t308 - qJD(2) * t371 - t302 * t409 + t303 * t411 - t370 * t391;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t428 - mrSges(2,2) * t424 + t413 * (mrSges(3,2) * t397 - mrSges(3,3) * t383 + Ifges(3,1) * t401 + Ifges(3,4) * t402 + Ifges(3,5) * qJDD(2) - qJ(3) * t301 - qJD(2) * t394 + t434 * t299 - t410 * t300) + t416 * (-mrSges(3,1) * t397 + mrSges(3,3) * t384 + Ifges(3,4) * t401 + Ifges(3,2) * t402 + Ifges(3,6) * qJDD(2) - pkin(2) * t307 + qJ(3) * t427 + qJD(2) * t395 + t410 * t299 + t434 * t300) + pkin(1) * (-m(3) * t397 + mrSges(3,1) * t402 - t401 * mrSges(3,2) + (-t404 * t413 + t405 * t416) * qJD(1) - t307) + pkin(6) * (t416 * (m(3) * t384 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t402 - qJD(2) * t404 + t400 * t431 + t427) - t413 * (m(3) * t383 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t401 + qJD(2) * t405 - t400 * t432 + t301)); Ifges(3,5) * t401 + Ifges(3,6) * t402 + mrSges(3,1) * t383 - mrSges(3,2) * t384 + Ifges(4,5) * t377 - Ifges(4,6) * t376 + t392 * t371 + t391 * t372 + mrSges(4,1) * t345 - mrSges(4,2) * t346 + t409 * t303 + t411 * t302 - pkin(3) * t324 + qJ(4) * t426 + pkin(2) * t301 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t413 * t394 - t416 * t395) * qJD(1); t307; t324; t420;];
tauJ = t1;
