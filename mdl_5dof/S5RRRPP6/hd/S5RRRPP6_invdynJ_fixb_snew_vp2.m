% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:18
% EndTime: 2019-12-31 21:00:21
% DurationCPUTime: 1.98s
% Computational Cost: add. (13289->251), mult. (27265->309), div. (0->0), fcn. (17961->8), ass. (0->101)
t430 = Ifges(5,1) + Ifges(6,1);
t421 = Ifges(5,4) - Ifges(6,5);
t420 = Ifges(5,5) + Ifges(6,4);
t429 = -Ifges(5,2) - Ifges(6,3);
t419 = Ifges(5,6) - Ifges(6,6);
t428 = -Ifges(6,2) - Ifges(5,3);
t398 = qJD(1) ^ 2;
t393 = sin(qJ(1));
t396 = cos(qJ(1));
t406 = g(1) * t393 - t396 * g(2);
t372 = -qJDD(1) * pkin(1) - pkin(6) * t398 - t406;
t392 = sin(qJ(2));
t395 = cos(qJ(2));
t411 = qJD(1) * qJD(2);
t407 = t395 * t411;
t382 = qJDD(1) * t392 + t407;
t408 = t392 * t411;
t383 = qJDD(1) * t395 - t408;
t337 = (-t382 - t407) * pkin(7) + (-t383 + t408) * pkin(2) + t372;
t402 = -g(1) * t396 - g(2) * t393;
t373 = -pkin(1) * t398 + qJDD(1) * pkin(6) + t402;
t364 = -g(3) * t392 + t395 * t373;
t381 = (-t395 * pkin(2) - t392 * pkin(7)) * qJD(1);
t397 = qJD(2) ^ 2;
t413 = qJD(1) * t395;
t342 = -pkin(2) * t397 + qJDD(2) * pkin(7) + t381 * t413 + t364;
t391 = sin(qJ(3));
t394 = cos(qJ(3));
t316 = t394 * t337 - t342 * t391;
t412 = t392 * qJD(1);
t378 = qJD(2) * t394 - t391 * t412;
t356 = qJD(3) * t378 + qJDD(2) * t391 + t382 * t394;
t377 = qJDD(3) - t383;
t379 = qJD(2) * t391 + t394 * t412;
t387 = qJD(3) - t413;
t312 = (t378 * t387 - t356) * qJ(4) + (t378 * t379 + t377) * pkin(3) + t316;
t317 = t391 * t337 + t394 * t342;
t355 = -qJD(3) * t379 + qJDD(2) * t394 - t382 * t391;
t361 = pkin(3) * t387 - qJ(4) * t379;
t376 = t378 ^ 2;
t314 = -pkin(3) * t376 + qJ(4) * t355 - t361 * t387 + t317;
t390 = sin(pkin(8));
t418 = cos(pkin(8));
t357 = -t418 * t378 + t379 * t390;
t423 = -2 * qJD(4);
t308 = t390 * t312 + t418 * t314 + t357 * t423;
t327 = -t418 * t355 + t356 * t390;
t358 = t390 * t378 + t418 * t379;
t344 = mrSges(5,1) * t387 - mrSges(5,3) * t358;
t332 = pkin(4) * t357 - qJ(5) * t358;
t386 = t387 ^ 2;
t305 = -pkin(4) * t386 + qJ(5) * t377 + 0.2e1 * qJD(5) * t387 - t332 * t357 + t308;
t345 = -mrSges(6,1) * t387 + mrSges(6,2) * t358;
t409 = m(6) * t305 + t377 * mrSges(6,3) + t387 * t345;
t333 = mrSges(6,1) * t357 - mrSges(6,3) * t358;
t414 = -mrSges(5,1) * t357 - mrSges(5,2) * t358 - t333;
t422 = -mrSges(5,3) - mrSges(6,2);
t297 = m(5) * t308 - mrSges(5,2) * t377 + t422 * t327 - t344 * t387 + t414 * t357 + t409;
t401 = t418 * t312 - t390 * t314;
t307 = t358 * t423 + t401;
t328 = t390 * t355 + t418 * t356;
t343 = -mrSges(5,2) * t387 - mrSges(5,3) * t357;
t306 = -t377 * pkin(4) - t386 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t332) * t358 - t401;
t346 = -mrSges(6,2) * t357 + mrSges(6,3) * t387;
t403 = -m(6) * t306 + t377 * mrSges(6,1) + t387 * t346;
t299 = m(5) * t307 + mrSges(5,1) * t377 + t422 * t328 + t343 * t387 + t414 * t358 + t403;
t294 = t390 * t297 + t418 * t299;
t302 = mrSges(6,2) * t328 + t333 * t358 - t403;
t349 = Ifges(4,4) * t379 + Ifges(4,2) * t378 + Ifges(4,6) * t387;
t350 = Ifges(4,1) * t379 + Ifges(4,4) * t378 + Ifges(4,5) * t387;
t415 = -t421 * t357 + t430 * t358 + t420 * t387;
t416 = t429 * t357 + t421 * t358 + t419 * t387;
t427 = -t419 * t327 + t420 * t328 + mrSges(4,1) * t316 + mrSges(5,1) * t307 - mrSges(6,1) * t306 - mrSges(4,2) * t317 - mrSges(5,2) * t308 + mrSges(6,3) * t305 + Ifges(4,5) * t356 + Ifges(4,6) * t355 + pkin(3) * t294 - pkin(4) * t302 + qJ(5) * (-mrSges(6,2) * t327 - t333 * t357 + t409) + t379 * t349 - t378 * t350 + t416 * t358 + t415 * t357;
t426 = (Ifges(4,3) - t428) * t377;
t417 = t419 * t357 - t420 * t358 + t428 * t387;
t363 = -t395 * g(3) - t392 * t373;
t359 = -mrSges(4,1) * t378 + mrSges(4,2) * t379;
t360 = -mrSges(4,2) * t387 + mrSges(4,3) * t378;
t290 = m(4) * t316 + mrSges(4,1) * t377 - mrSges(4,3) * t356 - t359 * t379 + t360 * t387 + t294;
t362 = mrSges(4,1) * t387 - mrSges(4,3) * t379;
t404 = t418 * t297 - t299 * t390;
t291 = m(4) * t317 - mrSges(4,2) * t377 + mrSges(4,3) * t355 + t359 * t378 - t362 * t387 + t404;
t405 = -t290 * t391 + t394 * t291;
t341 = -qJDD(2) * pkin(2) - pkin(7) * t397 + t381 * t412 - t363;
t315 = -pkin(3) * t355 - qJ(4) * t376 + t379 * t361 + qJDD(4) + t341;
t310 = -0.2e1 * qJD(5) * t358 + (t357 * t387 - t328) * qJ(5) + (t358 * t387 + t327) * pkin(4) + t315;
t303 = m(6) * t310 + t327 * mrSges(6,1) - t328 * mrSges(6,3) - t358 * t345 + t357 * t346;
t288 = t290 * t394 + t291 * t391;
t300 = m(5) * t315 + t327 * mrSges(5,1) + t328 * mrSges(5,2) + t357 * t343 + t358 * t344 + t303;
t400 = -m(4) * t341 + t355 * mrSges(4,1) - t356 * mrSges(4,2) + t378 * t360 - t379 * t362 - t300;
t385 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t413;
t384 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t412;
t380 = (-t395 * mrSges(3,1) + t392 * mrSges(3,2)) * qJD(1);
t371 = Ifges(3,5) * qJD(2) + (t392 * Ifges(3,1) + t395 * Ifges(3,4)) * qJD(1);
t370 = Ifges(3,6) * qJD(2) + (t392 * Ifges(3,4) + t395 * Ifges(3,2)) * qJD(1);
t348 = Ifges(4,5) * t379 + Ifges(4,6) * t378 + Ifges(4,3) * t387;
t293 = mrSges(5,2) * t315 + mrSges(6,2) * t306 - mrSges(5,3) * t307 - mrSges(6,3) * t310 - qJ(5) * t303 - t421 * t327 + t430 * t328 + t417 * t357 + t420 * t377 - t416 * t387;
t292 = -mrSges(5,1) * t315 - mrSges(6,1) * t310 + mrSges(6,2) * t305 + mrSges(5,3) * t308 - pkin(4) * t303 + t429 * t327 + t421 * t328 + t417 * t358 + t419 * t377 + t415 * t387;
t287 = mrSges(4,2) * t341 - mrSges(4,3) * t316 + Ifges(4,1) * t356 + Ifges(4,4) * t355 + Ifges(4,5) * t377 - qJ(4) * t294 - t390 * t292 + t418 * t293 + t378 * t348 - t387 * t349;
t286 = -mrSges(4,1) * t341 + mrSges(4,3) * t317 + Ifges(4,4) * t356 + Ifges(4,2) * t355 + Ifges(4,6) * t377 - pkin(3) * t300 + qJ(4) * t404 + t418 * t292 + t390 * t293 - t379 * t348 + t387 * t350;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t406 - mrSges(2,2) * t402 + t392 * (mrSges(3,2) * t372 - mrSges(3,3) * t363 + Ifges(3,1) * t382 + Ifges(3,4) * t383 + Ifges(3,5) * qJDD(2) - pkin(7) * t288 - qJD(2) * t370 - t391 * t286 + t394 * t287) + t395 * (-mrSges(3,1) * t372 + mrSges(3,3) * t364 + Ifges(3,4) * t382 + Ifges(3,2) * t383 + Ifges(3,6) * qJDD(2) - pkin(2) * t288 + qJD(2) * t371 - t427) + pkin(1) * (-m(3) * t372 + mrSges(3,1) * t383 - mrSges(3,2) * t382 + (-t384 * t392 + t385 * t395) * qJD(1) - t288) + pkin(6) * (t395 * (m(3) * t364 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t383 - qJD(2) * t384 + t380 * t413 + t405) - t392 * (m(3) * t363 + qJDD(2) * mrSges(3,1) - t382 * mrSges(3,3) + qJD(2) * t385 - t380 * t412 + t400)) - t395 * t426; Ifges(3,5) * t382 + Ifges(3,6) * t383 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t363 - mrSges(3,2) * t364 + t391 * t287 + t394 * t286 + pkin(2) * t400 + pkin(7) * t405 + (t392 * t370 - t395 * t371) * qJD(1); t426 + t427; t300; t302;];
tauJ = t1;
