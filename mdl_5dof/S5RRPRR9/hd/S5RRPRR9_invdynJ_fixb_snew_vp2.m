% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR9
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:28
% EndTime: 2019-12-31 20:20:32
% DurationCPUTime: 2.70s
% Computational Cost: add. (24559->276), mult. (55817->355), div. (0->0), fcn. (38588->10), ass. (0->109)
t415 = sin(qJ(2));
t419 = cos(qJ(2));
t434 = qJD(1) * qJD(2);
t403 = qJDD(1) * t415 + t419 * t434;
t422 = qJD(1) ^ 2;
t416 = sin(qJ(1));
t420 = cos(qJ(1));
t428 = -g(1) * t420 - g(2) * t416;
t400 = -pkin(1) * t422 + qJDD(1) * pkin(6) + t428;
t438 = t400 * t415;
t439 = pkin(2) * t422;
t363 = qJDD(2) * pkin(2) - qJ(3) * t403 - t438 + (qJ(3) * t434 + t415 * t439 - g(3)) * t419;
t385 = -g(3) * t415 + t419 * t400;
t404 = qJDD(1) * t419 - t415 * t434;
t437 = qJD(1) * t415;
t405 = qJD(2) * pkin(2) - qJ(3) * t437;
t410 = t419 ^ 2;
t364 = qJ(3) * t404 - qJD(2) * t405 - t410 * t439 + t385;
t411 = sin(pkin(9));
t412 = cos(pkin(9));
t394 = (t411 * t419 + t412 * t415) * qJD(1);
t346 = -0.2e1 * qJD(3) * t394 + t363 * t412 - t411 * t364;
t436 = qJD(1) * t419;
t393 = -t411 * t437 + t412 * t436;
t347 = 0.2e1 * qJD(3) * t393 + t411 * t363 + t412 * t364;
t374 = -mrSges(4,1) * t393 + mrSges(4,2) * t394;
t379 = -t411 * t403 + t404 * t412;
t387 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t394;
t376 = -pkin(3) * t393 - pkin(7) * t394;
t421 = qJD(2) ^ 2;
t338 = -pkin(3) * t421 + qJDD(2) * pkin(7) + t376 * t393 + t347;
t433 = g(1) * t416 - t420 * g(2);
t427 = -qJDD(1) * pkin(1) - t433;
t367 = -pkin(2) * t404 + qJDD(3) + t405 * t437 + (-qJ(3) * t410 - pkin(6)) * t422 + t427;
t380 = t403 * t412 + t404 * t411;
t341 = (-qJD(2) * t393 - t380) * pkin(7) + (qJD(2) * t394 - t379) * pkin(3) + t367;
t414 = sin(qJ(4));
t418 = cos(qJ(4));
t328 = -t338 * t414 + t418 * t341;
t382 = qJD(2) * t418 - t394 * t414;
t357 = qJD(4) * t382 + qJDD(2) * t414 + t380 * t418;
t378 = qJDD(4) - t379;
t383 = qJD(2) * t414 + t394 * t418;
t392 = qJD(4) - t393;
t325 = (t382 * t392 - t357) * pkin(8) + (t382 * t383 + t378) * pkin(4) + t328;
t329 = t418 * t338 + t414 * t341;
t356 = -qJD(4) * t383 + qJDD(2) * t418 - t380 * t414;
t370 = pkin(4) * t392 - pkin(8) * t383;
t381 = t382 ^ 2;
t326 = -pkin(4) * t381 + pkin(8) * t356 - t370 * t392 + t329;
t413 = sin(qJ(5));
t417 = cos(qJ(5));
t323 = t325 * t417 - t326 * t413;
t361 = t382 * t417 - t383 * t413;
t334 = qJD(5) * t361 + t356 * t413 + t357 * t417;
t362 = t382 * t413 + t383 * t417;
t348 = -mrSges(6,1) * t361 + mrSges(6,2) * t362;
t388 = qJD(5) + t392;
t349 = -mrSges(6,2) * t388 + mrSges(6,3) * t361;
t377 = qJDD(5) + t378;
t320 = m(6) * t323 + mrSges(6,1) * t377 - mrSges(6,3) * t334 - t348 * t362 + t349 * t388;
t324 = t325 * t413 + t326 * t417;
t333 = -qJD(5) * t362 + t356 * t417 - t357 * t413;
t350 = mrSges(6,1) * t388 - mrSges(6,3) * t362;
t321 = m(6) * t324 - mrSges(6,2) * t377 + mrSges(6,3) * t333 + t348 * t361 - t350 * t388;
t312 = t320 * t417 + t321 * t413;
t365 = -mrSges(5,1) * t382 + mrSges(5,2) * t383;
t368 = -mrSges(5,2) * t392 + mrSges(5,3) * t382;
t310 = m(5) * t328 + mrSges(5,1) * t378 - mrSges(5,3) * t357 - t365 * t383 + t368 * t392 + t312;
t369 = mrSges(5,1) * t392 - mrSges(5,3) * t383;
t430 = -t320 * t413 + t321 * t417;
t311 = m(5) * t329 - mrSges(5,2) * t378 + mrSges(5,3) * t356 + t365 * t382 - t369 * t392 + t430;
t431 = -t310 * t414 + t311 * t418;
t304 = m(4) * t347 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t379 - qJD(2) * t387 + t374 * t393 + t431;
t386 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t393;
t337 = -qJDD(2) * pkin(3) - pkin(7) * t421 + t394 * t376 - t346;
t327 = -pkin(4) * t356 - pkin(8) * t381 + t370 * t383 + t337;
t426 = m(6) * t327 - t333 * mrSges(6,1) + mrSges(6,2) * t334 - t361 * t349 + t350 * t362;
t424 = -m(5) * t337 + t356 * mrSges(5,1) - mrSges(5,2) * t357 + t382 * t368 - t369 * t383 - t426;
t316 = m(4) * t346 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t380 + qJD(2) * t386 - t374 * t394 + t424;
t300 = t304 * t411 + t316 * t412;
t306 = t310 * t418 + t311 * t414;
t432 = t304 * t412 - t316 * t411;
t343 = Ifges(6,4) * t362 + Ifges(6,2) * t361 + Ifges(6,6) * t388;
t344 = Ifges(6,1) * t362 + Ifges(6,4) * t361 + Ifges(6,5) * t388;
t425 = -mrSges(6,1) * t323 + mrSges(6,2) * t324 - Ifges(6,5) * t334 - Ifges(6,6) * t333 - Ifges(6,3) * t377 - t362 * t343 + t344 * t361;
t305 = m(4) * t367 - mrSges(4,1) * t379 + mrSges(4,2) * t380 - t386 * t393 + t387 * t394 + t306;
t352 = Ifges(5,4) * t383 + Ifges(5,2) * t382 + Ifges(5,6) * t392;
t353 = Ifges(5,1) * t383 + Ifges(5,4) * t382 + Ifges(5,5) * t392;
t423 = mrSges(5,1) * t328 - mrSges(5,2) * t329 + Ifges(5,5) * t357 + Ifges(5,6) * t356 + Ifges(5,3) * t378 + pkin(4) * t312 + t352 * t383 - t353 * t382 - t425;
t407 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t436;
t406 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t437;
t402 = (-mrSges(3,1) * t419 + mrSges(3,2) * t415) * qJD(1);
t399 = -pkin(6) * t422 + t427;
t397 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t415 + Ifges(3,4) * t419) * qJD(1);
t396 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t415 + Ifges(3,2) * t419) * qJD(1);
t384 = -g(3) * t419 - t438;
t373 = Ifges(4,1) * t394 + Ifges(4,4) * t393 + Ifges(4,5) * qJD(2);
t372 = Ifges(4,4) * t394 + Ifges(4,2) * t393 + Ifges(4,6) * qJD(2);
t371 = Ifges(4,5) * t394 + Ifges(4,6) * t393 + Ifges(4,3) * qJD(2);
t351 = Ifges(5,5) * t383 + Ifges(5,6) * t382 + Ifges(5,3) * t392;
t342 = Ifges(6,5) * t362 + Ifges(6,6) * t361 + Ifges(6,3) * t388;
t314 = mrSges(6,2) * t327 - mrSges(6,3) * t323 + Ifges(6,1) * t334 + Ifges(6,4) * t333 + Ifges(6,5) * t377 + t342 * t361 - t343 * t388;
t313 = -mrSges(6,1) * t327 + mrSges(6,3) * t324 + Ifges(6,4) * t334 + Ifges(6,2) * t333 + Ifges(6,6) * t377 - t342 * t362 + t344 * t388;
t301 = mrSges(5,2) * t337 - mrSges(5,3) * t328 + Ifges(5,1) * t357 + Ifges(5,4) * t356 + Ifges(5,5) * t378 - pkin(8) * t312 - t313 * t413 + t314 * t417 + t351 * t382 - t352 * t392;
t299 = -mrSges(5,1) * t337 + mrSges(5,3) * t329 + Ifges(5,4) * t357 + Ifges(5,2) * t356 + Ifges(5,6) * t378 - pkin(4) * t426 + pkin(8) * t430 + t417 * t313 + t413 * t314 - t383 * t351 + t392 * t353;
t298 = -mrSges(4,1) * t367 + mrSges(4,3) * t347 + Ifges(4,4) * t380 + Ifges(4,2) * t379 + Ifges(4,6) * qJDD(2) - pkin(3) * t306 + qJD(2) * t373 - t371 * t394 - t423;
t297 = mrSges(4,2) * t367 - mrSges(4,3) * t346 + Ifges(4,1) * t380 + Ifges(4,4) * t379 + Ifges(4,5) * qJDD(2) - pkin(7) * t306 - qJD(2) * t372 - t299 * t414 + t301 * t418 + t371 * t393;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t433 - mrSges(2,2) * t428 + t415 * (mrSges(3,2) * t399 - mrSges(3,3) * t384 + Ifges(3,1) * t403 + Ifges(3,4) * t404 + Ifges(3,5) * qJDD(2) - qJ(3) * t300 - qJD(2) * t396 + t412 * t297 - t411 * t298) + t419 * (-mrSges(3,1) * t399 + mrSges(3,3) * t385 + Ifges(3,4) * t403 + Ifges(3,2) * t404 + Ifges(3,6) * qJDD(2) - pkin(2) * t305 + qJ(3) * t432 + qJD(2) * t397 + t411 * t297 + t412 * t298) + pkin(1) * (-m(3) * t399 + mrSges(3,1) * t404 - mrSges(3,2) * t403 + (-t406 * t415 + t407 * t419) * qJD(1) - t305) + pkin(6) * (t419 * (m(3) * t385 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t404 - qJD(2) * t406 + t402 * t436 + t432) - t415 * (m(3) * t384 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t403 + qJD(2) * t407 - t402 * t437 + t300)); Ifges(3,5) * t403 + Ifges(3,6) * t404 + mrSges(3,1) * t384 - mrSges(3,2) * t385 + Ifges(4,5) * t380 + Ifges(4,6) * t379 + t394 * t372 - t393 * t373 + mrSges(4,1) * t346 - mrSges(4,2) * t347 + t414 * t301 + t418 * t299 + pkin(3) * t424 + pkin(7) * t431 + pkin(2) * t300 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t415 * t396 - t419 * t397) * qJD(1); t305; t423; -t425;];
tauJ = t1;
