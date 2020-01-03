% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP9
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:04:08
% EndTime: 2019-12-31 22:04:12
% DurationCPUTime: 2.06s
% Computational Cost: add. (14592->251), mult. (28946->307), div. (0->0), fcn. (19277->8), ass. (0->101)
t436 = Ifges(5,1) + Ifges(6,1);
t427 = Ifges(5,4) - Ifges(6,5);
t434 = Ifges(6,4) + Ifges(5,5);
t435 = Ifges(5,2) + Ifges(6,3);
t433 = Ifges(5,6) - Ifges(6,6);
t432 = -Ifges(5,3) - Ifges(6,2);
t403 = sin(qJ(3));
t406 = cos(qJ(3));
t404 = sin(qJ(2));
t424 = qJD(1) * t404;
t388 = qJD(2) * t406 - t403 * t424;
t389 = qJD(2) * t403 + t406 * t424;
t402 = sin(qJ(4));
t429 = cos(qJ(4));
t364 = -t429 * t388 + t389 * t402;
t365 = t402 * t388 + t429 * t389;
t407 = cos(qJ(2));
t423 = qJD(1) * t407;
t398 = qJD(3) - t423;
t397 = qJD(4) + t398;
t431 = t435 * t364 - t427 * t365 - t433 * t397;
t430 = -t427 * t364 + t436 * t365 + t434 * t397;
t428 = -mrSges(5,3) - mrSges(6,2);
t410 = qJD(1) ^ 2;
t405 = sin(qJ(1));
t408 = cos(qJ(1));
t419 = g(1) * t405 - t408 * g(2);
t383 = -qJDD(1) * pkin(1) - pkin(6) * t410 - t419;
t422 = qJD(1) * qJD(2);
t420 = t407 * t422;
t392 = qJDD(1) * t404 + t420;
t399 = t404 * t422;
t393 = qJDD(1) * t407 - t399;
t345 = (-t392 - t420) * pkin(7) + (-t393 + t399) * pkin(2) + t383;
t415 = -g(1) * t408 - g(2) * t405;
t384 = -pkin(1) * t410 + qJDD(1) * pkin(6) + t415;
t370 = -g(3) * t404 + t407 * t384;
t391 = (-pkin(2) * t407 - pkin(7) * t404) * qJD(1);
t409 = qJD(2) ^ 2;
t350 = -pkin(2) * t409 + qJDD(2) * pkin(7) + t391 * t423 + t370;
t329 = t406 * t345 - t350 * t403;
t363 = qJD(3) * t388 + qJDD(2) * t403 + t392 * t406;
t387 = qJDD(3) - t393;
t314 = (t388 * t398 - t363) * pkin(8) + (t388 * t389 + t387) * pkin(3) + t329;
t330 = t403 * t345 + t406 * t350;
t362 = -qJD(3) * t389 + qJDD(2) * t406 - t392 * t403;
t371 = pkin(3) * t398 - pkin(8) * t389;
t386 = t388 ^ 2;
t316 = -pkin(3) * t386 + pkin(8) * t362 - t371 * t398 + t330;
t312 = t402 * t314 + t429 * t316;
t327 = qJD(4) * t365 - t429 * t362 + t363 * t402;
t353 = mrSges(5,1) * t397 - mrSges(5,3) * t365;
t385 = qJDD(4) + t387;
t340 = pkin(4) * t364 - qJ(5) * t365;
t396 = t397 ^ 2;
t306 = -pkin(4) * t396 + qJ(5) * t385 + 0.2e1 * qJD(5) * t397 - t340 * t364 + t312;
t354 = -mrSges(6,1) * t397 + mrSges(6,2) * t365;
t421 = m(6) * t306 + t385 * mrSges(6,3) + t397 * t354;
t341 = mrSges(6,1) * t364 - mrSges(6,3) * t365;
t425 = -mrSges(5,1) * t364 - mrSges(5,2) * t365 - t341;
t296 = m(5) * t312 - mrSges(5,2) * t385 + t428 * t327 - t353 * t397 + t425 * t364 + t421;
t311 = t429 * t314 - t402 * t316;
t328 = -t364 * qJD(4) + t402 * t362 + t429 * t363;
t352 = -mrSges(5,2) * t397 - mrSges(5,3) * t364;
t307 = -t385 * pkin(4) - t396 * qJ(5) + t365 * t340 + qJDD(5) - t311;
t351 = -mrSges(6,2) * t364 + mrSges(6,3) * t397;
t416 = -m(6) * t307 + t385 * mrSges(6,1) + t397 * t351;
t298 = m(5) * t311 + mrSges(5,1) * t385 + t428 * t328 + t352 * t397 + t425 * t365 + t416;
t293 = t402 * t296 + t429 * t298;
t426 = t433 * t364 - t434 * t365 + t432 * t397;
t369 = -t407 * g(3) - t404 * t384;
t366 = -mrSges(4,1) * t388 + mrSges(4,2) * t389;
t367 = -mrSges(4,2) * t398 + mrSges(4,3) * t388;
t291 = m(4) * t329 + mrSges(4,1) * t387 - mrSges(4,3) * t363 - t366 * t389 + t367 * t398 + t293;
t368 = mrSges(4,1) * t398 - mrSges(4,3) * t389;
t417 = t429 * t296 - t298 * t402;
t292 = m(4) * t330 - mrSges(4,2) * t387 + mrSges(4,3) * t362 + t366 * t388 - t368 * t398 + t417;
t418 = -t291 * t403 + t406 * t292;
t349 = -qJDD(2) * pkin(2) - pkin(7) * t409 + t391 * t424 - t369;
t317 = -pkin(3) * t362 - pkin(8) * t386 + t389 * t371 + t349;
t309 = -0.2e1 * qJD(5) * t365 + (t364 * t397 - t328) * qJ(5) + (t365 * t397 + t327) * pkin(4) + t317;
t299 = m(6) * t309 + t327 * mrSges(6,1) - t328 * mrSges(6,3) + t364 * t351 - t365 * t354;
t287 = t291 * t406 + t292 * t403;
t414 = m(5) * t317 + t327 * mrSges(5,1) + t328 * mrSges(5,2) + t364 * t352 + t365 * t353 + t299;
t303 = mrSges(6,2) * t328 + t341 * t365 - t416;
t413 = -mrSges(5,1) * t311 + mrSges(6,1) * t307 + mrSges(5,2) * t312 - mrSges(6,3) * t306 + pkin(4) * t303 - qJ(5) * t421 + t432 * t385 + t431 * t365 + (qJ(5) * t341 - t430) * t364 - t434 * t328 + (qJ(5) * mrSges(6,2) + t433) * t327;
t412 = -m(4) * t349 + t362 * mrSges(4,1) - t363 * mrSges(4,2) + t388 * t367 - t389 * t368 - t414;
t357 = Ifges(4,4) * t389 + Ifges(4,2) * t388 + Ifges(4,6) * t398;
t358 = Ifges(4,1) * t389 + Ifges(4,4) * t388 + Ifges(4,5) * t398;
t411 = mrSges(4,1) * t329 - mrSges(4,2) * t330 + Ifges(4,5) * t363 + Ifges(4,6) * t362 + Ifges(4,3) * t387 + pkin(3) * t293 + t389 * t357 - t388 * t358 - t413;
t395 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t423;
t394 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t424;
t390 = (-mrSges(3,1) * t407 + mrSges(3,2) * t404) * qJD(1);
t382 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t404 + Ifges(3,4) * t407) * qJD(1);
t381 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t404 + Ifges(3,2) * t407) * qJD(1);
t356 = Ifges(4,5) * t389 + Ifges(4,6) * t388 + Ifges(4,3) * t398;
t289 = mrSges(5,2) * t317 + mrSges(6,2) * t307 - mrSges(5,3) * t311 - mrSges(6,3) * t309 - qJ(5) * t299 - t427 * t327 + t436 * t328 + t426 * t364 + t434 * t385 + t431 * t397;
t288 = -mrSges(5,1) * t317 - mrSges(6,1) * t309 + mrSges(6,2) * t306 + mrSges(5,3) * t312 - pkin(4) * t299 - t435 * t327 + t427 * t328 + t426 * t365 + t433 * t385 + t430 * t397;
t286 = mrSges(4,2) * t349 - mrSges(4,3) * t329 + Ifges(4,1) * t363 + Ifges(4,4) * t362 + Ifges(4,5) * t387 - pkin(8) * t293 - t402 * t288 + t429 * t289 + t388 * t356 - t398 * t357;
t285 = -mrSges(4,1) * t349 + mrSges(4,3) * t330 + Ifges(4,4) * t363 + Ifges(4,2) * t362 + Ifges(4,6) * t387 - pkin(3) * t414 + pkin(8) * t417 + t429 * t288 + t402 * t289 - t389 * t356 + t398 * t358;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t419 - mrSges(2,2) * t415 + t404 * (mrSges(3,2) * t383 - mrSges(3,3) * t369 + Ifges(3,1) * t392 + Ifges(3,4) * t393 + Ifges(3,5) * qJDD(2) - pkin(7) * t287 - qJD(2) * t381 - t403 * t285 + t406 * t286) + t407 * (-mrSges(3,1) * t383 + mrSges(3,3) * t370 + Ifges(3,4) * t392 + Ifges(3,2) * t393 + Ifges(3,6) * qJDD(2) - pkin(2) * t287 + qJD(2) * t382 - t411) + pkin(1) * (-m(3) * t383 + mrSges(3,1) * t393 - mrSges(3,2) * t392 + (-t394 * t404 + t395 * t407) * qJD(1) - t287) + pkin(6) * (t407 * (m(3) * t370 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t393 - qJD(2) * t394 + t390 * t423 + t418) - t404 * (m(3) * t369 + qJDD(2) * mrSges(3,1) - t392 * mrSges(3,3) + qJD(2) * t395 - t390 * t424 + t412)); Ifges(3,5) * t392 + Ifges(3,6) * t393 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t369 - mrSges(3,2) * t370 + t403 * t286 + t406 * t285 + pkin(2) * t412 + pkin(7) * t418 + (t404 * t381 - t407 * t382) * qJD(1); t411; -t413; t303;];
tauJ = t1;
