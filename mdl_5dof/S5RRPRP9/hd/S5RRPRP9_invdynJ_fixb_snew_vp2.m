% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:59
% EndTime: 2019-12-31 20:06:02
% DurationCPUTime: 1.95s
% Computational Cost: add. (11940->249), mult. (25619->308), div. (0->0), fcn. (16871->8), ass. (0->98)
t428 = Ifges(5,1) + Ifges(6,1);
t418 = Ifges(5,4) - Ifges(6,5);
t424 = -Ifges(5,5) - Ifges(6,4);
t427 = Ifges(5,2) + Ifges(6,3);
t416 = Ifges(5,6) - Ifges(6,6);
t391 = sin(pkin(8));
t392 = cos(pkin(8));
t394 = sin(qJ(2));
t410 = t394 * qJD(1);
t377 = qJD(2) * t392 - t391 * t410;
t378 = qJD(2) * t391 + t392 * t410;
t393 = sin(qJ(4));
t420 = cos(qJ(4));
t351 = -t377 * t420 + t378 * t393;
t396 = cos(qJ(2));
t409 = qJD(1) * qJD(2);
t407 = t396 * t409;
t382 = qJDD(1) * t394 + t407;
t360 = qJDD(2) * t392 - t382 * t391;
t361 = qJDD(2) * t391 + t382 * t392;
t321 = -t351 * qJD(4) + t393 * t360 + t361 * t420;
t352 = t393 * t377 + t378 * t420;
t332 = mrSges(6,1) * t351 - mrSges(6,3) * t352;
t399 = qJD(1) ^ 2;
t395 = sin(qJ(1));
t397 = cos(qJ(1));
t406 = g(1) * t395 - t397 * g(2);
t372 = -qJDD(1) * pkin(1) - pkin(6) * t399 - t406;
t388 = t394 * t409;
t383 = qJDD(1) * t396 - t388;
t336 = (-t382 - t407) * qJ(3) + (-t383 + t388) * pkin(2) + t372;
t402 = -g(1) * t397 - g(2) * t395;
t373 = -pkin(1) * t399 + qJDD(1) * pkin(6) + t402;
t356 = -g(3) * t394 + t396 * t373;
t380 = (-pkin(2) * t396 - qJ(3) * t394) * qJD(1);
t398 = qJD(2) ^ 2;
t411 = qJD(1) * t396;
t341 = -pkin(2) * t398 + qJDD(2) * qJ(3) + t380 * t411 + t356;
t314 = -0.2e1 * qJD(3) * t378 + t392 * t336 - t341 * t391;
t311 = (-t377 * t411 - t361) * pkin(7) + (t377 * t378 - t383) * pkin(3) + t314;
t315 = 0.2e1 * qJD(3) * t377 + t391 * t336 + t392 * t341;
t362 = -pkin(3) * t411 - pkin(7) * t378;
t376 = t377 ^ 2;
t313 = -pkin(3) * t376 + pkin(7) * t360 + t362 * t411 + t315;
t308 = t311 * t420 - t393 * t313;
t331 = pkin(4) * t351 - qJ(5) * t352;
t379 = qJDD(4) - t383;
t387 = qJD(4) - t411;
t386 = t387 ^ 2;
t305 = -t379 * pkin(4) - t386 * qJ(5) + t352 * t331 + qJDD(5) - t308;
t345 = -mrSges(6,2) * t351 + mrSges(6,3) * t387;
t403 = -m(6) * t305 + t379 * mrSges(6,1) + t387 * t345;
t301 = mrSges(6,2) * t321 + t332 * t352 - t403;
t309 = t393 * t311 + t420 * t313;
t304 = -pkin(4) * t386 + qJ(5) * t379 + 0.2e1 * qJD(5) * t387 - t331 * t351 + t309;
t320 = qJD(4) * t352 - t360 * t420 + t361 * t393;
t344 = -mrSges(6,1) * t387 + mrSges(6,2) * t352;
t408 = m(6) * t304 + t379 * mrSges(6,3) + t387 * t344;
t413 = -t418 * t351 + t428 * t352 - t424 * t387;
t414 = t427 * t351 - t418 * t352 - t416 * t387;
t423 = -Ifges(5,3) - Ifges(6,2);
t426 = -t424 * t321 - t416 * t320 - t423 * t379 + mrSges(5,1) * t308 - mrSges(6,1) * t305 - mrSges(5,2) * t309 + mrSges(6,3) * t304 - pkin(4) * t301 + qJ(5) * (-mrSges(6,2) * t320 - t332 * t351 + t408) - t414 * t352 + t413 * t351;
t419 = -mrSges(5,3) - mrSges(6,2);
t343 = mrSges(5,1) * t387 - mrSges(5,3) * t352;
t412 = -mrSges(5,1) * t351 - mrSges(5,2) * t352 - t332;
t296 = m(5) * t309 - mrSges(5,2) * t379 + t320 * t419 - t343 * t387 + t412 * t351 + t408;
t342 = -mrSges(5,2) * t387 - mrSges(5,3) * t351;
t298 = m(5) * t308 + mrSges(5,1) * t379 + t321 * t419 + t342 * t387 + t412 * t352 + t403;
t293 = t393 * t296 + t420 * t298;
t415 = t416 * t351 + t424 * t352 + t423 * t387;
t355 = -t396 * g(3) - t394 * t373;
t353 = -mrSges(4,1) * t377 + mrSges(4,2) * t378;
t358 = mrSges(4,2) * t411 + mrSges(4,3) * t377;
t289 = m(4) * t314 - mrSges(4,1) * t383 - mrSges(4,3) * t361 - t353 * t378 - t358 * t411 + t293;
t359 = -mrSges(4,1) * t411 - mrSges(4,3) * t378;
t404 = t420 * t296 - t298 * t393;
t290 = m(4) * t315 + mrSges(4,2) * t383 + mrSges(4,3) * t360 + t353 * t377 + t359 * t411 + t404;
t405 = -t289 * t391 + t392 * t290;
t340 = -qJDD(2) * pkin(2) - qJ(3) * t398 + t380 * t410 + qJDD(3) - t355;
t316 = -pkin(3) * t360 - pkin(7) * t376 + t378 * t362 + t340;
t307 = -0.2e1 * qJD(5) * t352 + (t351 * t387 - t321) * qJ(5) + (t352 * t387 + t320) * pkin(4) + t316;
t302 = m(6) * t307 + t320 * mrSges(6,1) - t321 * mrSges(6,3) - t352 * t344 + t351 * t345;
t287 = t289 * t392 + t290 * t391;
t401 = m(5) * t316 + t320 * mrSges(5,1) + t321 * mrSges(5,2) + t351 * t342 + t352 * t343 + t302;
t299 = m(4) * t340 - t360 * mrSges(4,1) + t361 * mrSges(4,2) - t377 * t358 + t378 * t359 + t401;
t385 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t411;
t384 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t410;
t381 = (-t396 * mrSges(3,1) + t394 * mrSges(3,2)) * qJD(1);
t371 = Ifges(3,5) * qJD(2) + (t394 * Ifges(3,1) + t396 * Ifges(3,4)) * qJD(1);
t370 = Ifges(3,6) * qJD(2) + (t394 * Ifges(3,4) + Ifges(3,2) * t396) * qJD(1);
t348 = Ifges(4,1) * t378 + Ifges(4,4) * t377 - Ifges(4,5) * t411;
t347 = Ifges(4,4) * t378 + Ifges(4,2) * t377 - Ifges(4,6) * t411;
t346 = Ifges(4,5) * t378 + Ifges(4,6) * t377 - Ifges(4,3) * t411;
t292 = mrSges(5,2) * t316 + mrSges(6,2) * t305 - mrSges(5,3) * t308 - mrSges(6,3) * t307 - qJ(5) * t302 - t418 * t320 + t428 * t321 + t415 * t351 - t424 * t379 + t414 * t387;
t291 = -mrSges(5,1) * t316 - mrSges(6,1) * t307 + mrSges(6,2) * t304 + mrSges(5,3) * t309 - pkin(4) * t302 - t427 * t320 + t418 * t321 + t415 * t352 + t416 * t379 + t413 * t387;
t286 = mrSges(4,2) * t340 - mrSges(4,3) * t314 + Ifges(4,1) * t361 + Ifges(4,4) * t360 - Ifges(4,5) * t383 - pkin(7) * t293 - t393 * t291 + t292 * t420 + t377 * t346 + t347 * t411;
t285 = -mrSges(4,1) * t340 + mrSges(4,3) * t315 + Ifges(4,4) * t361 + Ifges(4,2) * t360 - Ifges(4,6) * t383 - pkin(3) * t401 + pkin(7) * t404 + t291 * t420 + t393 * t292 - t378 * t346 - t348 * t411;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t406 - mrSges(2,2) * t402 + t394 * (mrSges(3,2) * t372 - mrSges(3,3) * t355 + Ifges(3,1) * t382 + Ifges(3,4) * t383 + Ifges(3,5) * qJDD(2) - qJ(3) * t287 - qJD(2) * t370 - t391 * t285 + t392 * t286) + t396 * (Ifges(3,6) * qJDD(2) + Ifges(3,4) * t382 + t377 * t348 - t378 * t347 - Ifges(4,6) * t360 - Ifges(4,5) * t361 + qJD(2) * t371 - mrSges(3,1) * t372 + mrSges(3,3) * t356 + (Ifges(4,3) + Ifges(3,2)) * t383 - mrSges(4,1) * t314 + mrSges(4,2) * t315 - pkin(3) * t293 - pkin(2) * t287 - t426) + pkin(1) * (-m(3) * t372 + mrSges(3,1) * t383 - mrSges(3,2) * t382 + (-t384 * t394 + t385 * t396) * qJD(1) - t287) + pkin(6) * (t396 * (m(3) * t356 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t383 - qJD(2) * t384 + t381 * t411 + t405) - t394 * (m(3) * t355 + qJDD(2) * mrSges(3,1) - t382 * mrSges(3,3) + qJD(2) * t385 - t381 * t410 - t299)); Ifges(3,5) * t382 + Ifges(3,6) * t383 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t355 - mrSges(3,2) * t356 + t391 * t286 + t392 * t285 - pkin(2) * t299 + qJ(3) * t405 + (t394 * t370 - t396 * t371) * qJD(1); t299; t426; t301;];
tauJ = t1;
