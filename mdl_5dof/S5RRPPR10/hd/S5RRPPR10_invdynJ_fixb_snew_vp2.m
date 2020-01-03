% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:16
% EndTime: 2019-12-31 19:43:19
% DurationCPUTime: 1.68s
% Computational Cost: add. (8149->258), mult. (17819->313), div. (0->0), fcn. (11018->8), ass. (0->103)
t454 = -2 * qJD(3);
t453 = Ifges(4,1) + Ifges(5,1);
t445 = Ifges(4,4) - Ifges(5,5);
t452 = Ifges(4,5) + Ifges(5,4);
t451 = Ifges(4,2) + Ifges(5,3);
t450 = Ifges(5,2) + Ifges(4,3);
t449 = Ifges(4,6) - Ifges(5,6);
t419 = qJD(1) ^ 2;
t414 = sin(qJ(1));
t417 = cos(qJ(1));
t430 = g(1) * t414 - t417 * g(2);
t387 = -qJDD(1) * pkin(1) - pkin(6) * t419 - t430;
t413 = sin(qJ(2));
t416 = cos(qJ(2));
t433 = qJD(1) * qJD(2);
t431 = t416 * t433;
t399 = qJDD(1) * t413 + t431;
t405 = t413 * t433;
t400 = qJDD(1) * t416 - t405;
t346 = (-t399 - t431) * qJ(3) + (-t400 + t405) * pkin(2) + t387;
t427 = -g(1) * t417 - g(2) * t414;
t388 = -pkin(1) * t419 + qJDD(1) * pkin(6) + t427;
t369 = -g(3) * t413 + t416 * t388;
t397 = (-pkin(2) * t416 - qJ(3) * t413) * qJD(1);
t418 = qJD(2) ^ 2;
t435 = qJD(1) * t416;
t350 = -pkin(2) * t418 + qJDD(2) * qJ(3) + t397 * t435 + t369;
t411 = sin(pkin(8));
t436 = qJD(1) * t413;
t442 = cos(pkin(8));
t393 = t411 * qJD(2) + t442 * t436;
t330 = t442 * t346 - t411 * t350 + t393 * t454;
t377 = t411 * qJDD(2) + t442 * t399;
t368 = -t416 * g(3) - t413 * t388;
t421 = qJDD(2) * pkin(2) + qJ(3) * t418 - t397 * t436 - qJDD(3) + t368;
t392 = -t442 * qJD(2) + t411 * t436;
t432 = t392 * t435;
t448 = -(t377 + t432) * qJ(4) - t421;
t447 = -2 * qJD(4);
t446 = -mrSges(4,3) - mrSges(5,2);
t441 = t416 ^ 2 * t419;
t440 = t392 * t451 - t393 * t445 + t435 * t449;
t439 = t392 * t449 - t393 * t452 + t435 * t450;
t438 = t392 * t445 - t393 * t453 + t435 * t452;
t366 = mrSges(5,1) * t392 - mrSges(5,3) * t393;
t437 = -mrSges(4,1) * t392 - mrSges(4,2) * t393 - t366;
t331 = t411 * t346 + t442 * t350 + t392 * t454;
t374 = -mrSges(4,1) * t435 - mrSges(4,3) * t393;
t375 = mrSges(5,1) * t435 + mrSges(5,2) * t393;
t376 = -t442 * qJDD(2) + t399 * t411;
t365 = pkin(3) * t392 - qJ(4) * t393;
t327 = -pkin(3) * t441 - t400 * qJ(4) - t392 * t365 + t435 * t447 + t331;
t328 = t400 * pkin(3) - qJ(4) * t441 + t393 * t365 + qJDD(4) - t330;
t322 = (-t377 + t432) * pkin(7) + (t392 * t393 + t400) * pkin(4) + t328;
t378 = pkin(4) * t435 - pkin(7) * t393;
t390 = t392 ^ 2;
t323 = -pkin(4) * t390 + pkin(7) * t376 - t378 * t435 + t327;
t412 = sin(qJ(5));
t415 = cos(qJ(5));
t320 = t322 * t415 - t323 * t412;
t363 = t392 * t415 - t393 * t412;
t338 = qJD(5) * t363 + t376 * t412 + t377 * t415;
t364 = t392 * t412 + t393 * t415;
t344 = -mrSges(6,1) * t363 + mrSges(6,2) * t364;
t403 = qJD(5) + t435;
t351 = -mrSges(6,2) * t403 + mrSges(6,3) * t363;
t396 = qJDD(5) + t400;
t317 = m(6) * t320 + mrSges(6,1) * t396 - mrSges(6,3) * t338 - t344 * t364 + t351 * t403;
t321 = t322 * t412 + t323 * t415;
t337 = -qJD(5) * t364 + t376 * t415 - t377 * t412;
t352 = mrSges(6,1) * t403 - mrSges(6,3) * t364;
t318 = m(6) * t321 - mrSges(6,2) * t396 + mrSges(6,3) * t337 + t344 * t363 - t352 * t403;
t428 = -t317 * t412 + t415 * t318;
t423 = m(5) * t327 - t400 * mrSges(5,3) + t428;
t308 = m(4) * t331 + mrSges(4,2) * t400 + t437 * t392 + t446 * t376 + (t374 - t375) * t435 + t423;
t372 = -mrSges(5,2) * t392 - mrSges(5,3) * t435;
t373 = mrSges(4,2) * t435 - mrSges(4,3) * t392;
t311 = t317 * t415 + t318 * t412;
t422 = -m(5) * t328 - t400 * mrSges(5,1) - t311;
t309 = m(4) * t330 - mrSges(4,1) * t400 + t437 * t393 + t446 * t377 + (-t372 - t373) * t435 + t422;
t429 = t442 * t308 - t309 * t411;
t325 = -pkin(7) * t390 + (-pkin(3) - pkin(4)) * t376 + (pkin(3) * t435 + (2 * qJD(4)) + t378) * t393 - t448;
t424 = -m(6) * t325 + t337 * mrSges(6,1) - t338 * mrSges(6,2) + t363 * t351 - t364 * t352;
t306 = t411 * t308 + t442 * t309;
t340 = Ifges(6,4) * t364 + Ifges(6,2) * t363 + Ifges(6,6) * t403;
t341 = Ifges(6,1) * t364 + Ifges(6,4) * t363 + Ifges(6,5) * t403;
t420 = mrSges(6,1) * t320 - mrSges(6,2) * t321 + Ifges(6,5) * t338 + Ifges(6,6) * t337 + Ifges(6,3) * t396 + t364 * t340 - t363 * t341;
t329 = t393 * t447 + (-t393 * t435 + t376) * pkin(3) + t448;
t315 = m(5) * t329 + t376 * mrSges(5,1) - t377 * mrSges(5,3) + t392 * t372 - t393 * t375 + t424;
t314 = -m(4) * t421 + t376 * mrSges(4,1) + t377 * mrSges(4,2) + t392 * t373 + t393 * t374 + t315;
t402 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t435;
t401 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t436;
t398 = (-mrSges(3,1) * t416 + mrSges(3,2) * t413) * qJD(1);
t386 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t413 + Ifges(3,4) * t416) * qJD(1);
t385 = Ifges(3,6) * qJD(2) + (t413 * Ifges(3,4) + Ifges(3,2) * t416) * qJD(1);
t384 = Ifges(3,3) * qJD(2) + (t413 * Ifges(3,5) + t416 * Ifges(3,6)) * qJD(1);
t339 = Ifges(6,5) * t364 + Ifges(6,6) * t363 + Ifges(6,3) * t403;
t313 = mrSges(6,2) * t325 - mrSges(6,3) * t320 + Ifges(6,1) * t338 + Ifges(6,4) * t337 + Ifges(6,5) * t396 + t339 * t363 - t340 * t403;
t312 = -mrSges(6,1) * t325 + mrSges(6,3) * t321 + Ifges(6,4) * t338 + Ifges(6,2) * t337 + Ifges(6,6) * t396 - t339 * t364 + t341 * t403;
t310 = mrSges(5,2) * t377 + t366 * t393 + t372 * t435 - t422;
t305 = -mrSges(4,2) * t421 + mrSges(5,2) * t328 - mrSges(4,3) * t330 - mrSges(5,3) * t329 - pkin(7) * t311 - qJ(4) * t315 - t412 * t312 + t415 * t313 - t445 * t376 + t377 * t453 + t439 * t392 - t400 * t452 - t440 * t435;
t304 = mrSges(4,1) * t421 - mrSges(5,1) * t329 + mrSges(5,2) * t327 + mrSges(4,3) * t331 - pkin(3) * t315 - pkin(4) * t424 - pkin(7) * t428 - t415 * t312 - t412 * t313 - t376 * t451 + t445 * t377 + t439 * t393 - t400 * t449 + t438 * t435;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t430 - mrSges(2,2) * t427 + t413 * (mrSges(3,2) * t387 - mrSges(3,3) * t368 + Ifges(3,1) * t399 + Ifges(3,4) * t400 + Ifges(3,5) * qJDD(2) - qJ(3) * t306 - qJD(2) * t385 - t411 * t304 + t442 * t305 + t384 * t435) + t416 * (t420 + Ifges(3,6) * qJDD(2) - qJ(4) * (-t375 * t435 + t423) + (Ifges(3,2) + t450) * t400 + t440 * t393 + (qJ(4) * t366 + t438) * t392 - t452 * t377 + (qJ(4) * mrSges(5,2) + t449) * t376 - pkin(2) * t306 - t384 * t436 + pkin(3) * t310 + pkin(4) * t311 - mrSges(5,3) * t327 + mrSges(5,1) * t328 - mrSges(4,1) * t330 + mrSges(4,2) * t331 + mrSges(3,3) * t369 + qJD(2) * t386 - mrSges(3,1) * t387 + Ifges(3,4) * t399) + pkin(1) * (-m(3) * t387 + t400 * mrSges(3,1) - t399 * mrSges(3,2) + (-t401 * t413 + t402 * t416) * qJD(1) - t306) + pkin(6) * (t416 * (m(3) * t369 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t400 - qJD(2) * t401 + t398 * t435 + t429) - t413 * (m(3) * t368 + qJDD(2) * mrSges(3,1) - t399 * mrSges(3,3) + qJD(2) * t402 - t398 * t436 - t314)); Ifges(3,5) * t399 + Ifges(3,6) * t400 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t368 - mrSges(3,2) * t369 + t411 * t305 + t442 * t304 - pkin(2) * t314 + qJ(3) * t429 + (t413 * t385 - t416 * t386) * qJD(1); t314; t310; t420;];
tauJ = t1;
