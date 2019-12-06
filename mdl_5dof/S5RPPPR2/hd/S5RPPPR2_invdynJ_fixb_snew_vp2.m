% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:17
% EndTime: 2019-12-05 17:31:20
% DurationCPUTime: 2.11s
% Computational Cost: add. (8382->231), mult. (23494->327), div. (0->0), fcn. (15854->10), ass. (0->111)
t398 = sin(pkin(7));
t401 = cos(pkin(7));
t406 = qJD(1) ^ 2;
t403 = sin(qJ(1));
t405 = cos(qJ(1));
t424 = t403 * g(2) - g(3) * t405;
t444 = -pkin(1) * t406 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t424;
t359 = -g(1) * t398 + t444 * t401;
t416 = -pkin(2) * t401 - qJ(3) * t398;
t382 = t416 * qJD(1);
t434 = qJD(1) * t401;
t350 = t382 * t434 + t359;
t397 = sin(pkin(8));
t400 = cos(pkin(8));
t419 = g(2) * t405 + g(3) * t403;
t409 = -qJ(2) * t406 + qJDD(2) - t419;
t435 = qJD(1) * t398;
t443 = (-pkin(1) + t416) * qJDD(1) + t409 - 0.2e1 * qJD(3) * t435;
t330 = -t397 * t350 + t443 * t400;
t358 = -t401 * g(1) - t444 * t398;
t396 = sin(pkin(9));
t399 = cos(pkin(9));
t438 = t398 * t400;
t410 = t396 * t438 + t399 * t401;
t371 = t410 * qJD(1);
t369 = t410 * qJDD(1);
t442 = 2 * qJD(4);
t441 = Ifges(4,4) * t400;
t395 = t401 ^ 2;
t440 = t395 * t406;
t439 = t397 * t398;
t437 = t401 * t406;
t331 = t400 * t350 + t443 * t397;
t375 = (pkin(3) * t397 - qJ(4) * t400) * t435;
t429 = t397 * t435;
t431 = qJDD(1) * t401;
t328 = -pkin(3) * t440 - qJ(4) * t431 - t375 * t429 + t331;
t349 = t382 * t435 + qJDD(3) - t358;
t336 = ((-qJDD(1) * t400 - t397 * t437) * qJ(4) + (qJDD(1) * t397 - t400 * t437) * pkin(3)) * t398 + t349;
t324 = t399 * t328 + t396 * t336 - t371 * t442;
t428 = t400 * t435;
t372 = -t396 * t434 + t399 * t428;
t351 = mrSges(5,1) * t371 + mrSges(5,2) * t372;
t357 = mrSges(5,1) * t429 - mrSges(5,3) * t372;
t352 = pkin(4) * t371 - pkin(6) * t372;
t432 = qJDD(1) * t398;
t425 = t397 * t432;
t394 = t398 ^ 2;
t430 = t397 ^ 2 * t394 * t406;
t322 = -pkin(4) * t430 + pkin(6) * t425 - t352 * t371 + t324;
t327 = pkin(3) * t431 - qJ(4) * t440 + t375 * t428 + qJDD(4) - t330;
t370 = (-t396 * t401 + t399 * t438) * qJDD(1);
t325 = (t371 * t429 - t370) * pkin(6) + (t372 * t429 + t369) * pkin(4) + t327;
t402 = sin(qJ(5));
t404 = cos(qJ(5));
t319 = -t322 * t402 + t325 * t404;
t353 = -t372 * t402 + t404 * t429;
t354 = t372 * t404 + t402 * t429;
t338 = -mrSges(6,1) * t353 + mrSges(6,2) * t354;
t340 = qJD(5) * t353 + t370 * t404 + t402 * t425;
t368 = qJD(5) + t371;
t341 = -mrSges(6,2) * t368 + mrSges(6,3) * t353;
t367 = qJDD(5) + t369;
t317 = m(6) * t319 + mrSges(6,1) * t367 - mrSges(6,3) * t340 - t338 * t354 + t341 * t368;
t320 = t322 * t404 + t325 * t402;
t339 = -qJD(5) * t354 - t370 * t402 + t404 * t425;
t342 = mrSges(6,1) * t368 - mrSges(6,3) * t354;
t318 = m(6) * t320 - mrSges(6,2) * t367 + mrSges(6,3) * t339 + t338 * t353 - t342 * t368;
t422 = -t317 * t402 + t404 * t318;
t309 = m(5) * t324 - mrSges(5,3) * t369 - t351 * t371 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t357) * t439 + t422;
t415 = t328 * t396 - t336 * t399;
t323 = -0.2e1 * qJD(4) * t372 - t415;
t356 = -mrSges(5,2) * t429 - mrSges(5,3) * t371;
t321 = -pkin(4) * t425 - pkin(6) * t430 + (t442 + t352) * t372 + t415;
t408 = -m(6) * t321 + t339 * mrSges(6,1) - mrSges(6,2) * t340 + t353 * t341 - t342 * t354;
t315 = m(5) * t323 - mrSges(5,3) * t370 - t351 * t372 + (mrSges(5,1) * qJDD(1) + qJD(1) * t356) * t439 + t408;
t305 = t396 * t309 + t399 * t315;
t423 = t399 * t309 - t315 * t396;
t420 = m(4) * t349 + t305;
t418 = -mrSges(3,1) * t401 + mrSges(3,2) * t398;
t417 = mrSges(4,1) * t397 + mrSges(4,2) * t400;
t376 = t417 * t435;
t412 = -mrSges(4,1) * t401 - mrSges(4,3) * t438;
t380 = t412 * qJD(1);
t411 = mrSges(4,2) * t401 - mrSges(4,3) * t439;
t304 = m(4) * t331 + t411 * qJDD(1) + (-t376 * t439 + t380 * t401) * qJD(1) + t423;
t311 = t317 * t404 + t318 * t402;
t310 = m(5) * t327 + t369 * mrSges(5,1) + mrSges(5,2) * t370 + t371 * t356 + t357 * t372 + t311;
t379 = t411 * qJD(1);
t306 = m(4) * t330 + t412 * qJDD(1) + (-t376 * t438 - t379 * t401) * qJD(1) - t310;
t301 = t304 * t397 + t306 * t400;
t414 = t379 * t397 + t380 * t400;
t413 = -Ifges(4,5) * t400 + Ifges(4,6) * t397 + Ifges(3,4);
t333 = Ifges(6,4) * t354 + Ifges(6,2) * t353 + Ifges(6,6) * t368;
t334 = Ifges(6,1) * t354 + Ifges(6,4) * t353 + Ifges(6,5) * t368;
t407 = mrSges(6,1) * t319 - mrSges(6,2) * t320 + Ifges(6,5) * t340 + Ifges(6,6) * t339 + Ifges(6,3) * t367 + t354 * t333 - t353 * t334;
t384 = (Ifges(3,5) * t398 + Ifges(3,6) * t401) * qJD(1);
t383 = t418 * qJD(1);
t378 = -qJDD(1) * pkin(1) + t409;
t364 = (-Ifges(4,5) * t401 + (Ifges(4,1) * t400 - Ifges(4,4) * t397) * t398) * qJD(1);
t363 = (-Ifges(4,6) * t401 + (-Ifges(4,2) * t397 + t441) * t398) * qJD(1);
t345 = Ifges(5,1) * t372 - Ifges(5,4) * t371 + Ifges(5,5) * t429;
t344 = Ifges(5,4) * t372 - Ifges(5,2) * t371 + Ifges(5,6) * t429;
t343 = Ifges(5,5) * t372 - Ifges(5,6) * t371 + Ifges(5,3) * t429;
t332 = Ifges(6,5) * t354 + Ifges(6,6) * t353 + Ifges(6,3) * t368;
t313 = mrSges(6,2) * t321 - mrSges(6,3) * t319 + Ifges(6,1) * t340 + Ifges(6,4) * t339 + Ifges(6,5) * t367 + t332 * t353 - t333 * t368;
t312 = -mrSges(6,1) * t321 + mrSges(6,3) * t320 + Ifges(6,4) * t340 + Ifges(6,2) * t339 + Ifges(6,6) * t367 - t332 * t354 + t334 * t368;
t303 = -mrSges(5,1) * t327 + mrSges(5,3) * t324 + Ifges(5,4) * t370 - Ifges(5,2) * t369 - pkin(4) * t311 - t372 * t343 + (Ifges(5,6) * qJDD(1) + qJD(1) * t345) * t439 - t407;
t302 = mrSges(5,2) * t327 - mrSges(5,3) * t323 + Ifges(5,1) * t370 - Ifges(5,4) * t369 - pkin(6) * t311 - t312 * t402 + t313 * t404 - t343 * t371 + (Ifges(5,5) * qJDD(1) - qJD(1) * t344) * t439;
t300 = m(3) * t378 + t418 * qJDD(1) + (-t394 - t395) * t406 * mrSges(3,3) + t301;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t419 - mrSges(2,2) * t424 + t398 * (t384 * t434 + mrSges(3,2) * t378 - mrSges(3,3) * t358 + t400 * (mrSges(4,2) * t349 - mrSges(4,3) * t330 - qJ(4) * t305 + t399 * t302 - t396 * t303 + t363 * t434) - t397 * (-mrSges(4,1) * t349 - mrSges(5,1) * t323 + mrSges(5,2) * t324 + mrSges(4,3) * t331 - Ifges(5,5) * t370 + Ifges(5,6) * t369 - pkin(3) * t305 - pkin(4) * t408 - pkin(6) * t422 - t404 * t312 - t402 * t313 - t372 * t344 - t371 * t345 - t364 * t434) - qJ(3) * t301 + (t401 * t413 + (Ifges(4,1) * t400 ^ 2 + Ifges(3,1) + (-0.2e1 * t441 + (Ifges(4,2) + Ifges(5,3)) * t397) * t397) * t398) * qJDD(1)) + t401 * (-mrSges(3,1) * t378 + mrSges(3,3) * t359 - mrSges(4,1) * t330 + mrSges(4,2) * t331 - t396 * t302 - t399 * t303 + pkin(3) * t310 - qJ(4) * t423 - pkin(2) * t301 + (Ifges(3,2) + Ifges(4,3)) * t431 + (t413 * qJDD(1) + (-t363 * t400 - t364 * t397 - t384) * qJD(1)) * t398) - pkin(1) * t300 + qJ(2) * ((m(3) * t359 + t400 * t304 - t397 * t306 + (qJDD(1) * mrSges(3,3) + qJD(1) * t383) * t401) * t401 + (-m(3) * t358 + (mrSges(3,3) + t417) * t432 + (t383 + t414) * t435 + t420) * t398); t300; (t414 * qJD(1) + t417 * qJDD(1)) * t398 + t420; t310; t407;];
tauJ = t1;
