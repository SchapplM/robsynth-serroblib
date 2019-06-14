% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:40:30
% EndTime: 2019-05-05 15:40:33
% DurationCPUTime: 1.53s
% Computational Cost: add. (16051->238), mult. (29195->296), div. (0->0), fcn. (15973->10), ass. (0->102)
t410 = -pkin(1) - pkin(2);
t392 = qJD(1) ^ 2;
t386 = sin(qJ(1));
t390 = cos(qJ(1));
t401 = -t390 * g(1) - t386 * g(2);
t399 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t401;
t350 = t410 * t392 + t399;
t405 = t386 * g(1) - t390 * g(2);
t398 = -t392 * qJ(2) + qJDD(2) - t405;
t351 = t410 * qJDD(1) + t398;
t381 = sin(pkin(10));
t382 = cos(pkin(10));
t330 = -t381 * t350 + t382 * t351;
t328 = qJDD(1) * pkin(3) - t392 * pkin(7) - t330;
t385 = sin(qJ(4));
t389 = cos(qJ(4));
t408 = qJD(1) * qJD(4);
t406 = t389 * t408;
t367 = -t385 * qJDD(1) - t406;
t407 = t385 * t408;
t368 = -t389 * qJDD(1) + t407;
t315 = (-t367 + t406) * pkin(8) + (-t368 - t407) * pkin(4) + t328;
t331 = t382 * t350 + t381 * t351;
t329 = -t392 * pkin(3) - qJDD(1) * pkin(7) + t331;
t379 = g(3) + qJDD(3);
t325 = t389 * t329 + t385 * t379;
t366 = (t389 * pkin(4) + t385 * pkin(8)) * qJD(1);
t375 = t389 * qJD(1);
t391 = qJD(4) ^ 2;
t318 = -t391 * pkin(4) + qJDD(4) * pkin(8) - t366 * t375 + t325;
t384 = sin(qJ(5));
t388 = cos(qJ(5));
t304 = t388 * t315 - t384 * t318;
t409 = t385 * qJD(1);
t363 = t388 * qJD(4) + t384 * t409;
t340 = t363 * qJD(5) + t384 * qJDD(4) + t388 * t367;
t362 = qJDD(5) - t368;
t364 = t384 * qJD(4) - t388 * t409;
t372 = t375 + qJD(5);
t302 = (t363 * t372 - t340) * pkin(9) + (t363 * t364 + t362) * pkin(5) + t304;
t305 = t384 * t315 + t388 * t318;
t339 = -t364 * qJD(5) + t388 * qJDD(4) - t384 * t367;
t349 = t372 * pkin(5) - t364 * pkin(9);
t361 = t363 ^ 2;
t303 = -t361 * pkin(5) + t339 * pkin(9) - t372 * t349 + t305;
t383 = sin(qJ(6));
t387 = cos(qJ(6));
t300 = t387 * t302 - t383 * t303;
t341 = t387 * t363 - t383 * t364;
t312 = t341 * qJD(6) + t383 * t339 + t387 * t340;
t342 = t383 * t363 + t387 * t364;
t323 = -t341 * mrSges(7,1) + t342 * mrSges(7,2);
t371 = qJD(6) + t372;
t332 = -t371 * mrSges(7,2) + t341 * mrSges(7,3);
t358 = qJDD(6) + t362;
t297 = m(7) * t300 + t358 * mrSges(7,1) - t312 * mrSges(7,3) - t342 * t323 + t371 * t332;
t301 = t383 * t302 + t387 * t303;
t311 = -t342 * qJD(6) + t387 * t339 - t383 * t340;
t333 = t371 * mrSges(7,1) - t342 * mrSges(7,3);
t298 = m(7) * t301 - t358 * mrSges(7,2) + t311 * mrSges(7,3) + t341 * t323 - t371 * t333;
t290 = t387 * t297 + t383 * t298;
t365 = (t389 * mrSges(5,1) - t385 * mrSges(5,2)) * qJD(1);
t369 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t409;
t343 = -t363 * mrSges(6,1) + t364 * mrSges(6,2);
t347 = -t372 * mrSges(6,2) + t363 * mrSges(6,3);
t288 = m(6) * t304 + t362 * mrSges(6,1) - t340 * mrSges(6,3) - t364 * t343 + t372 * t347 + t290;
t348 = t372 * mrSges(6,1) - t364 * mrSges(6,3);
t402 = -t383 * t297 + t387 * t298;
t289 = m(6) * t305 - t362 * mrSges(6,2) + t339 * mrSges(6,3) + t363 * t343 - t372 * t348 + t402;
t403 = -t384 * t288 + t388 * t289;
t285 = m(5) * t325 - qJDD(4) * mrSges(5,2) + t368 * mrSges(5,3) - qJD(4) * t369 - t365 * t375 + t403;
t324 = -t385 * t329 + t389 * t379;
t370 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t375;
t317 = -qJDD(4) * pkin(4) - t391 * pkin(8) - t366 * t409 - t324;
t306 = -t339 * pkin(5) - t361 * pkin(9) + t364 * t349 + t317;
t397 = m(7) * t306 - t311 * mrSges(7,1) + t312 * mrSges(7,2) - t341 * t332 + t342 * t333;
t394 = -m(6) * t317 + t339 * mrSges(6,1) - t340 * mrSges(6,2) + t363 * t347 - t364 * t348 - t397;
t293 = m(5) * t324 + qJDD(4) * mrSges(5,1) - t367 * mrSges(5,3) + qJD(4) * t370 + t365 * t409 + t394;
t404 = t389 * t285 - t385 * t293;
t280 = m(4) * t331 - t392 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t404;
t286 = t388 * t288 + t384 * t289;
t395 = -m(5) * t328 + t368 * mrSges(5,1) - t367 * mrSges(5,2) + t369 * t409 - t370 * t375 - t286;
t283 = m(4) * t330 - qJDD(1) * mrSges(4,1) - t392 * mrSges(4,2) + t395;
t400 = t381 * t280 + t382 * t283;
t320 = Ifges(7,4) * t342 + Ifges(7,2) * t341 + Ifges(7,6) * t371;
t321 = Ifges(7,1) * t342 + Ifges(7,4) * t341 + Ifges(7,5) * t371;
t396 = -mrSges(7,1) * t300 + mrSges(7,2) * t301 - Ifges(7,5) * t312 - Ifges(7,6) * t311 - Ifges(7,3) * t358 - t342 * t320 + t341 * t321;
t335 = Ifges(6,4) * t364 + Ifges(6,2) * t363 + Ifges(6,6) * t372;
t336 = Ifges(6,1) * t364 + Ifges(6,4) * t363 + Ifges(6,5) * t372;
t393 = mrSges(6,1) * t304 - mrSges(6,2) * t305 + Ifges(6,5) * t340 + Ifges(6,6) * t339 + Ifges(6,3) * t362 + pkin(5) * t290 + t364 * t335 - t363 * t336 - t396;
t357 = (Ifges(5,5) * qJD(4)) + (-t385 * Ifges(5,1) - t389 * Ifges(5,4)) * qJD(1);
t356 = (Ifges(5,6) * qJD(4)) + (-t385 * Ifges(5,4) - t389 * Ifges(5,2)) * qJD(1);
t353 = -qJDD(1) * pkin(1) + t398;
t352 = -t392 * pkin(1) + t399;
t334 = Ifges(6,5) * t364 + Ifges(6,6) * t363 + Ifges(6,3) * t372;
t319 = Ifges(7,5) * t342 + Ifges(7,6) * t341 + Ifges(7,3) * t371;
t292 = mrSges(7,2) * t306 - mrSges(7,3) * t300 + Ifges(7,1) * t312 + Ifges(7,4) * t311 + Ifges(7,5) * t358 + t341 * t319 - t371 * t320;
t291 = -mrSges(7,1) * t306 + mrSges(7,3) * t301 + Ifges(7,4) * t312 + Ifges(7,2) * t311 + Ifges(7,6) * t358 - t342 * t319 + t371 * t321;
t282 = mrSges(6,2) * t317 - mrSges(6,3) * t304 + Ifges(6,1) * t340 + Ifges(6,4) * t339 + Ifges(6,5) * t362 - pkin(9) * t290 - t383 * t291 + t387 * t292 + t363 * t334 - t372 * t335;
t281 = -mrSges(6,1) * t317 + mrSges(6,3) * t305 + Ifges(6,4) * t340 + Ifges(6,2) * t339 + Ifges(6,6) * t362 - pkin(5) * t397 + pkin(9) * t402 + t387 * t291 + t383 * t292 - t364 * t334 + t372 * t336;
t279 = m(3) * t353 - qJDD(1) * mrSges(3,1) - t392 * mrSges(3,3) + t400;
t1 = [-pkin(1) * t279 + qJ(2) * (m(3) * t352 - t392 * mrSges(3,1) + t382 * t280 - t381 * t283) + mrSges(2,1) * t405 - mrSges(2,2) * t401 - pkin(2) * t400 - mrSges(3,1) * t353 + mrSges(3,3) * t352 - t389 * (-mrSges(5,1) * t328 + mrSges(5,3) * t325 + Ifges(5,4) * t367 + Ifges(5,2) * t368 + Ifges(5,6) * qJDD(4) - pkin(4) * t286 + qJD(4) * t357 - t393) - pkin(3) * t395 - pkin(7) * t404 - t385 * (mrSges(5,2) * t328 - mrSges(5,3) * t324 + Ifges(5,1) * t367 + Ifges(5,4) * t368 + Ifges(5,5) * qJDD(4) - pkin(8) * t286 - qJD(4) * t356 - t384 * t281 + t388 * t282) - mrSges(4,1) * t330 + mrSges(4,2) * t331 + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1); t279; m(4) * t379 + t385 * t285 + t389 * t293; Ifges(5,5) * t367 + Ifges(5,6) * t368 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t324 - mrSges(5,2) * t325 + t384 * t282 + t388 * t281 + pkin(4) * t394 + pkin(8) * t403 + (-t385 * t356 + t389 * t357) * qJD(1); t393; -t396;];
tauJ  = t1;
