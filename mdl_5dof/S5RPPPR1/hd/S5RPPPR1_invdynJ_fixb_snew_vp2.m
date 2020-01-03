% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPPR1
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:23
% EndTime: 2020-01-03 11:20:24
% DurationCPUTime: 1.14s
% Computational Cost: add. (4605->176), mult. (10702->256), div. (0->0), fcn. (6478->10), ass. (0->94)
t364 = sin(qJ(1));
t366 = cos(qJ(1));
t382 = -g(2) * t366 - g(3) * t364;
t341 = qJDD(1) * pkin(1) + t382;
t367 = qJD(1) ^ 2;
t387 = -g(2) * t364 + t366 * g(3);
t342 = -pkin(1) * t367 + t387;
t359 = sin(pkin(7));
t362 = cos(pkin(7));
t324 = t359 * t341 + t362 * t342;
t403 = -pkin(2) * t367 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t324;
t323 = t341 * t362 - t359 * t342;
t372 = -qJ(3) * t367 + qJDD(3) - t323;
t358 = sin(pkin(8));
t361 = cos(pkin(8));
t380 = -t361 * pkin(3) - t358 * qJ(4);
t393 = t358 * qJD(1);
t402 = (-pkin(2) + t380) * qJDD(1) + t372 - 0.2e1 * qJD(4) * t393;
t356 = -g(1) + qJDD(2);
t309 = t356 * t361 - t403 * t358;
t360 = cos(pkin(9));
t401 = Ifges(5,4) * t360;
t357 = sin(pkin(9));
t400 = Ifges(5,2) * t357;
t399 = t358 * mrSges(4,2);
t355 = t358 ^ 2;
t398 = t355 * t367;
t397 = t357 * t358;
t396 = t358 * t360;
t310 = t358 * t356 + t403 * t361;
t338 = t380 * qJD(1);
t394 = qJD(1) * t361;
t303 = t338 * t394 + t310;
t378 = -t361 * pkin(4) - pkin(6) * t396;
t395 = t402 * t360;
t295 = t378 * qJDD(1) + (-t303 + (-pkin(4) * t355 * t360 + pkin(6) * t358 * t361) * t367) * t357 + t395;
t298 = t360 * t303 + t402 * t357;
t337 = t378 * qJD(1);
t389 = t357 ^ 2 * t398;
t391 = qJDD(1) * t358;
t296 = -pkin(6) * t357 * t391 - pkin(4) * t389 + t337 * t394 + t298;
t363 = sin(qJ(5));
t365 = cos(qJ(5));
t293 = t295 * t365 - t296 * t363;
t374 = t358 * (-t357 * t365 - t360 * t363);
t328 = qJD(1) * t374;
t373 = t358 * (-t357 * t363 + t360 * t365);
t329 = qJD(1) * t373;
t313 = -mrSges(6,1) * t328 + mrSges(6,2) * t329;
t316 = qJD(5) * t328 + qJDD(1) * t373;
t346 = qJD(5) - t394;
t321 = -mrSges(6,2) * t346 + mrSges(6,3) * t328;
t390 = qJDD(1) * t361;
t345 = qJDD(5) - t390;
t291 = m(6) * t293 + mrSges(6,1) * t345 - mrSges(6,3) * t316 - t313 * t329 + t321 * t346;
t294 = t295 * t363 + t296 * t365;
t315 = -qJD(5) * t329 + qJDD(1) * t374;
t322 = mrSges(6,1) * t346 - mrSges(6,3) * t329;
t292 = m(6) * t294 - mrSges(6,2) * t345 + mrSges(6,3) * t315 + t313 * t328 - t322 * t346;
t284 = t365 * t291 + t363 * t292;
t297 = -t303 * t357 + t395;
t381 = t357 * mrSges(5,1) + t360 * mrSges(5,2);
t330 = t381 * t393;
t375 = t361 * mrSges(5,2) - mrSges(5,3) * t397;
t332 = t375 * qJD(1);
t376 = -t361 * mrSges(5,1) - mrSges(5,3) * t396;
t282 = m(5) * t297 + t376 * qJDD(1) + (-t330 * t396 - t332 * t361) * qJD(1) + t284;
t333 = t376 * qJD(1);
t385 = -t291 * t363 + t365 * t292;
t283 = m(5) * t298 + t375 * qJDD(1) + (-t330 * t397 + t333 * t361) * qJD(1) + t385;
t339 = (-t361 * mrSges(4,1) + t399) * qJD(1);
t279 = m(4) * t310 - t282 * t357 + t283 * t360 + (qJDD(1) * mrSges(4,3) + qJD(1) * t339) * t361;
t302 = t338 * t393 + qJDD(4) - t309;
t300 = -pkin(6) * t389 + (pkin(4) * qJDD(1) * t357 + qJD(1) * t337 * t360) * t358 + t302;
t371 = -m(6) * t300 + mrSges(6,1) * t315 - t316 * mrSges(6,2) + t321 * t328 - t329 * t322;
t370 = m(5) * t302 - t371;
t379 = t332 * t357 + t333 * t360;
t287 = m(4) * t309 + ((-mrSges(4,3) - t381) * qJDD(1) + (-t339 - t379) * qJD(1)) * t358 - t370;
t386 = t361 * t279 - t287 * t358;
t281 = t282 * t360 + t283 * t357;
t377 = -Ifges(5,5) * t360 + Ifges(5,6) * t357 + Ifges(4,4);
t319 = -qJDD(1) * pkin(2) + t372;
t369 = -m(4) * t319 + mrSges(4,1) * t390 - t281 + (t361 ^ 2 * t367 + t398) * mrSges(4,3);
t305 = Ifges(6,4) * t329 + Ifges(6,2) * t328 + Ifges(6,6) * t346;
t306 = Ifges(6,1) * t329 + Ifges(6,4) * t328 + Ifges(6,5) * t346;
t368 = mrSges(6,1) * t293 - mrSges(6,2) * t294 + Ifges(6,5) * t316 + Ifges(6,6) * t315 + Ifges(6,3) * t345 + t329 * t305 - t328 * t306;
t340 = (Ifges(4,5) * t358 + Ifges(4,6) * t361) * qJD(1);
t327 = (-Ifges(5,5) * t361 + (Ifges(5,1) * t360 - Ifges(5,4) * t357) * t358) * qJD(1);
t326 = (-Ifges(5,6) * t361 + (-t400 + t401) * t358) * qJD(1);
t304 = Ifges(6,5) * t329 + Ifges(6,6) * t328 + Ifges(6,3) * t346;
t286 = mrSges(6,2) * t300 - mrSges(6,3) * t293 + Ifges(6,1) * t316 + Ifges(6,4) * t315 + Ifges(6,5) * t345 + t304 * t328 - t305 * t346;
t285 = -mrSges(6,1) * t300 + mrSges(6,3) * t294 + Ifges(6,4) * t316 + Ifges(6,2) * t315 + Ifges(6,6) * t345 - t304 * t329 + t306 * t346;
t280 = mrSges(4,2) * t391 - t369;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t382 - mrSges(2,2) * t387 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t323 - mrSges(3,2) * t324 + t358 * (t340 * t394 + mrSges(4,2) * t319 - mrSges(4,3) * t309 + t360 * (mrSges(5,2) * t302 - mrSges(5,3) * t297 - pkin(6) * t284 - t363 * t285 + t365 * t286 + t326 * t394) - t357 * (-mrSges(5,1) * t302 + mrSges(5,3) * t298 + pkin(4) * t371 + pkin(6) * t385 + t365 * t285 + t363 * t286 - t327 * t394) - qJ(4) * t281 + (t377 * t361 + (Ifges(5,1) * t360 ^ 2 + Ifges(4,1) + (t400 - 0.2e1 * t401) * t357) * t358) * qJDD(1)) + t361 * (-mrSges(4,1) * t319 - mrSges(5,1) * t297 + mrSges(5,2) * t298 + mrSges(4,3) * t310 - pkin(3) * t281 - pkin(4) * t284 + (Ifges(4,2) + Ifges(5,3)) * t390 + (t377 * qJDD(1) + (-t326 * t360 - t327 * t357 - t340) * qJD(1)) * t358 - t368) - pkin(2) * t280 + qJ(3) * t386 + pkin(1) * (t359 * (m(3) * t324 - mrSges(3,1) * t367 - qJDD(1) * mrSges(3,2) + t386) + t362 * (m(3) * t323 - mrSges(3,2) * t367 + (mrSges(3,1) - t399) * qJDD(1) + t369)); m(3) * t356 + t279 * t358 + t287 * t361; t280; (t379 * qJD(1) + t381 * qJDD(1)) * t358 + t370; t368;];
tauJ = t1;
