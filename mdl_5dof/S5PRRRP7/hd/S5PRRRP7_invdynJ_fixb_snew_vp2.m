% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:26
% EndTime: 2019-12-05 16:54:27
% DurationCPUTime: 0.94s
% Computational Cost: add. (4020->194), mult. (7626->238), div. (0->0), fcn. (4960->10), ass. (0->86)
t391 = Ifges(5,1) + Ifges(6,1);
t382 = Ifges(5,4) + Ifges(6,4);
t381 = Ifges(5,5) + Ifges(6,5);
t390 = Ifges(5,2) + Ifges(6,2);
t380 = Ifges(5,6) + Ifges(6,6);
t352 = sin(qJ(4));
t355 = cos(qJ(4));
t353 = sin(qJ(3));
t373 = qJD(2) * t353;
t335 = qJD(3) * t355 - t352 * t373;
t356 = cos(qJ(3));
t371 = qJD(2) * qJD(3);
t367 = t356 * t371;
t339 = qJDD(2) * t353 + t367;
t312 = qJD(4) * t335 + qJDD(3) * t352 + t339 * t355;
t336 = qJD(3) * t352 + t355 * t373;
t314 = -mrSges(6,1) * t335 + mrSges(6,2) * t336;
t348 = sin(pkin(9));
t350 = cos(pkin(9));
t342 = -g(1) * t350 - g(2) * t348;
t354 = sin(qJ(2));
t357 = cos(qJ(2));
t341 = g(1) * t348 - g(2) * t350;
t347 = -g(3) + qJDD(1);
t349 = sin(pkin(5));
t351 = cos(pkin(5));
t388 = t341 * t351 + t347 * t349;
t299 = t357 * t342 + t388 * t354;
t359 = qJD(2) ^ 2;
t296 = -pkin(2) * t359 + qJDD(2) * pkin(7) + t299;
t322 = -t341 * t349 + t347 * t351;
t292 = t356 * t296 + t353 * t322;
t338 = (-t356 * pkin(3) - t353 * pkin(8)) * qJD(2);
t358 = qJD(3) ^ 2;
t372 = qJD(2) * t356;
t287 = -pkin(3) * t358 + qJDD(3) * pkin(8) + t338 * t372 + t292;
t298 = -t354 * t342 + t388 * t357;
t295 = -qJDD(2) * pkin(2) - pkin(7) * t359 - t298;
t368 = t353 * t371;
t340 = qJDD(2) * t356 - t368;
t290 = (-t339 - t367) * pkin(8) + (-t340 + t368) * pkin(3) + t295;
t282 = -t287 * t352 + t355 * t290;
t332 = qJDD(4) - t340;
t346 = qJD(4) - t372;
t279 = -0.2e1 * qJD(5) * t336 + (t335 * t346 - t312) * qJ(5) + (t335 * t336 + t332) * pkin(4) + t282;
t317 = -mrSges(6,2) * t346 + mrSges(6,3) * t335;
t370 = m(6) * t279 + t332 * mrSges(6,1) + t346 * t317;
t276 = -mrSges(6,3) * t312 - t314 * t336 + t370;
t283 = t355 * t287 + t352 * t290;
t311 = -qJD(4) * t336 + qJDD(3) * t355 - t339 * t352;
t319 = pkin(4) * t346 - qJ(5) * t336;
t331 = t335 ^ 2;
t281 = -pkin(4) * t331 + qJ(5) * t311 + 0.2e1 * qJD(5) * t335 - t319 * t346 + t283;
t376 = -t390 * t335 - t382 * t336 - t380 * t346;
t384 = t382 * t335 + t391 * t336 + t381 * t346;
t386 = Ifges(5,3) + Ifges(6,3);
t389 = mrSges(5,1) * t282 + mrSges(6,1) * t279 - mrSges(5,2) * t283 - mrSges(6,2) * t281 + pkin(4) * t276 + t380 * t311 + t381 * t312 + t386 * t332 - t384 * t335 - t376 * t336;
t383 = -mrSges(5,2) - mrSges(6,2);
t377 = -t380 * t335 - t381 * t336 - t386 * t346;
t320 = mrSges(6,1) * t346 - mrSges(6,3) * t336;
t374 = -mrSges(5,1) * t346 + mrSges(5,3) * t336 - t320;
t369 = m(6) * t281 + t311 * mrSges(6,3) + t335 * t314;
t337 = (-t356 * mrSges(4,1) + t353 * mrSges(4,2)) * qJD(2);
t343 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t373;
t315 = -mrSges(5,1) * t335 + mrSges(5,2) * t336;
t318 = -mrSges(5,2) * t346 + mrSges(5,3) * t335;
t272 = m(5) * t282 + mrSges(5,1) * t332 + t318 * t346 + (-t314 - t315) * t336 + (-mrSges(5,3) - mrSges(6,3)) * t312 + t370;
t274 = m(5) * t283 + mrSges(5,3) * t311 + t315 * t335 + t383 * t332 + t374 * t346 + t369;
t365 = -t272 * t352 + t355 * t274;
t269 = m(4) * t292 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t340 - qJD(3) * t343 + t337 * t372 + t365;
t291 = -t353 * t296 + t322 * t356;
t344 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t372;
t286 = -qJDD(3) * pkin(3) - pkin(8) * t358 + t338 * t373 - t291;
t284 = -pkin(4) * t311 - qJ(5) * t331 + t319 * t336 + qJDD(5) + t286;
t364 = -m(6) * t284 + t311 * mrSges(6,1) + t335 * t317;
t360 = -m(5) * t286 + t311 * mrSges(5,1) + t383 * t312 + t335 * t318 + t374 * t336 + t364;
t275 = m(4) * t291 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t339 + qJD(3) * t344 - t337 * t373 + t360;
t366 = t356 * t269 - t275 * t353;
t271 = t272 * t355 + t274 * t352;
t361 = -m(4) * t295 + t340 * mrSges(4,1) - mrSges(4,2) * t339 - t343 * t373 + t344 * t372 - t271;
t327 = Ifges(4,5) * qJD(3) + (t353 * Ifges(4,1) + t356 * Ifges(4,4)) * qJD(2);
t326 = Ifges(4,6) * qJD(3) + (t353 * Ifges(4,4) + t356 * Ifges(4,2)) * qJD(2);
t277 = mrSges(6,2) * t312 + t320 * t336 - t364;
t270 = mrSges(5,2) * t286 + mrSges(6,2) * t284 - mrSges(5,3) * t282 - mrSges(6,3) * t279 - qJ(5) * t276 + t382 * t311 + t391 * t312 + t381 * t332 - t377 * t335 + t376 * t346;
t267 = -mrSges(5,1) * t286 + mrSges(5,3) * t283 - mrSges(6,1) * t284 + mrSges(6,3) * t281 - pkin(4) * t277 + qJ(5) * t369 + (-qJ(5) * t320 + t384) * t346 + t377 * t336 + (-qJ(5) * mrSges(6,2) + t380) * t332 + t382 * t312 + t390 * t311;
t1 = [m(2) * t347 + t351 * (m(3) * t322 + t269 * t353 + t275 * t356) + (t354 * (m(3) * t299 - mrSges(3,1) * t359 - qJDD(2) * mrSges(3,2) + t366) + t357 * (m(3) * t298 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t359 + t361)) * t349; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t298 - mrSges(3,2) * t299 + t353 * (mrSges(4,2) * t295 - mrSges(4,3) * t291 + Ifges(4,1) * t339 + Ifges(4,4) * t340 + Ifges(4,5) * qJDD(3) - pkin(8) * t271 - qJD(3) * t326 - t267 * t352 + t270 * t355) + t356 * (-mrSges(4,1) * t295 + mrSges(4,3) * t292 + Ifges(4,4) * t339 + Ifges(4,2) * t340 + Ifges(4,6) * qJDD(3) - pkin(3) * t271 + qJD(3) * t327 - t389) + pkin(2) * t361 + pkin(7) * t366; Ifges(4,5) * t339 + Ifges(4,6) * t340 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t291 - mrSges(4,2) * t292 + t352 * t270 + t355 * t267 + pkin(3) * t360 + pkin(8) * t365 + (t353 * t326 - t356 * t327) * qJD(2); t389; t277;];
tauJ = t1;
