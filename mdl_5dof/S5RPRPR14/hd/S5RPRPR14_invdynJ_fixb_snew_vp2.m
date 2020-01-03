% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:41
% EndTime: 2019-12-31 18:34:42
% DurationCPUTime: 0.99s
% Computational Cost: add. (6137->216), mult. (13150->277), div. (0->0), fcn. (8009->8), ass. (0->87)
t356 = sin(qJ(1));
t359 = cos(qJ(1));
t370 = -t359 * g(1) - t356 * g(2);
t365 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t370;
t352 = sin(pkin(8));
t353 = cos(pkin(8));
t355 = sin(qJ(3));
t358 = cos(qJ(3));
t334 = (t358 * t352 + t355 * t353) * qJD(1);
t380 = 2 * qJD(4);
t379 = -pkin(1) - pkin(6);
t361 = qJD(1) ^ 2;
t373 = t356 * g(1) - t359 * g(2);
t364 = -t361 * qJ(2) + qJDD(2) - t373;
t329 = t379 * qJDD(1) + t364;
t320 = t355 * g(3) + t358 * t329;
t376 = qJD(1) * qJD(3);
t374 = t355 * t376;
t342 = t358 * qJDD(1) - t374;
t303 = (-t342 - t374) * qJ(4) + (-t355 * t358 * t361 + qJDD(3)) * pkin(3) + t320;
t321 = -t358 * g(3) + t355 * t329;
t341 = -t355 * qJDD(1) - t358 * t376;
t377 = t358 * qJD(1);
t344 = qJD(3) * pkin(3) - qJ(4) * t377;
t351 = t355 ^ 2;
t304 = -t351 * t361 * pkin(3) + t341 * qJ(4) - qJD(3) * t344 + t321;
t293 = t352 * t303 + t353 * t304 - t334 * t380;
t378 = qJD(1) * t355;
t335 = -t352 * t378 + t353 * t377;
t314 = t334 * mrSges(5,1) + t335 * mrSges(5,2);
t318 = t353 * t341 - t352 * t342;
t328 = qJD(3) * mrSges(5,1) - t335 * mrSges(5,3);
t315 = t334 * pkin(4) - t335 * pkin(7);
t360 = qJD(3) ^ 2;
t290 = -t360 * pkin(4) + qJDD(3) * pkin(7) - t334 * t315 + t293;
t306 = -t341 * pkin(3) + qJDD(4) + t344 * t377 + (-qJ(4) * t351 + t379) * t361 + t365;
t319 = t352 * t341 + t353 * t342;
t291 = (qJD(3) * t334 - t319) * pkin(7) + (qJD(3) * t335 - t318) * pkin(4) + t306;
t354 = sin(qJ(5));
t357 = cos(qJ(5));
t287 = -t354 * t290 + t357 * t291;
t322 = t357 * qJD(3) - t354 * t335;
t300 = t322 * qJD(5) + t354 * qJDD(3) + t357 * t319;
t323 = t354 * qJD(3) + t357 * t335;
t307 = -t322 * mrSges(6,1) + t323 * mrSges(6,2);
t332 = qJD(5) + t334;
t308 = -t332 * mrSges(6,2) + t322 * mrSges(6,3);
t317 = qJDD(5) - t318;
t285 = m(6) * t287 + t317 * mrSges(6,1) - t300 * mrSges(6,3) - t323 * t307 + t332 * t308;
t288 = t357 * t290 + t354 * t291;
t299 = -t323 * qJD(5) + t357 * qJDD(3) - t354 * t319;
t309 = t332 * mrSges(6,1) - t323 * mrSges(6,3);
t286 = m(6) * t288 - t317 * mrSges(6,2) + t299 * mrSges(6,3) + t322 * t307 - t332 * t309;
t371 = -t354 * t285 + t357 * t286;
t275 = m(5) * t293 - qJDD(3) * mrSges(5,2) + t318 * mrSges(5,3) - qJD(3) * t328 - t334 * t314 + t371;
t368 = -t353 * t303 + t352 * t304;
t292 = -0.2e1 * qJD(4) * t335 - t368;
t327 = -qJD(3) * mrSges(5,2) - t334 * mrSges(5,3);
t289 = -qJDD(3) * pkin(4) - t360 * pkin(7) + (t380 + t315) * t335 + t368;
t363 = -m(6) * t289 + t299 * mrSges(6,1) - t300 * mrSges(6,2) + t322 * t308 - t323 * t309;
t281 = m(5) * t292 + qJDD(3) * mrSges(5,1) - t319 * mrSges(5,3) + qJD(3) * t327 - t335 * t314 + t363;
t272 = t352 * t275 + t353 * t281;
t277 = t357 * t285 + t354 * t286;
t372 = t353 * t275 - t352 * t281;
t340 = (t355 * mrSges(4,1) + t358 * mrSges(4,2)) * qJD(1);
t343 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t378;
t345 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t377;
t369 = t358 * (m(4) * t320 + qJDD(3) * mrSges(4,1) - t342 * mrSges(4,3) + qJD(3) * t343 - t340 * t377 + t272) + t355 * (m(4) * t321 - qJDD(3) * mrSges(4,2) + t341 * mrSges(4,3) - qJD(3) * t345 - t340 * t378 + t372);
t276 = m(5) * t306 - t318 * mrSges(5,1) + t319 * mrSges(5,2) + t334 * t327 + t335 * t328 + t277;
t295 = Ifges(6,4) * t323 + Ifges(6,2) * t322 + Ifges(6,6) * t332;
t296 = Ifges(6,1) * t323 + Ifges(6,4) * t322 + Ifges(6,5) * t332;
t362 = mrSges(6,1) * t287 - mrSges(6,2) * t288 + Ifges(6,5) * t300 + Ifges(6,6) * t299 + Ifges(6,3) * t317 + t323 * t295 - t322 * t296;
t338 = (Ifges(4,5) * qJD(3)) + (t358 * Ifges(4,1) - t355 * Ifges(4,4)) * qJD(1);
t337 = (Ifges(4,6) * qJD(3)) + (t358 * Ifges(4,4) - t355 * Ifges(4,2)) * qJD(1);
t333 = -qJDD(1) * pkin(1) + t364;
t330 = t361 * pkin(1) - t365;
t326 = t379 * t361 + t365;
t312 = Ifges(5,1) * t335 - Ifges(5,4) * t334 + (Ifges(5,5) * qJD(3));
t311 = Ifges(5,4) * t335 - Ifges(5,2) * t334 + (Ifges(5,6) * qJD(3));
t310 = Ifges(5,5) * t335 - Ifges(5,6) * t334 + (Ifges(5,3) * qJD(3));
t294 = Ifges(6,5) * t323 + Ifges(6,6) * t322 + Ifges(6,3) * t332;
t279 = mrSges(6,2) * t289 - mrSges(6,3) * t287 + Ifges(6,1) * t300 + Ifges(6,4) * t299 + Ifges(6,5) * t317 + t322 * t294 - t332 * t295;
t278 = -mrSges(6,1) * t289 + mrSges(6,3) * t288 + Ifges(6,4) * t300 + Ifges(6,2) * t299 + Ifges(6,6) * t317 - t323 * t294 + t332 * t296;
t269 = -mrSges(5,1) * t306 + mrSges(5,3) * t293 + Ifges(5,4) * t319 + Ifges(5,2) * t318 + Ifges(5,6) * qJDD(3) - pkin(4) * t277 + qJD(3) * t312 - t335 * t310 - t362;
t268 = mrSges(5,2) * t306 - mrSges(5,3) * t292 + Ifges(5,1) * t319 + Ifges(5,4) * t318 + Ifges(5,5) * qJDD(3) - pkin(7) * t277 - qJD(3) * t311 - t354 * t278 + t357 * t279 - t334 * t310;
t267 = m(3) * t333 + qJDD(1) * mrSges(3,2) - (t361 * mrSges(3,3)) + t369;
t1 = [mrSges(2,1) * t373 - mrSges(2,2) * t370 + mrSges(3,2) * t333 - mrSges(3,3) * t330 + t358 * (mrSges(4,2) * t326 - mrSges(4,3) * t320 + Ifges(4,1) * t342 + Ifges(4,4) * t341 + Ifges(4,5) * qJDD(3) - qJ(4) * t272 - qJD(3) * t337 + t353 * t268 - t352 * t269) - t355 * (-mrSges(4,1) * t326 + mrSges(4,3) * t321 + Ifges(4,4) * t342 + Ifges(4,2) * t341 + Ifges(4,6) * qJDD(3) - pkin(3) * t276 + qJ(4) * t372 + qJD(3) * t338 + t352 * t268 + t353 * t269) - pkin(6) * t369 - pkin(1) * t267 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t330 + m(4) * t326 - t341 * mrSges(4,1) + t361 * mrSges(3,2) + t342 * mrSges(4,2) + t276 + qJDD(1) * mrSges(3,3) + (t343 * t355 + t345 * t358) * qJD(1)) * qJ(2); t267; Ifges(4,5) * t342 + Ifges(4,6) * t341 + mrSges(4,1) * t320 - mrSges(4,2) * t321 + Ifges(5,5) * t319 + Ifges(5,6) * t318 + t335 * t311 + t334 * t312 + mrSges(5,1) * t292 - mrSges(5,2) * t293 + t354 * t279 + t357 * t278 + pkin(4) * t363 + pkin(7) * t371 + pkin(3) * t272 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t358 * t337 + t355 * t338) * qJD(1); t276; t362;];
tauJ = t1;
