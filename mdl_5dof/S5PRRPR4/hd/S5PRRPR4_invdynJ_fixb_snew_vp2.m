% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:54
% EndTime: 2019-12-05 16:21:56
% DurationCPUTime: 0.96s
% Computational Cost: add. (7216->211), mult. (15899->277), div. (0->0), fcn. (10650->10), ass. (0->86)
t351 = sin(pkin(8));
t353 = cos(pkin(8));
t340 = -t353 * g(1) - t351 * g(2);
t349 = -g(3) + qJDD(1);
t356 = sin(qJ(2));
t359 = cos(qJ(2));
t320 = t359 * t340 + t356 * t349;
t360 = qJD(2) ^ 2;
t315 = -t360 * pkin(2) + qJDD(2) * pkin(6) + t320;
t339 = -t351 * g(1) + t353 * g(2);
t355 = sin(qJ(3));
t358 = cos(qJ(3));
t301 = -t355 * t315 + t358 * t339;
t369 = qJD(2) * qJD(3);
t368 = t358 * t369;
t337 = t355 * qJDD(2) + t368;
t296 = (-t337 + t368) * qJ(4) + (t355 * t358 * t360 + qJDD(3)) * pkin(3) + t301;
t302 = t358 * t315 + t355 * t339;
t338 = t358 * qJDD(2) - t355 * t369;
t370 = t355 * qJD(2);
t341 = qJD(3) * pkin(3) - qJ(4) * t370;
t348 = t358 ^ 2;
t297 = -t348 * t360 * pkin(3) + t338 * qJ(4) - qJD(3) * t341 + t302;
t350 = sin(pkin(9));
t352 = cos(pkin(9));
t325 = (t358 * t350 + t355 * t352) * qJD(2);
t276 = -0.2e1 * qJD(4) * t325 + t352 * t296 - t350 * t297;
t312 = t352 * t337 + t350 * t338;
t324 = (-t355 * t350 + t358 * t352) * qJD(2);
t274 = (qJD(3) * t324 - t312) * pkin(7) + (t324 * t325 + qJDD(3)) * pkin(4) + t276;
t277 = 0.2e1 * qJD(4) * t324 + t350 * t296 + t352 * t297;
t311 = -t350 * t337 + t352 * t338;
t318 = qJD(3) * pkin(4) - t325 * pkin(7);
t323 = t324 ^ 2;
t275 = -t323 * pkin(4) + t311 * pkin(7) - qJD(3) * t318 + t277;
t354 = sin(qJ(5));
t357 = cos(qJ(5));
t272 = t357 * t274 - t354 * t275;
t306 = t357 * t324 - t354 * t325;
t286 = t306 * qJD(5) + t354 * t311 + t357 * t312;
t307 = t354 * t324 + t357 * t325;
t292 = -t306 * mrSges(6,1) + t307 * mrSges(6,2);
t347 = qJD(3) + qJD(5);
t299 = -t347 * mrSges(6,2) + t306 * mrSges(6,3);
t346 = qJDD(3) + qJDD(5);
t268 = m(6) * t272 + t346 * mrSges(6,1) - t286 * mrSges(6,3) - t307 * t292 + t347 * t299;
t273 = t354 * t274 + t357 * t275;
t285 = -t307 * qJD(5) + t357 * t311 - t354 * t312;
t300 = t347 * mrSges(6,1) - t307 * mrSges(6,3);
t269 = m(6) * t273 - t346 * mrSges(6,2) + t285 * mrSges(6,3) + t306 * t292 - t347 * t300;
t262 = t357 * t268 + t354 * t269;
t309 = -t324 * mrSges(5,1) + t325 * mrSges(5,2);
t316 = -qJD(3) * mrSges(5,2) + t324 * mrSges(5,3);
t260 = m(5) * t276 + qJDD(3) * mrSges(5,1) - t312 * mrSges(5,3) + qJD(3) * t316 - t325 * t309 + t262;
t317 = qJD(3) * mrSges(5,1) - t325 * mrSges(5,3);
t365 = -t354 * t268 + t357 * t269;
t261 = m(5) * t277 - qJDD(3) * mrSges(5,2) + t311 * mrSges(5,3) - qJD(3) * t317 + t324 * t309 + t365;
t256 = t352 * t260 + t350 * t261;
t371 = qJD(2) * t358;
t336 = (-t358 * mrSges(4,1) + t355 * mrSges(4,2)) * qJD(2);
t342 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t370;
t343 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t371;
t366 = -t350 * t260 + t352 * t261;
t367 = -t355 * (m(4) * t301 + qJDD(3) * mrSges(4,1) - t337 * mrSges(4,3) + qJD(3) * t343 - t336 * t370 + t256) + t358 * (m(4) * t302 - qJDD(3) * mrSges(4,2) + t338 * mrSges(4,3) - qJD(3) * t342 + t336 * t371 + t366);
t319 = -t356 * t340 + t359 * t349;
t363 = -qJDD(2) * pkin(2) - t319;
t298 = -t338 * pkin(3) + qJDD(4) + t341 * t370 + (-qJ(4) * t348 - pkin(6)) * t360 + t363;
t279 = -t311 * pkin(4) - t323 * pkin(7) + t325 * t318 + t298;
t364 = m(6) * t279 - t285 * mrSges(6,1) + t286 * mrSges(6,2) - t306 * t299 + t307 * t300;
t288 = Ifges(6,4) * t307 + Ifges(6,2) * t306 + Ifges(6,6) * t347;
t289 = Ifges(6,1) * t307 + Ifges(6,4) * t306 + Ifges(6,5) * t347;
t362 = mrSges(6,1) * t272 - mrSges(6,2) * t273 + Ifges(6,5) * t286 + Ifges(6,6) * t285 + Ifges(6,3) * t346 + t307 * t288 - t306 * t289;
t270 = m(5) * t298 - t311 * mrSges(5,1) + t312 * mrSges(5,2) - t324 * t316 + t325 * t317 + t364;
t314 = -t360 * pkin(6) + t363;
t361 = -m(4) * t314 + t338 * mrSges(4,1) - t337 * mrSges(4,2) - t342 * t370 + t343 * t371 - t270;
t328 = Ifges(4,5) * qJD(3) + (t355 * Ifges(4,1) + t358 * Ifges(4,4)) * qJD(2);
t327 = Ifges(4,6) * qJD(3) + (t355 * Ifges(4,4) + t358 * Ifges(4,2)) * qJD(2);
t305 = Ifges(5,1) * t325 + Ifges(5,4) * t324 + Ifges(5,5) * qJD(3);
t304 = Ifges(5,4) * t325 + Ifges(5,2) * t324 + Ifges(5,6) * qJD(3);
t303 = Ifges(5,5) * t325 + Ifges(5,6) * t324 + Ifges(5,3) * qJD(3);
t287 = Ifges(6,5) * t307 + Ifges(6,6) * t306 + Ifges(6,3) * t347;
t264 = mrSges(6,2) * t279 - mrSges(6,3) * t272 + Ifges(6,1) * t286 + Ifges(6,4) * t285 + Ifges(6,5) * t346 + t306 * t287 - t347 * t288;
t263 = -mrSges(6,1) * t279 + mrSges(6,3) * t273 + Ifges(6,4) * t286 + Ifges(6,2) * t285 + Ifges(6,6) * t346 - t307 * t287 + t347 * t289;
t253 = mrSges(5,2) * t298 - mrSges(5,3) * t276 + Ifges(5,1) * t312 + Ifges(5,4) * t311 + Ifges(5,5) * qJDD(3) - pkin(7) * t262 - qJD(3) * t304 - t354 * t263 + t357 * t264 + t324 * t303;
t252 = -mrSges(5,1) * t298 + mrSges(5,3) * t277 + Ifges(5,4) * t312 + Ifges(5,2) * t311 + Ifges(5,6) * qJDD(3) - pkin(4) * t364 + pkin(7) * t365 + qJD(3) * t305 + t357 * t263 + t354 * t264 - t325 * t303;
t1 = [m(2) * t349 + t356 * (m(3) * t320 - t360 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t367) + t359 * (m(3) * t319 + qJDD(2) * mrSges(3,1) - t360 * mrSges(3,2) + t361); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t319 - mrSges(3,2) * t320 + t355 * (mrSges(4,2) * t314 - mrSges(4,3) * t301 + Ifges(4,1) * t337 + Ifges(4,4) * t338 + Ifges(4,5) * qJDD(3) - qJ(4) * t256 - qJD(3) * t327 - t350 * t252 + t352 * t253) + t358 * (-mrSges(4,1) * t314 + mrSges(4,3) * t302 + Ifges(4,4) * t337 + Ifges(4,2) * t338 + Ifges(4,6) * qJDD(3) - pkin(3) * t270 + qJ(4) * t366 + qJD(3) * t328 + t352 * t252 + t350 * t253) + pkin(2) * t361 + pkin(6) * t367; t362 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t355 * t327 - t358 * t328) * qJD(2) + Ifges(4,5) * t337 + Ifges(4,6) * t338 - t324 * t305 + t325 * t304 + Ifges(5,6) * t311 + Ifges(5,5) * t312 + mrSges(4,1) * t301 - mrSges(4,2) * t302 - mrSges(5,2) * t277 + mrSges(5,1) * t276 + pkin(3) * t256 + pkin(4) * t262; t270; t362;];
tauJ = t1;
