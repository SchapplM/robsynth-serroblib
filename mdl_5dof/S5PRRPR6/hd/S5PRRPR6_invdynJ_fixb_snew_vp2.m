% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:42
% EndTime: 2019-12-05 16:30:44
% DurationCPUTime: 1.26s
% Computational Cost: add. (8364->215), mult. (17093->281), div. (0->0), fcn. (11745->12), ass. (0->92)
t354 = sin(pkin(9));
t357 = cos(pkin(9));
t345 = t354 * g(1) - t357 * g(2);
t352 = -g(3) + qJDD(1);
t355 = sin(pkin(5));
t358 = cos(pkin(5));
t380 = t345 * t358 + t352 * t355;
t346 = -t357 * g(1) - t354 * g(2);
t361 = sin(qJ(2));
t364 = cos(qJ(2));
t309 = -t361 * t346 + t380 * t364;
t310 = t364 * t346 + t380 * t361;
t366 = qJD(2) ^ 2;
t306 = -t366 * pkin(2) + qJDD(2) * pkin(7) + t310;
t325 = -t355 * t345 + t358 * t352;
t360 = sin(qJ(3));
t363 = cos(qJ(3));
t301 = t363 * t306 + t360 * t325;
t341 = (-t363 * pkin(3) - t360 * qJ(4)) * qJD(2);
t365 = qJD(3) ^ 2;
t376 = t363 * qJD(2);
t289 = -t365 * pkin(3) + qJDD(3) * qJ(4) + t341 * t376 + t301;
t305 = -qJDD(2) * pkin(2) - t366 * pkin(7) - t309;
t375 = qJD(2) * qJD(3);
t374 = t363 * t375;
t343 = t360 * qJDD(2) + t374;
t351 = t360 * t375;
t344 = t363 * qJDD(2) - t351;
t293 = (-t343 - t374) * qJ(4) + (-t344 + t351) * pkin(3) + t305;
t353 = sin(pkin(10));
t356 = cos(pkin(10));
t377 = t360 * qJD(2);
t336 = t353 * qJD(3) + t356 * t377;
t284 = -0.2e1 * qJD(4) * t336 - t353 * t289 + t356 * t293;
t323 = t353 * qJDD(3) + t356 * t343;
t335 = t356 * qJD(3) - t353 * t377;
t282 = (-t335 * t376 - t323) * pkin(8) + (t335 * t336 - t344) * pkin(4) + t284;
t285 = 0.2e1 * qJD(4) * t335 + t356 * t289 + t353 * t293;
t322 = t356 * qJDD(3) - t353 * t343;
t324 = -pkin(4) * t376 - t336 * pkin(8);
t334 = t335 ^ 2;
t283 = -t334 * pkin(4) + t322 * pkin(8) + t324 * t376 + t285;
t359 = sin(qJ(5));
t362 = cos(qJ(5));
t280 = t362 * t282 - t359 * t283;
t315 = t362 * t335 - t359 * t336;
t295 = t315 * qJD(5) + t359 * t322 + t362 * t323;
t316 = t359 * t335 + t362 * t336;
t302 = -t315 * mrSges(6,1) + t316 * mrSges(6,2);
t350 = qJD(5) - t376;
t307 = -t350 * mrSges(6,2) + t315 * mrSges(6,3);
t338 = qJDD(5) - t344;
t277 = m(6) * t280 + t338 * mrSges(6,1) - t295 * mrSges(6,3) - t316 * t302 + t350 * t307;
t281 = t359 * t282 + t362 * t283;
t294 = -t316 * qJD(5) + t362 * t322 - t359 * t323;
t308 = t350 * mrSges(6,1) - t316 * mrSges(6,3);
t278 = m(6) * t281 - t338 * mrSges(6,2) + t294 * mrSges(6,3) + t315 * t302 - t350 * t308;
t270 = t362 * t277 + t359 * t278;
t342 = (-t363 * mrSges(4,1) + t360 * mrSges(4,2)) * qJD(2);
t347 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t377;
t317 = -t335 * mrSges(5,1) + t336 * mrSges(5,2);
t320 = mrSges(5,2) * t376 + t335 * mrSges(5,3);
t268 = m(5) * t284 - t344 * mrSges(5,1) - t323 * mrSges(5,3) - t336 * t317 - t320 * t376 + t270;
t321 = -mrSges(5,1) * t376 - t336 * mrSges(5,3);
t371 = -t359 * t277 + t362 * t278;
t269 = m(5) * t285 + t344 * mrSges(5,2) + t322 * mrSges(5,3) + t335 * t317 + t321 * t376 + t371;
t372 = -t353 * t268 + t356 * t269;
t265 = m(4) * t301 - qJDD(3) * mrSges(4,2) + t344 * mrSges(4,3) - qJD(3) * t347 + t342 * t376 + t372;
t300 = -t360 * t306 + t363 * t325;
t288 = -qJDD(3) * pkin(3) - t365 * qJ(4) + t341 * t377 + qJDD(4) - t300;
t286 = -t322 * pkin(4) - t334 * pkin(8) + t336 * t324 + t288;
t369 = m(6) * t286 - t294 * mrSges(6,1) + t295 * mrSges(6,2) - t315 * t307 + t316 * t308;
t279 = m(5) * t288 - t322 * mrSges(5,1) + t323 * mrSges(5,2) - t335 * t320 + t336 * t321 + t369;
t348 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t376;
t273 = m(4) * t300 + qJDD(3) * mrSges(4,1) - t343 * mrSges(4,3) + qJD(3) * t348 - t342 * t377 - t279;
t373 = t363 * t265 - t360 * t273;
t266 = t356 * t268 + t353 * t269;
t368 = -m(4) * t305 + t344 * mrSges(4,1) - t343 * mrSges(4,2) - t347 * t377 + t348 * t376 - t266;
t297 = Ifges(6,4) * t316 + Ifges(6,2) * t315 + Ifges(6,6) * t350;
t298 = Ifges(6,1) * t316 + Ifges(6,4) * t315 + Ifges(6,5) * t350;
t367 = mrSges(6,1) * t280 - mrSges(6,2) * t281 + Ifges(6,5) * t295 + Ifges(6,6) * t294 + Ifges(6,3) * t338 + t316 * t297 - t315 * t298;
t332 = Ifges(4,5) * qJD(3) + (t360 * Ifges(4,1) + t363 * Ifges(4,4)) * qJD(2);
t331 = Ifges(4,6) * qJD(3) + (t360 * Ifges(4,4) + Ifges(4,2) * t363) * qJD(2);
t313 = Ifges(5,1) * t336 + Ifges(5,4) * t335 - Ifges(5,5) * t376;
t312 = Ifges(5,4) * t336 + Ifges(5,2) * t335 - Ifges(5,6) * t376;
t311 = Ifges(5,5) * t336 + Ifges(5,6) * t335 - Ifges(5,3) * t376;
t296 = Ifges(6,5) * t316 + Ifges(6,6) * t315 + Ifges(6,3) * t350;
t272 = mrSges(6,2) * t286 - mrSges(6,3) * t280 + Ifges(6,1) * t295 + Ifges(6,4) * t294 + Ifges(6,5) * t338 + t315 * t296 - t350 * t297;
t271 = -mrSges(6,1) * t286 + mrSges(6,3) * t281 + Ifges(6,4) * t295 + Ifges(6,2) * t294 + Ifges(6,6) * t338 - t316 * t296 + t350 * t298;
t263 = mrSges(5,2) * t288 - mrSges(5,3) * t284 + Ifges(5,1) * t323 + Ifges(5,4) * t322 - Ifges(5,5) * t344 - pkin(8) * t270 - t359 * t271 + t362 * t272 + t335 * t311 + t312 * t376;
t262 = -mrSges(5,1) * t288 + mrSges(5,3) * t285 + Ifges(5,4) * t323 + Ifges(5,2) * t322 - Ifges(5,6) * t344 - pkin(4) * t369 + pkin(8) * t371 + t362 * t271 + t359 * t272 - t336 * t311 - t313 * t376;
t1 = [m(2) * t352 + t358 * (m(3) * t325 + t360 * t265 + t363 * t273) + (t361 * (m(3) * t310 - t366 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t373) + t364 * (m(3) * t309 + qJDD(2) * mrSges(3,1) - t366 * mrSges(3,2) + t368)) * t355; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t309 - mrSges(3,2) * t310 + t360 * (mrSges(4,2) * t305 - mrSges(4,3) * t300 + Ifges(4,1) * t343 + Ifges(4,4) * t344 + Ifges(4,5) * qJDD(3) - qJ(4) * t266 - qJD(3) * t331 - t353 * t262 + t356 * t263) + t363 * (-mrSges(4,1) * t305 - mrSges(5,1) * t284 + mrSges(5,2) * t285 + mrSges(4,3) * t301 + Ifges(4,4) * t343 - Ifges(5,5) * t323 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t322 - pkin(3) * t266 - pkin(4) * t270 + qJD(3) * t332 - t336 * t312 + t335 * t313 - t367 + (Ifges(4,2) + Ifges(5,3)) * t344) + pkin(2) * t368 + pkin(7) * t373; Ifges(4,5) * t343 + Ifges(4,6) * t344 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t300 - mrSges(4,2) * t301 + t353 * t263 + t356 * t262 - pkin(3) * t279 + qJ(4) * t372 + (t360 * t331 - t363 * t332) * qJD(2); t279; t367;];
tauJ = t1;
