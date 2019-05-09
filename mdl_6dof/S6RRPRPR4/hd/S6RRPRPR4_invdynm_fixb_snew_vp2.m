% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 13:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:30:51
% EndTime: 2019-05-06 13:32:44
% DurationCPUTime: 70.96s
% Computational Cost: add. (1124220->396), mult. (2986482->519), div. (0->0), fcn. (2385921->14), ass. (0->165)
t391 = -2 * qJD(3);
t347 = sin(pkin(11));
t350 = cos(pkin(11));
t354 = sin(qJ(2));
t358 = cos(qJ(2));
t348 = sin(pkin(6));
t383 = qJD(1) * t348;
t324 = (t347 * t354 - t350 * t358) * t383;
t381 = qJD(1) * qJD(2);
t333 = (qJDD(1) * t354 + t358 * t381) * t348;
t351 = cos(pkin(6));
t340 = qJDD(1) * t351 + qJDD(2);
t341 = qJD(1) * t351 + qJD(2);
t355 = sin(qJ(1));
t359 = cos(qJ(1));
t337 = t355 * g(1) - g(2) * t359;
t360 = qJD(1) ^ 2;
t389 = pkin(8) * t348;
t330 = qJDD(1) * pkin(1) + t360 * t389 + t337;
t338 = -g(1) * t359 - g(2) * t355;
t331 = -pkin(1) * t360 + qJDD(1) * t389 + t338;
t384 = t351 * t358;
t372 = t330 * t384 - t354 * t331;
t388 = t348 ^ 2 * t360;
t269 = t340 * pkin(2) - t333 * qJ(3) + (pkin(2) * t354 * t388 + (qJ(3) * qJD(1) * t341 - g(3)) * t348) * t358 + t372;
t385 = t351 * t354;
t387 = t348 * t354;
t298 = -g(3) * t387 + t330 * t385 + t358 * t331;
t379 = t354 * t383;
t327 = pkin(2) * t341 - qJ(3) * t379;
t334 = (qJDD(1) * t358 - t354 * t381) * t348;
t380 = t358 ^ 2 * t388;
t272 = -pkin(2) * t380 + qJ(3) * t334 - t327 * t341 + t298;
t325 = (t347 * t358 + t350 * t354) * t383;
t244 = t350 * t269 - t347 * t272 + t325 * t391;
t390 = 2 * qJD(5);
t386 = t348 * t358;
t245 = t347 * t269 + t350 * t272 + t324 * t391;
t299 = mrSges(4,1) * t324 + mrSges(4,2) * t325;
t306 = -t333 * t347 + t334 * t350;
t312 = mrSges(4,1) * t341 - mrSges(4,3) * t325;
t300 = pkin(3) * t324 - pkin(9) * t325;
t339 = t341 ^ 2;
t237 = -pkin(3) * t339 + pkin(9) * t340 - t300 * t324 + t245;
t316 = -t351 * g(3) - t348 * t330;
t283 = -t334 * pkin(2) - qJ(3) * t380 + t327 * t379 + qJDD(3) + t316;
t307 = t333 * t350 + t334 * t347;
t248 = (t324 * t341 - t307) * pkin(9) + (t325 * t341 - t306) * pkin(3) + t283;
t353 = sin(qJ(4));
t357 = cos(qJ(4));
t229 = -t353 * t237 + t357 * t248;
t309 = -t325 * t353 + t341 * t357;
t280 = qJD(4) * t309 + t307 * t357 + t340 * t353;
t305 = qJDD(4) - t306;
t310 = t325 * t357 + t341 * t353;
t323 = qJD(4) + t324;
t226 = (t309 * t323 - t280) * qJ(5) + (t309 * t310 + t305) * pkin(4) + t229;
t230 = t357 * t237 + t353 * t248;
t279 = -qJD(4) * t310 - t307 * t353 + t340 * t357;
t291 = pkin(4) * t323 - qJ(5) * t310;
t308 = t309 ^ 2;
t228 = -pkin(4) * t308 + qJ(5) * t279 - t291 * t323 + t230;
t346 = sin(pkin(12));
t349 = cos(pkin(12));
t285 = t309 * t349 - t310 * t346;
t223 = t346 * t226 + t349 * t228 + t285 * t390;
t254 = t279 * t349 - t280 * t346;
t286 = t309 * t346 + t310 * t349;
t262 = -mrSges(6,1) * t285 + mrSges(6,2) * t286;
t271 = mrSges(6,1) * t323 - mrSges(6,3) * t286;
t263 = -pkin(5) * t285 - pkin(10) * t286;
t322 = t323 ^ 2;
t220 = -pkin(5) * t322 + pkin(10) * t305 + t263 * t285 + t223;
t236 = -t340 * pkin(3) - t339 * pkin(9) + t325 * t300 - t244;
t231 = -t279 * pkin(4) - t308 * qJ(5) + t310 * t291 + qJDD(5) + t236;
t255 = t279 * t346 + t280 * t349;
t224 = (-t285 * t323 - t255) * pkin(10) + (t286 * t323 - t254) * pkin(5) + t231;
t352 = sin(qJ(6));
t356 = cos(qJ(6));
t217 = -t220 * t352 + t224 * t356;
t265 = -t286 * t352 + t323 * t356;
t234 = qJD(6) * t265 + t255 * t356 + t305 * t352;
t266 = t286 * t356 + t323 * t352;
t249 = -mrSges(7,1) * t265 + mrSges(7,2) * t266;
t253 = qJDD(6) - t254;
t284 = qJD(6) - t285;
t256 = -mrSges(7,2) * t284 + mrSges(7,3) * t265;
t213 = m(7) * t217 + mrSges(7,1) * t253 - mrSges(7,3) * t234 - t249 * t266 + t256 * t284;
t218 = t220 * t356 + t224 * t352;
t233 = -qJD(6) * t266 - t255 * t352 + t305 * t356;
t257 = mrSges(7,1) * t284 - mrSges(7,3) * t266;
t214 = m(7) * t218 - mrSges(7,2) * t253 + mrSges(7,3) * t233 + t249 * t265 - t257 * t284;
t374 = -t213 * t352 + t356 * t214;
t199 = m(6) * t223 - mrSges(6,2) * t305 + mrSges(6,3) * t254 + t262 * t285 - t271 * t323 + t374;
t371 = -t349 * t226 + t346 * t228;
t222 = -0.2e1 * qJD(5) * t286 - t371;
t270 = -mrSges(6,2) * t323 + mrSges(6,3) * t285;
t219 = -t305 * pkin(5) - t322 * pkin(10) + (t390 + t263) * t286 + t371;
t369 = -m(7) * t219 + t233 * mrSges(7,1) - mrSges(7,2) * t234 + t265 * t256 - t257 * t266;
t209 = m(6) * t222 + mrSges(6,1) * t305 - mrSges(6,3) * t255 - t262 * t286 + t270 * t323 + t369;
t194 = t346 * t199 + t349 * t209;
t288 = -mrSges(5,1) * t309 + mrSges(5,2) * t310;
t290 = -mrSges(5,2) * t323 + mrSges(5,3) * t309;
t192 = m(5) * t229 + mrSges(5,1) * t305 - mrSges(5,3) * t280 - t288 * t310 + t290 * t323 + t194;
t292 = mrSges(5,1) * t323 - mrSges(5,3) * t310;
t375 = t349 * t199 - t209 * t346;
t193 = m(5) * t230 - mrSges(5,2) * t305 + mrSges(5,3) * t279 + t288 * t309 - t292 * t323 + t375;
t376 = -t192 * t353 + t357 * t193;
t183 = m(4) * t245 - mrSges(4,2) * t340 + mrSges(4,3) * t306 - t299 * t324 - t312 * t341 + t376;
t311 = -mrSges(4,2) * t341 - mrSges(4,3) * t324;
t202 = t356 * t213 + t352 * t214;
t366 = m(6) * t231 - t254 * mrSges(6,1) + mrSges(6,2) * t255 - t285 * t270 + t271 * t286 + t202;
t362 = -m(5) * t236 + t279 * mrSges(5,1) - mrSges(5,2) * t280 + t309 * t290 - t292 * t310 - t366;
t196 = m(4) * t244 + mrSges(4,1) * t340 - mrSges(4,3) * t307 - t299 * t325 + t311 * t341 + t362;
t180 = t347 * t183 + t350 * t196;
t186 = t357 * t192 + t353 * t193;
t378 = t358 * t383;
t297 = -g(3) * t386 + t372;
t329 = -mrSges(3,2) * t341 + mrSges(3,3) * t378;
t332 = (-mrSges(3,1) * t358 + mrSges(3,2) * t354) * t383;
t178 = m(3) * t297 + mrSges(3,1) * t340 - mrSges(3,3) * t333 + t329 * t341 - t332 * t379 + t180;
t328 = mrSges(3,1) * t341 - mrSges(3,3) * t379;
t377 = t350 * t183 - t196 * t347;
t179 = m(3) * t298 - mrSges(3,2) * t340 + mrSges(3,3) * t334 - t328 * t341 + t332 * t378 + t377;
t169 = -t178 * t354 + t358 * t179;
t367 = m(4) * t283 - t306 * mrSges(4,1) + t307 * mrSges(4,2) + t324 * t311 + t325 * t312 + t186;
t184 = m(3) * t316 - t334 * mrSges(3,1) + t333 * mrSges(3,2) + (t328 * t354 - t329 * t358) * t383 + t367;
t166 = t178 * t384 + t179 * t385 - t184 * t348;
t238 = Ifges(7,5) * t266 + Ifges(7,6) * t265 + Ifges(7,3) * t284;
t240 = Ifges(7,1) * t266 + Ifges(7,4) * t265 + Ifges(7,5) * t284;
t206 = -mrSges(7,1) * t219 + mrSges(7,3) * t218 + Ifges(7,4) * t234 + Ifges(7,2) * t233 + Ifges(7,6) * t253 - t238 * t266 + t240 * t284;
t239 = Ifges(7,4) * t266 + Ifges(7,2) * t265 + Ifges(7,6) * t284;
t207 = mrSges(7,2) * t219 - mrSges(7,3) * t217 + Ifges(7,1) * t234 + Ifges(7,4) * t233 + Ifges(7,5) * t253 + t238 * t265 - t239 * t284;
t258 = Ifges(6,5) * t286 + Ifges(6,6) * t285 + Ifges(6,3) * t323;
t259 = Ifges(6,4) * t286 + Ifges(6,2) * t285 + Ifges(6,6) * t323;
t187 = mrSges(6,2) * t231 - mrSges(6,3) * t222 + Ifges(6,1) * t255 + Ifges(6,4) * t254 + Ifges(6,5) * t305 - pkin(10) * t202 - t206 * t352 + t207 * t356 + t258 * t285 - t259 * t323;
t260 = Ifges(6,1) * t286 + Ifges(6,4) * t285 + Ifges(6,5) * t323;
t363 = mrSges(7,1) * t217 - mrSges(7,2) * t218 + Ifges(7,5) * t234 + Ifges(7,6) * t233 + Ifges(7,3) * t253 + t239 * t266 - t240 * t265;
t188 = -mrSges(6,1) * t231 + mrSges(6,3) * t223 + Ifges(6,4) * t255 + Ifges(6,2) * t254 + Ifges(6,6) * t305 - pkin(5) * t202 - t258 * t286 + t260 * t323 - t363;
t273 = Ifges(5,5) * t310 + Ifges(5,6) * t309 + Ifges(5,3) * t323;
t275 = Ifges(5,1) * t310 + Ifges(5,4) * t309 + Ifges(5,5) * t323;
t172 = -mrSges(5,1) * t236 + mrSges(5,3) * t230 + Ifges(5,4) * t280 + Ifges(5,2) * t279 + Ifges(5,6) * t305 - pkin(4) * t366 + qJ(5) * t375 + t346 * t187 + t349 * t188 - t310 * t273 + t323 * t275;
t274 = Ifges(5,4) * t310 + Ifges(5,2) * t309 + Ifges(5,6) * t323;
t174 = mrSges(5,2) * t236 - mrSges(5,3) * t229 + Ifges(5,1) * t280 + Ifges(5,4) * t279 + Ifges(5,5) * t305 - qJ(5) * t194 + t187 * t349 - t188 * t346 + t273 * t309 - t274 * t323;
t293 = Ifges(4,5) * t325 - Ifges(4,6) * t324 + Ifges(4,3) * t341;
t294 = Ifges(4,4) * t325 - Ifges(4,2) * t324 + Ifges(4,6) * t341;
t162 = mrSges(4,2) * t283 - mrSges(4,3) * t244 + Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * t340 - pkin(9) * t186 - t172 * t353 + t174 * t357 - t293 * t324 - t294 * t341;
t295 = Ifges(4,1) * t325 - Ifges(4,4) * t324 + Ifges(4,5) * t341;
t364 = -mrSges(6,1) * t222 + mrSges(6,2) * t223 - Ifges(6,5) * t255 - Ifges(6,6) * t254 - Ifges(6,3) * t305 - pkin(5) * t369 - pkin(10) * t374 - t356 * t206 - t352 * t207 - t286 * t259 + t285 * t260;
t361 = mrSges(5,1) * t229 - mrSges(5,2) * t230 + Ifges(5,5) * t280 + Ifges(5,6) * t279 + Ifges(5,3) * t305 + pkin(4) * t194 + t310 * t274 - t309 * t275 - t364;
t170 = -mrSges(4,1) * t283 + mrSges(4,3) * t245 + Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * t340 - pkin(3) * t186 - t325 * t293 + t341 * t295 - t361;
t313 = Ifges(3,3) * t341 + (Ifges(3,5) * t354 + Ifges(3,6) * t358) * t383;
t315 = Ifges(3,5) * t341 + (Ifges(3,1) * t354 + Ifges(3,4) * t358) * t383;
t157 = -mrSges(3,1) * t316 + mrSges(3,3) * t298 + Ifges(3,4) * t333 + Ifges(3,2) * t334 + Ifges(3,6) * t340 - pkin(2) * t367 + qJ(3) * t377 + t347 * t162 + t350 * t170 - t313 * t379 + t341 * t315;
t314 = Ifges(3,6) * t341 + (Ifges(3,4) * t354 + Ifges(3,2) * t358) * t383;
t159 = mrSges(3,2) * t316 - mrSges(3,3) * t297 + Ifges(3,1) * t333 + Ifges(3,4) * t334 + Ifges(3,5) * t340 - qJ(3) * t180 + t162 * t350 - t170 * t347 + t313 * t378 - t314 * t341;
t365 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * t340 + pkin(3) * t362 + pkin(9) * t376 + t357 * t172 + t353 * t174 + t325 * t294 + t324 * t295;
t161 = pkin(2) * t180 + Ifges(3,3) * t340 + Ifges(3,5) * t333 + Ifges(3,6) * t334 + mrSges(3,1) * t297 - mrSges(3,2) * t298 + (t314 * t354 - t315 * t358) * t383 + t365;
t368 = mrSges(2,1) * t337 - mrSges(2,2) * t338 + Ifges(2,3) * qJDD(1) + pkin(1) * t166 + t157 * t386 + t159 * t387 + t351 * t161 + t169 * t389;
t167 = m(2) * t338 - mrSges(2,1) * t360 - qJDD(1) * mrSges(2,2) + t169;
t165 = t351 * t184 + (t178 * t358 + t179 * t354) * t348;
t163 = m(2) * t337 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t360 + t166;
t155 = -mrSges(2,2) * g(3) - mrSges(2,3) * t337 + Ifges(2,5) * qJDD(1) - t360 * Ifges(2,6) - t354 * t157 + t358 * t159 + (-t165 * t348 - t166 * t351) * pkin(8);
t154 = mrSges(2,1) * g(3) + mrSges(2,3) * t338 + t360 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t165 - t348 * t161 + (pkin(8) * t169 + t157 * t358 + t159 * t354) * t351;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t359 * t155 - t355 * t154 - pkin(7) * (t163 * t359 + t167 * t355), t155, t159, t162, t174, t187, t207; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t355 * t155 + t359 * t154 + pkin(7) * (-t163 * t355 + t167 * t359), t154, t157, t170, t172, t188, t206; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t368, t368, t161, t365, t361, -t364, t363;];
m_new  = t1;
