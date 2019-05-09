% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR9
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
% Datum: 2019-05-06 15:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:12:40
% EndTime: 2019-05-06 15:14:35
% DurationCPUTime: 78.42s
% Computational Cost: add. (1399310->399), mult. (3203209->519), div. (0->0), fcn. (2606059->14), ass. (0->163)
t381 = cos(qJ(4));
t344 = sin(pkin(6));
t380 = pkin(8) * t344;
t347 = cos(pkin(6));
t379 = t347 * g(3);
t350 = sin(qJ(2));
t378 = t344 * t350;
t353 = cos(qJ(2));
t377 = t344 * t353;
t376 = t347 * t350;
t375 = t347 * t353;
t373 = qJD(1) * t344;
t326 = (-pkin(2) * t353 - qJ(3) * t350) * t373;
t338 = qJD(1) * t347 + qJD(2);
t336 = t338 ^ 2;
t337 = qJDD(1) * t347 + qJDD(2);
t372 = qJD(1) * t353;
t351 = sin(qJ(1));
t354 = cos(qJ(1));
t333 = g(1) * t351 - g(2) * t354;
t355 = qJD(1) ^ 2;
t324 = qJDD(1) * pkin(1) + t355 * t380 + t333;
t334 = -g(1) * t354 - g(2) * t351;
t371 = qJDD(1) * t344;
t325 = -pkin(1) * t355 + pkin(8) * t371 + t334;
t374 = t324 * t376 + t325 * t353;
t279 = -t336 * pkin(2) + t337 * qJ(3) + (-g(3) * t350 + t326 * t372) * t344 + t374;
t328 = (qJD(2) * t372 + qJDD(1) * t350) * t344;
t370 = t350 * t373;
t329 = -qJD(2) * t370 + t353 * t371;
t280 = -t329 * pkin(2) - t379 - t328 * qJ(3) + (-t324 + (pkin(2) * t350 - qJ(3) * t353) * t338 * qJD(1)) * t344;
t343 = sin(pkin(11));
t346 = cos(pkin(11));
t317 = t338 * t343 + t346 * t370;
t242 = -0.2e1 * qJD(3) * t317 - t343 * t279 + t280 * t346;
t304 = t328 * t346 + t337 * t343;
t316 = t338 * t346 - t343 * t370;
t369 = t344 * t372;
t233 = (-t316 * t369 - t304) * pkin(9) + (t316 * t317 - t329) * pkin(3) + t242;
t243 = 0.2e1 * qJD(3) * t316 + t279 * t346 + t280 * t343;
t303 = -t328 * t343 + t337 * t346;
t305 = -pkin(3) * t369 - pkin(9) * t317;
t315 = t316 ^ 2;
t240 = -pkin(3) * t315 + pkin(9) * t303 + t305 * t369 + t243;
t349 = sin(qJ(4));
t222 = t233 * t349 + t240 * t381;
t297 = t316 * t349 + t317 * t381;
t264 = qJD(4) * t297 - t303 * t381 + t304 * t349;
t296 = -t316 * t381 + t317 * t349;
t274 = mrSges(5,1) * t296 + mrSges(5,2) * t297;
t332 = qJD(4) - t369;
t287 = mrSges(5,1) * t332 - mrSges(5,3) * t297;
t321 = qJDD(4) - t329;
t273 = pkin(4) * t296 - qJ(5) * t297;
t331 = t332 ^ 2;
t219 = -pkin(4) * t331 + qJ(5) * t321 - t273 * t296 + t222;
t293 = -g(3) * t377 + t324 * t375 - t325 * t350;
t278 = -t337 * pkin(2) - t336 * qJ(3) + t326 * t370 + qJDD(3) - t293;
t244 = -t303 * pkin(3) - t315 * pkin(9) + t305 * t317 + t278;
t265 = -qJD(4) * t296 + t303 * t349 + t304 * t381;
t225 = (t296 * t332 - t265) * qJ(5) + (t297 * t332 + t264) * pkin(4) + t244;
t342 = sin(pkin(12));
t345 = cos(pkin(12));
t285 = t297 * t345 + t332 * t342;
t214 = -0.2e1 * qJD(5) * t285 - t342 * t219 + t225 * t345;
t255 = t265 * t345 + t321 * t342;
t284 = -t297 * t342 + t332 * t345;
t212 = (t284 * t296 - t255) * pkin(10) + (t284 * t285 + t264) * pkin(5) + t214;
t215 = 0.2e1 * qJD(5) * t284 + t219 * t345 + t225 * t342;
t254 = -t265 * t342 + t321 * t345;
t268 = pkin(5) * t296 - pkin(10) * t285;
t283 = t284 ^ 2;
t213 = -pkin(5) * t283 + pkin(10) * t254 - t268 * t296 + t215;
t348 = sin(qJ(6));
t352 = cos(qJ(6));
t210 = t212 * t352 - t213 * t348;
t256 = t284 * t352 - t285 * t348;
t231 = qJD(6) * t256 + t254 * t348 + t255 * t352;
t257 = t284 * t348 + t285 * t352;
t241 = -mrSges(7,1) * t256 + mrSges(7,2) * t257;
t295 = qJD(6) + t296;
t245 = -mrSges(7,2) * t295 + mrSges(7,3) * t256;
t263 = qJDD(6) + t264;
t205 = m(7) * t210 + mrSges(7,1) * t263 - mrSges(7,3) * t231 - t241 * t257 + t245 * t295;
t211 = t212 * t348 + t213 * t352;
t230 = -qJD(6) * t257 + t254 * t352 - t255 * t348;
t246 = mrSges(7,1) * t295 - mrSges(7,3) * t257;
t206 = m(7) * t211 - mrSges(7,2) * t263 + mrSges(7,3) * t230 + t241 * t256 - t246 * t295;
t197 = t205 * t352 + t206 * t348;
t258 = -mrSges(6,1) * t284 + mrSges(6,2) * t285;
t266 = -mrSges(6,2) * t296 + mrSges(6,3) * t284;
t195 = m(6) * t214 + mrSges(6,1) * t264 - mrSges(6,3) * t255 - t258 * t285 + t266 * t296 + t197;
t267 = mrSges(6,1) * t296 - mrSges(6,3) * t285;
t365 = -t205 * t348 + t206 * t352;
t196 = m(6) * t215 - mrSges(6,2) * t264 + mrSges(6,3) * t254 + t258 * t284 - t267 * t296 + t365;
t366 = -t195 * t342 + t196 * t345;
t188 = m(5) * t222 - mrSges(5,2) * t321 - mrSges(5,3) * t264 - t274 * t296 - t287 * t332 + t366;
t221 = t233 * t381 - t240 * t349;
t286 = -mrSges(5,2) * t332 - mrSges(5,3) * t296;
t218 = -pkin(4) * t321 - qJ(5) * t331 + t273 * t297 + qJDD(5) - t221;
t216 = -pkin(5) * t254 - pkin(10) * t283 + t268 * t285 + t218;
t364 = m(7) * t216 - mrSges(7,1) * t230 + mrSges(7,2) * t231 - t245 * t256 + t246 * t257;
t359 = -m(6) * t218 + mrSges(6,1) * t254 - mrSges(6,2) * t255 + t266 * t284 - t267 * t285 - t364;
t201 = m(5) * t221 + mrSges(5,1) * t321 - mrSges(5,3) * t265 - t274 * t297 + t286 * t332 + t359;
t179 = t188 * t349 + t201 * t381;
t298 = -mrSges(4,1) * t316 + mrSges(4,2) * t317;
t301 = mrSges(4,2) * t369 + mrSges(4,3) * t316;
t177 = m(4) * t242 - mrSges(4,1) * t329 - mrSges(4,3) * t304 - t298 * t317 - t301 * t369 + t179;
t302 = -mrSges(4,1) * t369 - mrSges(4,3) * t317;
t367 = t188 * t381 - t201 * t349;
t178 = m(4) * t243 + mrSges(4,2) * t329 + mrSges(4,3) * t303 + t298 * t316 + t302 * t369 + t367;
t172 = t177 * t346 + t178 * t343;
t190 = t195 * t345 + t196 * t342;
t294 = -g(3) * t378 + t374;
t322 = mrSges(3,1) * t338 - mrSges(3,3) * t370;
t327 = (-mrSges(3,1) * t353 + mrSges(3,2) * t350) * t373;
t368 = -t177 * t343 + t178 * t346;
t170 = m(3) * t294 - mrSges(3,2) * t337 + mrSges(3,3) * t329 - t322 * t338 + t327 * t369 + t368;
t323 = -mrSges(3,2) * t338 + mrSges(3,3) * t369;
t361 = m(5) * t244 + mrSges(5,1) * t264 + mrSges(5,2) * t265 + t286 * t296 + t287 * t297 + t190;
t358 = -m(4) * t278 + mrSges(4,1) * t303 - mrSges(4,2) * t304 + t301 * t316 - t302 * t317 - t361;
t185 = m(3) * t293 + mrSges(3,1) * t337 - mrSges(3,3) * t328 + t323 * t338 - t327 * t370 + t358;
t166 = t170 * t353 - t185 * t350;
t309 = -t344 * t324 - t379;
t171 = m(3) * t309 - t329 * mrSges(3,1) + t328 * mrSges(3,2) + (t322 * t350 - t323 * t353) * t373 + t172;
t163 = t170 * t376 - t171 * t344 + t185 * t375;
t234 = Ifges(7,5) * t257 + Ifges(7,6) * t256 + Ifges(7,3) * t295;
t236 = Ifges(7,1) * t257 + Ifges(7,4) * t256 + Ifges(7,5) * t295;
t198 = -mrSges(7,1) * t216 + mrSges(7,3) * t211 + Ifges(7,4) * t231 + Ifges(7,2) * t230 + Ifges(7,6) * t263 - t234 * t257 + t236 * t295;
t235 = Ifges(7,4) * t257 + Ifges(7,2) * t256 + Ifges(7,6) * t295;
t199 = mrSges(7,2) * t216 - mrSges(7,3) * t210 + Ifges(7,1) * t231 + Ifges(7,4) * t230 + Ifges(7,5) * t263 + t234 * t256 - t235 * t295;
t247 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * t296;
t249 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * t296;
t181 = -mrSges(6,1) * t218 + mrSges(6,3) * t215 + Ifges(6,4) * t255 + Ifges(6,2) * t254 + Ifges(6,6) * t264 - pkin(5) * t364 + pkin(10) * t365 + t352 * t198 + t348 * t199 - t285 * t247 + t296 * t249;
t248 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * t296;
t183 = mrSges(6,2) * t218 - mrSges(6,3) * t214 + Ifges(6,1) * t255 + Ifges(6,4) * t254 + Ifges(6,5) * t264 - pkin(10) * t197 - t198 * t348 + t199 * t352 + t247 * t284 - t248 * t296;
t269 = Ifges(5,5) * t297 - Ifges(5,6) * t296 + Ifges(5,3) * t332;
t270 = Ifges(5,4) * t297 - Ifges(5,2) * t296 + Ifges(5,6) * t332;
t167 = mrSges(5,2) * t244 - mrSges(5,3) * t221 + Ifges(5,1) * t265 - Ifges(5,4) * t264 + Ifges(5,5) * t321 - qJ(5) * t190 - t181 * t342 + t183 * t345 - t269 * t296 - t270 * t332;
t271 = Ifges(5,1) * t297 - Ifges(5,4) * t296 + Ifges(5,5) * t332;
t362 = -mrSges(7,1) * t210 + mrSges(7,2) * t211 - Ifges(7,5) * t231 - Ifges(7,6) * t230 - Ifges(7,3) * t263 - t235 * t257 + t256 * t236;
t357 = -mrSges(6,1) * t214 + mrSges(6,2) * t215 - Ifges(6,5) * t255 - Ifges(6,6) * t254 - pkin(5) * t197 - t285 * t248 + t284 * t249 + t362;
t173 = (-Ifges(5,2) - Ifges(6,3)) * t264 + t332 * t271 + Ifges(5,6) * t321 - t297 * t269 + Ifges(5,4) * t265 - mrSges(5,1) * t244 + mrSges(5,3) * t222 - pkin(4) * t190 + t357;
t288 = Ifges(4,5) * t317 + Ifges(4,6) * t316 - Ifges(4,3) * t369;
t290 = Ifges(4,1) * t317 + Ifges(4,4) * t316 - Ifges(4,5) * t369;
t156 = -mrSges(4,1) * t278 + mrSges(4,3) * t243 + Ifges(4,4) * t304 + Ifges(4,2) * t303 - Ifges(4,6) * t329 - pkin(3) * t361 + pkin(9) * t367 + t349 * t167 + t173 * t381 - t317 * t288 - t290 * t369;
t289 = Ifges(4,4) * t317 + Ifges(4,2) * t316 - Ifges(4,6) * t369;
t159 = mrSges(4,2) * t278 - mrSges(4,3) * t242 + Ifges(4,1) * t304 + Ifges(4,4) * t303 - Ifges(4,5) * t329 - pkin(9) * t179 + t167 * t381 - t173 * t349 + t288 * t316 + t289 * t369;
t307 = Ifges(3,6) * t338 + (Ifges(3,4) * t350 + Ifges(3,2) * t353) * t373;
t308 = Ifges(3,5) * t338 + (Ifges(3,1) * t350 + Ifges(3,4) * t353) * t373;
t153 = Ifges(3,5) * t328 + Ifges(3,6) * t329 + Ifges(3,3) * t337 + mrSges(3,1) * t293 - mrSges(3,2) * t294 + t343 * t159 + t346 * t156 + pkin(2) * t358 + qJ(3) * t368 + (t307 * t350 - t308 * t353) * t373;
t306 = Ifges(3,3) * t338 + (Ifges(3,5) * t350 + Ifges(3,6) * t353) * t373;
t155 = mrSges(3,2) * t309 - mrSges(3,3) * t293 + Ifges(3,1) * t328 + Ifges(3,4) * t329 + Ifges(3,5) * t337 - qJ(3) * t172 - t156 * t343 + t159 * t346 + t306 * t369 - t307 * t338;
t360 = -mrSges(5,1) * t221 + mrSges(5,2) * t222 - Ifges(5,5) * t265 + Ifges(5,6) * t264 - Ifges(5,3) * t321 - pkin(4) * t359 - qJ(5) * t366 - t181 * t345 - t183 * t342 - t270 * t297 - t296 * t271;
t356 = -mrSges(4,1) * t242 + mrSges(4,2) * t243 - Ifges(4,5) * t304 - Ifges(4,6) * t303 - pkin(3) * t179 - t317 * t289 + t316 * t290 + t360;
t158 = -t306 * t370 + (Ifges(3,2) + Ifges(4,3)) * t329 + Ifges(3,6) * t337 + t338 * t308 + Ifges(3,4) * t328 - mrSges(3,1) * t309 + mrSges(3,3) * t294 - pkin(2) * t172 + t356;
t363 = mrSges(2,1) * t333 - mrSges(2,2) * t334 + Ifges(2,3) * qJDD(1) + pkin(1) * t163 + t153 * t347 + t155 * t378 + t158 * t377 + t166 * t380;
t164 = m(2) * t334 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t166;
t162 = t347 * t171 + (t170 * t350 + t185 * t353) * t344;
t160 = m(2) * t333 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t355 + t163;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t333 + Ifges(2,5) * qJDD(1) - t355 * Ifges(2,6) + t353 * t155 - t350 * t158 + (-t162 * t344 - t163 * t347) * pkin(8);
t150 = mrSges(2,1) * g(3) + mrSges(2,3) * t334 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t344 * t153 + (pkin(8) * t166 + t155 * t350 + t158 * t353) * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t151 - t351 * t150 - pkin(7) * (t160 * t354 + t164 * t351), t151, t155, t159, t167, t183, t199; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t151 + t354 * t150 + pkin(7) * (-t160 * t351 + t164 * t354), t150, t158, t156, t173, t181, t198; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t363, t363, t153, -Ifges(4,3) * t329 - t356, -t360, Ifges(6,3) * t264 - t357, -t362;];
m_new  = t1;
