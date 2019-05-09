% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 11:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 11:07:39
% EndTime: 2019-05-08 11:10:44
% DurationCPUTime: 89.05s
% Computational Cost: add. (1684091->399), mult. (3558843->515), div. (0->0), fcn. (2923669->14), ass. (0->167)
t344 = sin(pkin(6));
t383 = pkin(8) * t344;
t345 = cos(pkin(6));
t382 = t345 * g(3);
t350 = sin(qJ(2));
t381 = t344 * t350;
t356 = cos(qJ(2));
t380 = t344 * t356;
t379 = t345 * t350;
t378 = t345 * t356;
t376 = qJD(1) * t344;
t328 = (-pkin(2) * t356 - pkin(9) * t350) * t376;
t340 = qJD(1) * t345 + qJD(2);
t338 = t340 ^ 2;
t339 = qJDD(1) * t345 + qJDD(2);
t375 = qJD(1) * t356;
t351 = sin(qJ(1));
t357 = cos(qJ(1));
t335 = t351 * g(1) - g(2) * t357;
t358 = qJD(1) ^ 2;
t325 = qJDD(1) * pkin(1) + t358 * t383 + t335;
t336 = -g(1) * t357 - g(2) * t351;
t374 = qJDD(1) * t344;
t326 = -pkin(1) * t358 + pkin(8) * t374 + t336;
t377 = t325 * t379 + t356 * t326;
t281 = -t338 * pkin(2) + t339 * pkin(9) + (-g(3) * t350 + t328 * t375) * t344 + t377;
t329 = (qJD(2) * t375 + qJDD(1) * t350) * t344;
t373 = t350 * t376;
t330 = -qJD(2) * t373 + t356 * t374;
t282 = -t330 * pkin(2) - t329 * pkin(9) - t382 + (-t325 + (pkin(2) * t350 - pkin(9) * t356) * t340 * qJD(1)) * t344;
t349 = sin(qJ(3));
t355 = cos(qJ(3));
t249 = -t349 * t281 + t355 * t282;
t317 = t340 * t355 - t349 * t373;
t296 = qJD(3) * t317 + t329 * t355 + t339 * t349;
t318 = t340 * t349 + t355 * t373;
t322 = qJDD(3) - t330;
t372 = t344 * t375;
t334 = qJD(3) - t372;
t234 = (t317 * t334 - t296) * pkin(10) + (t317 * t318 + t322) * pkin(3) + t249;
t250 = t355 * t281 + t349 * t282;
t295 = -qJD(3) * t318 - t329 * t349 + t339 * t355;
t307 = pkin(3) * t334 - pkin(10) * t318;
t316 = t317 ^ 2;
t241 = -pkin(3) * t316 + pkin(10) * t295 - t307 * t334 + t250;
t348 = sin(qJ(4));
t354 = cos(qJ(4));
t231 = t348 * t234 + t354 * t241;
t303 = t317 * t348 + t318 * t354;
t261 = -t303 * qJD(4) + t295 * t354 - t348 * t296;
t302 = t317 * t354 - t348 * t318;
t275 = -mrSges(5,1) * t302 + mrSges(5,2) * t303;
t332 = qJD(4) + t334;
t287 = mrSges(5,1) * t332 - mrSges(5,3) * t303;
t321 = qJDD(4) + t322;
t276 = -pkin(4) * t302 - pkin(11) * t303;
t331 = t332 ^ 2;
t220 = -pkin(4) * t331 + pkin(11) * t321 + t276 * t302 + t231;
t297 = -g(3) * t380 + t325 * t378 - t350 * t326;
t280 = -t339 * pkin(2) - t338 * pkin(9) + t328 * t373 - t297;
t246 = -t295 * pkin(3) - t316 * pkin(10) + t318 * t307 + t280;
t262 = qJD(4) * t302 + t295 * t348 + t296 * t354;
t223 = (-t302 * t332 - t262) * pkin(11) + (t303 * t332 - t261) * pkin(4) + t246;
t347 = sin(qJ(5));
t353 = cos(qJ(5));
t215 = -t347 * t220 + t353 * t223;
t284 = -t303 * t347 + t332 * t353;
t245 = qJD(5) * t284 + t262 * t353 + t321 * t347;
t260 = qJDD(5) - t261;
t285 = t303 * t353 + t332 * t347;
t301 = qJD(5) - t302;
t213 = (t284 * t301 - t245) * pkin(12) + (t284 * t285 + t260) * pkin(5) + t215;
t216 = t353 * t220 + t347 * t223;
t244 = -qJD(5) * t285 - t262 * t347 + t321 * t353;
t270 = pkin(5) * t301 - pkin(12) * t285;
t283 = t284 ^ 2;
t214 = -pkin(5) * t283 + pkin(12) * t244 - t270 * t301 + t216;
t346 = sin(qJ(6));
t352 = cos(qJ(6));
t211 = t213 * t352 - t214 * t346;
t265 = t284 * t352 - t285 * t346;
t228 = qJD(6) * t265 + t244 * t346 + t245 * t352;
t266 = t284 * t346 + t285 * t352;
t242 = -mrSges(7,1) * t265 + mrSges(7,2) * t266;
t299 = qJD(6) + t301;
t247 = -mrSges(7,2) * t299 + mrSges(7,3) * t265;
t256 = qJDD(6) + t260;
t206 = m(7) * t211 + mrSges(7,1) * t256 - mrSges(7,3) * t228 - t242 * t266 + t247 * t299;
t212 = t213 * t346 + t214 * t352;
t227 = -qJD(6) * t266 + t244 * t352 - t245 * t346;
t248 = mrSges(7,1) * t299 - mrSges(7,3) * t266;
t207 = m(7) * t212 - mrSges(7,2) * t256 + mrSges(7,3) * t227 + t242 * t265 - t248 * t299;
t198 = t352 * t206 + t346 * t207;
t267 = -mrSges(6,1) * t284 + mrSges(6,2) * t285;
t268 = -mrSges(6,2) * t301 + mrSges(6,3) * t284;
t196 = m(6) * t215 + mrSges(6,1) * t260 - mrSges(6,3) * t245 - t267 * t285 + t268 * t301 + t198;
t269 = mrSges(6,1) * t301 - mrSges(6,3) * t285;
t368 = -t206 * t346 + t352 * t207;
t197 = m(6) * t216 - mrSges(6,2) * t260 + mrSges(6,3) * t244 + t267 * t284 - t269 * t301 + t368;
t369 = -t196 * t347 + t353 * t197;
t189 = m(5) * t231 - mrSges(5,2) * t321 + mrSges(5,3) * t261 + t275 * t302 - t287 * t332 + t369;
t230 = t234 * t354 - t348 * t241;
t286 = -mrSges(5,2) * t332 + mrSges(5,3) * t302;
t219 = -pkin(4) * t321 - pkin(11) * t331 + t303 * t276 - t230;
t217 = -pkin(5) * t244 - pkin(12) * t283 + t270 * t285 + t219;
t367 = m(7) * t217 - t227 * mrSges(7,1) + mrSges(7,2) * t228 - t265 * t247 + t248 * t266;
t362 = -m(6) * t219 + t244 * mrSges(6,1) - mrSges(6,2) * t245 + t284 * t268 - t269 * t285 - t367;
t202 = m(5) * t230 + mrSges(5,1) * t321 - mrSges(5,3) * t262 - t275 * t303 + t286 * t332 + t362;
t180 = t348 * t189 + t354 * t202;
t304 = -mrSges(4,1) * t317 + mrSges(4,2) * t318;
t305 = -mrSges(4,2) * t334 + mrSges(4,3) * t317;
t178 = m(4) * t249 + mrSges(4,1) * t322 - mrSges(4,3) * t296 - t304 * t318 + t305 * t334 + t180;
t306 = mrSges(4,1) * t334 - mrSges(4,3) * t318;
t370 = t354 * t189 - t202 * t348;
t179 = m(4) * t250 - mrSges(4,2) * t322 + mrSges(4,3) * t295 + t304 * t317 - t306 * t334 + t370;
t173 = t355 * t178 + t349 * t179;
t191 = t353 * t196 + t347 * t197;
t298 = -g(3) * t381 + t377;
t323 = mrSges(3,1) * t340 - mrSges(3,3) * t373;
t327 = (-mrSges(3,1) * t356 + mrSges(3,2) * t350) * t376;
t371 = -t178 * t349 + t355 * t179;
t171 = m(3) * t298 - mrSges(3,2) * t339 + mrSges(3,3) * t330 - t323 * t340 + t327 * t372 + t371;
t324 = -mrSges(3,2) * t340 + mrSges(3,3) * t372;
t364 = m(5) * t246 - t261 * mrSges(5,1) + mrSges(5,2) * t262 - t302 * t286 + t287 * t303 + t191;
t361 = -m(4) * t280 + t295 * mrSges(4,1) - mrSges(4,2) * t296 + t317 * t305 - t306 * t318 - t364;
t186 = m(3) * t297 + mrSges(3,1) * t339 - mrSges(3,3) * t329 + t324 * t340 - t327 * t373 + t361;
t167 = t356 * t171 - t186 * t350;
t311 = -t344 * t325 - t382;
t172 = m(3) * t311 - t330 * mrSges(3,1) + t329 * mrSges(3,2) + (t323 * t350 - t324 * t356) * t376 + t173;
t164 = t171 * t379 - t172 * t344 + t186 * t378;
t235 = Ifges(7,5) * t266 + Ifges(7,6) * t265 + Ifges(7,3) * t299;
t237 = Ifges(7,1) * t266 + Ifges(7,4) * t265 + Ifges(7,5) * t299;
t199 = -mrSges(7,1) * t217 + mrSges(7,3) * t212 + Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t256 - t235 * t266 + t237 * t299;
t236 = Ifges(7,4) * t266 + Ifges(7,2) * t265 + Ifges(7,6) * t299;
t200 = mrSges(7,2) * t217 - mrSges(7,3) * t211 + Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t256 + t235 * t265 - t236 * t299;
t251 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * t301;
t253 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * t301;
t182 = -mrSges(6,1) * t219 + mrSges(6,3) * t216 + Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t260 - pkin(5) * t367 + pkin(12) * t368 + t352 * t199 + t346 * t200 - t285 * t251 + t301 * t253;
t252 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * t301;
t184 = mrSges(6,2) * t219 - mrSges(6,3) * t215 + Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t260 - pkin(12) * t198 - t199 * t346 + t200 * t352 + t251 * t284 - t252 * t301;
t271 = Ifges(5,5) * t303 + Ifges(5,6) * t302 + Ifges(5,3) * t332;
t272 = Ifges(5,4) * t303 + Ifges(5,2) * t302 + Ifges(5,6) * t332;
t168 = mrSges(5,2) * t246 - mrSges(5,3) * t230 + Ifges(5,1) * t262 + Ifges(5,4) * t261 + Ifges(5,5) * t321 - pkin(11) * t191 - t182 * t347 + t184 * t353 + t271 * t302 - t272 * t332;
t273 = Ifges(5,1) * t303 + Ifges(5,4) * t302 + Ifges(5,5) * t332;
t365 = -mrSges(7,1) * t211 + mrSges(7,2) * t212 - Ifges(7,5) * t228 - Ifges(7,6) * t227 - Ifges(7,3) * t256 - t266 * t236 + t265 * t237;
t360 = mrSges(6,1) * t215 - mrSges(6,2) * t216 + Ifges(6,5) * t245 + Ifges(6,6) * t244 + Ifges(6,3) * t260 + pkin(5) * t198 + t285 * t252 - t284 * t253 - t365;
t174 = -mrSges(5,1) * t246 + mrSges(5,3) * t231 + Ifges(5,4) * t262 + Ifges(5,2) * t261 + Ifges(5,6) * t321 - pkin(4) * t191 - t303 * t271 + t332 * t273 - t360;
t289 = Ifges(4,5) * t318 + Ifges(4,6) * t317 + Ifges(4,3) * t334;
t291 = Ifges(4,1) * t318 + Ifges(4,4) * t317 + Ifges(4,5) * t334;
t157 = -mrSges(4,1) * t280 + mrSges(4,3) * t250 + Ifges(4,4) * t296 + Ifges(4,2) * t295 + Ifges(4,6) * t322 - pkin(3) * t364 + pkin(10) * t370 + t348 * t168 + t354 * t174 - t318 * t289 + t334 * t291;
t290 = Ifges(4,4) * t318 + Ifges(4,2) * t317 + Ifges(4,6) * t334;
t160 = mrSges(4,2) * t280 - mrSges(4,3) * t249 + Ifges(4,1) * t296 + Ifges(4,4) * t295 + Ifges(4,5) * t322 - pkin(10) * t180 + t168 * t354 - t174 * t348 + t289 * t317 - t290 * t334;
t309 = Ifges(3,6) * t340 + (Ifges(3,4) * t350 + Ifges(3,2) * t356) * t376;
t310 = Ifges(3,5) * t340 + (Ifges(3,1) * t350 + Ifges(3,4) * t356) * t376;
t154 = Ifges(3,5) * t329 + Ifges(3,6) * t330 + Ifges(3,3) * t339 + mrSges(3,1) * t297 - mrSges(3,2) * t298 + t349 * t160 + t355 * t157 + pkin(2) * t361 + pkin(9) * t371 + (t309 * t350 - t310 * t356) * t376;
t308 = Ifges(3,3) * t340 + (Ifges(3,5) * t350 + Ifges(3,6) * t356) * t376;
t156 = mrSges(3,2) * t311 - mrSges(3,3) * t297 + Ifges(3,1) * t329 + Ifges(3,4) * t330 + Ifges(3,5) * t339 - pkin(9) * t173 - t157 * t349 + t160 * t355 + t308 * t372 - t309 * t340;
t363 = -mrSges(5,1) * t230 + mrSges(5,2) * t231 - Ifges(5,5) * t262 - Ifges(5,6) * t261 - Ifges(5,3) * t321 - pkin(4) * t362 - pkin(11) * t369 - t353 * t182 - t347 * t184 - t303 * t272 + t302 * t273;
t359 = mrSges(4,1) * t249 - mrSges(4,2) * t250 + Ifges(4,5) * t296 + Ifges(4,6) * t295 + Ifges(4,3) * t322 + pkin(3) * t180 + t318 * t290 - t317 * t291 - t363;
t159 = -mrSges(3,1) * t311 + mrSges(3,3) * t298 + Ifges(3,4) * t329 + Ifges(3,2) * t330 + Ifges(3,6) * t339 - pkin(2) * t173 - t308 * t373 + t340 * t310 - t359;
t366 = mrSges(2,1) * t335 - mrSges(2,2) * t336 + Ifges(2,3) * qJDD(1) + pkin(1) * t164 + t345 * t154 + t156 * t381 + t159 * t380 + t167 * t383;
t165 = m(2) * t336 - mrSges(2,1) * t358 - qJDD(1) * mrSges(2,2) + t167;
t163 = t345 * t172 + (t171 * t350 + t186 * t356) * t344;
t161 = m(2) * t335 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t358 + t164;
t152 = -mrSges(2,2) * g(3) - mrSges(2,3) * t335 + Ifges(2,5) * qJDD(1) - t358 * Ifges(2,6) + t356 * t156 - t350 * t159 + (-t163 * t344 - t164 * t345) * pkin(8);
t151 = mrSges(2,1) * g(3) + mrSges(2,3) * t336 + t358 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t163 - t344 * t154 + (pkin(8) * t167 + t156 * t350 + t159 * t356) * t345;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t357 * t152 - t351 * t151 - pkin(7) * (t161 * t357 + t165 * t351), t152, t156, t160, t168, t184, t200; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t152 + t357 * t151 + pkin(7) * (-t161 * t351 + t165 * t357), t151, t159, t157, t174, t182, t199; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t154, t359, -t363, t360, -t365;];
m_new  = t1;
