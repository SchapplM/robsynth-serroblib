% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 23:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:37:19
% EndTime: 2019-05-06 23:39:34
% DurationCPUTime: 81.47s
% Computational Cost: add. (1455040->399), mult. (3288211->517), div. (0->0), fcn. (2689739->14), ass. (0->165)
t345 = sin(pkin(6));
t383 = pkin(8) * t345;
t347 = cos(pkin(6));
t382 = t347 * g(3);
t351 = sin(qJ(2));
t381 = t345 * t351;
t356 = cos(qJ(2));
t380 = t345 * t356;
t379 = t347 * t351;
t378 = t347 * t356;
t376 = qJD(1) * t345;
t328 = (-pkin(2) * t356 - qJ(3) * t351) * t376;
t340 = t347 * qJD(1) + qJD(2);
t338 = t340 ^ 2;
t339 = t347 * qJDD(1) + qJDD(2);
t375 = qJD(1) * t356;
t352 = sin(qJ(1));
t357 = cos(qJ(1));
t335 = t352 * g(1) - t357 * g(2);
t358 = qJD(1) ^ 2;
t326 = qJDD(1) * pkin(1) + t358 * t383 + t335;
t336 = -t357 * g(1) - t352 * g(2);
t374 = qJDD(1) * t345;
t327 = -t358 * pkin(1) + pkin(8) * t374 + t336;
t377 = t326 * t379 + t356 * t327;
t281 = -t338 * pkin(2) + t339 * qJ(3) + (-g(3) * t351 + t328 * t375) * t345 + t377;
t330 = (qJD(2) * t375 + qJDD(1) * t351) * t345;
t373 = t351 * t376;
t331 = -qJD(2) * t373 + t356 * t374;
t282 = -t331 * pkin(2) - t382 - t330 * qJ(3) + (-t326 + (pkin(2) * t351 - qJ(3) * t356) * t340 * qJD(1)) * t345;
t344 = sin(pkin(12));
t346 = cos(pkin(12));
t319 = t344 * t340 + t346 * t373;
t243 = -0.2e1 * qJD(3) * t319 - t344 * t281 + t346 * t282;
t306 = t346 * t330 + t344 * t339;
t318 = t346 * t340 - t344 * t373;
t372 = t345 * t375;
t234 = (-t318 * t372 - t306) * pkin(9) + (t318 * t319 - t331) * pkin(3) + t243;
t244 = 0.2e1 * qJD(3) * t318 + t346 * t281 + t344 * t282;
t305 = -t344 * t330 + t346 * t339;
t307 = -pkin(3) * t372 - t319 * pkin(9);
t317 = t318 ^ 2;
t241 = -t317 * pkin(3) + t305 * pkin(9) + t307 * t372 + t244;
t350 = sin(qJ(4));
t355 = cos(qJ(4));
t223 = t350 * t234 + t355 * t241;
t299 = t350 * t318 + t355 * t319;
t266 = -t299 * qJD(4) + t355 * t305 - t350 * t306;
t298 = t355 * t318 - t350 * t319;
t275 = -t298 * mrSges(5,1) + t299 * mrSges(5,2);
t334 = qJD(4) - t372;
t287 = t334 * mrSges(5,1) - t299 * mrSges(5,3);
t323 = qJDD(4) - t331;
t276 = -t298 * pkin(4) - t299 * pkin(10);
t333 = t334 ^ 2;
t220 = -t333 * pkin(4) + t323 * pkin(10) + t298 * t276 + t223;
t294 = -g(3) * t380 + t326 * t378 - t351 * t327;
t280 = -t339 * pkin(2) - t338 * qJ(3) + t328 * t373 + qJDD(3) - t294;
t248 = -t305 * pkin(3) - t317 * pkin(9) + t319 * t307 + t280;
t267 = t298 * qJD(4) + t350 * t305 + t355 * t306;
t231 = (-t298 * t334 - t267) * pkin(10) + (t299 * t334 - t266) * pkin(4) + t248;
t349 = sin(qJ(5));
t354 = cos(qJ(5));
t215 = -t349 * t220 + t354 * t231;
t284 = -t349 * t299 + t354 * t334;
t247 = t284 * qJD(5) + t354 * t267 + t349 * t323;
t265 = qJDD(5) - t266;
t285 = t354 * t299 + t349 * t334;
t297 = qJD(5) - t298;
t213 = (t284 * t297 - t247) * pkin(11) + (t284 * t285 + t265) * pkin(5) + t215;
t216 = t354 * t220 + t349 * t231;
t246 = -t285 * qJD(5) - t349 * t267 + t354 * t323;
t270 = t297 * pkin(5) - t285 * pkin(11);
t283 = t284 ^ 2;
t214 = -t283 * pkin(5) + t246 * pkin(11) - t297 * t270 + t216;
t348 = sin(qJ(6));
t353 = cos(qJ(6));
t211 = t353 * t213 - t348 * t214;
t257 = t353 * t284 - t348 * t285;
t228 = t257 * qJD(6) + t348 * t246 + t353 * t247;
t258 = t348 * t284 + t353 * t285;
t242 = -t257 * mrSges(7,1) + t258 * mrSges(7,2);
t293 = qJD(6) + t297;
t249 = -t293 * mrSges(7,2) + t257 * mrSges(7,3);
t261 = qJDD(6) + t265;
t206 = m(7) * t211 + t261 * mrSges(7,1) - t228 * mrSges(7,3) - t258 * t242 + t293 * t249;
t212 = t348 * t213 + t353 * t214;
t227 = -t258 * qJD(6) + t353 * t246 - t348 * t247;
t250 = t293 * mrSges(7,1) - t258 * mrSges(7,3);
t207 = m(7) * t212 - t261 * mrSges(7,2) + t227 * mrSges(7,3) + t257 * t242 - t293 * t250;
t198 = t353 * t206 + t348 * t207;
t259 = -t284 * mrSges(6,1) + t285 * mrSges(6,2);
t268 = -t297 * mrSges(6,2) + t284 * mrSges(6,3);
t196 = m(6) * t215 + t265 * mrSges(6,1) - t247 * mrSges(6,3) - t285 * t259 + t297 * t268 + t198;
t269 = t297 * mrSges(6,1) - t285 * mrSges(6,3);
t368 = -t348 * t206 + t353 * t207;
t197 = m(6) * t216 - t265 * mrSges(6,2) + t246 * mrSges(6,3) + t284 * t259 - t297 * t269 + t368;
t369 = -t349 * t196 + t354 * t197;
t189 = m(5) * t223 - t323 * mrSges(5,2) + t266 * mrSges(5,3) + t298 * t275 - t334 * t287 + t369;
t222 = t355 * t234 - t350 * t241;
t286 = -t334 * mrSges(5,2) + t298 * mrSges(5,3);
t219 = -t323 * pkin(4) - t333 * pkin(10) + t299 * t276 - t222;
t217 = -t246 * pkin(5) - t283 * pkin(11) + t285 * t270 + t219;
t367 = m(7) * t217 - t227 * mrSges(7,1) + t228 * mrSges(7,2) - t257 * t249 + t258 * t250;
t362 = -m(6) * t219 + t246 * mrSges(6,1) - t247 * mrSges(6,2) + t284 * t268 - t285 * t269 - t367;
t202 = m(5) * t222 + t323 * mrSges(5,1) - t267 * mrSges(5,3) - t299 * t275 + t334 * t286 + t362;
t180 = t350 * t189 + t355 * t202;
t300 = -t318 * mrSges(4,1) + t319 * mrSges(4,2);
t303 = mrSges(4,2) * t372 + t318 * mrSges(4,3);
t178 = m(4) * t243 - t331 * mrSges(4,1) - t306 * mrSges(4,3) - t319 * t300 - t303 * t372 + t180;
t304 = -mrSges(4,1) * t372 - t319 * mrSges(4,3);
t370 = t355 * t189 - t350 * t202;
t179 = m(4) * t244 + t331 * mrSges(4,2) + t305 * mrSges(4,3) + t318 * t300 + t304 * t372 + t370;
t173 = t346 * t178 + t344 * t179;
t191 = t354 * t196 + t349 * t197;
t295 = -g(3) * t381 + t377;
t324 = t340 * mrSges(3,1) - mrSges(3,3) * t373;
t329 = (-mrSges(3,1) * t356 + mrSges(3,2) * t351) * t376;
t371 = -t344 * t178 + t346 * t179;
t171 = m(3) * t295 - t339 * mrSges(3,2) + t331 * mrSges(3,3) - t340 * t324 + t329 * t372 + t371;
t325 = -t340 * mrSges(3,2) + mrSges(3,3) * t372;
t364 = m(5) * t248 - t266 * mrSges(5,1) + t267 * mrSges(5,2) - t298 * t286 + t299 * t287 + t191;
t361 = -m(4) * t280 + t305 * mrSges(4,1) - t306 * mrSges(4,2) + t318 * t303 - t319 * t304 - t364;
t186 = m(3) * t294 + t339 * mrSges(3,1) - t330 * mrSges(3,3) + t340 * t325 - t329 * t373 + t361;
t167 = t356 * t171 - t351 * t186;
t311 = -t345 * t326 - t382;
t172 = m(3) * t311 - t331 * mrSges(3,1) + t330 * mrSges(3,2) + (t324 * t351 - t325 * t356) * t376 + t173;
t164 = t171 * t379 - t345 * t172 + t186 * t378;
t235 = Ifges(7,5) * t258 + Ifges(7,6) * t257 + Ifges(7,3) * t293;
t237 = Ifges(7,1) * t258 + Ifges(7,4) * t257 + Ifges(7,5) * t293;
t199 = -mrSges(7,1) * t217 + mrSges(7,3) * t212 + Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t261 - t258 * t235 + t293 * t237;
t236 = Ifges(7,4) * t258 + Ifges(7,2) * t257 + Ifges(7,6) * t293;
t200 = mrSges(7,2) * t217 - mrSges(7,3) * t211 + Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t261 + t257 * t235 - t293 * t236;
t251 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * t297;
t253 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * t297;
t182 = -mrSges(6,1) * t219 + mrSges(6,3) * t216 + Ifges(6,4) * t247 + Ifges(6,2) * t246 + Ifges(6,6) * t265 - pkin(5) * t367 + pkin(11) * t368 + t353 * t199 + t348 * t200 - t285 * t251 + t297 * t253;
t252 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * t297;
t184 = mrSges(6,2) * t219 - mrSges(6,3) * t215 + Ifges(6,1) * t247 + Ifges(6,4) * t246 + Ifges(6,5) * t265 - pkin(11) * t198 - t348 * t199 + t353 * t200 + t284 * t251 - t297 * t252;
t271 = Ifges(5,5) * t299 + Ifges(5,6) * t298 + Ifges(5,3) * t334;
t272 = Ifges(5,4) * t299 + Ifges(5,2) * t298 + Ifges(5,6) * t334;
t168 = mrSges(5,2) * t248 - mrSges(5,3) * t222 + Ifges(5,1) * t267 + Ifges(5,4) * t266 + Ifges(5,5) * t323 - pkin(10) * t191 - t349 * t182 + t354 * t184 + t298 * t271 - t334 * t272;
t273 = Ifges(5,1) * t299 + Ifges(5,4) * t298 + Ifges(5,5) * t334;
t365 = -mrSges(7,1) * t211 + mrSges(7,2) * t212 - Ifges(7,5) * t228 - Ifges(7,6) * t227 - Ifges(7,3) * t261 - t258 * t236 + t257 * t237;
t360 = mrSges(6,1) * t215 - mrSges(6,2) * t216 + Ifges(6,5) * t247 + Ifges(6,6) * t246 + Ifges(6,3) * t265 + pkin(5) * t198 + t285 * t252 - t284 * t253 - t365;
t174 = -mrSges(5,1) * t248 + mrSges(5,3) * t223 + Ifges(5,4) * t267 + Ifges(5,2) * t266 + Ifges(5,6) * t323 - pkin(4) * t191 - t299 * t271 + t334 * t273 - t360;
t288 = Ifges(4,5) * t319 + Ifges(4,6) * t318 - Ifges(4,3) * t372;
t290 = Ifges(4,1) * t319 + Ifges(4,4) * t318 - Ifges(4,5) * t372;
t157 = -mrSges(4,1) * t280 + mrSges(4,3) * t244 + Ifges(4,4) * t306 + Ifges(4,2) * t305 - Ifges(4,6) * t331 - pkin(3) * t364 + pkin(9) * t370 + t350 * t168 + t355 * t174 - t319 * t288 - t290 * t372;
t289 = Ifges(4,4) * t319 + Ifges(4,2) * t318 - Ifges(4,6) * t372;
t160 = mrSges(4,2) * t280 - mrSges(4,3) * t243 + Ifges(4,1) * t306 + Ifges(4,4) * t305 - Ifges(4,5) * t331 - pkin(9) * t180 + t355 * t168 - t350 * t174 + t318 * t288 + t289 * t372;
t309 = Ifges(3,6) * t340 + (Ifges(3,4) * t351 + Ifges(3,2) * t356) * t376;
t310 = Ifges(3,5) * t340 + (Ifges(3,1) * t351 + Ifges(3,4) * t356) * t376;
t154 = Ifges(3,5) * t330 + Ifges(3,6) * t331 + Ifges(3,3) * t339 + mrSges(3,1) * t294 - mrSges(3,2) * t295 + t344 * t160 + t346 * t157 + pkin(2) * t361 + qJ(3) * t371 + (t309 * t351 - t310 * t356) * t376;
t308 = Ifges(3,3) * t340 + (Ifges(3,5) * t351 + Ifges(3,6) * t356) * t376;
t156 = mrSges(3,2) * t311 - mrSges(3,3) * t294 + Ifges(3,1) * t330 + Ifges(3,4) * t331 + Ifges(3,5) * t339 - qJ(3) * t173 - t344 * t157 + t346 * t160 + t308 * t372 - t340 * t309;
t363 = -mrSges(5,1) * t222 + mrSges(5,2) * t223 - Ifges(5,5) * t267 - Ifges(5,6) * t266 - Ifges(5,3) * t323 - pkin(4) * t362 - pkin(10) * t369 - t354 * t182 - t349 * t184 - t299 * t272 + t298 * t273;
t359 = -mrSges(4,1) * t243 + mrSges(4,2) * t244 - Ifges(4,5) * t306 - Ifges(4,6) * t305 - pkin(3) * t180 - t319 * t289 + t318 * t290 + t363;
t159 = (Ifges(3,2) + Ifges(4,3)) * t331 - t308 * t373 + t340 * t310 + Ifges(3,6) * t339 + Ifges(3,4) * t330 - mrSges(3,1) * t311 + mrSges(3,3) * t295 - pkin(2) * t173 + t359;
t366 = mrSges(2,1) * t335 - mrSges(2,2) * t336 + Ifges(2,3) * qJDD(1) + pkin(1) * t164 + t347 * t154 + t156 * t381 + t159 * t380 + t167 * t383;
t165 = m(2) * t336 - t358 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t167;
t163 = t347 * t172 + (t171 * t351 + t186 * t356) * t345;
t161 = m(2) * t335 + qJDD(1) * mrSges(2,1) - t358 * mrSges(2,2) + t164;
t152 = -mrSges(2,2) * g(3) - mrSges(2,3) * t335 + Ifges(2,5) * qJDD(1) - t358 * Ifges(2,6) + t356 * t156 - t351 * t159 + (-t163 * t345 - t164 * t347) * pkin(8);
t151 = mrSges(2,1) * g(3) + mrSges(2,3) * t336 + t358 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t163 - t345 * t154 + (pkin(8) * t167 + t156 * t351 + t159 * t356) * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t357 * t152 - t352 * t151 - pkin(7) * (t357 * t161 + t352 * t165), t152, t156, t160, t168, t184, t200; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t352 * t152 + t357 * t151 + pkin(7) * (-t352 * t161 + t357 * t165), t151, t159, t157, t174, t182, t199; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t154, -Ifges(4,3) * t331 - t359, -t363, t360, -t365;];
m_new  = t1;
