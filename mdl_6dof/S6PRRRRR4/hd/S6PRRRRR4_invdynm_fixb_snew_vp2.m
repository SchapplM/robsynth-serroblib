% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRR4
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:24:12
% EndTime: 2019-05-05 11:25:33
% DurationCPUTime: 68.57s
% Computational Cost: add. (1335267->352), mult. (2757201->466), div. (0->0), fcn. (2246135->16), ass. (0->158)
t314 = sin(pkin(13));
t317 = cos(pkin(13));
t305 = t314 * g(1) - t317 * g(2);
t306 = -t317 * g(1) - t314 * g(2);
t313 = -g(3) + qJDD(1);
t324 = sin(qJ(2));
t319 = cos(pkin(6));
t329 = cos(qJ(2));
t348 = t319 * t329;
t316 = sin(pkin(6));
t351 = t316 * t329;
t274 = t305 * t348 - t324 * t306 + t313 * t351;
t315 = sin(pkin(7));
t330 = qJD(2) ^ 2;
t271 = t330 * t315 * pkin(9) + qJDD(2) * pkin(2) + t274;
t349 = t319 * t324;
t352 = t316 * t324;
t275 = t305 * t349 + t329 * t306 + t313 * t352;
t345 = qJDD(2) * t315;
t272 = -t330 * pkin(2) + pkin(9) * t345 + t275;
t289 = -t316 * t305 + t319 * t313;
t318 = cos(pkin(7));
t323 = sin(qJ(3));
t328 = cos(qJ(3));
t242 = -t323 * t272 + (t271 * t318 + t289 * t315) * t328;
t350 = t318 * t323;
t353 = t315 * t323;
t243 = t271 * t350 + t328 * t272 + t289 * t353;
t312 = t318 * qJD(2) + qJD(3);
t347 = qJD(2) * t315;
t344 = t323 * t347;
t292 = t312 * mrSges(4,1) - mrSges(4,3) * t344;
t294 = (-mrSges(4,1) * t328 + mrSges(4,2) * t323) * t347;
t297 = -qJD(3) * t344 + t328 * t345;
t311 = t318 * qJDD(2) + qJDD(3);
t295 = (-pkin(3) * t328 - pkin(10) * t323) * t347;
t310 = t312 ^ 2;
t346 = qJD(2) * t328;
t343 = t315 * t346;
t230 = -t310 * pkin(3) + t311 * pkin(10) + t295 * t343 + t243;
t284 = t318 * t289;
t296 = (qJD(3) * t346 + qJDD(2) * t323) * t315;
t240 = -t297 * pkin(3) - t296 * pkin(10) + t284 + (-t271 + (pkin(3) * t323 - pkin(10) * t328) * t312 * qJD(2)) * t315;
t322 = sin(qJ(4));
t327 = cos(qJ(4));
t218 = -t322 * t230 + t327 * t240;
t287 = t327 * t312 - t322 * t344;
t264 = t287 * qJD(4) + t327 * t296 + t322 * t311;
t288 = t322 * t312 + t327 * t344;
t291 = qJDD(4) - t297;
t304 = qJD(4) - t343;
t215 = (t287 * t304 - t264) * pkin(11) + (t287 * t288 + t291) * pkin(4) + t218;
t219 = t327 * t230 + t322 * t240;
t263 = -t288 * qJD(4) - t322 * t296 + t327 * t311;
t278 = t304 * pkin(4) - t288 * pkin(11);
t286 = t287 ^ 2;
t217 = -t286 * pkin(4) + t263 * pkin(11) - t304 * t278 + t219;
t321 = sin(qJ(5));
t326 = cos(qJ(5));
t212 = t321 * t215 + t326 * t217;
t270 = t321 * t287 + t326 * t288;
t235 = -t270 * qJD(5) + t326 * t263 - t321 * t264;
t269 = t326 * t287 - t321 * t288;
t250 = -t269 * mrSges(6,1) + t270 * mrSges(6,2);
t303 = qJD(5) + t304;
t256 = t303 * mrSges(6,1) - t270 * mrSges(6,3);
t290 = qJDD(5) + t291;
t251 = -t269 * pkin(5) - t270 * pkin(12);
t301 = t303 ^ 2;
t209 = -t301 * pkin(5) + t290 * pkin(12) + t269 * t251 + t212;
t229 = -t311 * pkin(3) - t310 * pkin(10) + t295 * t344 - t242;
t220 = -t263 * pkin(4) - t286 * pkin(11) + t288 * t278 + t229;
t236 = t269 * qJD(5) + t321 * t263 + t326 * t264;
t213 = (-t269 * t303 - t236) * pkin(12) + (t270 * t303 - t235) * pkin(5) + t220;
t320 = sin(qJ(6));
t325 = cos(qJ(6));
t206 = -t320 * t209 + t325 * t213;
t253 = -t320 * t270 + t325 * t303;
t223 = t253 * qJD(6) + t325 * t236 + t320 * t290;
t234 = qJDD(6) - t235;
t254 = t325 * t270 + t320 * t303;
t241 = -t253 * mrSges(7,1) + t254 * mrSges(7,2);
t266 = qJD(6) - t269;
t244 = -t266 * mrSges(7,2) + t253 * mrSges(7,3);
t202 = m(7) * t206 + t234 * mrSges(7,1) - t223 * mrSges(7,3) - t254 * t241 + t266 * t244;
t207 = t325 * t209 + t320 * t213;
t222 = -t254 * qJD(6) - t320 * t236 + t325 * t290;
t245 = t266 * mrSges(7,1) - t254 * mrSges(7,3);
t203 = m(7) * t207 - t234 * mrSges(7,2) + t222 * mrSges(7,3) + t253 * t241 - t266 * t245;
t340 = -t320 * t202 + t325 * t203;
t189 = m(6) * t212 - t290 * mrSges(6,2) + t235 * mrSges(6,3) + t269 * t250 - t303 * t256 + t340;
t211 = t326 * t215 - t321 * t217;
t255 = -t303 * mrSges(6,2) + t269 * mrSges(6,3);
t208 = -t290 * pkin(5) - t301 * pkin(12) + t270 * t251 - t211;
t336 = -m(7) * t208 + t222 * mrSges(7,1) - t223 * mrSges(7,2) + t253 * t244 - t254 * t245;
t198 = m(6) * t211 + t290 * mrSges(6,1) - t236 * mrSges(6,3) - t270 * t250 + t303 * t255 + t336;
t183 = t321 * t189 + t326 * t198;
t273 = -t287 * mrSges(5,1) + t288 * mrSges(5,2);
t276 = -t304 * mrSges(5,2) + t287 * mrSges(5,3);
t181 = m(5) * t218 + t291 * mrSges(5,1) - t264 * mrSges(5,3) - t288 * t273 + t304 * t276 + t183;
t277 = t304 * mrSges(5,1) - t288 * mrSges(5,3);
t341 = t326 * t189 - t321 * t198;
t182 = m(5) * t219 - t291 * mrSges(5,2) + t263 * mrSges(5,3) + t287 * t273 - t304 * t277 + t341;
t342 = -t322 * t181 + t327 * t182;
t172 = m(4) * t243 - t311 * mrSges(4,2) + t297 * mrSges(4,3) - t312 * t292 + t294 * t343 + t342;
t175 = t327 * t181 + t322 * t182;
t252 = -t315 * t271 + t284;
t293 = -t312 * mrSges(4,2) + mrSges(4,3) * t343;
t174 = m(4) * t252 - t297 * mrSges(4,1) + t296 * mrSges(4,2) + (t292 * t323 - t293 * t328) * t347 + t175;
t191 = t325 * t202 + t320 * t203;
t335 = m(6) * t220 - t235 * mrSges(6,1) + t236 * mrSges(6,2) - t269 * t255 + t270 * t256 + t191;
t332 = -m(5) * t229 + t263 * mrSges(5,1) - t264 * mrSges(5,2) + t287 * t276 - t288 * t277 - t335;
t186 = m(4) * t242 + t311 * mrSges(4,1) - t296 * mrSges(4,3) + t312 * t293 - t294 * t344 + t332;
t354 = t186 * t328;
t162 = t172 * t350 - t315 * t174 + t318 * t354;
t159 = m(3) * t274 + qJDD(2) * mrSges(3,1) - t330 * mrSges(3,2) + t162;
t168 = t328 * t172 - t323 * t186;
t167 = m(3) * t275 - t330 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t168;
t156 = -t324 * t159 + t329 * t167;
t355 = pkin(8) * t156;
t161 = t172 * t353 + t318 * t174 + t315 * t354;
t160 = m(3) * t289 + t161;
t151 = t159 * t348 - t316 * t160 + t167 * t349;
t224 = Ifges(7,5) * t254 + Ifges(7,6) * t253 + Ifges(7,3) * t266;
t226 = Ifges(7,1) * t254 + Ifges(7,4) * t253 + Ifges(7,5) * t266;
t195 = -mrSges(7,1) * t208 + mrSges(7,3) * t207 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t234 - t254 * t224 + t266 * t226;
t225 = Ifges(7,4) * t254 + Ifges(7,2) * t253 + Ifges(7,6) * t266;
t196 = mrSges(7,2) * t208 - mrSges(7,3) * t206 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t234 + t253 * t224 - t266 * t225;
t246 = Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t303;
t247 = Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t303;
t176 = mrSges(6,2) * t220 - mrSges(6,3) * t211 + Ifges(6,1) * t236 + Ifges(6,4) * t235 + Ifges(6,5) * t290 - pkin(12) * t191 - t320 * t195 + t325 * t196 + t269 * t246 - t303 * t247;
t248 = Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t303;
t333 = mrSges(7,1) * t206 - mrSges(7,2) * t207 + Ifges(7,5) * t223 + Ifges(7,6) * t222 + Ifges(7,3) * t234 + t254 * t225 - t253 * t226;
t177 = -mrSges(6,1) * t220 + mrSges(6,3) * t212 + Ifges(6,4) * t236 + Ifges(6,2) * t235 + Ifges(6,6) * t290 - pkin(5) * t191 - t270 * t246 + t303 * t248 - t333;
t257 = Ifges(5,5) * t288 + Ifges(5,6) * t287 + Ifges(5,3) * t304;
t259 = Ifges(5,1) * t288 + Ifges(5,4) * t287 + Ifges(5,5) * t304;
t163 = -mrSges(5,1) * t229 + mrSges(5,3) * t219 + Ifges(5,4) * t264 + Ifges(5,2) * t263 + Ifges(5,6) * t291 - pkin(4) * t335 + pkin(11) * t341 + t321 * t176 + t326 * t177 - t288 * t257 + t304 * t259;
t258 = Ifges(5,4) * t288 + Ifges(5,2) * t287 + Ifges(5,6) * t304;
t164 = mrSges(5,2) * t229 - mrSges(5,3) * t218 + Ifges(5,1) * t264 + Ifges(5,4) * t263 + Ifges(5,5) * t291 - pkin(11) * t183 + t326 * t176 - t321 * t177 + t287 * t257 - t304 * t258;
t280 = Ifges(4,3) * t312 + (Ifges(4,5) * t323 + Ifges(4,6) * t328) * t347;
t281 = Ifges(4,6) * t312 + (Ifges(4,4) * t323 + Ifges(4,2) * t328) * t347;
t153 = mrSges(4,2) * t252 - mrSges(4,3) * t242 + Ifges(4,1) * t296 + Ifges(4,4) * t297 + Ifges(4,5) * t311 - pkin(10) * t175 - t322 * t163 + t327 * t164 + t280 * t343 - t312 * t281;
t282 = Ifges(4,5) * t312 + (Ifges(4,1) * t323 + Ifges(4,4) * t328) * t347;
t334 = -mrSges(6,1) * t211 + mrSges(6,2) * t212 - Ifges(6,5) * t236 - Ifges(6,6) * t235 - Ifges(6,3) * t290 - pkin(5) * t336 - pkin(12) * t340 - t325 * t195 - t320 * t196 - t270 * t247 + t269 * t248;
t331 = mrSges(5,1) * t218 - mrSges(5,2) * t219 + Ifges(5,5) * t264 + Ifges(5,6) * t263 + Ifges(5,3) * t291 + pkin(4) * t183 + t288 * t258 - t287 * t259 - t334;
t157 = -mrSges(4,1) * t252 + mrSges(4,3) * t243 + Ifges(4,4) * t296 + Ifges(4,2) * t297 + Ifges(4,6) * t311 - pkin(3) * t175 - t280 * t344 + t312 * t282 - t331;
t338 = pkin(9) * t168 + t153 * t323 + t157 * t328;
t152 = Ifges(4,5) * t296 + Ifges(4,6) * t297 + Ifges(4,3) * t311 + mrSges(4,1) * t242 - mrSges(4,2) * t243 + t322 * t164 + t327 * t163 + pkin(3) * t332 + pkin(10) * t342 + (t281 * t323 - t282 * t328) * t347;
t143 = mrSges(3,1) * t274 - mrSges(3,2) * t275 + Ifges(3,3) * qJDD(2) + pkin(2) * t162 + t318 * t152 + t338 * t315;
t145 = -mrSges(3,1) * t289 + mrSges(3,3) * t275 + t330 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t161 - t315 * t152 + t338 * t318;
t147 = mrSges(3,2) * t289 - mrSges(3,3) * t274 + Ifges(3,5) * qJDD(2) - t330 * Ifges(3,6) + t328 * t153 - t323 * t157 + (-t161 * t315 - t162 * t318) * pkin(9);
t337 = mrSges(2,1) * t305 - mrSges(2,2) * t306 + pkin(1) * t151 + t319 * t143 + t145 * t351 + t147 * t352 + t316 * t355;
t154 = m(2) * t306 + t156;
t150 = t319 * t160 + (t159 * t329 + t167 * t324) * t316;
t148 = m(2) * t305 + t151;
t141 = mrSges(2,2) * t313 - mrSges(2,3) * t305 - t324 * t145 + t329 * t147 + (-t150 * t316 - t151 * t319) * pkin(8);
t140 = -mrSges(2,1) * t313 + mrSges(2,3) * t306 - pkin(1) * t150 - t316 * t143 + (t145 * t329 + t147 * t324 + t355) * t319;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t317 * t141 - t314 * t140 - qJ(1) * (t317 * t148 + t314 * t154), t141, t147, t153, t164, t176, t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t314 * t141 + t317 * t140 + qJ(1) * (-t314 * t148 + t317 * t154), t140, t145, t157, t163, t177, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t337, t337, t143, t152, t331, -t334, t333;];
m_new  = t1;
