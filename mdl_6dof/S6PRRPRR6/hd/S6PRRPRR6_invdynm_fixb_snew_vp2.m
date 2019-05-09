% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:40:06
% EndTime: 2019-05-05 05:41:13
% DurationCPUTime: 63.79s
% Computational Cost: add. (1199136->352), mult. (2602187->468), div. (0->0), fcn. (2109010->16), ass. (0->156)
t315 = sin(pkin(12));
t319 = cos(pkin(12));
t305 = t315 * g(1) - t319 * g(2);
t306 = -t319 * g(1) - t315 * g(2);
t313 = -g(3) + qJDD(1);
t325 = sin(qJ(2));
t321 = cos(pkin(6));
t329 = cos(qJ(2));
t348 = t321 * t329;
t317 = sin(pkin(6));
t351 = t317 * t329;
t271 = t305 * t348 - t325 * t306 + t313 * t351;
t316 = sin(pkin(7));
t330 = qJD(2) ^ 2;
t269 = t330 * t316 * pkin(9) + qJDD(2) * pkin(2) + t271;
t349 = t321 * t325;
t352 = t317 * t325;
t272 = t305 * t349 + t329 * t306 + t313 * t352;
t345 = qJDD(2) * t316;
t270 = -t330 * pkin(2) + pkin(9) * t345 + t272;
t290 = -t317 * t305 + t321 * t313;
t320 = cos(pkin(7));
t324 = sin(qJ(3));
t328 = cos(qJ(3));
t236 = -t324 * t270 + (t269 * t320 + t290 * t316) * t328;
t350 = t320 * t324;
t353 = t316 * t324;
t237 = t269 * t350 + t328 * t270 + t290 * t353;
t312 = t320 * qJD(2) + qJD(3);
t347 = qJD(2) * t316;
t344 = t324 * t347;
t293 = t312 * mrSges(4,1) - mrSges(4,3) * t344;
t296 = (-mrSges(4,1) * t328 + mrSges(4,2) * t324) * t347;
t298 = -qJD(3) * t344 + t328 * t345;
t311 = t320 * qJDD(2) + qJDD(3);
t295 = (-pkin(3) * t328 - qJ(4) * t324) * t347;
t310 = t312 ^ 2;
t346 = qJD(2) * t328;
t343 = t316 * t346;
t230 = -t310 * pkin(3) + t311 * qJ(4) + t295 * t343 + t237;
t286 = t320 * t290;
t297 = (qJD(3) * t346 + qJDD(2) * t324) * t316;
t234 = -t298 * pkin(3) - t297 * qJ(4) + t286 + (-t269 + (pkin(3) * t324 - qJ(4) * t328) * t312 * qJD(2)) * t316;
t314 = sin(pkin(13));
t318 = cos(pkin(13));
t289 = t314 * t312 + t318 * t344;
t218 = -0.2e1 * qJD(4) * t289 - t314 * t230 + t318 * t234;
t277 = t318 * t297 + t314 * t311;
t288 = t318 * t312 - t314 * t344;
t215 = (-t288 * t343 - t277) * pkin(10) + (t288 * t289 - t298) * pkin(4) + t218;
t219 = 0.2e1 * qJD(4) * t288 + t318 * t230 + t314 * t234;
t276 = -t314 * t297 + t318 * t311;
t278 = -pkin(4) * t343 - t289 * pkin(10);
t287 = t288 ^ 2;
t217 = -t287 * pkin(4) + t276 * pkin(10) + t278 * t343 + t219;
t323 = sin(qJ(5));
t327 = cos(qJ(5));
t212 = t323 * t215 + t327 * t217;
t265 = t323 * t288 + t327 * t289;
t242 = -t265 * qJD(5) + t327 * t276 - t323 * t277;
t264 = t327 * t288 - t323 * t289;
t250 = -t264 * mrSges(6,1) + t265 * mrSges(6,2);
t304 = qJD(5) - t343;
t256 = t304 * mrSges(6,1) - t265 * mrSges(6,3);
t292 = qJDD(5) - t298;
t251 = -t264 * pkin(5) - t265 * pkin(11);
t303 = t304 ^ 2;
t209 = -t303 * pkin(5) + t292 * pkin(11) + t264 * t251 + t212;
t229 = -t311 * pkin(3) - t310 * qJ(4) + t295 * t344 + qJDD(4) - t236;
t220 = -t276 * pkin(4) - t287 * pkin(10) + t289 * t278 + t229;
t243 = t264 * qJD(5) + t323 * t276 + t327 * t277;
t213 = (-t264 * t304 - t243) * pkin(11) + (t265 * t304 - t242) * pkin(5) + t220;
t322 = sin(qJ(6));
t326 = cos(qJ(6));
t206 = -t322 * t209 + t326 * t213;
t253 = -t322 * t265 + t326 * t304;
t223 = t253 * qJD(6) + t326 * t243 + t322 * t292;
t254 = t326 * t265 + t322 * t304;
t235 = -t253 * mrSges(7,1) + t254 * mrSges(7,2);
t241 = qJDD(6) - t242;
t263 = qJD(6) - t264;
t244 = -t263 * mrSges(7,2) + t253 * mrSges(7,3);
t202 = m(7) * t206 + t241 * mrSges(7,1) - t223 * mrSges(7,3) - t254 * t235 + t263 * t244;
t207 = t326 * t209 + t322 * t213;
t222 = -t254 * qJD(6) - t322 * t243 + t326 * t292;
t245 = t263 * mrSges(7,1) - t254 * mrSges(7,3);
t203 = m(7) * t207 - t241 * mrSges(7,2) + t222 * mrSges(7,3) + t253 * t235 - t263 * t245;
t340 = -t322 * t202 + t326 * t203;
t189 = m(6) * t212 - t292 * mrSges(6,2) + t242 * mrSges(6,3) + t264 * t250 - t304 * t256 + t340;
t211 = t327 * t215 - t323 * t217;
t255 = -t304 * mrSges(6,2) + t264 * mrSges(6,3);
t208 = -t292 * pkin(5) - t303 * pkin(11) + t265 * t251 - t211;
t336 = -m(7) * t208 + t222 * mrSges(7,1) - t223 * mrSges(7,2) + t253 * t244 - t254 * t245;
t198 = m(6) * t211 + t292 * mrSges(6,1) - t243 * mrSges(6,3) - t265 * t250 + t304 * t255 + t336;
t183 = t323 * t189 + t327 * t198;
t268 = -t288 * mrSges(5,1) + t289 * mrSges(5,2);
t274 = mrSges(5,2) * t343 + t288 * mrSges(5,3);
t181 = m(5) * t218 - t298 * mrSges(5,1) - t277 * mrSges(5,3) - t289 * t268 - t274 * t343 + t183;
t275 = -mrSges(5,1) * t343 - t289 * mrSges(5,3);
t341 = t327 * t189 - t323 * t198;
t182 = m(5) * t219 + t298 * mrSges(5,2) + t276 * mrSges(5,3) + t288 * t268 + t275 * t343 + t341;
t342 = -t314 * t181 + t318 * t182;
t172 = m(4) * t237 - t311 * mrSges(4,2) + t298 * mrSges(4,3) - t312 * t293 + t296 * t343 + t342;
t175 = t318 * t181 + t314 * t182;
t252 = -t316 * t269 + t286;
t294 = -t312 * mrSges(4,2) + mrSges(4,3) * t343;
t174 = m(4) * t252 - t298 * mrSges(4,1) + t297 * mrSges(4,2) + (t293 * t324 - t294 * t328) * t347 + t175;
t191 = t326 * t202 + t322 * t203;
t335 = m(6) * t220 - t242 * mrSges(6,1) + t243 * mrSges(6,2) - t264 * t255 + t265 * t256 + t191;
t332 = -m(5) * t229 + t276 * mrSges(5,1) - t277 * mrSges(5,2) + t288 * t274 - t289 * t275 - t335;
t186 = m(4) * t236 + t311 * mrSges(4,1) - t297 * mrSges(4,3) + t312 * t294 - t296 * t344 + t332;
t354 = t186 * t328;
t162 = t172 * t350 - t316 * t174 + t320 * t354;
t159 = m(3) * t271 + qJDD(2) * mrSges(3,1) - t330 * mrSges(3,2) + t162;
t168 = t328 * t172 - t324 * t186;
t167 = m(3) * t272 - t330 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t168;
t156 = -t325 * t159 + t329 * t167;
t355 = pkin(8) * t156;
t161 = t172 * t353 + t320 * t174 + t316 * t354;
t160 = m(3) * t290 + t161;
t151 = t159 * t348 - t317 * t160 + t167 * t349;
t224 = Ifges(7,5) * t254 + Ifges(7,6) * t253 + Ifges(7,3) * t263;
t226 = Ifges(7,1) * t254 + Ifges(7,4) * t253 + Ifges(7,5) * t263;
t195 = -mrSges(7,1) * t208 + mrSges(7,3) * t207 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t241 - t254 * t224 + t263 * t226;
t225 = Ifges(7,4) * t254 + Ifges(7,2) * t253 + Ifges(7,6) * t263;
t196 = mrSges(7,2) * t208 - mrSges(7,3) * t206 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t241 + t253 * t224 - t263 * t225;
t246 = Ifges(6,5) * t265 + Ifges(6,6) * t264 + Ifges(6,3) * t304;
t247 = Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * t304;
t176 = mrSges(6,2) * t220 - mrSges(6,3) * t211 + Ifges(6,1) * t243 + Ifges(6,4) * t242 + Ifges(6,5) * t292 - pkin(11) * t191 - t322 * t195 + t326 * t196 + t264 * t246 - t304 * t247;
t248 = Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * t304;
t333 = mrSges(7,1) * t206 - mrSges(7,2) * t207 + Ifges(7,5) * t223 + Ifges(7,6) * t222 + Ifges(7,3) * t241 + t254 * t225 - t253 * t226;
t177 = -mrSges(6,1) * t220 + mrSges(6,3) * t212 + Ifges(6,4) * t243 + Ifges(6,2) * t242 + Ifges(6,6) * t292 - pkin(5) * t191 - t265 * t246 + t304 * t248 - t333;
t257 = Ifges(5,5) * t289 + Ifges(5,6) * t288 - Ifges(5,3) * t343;
t259 = Ifges(5,1) * t289 + Ifges(5,4) * t288 - Ifges(5,5) * t343;
t163 = -mrSges(5,1) * t229 + mrSges(5,3) * t219 + Ifges(5,4) * t277 + Ifges(5,2) * t276 - Ifges(5,6) * t298 - pkin(4) * t335 + pkin(10) * t341 + t323 * t176 + t327 * t177 - t289 * t257 - t259 * t343;
t258 = Ifges(5,4) * t289 + Ifges(5,2) * t288 - Ifges(5,6) * t343;
t164 = mrSges(5,2) * t229 - mrSges(5,3) * t218 + Ifges(5,1) * t277 + Ifges(5,4) * t276 - Ifges(5,5) * t298 - pkin(10) * t183 + t327 * t176 - t323 * t177 + t288 * t257 + t258 * t343;
t280 = Ifges(4,3) * t312 + (Ifges(4,5) * t324 + Ifges(4,6) * t328) * t347;
t281 = Ifges(4,6) * t312 + (Ifges(4,4) * t324 + Ifges(4,2) * t328) * t347;
t153 = mrSges(4,2) * t252 - mrSges(4,3) * t236 + Ifges(4,1) * t297 + Ifges(4,4) * t298 + Ifges(4,5) * t311 - qJ(4) * t175 - t314 * t163 + t318 * t164 + t280 * t343 - t312 * t281;
t282 = Ifges(4,5) * t312 + (Ifges(4,1) * t324 + Ifges(4,4) * t328) * t347;
t334 = -mrSges(6,1) * t211 + mrSges(6,2) * t212 - Ifges(6,5) * t243 - Ifges(6,6) * t242 - Ifges(6,3) * t292 - pkin(5) * t336 - pkin(11) * t340 - t326 * t195 - t322 * t196 - t265 * t247 + t264 * t248;
t331 = -mrSges(5,1) * t218 + mrSges(5,2) * t219 - Ifges(5,5) * t277 - Ifges(5,6) * t276 - pkin(4) * t183 - t289 * t258 + t288 * t259 + t334;
t157 = t331 - pkin(3) * t175 - t280 * t344 + (Ifges(4,2) + Ifges(5,3)) * t298 + mrSges(4,3) * t237 - mrSges(4,1) * t252 + Ifges(4,4) * t297 + Ifges(4,6) * t311 + t312 * t282;
t338 = pkin(9) * t168 + t153 * t324 + t157 * t328;
t152 = Ifges(4,5) * t297 + Ifges(4,6) * t298 + Ifges(4,3) * t311 + mrSges(4,1) * t236 - mrSges(4,2) * t237 + t314 * t164 + t318 * t163 + pkin(3) * t332 + qJ(4) * t342 + (t281 * t324 - t282 * t328) * t347;
t143 = mrSges(3,1) * t271 - mrSges(3,2) * t272 + Ifges(3,3) * qJDD(2) + pkin(2) * t162 + t320 * t152 + t338 * t316;
t145 = -mrSges(3,1) * t290 + mrSges(3,3) * t272 + t330 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t161 - t316 * t152 + t338 * t320;
t147 = mrSges(3,2) * t290 - mrSges(3,3) * t271 + Ifges(3,5) * qJDD(2) - t330 * Ifges(3,6) + t328 * t153 - t324 * t157 + (-t161 * t316 - t162 * t320) * pkin(9);
t337 = mrSges(2,1) * t305 - mrSges(2,2) * t306 + pkin(1) * t151 + t321 * t143 + t145 * t351 + t147 * t352 + t317 * t355;
t154 = m(2) * t306 + t156;
t150 = t321 * t160 + (t159 * t329 + t167 * t325) * t317;
t148 = m(2) * t305 + t151;
t141 = mrSges(2,2) * t313 - mrSges(2,3) * t305 - t325 * t145 + t329 * t147 + (-t150 * t317 - t151 * t321) * pkin(8);
t140 = -mrSges(2,1) * t313 + mrSges(2,3) * t306 - pkin(1) * t150 - t317 * t143 + (t145 * t329 + t147 * t325 + t355) * t321;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t319 * t141 - t315 * t140 - qJ(1) * (t319 * t148 + t315 * t154), t141, t147, t153, t164, t176, t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t315 * t141 + t319 * t140 + qJ(1) * (-t315 * t148 + t319 * t154), t140, t145, t157, t163, t177, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t337, t337, t143, t152, -Ifges(5,3) * t298 - t331, -t334, t333;];
m_new  = t1;
