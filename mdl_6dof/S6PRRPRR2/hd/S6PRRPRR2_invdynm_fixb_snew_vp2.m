% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:28:28
% EndTime: 2019-05-05 04:29:08
% DurationCPUTime: 29.55s
% Computational Cost: add. (539235->340), mult. (1167174->440), div. (0->0), fcn. (850034->14), ass. (0->144)
t308 = sin(pkin(11));
t311 = cos(pkin(11));
t297 = g(1) * t308 - g(2) * t311;
t298 = -g(1) * t311 - g(2) * t308;
t306 = -g(3) + qJDD(1);
t320 = cos(qJ(2));
t312 = cos(pkin(6));
t316 = sin(qJ(2));
t345 = t312 * t316;
t309 = sin(pkin(6));
t346 = t309 * t316;
t257 = t297 * t345 + t320 * t298 + t306 * t346;
t322 = qJD(2) ^ 2;
t252 = -pkin(2) * t322 + qJDD(2) * pkin(8) + t257;
t276 = -t297 * t309 + t306 * t312;
t315 = sin(qJ(3));
t319 = cos(qJ(3));
t233 = -t315 * t252 + t319 * t276;
t341 = qJD(2) * qJD(3);
t340 = t319 * t341;
t294 = qJDD(2) * t315 + t340;
t228 = (-t294 + t340) * qJ(4) + (t315 * t319 * t322 + qJDD(3)) * pkin(3) + t233;
t234 = t319 * t252 + t315 * t276;
t295 = qJDD(2) * t319 - t315 * t341;
t344 = qJD(2) * t315;
t299 = qJD(3) * pkin(3) - qJ(4) * t344;
t305 = t319 ^ 2;
t229 = -pkin(3) * t305 * t322 + qJ(4) * t295 - qJD(3) * t299 + t234;
t307 = sin(pkin(12));
t310 = cos(pkin(12));
t343 = qJD(2) * t319;
t281 = -t307 * t344 + t310 * t343;
t210 = 0.2e1 * qJD(4) * t281 + t307 * t228 + t310 * t229;
t282 = (t307 * t319 + t310 * t315) * qJD(2);
t259 = -mrSges(5,1) * t281 + mrSges(5,2) * t282;
t267 = -t307 * t294 + t295 * t310;
t275 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t282;
t261 = -pkin(4) * t281 - pkin(9) * t282;
t321 = qJD(3) ^ 2;
t207 = -pkin(4) * t321 + qJDD(3) * pkin(9) + t261 * t281 + t210;
t256 = -t316 * t298 + (t297 * t312 + t306 * t309) * t320;
t329 = -qJDD(2) * pkin(2) - t256;
t230 = -t295 * pkin(3) + qJDD(4) + t299 * t344 + (-qJ(4) * t305 - pkin(8)) * t322 + t329;
t268 = t294 * t310 + t295 * t307;
t219 = (-qJD(3) * t281 - t268) * pkin(9) + (qJD(3) * t282 - t267) * pkin(4) + t230;
t314 = sin(qJ(5));
t318 = cos(qJ(5));
t202 = -t314 * t207 + t318 * t219;
t270 = qJD(3) * t318 - t282 * t314;
t241 = qJD(5) * t270 + qJDD(3) * t314 + t268 * t318;
t266 = qJDD(5) - t267;
t271 = qJD(3) * t314 + t282 * t318;
t280 = qJD(5) - t281;
t200 = (t270 * t280 - t241) * pkin(10) + (t270 * t271 + t266) * pkin(5) + t202;
t203 = t318 * t207 + t314 * t219;
t240 = -qJD(5) * t271 + qJDD(3) * t318 - t268 * t314;
t250 = pkin(5) * t280 - pkin(10) * t271;
t269 = t270 ^ 2;
t201 = -pkin(5) * t269 + pkin(10) * t240 - t250 * t280 + t203;
t313 = sin(qJ(6));
t317 = cos(qJ(6));
t198 = t200 * t317 - t201 * t313;
t242 = t270 * t317 - t271 * t313;
t215 = qJD(6) * t242 + t240 * t313 + t241 * t317;
t243 = t270 * t313 + t271 * t317;
t224 = -mrSges(7,1) * t242 + mrSges(7,2) * t243;
t277 = qJD(6) + t280;
t231 = -mrSges(7,2) * t277 + mrSges(7,3) * t242;
t262 = qJDD(6) + t266;
t193 = m(7) * t198 + mrSges(7,1) * t262 - mrSges(7,3) * t215 - t224 * t243 + t231 * t277;
t199 = t200 * t313 + t201 * t317;
t214 = -qJD(6) * t243 + t240 * t317 - t241 * t313;
t232 = mrSges(7,1) * t277 - mrSges(7,3) * t243;
t194 = m(7) * t199 - mrSges(7,2) * t262 + mrSges(7,3) * t214 + t224 * t242 - t232 * t277;
t185 = t317 * t193 + t313 * t194;
t245 = -mrSges(6,1) * t270 + mrSges(6,2) * t271;
t248 = -mrSges(6,2) * t280 + mrSges(6,3) * t270;
t183 = m(6) * t202 + mrSges(6,1) * t266 - mrSges(6,3) * t241 - t245 * t271 + t248 * t280 + t185;
t249 = mrSges(6,1) * t280 - mrSges(6,3) * t271;
t336 = -t193 * t313 + t317 * t194;
t184 = m(6) * t203 - mrSges(6,2) * t266 + mrSges(6,3) * t240 + t245 * t270 - t249 * t280 + t336;
t337 = -t183 * t314 + t318 * t184;
t176 = m(5) * t210 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t267 - qJD(3) * t275 + t259 * t281 + t337;
t209 = -0.2e1 * qJD(4) * t282 + t228 * t310 - t307 * t229;
t274 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t281;
t206 = -qJDD(3) * pkin(4) - pkin(9) * t321 + t282 * t261 - t209;
t204 = -pkin(5) * t240 - pkin(10) * t269 + t250 * t271 + t206;
t331 = m(7) * t204 - t214 * mrSges(7,1) + mrSges(7,2) * t215 - t242 * t231 + t232 * t243;
t326 = -m(6) * t206 + t240 * mrSges(6,1) - mrSges(6,2) * t241 + t270 * t248 - t249 * t271 - t331;
t189 = m(5) * t209 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t268 + qJD(3) * t274 - t259 * t282 + t326;
t167 = t307 * t176 + t310 * t189;
t293 = (-mrSges(4,1) * t319 + mrSges(4,2) * t315) * qJD(2);
t301 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t343;
t165 = m(4) * t233 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t294 + qJD(3) * t301 - t293 * t344 + t167;
t300 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t344;
t338 = t310 * t176 - t189 * t307;
t166 = m(4) * t234 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t295 - qJD(3) * t300 + t293 * t343 + t338;
t160 = t319 * t165 + t315 * t166;
t285 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t315 + Ifges(4,2) * t319) * qJD(2);
t286 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t315 + Ifges(4,4) * t319) * qJD(2);
t220 = Ifges(7,5) * t243 + Ifges(7,6) * t242 + Ifges(7,3) * t277;
t222 = Ifges(7,1) * t243 + Ifges(7,4) * t242 + Ifges(7,5) * t277;
t186 = -mrSges(7,1) * t204 + mrSges(7,3) * t199 + Ifges(7,4) * t215 + Ifges(7,2) * t214 + Ifges(7,6) * t262 - t220 * t243 + t222 * t277;
t221 = Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t277;
t187 = mrSges(7,2) * t204 - mrSges(7,3) * t198 + Ifges(7,1) * t215 + Ifges(7,4) * t214 + Ifges(7,5) * t262 + t220 * t242 - t221 * t277;
t235 = Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t280;
t237 = Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t280;
t169 = -mrSges(6,1) * t206 + mrSges(6,3) * t203 + Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t266 - pkin(5) * t331 + pkin(10) * t336 + t317 * t186 + t313 * t187 - t271 * t235 + t280 * t237;
t236 = Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t280;
t171 = mrSges(6,2) * t206 - mrSges(6,3) * t202 + Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t266 - pkin(10) * t185 - t186 * t313 + t187 * t317 + t235 * t270 - t236 * t280;
t254 = Ifges(5,4) * t282 + Ifges(5,2) * t281 + Ifges(5,6) * qJD(3);
t255 = Ifges(5,1) * t282 + Ifges(5,4) * t281 + Ifges(5,5) * qJD(3);
t327 = -mrSges(5,1) * t209 + mrSges(5,2) * t210 - Ifges(5,5) * t268 - Ifges(5,6) * t267 - Ifges(5,3) * qJDD(3) - pkin(4) * t326 - pkin(9) * t337 - t318 * t169 - t314 * t171 - t282 * t254 + t281 * t255;
t350 = mrSges(4,1) * t233 - mrSges(4,2) * t234 + Ifges(4,5) * t294 + Ifges(4,6) * t295 + Ifges(4,3) * qJDD(3) + pkin(3) * t167 + (t285 * t315 - t286 * t319) * qJD(2) - t327;
t146 = -mrSges(3,1) * t276 + mrSges(3,3) * t257 + t322 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t160 - t350;
t339 = -t165 * t315 + t319 * t166;
t158 = m(3) * t257 - mrSges(3,1) * t322 - qJDD(2) * mrSges(3,2) + t339;
t251 = -t322 * pkin(8) + t329;
t178 = t318 * t183 + t314 * t184;
t328 = m(5) * t230 - t267 * mrSges(5,1) + mrSges(5,2) * t268 - t281 * t274 + t275 * t282 + t178;
t325 = -m(4) * t251 + t295 * mrSges(4,1) - mrSges(4,2) * t294 - t300 * t344 + t301 * t343 - t328;
t173 = m(3) * t256 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t322 + t325;
t154 = t320 * t158 - t173 * t316;
t351 = pkin(7) * t154 + t146 * t320;
t347 = t173 * t320;
t159 = m(3) * t276 + t160;
t151 = t158 * t345 - t159 * t309 + t312 * t347;
t253 = Ifges(5,5) * t282 + Ifges(5,6) * t281 + Ifges(5,3) * qJD(3);
t155 = mrSges(5,2) * t230 - mrSges(5,3) * t209 + Ifges(5,1) * t268 + Ifges(5,4) * t267 + Ifges(5,5) * qJDD(3) - pkin(9) * t178 - qJD(3) * t254 - t169 * t314 + t171 * t318 + t253 * t281;
t330 = -mrSges(7,1) * t198 + mrSges(7,2) * t199 - Ifges(7,5) * t215 - Ifges(7,6) * t214 - Ifges(7,3) * t262 - t243 * t221 + t242 * t222;
t323 = mrSges(6,1) * t202 - mrSges(6,2) * t203 + Ifges(6,5) * t241 + Ifges(6,6) * t240 + Ifges(6,3) * t266 + pkin(5) * t185 + t271 * t236 - t270 * t237 - t330;
t161 = -mrSges(5,1) * t230 + mrSges(5,3) * t210 + Ifges(5,4) * t268 + Ifges(5,2) * t267 + Ifges(5,6) * qJDD(3) - pkin(4) * t178 + qJD(3) * t255 - t282 * t253 - t323;
t284 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t315 + Ifges(4,6) * t319) * qJD(2);
t144 = -mrSges(4,1) * t251 + mrSges(4,3) * t234 + Ifges(4,4) * t294 + Ifges(4,2) * t295 + Ifges(4,6) * qJDD(3) - pkin(3) * t328 + qJ(4) * t338 + qJD(3) * t286 + t307 * t155 + t310 * t161 - t284 * t344;
t147 = mrSges(4,2) * t251 - mrSges(4,3) * t233 + Ifges(4,1) * t294 + Ifges(4,4) * t295 + Ifges(4,5) * qJDD(3) - qJ(4) * t167 - qJD(3) * t285 + t155 * t310 - t161 * t307 + t284 * t343;
t141 = mrSges(3,1) * t256 - mrSges(3,2) * t257 + Ifges(3,3) * qJDD(2) + pkin(2) * t325 + pkin(8) * t339 + t319 * t144 + t315 * t147;
t143 = mrSges(3,2) * t276 - mrSges(3,3) * t256 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t322 - pkin(8) * t160 - t144 * t315 + t147 * t319;
t332 = mrSges(2,1) * t297 - mrSges(2,2) * t298 + pkin(1) * t151 + t312 * t141 + t143 * t346 + t351 * t309;
t152 = m(2) * t298 + t154;
t150 = t312 * t159 + (t158 * t316 + t347) * t309;
t148 = m(2) * t297 + t151;
t139 = mrSges(2,2) * t306 - mrSges(2,3) * t297 + t320 * t143 - t316 * t146 + (-t150 * t309 - t151 * t312) * pkin(7);
t138 = -mrSges(2,1) * t306 + mrSges(2,3) * t298 - pkin(1) * t150 - t309 * t141 + (t143 * t316 + t351) * t312;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t311 * t139 - t308 * t138 - qJ(1) * (t148 * t311 + t152 * t308), t139, t143, t147, t155, t171, t187; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t308 * t139 + t311 * t138 + qJ(1) * (-t148 * t308 + t152 * t311), t138, t146, t144, t161, t169, t186; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t332, t332, t141, t350, -t327, t323, -t330;];
m_new  = t1;
