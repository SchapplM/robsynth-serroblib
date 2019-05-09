% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPPR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-05-05 02:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:44:12
% EndTime: 2019-05-05 02:44:30
% DurationCPUTime: 11.22s
% Computational Cost: add. (186772->344), mult. (405461->422), div. (0->0), fcn. (276576->12), ass. (0->140)
t309 = sin(pkin(10));
t311 = cos(pkin(10));
t296 = g(1) * t309 - g(2) * t311;
t297 = -g(1) * t311 - g(2) * t309;
t307 = -g(3) + qJDD(1);
t318 = cos(qJ(2));
t312 = cos(pkin(6));
t315 = sin(qJ(2));
t350 = t312 * t315;
t310 = sin(pkin(6));
t351 = t310 * t315;
t244 = t296 * t350 + t297 * t318 + t307 * t351;
t320 = qJD(2) ^ 2;
t236 = -pkin(2) * t320 + qJDD(2) * pkin(8) + t244;
t273 = -t296 * t310 + t307 * t312;
t314 = sin(qJ(3));
t317 = cos(qJ(3));
t216 = -t314 * t236 + t273 * t317;
t342 = qJD(2) * qJD(3);
t341 = t317 * t342;
t293 = qJDD(2) * t314 + t341;
t212 = (-t293 + t341) * qJ(4) + (t314 * t317 * t320 + qJDD(3)) * pkin(3) + t216;
t217 = t236 * t317 + t273 * t314;
t294 = qJDD(2) * t317 - t314 * t342;
t346 = qJD(2) * t314;
t298 = qJD(3) * pkin(3) - qJ(4) * t346;
t306 = t317 ^ 2;
t213 = -pkin(3) * t306 * t320 + qJ(4) * t294 - qJD(3) * t298 + t217;
t308 = sin(pkin(11));
t354 = cos(pkin(11));
t281 = (t308 * t317 + t314 * t354) * qJD(2);
t361 = -2 * qJD(4);
t205 = t212 * t354 - t213 * t308 + t281 * t361;
t345 = qJD(2) * t317;
t280 = t308 * t346 - t345 * t354;
t276 = t280 * t361;
t349 = t212 * t308 + t213 * t354;
t206 = t276 + t349;
t240 = Ifges(5,4) * t281 - Ifges(5,2) * t280 + Ifges(5,6) * qJD(3);
t250 = -mrSges(6,2) * t280 - mrSges(6,3) * t281;
t261 = t293 * t308 - t294 * t354;
t262 = t293 * t354 + t294 * t308;
t270 = mrSges(6,1) * t280 - qJD(3) * mrSges(6,3);
t248 = pkin(4) * t280 - qJ(5) * t281;
t319 = qJD(3) ^ 2;
t202 = -qJDD(3) * pkin(4) - t319 * qJ(5) + t248 * t281 + qJDD(5) - t205;
t344 = qJD(3) * t280;
t197 = (t280 * t281 - qJDD(3)) * pkin(9) + (t262 + t344) * pkin(5) + t202;
t272 = pkin(5) * t281 - qJD(3) * pkin(9);
t279 = t280 ^ 2;
t243 = -t315 * t297 + (t296 * t312 + t307 * t310) * t318;
t330 = -qJDD(2) * pkin(2) - t243;
t215 = -t294 * pkin(3) + qJDD(4) + t298 * t346 + (-qJ(4) * t306 - pkin(8)) * t320 + t330;
t357 = -2 * qJD(5);
t322 = (-t262 + t344) * qJ(5) + t215 + (pkin(4) * qJD(3) + t357) * t281;
t203 = -t279 * pkin(5) - t281 * t272 + (pkin(4) + pkin(9)) * t261 + t322;
t313 = sin(qJ(6));
t316 = cos(qJ(6));
t194 = t197 * t316 - t203 * t313;
t263 = -qJD(3) * t313 + t280 * t316;
t226 = qJD(6) * t263 + qJDD(3) * t316 + t261 * t313;
t264 = qJD(3) * t316 + t280 * t313;
t229 = -mrSges(7,1) * t263 + mrSges(7,2) * t264;
t278 = qJD(6) + t281;
t233 = -mrSges(7,2) * t278 + mrSges(7,3) * t263;
t260 = qJDD(6) + t262;
t191 = m(7) * t194 + mrSges(7,1) * t260 - mrSges(7,3) * t226 - t229 * t264 + t233 * t278;
t195 = t197 * t313 + t203 * t316;
t225 = -qJD(6) * t264 - qJDD(3) * t313 + t261 * t316;
t234 = mrSges(7,1) * t278 - mrSges(7,3) * t264;
t192 = m(7) * t195 - mrSges(7,2) * t260 + mrSges(7,3) * t225 + t229 * t263 - t234 * t278;
t179 = t316 * t191 + t313 * t192;
t334 = t319 * pkin(4) - qJDD(3) * qJ(5) - t349;
t199 = -t261 * pkin(5) - t279 * pkin(9) - t280 * t248 + t276 + ((2 * qJD(5)) + t272) * qJD(3) - t334;
t218 = Ifges(7,5) * t264 + Ifges(7,6) * t263 + Ifges(7,3) * t278;
t220 = Ifges(7,1) * t264 + Ifges(7,4) * t263 + Ifges(7,5) * t278;
t184 = -mrSges(7,1) * t199 + mrSges(7,3) * t195 + Ifges(7,4) * t226 + Ifges(7,2) * t225 + Ifges(7,6) * t260 - t218 * t264 + t220 * t278;
t219 = Ifges(7,4) * t264 + Ifges(7,2) * t263 + Ifges(7,6) * t278;
t185 = mrSges(7,2) * t199 - mrSges(7,3) * t194 + Ifges(7,1) * t226 + Ifges(7,4) * t225 + Ifges(7,5) * t260 + t218 * t263 - t219 * t278;
t200 = qJD(3) * t357 + ((2 * qJD(4)) + t248) * t280 + t334;
t237 = Ifges(6,5) * qJD(3) - Ifges(6,6) * t281 + Ifges(6,3) * t280;
t328 = -mrSges(6,2) * t202 + mrSges(6,3) * t200 - Ifges(6,1) * qJDD(3) + Ifges(6,4) * t262 - Ifges(6,5) * t261 + pkin(9) * t179 + t313 * t184 - t185 * t316 + t237 * t281;
t196 = -m(7) * t199 + t225 * mrSges(7,1) - mrSges(7,2) * t226 + t263 * t233 - t234 * t264;
t271 = mrSges(6,1) * t281 + qJD(3) * mrSges(6,2);
t329 = -m(6) * t200 + qJDD(3) * mrSges(6,3) + qJD(3) * t271 - t196;
t332 = -m(6) * t202 - mrSges(6,1) * t262 - t250 * t281 - t179;
t239 = Ifges(6,4) * qJD(3) - Ifges(6,2) * t281 + Ifges(6,6) * t280;
t347 = Ifges(5,1) * t281 - Ifges(5,4) * t280 + Ifges(5,5) * qJD(3) - t239;
t362 = -mrSges(5,2) * t206 + pkin(4) * (-qJDD(3) * mrSges(6,2) - qJD(3) * t270 + t332) + qJ(5) * (-t261 * mrSges(6,1) - t280 * t250 + t329) + mrSges(5,1) * t205 + t281 * t240 - Ifges(5,6) * t261 + Ifges(5,5) * t262 + Ifges(5,3) * qJDD(3) - t328 + t347 * t280;
t249 = mrSges(5,1) * t280 + mrSges(5,2) * t281;
t268 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t280;
t173 = m(5) * t205 - t262 * mrSges(5,3) - t281 * t249 + (mrSges(5,1) - mrSges(6,2)) * qJDD(3) + (t268 - t270) * qJD(3) + t332;
t269 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t281;
t186 = m(5) * t206 - qJDD(3) * mrSges(5,2) - qJD(3) * t269 + (-t249 - t250) * t280 + (-mrSges(5,3) - mrSges(6,1)) * t261 + t329;
t171 = t173 * t354 + t186 * t308;
t284 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t314 + Ifges(4,2) * t317) * qJD(2);
t285 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t314 + Ifges(4,4) * t317) * qJD(2);
t360 = mrSges(4,1) * t216 - mrSges(4,2) * t217 + Ifges(4,5) * t293 + Ifges(4,6) * t294 + Ifges(4,3) * qJDD(3) + pkin(3) * t171 + (t284 * t314 - t285 * t317) * qJD(2) + t362;
t292 = (-mrSges(4,1) * t317 + mrSges(4,2) * t314) * qJD(2);
t300 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t345;
t169 = m(4) * t216 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t293 + qJD(3) * t300 - t292 * t346 + t171;
t299 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t346;
t338 = -t173 * t308 + t186 * t354;
t170 = m(4) * t217 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t294 - qJD(3) * t299 + t292 * t345 + t338;
t164 = t169 * t317 + t170 * t314;
t150 = -mrSges(3,1) * t273 + mrSges(3,3) * t244 + t320 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t164 - t360;
t339 = -t169 * t314 + t170 * t317;
t162 = m(3) * t244 - mrSges(3,1) * t320 - qJDD(2) * mrSges(3,2) + t339;
t235 = -t320 * pkin(8) + t330;
t180 = -t191 * t313 + t192 * t316;
t208 = t261 * pkin(4) + t322;
t178 = m(6) * t208 - mrSges(6,2) * t261 - mrSges(6,3) * t262 - t270 * t280 - t271 * t281 + t180;
t326 = m(5) * t215 + mrSges(5,1) * t261 + mrSges(5,2) * t262 + t268 * t280 + t269 * t281 + t178;
t323 = -m(4) * t235 + mrSges(4,1) * t294 - mrSges(4,2) * t293 - t299 * t346 + t300 * t345 - t326;
t175 = m(3) * t243 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t320 + t323;
t158 = t162 * t318 - t175 * t315;
t359 = pkin(7) * t158 + t150 * t318;
t355 = Ifges(5,4) + Ifges(6,6);
t352 = t175 * t318;
t241 = Ifges(6,1) * qJD(3) - Ifges(6,4) * t281 + Ifges(6,5) * t280;
t348 = -Ifges(5,5) * t281 + Ifges(5,6) * t280 - Ifges(5,3) * qJD(3) - t241;
t163 = m(3) * t273 + t164;
t155 = t162 * t350 - t163 * t310 + t312 * t352;
t327 = -mrSges(6,1) * t200 + mrSges(6,2) * t208 - pkin(5) * t196 - pkin(9) * t180 - t316 * t184 - t313 * t185;
t159 = -mrSges(5,1) * t215 + mrSges(5,3) * t206 - pkin(4) * t178 + t348 * t281 + t355 * t262 + (-Ifges(5,2) - Ifges(6,3)) * t261 + (Ifges(5,6) - Ifges(6,5)) * qJDD(3) + t347 * qJD(3) + t327;
t331 = mrSges(7,1) * t194 - mrSges(7,2) * t195 + Ifges(7,5) * t226 + Ifges(7,6) * t225 + Ifges(7,3) * t260 + t219 * t264 - t263 * t220;
t325 = mrSges(6,1) * t202 - mrSges(6,3) * t208 + pkin(5) * t179 + t331;
t165 = t325 + t348 * t280 + (Ifges(5,1) + Ifges(6,2)) * t262 - t355 * t261 + (Ifges(5,5) - Ifges(6,4)) * qJDD(3) + (-t240 + t237) * qJD(3) + mrSges(5,2) * t215 - mrSges(5,3) * t205 - qJ(5) * t178;
t283 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t314 + Ifges(4,6) * t317) * qJD(2);
t148 = -mrSges(4,1) * t235 + mrSges(4,3) * t217 + Ifges(4,4) * t293 + Ifges(4,2) * t294 + Ifges(4,6) * qJDD(3) - pkin(3) * t326 + qJ(4) * t338 + qJD(3) * t285 + t159 * t354 + t308 * t165 - t283 * t346;
t151 = mrSges(4,2) * t235 - mrSges(4,3) * t216 + Ifges(4,1) * t293 + Ifges(4,4) * t294 + Ifges(4,5) * qJDD(3) - qJ(4) * t171 - qJD(3) * t284 - t159 * t308 + t165 * t354 + t283 * t345;
t145 = mrSges(3,1) * t243 - mrSges(3,2) * t244 + Ifges(3,3) * qJDD(2) + pkin(2) * t323 + pkin(8) * t339 + t317 * t148 + t314 * t151;
t147 = mrSges(3,2) * t273 - mrSges(3,3) * t243 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t320 - pkin(8) * t164 - t148 * t314 + t151 * t317;
t333 = mrSges(2,1) * t296 - mrSges(2,2) * t297 + pkin(1) * t155 + t312 * t145 + t147 * t351 + t310 * t359;
t156 = m(2) * t297 + t158;
t154 = t312 * t163 + (t162 * t315 + t352) * t310;
t152 = m(2) * t296 + t155;
t143 = mrSges(2,2) * t307 - mrSges(2,3) * t296 + t318 * t147 - t315 * t150 + (-t154 * t310 - t155 * t312) * pkin(7);
t142 = -mrSges(2,1) * t307 + mrSges(2,3) * t297 - pkin(1) * t154 - t310 * t145 + (t147 * t315 + t359) * t312;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t311 * t143 - t309 * t142 - qJ(1) * (t152 * t311 + t156 * t309), t143, t147, t151, t165, -t280 * t239 - t328, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t309 * t143 + t311 * t142 + qJ(1) * (-t152 * t309 + t156 * t311), t142, t150, t148, t159, Ifges(6,4) * qJDD(3) - Ifges(6,2) * t262 + Ifges(6,6) * t261 - qJD(3) * t237 + t280 * t241 - t325, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t333, t333, t145, t360, t362, Ifges(6,5) * qJDD(3) - Ifges(6,6) * t262 + Ifges(6,3) * t261 + qJD(3) * t239 + t281 * t241 - t327, t331;];
m_new  = t1;
