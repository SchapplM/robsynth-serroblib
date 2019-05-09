% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRR3
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
% Datum: 2019-05-05 04:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:45:16
% EndTime: 2019-05-05 04:46:22
% DurationCPUTime: 59.14s
% Computational Cost: add. (1052491->349), mult. (2556952->467), div. (0->0), fcn. (2038362->16), ass. (0->158)
t364 = -2 * qJD(4);
t319 = sin(pkin(13));
t323 = cos(pkin(13));
t329 = sin(qJ(3));
t333 = cos(qJ(3));
t321 = sin(pkin(7));
t352 = qJD(2) * t321;
t295 = (t319 * t329 - t323 * t333) * t352;
t320 = sin(pkin(12));
t324 = cos(pkin(12));
t308 = t320 * g(1) - t324 * g(2);
t309 = -t324 * g(1) - t320 * g(2);
t318 = -g(3) + qJDD(1);
t330 = sin(qJ(2));
t326 = cos(pkin(6));
t334 = cos(qJ(2));
t353 = t326 * t334;
t322 = sin(pkin(6));
t357 = t322 * t334;
t277 = t308 * t353 - t330 * t309 + t318 * t357;
t335 = qJD(2) ^ 2;
t362 = pkin(9) * t321;
t267 = qJDD(2) * pkin(2) + t335 * t362 + t277;
t354 = t326 * t330;
t358 = t322 * t330;
t278 = t308 * t354 + t334 * t309 + t318 * t358;
t268 = -t335 * pkin(2) + qJDD(2) * t362 + t278;
t297 = -t322 * t308 + t326 * t318;
t325 = cos(pkin(7));
t355 = t325 * t333;
t359 = t321 * t333;
t238 = t267 * t355 - t329 * t268 + t297 * t359;
t350 = qJD(2) * qJD(3);
t302 = (qJDD(2) * t329 + t333 * t350) * t321;
t314 = t325 * qJDD(2) + qJDD(3);
t315 = t325 * qJD(2) + qJD(3);
t347 = t333 * t352;
t361 = t321 ^ 2 * t335;
t230 = (t315 * t347 - t302) * qJ(4) + (t329 * t333 * t361 + t314) * pkin(3) + t238;
t356 = t325 * t329;
t360 = t321 * t329;
t239 = t267 * t356 + t333 * t268 + t297 * t360;
t348 = t329 * t352;
t298 = t315 * pkin(3) - qJ(4) * t348;
t303 = (qJDD(2) * t333 - t329 * t350) * t321;
t349 = t333 ^ 2 * t361;
t231 = -pkin(3) * t349 + t303 * qJ(4) - t315 * t298 + t239;
t296 = (t319 * t333 + t323 * t329) * t352;
t220 = t323 * t230 - t319 * t231 + t296 * t364;
t221 = t319 * t230 + t323 * t231 + t295 * t364;
t269 = t295 * mrSges(5,1) + t296 * mrSges(5,2);
t275 = -t319 * t302 + t323 * t303;
t283 = t315 * mrSges(5,1) - t296 * mrSges(5,3);
t270 = t295 * pkin(4) - t296 * pkin(10);
t313 = t315 ^ 2;
t218 = -t313 * pkin(4) + t314 * pkin(10) - t295 * t270 + t221;
t252 = -t321 * t267 + t325 * t297;
t240 = -t303 * pkin(3) - qJ(4) * t349 + t298 * t348 + qJDD(4) + t252;
t276 = t323 * t302 + t319 * t303;
t223 = (t295 * t315 - t276) * pkin(10) + (t296 * t315 - t275) * pkin(4) + t240;
t328 = sin(qJ(5));
t332 = cos(qJ(5));
t215 = t332 * t218 + t328 * t223;
t280 = -t328 * t296 + t332 * t315;
t281 = t332 * t296 + t328 * t315;
t254 = -t280 * pkin(5) - t281 * pkin(11);
t274 = qJDD(5) - t275;
t294 = qJD(5) + t295;
t293 = t294 ^ 2;
t212 = -t293 * pkin(5) + t274 * pkin(11) + t280 * t254 + t215;
t217 = -t314 * pkin(4) - t313 * pkin(10) + t296 * t270 - t220;
t249 = -t281 * qJD(5) - t328 * t276 + t332 * t314;
t250 = t280 * qJD(5) + t332 * t276 + t328 * t314;
t213 = (-t280 * t294 - t250) * pkin(11) + (t281 * t294 - t249) * pkin(5) + t217;
t327 = sin(qJ(6));
t331 = cos(qJ(6));
t208 = -t327 * t212 + t331 * t213;
t256 = -t327 * t281 + t331 * t294;
t226 = t256 * qJD(6) + t331 * t250 + t327 * t274;
t257 = t331 * t281 + t327 * t294;
t236 = -t256 * mrSges(7,1) + t257 * mrSges(7,2);
t279 = qJD(6) - t280;
t241 = -t279 * mrSges(7,2) + t256 * mrSges(7,3);
t248 = qJDD(6) - t249;
t206 = m(7) * t208 + t248 * mrSges(7,1) - t226 * mrSges(7,3) - t257 * t236 + t279 * t241;
t209 = t331 * t212 + t327 * t213;
t225 = -t257 * qJD(6) - t327 * t250 + t331 * t274;
t242 = t279 * mrSges(7,1) - t257 * mrSges(7,3);
t207 = m(7) * t209 - t248 * mrSges(7,2) + t225 * mrSges(7,3) + t256 * t236 - t279 * t242;
t200 = -t327 * t206 + t331 * t207;
t253 = -t280 * mrSges(6,1) + t281 * mrSges(6,2);
t259 = t294 * mrSges(6,1) - t281 * mrSges(6,3);
t197 = m(6) * t215 - t274 * mrSges(6,2) + t249 * mrSges(6,3) + t280 * t253 - t294 * t259 + t200;
t214 = -t328 * t218 + t332 * t223;
t211 = -t274 * pkin(5) - t293 * pkin(11) + t281 * t254 - t214;
t210 = -m(7) * t211 + t225 * mrSges(7,1) - t226 * mrSges(7,2) + t256 * t241 - t257 * t242;
t258 = -t294 * mrSges(6,2) + t280 * mrSges(6,3);
t204 = m(6) * t214 + t274 * mrSges(6,1) - t250 * mrSges(6,3) - t281 * t253 + t294 * t258 + t210;
t345 = t332 * t197 - t328 * t204;
t188 = m(5) * t221 - t314 * mrSges(5,2) + t275 * mrSges(5,3) - t295 * t269 - t315 * t283 + t345;
t282 = -t315 * mrSges(5,2) - t295 * mrSges(5,3);
t199 = t331 * t206 + t327 * t207;
t338 = -m(6) * t217 + t249 * mrSges(6,1) - t250 * mrSges(6,2) + t280 * t258 - t281 * t259 - t199;
t194 = m(5) * t220 + t314 * mrSges(5,1) - t276 * mrSges(5,3) - t296 * t269 + t315 * t282 + t338;
t181 = t319 * t188 + t323 * t194;
t300 = -t315 * mrSges(4,2) + mrSges(4,3) * t347;
t301 = (-mrSges(4,1) * t333 + mrSges(4,2) * t329) * t352;
t179 = m(4) * t238 + t314 * mrSges(4,1) - t302 * mrSges(4,3) + t315 * t300 - t301 * t348 + t181;
t299 = t315 * mrSges(4,1) - mrSges(4,3) * t348;
t346 = t323 * t188 - t319 * t194;
t180 = m(4) * t239 - t314 * mrSges(4,2) + t303 * mrSges(4,3) - t315 * t299 + t301 * t347 + t346;
t192 = t328 * t197 + t332 * t204;
t340 = m(5) * t240 - t275 * mrSges(5,1) + t276 * mrSges(5,2) + t295 * t282 + t296 * t283 + t192;
t190 = m(4) * t252 - t303 * mrSges(4,1) + t302 * mrSges(4,2) + (t299 * t329 - t300 * t333) * t352 + t340;
t167 = t179 * t355 + t180 * t356 - t321 * t190;
t164 = m(3) * t277 + qJDD(2) * mrSges(3,1) - t335 * mrSges(3,2) + t167;
t172 = -t329 * t179 + t333 * t180;
t171 = m(3) * t278 - t335 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t172;
t161 = -t330 * t164 + t334 * t171;
t363 = pkin(8) * t161;
t166 = t179 * t359 + t180 * t360 + t325 * t190;
t165 = m(3) * t297 + t166;
t156 = t164 * t353 - t322 * t165 + t171 * t354;
t232 = Ifges(7,5) * t257 + Ifges(7,6) * t256 + Ifges(7,3) * t279;
t234 = Ifges(7,1) * t257 + Ifges(7,4) * t256 + Ifges(7,5) * t279;
t201 = -mrSges(7,1) * t211 + mrSges(7,3) * t209 + Ifges(7,4) * t226 + Ifges(7,2) * t225 + Ifges(7,6) * t248 - t257 * t232 + t279 * t234;
t233 = Ifges(7,4) * t257 + Ifges(7,2) * t256 + Ifges(7,6) * t279;
t202 = mrSges(7,2) * t211 - mrSges(7,3) * t208 + Ifges(7,1) * t226 + Ifges(7,4) * t225 + Ifges(7,5) * t248 + t256 * t232 - t279 * t233;
t243 = Ifges(6,5) * t281 + Ifges(6,6) * t280 + Ifges(6,3) * t294;
t244 = Ifges(6,4) * t281 + Ifges(6,2) * t280 + Ifges(6,6) * t294;
t183 = mrSges(6,2) * t217 - mrSges(6,3) * t214 + Ifges(6,1) * t250 + Ifges(6,4) * t249 + Ifges(6,5) * t274 - pkin(11) * t199 - t327 * t201 + t331 * t202 + t280 * t243 - t294 * t244;
t245 = Ifges(6,1) * t281 + Ifges(6,4) * t280 + Ifges(6,5) * t294;
t337 = mrSges(7,1) * t208 - mrSges(7,2) * t209 + Ifges(7,5) * t226 + Ifges(7,6) * t225 + Ifges(7,3) * t248 + t257 * t233 - t256 * t234;
t185 = -mrSges(6,1) * t217 + mrSges(6,3) * t215 + Ifges(6,4) * t250 + Ifges(6,2) * t249 + Ifges(6,6) * t274 - pkin(5) * t199 - t281 * t243 + t294 * t245 - t337;
t260 = Ifges(5,5) * t296 - Ifges(5,6) * t295 + Ifges(5,3) * t315;
t261 = Ifges(5,4) * t296 - Ifges(5,2) * t295 + Ifges(5,6) * t315;
t168 = mrSges(5,2) * t240 - mrSges(5,3) * t220 + Ifges(5,1) * t276 + Ifges(5,4) * t275 + Ifges(5,5) * t314 - pkin(10) * t192 + t332 * t183 - t328 * t185 - t295 * t260 - t315 * t261;
t262 = Ifges(5,1) * t296 - Ifges(5,4) * t295 + Ifges(5,5) * t315;
t336 = mrSges(6,1) * t214 - mrSges(6,2) * t215 + Ifges(6,5) * t250 + Ifges(6,6) * t249 + Ifges(6,3) * t274 + pkin(5) * t210 + pkin(11) * t200 + t331 * t201 + t327 * t202 + t281 * t244 - t280 * t245;
t173 = -mrSges(5,1) * t240 + mrSges(5,3) * t221 + Ifges(5,4) * t276 + Ifges(5,2) * t275 + Ifges(5,6) * t314 - pkin(4) * t192 - t296 * t260 + t315 * t262 - t336;
t286 = Ifges(4,3) * t315 + (Ifges(4,5) * t329 + Ifges(4,6) * t333) * t352;
t288 = Ifges(4,5) * t315 + (Ifges(4,1) * t329 + Ifges(4,4) * t333) * t352;
t157 = -mrSges(4,1) * t252 + mrSges(4,3) * t239 + Ifges(4,4) * t302 + Ifges(4,2) * t303 + Ifges(4,6) * t314 - pkin(3) * t340 + qJ(4) * t346 + t319 * t168 + t323 * t173 - t286 * t348 + t315 * t288;
t287 = Ifges(4,6) * t315 + (Ifges(4,4) * t329 + Ifges(4,2) * t333) * t352;
t158 = mrSges(4,2) * t252 - mrSges(4,3) * t238 + Ifges(4,1) * t302 + Ifges(4,4) * t303 + Ifges(4,5) * t314 - qJ(4) * t181 + t323 * t168 - t319 * t173 + t286 * t347 - t315 * t287;
t342 = pkin(9) * t172 + t157 * t333 + t158 * t329;
t339 = mrSges(5,1) * t220 - mrSges(5,2) * t221 + Ifges(5,5) * t276 + Ifges(5,6) * t275 + Ifges(5,3) * t314 + pkin(4) * t338 + pkin(10) * t345 + t328 * t183 + t332 * t185 + t296 * t261 + t295 * t262;
t162 = Ifges(4,3) * t314 + Ifges(4,5) * t302 + Ifges(4,6) * t303 + mrSges(4,1) * t238 - mrSges(4,2) * t239 + pkin(3) * t181 + t339 + (t287 * t329 - t288 * t333) * t352;
t148 = mrSges(3,1) * t277 - mrSges(3,2) * t278 + Ifges(3,3) * qJDD(2) + pkin(2) * t167 + t325 * t162 + t342 * t321;
t150 = -mrSges(3,1) * t297 + mrSges(3,3) * t278 + t335 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t166 - t321 * t162 + t342 * t325;
t152 = mrSges(3,2) * t297 - mrSges(3,3) * t277 + Ifges(3,5) * qJDD(2) - t335 * Ifges(3,6) - t329 * t157 + t333 * t158 + (-t166 * t321 - t167 * t325) * pkin(9);
t341 = mrSges(2,1) * t308 - mrSges(2,2) * t309 + pkin(1) * t156 + t326 * t148 + t150 * t357 + t152 * t358 + t322 * t363;
t159 = m(2) * t309 + t161;
t155 = t326 * t165 + (t164 * t334 + t171 * t330) * t322;
t153 = m(2) * t308 + t156;
t146 = mrSges(2,2) * t318 - mrSges(2,3) * t308 - t330 * t150 + t334 * t152 + (-t155 * t322 - t156 * t326) * pkin(8);
t145 = -mrSges(2,1) * t318 + mrSges(2,3) * t309 - pkin(1) * t155 - t322 * t148 + (t150 * t334 + t152 * t330 + t363) * t326;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t324 * t146 - t320 * t145 - qJ(1) * (t324 * t153 + t320 * t159), t146, t152, t158, t168, t183, t202; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t320 * t146 + t324 * t145 + qJ(1) * (-t320 * t153 + t324 * t159), t145, t150, t157, t173, t185, t201; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t341, t341, t148, t162, t339, t336, t337;];
m_new  = t1;
