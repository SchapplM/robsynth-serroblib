% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 09:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:05:32
% EndTime: 2019-05-05 09:06:16
% DurationCPUTime: 25.04s
% Computational Cost: add. (462530->356), mult. (956493->447), div. (0->0), fcn. (744400->14), ass. (0->152)
t313 = sin(pkin(7));
t320 = sin(qJ(3));
t323 = cos(qJ(3));
t344 = qJD(2) * qJD(3);
t294 = (-qJDD(2) * t323 + t320 * t344) * t313;
t312 = sin(pkin(12));
t315 = cos(pkin(12));
t304 = g(1) * t312 - g(2) * t315;
t305 = -g(1) * t315 - g(2) * t312;
t311 = -g(3) + qJDD(1);
t321 = sin(qJ(2));
t317 = cos(pkin(6));
t324 = cos(qJ(2));
t348 = t317 * t324;
t314 = sin(pkin(6));
t351 = t314 * t324;
t264 = t304 * t348 - t305 * t321 + t311 * t351;
t325 = qJD(2) ^ 2;
t357 = pkin(9) * t313;
t258 = qJDD(2) * pkin(2) + t325 * t357 + t264;
t349 = t317 * t321;
t352 = t314 * t321;
t265 = t304 * t349 + t324 * t305 + t311 * t352;
t259 = -pkin(2) * t325 + qJDD(2) * t357 + t265;
t285 = -t304 * t314 + t311 * t317;
t316 = cos(pkin(7));
t220 = -t320 * t259 + (t258 * t316 + t285 * t313) * t323;
t350 = t316 * t320;
t353 = t313 * t320;
t221 = t258 * t350 + t323 * t259 + t285 * t353;
t345 = qJD(2) * t313;
t292 = (-pkin(3) * t323 - pkin(10) * t320) * t345;
t310 = qJD(2) * t316 + qJD(3);
t308 = t310 ^ 2;
t309 = qJDD(2) * t316 + qJDD(3);
t342 = t323 * t345;
t216 = -pkin(3) * t308 + pkin(10) * t309 + t292 * t342 + t221;
t279 = t316 * t285;
t293 = (qJDD(2) * t320 + t323 * t344) * t313;
t218 = pkin(3) * t294 - pkin(10) * t293 + t279 + (-t258 + (pkin(3) * t320 - pkin(10) * t323) * t310 * qJD(2)) * t313;
t319 = sin(qJ(4));
t359 = cos(qJ(4));
t211 = -t319 * t216 + t218 * t359;
t212 = t359 * t216 + t319 * t218;
t343 = t320 * t345;
t283 = -t310 * t359 + t319 * t343;
t284 = t319 * t310 + t343 * t359;
t302 = -qJD(4) + t342;
t242 = Ifges(5,4) * t284 - Ifges(5,2) * t283 - Ifges(5,6) * t302;
t253 = qJD(4) * t284 + t293 * t319 - t309 * t359;
t254 = -t283 * qJD(4) + t293 * t359 + t319 * t309;
t262 = -mrSges(6,2) * t283 - mrSges(6,3) * t284;
t270 = mrSges(6,1) * t283 + mrSges(6,3) * t302;
t288 = qJDD(4) + t294;
t260 = pkin(4) * t283 - qJ(5) * t284;
t301 = t302 ^ 2;
t209 = -t288 * pkin(4) - t301 * qJ(5) + t284 * t260 + qJDD(5) - t211;
t354 = t283 * t302;
t203 = (t283 * t284 - t288) * pkin(11) + (t254 - t354) * pkin(5) + t209;
t272 = pkin(5) * t284 + pkin(11) * t302;
t282 = t283 ^ 2;
t215 = -pkin(3) * t309 - pkin(10) * t308 + t292 * t343 - t220;
t360 = -2 * qJD(5);
t327 = (-t254 - t354) * qJ(5) + t215 + (-pkin(4) * t302 + t360) * t284;
t206 = -pkin(5) * t282 - t272 * t284 + (pkin(4) + pkin(11)) * t253 + t327;
t318 = sin(qJ(6));
t322 = cos(qJ(6));
t200 = t203 * t322 - t206 * t318;
t266 = t283 * t322 + t302 * t318;
t226 = qJD(6) * t266 + t253 * t318 + t288 * t322;
t267 = t283 * t318 - t302 * t322;
t231 = -mrSges(7,1) * t266 + mrSges(7,2) * t267;
t281 = qJD(6) + t284;
t236 = -mrSges(7,2) * t281 + mrSges(7,3) * t266;
t250 = qJDD(6) + t254;
t197 = m(7) * t200 + mrSges(7,1) * t250 - mrSges(7,3) * t226 - t231 * t267 + t236 * t281;
t201 = t203 * t318 + t206 * t322;
t225 = -qJD(6) * t267 + t253 * t322 - t288 * t318;
t237 = mrSges(7,1) * t281 - mrSges(7,3) * t267;
t198 = m(7) * t201 - mrSges(7,2) * t250 + mrSges(7,3) * t225 + t231 * t266 - t237 * t281;
t186 = t197 * t322 + t198 * t318;
t334 = -pkin(4) * t301 + qJ(5) * t288 - t260 * t283 + t212;
t205 = -pkin(5) * t253 - pkin(11) * t282 + (t360 - t272) * t302 + t334;
t227 = Ifges(7,5) * t267 + Ifges(7,6) * t266 + Ifges(7,3) * t281;
t229 = Ifges(7,1) * t267 + Ifges(7,4) * t266 + Ifges(7,5) * t281;
t189 = -mrSges(7,1) * t205 + mrSges(7,3) * t201 + Ifges(7,4) * t226 + Ifges(7,2) * t225 + Ifges(7,6) * t250 - t227 * t267 + t229 * t281;
t228 = Ifges(7,4) * t267 + Ifges(7,2) * t266 + Ifges(7,6) * t281;
t190 = mrSges(7,2) * t205 - mrSges(7,3) * t200 + Ifges(7,1) * t226 + Ifges(7,4) * t225 + Ifges(7,5) * t250 + t227 * t266 - t228 * t281;
t207 = 0.2e1 * qJD(5) * t302 - t334;
t239 = -Ifges(6,5) * t302 - Ifges(6,6) * t284 + Ifges(6,3) * t283;
t331 = -mrSges(6,2) * t209 + mrSges(6,3) * t207 - Ifges(6,1) * t288 + Ifges(6,4) * t254 - Ifges(6,5) * t253 + pkin(11) * t186 + t318 * t189 - t322 * t190 + t284 * t239;
t202 = -m(7) * t205 + mrSges(7,1) * t225 - t226 * mrSges(7,2) + t236 * t266 - t267 * t237;
t271 = mrSges(6,1) * t284 - mrSges(6,2) * t302;
t332 = -m(6) * t207 + t288 * mrSges(6,3) - t302 * t271 - t202;
t335 = -m(6) * t209 - t254 * mrSges(6,1) - t284 * t262 - t186;
t241 = -Ifges(6,4) * t302 - Ifges(6,2) * t284 + Ifges(6,6) * t283;
t346 = -Ifges(5,1) * t284 + Ifges(5,4) * t283 + Ifges(5,5) * t302 + t241;
t361 = -t283 * t346 + mrSges(5,1) * t211 - mrSges(5,2) * t212 + Ifges(5,5) * t254 - Ifges(5,6) * t253 + Ifges(5,3) * t288 + pkin(4) * (-mrSges(6,2) * t288 + t270 * t302 + t335) + qJ(5) * (-mrSges(6,1) * t253 - t262 * t283 + t332) + t284 * t242 - t331;
t289 = mrSges(4,1) * t310 - mrSges(4,3) * t343;
t291 = (-mrSges(4,1) * t323 + mrSges(4,2) * t320) * t345;
t261 = mrSges(5,1) * t283 + mrSges(5,2) * t284;
t268 = mrSges(5,2) * t302 - mrSges(5,3) * t283;
t183 = m(5) * t211 - mrSges(5,3) * t254 - t261 * t284 + (-t268 + t270) * t302 + (mrSges(5,1) - mrSges(6,2)) * t288 + t335;
t269 = -mrSges(5,1) * t302 - mrSges(5,3) * t284;
t193 = m(5) * t212 - mrSges(5,2) * t288 + t269 * t302 + (-t261 - t262) * t283 + (-mrSges(5,3) - mrSges(6,1)) * t253 + t332;
t341 = -t319 * t183 + t359 * t193;
t175 = m(4) * t221 - mrSges(4,2) * t309 - mrSges(4,3) * t294 - t289 * t310 + t291 * t342 + t341;
t178 = t359 * t183 + t319 * t193;
t233 = -t258 * t313 + t279;
t290 = -mrSges(4,2) * t310 + mrSges(4,3) * t342;
t177 = m(4) * t233 + mrSges(4,1) * t294 + mrSges(4,2) * t293 + (t289 * t320 - t290 * t323) * t345 + t178;
t187 = -t318 * t197 + t322 * t198;
t210 = pkin(4) * t253 + t327;
t339 = -m(6) * t210 + t253 * mrSges(6,2) + t283 * t270 - t187;
t328 = -m(5) * t215 - t253 * mrSges(5,1) - t283 * t268 + (-t269 + t271) * t284 + (-mrSges(5,2) + mrSges(6,3)) * t254 + t339;
t181 = m(4) * t220 + mrSges(4,1) * t309 - mrSges(4,3) * t293 + t290 * t310 - t291 * t343 + t328;
t355 = t181 * t323;
t165 = t175 * t350 - t177 * t313 + t316 * t355;
t162 = m(3) * t264 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t325 + t165;
t171 = t323 * t175 - t181 * t320;
t170 = m(3) * t265 - mrSges(3,1) * t325 - qJDD(2) * mrSges(3,2) + t171;
t159 = -t162 * t321 + t324 * t170;
t358 = pkin(8) * t159;
t356 = Ifges(5,4) + Ifges(6,6);
t243 = -Ifges(6,1) * t302 - Ifges(6,4) * t284 + Ifges(6,5) * t283;
t347 = -Ifges(5,5) * t284 + Ifges(5,6) * t283 + Ifges(5,3) * t302 - t243;
t164 = t175 * t353 + t316 * t177 + t313 * t355;
t163 = m(3) * t285 + t164;
t154 = t162 * t348 - t163 * t314 + t170 * t349;
t185 = -t254 * mrSges(6,3) - t284 * t271 - t339;
t330 = -mrSges(6,1) * t207 + mrSges(6,2) * t210 - pkin(5) * t202 - pkin(11) * t187 - t322 * t189 - t318 * t190;
t166 = -mrSges(5,1) * t215 + mrSges(5,3) * t212 - pkin(4) * t185 + t346 * t302 + (Ifges(5,6) - Ifges(6,5)) * t288 + t347 * t284 + t356 * t254 + (-Ifges(5,2) - Ifges(6,3)) * t253 + t330;
t333 = mrSges(7,1) * t200 - mrSges(7,2) * t201 + Ifges(7,5) * t226 + Ifges(7,6) * t225 + Ifges(7,3) * t250 + t267 * t228 - t266 * t229;
t329 = mrSges(6,1) * t209 - mrSges(6,3) * t210 + pkin(5) * t186 + t333;
t167 = (t242 - t239) * t302 + (Ifges(5,5) - Ifges(6,4)) * t288 + t347 * t283 + (Ifges(5,1) + Ifges(6,2)) * t254 - t356 * t253 + t329 - mrSges(5,3) * t211 + mrSges(5,2) * t215 - qJ(5) * t185;
t275 = Ifges(4,3) * t310 + (Ifges(4,5) * t320 + Ifges(4,6) * t323) * t345;
t276 = Ifges(4,6) * t310 + (Ifges(4,4) * t320 + Ifges(4,2) * t323) * t345;
t156 = mrSges(4,2) * t233 - mrSges(4,3) * t220 + Ifges(4,1) * t293 - Ifges(4,4) * t294 + Ifges(4,5) * t309 - pkin(10) * t178 - t319 * t166 + t167 * t359 + t275 * t342 - t310 * t276;
t277 = Ifges(4,5) * t310 + (Ifges(4,1) * t320 + Ifges(4,4) * t323) * t345;
t160 = -mrSges(4,1) * t233 + mrSges(4,3) * t221 + Ifges(4,4) * t293 - Ifges(4,2) * t294 + Ifges(4,6) * t309 - pkin(3) * t178 - t275 * t343 + t310 * t277 - t361;
t337 = pkin(9) * t171 + t156 * t320 + t160 * t323;
t155 = Ifges(4,5) * t293 - Ifges(4,6) * t294 + Ifges(4,3) * t309 + mrSges(4,1) * t220 - mrSges(4,2) * t221 + t319 * t167 + t359 * t166 + pkin(3) * t328 + pkin(10) * t341 + (t276 * t320 - t277 * t323) * t345;
t146 = mrSges(3,1) * t264 - mrSges(3,2) * t265 + Ifges(3,3) * qJDD(2) + pkin(2) * t165 + t155 * t316 + t313 * t337;
t148 = -mrSges(3,1) * t285 + mrSges(3,3) * t265 + Ifges(3,5) * t325 + Ifges(3,6) * qJDD(2) - pkin(2) * t164 - t155 * t313 + t316 * t337;
t150 = mrSges(3,2) * t285 - mrSges(3,3) * t264 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t325 + t156 * t323 - t160 * t320 + (-t164 * t313 - t165 * t316) * pkin(9);
t336 = mrSges(2,1) * t304 - mrSges(2,2) * t305 + pkin(1) * t154 + t317 * t146 + t148 * t351 + t150 * t352 + t314 * t358;
t157 = m(2) * t305 + t159;
t153 = t163 * t317 + (t162 * t324 + t170 * t321) * t314;
t151 = m(2) * t304 + t154;
t144 = mrSges(2,2) * t311 - mrSges(2,3) * t304 - t148 * t321 + t150 * t324 + (-t153 * t314 - t154 * t317) * pkin(8);
t143 = -mrSges(2,1) * t311 + mrSges(2,3) * t305 - pkin(1) * t153 - t146 * t314 + (t148 * t324 + t150 * t321 + t358) * t317;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t315 * t144 - t312 * t143 - qJ(1) * (t151 * t315 + t157 * t312), t144, t150, t156, t167, -t283 * t241 - t331, t190; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t312 * t144 + t315 * t143 + qJ(1) * (-t151 * t312 + t157 * t315), t143, t148, t160, t166, Ifges(6,4) * t288 - Ifges(6,2) * t254 + Ifges(6,6) * t253 + t302 * t239 + t283 * t243 - t329, t189; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t336, t336, t146, t155, t361, Ifges(6,5) * t288 - Ifges(6,6) * t254 + Ifges(6,3) * t253 - t302 * t241 + t284 * t243 - t330, t333;];
m_new  = t1;
