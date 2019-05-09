% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 08:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:00:18
% EndTime: 2019-05-05 08:01:28
% DurationCPUTime: 65.26s
% Computational Cost: add. (1256969->352), mult. (2660385->468), div. (0->0), fcn. (2152643->16), ass. (0->158)
t314 = sin(pkin(7));
t322 = sin(qJ(3));
t326 = cos(qJ(3));
t345 = qJD(2) * qJD(3);
t297 = (-qJDD(2) * t326 + t322 * t345) * t314;
t313 = sin(pkin(12));
t317 = cos(pkin(12));
t304 = t313 * g(1) - t317 * g(2);
t305 = -t317 * g(1) - t313 * g(2);
t311 = -g(3) + qJDD(1);
t323 = sin(qJ(2));
t319 = cos(pkin(6));
t327 = cos(qJ(2));
t347 = t319 * t327;
t315 = sin(pkin(6));
t350 = t315 * t327;
t275 = t304 * t347 - t323 * t305 + t311 * t350;
t328 = qJD(2) ^ 2;
t354 = pkin(9) * t314;
t272 = qJDD(2) * pkin(2) + t328 * t354 + t275;
t348 = t319 * t323;
t351 = t315 * t323;
t276 = t304 * t348 + t327 * t305 + t311 * t351;
t273 = -t328 * pkin(2) + qJDD(2) * t354 + t276;
t289 = -t315 * t304 + t319 * t311;
t318 = cos(pkin(7));
t236 = -t322 * t273 + (t272 * t318 + t289 * t314) * t326;
t356 = 2 * qJD(5);
t349 = t318 * t322;
t352 = t314 * t322;
t237 = t272 * t349 + t326 * t273 + t289 * t352;
t310 = t318 * qJD(2) + qJD(3);
t346 = qJD(2) * t314;
t344 = t322 * t346;
t292 = t310 * mrSges(4,1) - mrSges(4,3) * t344;
t294 = (-mrSges(4,1) * t326 + mrSges(4,2) * t322) * t346;
t309 = t318 * qJDD(2) + qJDD(3);
t295 = (-pkin(3) * t326 - pkin(10) * t322) * t346;
t308 = t310 ^ 2;
t343 = t326 * t346;
t230 = -t308 * pkin(3) + t309 * pkin(10) + t295 * t343 + t237;
t285 = t318 * t289;
t296 = (qJDD(2) * t322 + t326 * t345) * t314;
t234 = t297 * pkin(3) - t296 * pkin(10) + t285 + (-t272 + (pkin(3) * t322 - pkin(10) * t326) * t310 * qJD(2)) * t314;
t321 = sin(qJ(4));
t325 = cos(qJ(4));
t218 = -t321 * t230 + t325 * t234;
t287 = t325 * t310 - t321 * t344;
t264 = t287 * qJD(4) + t325 * t296 + t321 * t309;
t288 = t321 * t310 + t325 * t344;
t291 = qJDD(4) + t297;
t303 = qJD(4) - t343;
t215 = (t287 * t303 - t264) * qJ(5) + (t287 * t288 + t291) * pkin(4) + t218;
t219 = t325 * t230 + t321 * t234;
t263 = -t288 * qJD(4) - t321 * t296 + t325 * t309;
t278 = t303 * pkin(4) - t288 * qJ(5);
t286 = t287 ^ 2;
t217 = -t286 * pkin(4) + t263 * qJ(5) - t303 * t278 + t219;
t312 = sin(pkin(13));
t316 = cos(pkin(13));
t270 = t316 * t287 - t312 * t288;
t212 = t312 * t215 + t316 * t217 + t270 * t356;
t244 = t316 * t263 - t312 * t264;
t271 = t312 * t287 + t316 * t288;
t250 = -t270 * mrSges(6,1) + t271 * mrSges(6,2);
t256 = t303 * mrSges(6,1) - t271 * mrSges(6,3);
t251 = -t270 * pkin(5) - t271 * pkin(11);
t302 = t303 ^ 2;
t209 = -t302 * pkin(5) + t291 * pkin(11) + t270 * t251 + t212;
t229 = -t309 * pkin(3) - t308 * pkin(10) + t295 * t344 - t236;
t220 = -t263 * pkin(4) - t286 * qJ(5) + t288 * t278 + qJDD(5) + t229;
t245 = t312 * t263 + t316 * t264;
t213 = (-t270 * t303 - t245) * pkin(11) + (t271 * t303 - t244) * pkin(5) + t220;
t320 = sin(qJ(6));
t324 = cos(qJ(6));
t206 = -t320 * t209 + t324 * t213;
t253 = -t320 * t271 + t324 * t303;
t223 = t253 * qJD(6) + t324 * t245 + t320 * t291;
t254 = t324 * t271 + t320 * t303;
t235 = -t253 * mrSges(7,1) + t254 * mrSges(7,2);
t267 = qJD(6) - t270;
t238 = -t267 * mrSges(7,2) + t253 * mrSges(7,3);
t243 = qJDD(6) - t244;
t202 = m(7) * t206 + t243 * mrSges(7,1) - t223 * mrSges(7,3) - t254 * t235 + t267 * t238;
t207 = t324 * t209 + t320 * t213;
t222 = -t254 * qJD(6) - t320 * t245 + t324 * t291;
t239 = t267 * mrSges(7,1) - t254 * mrSges(7,3);
t203 = m(7) * t207 - t243 * mrSges(7,2) + t222 * mrSges(7,3) + t253 * t235 - t267 * t239;
t340 = -t320 * t202 + t324 * t203;
t189 = m(6) * t212 - t291 * mrSges(6,2) + t244 * mrSges(6,3) + t270 * t250 - t303 * t256 + t340;
t339 = -t316 * t215 + t312 * t217;
t211 = -0.2e1 * qJD(5) * t271 - t339;
t255 = -t303 * mrSges(6,2) + t270 * mrSges(6,3);
t208 = -t291 * pkin(5) - t302 * pkin(11) + (t356 + t251) * t271 + t339;
t334 = -m(7) * t208 + t222 * mrSges(7,1) - t223 * mrSges(7,2) + t253 * t238 - t254 * t239;
t198 = m(6) * t211 + t291 * mrSges(6,1) - t245 * mrSges(6,3) - t271 * t250 + t303 * t255 + t334;
t183 = t312 * t189 + t316 * t198;
t274 = -t287 * mrSges(5,1) + t288 * mrSges(5,2);
t277 = -t303 * mrSges(5,2) + t287 * mrSges(5,3);
t181 = m(5) * t218 + t291 * mrSges(5,1) - t264 * mrSges(5,3) - t288 * t274 + t303 * t277 + t183;
t279 = t303 * mrSges(5,1) - t288 * mrSges(5,3);
t341 = t316 * t189 - t312 * t198;
t182 = m(5) * t219 - t291 * mrSges(5,2) + t263 * mrSges(5,3) + t287 * t274 - t303 * t279 + t341;
t342 = -t321 * t181 + t325 * t182;
t172 = m(4) * t237 - t309 * mrSges(4,2) - t297 * mrSges(4,3) - t310 * t292 + t294 * t343 + t342;
t175 = t325 * t181 + t321 * t182;
t252 = -t314 * t272 + t285;
t293 = -t310 * mrSges(4,2) + mrSges(4,3) * t343;
t174 = m(4) * t252 + t297 * mrSges(4,1) + t296 * mrSges(4,2) + (t292 * t322 - t293 * t326) * t346 + t175;
t191 = t324 * t202 + t320 * t203;
t333 = m(6) * t220 - t244 * mrSges(6,1) + t245 * mrSges(6,2) - t270 * t255 + t271 * t256 + t191;
t330 = -m(5) * t229 + t263 * mrSges(5,1) - t264 * mrSges(5,2) + t287 * t277 - t288 * t279 - t333;
t186 = m(4) * t236 + t309 * mrSges(4,1) - t296 * mrSges(4,3) + t310 * t293 - t294 * t344 + t330;
t353 = t186 * t326;
t162 = t172 * t349 - t314 * t174 + t318 * t353;
t159 = m(3) * t275 + qJDD(2) * mrSges(3,1) - t328 * mrSges(3,2) + t162;
t168 = t326 * t172 - t322 * t186;
t167 = m(3) * t276 - t328 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t168;
t156 = -t323 * t159 + t327 * t167;
t355 = pkin(8) * t156;
t161 = t172 * t352 + t318 * t174 + t314 * t353;
t160 = m(3) * t289 + t161;
t151 = t159 * t347 - t315 * t160 + t167 * t348;
t224 = Ifges(7,5) * t254 + Ifges(7,6) * t253 + Ifges(7,3) * t267;
t226 = Ifges(7,1) * t254 + Ifges(7,4) * t253 + Ifges(7,5) * t267;
t195 = -mrSges(7,1) * t208 + mrSges(7,3) * t207 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t243 - t254 * t224 + t267 * t226;
t225 = Ifges(7,4) * t254 + Ifges(7,2) * t253 + Ifges(7,6) * t267;
t196 = mrSges(7,2) * t208 - mrSges(7,3) * t206 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t243 + t253 * t224 - t267 * t225;
t246 = Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t303;
t247 = Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t303;
t176 = mrSges(6,2) * t220 - mrSges(6,3) * t211 + Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t291 - pkin(11) * t191 - t320 * t195 + t324 * t196 + t270 * t246 - t303 * t247;
t248 = Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t303;
t331 = mrSges(7,1) * t206 - mrSges(7,2) * t207 + Ifges(7,5) * t223 + Ifges(7,6) * t222 + Ifges(7,3) * t243 + t254 * t225 - t253 * t226;
t177 = -mrSges(6,1) * t220 + mrSges(6,3) * t212 + Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t291 - pkin(5) * t191 - t271 * t246 + t303 * t248 - t331;
t257 = Ifges(5,5) * t288 + Ifges(5,6) * t287 + Ifges(5,3) * t303;
t259 = Ifges(5,1) * t288 + Ifges(5,4) * t287 + Ifges(5,5) * t303;
t163 = -mrSges(5,1) * t229 + mrSges(5,3) * t219 + Ifges(5,4) * t264 + Ifges(5,2) * t263 + Ifges(5,6) * t291 - pkin(4) * t333 + qJ(5) * t341 + t312 * t176 + t316 * t177 - t288 * t257 + t303 * t259;
t258 = Ifges(5,4) * t288 + Ifges(5,2) * t287 + Ifges(5,6) * t303;
t164 = mrSges(5,2) * t229 - mrSges(5,3) * t218 + Ifges(5,1) * t264 + Ifges(5,4) * t263 + Ifges(5,5) * t291 - qJ(5) * t183 + t316 * t176 - t312 * t177 + t287 * t257 - t303 * t258;
t281 = Ifges(4,3) * t310 + (Ifges(4,5) * t322 + Ifges(4,6) * t326) * t346;
t282 = Ifges(4,6) * t310 + (Ifges(4,4) * t322 + Ifges(4,2) * t326) * t346;
t153 = mrSges(4,2) * t252 - mrSges(4,3) * t236 + Ifges(4,1) * t296 - Ifges(4,4) * t297 + Ifges(4,5) * t309 - pkin(10) * t175 - t321 * t163 + t325 * t164 + t281 * t343 - t310 * t282;
t283 = Ifges(4,5) * t310 + (Ifges(4,1) * t322 + Ifges(4,4) * t326) * t346;
t332 = -mrSges(6,1) * t211 + mrSges(6,2) * t212 - Ifges(6,5) * t245 - Ifges(6,6) * t244 - Ifges(6,3) * t291 - pkin(5) * t334 - pkin(11) * t340 - t324 * t195 - t320 * t196 - t271 * t247 + t270 * t248;
t329 = mrSges(5,1) * t218 - mrSges(5,2) * t219 + Ifges(5,5) * t264 + Ifges(5,6) * t263 + Ifges(5,3) * t291 + pkin(4) * t183 + t288 * t258 - t287 * t259 - t332;
t157 = -mrSges(4,1) * t252 + mrSges(4,3) * t237 + Ifges(4,4) * t296 - Ifges(4,2) * t297 + Ifges(4,6) * t309 - pkin(3) * t175 - t281 * t344 + t310 * t283 - t329;
t336 = pkin(9) * t168 + t153 * t322 + t157 * t326;
t152 = Ifges(4,5) * t296 - Ifges(4,6) * t297 + Ifges(4,3) * t309 + mrSges(4,1) * t236 - mrSges(4,2) * t237 + t321 * t164 + t325 * t163 + pkin(3) * t330 + pkin(10) * t342 + (t282 * t322 - t283 * t326) * t346;
t143 = mrSges(3,1) * t275 - mrSges(3,2) * t276 + Ifges(3,3) * qJDD(2) + pkin(2) * t162 + t318 * t152 + t336 * t314;
t145 = -mrSges(3,1) * t289 + mrSges(3,3) * t276 + t328 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t161 - t314 * t152 + t336 * t318;
t147 = mrSges(3,2) * t289 - mrSges(3,3) * t275 + Ifges(3,5) * qJDD(2) - t328 * Ifges(3,6) + t326 * t153 - t322 * t157 + (-t161 * t314 - t162 * t318) * pkin(9);
t335 = mrSges(2,1) * t304 - mrSges(2,2) * t305 + pkin(1) * t151 + t319 * t143 + t145 * t350 + t147 * t351 + t315 * t355;
t154 = m(2) * t305 + t156;
t150 = t319 * t160 + (t159 * t327 + t167 * t323) * t315;
t148 = m(2) * t304 + t151;
t141 = mrSges(2,2) * t311 - mrSges(2,3) * t304 - t323 * t145 + t327 * t147 + (-t150 * t315 - t151 * t319) * pkin(8);
t140 = -mrSges(2,1) * t311 + mrSges(2,3) * t305 - pkin(1) * t150 - t315 * t143 + (t145 * t327 + t147 * t323 + t355) * t319;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t317 * t141 - t313 * t140 - qJ(1) * (t317 * t148 + t313 * t154), t141, t147, t153, t164, t176, t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t313 * t141 + t317 * t140 + qJ(1) * (-t313 * t148 + t317 * t154), t140, t145, t157, t163, t177, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t335, t335, t143, t152, t329, -t332, t331;];
m_new  = t1;
