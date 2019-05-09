% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRR5
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
% Datum: 2019-05-05 12:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:56:04
% EndTime: 2019-05-05 11:57:31
% DurationCPUTime: 70.80s
% Computational Cost: add. (1376394->352), mult. (2830566->466), div. (0->0), fcn. (2293492->16), ass. (0->158)
t313 = sin(pkin(7));
t321 = sin(qJ(3));
t326 = cos(qJ(3));
t342 = qJD(2) * qJD(3);
t296 = (-qJDD(2) * t326 + t321 * t342) * t313;
t312 = sin(pkin(13));
t315 = cos(pkin(13));
t304 = g(1) * t312 - g(2) * t315;
t305 = -g(1) * t315 - g(2) * t312;
t311 = -g(3) + qJDD(1);
t322 = sin(qJ(2));
t317 = cos(pkin(6));
t327 = cos(qJ(2));
t344 = t317 * t327;
t314 = sin(pkin(6));
t347 = t314 * t327;
t270 = t304 * t344 - t305 * t322 + t311 * t347;
t328 = qJD(2) ^ 2;
t351 = pkin(9) * t313;
t266 = qJDD(2) * pkin(2) + t328 * t351 + t270;
t345 = t317 * t322;
t348 = t314 * t322;
t271 = t304 * t345 + t305 * t327 + t311 * t348;
t267 = -pkin(2) * t328 + qJDD(2) * t351 + t271;
t288 = -t304 * t314 + t311 * t317;
t316 = cos(pkin(7));
t234 = -t321 * t267 + (t266 * t316 + t288 * t313) * t326;
t346 = t316 * t321;
t349 = t313 * t321;
t235 = t266 * t346 + t267 * t326 + t288 * t349;
t310 = qJD(2) * t316 + qJD(3);
t343 = qJD(2) * t313;
t341 = t321 * t343;
t291 = mrSges(4,1) * t310 - mrSges(4,3) * t341;
t293 = (-mrSges(4,1) * t326 + mrSges(4,2) * t321) * t343;
t309 = qJDD(2) * t316 + qJDD(3);
t294 = (-pkin(3) * t326 - pkin(10) * t321) * t343;
t308 = t310 ^ 2;
t340 = t326 * t343;
t230 = -pkin(3) * t308 + pkin(10) * t309 + t294 * t340 + t235;
t282 = t316 * t288;
t295 = (qJDD(2) * t321 + t326 * t342) * t313;
t233 = t296 * pkin(3) - t295 * pkin(10) + t282 + (-t266 + (pkin(3) * t321 - pkin(10) * t326) * t310 * qJD(2)) * t313;
t320 = sin(qJ(4));
t325 = cos(qJ(4));
t221 = t230 * t325 + t233 * t320;
t286 = t310 * t325 - t320 * t341;
t287 = t310 * t320 + t325 * t341;
t269 = -pkin(4) * t286 - pkin(11) * t287;
t290 = qJDD(4) + t296;
t303 = qJD(4) - t340;
t301 = t303 ^ 2;
t211 = -pkin(4) * t301 + pkin(11) * t290 + t269 * t286 + t221;
t229 = -t309 * pkin(3) - t308 * pkin(10) + t294 * t341 - t234;
t261 = -qJD(4) * t287 - t295 * t320 + t309 * t325;
t262 = qJD(4) * t286 + t295 * t325 + t309 * t320;
t214 = (-t286 * t303 - t262) * pkin(11) + (t287 * t303 - t261) * pkin(4) + t229;
t319 = sin(qJ(5));
t324 = cos(qJ(5));
t206 = -t319 * t211 + t214 * t324;
t273 = -t287 * t319 + t303 * t324;
t238 = qJD(5) * t273 + t262 * t324 + t290 * t319;
t259 = qJDD(5) - t261;
t274 = t287 * t324 + t303 * t319;
t285 = qJD(5) - t286;
t204 = (t273 * t285 - t238) * pkin(12) + (t273 * t274 + t259) * pkin(5) + t206;
t207 = t211 * t324 + t214 * t319;
t237 = -qJD(5) * t274 - t262 * t319 + t290 * t324;
t252 = pkin(5) * t285 - pkin(12) * t274;
t272 = t273 ^ 2;
t205 = -pkin(5) * t272 + pkin(12) * t237 - t252 * t285 + t207;
t318 = sin(qJ(6));
t323 = cos(qJ(6));
t202 = t204 * t323 - t205 * t318;
t245 = t273 * t323 - t274 * t318;
t219 = qJD(6) * t245 + t237 * t318 + t238 * t323;
t246 = t273 * t318 + t274 * t323;
t231 = -mrSges(7,1) * t245 + mrSges(7,2) * t246;
t283 = qJD(6) + t285;
t239 = -mrSges(7,2) * t283 + mrSges(7,3) * t245;
t254 = qJDD(6) + t259;
t198 = m(7) * t202 + mrSges(7,1) * t254 - mrSges(7,3) * t219 - t231 * t246 + t239 * t283;
t203 = t204 * t318 + t205 * t323;
t218 = -qJD(6) * t246 + t237 * t323 - t238 * t318;
t240 = mrSges(7,1) * t283 - mrSges(7,3) * t246;
t199 = m(7) * t203 - mrSges(7,2) * t254 + mrSges(7,3) * t218 + t231 * t245 - t240 * t283;
t190 = t198 * t323 + t199 * t318;
t247 = -mrSges(6,1) * t273 + mrSges(6,2) * t274;
t250 = -mrSges(6,2) * t285 + mrSges(6,3) * t273;
t188 = m(6) * t206 + mrSges(6,1) * t259 - mrSges(6,3) * t238 - t247 * t274 + t250 * t285 + t190;
t251 = mrSges(6,1) * t285 - mrSges(6,3) * t274;
t338 = -t198 * t318 + t199 * t323;
t189 = m(6) * t207 - mrSges(6,2) * t259 + mrSges(6,3) * t237 + t247 * t273 - t251 * t285 + t338;
t186 = -t188 * t319 + t189 * t324;
t268 = -mrSges(5,1) * t286 + mrSges(5,2) * t287;
t276 = mrSges(5,1) * t303 - mrSges(5,3) * t287;
t184 = m(5) * t221 - mrSges(5,2) * t290 + mrSges(5,3) * t261 + t268 * t286 - t276 * t303 + t186;
t220 = -t230 * t320 + t233 * t325;
t210 = -pkin(4) * t290 - pkin(11) * t301 + t269 * t287 - t220;
t208 = -pkin(5) * t237 - pkin(12) * t272 + t252 * t274 + t210;
t333 = m(7) * t208 - mrSges(7,1) * t218 + mrSges(7,2) * t219 - t239 * t245 + t240 * t246;
t200 = -m(6) * t210 + mrSges(6,1) * t237 - mrSges(6,2) * t238 + t250 * t273 - t251 * t274 - t333;
t275 = -mrSges(5,2) * t303 + mrSges(5,3) * t286;
t194 = m(5) * t220 + mrSges(5,1) * t290 - mrSges(5,3) * t262 - t268 * t287 + t275 * t303 + t200;
t339 = t184 * t325 - t194 * t320;
t173 = m(4) * t235 - mrSges(4,2) * t309 - mrSges(4,3) * t296 - t291 * t310 + t293 * t340 + t339;
t176 = t184 * t320 + t194 * t325;
t248 = -t313 * t266 + t282;
t292 = -mrSges(4,2) * t310 + mrSges(4,3) * t340;
t175 = m(4) * t248 + t296 * mrSges(4,1) + t295 * mrSges(4,2) + (t291 * t321 - t292 * t326) * t343 + t176;
t185 = t188 * t324 + t189 * t319;
t331 = -m(5) * t229 + mrSges(5,1) * t261 - mrSges(5,2) * t262 + t275 * t286 - t276 * t287 - t185;
t181 = m(4) * t234 + mrSges(4,1) * t309 - mrSges(4,3) * t295 + t292 * t310 - t293 * t341 + t331;
t350 = t181 * t326;
t163 = t173 * t346 - t175 * t313 + t316 * t350;
t160 = m(3) * t270 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t328 + t163;
t168 = t173 * t326 - t181 * t321;
t167 = m(3) * t271 - mrSges(3,1) * t328 - qJDD(2) * mrSges(3,2) + t168;
t157 = -t160 * t322 + t167 * t327;
t352 = pkin(8) * t157;
t162 = t173 * t349 + t175 * t316 + t313 * t350;
t161 = m(3) * t288 + t162;
t152 = t160 * t344 - t161 * t314 + t167 * t345;
t223 = Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t283;
t225 = Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t283;
t191 = -mrSges(7,1) * t208 + mrSges(7,3) * t203 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t254 - t223 * t246 + t225 * t283;
t224 = Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t283;
t192 = mrSges(7,2) * t208 - mrSges(7,3) * t202 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t254 + t223 * t245 - t224 * t283;
t241 = Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * t285;
t243 = Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * t285;
t177 = -mrSges(6,1) * t210 + mrSges(6,3) * t207 + Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t259 - pkin(5) * t333 + pkin(12) * t338 + t323 * t191 + t318 * t192 - t274 * t241 + t285 * t243;
t242 = Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * t285;
t178 = mrSges(6,2) * t210 - mrSges(6,3) * t206 + Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t259 - pkin(12) * t190 - t191 * t318 + t192 * t323 + t241 * t273 - t242 * t285;
t255 = Ifges(5,5) * t287 + Ifges(5,6) * t286 + Ifges(5,3) * t303;
t256 = Ifges(5,4) * t287 + Ifges(5,2) * t286 + Ifges(5,6) * t303;
t164 = mrSges(5,2) * t229 - mrSges(5,3) * t220 + Ifges(5,1) * t262 + Ifges(5,4) * t261 + Ifges(5,5) * t290 - pkin(11) * t185 - t177 * t319 + t178 * t324 + t255 * t286 - t256 * t303;
t257 = Ifges(5,1) * t287 + Ifges(5,4) * t286 + Ifges(5,5) * t303;
t332 = -mrSges(7,1) * t202 + mrSges(7,2) * t203 - Ifges(7,5) * t219 - Ifges(7,6) * t218 - Ifges(7,3) * t254 - t224 * t246 + t225 * t245;
t329 = mrSges(6,1) * t206 - mrSges(6,2) * t207 + Ifges(6,5) * t238 + Ifges(6,6) * t237 + Ifges(6,3) * t259 + pkin(5) * t190 + t242 * t274 - t243 * t273 - t332;
t169 = -mrSges(5,1) * t229 + mrSges(5,3) * t221 + Ifges(5,4) * t262 + Ifges(5,2) * t261 + Ifges(5,6) * t290 - pkin(4) * t185 - t255 * t287 + t257 * t303 - t329;
t278 = Ifges(4,3) * t310 + (Ifges(4,5) * t321 + Ifges(4,6) * t326) * t343;
t279 = Ifges(4,6) * t310 + (Ifges(4,4) * t321 + Ifges(4,2) * t326) * t343;
t154 = mrSges(4,2) * t248 - mrSges(4,3) * t234 + Ifges(4,1) * t295 - Ifges(4,4) * t296 + Ifges(4,5) * t309 - pkin(10) * t176 + t164 * t325 - t169 * t320 + t278 * t340 - t279 * t310;
t280 = Ifges(4,5) * t310 + (Ifges(4,1) * t321 + Ifges(4,4) * t326) * t343;
t330 = mrSges(5,1) * t220 - mrSges(5,2) * t221 + Ifges(5,5) * t262 + Ifges(5,6) * t261 + Ifges(5,3) * t290 + pkin(4) * t200 + pkin(11) * t186 + t177 * t324 + t178 * t319 + t256 * t287 - t257 * t286;
t158 = -mrSges(4,1) * t248 + mrSges(4,3) * t235 + Ifges(4,4) * t295 - Ifges(4,2) * t296 + Ifges(4,6) * t309 - pkin(3) * t176 - t278 * t341 + t280 * t310 - t330;
t335 = pkin(9) * t168 + t154 * t321 + t158 * t326;
t153 = Ifges(4,5) * t295 - Ifges(4,6) * t296 + Ifges(4,3) * t309 + mrSges(4,1) * t234 - mrSges(4,2) * t235 + t320 * t164 + t325 * t169 + pkin(3) * t331 + pkin(10) * t339 + (t279 * t321 - t280 * t326) * t343;
t144 = mrSges(3,1) * t270 - mrSges(3,2) * t271 + Ifges(3,3) * qJDD(2) + pkin(2) * t163 + t316 * t153 + t313 * t335;
t146 = -mrSges(3,1) * t288 + mrSges(3,3) * t271 + t328 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t162 - t313 * t153 + t316 * t335;
t148 = mrSges(3,2) * t288 - mrSges(3,3) * t270 + Ifges(3,5) * qJDD(2) - t328 * Ifges(3,6) + t326 * t154 - t321 * t158 + (-t162 * t313 - t163 * t316) * pkin(9);
t334 = mrSges(2,1) * t304 - mrSges(2,2) * t305 + pkin(1) * t152 + t144 * t317 + t146 * t347 + t148 * t348 + t314 * t352;
t155 = m(2) * t305 + t157;
t151 = t317 * t161 + (t160 * t327 + t167 * t322) * t314;
t149 = m(2) * t304 + t152;
t142 = mrSges(2,2) * t311 - mrSges(2,3) * t304 - t322 * t146 + t327 * t148 + (-t151 * t314 - t152 * t317) * pkin(8);
t141 = -mrSges(2,1) * t311 + mrSges(2,3) * t305 - pkin(1) * t151 - t314 * t144 + (t146 * t327 + t148 * t322 + t352) * t317;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t315 * t142 - t312 * t141 - qJ(1) * (t149 * t315 + t155 * t312), t142, t148, t154, t164, t178, t192; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t312 * t142 + t315 * t141 + qJ(1) * (-t149 * t312 + t155 * t315), t141, t146, t158, t169, t177, t191; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t334, t334, t144, t153, t330, t329, -t332;];
m_new  = t1;
