% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPR7
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
% Datum: 2019-05-05 08:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:42:41
% EndTime: 2019-05-05 08:43:59
% DurationCPUTime: 68.43s
% Computational Cost: add. (1312724->352), mult. (2754186->468), div. (0->0), fcn. (2219098->16), ass. (0->156)
t312 = sin(pkin(7));
t320 = sin(qJ(3));
t323 = cos(qJ(3));
t339 = qJD(2) * qJD(3);
t293 = (-qJDD(2) * t323 + t320 * t339) * t312;
t311 = sin(pkin(12));
t315 = cos(pkin(12));
t302 = t311 * g(1) - t315 * g(2);
t303 = -t315 * g(1) - t311 * g(2);
t309 = -g(3) + qJDD(1);
t321 = sin(qJ(2));
t317 = cos(pkin(6));
t324 = cos(qJ(2));
t341 = t317 * t324;
t313 = sin(pkin(6));
t344 = t313 * t324;
t271 = t302 * t341 - t321 * t303 + t309 * t344;
t325 = qJD(2) ^ 2;
t348 = pkin(9) * t312;
t264 = qJDD(2) * pkin(2) + t325 * t348 + t271;
t342 = t317 * t321;
t345 = t313 * t321;
t272 = t302 * t342 + t324 * t303 + t309 * t345;
t265 = -t325 * pkin(2) + qJDD(2) * t348 + t272;
t286 = -t313 * t302 + t317 * t309;
t316 = cos(pkin(7));
t233 = -t320 * t265 + (t264 * t316 + t286 * t312) * t323;
t350 = cos(qJ(4));
t343 = t316 * t320;
t346 = t312 * t320;
t234 = t264 * t343 + t323 * t265 + t286 * t346;
t308 = t316 * qJD(2) + qJD(3);
t340 = qJD(2) * t312;
t338 = t320 * t340;
t288 = t308 * mrSges(4,1) - mrSges(4,3) * t338;
t290 = (-mrSges(4,1) * t323 + mrSges(4,2) * t320) * t340;
t307 = t316 * qJDD(2) + qJDD(3);
t291 = (-pkin(3) * t323 - pkin(10) * t320) * t340;
t306 = t308 ^ 2;
t337 = t323 * t340;
t229 = -t306 * pkin(3) + t307 * pkin(10) + t291 * t337 + t234;
t282 = t316 * t286;
t292 = (qJDD(2) * t320 + t323 * t339) * t312;
t232 = t293 * pkin(3) - t292 * pkin(10) + t282 + (-t264 + (pkin(3) * t320 - pkin(10) * t323) * t308 * qJD(2)) * t312;
t319 = sin(qJ(4));
t215 = t229 * t350 + t319 * t232;
t284 = -t308 * t350 + t319 * t338;
t285 = t319 * t308 + t338 * t350;
t266 = t284 * pkin(4) - t285 * qJ(5);
t287 = qJDD(4) + t293;
t301 = qJD(4) - t337;
t300 = t301 ^ 2;
t210 = -t300 * pkin(4) + t287 * qJ(5) - t284 * t266 + t215;
t228 = -t307 * pkin(3) - t306 * pkin(10) + t291 * t338 - t233;
t259 = t285 * qJD(4) + t319 * t292 - t307 * t350;
t260 = -t284 * qJD(4) + t292 * t350 + t319 * t307;
t213 = (t284 * t301 - t260) * qJ(5) + (t285 * t301 + t259) * pkin(4) + t228;
t310 = sin(pkin(13));
t314 = cos(pkin(13));
t274 = t314 * t285 + t310 * t301;
t205 = -0.2e1 * qJD(5) * t274 - t310 * t210 + t314 * t213;
t245 = t314 * t260 + t310 * t287;
t273 = -t310 * t285 + t314 * t301;
t203 = (t273 * t284 - t245) * pkin(11) + (t273 * t274 + t259) * pkin(5) + t205;
t206 = 0.2e1 * qJD(5) * t273 + t314 * t210 + t310 * t213;
t244 = -t310 * t260 + t314 * t287;
t251 = t284 * pkin(5) - t274 * pkin(11);
t270 = t273 ^ 2;
t204 = -t270 * pkin(5) + t244 * pkin(11) - t284 * t251 + t206;
t318 = sin(qJ(6));
t322 = cos(qJ(6));
t201 = t322 * t203 - t318 * t204;
t241 = t322 * t273 - t318 * t274;
t221 = t241 * qJD(6) + t318 * t244 + t322 * t245;
t242 = t318 * t273 + t322 * t274;
t230 = -t241 * mrSges(7,1) + t242 * mrSges(7,2);
t283 = qJD(6) + t284;
t235 = -t283 * mrSges(7,2) + t241 * mrSges(7,3);
t257 = qJDD(6) + t259;
t197 = m(7) * t201 + t257 * mrSges(7,1) - t221 * mrSges(7,3) - t242 * t230 + t283 * t235;
t202 = t318 * t203 + t322 * t204;
t220 = -t242 * qJD(6) + t322 * t244 - t318 * t245;
t236 = t283 * mrSges(7,1) - t242 * mrSges(7,3);
t198 = m(7) * t202 - t257 * mrSges(7,2) + t220 * mrSges(7,3) + t241 * t230 - t283 * t236;
t189 = t322 * t197 + t318 * t198;
t246 = -t273 * mrSges(6,1) + t274 * mrSges(6,2);
t249 = -t284 * mrSges(6,2) + t273 * mrSges(6,3);
t187 = m(6) * t205 + t259 * mrSges(6,1) - t245 * mrSges(6,3) - t274 * t246 + t284 * t249 + t189;
t250 = t284 * mrSges(6,1) - t274 * mrSges(6,3);
t335 = -t318 * t197 + t322 * t198;
t188 = m(6) * t206 - t259 * mrSges(6,2) + t244 * mrSges(6,3) + t273 * t246 - t284 * t250 + t335;
t185 = -t310 * t187 + t314 * t188;
t267 = t284 * mrSges(5,1) + t285 * mrSges(5,2);
t276 = t301 * mrSges(5,1) - t285 * mrSges(5,3);
t183 = m(5) * t215 - t287 * mrSges(5,2) - t259 * mrSges(5,3) - t284 * t267 - t301 * t276 + t185;
t214 = -t319 * t229 + t232 * t350;
t209 = -t287 * pkin(4) - t300 * qJ(5) + t285 * t266 + qJDD(5) - t214;
t207 = -t244 * pkin(5) - t270 * pkin(11) + t274 * t251 + t209;
t330 = m(7) * t207 - t220 * mrSges(7,1) + t221 * mrSges(7,2) - t241 * t235 + t242 * t236;
t199 = -m(6) * t209 + t244 * mrSges(6,1) - t245 * mrSges(6,2) + t273 * t249 - t274 * t250 - t330;
t275 = -t301 * mrSges(5,2) - t284 * mrSges(5,3);
t193 = m(5) * t214 + t287 * mrSges(5,1) - t260 * mrSges(5,3) - t285 * t267 + t301 * t275 + t199;
t336 = t183 * t350 - t319 * t193;
t172 = m(4) * t234 - t307 * mrSges(4,2) - t293 * mrSges(4,3) - t308 * t288 + t290 * t337 + t336;
t175 = t319 * t183 + t193 * t350;
t247 = -t312 * t264 + t282;
t289 = -t308 * mrSges(4,2) + mrSges(4,3) * t337;
t174 = m(4) * t247 + t293 * mrSges(4,1) + t292 * mrSges(4,2) + (t288 * t320 - t289 * t323) * t340 + t175;
t184 = t314 * t187 + t310 * t188;
t328 = -m(5) * t228 - t259 * mrSges(5,1) - t260 * mrSges(5,2) - t284 * t275 - t285 * t276 - t184;
t180 = m(4) * t233 + t307 * mrSges(4,1) - t292 * mrSges(4,3) + t308 * t289 - t290 * t338 + t328;
t347 = t180 * t323;
t162 = t172 * t343 - t312 * t174 + t316 * t347;
t159 = m(3) * t271 + qJDD(2) * mrSges(3,1) - t325 * mrSges(3,2) + t162;
t167 = t323 * t172 - t320 * t180;
t166 = m(3) * t272 - t325 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t167;
t156 = -t321 * t159 + t324 * t166;
t349 = pkin(8) * t156;
t161 = t172 * t346 + t316 * t174 + t312 * t347;
t160 = m(3) * t286 + t161;
t151 = t159 * t341 - t313 * t160 + t166 * t342;
t222 = Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t283;
t224 = Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t283;
t190 = -mrSges(7,1) * t207 + mrSges(7,3) * t202 + Ifges(7,4) * t221 + Ifges(7,2) * t220 + Ifges(7,6) * t257 - t242 * t222 + t283 * t224;
t223 = Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t283;
t191 = mrSges(7,2) * t207 - mrSges(7,3) * t201 + Ifges(7,1) * t221 + Ifges(7,4) * t220 + Ifges(7,5) * t257 + t241 * t222 - t283 * t223;
t237 = Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * t284;
t239 = Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * t284;
t176 = -mrSges(6,1) * t209 + mrSges(6,3) * t206 + Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t259 - pkin(5) * t330 + pkin(11) * t335 + t322 * t190 + t318 * t191 - t274 * t237 + t284 * t239;
t238 = Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * t284;
t177 = mrSges(6,2) * t209 - mrSges(6,3) * t205 + Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t259 - pkin(11) * t189 - t318 * t190 + t322 * t191 + t273 * t237 - t284 * t238;
t253 = Ifges(5,5) * t285 - Ifges(5,6) * t284 + Ifges(5,3) * t301;
t254 = Ifges(5,4) * t285 - Ifges(5,2) * t284 + Ifges(5,6) * t301;
t163 = mrSges(5,2) * t228 - mrSges(5,3) * t214 + Ifges(5,1) * t260 - Ifges(5,4) * t259 + Ifges(5,5) * t287 - qJ(5) * t184 - t310 * t176 + t314 * t177 - t284 * t253 - t301 * t254;
t255 = Ifges(5,1) * t285 - Ifges(5,4) * t284 + Ifges(5,5) * t301;
t329 = -mrSges(7,1) * t201 + mrSges(7,2) * t202 - Ifges(7,5) * t221 - Ifges(7,6) * t220 - Ifges(7,3) * t257 - t242 * t223 + t241 * t224;
t327 = -mrSges(6,1) * t205 + mrSges(6,2) * t206 - Ifges(6,5) * t245 - Ifges(6,6) * t244 - pkin(5) * t189 - t274 * t238 + t273 * t239 + t329;
t168 = (-Ifges(5,2) - Ifges(6,3)) * t259 + t327 + t301 * t255 - t285 * t253 + Ifges(5,6) * t287 + Ifges(5,4) * t260 - mrSges(5,1) * t228 + mrSges(5,3) * t215 - pkin(4) * t184;
t278 = Ifges(4,3) * t308 + (Ifges(4,5) * t320 + Ifges(4,6) * t323) * t340;
t279 = Ifges(4,6) * t308 + (Ifges(4,4) * t320 + Ifges(4,2) * t323) * t340;
t153 = mrSges(4,2) * t247 - mrSges(4,3) * t233 + Ifges(4,1) * t292 - Ifges(4,4) * t293 + Ifges(4,5) * t307 - pkin(10) * t175 + t163 * t350 - t319 * t168 + t278 * t337 - t308 * t279;
t280 = Ifges(4,5) * t308 + (Ifges(4,1) * t320 + Ifges(4,4) * t323) * t340;
t326 = mrSges(5,1) * t214 - mrSges(5,2) * t215 + Ifges(5,5) * t260 - Ifges(5,6) * t259 + Ifges(5,3) * t287 + pkin(4) * t199 + qJ(5) * t185 + t314 * t176 + t310 * t177 + t285 * t254 + t284 * t255;
t157 = -mrSges(4,1) * t247 + mrSges(4,3) * t234 + Ifges(4,4) * t292 - Ifges(4,2) * t293 + Ifges(4,6) * t307 - pkin(3) * t175 - t278 * t338 + t308 * t280 - t326;
t332 = pkin(9) * t167 + t153 * t320 + t157 * t323;
t152 = Ifges(4,5) * t292 - Ifges(4,6) * t293 + Ifges(4,3) * t307 + mrSges(4,1) * t233 - mrSges(4,2) * t234 + t319 * t163 + t350 * t168 + pkin(3) * t328 + pkin(10) * t336 + (t279 * t320 - t280 * t323) * t340;
t143 = mrSges(3,1) * t271 - mrSges(3,2) * t272 + Ifges(3,3) * qJDD(2) + pkin(2) * t162 + t316 * t152 + t312 * t332;
t145 = -mrSges(3,1) * t286 + mrSges(3,3) * t272 + t325 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t161 - t312 * t152 + t316 * t332;
t147 = mrSges(3,2) * t286 - mrSges(3,3) * t271 + Ifges(3,5) * qJDD(2) - t325 * Ifges(3,6) + t323 * t153 - t320 * t157 + (-t161 * t312 - t162 * t316) * pkin(9);
t331 = mrSges(2,1) * t302 - mrSges(2,2) * t303 + pkin(1) * t151 + t317 * t143 + t145 * t344 + t147 * t345 + t313 * t349;
t154 = m(2) * t303 + t156;
t150 = t317 * t160 + (t159 * t324 + t166 * t321) * t313;
t148 = m(2) * t302 + t151;
t141 = mrSges(2,2) * t309 - mrSges(2,3) * t302 - t321 * t145 + t324 * t147 + (-t150 * t313 - t151 * t317) * pkin(8);
t140 = -mrSges(2,1) * t309 + mrSges(2,3) * t303 - pkin(1) * t150 - t313 * t143 + (t145 * t324 + t147 * t321 + t349) * t317;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t315 * t141 - t311 * t140 - qJ(1) * (t315 * t148 + t311 * t154), t141, t147, t153, t163, t177, t191; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t311 * t141 + t315 * t140 + qJ(1) * (-t311 * t148 + t315 * t154), t140, t145, t157, t168, t176, t190; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t331, t331, t143, t152, t326, Ifges(6,3) * t259 - t327, -t329;];
m_new  = t1;
