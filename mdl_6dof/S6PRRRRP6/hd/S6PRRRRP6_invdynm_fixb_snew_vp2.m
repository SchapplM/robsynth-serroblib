% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:11:15
% EndTime: 2019-05-05 10:12:06
% DurationCPUTime: 30.86s
% Computational Cost: add. (592731->348), mult. (1213465->448), div. (0->0), fcn. (961700->14), ass. (0->147)
t305 = sin(pkin(7));
t312 = sin(qJ(3));
t315 = cos(qJ(3));
t333 = qJD(2) * qJD(3);
t288 = (-qJDD(2) * t315 + t312 * t333) * t305;
t304 = sin(pkin(12));
t307 = cos(pkin(12));
t296 = g(1) * t304 - g(2) * t307;
t297 = -g(1) * t307 - g(2) * t304;
t303 = -g(3) + qJDD(1);
t313 = sin(qJ(2));
t309 = cos(pkin(6));
t316 = cos(qJ(2));
t338 = t309 * t316;
t306 = sin(pkin(6));
t341 = t306 * t316;
t264 = t296 * t338 - t297 * t313 + t303 * t341;
t317 = qJD(2) ^ 2;
t346 = pkin(9) * t305;
t260 = qJDD(2) * pkin(2) + t317 * t346 + t264;
t339 = t309 * t313;
t342 = t306 * t313;
t265 = t296 * t339 + t297 * t316 + t303 * t342;
t261 = -pkin(2) * t317 + qJDD(2) * t346 + t265;
t281 = -t296 * t306 + t303 * t309;
t308 = cos(pkin(7));
t215 = -t312 * t261 + (t260 * t308 + t281 * t305) * t315;
t340 = t308 * t312;
t343 = t305 * t312;
t216 = t260 * t340 + t261 * t315 + t281 * t343;
t334 = qJD(2) * t305;
t286 = (-pkin(3) * t315 - pkin(10) * t312) * t334;
t302 = qJD(2) * t308 + qJD(3);
t300 = t302 ^ 2;
t301 = qJDD(2) * t308 + qJDD(3);
t330 = t315 * t334;
t211 = -pkin(3) * t300 + pkin(10) * t301 + t286 * t330 + t216;
t275 = t308 * t281;
t287 = (qJDD(2) * t312 + t315 * t333) * t305;
t213 = t288 * pkin(3) - t287 * pkin(10) + t275 + (-t260 + (pkin(3) * t312 - pkin(10) * t315) * t302 * qJD(2)) * t305;
t311 = sin(qJ(4));
t314 = cos(qJ(4));
t207 = t211 * t314 + t213 * t311;
t331 = t312 * t334;
t279 = t302 * t314 - t311 * t331;
t280 = t302 * t311 + t314 * t331;
t263 = -pkin(4) * t279 - pkin(11) * t280;
t282 = qJDD(4) + t288;
t295 = qJD(4) - t330;
t294 = t295 ^ 2;
t203 = -pkin(4) * t294 + pkin(11) * t282 + t263 * t279 + t207;
t210 = -t301 * pkin(3) - t300 * pkin(10) + t286 * t331 - t215;
t255 = -qJD(4) * t280 - t287 * t311 + t301 * t314;
t256 = qJD(4) * t279 + t287 * t314 + t301 * t311;
t205 = (-t279 * t295 - t256) * pkin(11) + (t280 * t295 - t255) * pkin(4) + t210;
t310 = sin(qJ(5));
t348 = cos(qJ(5));
t199 = -t203 * t310 + t205 * t348;
t200 = t203 * t348 + t205 * t310;
t267 = t280 * t348 + t295 * t310;
t224 = qJD(5) * t267 + t256 * t310 - t282 * t348;
t266 = t280 * t310 - t295 * t348;
t225 = -qJD(5) * t266 + t256 * t348 + t282 * t310;
t277 = qJD(5) - t279;
t226 = Ifges(7,5) * t267 + Ifges(7,6) * t277 + Ifges(7,3) * t266;
t229 = Ifges(6,4) * t267 - Ifges(6,2) * t266 + Ifges(6,6) * t277;
t231 = Ifges(6,1) * t267 - Ifges(6,4) * t266 + Ifges(6,5) * t277;
t237 = mrSges(7,1) * t266 - mrSges(7,3) * t267;
t253 = qJDD(5) - t255;
t236 = pkin(5) * t266 - qJ(6) * t267;
t276 = t277 ^ 2;
t195 = -pkin(5) * t276 + qJ(6) * t253 + 0.2e1 * qJD(6) * t277 - t236 * t266 + t200;
t197 = -pkin(5) * t253 - qJ(6) * t276 + t236 * t267 + qJDD(6) - t199;
t230 = Ifges(7,1) * t267 + Ifges(7,4) * t277 + Ifges(7,5) * t266;
t323 = mrSges(7,1) * t197 - mrSges(7,3) * t195 - Ifges(7,4) * t225 - Ifges(7,2) * t253 - Ifges(7,6) * t224 - t230 * t266;
t241 = -mrSges(7,2) * t266 + mrSges(7,3) * t277;
t328 = -m(7) * t197 + mrSges(7,1) * t253 + t241 * t277;
t244 = -mrSges(7,1) * t277 + mrSges(7,2) * t267;
t332 = m(7) * t195 + mrSges(7,3) * t253 + t244 * t277;
t349 = -(-t229 + t226) * t267 + mrSges(6,1) * t199 - mrSges(6,2) * t200 + Ifges(6,5) * t225 - Ifges(6,6) * t224 + Ifges(6,3) * t253 + pkin(5) * (-t225 * mrSges(7,2) - t267 * t237 + t328) + qJ(6) * (-t224 * mrSges(7,2) - t266 * t237 + t332) + t266 * t231 - t323;
t283 = mrSges(4,1) * t302 - mrSges(4,3) * t331;
t285 = (-mrSges(4,1) * t315 + mrSges(4,2) * t312) * t334;
t243 = mrSges(6,1) * t277 - mrSges(6,3) * t267;
t335 = -mrSges(6,1) * t266 - mrSges(6,2) * t267 - t237;
t345 = -mrSges(6,3) - mrSges(7,2);
t187 = m(6) * t200 - t253 * mrSges(6,2) + t224 * t345 - t277 * t243 + t266 * t335 + t332;
t242 = -mrSges(6,2) * t277 - mrSges(6,3) * t266;
t188 = m(6) * t199 + t253 * mrSges(6,1) + t225 * t345 + t277 * t242 + t267 * t335 + t328;
t183 = t187 * t348 - t188 * t310;
t262 = -mrSges(5,1) * t279 + mrSges(5,2) * t280;
t269 = mrSges(5,1) * t295 - mrSges(5,3) * t280;
t179 = m(5) * t207 - mrSges(5,2) * t282 + mrSges(5,3) * t255 + t262 * t279 - t269 * t295 + t183;
t206 = -t211 * t311 + t314 * t213;
t202 = -t282 * pkin(4) - t294 * pkin(11) + t263 * t280 - t206;
t198 = -0.2e1 * qJD(6) * t267 + (t266 * t277 - t225) * qJ(6) + (t267 * t277 + t224) * pkin(5) + t202;
t192 = m(7) * t198 + mrSges(7,1) * t224 - mrSges(7,3) * t225 + t241 * t266 - t244 * t267;
t189 = -m(6) * t202 - mrSges(6,1) * t224 - mrSges(6,2) * t225 - t242 * t266 - t243 * t267 - t192;
t268 = -mrSges(5,2) * t295 + mrSges(5,3) * t279;
t185 = m(5) * t206 + mrSges(5,1) * t282 - mrSges(5,3) * t256 - t262 * t280 + t268 * t295 + t189;
t329 = t179 * t314 - t185 * t311;
t170 = m(4) * t216 - mrSges(4,2) * t301 - mrSges(4,3) * t288 - t283 * t302 + t285 * t330 + t329;
t173 = t179 * t311 + t185 * t314;
t239 = -t305 * t260 + t275;
t284 = -mrSges(4,2) * t302 + mrSges(4,3) * t330;
t172 = m(4) * t239 + t288 * mrSges(4,1) + t287 * mrSges(4,2) + (t283 * t312 - t284 * t315) * t334 + t173;
t182 = t187 * t310 + t188 * t348;
t320 = -m(5) * t210 + mrSges(5,1) * t255 - mrSges(5,2) * t256 + t268 * t279 - t269 * t280 - t182;
t176 = m(4) * t215 + mrSges(4,1) * t301 - mrSges(4,3) * t287 + t284 * t302 - t285 * t331 + t320;
t344 = t176 * t315;
t160 = t170 * t340 - t172 * t305 + t308 * t344;
t157 = m(3) * t264 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t317 + t160;
t165 = t170 * t315 - t176 * t312;
t164 = m(3) * t265 - mrSges(3,1) * t317 - qJDD(2) * mrSges(3,2) + t165;
t154 = -t157 * t313 + t164 * t316;
t347 = pkin(8) * t154;
t228 = Ifges(7,4) * t267 + Ifges(7,2) * t277 + Ifges(7,6) * t266;
t337 = -Ifges(6,5) * t267 + Ifges(6,6) * t266 - Ifges(6,3) * t277 - t228;
t159 = t170 * t343 + t172 * t308 + t305 * t344;
t158 = m(3) * t281 + t159;
t149 = t157 * t338 - t158 * t306 + t164 * t339;
t327 = -mrSges(7,1) * t198 + mrSges(7,2) * t195;
t180 = -mrSges(6,1) * t202 + mrSges(6,3) * t200 - pkin(5) * t192 + (t230 + t231) * t277 + t337 * t267 + (Ifges(6,6) - Ifges(7,6)) * t253 + (Ifges(6,4) - Ifges(7,5)) * t225 + (-Ifges(6,2) - Ifges(7,3)) * t224 + t327;
t322 = mrSges(7,2) * t197 - mrSges(7,3) * t198 + Ifges(7,1) * t225 + Ifges(7,4) * t253 + Ifges(7,5) * t224 + t226 * t277;
t181 = mrSges(6,2) * t202 - mrSges(6,3) * t199 + Ifges(6,1) * t225 - Ifges(6,4) * t224 + Ifges(6,5) * t253 - qJ(6) * t192 - t277 * t229 + t266 * t337 + t322;
t249 = Ifges(5,5) * t280 + Ifges(5,6) * t279 + Ifges(5,3) * t295;
t250 = Ifges(5,4) * t280 + Ifges(5,2) * t279 + Ifges(5,6) * t295;
t161 = mrSges(5,2) * t210 - mrSges(5,3) * t206 + Ifges(5,1) * t256 + Ifges(5,4) * t255 + Ifges(5,5) * t282 - pkin(11) * t182 - t180 * t310 + t181 * t348 + t249 * t279 - t250 * t295;
t251 = Ifges(5,1) * t280 + Ifges(5,4) * t279 + Ifges(5,5) * t295;
t166 = -mrSges(5,1) * t210 + mrSges(5,3) * t207 + Ifges(5,4) * t256 + Ifges(5,2) * t255 + Ifges(5,6) * t282 - pkin(4) * t182 - t280 * t249 + t295 * t251 - t349;
t271 = Ifges(4,3) * t302 + (Ifges(4,5) * t312 + Ifges(4,6) * t315) * t334;
t272 = Ifges(4,6) * t302 + (Ifges(4,4) * t312 + Ifges(4,2) * t315) * t334;
t151 = mrSges(4,2) * t239 - mrSges(4,3) * t215 + Ifges(4,1) * t287 - Ifges(4,4) * t288 + Ifges(4,5) * t301 - pkin(10) * t173 + t161 * t314 - t166 * t311 + t271 * t330 - t272 * t302;
t273 = Ifges(4,5) * t302 + (Ifges(4,1) * t312 + Ifges(4,4) * t315) * t334;
t318 = mrSges(5,1) * t206 - mrSges(5,2) * t207 + Ifges(5,5) * t256 + Ifges(5,6) * t255 + Ifges(5,3) * t282 + pkin(4) * t189 + pkin(11) * t183 + t180 * t348 + t181 * t310 + t250 * t280 - t251 * t279;
t155 = -mrSges(4,1) * t239 + mrSges(4,3) * t216 + Ifges(4,4) * t287 - Ifges(4,2) * t288 + Ifges(4,6) * t301 - pkin(3) * t173 - t271 * t331 + t273 * t302 - t318;
t324 = pkin(9) * t165 + t151 * t312 + t155 * t315;
t150 = Ifges(4,5) * t287 - Ifges(4,6) * t288 + Ifges(4,3) * t301 + mrSges(4,1) * t215 - mrSges(4,2) * t216 + t311 * t161 + t314 * t166 + pkin(3) * t320 + pkin(10) * t329 + (t272 * t312 - t273 * t315) * t334;
t141 = mrSges(3,1) * t264 - mrSges(3,2) * t265 + Ifges(3,3) * qJDD(2) + pkin(2) * t160 + t308 * t150 + t324 * t305;
t143 = -mrSges(3,1) * t281 + mrSges(3,3) * t265 + t317 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t159 - t305 * t150 + t308 * t324;
t145 = mrSges(3,2) * t281 - mrSges(3,3) * t264 + Ifges(3,5) * qJDD(2) - t317 * Ifges(3,6) + t315 * t151 - t312 * t155 + (-t159 * t305 - t160 * t308) * pkin(9);
t321 = mrSges(2,1) * t296 - mrSges(2,2) * t297 + pkin(1) * t149 + t141 * t309 + t143 * t341 + t145 * t342 + t306 * t347;
t152 = m(2) * t297 + t154;
t148 = t309 * t158 + (t157 * t316 + t164 * t313) * t306;
t146 = m(2) * t296 + t149;
t139 = mrSges(2,2) * t303 - mrSges(2,3) * t296 - t313 * t143 + t316 * t145 + (-t148 * t306 - t149 * t309) * pkin(8);
t138 = -mrSges(2,1) * t303 + mrSges(2,3) * t297 - pkin(1) * t148 - t306 * t141 + (t143 * t316 + t145 * t313 + t347) * t309;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t307 * t139 - t304 * t138 - qJ(1) * (t146 * t307 + t152 * t304), t139, t145, t151, t161, t181, -t228 * t266 + t322; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t304 * t139 + t307 * t138 + qJ(1) * (-t146 * t304 + t152 * t307), t138, t143, t155, t166, t180, -t267 * t226 - t323; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t321, t321, t141, t150, t318, t349, Ifges(7,5) * t225 + Ifges(7,6) * t253 + Ifges(7,3) * t224 + t267 * t228 - t277 * t230 - t327;];
m_new  = t1;
