% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRP5
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
% Datum: 2019-05-05 10:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:00:57
% EndTime: 2019-05-05 10:01:52
% DurationCPUTime: 31.92s
% Computational Cost: add. (607633->347), mult. (1247141->448), div. (0->0), fcn. (989926->14), ass. (0->148)
t307 = sin(pkin(7));
t314 = sin(qJ(3));
t318 = cos(qJ(3));
t337 = qJD(2) * qJD(3);
t291 = (-qJDD(2) * t318 + t314 * t337) * t307;
t306 = sin(pkin(12));
t309 = cos(pkin(12));
t298 = g(1) * t306 - g(2) * t309;
t299 = -g(1) * t309 - g(2) * t306;
t305 = -g(3) + qJDD(1);
t315 = sin(qJ(2));
t311 = cos(pkin(6));
t319 = cos(qJ(2));
t341 = t311 * t319;
t308 = sin(pkin(6));
t344 = t308 * t319;
t268 = t298 * t341 - t299 * t315 + t305 * t344;
t320 = qJD(2) ^ 2;
t349 = pkin(9) * t307;
t264 = qJDD(2) * pkin(2) + t320 * t349 + t268;
t342 = t311 * t315;
t345 = t308 * t315;
t269 = t298 * t342 + t299 * t319 + t305 * t345;
t265 = -pkin(2) * t320 + qJDD(2) * t349 + t269;
t284 = -t298 * t308 + t305 * t311;
t310 = cos(pkin(7));
t216 = -t314 * t265 + (t264 * t310 + t284 * t307) * t318;
t304 = qJD(2) * t310 + qJD(3);
t313 = sin(qJ(4));
t317 = cos(qJ(4));
t338 = qJD(2) * t307;
t334 = t314 * t338;
t282 = t304 * t317 - t313 * t334;
t290 = (qJDD(2) * t314 + t318 * t337) * t307;
t303 = qJDD(2) * t310 + qJDD(3);
t260 = qJD(4) * t282 + t290 * t317 + t303 * t313;
t283 = t304 * t313 + t317 * t334;
t333 = t318 * t338;
t297 = qJD(4) - t333;
t312 = sin(qJ(5));
t316 = cos(qJ(5));
t271 = -t283 * t312 + t297 * t316;
t285 = qJDD(4) + t291;
t229 = qJD(5) * t271 + t260 * t316 + t285 * t312;
t272 = t283 * t316 + t297 * t312;
t240 = -mrSges(7,1) * t271 + mrSges(7,2) * t272;
t343 = t310 * t314;
t346 = t307 * t314;
t217 = t264 * t343 + t265 * t318 + t284 * t346;
t289 = (-pkin(3) * t318 - pkin(10) * t314) * t338;
t302 = t304 ^ 2;
t211 = -pkin(3) * t302 + pkin(10) * t303 + t289 * t333 + t217;
t280 = t310 * t284;
t213 = pkin(3) * t291 - pkin(10) * t290 + t280 + (-t264 + (pkin(3) * t314 - pkin(10) * t318) * t304 * qJD(2)) * t307;
t207 = t211 * t317 + t213 * t313;
t267 = -pkin(4) * t282 - pkin(11) * t283;
t296 = t297 ^ 2;
t202 = -pkin(4) * t296 + pkin(11) * t285 + t267 * t282 + t207;
t210 = -pkin(3) * t303 - pkin(10) * t302 + t289 * t334 - t216;
t259 = -qJD(4) * t283 - t290 * t313 + t303 * t317;
t205 = (-t282 * t297 - t260) * pkin(11) + (t283 * t297 - t259) * pkin(4) + t210;
t196 = -t202 * t312 + t205 * t316;
t257 = qJDD(5) - t259;
t281 = qJD(5) - t282;
t192 = -0.2e1 * qJD(6) * t272 + (t271 * t281 - t229) * qJ(6) + (t271 * t272 + t257) * pkin(5) + t196;
t244 = -mrSges(7,2) * t281 + mrSges(7,3) * t271;
t336 = m(7) * t192 + mrSges(7,1) * t257 + t244 * t281;
t189 = -mrSges(7,3) * t229 - t240 * t272 + t336;
t197 = t202 * t316 + t205 * t312;
t228 = -qJD(5) * t272 - t260 * t312 + t285 * t316;
t234 = Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t281;
t235 = Ifges(7,1) * t272 + Ifges(7,4) * t271 + Ifges(7,5) * t281;
t236 = Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t281;
t246 = pkin(5) * t281 - qJ(6) * t272;
t270 = t271 ^ 2;
t195 = -pkin(5) * t270 + qJ(6) * t228 + 0.2e1 * qJD(6) * t271 - t246 * t281 + t197;
t233 = Ifges(7,4) * t272 + Ifges(7,2) * t271 + Ifges(7,6) * t281;
t326 = -mrSges(7,1) * t192 + mrSges(7,2) * t195 - Ifges(7,5) * t229 - Ifges(7,6) * t228 - Ifges(7,3) * t257 - t233 * t272;
t351 = mrSges(6,1) * t196 - mrSges(6,2) * t197 + Ifges(6,5) * t229 + Ifges(6,6) * t228 + Ifges(6,3) * t257 + pkin(5) * t189 + t272 * t234 - (t236 + t235) * t271 - t326;
t286 = mrSges(4,1) * t304 - mrSges(4,3) * t334;
t288 = (-mrSges(4,1) * t318 + mrSges(4,2) * t314) * t338;
t241 = -mrSges(6,1) * t271 + mrSges(6,2) * t272;
t245 = -mrSges(6,2) * t281 + mrSges(6,3) * t271;
t183 = m(6) * t196 + mrSges(6,1) * t257 + t245 * t281 + (-t240 - t241) * t272 + (-mrSges(6,3) - mrSges(7,3)) * t229 + t336;
t335 = m(7) * t195 + mrSges(7,3) * t228 + t240 * t271;
t247 = mrSges(7,1) * t281 - mrSges(7,3) * t272;
t339 = -mrSges(6,1) * t281 + mrSges(6,3) * t272 - t247;
t348 = -mrSges(6,2) - mrSges(7,2);
t185 = m(6) * t197 + mrSges(6,3) * t228 + t241 * t271 + t257 * t348 + t281 * t339 + t335;
t182 = -t183 * t312 + t185 * t316;
t266 = -mrSges(5,1) * t282 + mrSges(5,2) * t283;
t274 = mrSges(5,1) * t297 - mrSges(5,3) * t283;
t179 = m(5) * t207 - mrSges(5,2) * t285 + mrSges(5,3) * t259 + t266 * t282 - t274 * t297 + t182;
t206 = -t211 * t313 + t213 * t317;
t201 = -pkin(4) * t285 - pkin(11) * t296 + t267 * t283 - t206;
t199 = -pkin(5) * t228 - qJ(6) * t270 + t246 * t272 + qJDD(6) + t201;
t331 = -m(7) * t199 + mrSges(7,1) * t228 + t244 * t271;
t188 = -m(6) * t201 + mrSges(6,1) * t228 + t229 * t348 + t245 * t271 + t272 * t339 + t331;
t273 = -mrSges(5,2) * t297 + mrSges(5,3) * t282;
t187 = m(5) * t206 + mrSges(5,1) * t285 - mrSges(5,3) * t260 - t266 * t283 + t273 * t297 + t188;
t332 = t179 * t317 - t187 * t313;
t169 = m(4) * t217 - mrSges(4,2) * t303 - mrSges(4,3) * t291 - t286 * t304 + t288 * t333 + t332;
t172 = t179 * t313 + t187 * t317;
t242 = -t264 * t307 + t280;
t287 = -mrSges(4,2) * t304 + mrSges(4,3) * t333;
t171 = m(4) * t242 + mrSges(4,1) * t291 + mrSges(4,2) * t290 + (t286 * t314 - t287 * t318) * t338 + t172;
t181 = t183 * t316 + t185 * t312;
t323 = -m(5) * t210 + mrSges(5,1) * t259 - mrSges(5,2) * t260 + t273 * t282 - t274 * t283 - t181;
t176 = m(4) * t216 + mrSges(4,1) * t303 - mrSges(4,3) * t290 + t287 * t304 - t288 * t334 + t323;
t347 = t176 * t318;
t159 = t169 * t343 - t171 * t307 + t310 * t347;
t156 = m(3) * t268 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t320 + t159;
t164 = t169 * t318 - t176 * t314;
t163 = m(3) * t269 - mrSges(3,1) * t320 - qJDD(2) * mrSges(3,2) + t164;
t153 = -t156 * t315 + t163 * t319;
t350 = pkin(8) * t153;
t158 = t169 * t346 + t171 * t310 + t307 * t347;
t157 = m(3) * t284 + t158;
t148 = t156 * t341 - t157 * t308 + t163 * t342;
t231 = Ifges(7,5) * t272 + Ifges(7,6) * t271 + Ifges(7,3) * t281;
t232 = Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t281;
t327 = -mrSges(7,1) * t199 + mrSges(7,3) * t195 + Ifges(7,4) * t229 + Ifges(7,2) * t228 + Ifges(7,6) * t257 + t235 * t281;
t173 = Ifges(6,4) * t229 + Ifges(6,2) * t228 + Ifges(6,6) * t257 + t281 * t236 - mrSges(6,1) * t201 + mrSges(6,3) * t197 - pkin(5) * (mrSges(7,2) * t229 - t331) + qJ(6) * (-mrSges(7,2) * t257 - t247 * t281 + t335) + (-pkin(5) * t247 - t231 - t232) * t272 + t327;
t325 = mrSges(7,2) * t199 - mrSges(7,3) * t192 + Ifges(7,1) * t229 + Ifges(7,4) * t228 + Ifges(7,5) * t257 + t231 * t271;
t180 = mrSges(6,2) * t201 - mrSges(6,3) * t196 + Ifges(6,1) * t229 + Ifges(6,4) * t228 + Ifges(6,5) * t257 - qJ(6) * t189 + t232 * t271 + (-t233 - t234) * t281 + t325;
t253 = Ifges(5,5) * t283 + Ifges(5,6) * t282 + Ifges(5,3) * t297;
t254 = Ifges(5,4) * t283 + Ifges(5,2) * t282 + Ifges(5,6) * t297;
t160 = mrSges(5,2) * t210 - mrSges(5,3) * t206 + Ifges(5,1) * t260 + Ifges(5,4) * t259 + Ifges(5,5) * t285 - pkin(11) * t181 - t173 * t312 + t180 * t316 + t253 * t282 - t254 * t297;
t255 = Ifges(5,1) * t283 + Ifges(5,4) * t282 + Ifges(5,5) * t297;
t165 = -mrSges(5,1) * t210 + mrSges(5,3) * t207 + Ifges(5,4) * t260 + Ifges(5,2) * t259 + Ifges(5,6) * t285 - pkin(4) * t181 - t283 * t253 + t297 * t255 - t351;
t276 = Ifges(4,3) * t304 + (Ifges(4,5) * t314 + Ifges(4,6) * t318) * t338;
t277 = Ifges(4,6) * t304 + (Ifges(4,4) * t314 + Ifges(4,2) * t318) * t338;
t150 = mrSges(4,2) * t242 - mrSges(4,3) * t216 + Ifges(4,1) * t290 - Ifges(4,4) * t291 + Ifges(4,5) * t303 - pkin(10) * t172 + t160 * t317 - t165 * t313 + t276 * t333 - t277 * t304;
t278 = Ifges(4,5) * t304 + (Ifges(4,1) * t314 + Ifges(4,4) * t318) * t338;
t321 = mrSges(5,1) * t206 - mrSges(5,2) * t207 + Ifges(5,5) * t260 + Ifges(5,6) * t259 + Ifges(5,3) * t285 + pkin(4) * t188 + pkin(11) * t182 + t173 * t316 + t180 * t312 + t254 * t283 - t255 * t282;
t154 = -mrSges(4,1) * t242 + mrSges(4,3) * t217 + Ifges(4,4) * t290 - Ifges(4,2) * t291 + Ifges(4,6) * t303 - pkin(3) * t172 - t276 * t334 + t278 * t304 - t321;
t328 = pkin(9) * t164 + t150 * t314 + t154 * t318;
t149 = Ifges(4,5) * t290 - Ifges(4,6) * t291 + Ifges(4,3) * t303 + mrSges(4,1) * t216 - mrSges(4,2) * t217 + t313 * t160 + t317 * t165 + pkin(3) * t323 + pkin(10) * t332 + (t277 * t314 - t278 * t318) * t338;
t140 = mrSges(3,1) * t268 - mrSges(3,2) * t269 + Ifges(3,3) * qJDD(2) + pkin(2) * t159 + t149 * t310 + t307 * t328;
t142 = -mrSges(3,1) * t284 + mrSges(3,3) * t269 + Ifges(3,5) * t320 + Ifges(3,6) * qJDD(2) - pkin(2) * t158 - t149 * t307 + t310 * t328;
t144 = mrSges(3,2) * t284 - mrSges(3,3) * t268 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t320 + t150 * t318 - t154 * t314 + (-t158 * t307 - t159 * t310) * pkin(9);
t324 = mrSges(2,1) * t298 - mrSges(2,2) * t299 + pkin(1) * t148 + t140 * t311 + t142 * t344 + t144 * t345 + t308 * t350;
t151 = m(2) * t299 + t153;
t147 = t157 * t311 + (t156 * t319 + t163 * t315) * t308;
t145 = m(2) * t298 + t148;
t138 = mrSges(2,2) * t305 - mrSges(2,3) * t298 - t142 * t315 + t144 * t319 + (-t147 * t308 - t148 * t311) * pkin(8);
t137 = -mrSges(2,1) * t305 + mrSges(2,3) * t299 - pkin(1) * t147 - t140 * t308 + (t142 * t319 + t144 * t315 + t350) * t311;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t138 - t306 * t137 - qJ(1) * (t145 * t309 + t151 * t306), t138, t144, t150, t160, t180, -t233 * t281 + t325; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t306 * t138 + t309 * t137 + qJ(1) * (-t145 * t306 + t151 * t309), t137, t142, t154, t165, t173, -t272 * t231 + t327; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t324, t324, t140, t149, t321, t351, -t271 * t235 - t326;];
m_new  = t1;
