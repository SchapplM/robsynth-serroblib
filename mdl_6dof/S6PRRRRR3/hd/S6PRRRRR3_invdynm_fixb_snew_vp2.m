% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:01:24
% EndTime: 2019-05-05 11:02:22
% DurationCPUTime: 38.38s
% Computational Cost: add. (737915->340), mult. (1455884->434), div. (0->0), fcn. (1073975->14), ass. (0->146)
t298 = sin(pkin(12));
t300 = cos(pkin(12));
t288 = g(1) * t298 - g(2) * t300;
t289 = -g(1) * t300 - g(2) * t298;
t297 = -g(3) + qJDD(1);
t311 = cos(qJ(2));
t301 = cos(pkin(6));
t306 = sin(qJ(2));
t331 = t301 * t306;
t299 = sin(pkin(6));
t332 = t299 * t306;
t248 = t288 * t331 + t311 * t289 + t297 * t332;
t313 = qJD(2) ^ 2;
t243 = -pkin(2) * t313 + qJDD(2) * pkin(8) + t248;
t265 = -t288 * t299 + t297 * t301;
t305 = sin(qJ(3));
t310 = cos(qJ(3));
t236 = t310 * t243 + t305 * t265;
t284 = (-pkin(3) * t310 - pkin(9) * t305) * qJD(2);
t312 = qJD(3) ^ 2;
t329 = qJD(2) * t310;
t225 = -pkin(3) * t312 + qJDD(3) * pkin(9) + t284 * t329 + t236;
t247 = -t306 * t289 + (t288 * t301 + t297 * t299) * t311;
t242 = -qJDD(2) * pkin(2) - t313 * pkin(8) - t247;
t328 = qJD(2) * qJD(3);
t327 = t310 * t328;
t285 = qJDD(2) * t305 + t327;
t296 = t305 * t328;
t286 = qJDD(2) * t310 - t296;
t230 = (-t285 - t327) * pkin(9) + (-t286 + t296) * pkin(3) + t242;
t304 = sin(qJ(4));
t309 = cos(qJ(4));
t208 = -t304 * t225 + t309 * t230;
t330 = qJD(2) * t305;
t281 = qJD(3) * t309 - t304 * t330;
t256 = qJD(4) * t281 + qJDD(3) * t304 + t285 * t309;
t278 = qJDD(4) - t286;
t282 = qJD(3) * t304 + t309 * t330;
t295 = qJD(4) - t329;
t204 = (t281 * t295 - t256) * pkin(10) + (t281 * t282 + t278) * pkin(4) + t208;
t209 = t309 * t225 + t304 * t230;
t255 = -qJD(4) * t282 + qJDD(3) * t309 - t285 * t304;
t264 = pkin(4) * t295 - pkin(10) * t282;
t277 = t281 ^ 2;
t206 = -pkin(4) * t277 + pkin(10) * t255 - t264 * t295 + t209;
t303 = sin(qJ(5));
t308 = cos(qJ(5));
t192 = t308 * t204 - t303 * t206;
t258 = t281 * t308 - t282 * t303;
t222 = qJD(5) * t258 + t255 * t303 + t256 * t308;
t259 = t281 * t303 + t282 * t308;
t274 = qJDD(5) + t278;
t294 = qJD(5) + t295;
t189 = (t258 * t294 - t222) * pkin(11) + (t258 * t259 + t274) * pkin(5) + t192;
t193 = t303 * t204 + t308 * t206;
t221 = -qJD(5) * t259 + t255 * t308 - t256 * t303;
t246 = pkin(5) * t294 - pkin(11) * t259;
t257 = t258 ^ 2;
t190 = -pkin(5) * t257 + pkin(11) * t221 - t246 * t294 + t193;
t302 = sin(qJ(6));
t307 = cos(qJ(6));
t187 = t189 * t307 - t190 * t302;
t237 = t258 * t307 - t259 * t302;
t201 = qJD(6) * t237 + t221 * t302 + t222 * t307;
t238 = t258 * t302 + t259 * t307;
t216 = -mrSges(7,1) * t237 + mrSges(7,2) * t238;
t290 = qJD(6) + t294;
t228 = -mrSges(7,2) * t290 + mrSges(7,3) * t237;
t269 = qJDD(6) + t274;
t183 = m(7) * t187 + mrSges(7,1) * t269 - mrSges(7,3) * t201 - t216 * t238 + t228 * t290;
t188 = t189 * t302 + t190 * t307;
t200 = -qJD(6) * t238 + t221 * t307 - t222 * t302;
t229 = mrSges(7,1) * t290 - mrSges(7,3) * t238;
t184 = m(7) * t188 - mrSges(7,2) * t269 + mrSges(7,3) * t200 + t216 * t237 - t229 * t290;
t175 = t307 * t183 + t302 * t184;
t239 = -mrSges(6,1) * t258 + mrSges(6,2) * t259;
t244 = -mrSges(6,2) * t294 + mrSges(6,3) * t258;
t172 = m(6) * t192 + mrSges(6,1) * t274 - mrSges(6,3) * t222 - t239 * t259 + t244 * t294 + t175;
t245 = mrSges(6,1) * t294 - mrSges(6,3) * t259;
t324 = -t183 * t302 + t307 * t184;
t173 = m(6) * t193 - mrSges(6,2) * t274 + mrSges(6,3) * t221 + t239 * t258 - t245 * t294 + t324;
t168 = t308 * t172 + t303 * t173;
t260 = -mrSges(5,1) * t281 + mrSges(5,2) * t282;
t262 = -mrSges(5,2) * t295 + mrSges(5,3) * t281;
t166 = m(5) * t208 + mrSges(5,1) * t278 - mrSges(5,3) * t256 - t260 * t282 + t262 * t295 + t168;
t263 = mrSges(5,1) * t295 - mrSges(5,3) * t282;
t325 = -t172 * t303 + t308 * t173;
t167 = m(5) * t209 - mrSges(5,2) * t278 + mrSges(5,3) * t255 + t260 * t281 - t263 * t295 + t325;
t162 = -t166 * t304 + t309 * t167;
t283 = (-mrSges(4,1) * t310 + mrSges(4,2) * t305) * qJD(2);
t291 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t330;
t160 = m(4) * t236 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t286 - qJD(3) * t291 + t283 * t329 + t162;
t235 = -t305 * t243 + t265 * t310;
t224 = -qJDD(3) * pkin(3) - pkin(9) * t312 + t284 * t330 - t235;
t210 = -pkin(4) * t255 - pkin(10) * t277 + t282 * t264 + t224;
t195 = -pkin(5) * t221 - pkin(11) * t257 + t246 * t259 + t210;
t323 = m(7) * t195 - t200 * mrSges(7,1) + t201 * mrSges(7,2) - t237 * t228 + t238 * t229;
t318 = m(6) * t210 - t221 * mrSges(6,1) + t222 * mrSges(6,2) - t258 * t244 + t259 * t245 + t323;
t185 = -m(5) * t224 + t255 * mrSges(5,1) - t256 * mrSges(5,2) + t281 * t262 - t282 * t263 - t318;
t292 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t329;
t179 = m(4) * t235 + qJDD(3) * mrSges(4,1) - t285 * mrSges(4,3) + qJD(3) * t292 - t283 * t330 + t185;
t155 = t305 * t160 + t310 * t179;
t211 = Ifges(7,5) * t238 + Ifges(7,6) * t237 + Ifges(7,3) * t290;
t213 = Ifges(7,1) * t238 + Ifges(7,4) * t237 + Ifges(7,5) * t290;
t176 = -mrSges(7,1) * t195 + mrSges(7,3) * t188 + Ifges(7,4) * t201 + Ifges(7,2) * t200 + Ifges(7,6) * t269 - t211 * t238 + t213 * t290;
t212 = Ifges(7,4) * t238 + Ifges(7,2) * t237 + Ifges(7,6) * t290;
t177 = mrSges(7,2) * t195 - mrSges(7,3) * t187 + Ifges(7,1) * t201 + Ifges(7,4) * t200 + Ifges(7,5) * t269 + t211 * t237 - t212 * t290;
t231 = Ifges(6,5) * t259 + Ifges(6,6) * t258 + Ifges(6,3) * t294;
t233 = Ifges(6,1) * t259 + Ifges(6,4) * t258 + Ifges(6,5) * t294;
t163 = -mrSges(6,1) * t210 + mrSges(6,3) * t193 + Ifges(6,4) * t222 + Ifges(6,2) * t221 + Ifges(6,6) * t274 - pkin(5) * t323 + pkin(11) * t324 + t307 * t176 + t302 * t177 - t259 * t231 + t294 * t233;
t232 = Ifges(6,4) * t259 + Ifges(6,2) * t258 + Ifges(6,6) * t294;
t164 = mrSges(6,2) * t210 - mrSges(6,3) * t192 + Ifges(6,1) * t222 + Ifges(6,4) * t221 + Ifges(6,5) * t274 - pkin(11) * t175 - t176 * t302 + t177 * t307 + t231 * t258 - t232 * t294;
t249 = Ifges(5,5) * t282 + Ifges(5,6) * t281 + Ifges(5,3) * t295;
t251 = Ifges(5,1) * t282 + Ifges(5,4) * t281 + Ifges(5,5) * t295;
t149 = -mrSges(5,1) * t224 + mrSges(5,3) * t209 + Ifges(5,4) * t256 + Ifges(5,2) * t255 + Ifges(5,6) * t278 - pkin(4) * t318 + pkin(10) * t325 + t308 * t163 + t303 * t164 - t282 * t249 + t295 * t251;
t250 = Ifges(5,4) * t282 + Ifges(5,2) * t281 + Ifges(5,6) * t295;
t150 = mrSges(5,2) * t224 - mrSges(5,3) * t208 + Ifges(5,1) * t256 + Ifges(5,4) * t255 + Ifges(5,5) * t278 - pkin(10) * t168 - t163 * t303 + t164 * t308 + t249 * t281 - t250 * t295;
t272 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t305 + Ifges(4,2) * t310) * qJD(2);
t273 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t305 + Ifges(4,4) * t310) * qJD(2);
t336 = mrSges(4,1) * t235 - mrSges(4,2) * t236 + Ifges(4,5) * t285 + Ifges(4,6) * t286 + Ifges(4,3) * qJDD(3) + pkin(3) * t185 + pkin(9) * t162 + t309 * t149 + t304 * t150 + (t272 * t305 - t273 * t310) * qJD(2);
t139 = -mrSges(3,1) * t265 + mrSges(3,3) * t248 + t313 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t155 - t336;
t326 = t310 * t160 - t179 * t305;
t153 = m(3) * t248 - mrSges(3,1) * t313 - qJDD(2) * mrSges(3,2) + t326;
t161 = t166 * t309 + t167 * t304;
t317 = -m(4) * t242 + t286 * mrSges(4,1) - mrSges(4,2) * t285 - t291 * t330 + t292 * t329 - t161;
t157 = m(3) * t247 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t313 + t317;
t147 = t311 * t153 - t157 * t306;
t337 = pkin(7) * t147 + t139 * t311;
t333 = t157 * t311;
t154 = m(3) * t265 + t155;
t144 = t153 * t331 - t154 * t299 + t301 * t333;
t271 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t305 + Ifges(4,6) * t310) * qJD(2);
t140 = mrSges(4,2) * t242 - mrSges(4,3) * t235 + Ifges(4,1) * t285 + Ifges(4,4) * t286 + Ifges(4,5) * qJDD(3) - pkin(9) * t161 - qJD(3) * t272 - t149 * t304 + t150 * t309 + t271 * t329;
t319 = -mrSges(7,1) * t187 + mrSges(7,2) * t188 - Ifges(7,5) * t201 - Ifges(7,6) * t200 - Ifges(7,3) * t269 - t238 * t212 + t237 * t213;
t316 = -mrSges(6,1) * t192 + mrSges(6,2) * t193 - Ifges(6,5) * t222 - Ifges(6,6) * t221 - Ifges(6,3) * t274 - pkin(5) * t175 - t259 * t232 + t258 * t233 + t319;
t314 = mrSges(5,1) * t208 - mrSges(5,2) * t209 + Ifges(5,5) * t256 + Ifges(5,6) * t255 + Ifges(5,3) * t278 + pkin(4) * t168 + t282 * t250 - t281 * t251 - t316;
t148 = -mrSges(4,1) * t242 + mrSges(4,3) * t236 + Ifges(4,4) * t285 + Ifges(4,2) * t286 + Ifges(4,6) * qJDD(3) - pkin(3) * t161 + qJD(3) * t273 - t271 * t330 - t314;
t135 = mrSges(3,1) * t247 - mrSges(3,2) * t248 + Ifges(3,3) * qJDD(2) + pkin(2) * t317 + pkin(8) * t326 + t305 * t140 + t310 * t148;
t137 = mrSges(3,2) * t265 - mrSges(3,3) * t247 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t313 - pkin(8) * t155 + t140 * t310 - t148 * t305;
t320 = mrSges(2,1) * t288 - mrSges(2,2) * t289 + pkin(1) * t144 + t301 * t135 + t137 * t332 + t337 * t299;
t145 = m(2) * t289 + t147;
t143 = t301 * t154 + (t153 * t306 + t333) * t299;
t141 = m(2) * t288 + t144;
t133 = mrSges(2,2) * t297 - mrSges(2,3) * t288 + t311 * t137 - t306 * t139 + (-t143 * t299 - t144 * t301) * pkin(7);
t132 = -mrSges(2,1) * t297 + mrSges(2,3) * t289 - pkin(1) * t143 - t299 * t135 + (t137 * t306 + t337) * t301;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t300 * t133 - t298 * t132 - qJ(1) * (t141 * t300 + t145 * t298), t133, t137, t140, t150, t164, t177; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t298 * t133 + t300 * t132 + qJ(1) * (-t141 * t298 + t145 * t300), t132, t139, t148, t149, t163, t176; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t320, t320, t135, t336, t314, -t316, -t319;];
m_new  = t1;
