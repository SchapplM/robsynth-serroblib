% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:12:56
% EndTime: 2019-05-05 07:13:38
% DurationCPUTime: 31.42s
% Computational Cost: add. (601540->341), mult. (1211433->440), div. (0->0), fcn. (884190->14), ass. (0->144)
t307 = sin(pkin(11));
t310 = cos(pkin(11));
t293 = t307 * g(1) - t310 * g(2);
t294 = -t310 * g(1) - t307 * g(2);
t305 = -g(3) + qJDD(1);
t318 = cos(qJ(2));
t311 = cos(pkin(6));
t315 = sin(qJ(2));
t340 = t311 * t315;
t308 = sin(pkin(6));
t341 = t308 * t315;
t262 = t293 * t340 + t318 * t294 + t305 * t341;
t319 = qJD(2) ^ 2;
t257 = -t319 * pkin(2) + qJDD(2) * pkin(8) + t262;
t275 = -t308 * t293 + t311 * t305;
t314 = sin(qJ(3));
t317 = cos(qJ(3));
t235 = -t314 * t257 + t317 * t275;
t337 = qJD(2) * qJD(3);
t336 = t317 * t337;
t290 = t314 * qJDD(2) + t336;
t226 = (-t290 + t336) * pkin(9) + (t314 * t317 * t319 + qJDD(3)) * pkin(3) + t235;
t236 = t317 * t257 + t314 * t275;
t291 = t317 * qJDD(2) - t314 * t337;
t339 = qJD(2) * t314;
t298 = qJD(3) * pkin(3) - pkin(9) * t339;
t304 = t317 ^ 2;
t227 = -t304 * t319 * pkin(3) + t291 * pkin(9) - qJD(3) * t298 + t236;
t313 = sin(qJ(4));
t345 = cos(qJ(4));
t211 = t313 * t226 + t227 * t345;
t282 = (t313 * t317 + t314 * t345) * qJD(2);
t251 = t282 * qJD(4) + t313 * t290 - t291 * t345;
t338 = qJD(2) * t317;
t281 = t313 * t339 - t338 * t345;
t265 = t281 * mrSges(5,1) + t282 * mrSges(5,2);
t303 = qJD(3) + qJD(4);
t274 = t303 * mrSges(5,1) - t282 * mrSges(5,3);
t302 = qJDD(3) + qJDD(4);
t264 = t281 * pkin(4) - t282 * qJ(5);
t301 = t303 ^ 2;
t205 = -t301 * pkin(4) + t302 * qJ(5) - t281 * t264 + t211;
t261 = -t315 * t294 + (t293 * t311 + t305 * t308) * t318;
t326 = -qJDD(2) * pkin(2) - t261;
t230 = -t291 * pkin(3) + t298 * t339 + (-pkin(9) * t304 - pkin(8)) * t319 + t326;
t252 = -t281 * qJD(4) + t290 * t345 + t313 * t291;
t208 = (t281 * t303 - t252) * qJ(5) + (t282 * t303 + t251) * pkin(4) + t230;
t306 = sin(pkin(12));
t309 = cos(pkin(12));
t270 = t309 * t282 + t306 * t303;
t200 = -0.2e1 * qJD(5) * t270 - t306 * t205 + t309 * t208;
t241 = t309 * t252 + t306 * t302;
t269 = -t306 * t282 + t309 * t303;
t198 = (t269 * t281 - t241) * pkin(10) + (t269 * t270 + t251) * pkin(5) + t200;
t201 = 0.2e1 * qJD(5) * t269 + t309 * t205 + t306 * t208;
t240 = -t306 * t252 + t309 * t302;
t255 = t281 * pkin(5) - t270 * pkin(10);
t268 = t269 ^ 2;
t199 = -t268 * pkin(5) + t240 * pkin(10) - t281 * t255 + t201;
t312 = sin(qJ(6));
t316 = cos(qJ(6));
t196 = t316 * t198 - t312 * t199;
t238 = t316 * t269 - t312 * t270;
t217 = t238 * qJD(6) + t312 * t240 + t316 * t241;
t239 = t312 * t269 + t316 * t270;
t222 = -t238 * mrSges(7,1) + t239 * mrSges(7,2);
t276 = qJD(6) + t281;
t228 = -t276 * mrSges(7,2) + t238 * mrSges(7,3);
t249 = qJDD(6) + t251;
t191 = m(7) * t196 + t249 * mrSges(7,1) - t217 * mrSges(7,3) - t239 * t222 + t276 * t228;
t197 = t312 * t198 + t316 * t199;
t216 = -t239 * qJD(6) + t316 * t240 - t312 * t241;
t229 = t276 * mrSges(7,1) - t239 * mrSges(7,3);
t192 = m(7) * t197 - t249 * mrSges(7,2) + t216 * mrSges(7,3) + t238 * t222 - t276 * t229;
t183 = t316 * t191 + t312 * t192;
t242 = -t269 * mrSges(6,1) + t270 * mrSges(6,2);
t253 = -t281 * mrSges(6,2) + t269 * mrSges(6,3);
t181 = m(6) * t200 + t251 * mrSges(6,1) - t241 * mrSges(6,3) - t270 * t242 + t281 * t253 + t183;
t254 = t281 * mrSges(6,1) - t270 * mrSges(6,3);
t332 = -t312 * t191 + t316 * t192;
t182 = m(6) * t201 - t251 * mrSges(6,2) + t240 * mrSges(6,3) + t269 * t242 - t281 * t254 + t332;
t333 = -t306 * t181 + t309 * t182;
t174 = m(5) * t211 - t302 * mrSges(5,2) - t251 * mrSges(5,3) - t281 * t265 - t303 * t274 + t333;
t210 = t226 * t345 - t313 * t227;
t273 = -t303 * mrSges(5,2) - t281 * mrSges(5,3);
t204 = -t302 * pkin(4) - t301 * qJ(5) + t282 * t264 + qJDD(5) - t210;
t202 = -t240 * pkin(5) - t268 * pkin(10) + t270 * t255 + t204;
t328 = m(7) * t202 - t216 * mrSges(7,1) + t217 * mrSges(7,2) - t238 * t228 + t239 * t229;
t323 = -m(6) * t204 + t240 * mrSges(6,1) - t241 * mrSges(6,2) + t269 * t253 - t270 * t254 - t328;
t187 = m(5) * t210 + t302 * mrSges(5,1) - t252 * mrSges(5,3) - t282 * t265 + t303 * t273 + t323;
t165 = t313 * t174 + t187 * t345;
t289 = (-mrSges(4,1) * t317 + mrSges(4,2) * t314) * qJD(2);
t296 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t338;
t163 = m(4) * t235 + qJDD(3) * mrSges(4,1) - t290 * mrSges(4,3) + qJD(3) * t296 - t289 * t339 + t165;
t295 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t339;
t334 = t174 * t345 - t313 * t187;
t164 = m(4) * t236 - qJDD(3) * mrSges(4,2) + t291 * mrSges(4,3) - qJD(3) * t295 + t289 * t338 + t334;
t158 = t317 * t163 + t314 * t164;
t279 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t314 + Ifges(4,2) * t317) * qJD(2);
t280 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t314 + Ifges(4,4) * t317) * qJD(2);
t218 = Ifges(7,5) * t239 + Ifges(7,6) * t238 + Ifges(7,3) * t276;
t220 = Ifges(7,1) * t239 + Ifges(7,4) * t238 + Ifges(7,5) * t276;
t184 = -mrSges(7,1) * t202 + mrSges(7,3) * t197 + Ifges(7,4) * t217 + Ifges(7,2) * t216 + Ifges(7,6) * t249 - t239 * t218 + t276 * t220;
t219 = Ifges(7,4) * t239 + Ifges(7,2) * t238 + Ifges(7,6) * t276;
t185 = mrSges(7,2) * t202 - mrSges(7,3) * t196 + Ifges(7,1) * t217 + Ifges(7,4) * t216 + Ifges(7,5) * t249 + t238 * t218 - t276 * t219;
t231 = Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t281;
t233 = Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t281;
t167 = -mrSges(6,1) * t204 + mrSges(6,3) * t201 + Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t251 - pkin(5) * t328 + pkin(10) * t332 + t316 * t184 + t312 * t185 - t270 * t231 + t281 * t233;
t232 = Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t281;
t169 = mrSges(6,2) * t204 - mrSges(6,3) * t200 + Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t251 - pkin(10) * t183 - t312 * t184 + t316 * t185 + t269 * t231 - t281 * t232;
t259 = Ifges(5,4) * t282 - Ifges(5,2) * t281 + Ifges(5,6) * t303;
t260 = Ifges(5,1) * t282 - Ifges(5,4) * t281 + Ifges(5,5) * t303;
t324 = -mrSges(5,1) * t210 + mrSges(5,2) * t211 - Ifges(5,5) * t252 + Ifges(5,6) * t251 - Ifges(5,3) * t302 - pkin(4) * t323 - qJ(5) * t333 - t309 * t167 - t306 * t169 - t282 * t259 - t281 * t260;
t346 = mrSges(4,1) * t235 - mrSges(4,2) * t236 + Ifges(4,5) * t290 + Ifges(4,6) * t291 + Ifges(4,3) * qJDD(3) + pkin(3) * t165 + (t314 * t279 - t317 * t280) * qJD(2) - t324;
t144 = -mrSges(3,1) * t275 + mrSges(3,3) * t262 + t319 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t158 - t346;
t335 = -t314 * t163 + t317 * t164;
t156 = m(3) * t262 - t319 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t335;
t256 = -t319 * pkin(8) + t326;
t176 = t309 * t181 + t306 * t182;
t325 = m(5) * t230 + t251 * mrSges(5,1) + t252 * mrSges(5,2) + t281 * t273 + t282 * t274 + t176;
t322 = -m(4) * t256 + t291 * mrSges(4,1) - t290 * mrSges(4,2) - t295 * t339 + t296 * t338 - t325;
t171 = m(3) * t261 + qJDD(2) * mrSges(3,1) - t319 * mrSges(3,2) + t322;
t152 = t318 * t156 - t315 * t171;
t347 = pkin(7) * t152 + t144 * t318;
t342 = t171 * t318;
t157 = m(3) * t275 + t158;
t149 = t156 * t340 - t308 * t157 + t311 * t342;
t258 = Ifges(5,5) * t282 - Ifges(5,6) * t281 + Ifges(5,3) * t303;
t153 = mrSges(5,2) * t230 - mrSges(5,3) * t210 + Ifges(5,1) * t252 - Ifges(5,4) * t251 + Ifges(5,5) * t302 - qJ(5) * t176 - t306 * t167 + t309 * t169 - t281 * t258 - t303 * t259;
t327 = -mrSges(7,1) * t196 + mrSges(7,2) * t197 - Ifges(7,5) * t217 - Ifges(7,6) * t216 - Ifges(7,3) * t249 - t239 * t219 + t238 * t220;
t321 = -mrSges(6,1) * t200 + mrSges(6,2) * t201 - Ifges(6,5) * t241 - Ifges(6,6) * t240 - pkin(5) * t183 - t270 * t232 + t269 * t233 + t327;
t159 = (-Ifges(5,2) - Ifges(6,3)) * t251 + Ifges(5,6) * t302 + t303 * t260 - t282 * t258 + t321 + Ifges(5,4) * t252 - mrSges(5,1) * t230 + mrSges(5,3) * t211 - pkin(4) * t176;
t278 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t314 + Ifges(4,6) * t317) * qJD(2);
t142 = -mrSges(4,1) * t256 + mrSges(4,3) * t236 + Ifges(4,4) * t290 + Ifges(4,2) * t291 + Ifges(4,6) * qJDD(3) - pkin(3) * t325 + pkin(9) * t334 + qJD(3) * t280 + t313 * t153 + t159 * t345 - t278 * t339;
t145 = mrSges(4,2) * t256 - mrSges(4,3) * t235 + Ifges(4,1) * t290 + Ifges(4,4) * t291 + Ifges(4,5) * qJDD(3) - pkin(9) * t165 - qJD(3) * t279 + t153 * t345 - t313 * t159 + t278 * t338;
t139 = mrSges(3,1) * t261 - mrSges(3,2) * t262 + Ifges(3,3) * qJDD(2) + pkin(2) * t322 + pkin(8) * t335 + t317 * t142 + t314 * t145;
t141 = mrSges(3,2) * t275 - mrSges(3,3) * t261 + Ifges(3,5) * qJDD(2) - t319 * Ifges(3,6) - pkin(8) * t158 - t314 * t142 + t317 * t145;
t329 = mrSges(2,1) * t293 - mrSges(2,2) * t294 + pkin(1) * t149 + t311 * t139 + t141 * t341 + t347 * t308;
t150 = m(2) * t294 + t152;
t148 = t311 * t157 + (t156 * t315 + t342) * t308;
t146 = m(2) * t293 + t149;
t137 = mrSges(2,2) * t305 - mrSges(2,3) * t293 + t318 * t141 - t315 * t144 + (-t148 * t308 - t149 * t311) * pkin(7);
t136 = -mrSges(2,1) * t305 + mrSges(2,3) * t294 - pkin(1) * t148 - t308 * t139 + (t141 * t315 + t347) * t311;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t310 * t137 - t307 * t136 - qJ(1) * (t310 * t146 + t307 * t150), t137, t141, t145, t153, t169, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t307 * t137 + t310 * t136 + qJ(1) * (-t307 * t146 + t310 * t150), t136, t144, t142, t159, t167, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t329, t329, t139, t346, -t324, Ifges(6,3) * t251 - t321, -t327;];
m_new  = t1;
