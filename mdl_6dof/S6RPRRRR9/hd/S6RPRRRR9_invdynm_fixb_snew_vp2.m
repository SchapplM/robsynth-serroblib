% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:26:28
% EndTime: 2019-05-06 04:26:59
% DurationCPUTime: 15.82s
% Computational Cost: add. (306670->341), mult. (608634->415), div. (0->0), fcn. (412755->10), ass. (0->139)
t303 = sin(qJ(1));
t308 = cos(qJ(1));
t284 = -t308 * g(1) - t303 * g(2);
t337 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t284;
t336 = (-pkin(1) - pkin(7));
t335 = mrSges(2,1) - mrSges(3,2);
t334 = -Ifges(3,4) + Ifges(2,5);
t333 = (Ifges(3,5) - Ifges(2,6));
t310 = qJD(1) ^ 2;
t251 = (t310 * t336) - t337;
t307 = cos(qJ(3));
t331 = qJD(1) * qJD(3);
t287 = t307 * t331;
t302 = sin(qJ(3));
t277 = -t302 * qJDD(1) - t287;
t329 = t302 * t331;
t278 = t307 * qJDD(1) - t329;
t224 = (-t278 + t329) * pkin(8) + (-t277 + t287) * pkin(3) + t251;
t283 = t303 * g(1) - t308 * g(2);
t323 = -t310 * qJ(2) + qJDD(2) - t283;
t252 = qJDD(1) * t336 + t323;
t245 = -t307 * g(3) + t302 * t252;
t276 = (pkin(3) * t302 - pkin(8) * t307) * qJD(1);
t289 = t302 * qJD(1);
t309 = qJD(3) ^ 2;
t227 = -t309 * pkin(3) + qJDD(3) * pkin(8) - t276 * t289 + t245;
t301 = sin(qJ(4));
t306 = cos(qJ(4));
t205 = t306 * t224 - t301 * t227;
t332 = qJD(1) * t307;
t273 = t306 * qJD(3) - t301 * t332;
t238 = t273 * qJD(4) + t301 * qJDD(3) + t306 * t278;
t272 = qJDD(4) - t277;
t274 = t301 * qJD(3) + t306 * t332;
t286 = t289 + qJD(4);
t195 = (t273 * t286 - t238) * pkin(9) + (t273 * t274 + t272) * pkin(4) + t205;
t206 = t301 * t224 + t306 * t227;
t237 = -t274 * qJD(4) + t306 * qJDD(3) - t301 * t278;
t250 = t286 * pkin(4) - t274 * pkin(9);
t271 = t273 ^ 2;
t197 = -t271 * pkin(4) + t237 * pkin(9) - t286 * t250 + t206;
t300 = sin(qJ(5));
t305 = cos(qJ(5));
t182 = t305 * t195 - t300 * t197;
t240 = t305 * t273 - t300 * t274;
t212 = t240 * qJD(5) + t300 * t237 + t305 * t238;
t241 = t300 * t273 + t305 * t274;
t265 = qJDD(5) + t272;
t285 = qJD(5) + t286;
t179 = (t240 * t285 - t212) * pkin(10) + (t240 * t241 + t265) * pkin(5) + t182;
t183 = t300 * t195 + t305 * t197;
t211 = -t241 * qJD(5) + t305 * t237 - t300 * t238;
t230 = t285 * pkin(5) - t241 * pkin(10);
t239 = t240 ^ 2;
t180 = -t239 * pkin(5) + t211 * pkin(10) - t285 * t230 + t183;
t299 = sin(qJ(6));
t304 = cos(qJ(6));
t177 = t304 * t179 - t299 * t180;
t219 = t304 * t240 - t299 * t241;
t191 = t219 * qJD(6) + t299 * t211 + t304 * t212;
t220 = t299 * t240 + t304 * t241;
t203 = -t219 * mrSges(7,1) + t220 * mrSges(7,2);
t280 = qJD(6) + t285;
t213 = -t280 * mrSges(7,2) + t219 * mrSges(7,3);
t261 = qJDD(6) + t265;
t174 = m(7) * t177 + t261 * mrSges(7,1) - t191 * mrSges(7,3) - t220 * t203 + t280 * t213;
t178 = t299 * t179 + t304 * t180;
t190 = -t220 * qJD(6) + t304 * t211 - t299 * t212;
t214 = t280 * mrSges(7,1) - t220 * mrSges(7,3);
t175 = m(7) * t178 - t261 * mrSges(7,2) + t190 * mrSges(7,3) + t219 * t203 - t280 * t214;
t165 = t304 * t174 + t299 * t175;
t221 = -t240 * mrSges(6,1) + t241 * mrSges(6,2);
t228 = -t285 * mrSges(6,2) + t240 * mrSges(6,3);
t162 = m(6) * t182 + t265 * mrSges(6,1) - t212 * mrSges(6,3) - t241 * t221 + t285 * t228 + t165;
t229 = t285 * mrSges(6,1) - t241 * mrSges(6,3);
t326 = -t299 * t174 + t304 * t175;
t163 = m(6) * t183 - t265 * mrSges(6,2) + t211 * mrSges(6,3) + t240 * t221 - t285 * t229 + t326;
t158 = t305 * t162 + t300 * t163;
t243 = -t273 * mrSges(5,1) + t274 * mrSges(5,2);
t246 = -t286 * mrSges(5,2) + t273 * mrSges(5,3);
t156 = m(5) * t205 + t272 * mrSges(5,1) - t238 * mrSges(5,3) - t274 * t243 + t286 * t246 + t158;
t247 = t286 * mrSges(5,1) - t274 * mrSges(5,3);
t327 = -t300 * t162 + t305 * t163;
t157 = m(5) * t206 - t272 * mrSges(5,2) + t237 * mrSges(5,3) + t273 * t243 - t286 * t247 + t327;
t149 = t306 * t156 + t301 * t157;
t275 = (mrSges(4,1) * t302 + mrSges(4,2) * t307) * qJD(1);
t282 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t332;
t328 = -t301 * t156 + t306 * t157;
t147 = m(4) * t245 - qJDD(3) * mrSges(4,2) + t277 * mrSges(4,3) - qJD(3) * t282 - t275 * t289 + t328;
t244 = t302 * g(3) + t307 * t252;
t281 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t289;
t226 = -qJDD(3) * pkin(3) - t309 * pkin(8) + t276 * t332 - t244;
t204 = -t237 * pkin(4) - t271 * pkin(9) + t274 * t250 + t226;
t185 = -t211 * pkin(5) - t239 * pkin(10) + t241 * t230 + t204;
t325 = m(7) * t185 - t190 * mrSges(7,1) + t191 * mrSges(7,2) - t219 * t213 + t220 * t214;
t316 = m(6) * t204 - t211 * mrSges(6,1) + t212 * mrSges(6,2) - t240 * t228 + t241 * t229 + t325;
t312 = -m(5) * t226 + t237 * mrSges(5,1) - t238 * mrSges(5,2) + t273 * t246 - t274 * t247 - t316;
t168 = m(4) * t244 + qJDD(3) * mrSges(4,1) - t278 * mrSges(4,3) + qJD(3) * t281 - t275 * t332 + t312;
t142 = t307 * t147 - t302 * t168;
t141 = t302 * t147 + t307 * t168;
t258 = -qJDD(1) * pkin(1) + t323;
t322 = -m(3) * t258 + (t310 * mrSges(3,3)) - t141;
t145 = -m(4) * t251 + t277 * mrSges(4,1) - t278 * mrSges(4,2) - t281 * t289 - t282 * t332 - t149;
t199 = Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t280;
t200 = Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t280;
t321 = -mrSges(7,1) * t177 + mrSges(7,2) * t178 - Ifges(7,5) * t191 - Ifges(7,6) * t190 - Ifges(7,3) * t261 - t220 * t199 + t219 * t200;
t198 = Ifges(7,5) * t220 + Ifges(7,6) * t219 + Ifges(7,3) * t280;
t166 = -mrSges(7,1) * t185 + mrSges(7,3) * t178 + Ifges(7,4) * t191 + Ifges(7,2) * t190 + Ifges(7,6) * t261 - t220 * t198 + t280 * t200;
t167 = mrSges(7,2) * t185 - mrSges(7,3) * t177 + Ifges(7,1) * t191 + Ifges(7,4) * t190 + Ifges(7,5) * t261 + t219 * t198 - t280 * t199;
t215 = Ifges(6,5) * t241 + Ifges(6,6) * t240 + Ifges(6,3) * t285;
t217 = Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t285;
t151 = -mrSges(6,1) * t204 + mrSges(6,3) * t183 + Ifges(6,4) * t212 + Ifges(6,2) * t211 + Ifges(6,6) * t265 - pkin(5) * t325 + pkin(10) * t326 + t304 * t166 + t299 * t167 - t241 * t215 + t285 * t217;
t216 = Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t285;
t152 = mrSges(6,2) * t204 - mrSges(6,3) * t182 + Ifges(6,1) * t212 + Ifges(6,4) * t211 + Ifges(6,5) * t265 - pkin(10) * t165 - t299 * t166 + t304 * t167 + t240 * t215 - t285 * t216;
t231 = Ifges(5,5) * t274 + Ifges(5,6) * t273 + Ifges(5,3) * t286;
t233 = Ifges(5,1) * t274 + Ifges(5,4) * t273 + Ifges(5,5) * t286;
t135 = -mrSges(5,1) * t226 + mrSges(5,3) * t206 + Ifges(5,4) * t238 + Ifges(5,2) * t237 + Ifges(5,6) * t272 - pkin(4) * t316 + pkin(9) * t327 + t305 * t151 + t300 * t152 - t274 * t231 + t286 * t233;
t232 = Ifges(5,4) * t274 + Ifges(5,2) * t273 + Ifges(5,6) * t286;
t137 = mrSges(5,2) * t226 - mrSges(5,3) * t205 + Ifges(5,1) * t238 + Ifges(5,4) * t237 + Ifges(5,5) * t272 - pkin(9) * t158 - t300 * t151 + t305 * t152 + t273 * t231 - t286 * t232;
t262 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t307 - Ifges(4,6) * t302) * qJD(1);
t263 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t307 - Ifges(4,2) * t302) * qJD(1);
t132 = mrSges(4,2) * t251 - mrSges(4,3) * t244 + Ifges(4,1) * t278 + Ifges(4,4) * t277 + Ifges(4,5) * qJDD(3) - pkin(8) * t149 - qJD(3) * t263 - t301 * t135 + t306 * t137 - t262 * t289;
t264 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t307 - Ifges(4,4) * t302) * qJD(1);
t313 = -mrSges(6,1) * t182 + mrSges(6,2) * t183 - Ifges(6,5) * t212 - Ifges(6,6) * t211 - Ifges(6,3) * t265 - pkin(5) * t165 - t241 * t216 + t240 * t217 + t321;
t311 = mrSges(5,1) * t205 - mrSges(5,2) * t206 + Ifges(5,5) * t238 + Ifges(5,6) * t237 + Ifges(5,3) * t272 + pkin(4) * t158 + t274 * t232 - t273 * t233 - t313;
t133 = -mrSges(4,1) * t251 + mrSges(4,3) * t245 + Ifges(4,4) * t278 + Ifges(4,2) * t277 + Ifges(4,6) * qJDD(3) - pkin(3) * t149 + qJD(3) * t264 - t262 * t332 - t311;
t255 = t310 * pkin(1) + t337;
t320 = mrSges(3,2) * t258 - mrSges(3,3) * t255 + Ifges(3,1) * qJDD(1) - pkin(7) * t141 + t307 * t132 - t302 * t133;
t319 = -mrSges(3,1) * t255 - pkin(2) * t145 - pkin(7) * t142 - t302 * t132 - t307 * t133;
t318 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t278 + Ifges(4,6) * t277 + Ifges(4,3) * qJDD(3) + pkin(3) * t312 + pkin(8) * t328 + t306 * t135 + t301 * t137 + t263 * t332 + t264 * t289;
t317 = -m(3) * t255 + t310 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t145;
t315 = -mrSges(2,2) * t284 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t322) + qJ(2) * t317 + mrSges(2,1) * t283 + Ifges(2,3) * qJDD(1) + t320;
t314 = mrSges(3,1) * t258 + pkin(2) * t141 + t318;
t143 = m(2) * t284 - t310 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t317;
t140 = -m(3) * g(3) + t142;
t138 = m(2) * t283 - t310 * mrSges(2,2) + qJDD(1) * t335 + t322;
t130 = t334 * qJDD(1) + (t333 * t310) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - qJ(2) * t140 + t314 - mrSges(2,3) * t283;
t129 = mrSges(2,3) * t284 - pkin(1) * t140 + g(3) * t335 - qJDD(1) * t333 + t310 * t334 + t319;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t308 * t130 - t303 * t129 - pkin(6) * (t308 * t138 + t303 * t143), t130, t320, t132, t137, t152, t167; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t303 * t130 + t308 * t129 + pkin(6) * (-t303 * t138 + t308 * t143), t129, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t310 * Ifges(3,5)) - t314, t133, t135, t151, t166; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t315, t315, mrSges(3,2) * g(3) + t310 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t319, t318, t311, -t313, -t321;];
m_new  = t1;
