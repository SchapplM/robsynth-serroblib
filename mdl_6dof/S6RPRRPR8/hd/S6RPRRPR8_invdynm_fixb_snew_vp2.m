% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:08:13
% EndTime: 2019-05-05 23:08:43
% DurationCPUTime: 15.20s
% Computational Cost: add. (289673->340), mult. (588269->417), div. (0->0), fcn. (392875->10), ass. (0->137)
t304 = sin(qJ(1));
t308 = cos(qJ(1));
t284 = -t308 * g(1) - t304 * g(2);
t337 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t284;
t336 = (-pkin(1) - pkin(7));
t335 = mrSges(2,1) - mrSges(3,2);
t334 = -Ifges(3,4) + Ifges(2,5);
t333 = (Ifges(3,5) - Ifges(2,6));
t310 = qJD(1) ^ 2;
t253 = (t310 * t336) - t337;
t307 = cos(qJ(3));
t331 = qJD(1) * qJD(3);
t287 = t307 * t331;
t303 = sin(qJ(3));
t278 = -t303 * qJDD(1) - t287;
t329 = t303 * t331;
t279 = qJDD(1) * t307 - t329;
t224 = (-t279 + t329) * pkin(8) + (-t278 + t287) * pkin(3) + t253;
t283 = t304 * g(1) - t308 * g(2);
t323 = -t310 * qJ(2) + qJDD(2) - t283;
t254 = qJDD(1) * t336 + t323;
t247 = -g(3) * t307 + t303 * t254;
t277 = (pkin(3) * t303 - pkin(8) * t307) * qJD(1);
t289 = t303 * qJD(1);
t309 = qJD(3) ^ 2;
t227 = -pkin(3) * t309 + qJDD(3) * pkin(8) - t277 * t289 + t247;
t302 = sin(qJ(4));
t306 = cos(qJ(4));
t205 = t306 * t224 - t302 * t227;
t332 = qJD(1) * t307;
t274 = qJD(3) * t306 - t302 * t332;
t240 = qJD(4) * t274 + qJDD(3) * t302 + t279 * t306;
t273 = qJDD(4) - t278;
t275 = qJD(3) * t302 + t306 * t332;
t286 = t289 + qJD(4);
t195 = (t274 * t286 - t240) * qJ(5) + (t274 * t275 + t273) * pkin(4) + t205;
t206 = t302 * t224 + t306 * t227;
t239 = -qJD(4) * t275 + qJDD(3) * t306 - t279 * t302;
t249 = pkin(4) * t286 - qJ(5) * t275;
t272 = t274 ^ 2;
t197 = -pkin(4) * t272 + qJ(5) * t239 - t249 * t286 + t206;
t299 = sin(pkin(10));
t300 = cos(pkin(10));
t243 = t274 * t299 + t275 * t300;
t182 = -0.2e1 * qJD(5) * t243 + t300 * t195 - t299 * t197;
t217 = t239 * t299 + t240 * t300;
t242 = t274 * t300 - t275 * t299;
t179 = (t242 * t286 - t217) * pkin(9) + (t242 * t243 + t273) * pkin(5) + t182;
t183 = 0.2e1 * qJD(5) * t242 + t299 * t195 + t300 * t197;
t216 = t239 * t300 - t240 * t299;
t230 = pkin(5) * t286 - pkin(9) * t243;
t241 = t242 ^ 2;
t180 = -pkin(5) * t241 + pkin(9) * t216 - t230 * t286 + t183;
t301 = sin(qJ(6));
t305 = cos(qJ(6));
t177 = t179 * t305 - t180 * t301;
t219 = t242 * t305 - t243 * t301;
t191 = qJD(6) * t219 + t216 * t301 + t217 * t305;
t220 = t242 * t301 + t243 * t305;
t203 = -mrSges(7,1) * t219 + mrSges(7,2) * t220;
t285 = qJD(6) + t286;
t208 = -mrSges(7,2) * t285 + mrSges(7,3) * t219;
t266 = qJDD(6) + t273;
t171 = m(7) * t177 + mrSges(7,1) * t266 - mrSges(7,3) * t191 - t203 * t220 + t208 * t285;
t178 = t179 * t301 + t180 * t305;
t190 = -qJD(6) * t220 + t216 * t305 - t217 * t301;
t209 = mrSges(7,1) * t285 - mrSges(7,3) * t220;
t172 = m(7) * t178 - mrSges(7,2) * t266 + mrSges(7,3) * t190 + t203 * t219 - t209 * t285;
t165 = t305 * t171 + t301 * t172;
t221 = -mrSges(6,1) * t242 + mrSges(6,2) * t243;
t228 = -mrSges(6,2) * t286 + mrSges(6,3) * t242;
t162 = m(6) * t182 + mrSges(6,1) * t273 - mrSges(6,3) * t217 - t221 * t243 + t228 * t286 + t165;
t229 = mrSges(6,1) * t286 - mrSges(6,3) * t243;
t326 = -t171 * t301 + t305 * t172;
t163 = m(6) * t183 - mrSges(6,2) * t273 + mrSges(6,3) * t216 + t221 * t242 - t229 * t286 + t326;
t158 = t300 * t162 + t299 * t163;
t245 = -mrSges(5,1) * t274 + mrSges(5,2) * t275;
t248 = -mrSges(5,2) * t286 + mrSges(5,3) * t274;
t156 = m(5) * t205 + mrSges(5,1) * t273 - mrSges(5,3) * t240 - t245 * t275 + t248 * t286 + t158;
t250 = mrSges(5,1) * t286 - mrSges(5,3) * t275;
t327 = -t162 * t299 + t300 * t163;
t157 = m(5) * t206 - mrSges(5,2) * t273 + mrSges(5,3) * t239 + t245 * t274 - t250 * t286 + t327;
t149 = t306 * t156 + t302 * t157;
t276 = (mrSges(4,1) * t303 + mrSges(4,2) * t307) * qJD(1);
t282 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t332;
t328 = -t156 * t302 + t306 * t157;
t147 = m(4) * t247 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t278 - qJD(3) * t282 - t276 * t289 + t328;
t246 = g(3) * t303 + t254 * t307;
t281 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t289;
t226 = -qJDD(3) * pkin(3) - pkin(8) * t309 + t277 * t332 - t246;
t204 = -pkin(4) * t239 - qJ(5) * t272 + t275 * t249 + qJDD(5) + t226;
t185 = -pkin(5) * t216 - pkin(9) * t241 + t230 * t243 + t204;
t325 = m(7) * t185 - t190 * mrSges(7,1) + t191 * mrSges(7,2) - t219 * t208 + t220 * t209;
t316 = m(6) * t204 - t216 * mrSges(6,1) + t217 * mrSges(6,2) - t242 * t228 + t243 * t229 + t325;
t312 = -m(5) * t226 + t239 * mrSges(5,1) - t240 * mrSges(5,2) + t274 * t248 - t275 * t250 - t316;
t173 = m(4) * t246 + qJDD(3) * mrSges(4,1) - t279 * mrSges(4,3) + qJD(3) * t281 - t276 * t332 + t312;
t142 = t307 * t147 - t173 * t303;
t141 = t303 * t147 + t307 * t173;
t259 = -qJDD(1) * pkin(1) + t323;
t322 = -m(3) * t259 + (t310 * mrSges(3,3)) - t141;
t145 = -m(4) * t253 + mrSges(4,1) * t278 - t279 * mrSges(4,2) - t281 * t289 - t282 * t332 - t149;
t199 = Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t285;
t200 = Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t285;
t321 = -mrSges(7,1) * t177 + mrSges(7,2) * t178 - Ifges(7,5) * t191 - Ifges(7,6) * t190 - Ifges(7,3) * t266 - t220 * t199 + t219 * t200;
t198 = Ifges(7,5) * t220 + Ifges(7,6) * t219 + Ifges(7,3) * t285;
t166 = -mrSges(7,1) * t185 + mrSges(7,3) * t178 + Ifges(7,4) * t191 + Ifges(7,2) * t190 + Ifges(7,6) * t266 - t198 * t220 + t200 * t285;
t167 = mrSges(7,2) * t185 - mrSges(7,3) * t177 + Ifges(7,1) * t191 + Ifges(7,4) * t190 + Ifges(7,5) * t266 + t198 * t219 - t199 * t285;
t213 = Ifges(6,5) * t243 + Ifges(6,6) * t242 + Ifges(6,3) * t286;
t215 = Ifges(6,1) * t243 + Ifges(6,4) * t242 + Ifges(6,5) * t286;
t151 = -mrSges(6,1) * t204 + mrSges(6,3) * t183 + Ifges(6,4) * t217 + Ifges(6,2) * t216 + Ifges(6,6) * t273 - pkin(5) * t325 + pkin(9) * t326 + t305 * t166 + t301 * t167 - t243 * t213 + t286 * t215;
t214 = Ifges(6,4) * t243 + Ifges(6,2) * t242 + Ifges(6,6) * t286;
t152 = mrSges(6,2) * t204 - mrSges(6,3) * t182 + Ifges(6,1) * t217 + Ifges(6,4) * t216 + Ifges(6,5) * t273 - pkin(9) * t165 - t166 * t301 + t167 * t305 + t213 * t242 - t214 * t286;
t231 = Ifges(5,5) * t275 + Ifges(5,6) * t274 + Ifges(5,3) * t286;
t233 = Ifges(5,1) * t275 + Ifges(5,4) * t274 + Ifges(5,5) * t286;
t135 = -mrSges(5,1) * t226 + mrSges(5,3) * t206 + Ifges(5,4) * t240 + Ifges(5,2) * t239 + Ifges(5,6) * t273 - pkin(4) * t316 + qJ(5) * t327 + t300 * t151 + t299 * t152 - t275 * t231 + t286 * t233;
t232 = Ifges(5,4) * t275 + Ifges(5,2) * t274 + Ifges(5,6) * t286;
t137 = mrSges(5,2) * t226 - mrSges(5,3) * t205 + Ifges(5,1) * t240 + Ifges(5,4) * t239 + Ifges(5,5) * t273 - qJ(5) * t158 - t151 * t299 + t152 * t300 + t231 * t274 - t232 * t286;
t262 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t307 - Ifges(4,6) * t303) * qJD(1);
t263 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t307 - Ifges(4,2) * t303) * qJD(1);
t132 = mrSges(4,2) * t253 - mrSges(4,3) * t246 + Ifges(4,1) * t279 + Ifges(4,4) * t278 + Ifges(4,5) * qJDD(3) - pkin(8) * t149 - qJD(3) * t263 - t135 * t302 + t137 * t306 - t262 * t289;
t264 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t307 - Ifges(4,4) * t303) * qJD(1);
t313 = -mrSges(6,1) * t182 + mrSges(6,2) * t183 - Ifges(6,5) * t217 - Ifges(6,6) * t216 - Ifges(6,3) * t273 - pkin(5) * t165 - t243 * t214 + t242 * t215 + t321;
t311 = mrSges(5,1) * t205 - mrSges(5,2) * t206 + Ifges(5,5) * t240 + Ifges(5,6) * t239 + Ifges(5,3) * t273 + pkin(4) * t158 + t275 * t232 - t274 * t233 - t313;
t133 = -mrSges(4,1) * t253 + mrSges(4,3) * t247 + Ifges(4,4) * t279 + Ifges(4,2) * t278 + Ifges(4,6) * qJDD(3) - pkin(3) * t149 + qJD(3) * t264 - t262 * t332 - t311;
t257 = t310 * pkin(1) + t337;
t320 = mrSges(3,2) * t259 - mrSges(3,3) * t257 + Ifges(3,1) * qJDD(1) - pkin(7) * t141 + t307 * t132 - t133 * t303;
t319 = -mrSges(3,1) * t257 - pkin(2) * t145 - pkin(7) * t142 - t303 * t132 - t307 * t133;
t318 = mrSges(4,1) * t246 - mrSges(4,2) * t247 + Ifges(4,5) * t279 + Ifges(4,6) * t278 + Ifges(4,3) * qJDD(3) + pkin(3) * t312 + pkin(8) * t328 + t306 * t135 + t302 * t137 + t263 * t332 + t264 * t289;
t317 = -m(3) * t257 + t310 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t145;
t315 = -mrSges(2,2) * t284 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t322) + qJ(2) * t317 + mrSges(2,1) * t283 + Ifges(2,3) * qJDD(1) + t320;
t314 = mrSges(3,1) * t259 + pkin(2) * t141 + t318;
t143 = m(2) * t284 - mrSges(2,1) * t310 - qJDD(1) * mrSges(2,2) + t317;
t140 = -m(3) * g(3) + t142;
t138 = m(2) * t283 - t310 * mrSges(2,2) + qJDD(1) * t335 + t322;
t130 = t334 * qJDD(1) + (t333 * t310) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - qJ(2) * t140 + t314 - mrSges(2,3) * t283;
t129 = mrSges(2,3) * t284 - pkin(1) * t140 + g(3) * t335 - qJDD(1) * t333 + t310 * t334 + t319;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t308 * t130 - t304 * t129 - pkin(6) * (t138 * t308 + t143 * t304), t130, t320, t132, t137, t152, t167; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t304 * t130 + t308 * t129 + pkin(6) * (-t138 * t304 + t143 * t308), t129, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t310 * Ifges(3,5)) - t314, t133, t135, t151, t166; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t315, t315, mrSges(3,2) * g(3) + t310 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t319, t318, t311, -t313, -t321;];
m_new  = t1;
