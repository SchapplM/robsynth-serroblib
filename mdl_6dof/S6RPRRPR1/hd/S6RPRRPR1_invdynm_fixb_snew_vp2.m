% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 21:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:55:29
% EndTime: 2019-05-05 21:55:55
% DurationCPUTime: 20.75s
% Computational Cost: add. (390185->342), mult. (836156->433), div. (0->0), fcn. (583663->12), ass. (0->139)
t304 = sin(qJ(1));
t308 = cos(qJ(1));
t283 = t304 * g(1) - t308 * g(2);
t274 = qJDD(1) * pkin(1) + t283;
t284 = -t308 * g(1) - t304 * g(2);
t309 = qJD(1) ^ 2;
t276 = -t309 * pkin(1) + t284;
t298 = sin(pkin(10));
t300 = cos(pkin(10));
t258 = t298 * t274 + t300 * t276;
t254 = -t309 * pkin(2) + qJDD(1) * pkin(7) + t258;
t296 = -g(3) + qJDD(2);
t303 = sin(qJ(3));
t307 = cos(qJ(3));
t237 = -t303 * t254 + t307 * t296;
t329 = qJD(1) * qJD(3);
t328 = t307 * t329;
t277 = t303 * qJDD(1) + t328;
t228 = (-t277 + t328) * pkin(8) + (t303 * t307 * t309 + qJDD(3)) * pkin(3) + t237;
t238 = t307 * t254 + t303 * t296;
t278 = t307 * qJDD(1) - t303 * t329;
t331 = qJD(1) * t303;
t282 = qJD(3) * pkin(3) - pkin(8) * t331;
t295 = t307 ^ 2;
t229 = -t295 * t309 * pkin(3) + t278 * pkin(8) - qJD(3) * t282 + t238;
t302 = sin(qJ(4));
t306 = cos(qJ(4));
t201 = t306 * t228 - t302 * t229;
t269 = (-t302 * t303 + t306 * t307) * qJD(1);
t240 = t269 * qJD(4) + t306 * t277 + t302 * t278;
t270 = (t302 * t307 + t303 * t306) * qJD(1);
t291 = qJDD(3) + qJDD(4);
t292 = qJD(3) + qJD(4);
t192 = (t269 * t292 - t240) * qJ(5) + (t269 * t270 + t291) * pkin(4) + t201;
t202 = t302 * t228 + t306 * t229;
t239 = -t270 * qJD(4) - t302 * t277 + t306 * t278;
t260 = t292 * pkin(4) - t270 * qJ(5);
t262 = t269 ^ 2;
t194 = -t262 * pkin(4) + t239 * qJ(5) - t292 * t260 + t202;
t297 = sin(pkin(11));
t299 = cos(pkin(11));
t251 = t299 * t269 - t297 * t270;
t332 = 2 * qJD(5);
t189 = t297 * t192 + t299 * t194 + t251 * t332;
t213 = t299 * t239 - t297 * t240;
t252 = t297 * t269 + t299 * t270;
t226 = -t251 * mrSges(6,1) + t252 * mrSges(6,2);
t242 = t292 * mrSges(6,1) - t252 * mrSges(6,3);
t227 = -t251 * pkin(5) - t252 * pkin(9);
t290 = t292 ^ 2;
t186 = -t290 * pkin(5) + t291 * pkin(9) + t251 * t227 + t189;
t257 = t300 * t274 - t298 * t276;
t320 = -qJDD(1) * pkin(2) - t257;
t230 = -t278 * pkin(3) + t282 * t331 + (-pkin(8) * t295 - pkin(7)) * t309 + t320;
t196 = -t239 * pkin(4) - t262 * qJ(5) + t270 * t260 + qJDD(5) + t230;
t214 = t297 * t239 + t299 * t240;
t190 = (-t251 * t292 - t214) * pkin(9) + (t252 * t292 - t213) * pkin(5) + t196;
t301 = sin(qJ(6));
t305 = cos(qJ(6));
t183 = -t301 * t186 + t305 * t190;
t235 = -t301 * t252 + t305 * t292;
t200 = t235 * qJD(6) + t305 * t214 + t301 * t291;
t212 = qJDD(6) - t213;
t236 = t305 * t252 + t301 * t292;
t215 = -t235 * mrSges(7,1) + t236 * mrSges(7,2);
t245 = qJD(6) - t251;
t216 = -t245 * mrSges(7,2) + t235 * mrSges(7,3);
t179 = m(7) * t183 + t212 * mrSges(7,1) - t200 * mrSges(7,3) - t236 * t215 + t245 * t216;
t184 = t305 * t186 + t301 * t190;
t199 = -t236 * qJD(6) - t301 * t214 + t305 * t291;
t217 = t245 * mrSges(7,1) - t236 * mrSges(7,3);
t180 = m(7) * t184 - t212 * mrSges(7,2) + t199 * mrSges(7,3) + t235 * t215 - t245 * t217;
t323 = -t301 * t179 + t305 * t180;
t166 = m(6) * t189 - t291 * mrSges(6,2) + t213 * mrSges(6,3) + t251 * t226 - t292 * t242 + t323;
t322 = -t299 * t192 + t297 * t194;
t188 = -0.2e1 * qJD(5) * t252 - t322;
t241 = -t292 * mrSges(6,2) + t251 * mrSges(6,3);
t185 = -t291 * pkin(5) - t290 * pkin(9) + (t332 + t227) * t252 + t322;
t317 = -m(7) * t185 + t199 * mrSges(7,1) - t200 * mrSges(7,2) + t235 * t216 - t236 * t217;
t175 = m(6) * t188 + t291 * mrSges(6,1) - t214 * mrSges(6,3) - t252 * t226 + t292 * t241 + t317;
t160 = t297 * t166 + t299 * t175;
t255 = -t269 * mrSges(5,1) + t270 * mrSges(5,2);
t259 = -t292 * mrSges(5,2) + t269 * mrSges(5,3);
t157 = m(5) * t201 + t291 * mrSges(5,1) - t240 * mrSges(5,3) - t270 * t255 + t292 * t259 + t160;
t261 = t292 * mrSges(5,1) - t270 * mrSges(5,3);
t324 = t299 * t166 - t297 * t175;
t158 = m(5) * t202 - t291 * mrSges(5,2) + t239 * mrSges(5,3) + t269 * t255 - t292 * t261 + t324;
t151 = t306 * t157 + t302 * t158;
t267 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t303 + Ifges(4,2) * t307) * qJD(1);
t268 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t303 + Ifges(4,4) * t307) * qJD(1);
t247 = Ifges(5,4) * t270 + Ifges(5,2) * t269 + Ifges(5,6) * t292;
t248 = Ifges(5,1) * t270 + Ifges(5,4) * t269 + Ifges(5,5) * t292;
t203 = Ifges(7,5) * t236 + Ifges(7,6) * t235 + Ifges(7,3) * t245;
t205 = Ifges(7,1) * t236 + Ifges(7,4) * t235 + Ifges(7,5) * t245;
t172 = -mrSges(7,1) * t185 + mrSges(7,3) * t184 + Ifges(7,4) * t200 + Ifges(7,2) * t199 + Ifges(7,6) * t212 - t236 * t203 + t245 * t205;
t204 = Ifges(7,4) * t236 + Ifges(7,2) * t235 + Ifges(7,6) * t245;
t173 = mrSges(7,2) * t185 - mrSges(7,3) * t183 + Ifges(7,1) * t200 + Ifges(7,4) * t199 + Ifges(7,5) * t212 + t235 * t203 - t245 * t204;
t219 = Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t292;
t220 = Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t292;
t315 = -mrSges(6,1) * t188 + mrSges(6,2) * t189 - Ifges(6,5) * t214 - Ifges(6,6) * t213 - Ifges(6,3) * t291 - pkin(5) * t317 - pkin(9) * t323 - t305 * t172 - t301 * t173 - t252 * t219 + t251 * t220;
t312 = -mrSges(5,1) * t201 + mrSges(5,2) * t202 - Ifges(5,5) * t240 - Ifges(5,6) * t239 - Ifges(5,3) * t291 - pkin(4) * t160 - t270 * t247 + t269 * t248 + t315;
t333 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t277 + Ifges(4,6) * t278 + Ifges(4,3) * qJDD(3) + pkin(3) * t151 + (t303 * t267 - t307 * t268) * qJD(1) - t312;
t275 = (-mrSges(4,1) * t307 + mrSges(4,2) * t303) * qJD(1);
t330 = qJD(1) * t307;
t281 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t330;
t149 = m(4) * t237 + qJDD(3) * mrSges(4,1) - t277 * mrSges(4,3) + qJD(3) * t281 - t275 * t331 + t151;
t280 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t331;
t325 = -t302 * t157 + t306 * t158;
t150 = m(4) * t238 - qJDD(3) * mrSges(4,2) + t278 * mrSges(4,3) - qJD(3) * t280 + t275 * t330 + t325;
t326 = -t303 * t149 + t307 * t150;
t142 = m(3) * t258 - t309 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t326;
t253 = -t309 * pkin(7) + t320;
t168 = t305 * t179 + t301 * t180;
t319 = m(6) * t196 - t213 * mrSges(6,1) + t214 * mrSges(6,2) - t251 * t241 + t252 * t242 + t168;
t314 = m(5) * t230 - t239 * mrSges(5,1) + t240 * mrSges(5,2) - t269 * t259 + t270 * t261 + t319;
t311 = -m(4) * t253 + t278 * mrSges(4,1) - t277 * mrSges(4,2) - t280 * t331 + t281 * t330 - t314;
t162 = m(3) * t257 + qJDD(1) * mrSges(3,1) - t309 * mrSges(3,2) + t311;
t138 = t298 * t142 + t300 * t162;
t144 = t307 * t149 + t303 * t150;
t327 = t300 * t142 - t298 * t162;
t218 = Ifges(6,5) * t252 + Ifges(6,6) * t251 + Ifges(6,3) * t292;
t152 = mrSges(6,2) * t196 - mrSges(6,3) * t188 + Ifges(6,1) * t214 + Ifges(6,4) * t213 + Ifges(6,5) * t291 - pkin(9) * t168 - t301 * t172 + t305 * t173 + t251 * t218 - t292 * t219;
t313 = mrSges(7,1) * t183 - mrSges(7,2) * t184 + Ifges(7,5) * t200 + Ifges(7,6) * t199 + Ifges(7,3) * t212 + t236 * t204 - t235 * t205;
t153 = -mrSges(6,1) * t196 + mrSges(6,3) * t189 + Ifges(6,4) * t214 + Ifges(6,2) * t213 + Ifges(6,6) * t291 - pkin(5) * t168 - t252 * t218 + t292 * t220 - t313;
t246 = Ifges(5,5) * t270 + Ifges(5,6) * t269 + Ifges(5,3) * t292;
t139 = -mrSges(5,1) * t230 + mrSges(5,3) * t202 + Ifges(5,4) * t240 + Ifges(5,2) * t239 + Ifges(5,6) * t291 - pkin(4) * t319 + qJ(5) * t324 + t297 * t152 + t299 * t153 - t270 * t246 + t292 * t248;
t145 = mrSges(5,2) * t230 - mrSges(5,3) * t201 + Ifges(5,1) * t240 + Ifges(5,4) * t239 + Ifges(5,5) * t291 - qJ(5) * t160 + t299 * t152 - t297 * t153 + t269 * t246 - t292 * t247;
t266 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t303 + Ifges(4,6) * t307) * qJD(1);
t131 = -mrSges(4,1) * t253 + mrSges(4,3) * t238 + Ifges(4,4) * t277 + Ifges(4,2) * t278 + Ifges(4,6) * qJDD(3) - pkin(3) * t314 + pkin(8) * t325 + qJD(3) * t268 + t306 * t139 + t302 * t145 - t266 * t331;
t133 = mrSges(4,2) * t253 - mrSges(4,3) * t237 + Ifges(4,1) * t277 + Ifges(4,4) * t278 + Ifges(4,5) * qJDD(3) - pkin(8) * t151 - qJD(3) * t267 - t302 * t139 + t306 * t145 + t266 * t330;
t318 = mrSges(3,1) * t257 - mrSges(3,2) * t258 + Ifges(3,3) * qJDD(1) + pkin(2) * t311 + pkin(7) * t326 + t307 * t131 + t303 * t133;
t316 = mrSges(2,1) * t283 - mrSges(2,2) * t284 + Ifges(2,3) * qJDD(1) + pkin(1) * t138 + t318;
t136 = m(2) * t284 - t309 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t327;
t135 = m(2) * t283 + qJDD(1) * mrSges(2,1) - t309 * mrSges(2,2) + t138;
t134 = -mrSges(3,1) * t296 + mrSges(3,3) * t258 + t309 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t144 - t333;
t129 = mrSges(3,2) * t296 - mrSges(3,3) * t257 + Ifges(3,5) * qJDD(1) - t309 * Ifges(3,6) - pkin(7) * t144 - t303 * t131 + t307 * t133;
t128 = -mrSges(2,2) * g(3) - mrSges(2,3) * t283 + Ifges(2,5) * qJDD(1) - t309 * Ifges(2,6) - qJ(2) * t138 + t300 * t129 - t298 * t134;
t127 = Ifges(2,6) * qJDD(1) + t309 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t284 + t298 * t129 + t300 * t134 - pkin(1) * (m(3) * t296 + t144) + qJ(2) * t327;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t308 * t128 - t304 * t127 - pkin(6) * (t308 * t135 + t304 * t136), t128, t129, t133, t145, t152, t173; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t304 * t128 + t308 * t127 + pkin(6) * (-t304 * t135 + t308 * t136), t127, t134, t131, t139, t153, t172; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t316, t316, t318, t333, -t312, -t315, t313;];
m_new  = t1;
