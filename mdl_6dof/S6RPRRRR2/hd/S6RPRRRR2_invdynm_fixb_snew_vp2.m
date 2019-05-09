% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 02:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:43:56
% EndTime: 2019-05-06 02:44:28
% DurationCPUTime: 20.46s
% Computational Cost: add. (408289->342), mult. (797586->430), div. (0->0), fcn. (551709->12), ass. (0->139)
t307 = sin(qJ(1));
t312 = cos(qJ(1));
t287 = t307 * g(1) - t312 * g(2);
t278 = qJDD(1) * pkin(1) + t287;
t288 = -t312 * g(1) - t307 * g(2);
t313 = qJD(1) ^ 2;
t280 = -t313 * pkin(1) + t288;
t301 = sin(pkin(11));
t302 = cos(pkin(11));
t257 = t301 * t278 + t302 * t280;
t252 = -t313 * pkin(2) + qJDD(1) * pkin(7) + t257;
t300 = -g(3) + qJDD(2);
t306 = sin(qJ(3));
t311 = cos(qJ(3));
t239 = -t306 * t252 + t311 * t300;
t332 = qJD(1) * qJD(3);
t331 = t311 * t332;
t281 = t306 * qJDD(1) + t331;
t221 = (-t281 + t331) * pkin(8) + (t306 * t311 * t313 + qJDD(3)) * pkin(3) + t239;
t240 = t311 * t252 + t306 * t300;
t282 = t311 * qJDD(1) - t306 * t332;
t334 = qJD(1) * t306;
t286 = qJD(3) * pkin(3) - pkin(8) * t334;
t299 = t311 ^ 2;
t222 = -t299 * t313 * pkin(3) + t282 * pkin(8) - qJD(3) * t286 + t240;
t305 = sin(qJ(4));
t310 = cos(qJ(4));
t207 = t305 * t221 + t310 * t222;
t273 = (t305 * t311 + t306 * t310) * qJD(1);
t241 = -t273 * qJD(4) - t305 * t281 + t310 * t282;
t333 = qJD(1) * t311;
t272 = -t305 * t334 + t310 * t333;
t253 = -t272 * mrSges(5,1) + t273 * mrSges(5,2);
t296 = qJD(3) + qJD(4);
t262 = t296 * mrSges(5,1) - t273 * mrSges(5,3);
t295 = qJDD(3) + qJDD(4);
t254 = -t272 * pkin(4) - t273 * pkin(9);
t294 = t296 ^ 2;
t200 = -t294 * pkin(4) + t295 * pkin(9) + t272 * t254 + t207;
t256 = t302 * t278 - t301 * t280;
t324 = -qJDD(1) * pkin(2) - t256;
t227 = -t282 * pkin(3) + t286 * t334 + (-pkin(8) * t299 - pkin(7)) * t313 + t324;
t242 = t272 * qJD(4) + t310 * t281 + t305 * t282;
t203 = (-t272 * t296 - t242) * pkin(9) + (t273 * t296 - t241) * pkin(4) + t227;
t304 = sin(qJ(5));
t309 = cos(qJ(5));
t190 = -t304 * t200 + t309 * t203;
t259 = -t304 * t273 + t309 * t296;
t215 = t259 * qJD(5) + t309 * t242 + t304 * t295;
t238 = qJDD(5) - t241;
t260 = t309 * t273 + t304 * t296;
t265 = qJD(5) - t272;
t188 = (t259 * t265 - t215) * pkin(10) + (t259 * t260 + t238) * pkin(5) + t190;
t191 = t309 * t200 + t304 * t203;
t214 = -t260 * qJD(5) - t304 * t242 + t309 * t295;
t245 = t265 * pkin(5) - t260 * pkin(10);
t258 = t259 ^ 2;
t189 = -t258 * pkin(5) + t214 * pkin(10) - t265 * t245 + t191;
t303 = sin(qJ(6));
t308 = cos(qJ(6));
t186 = t308 * t188 - t303 * t189;
t228 = t308 * t259 - t303 * t260;
t197 = t228 * qJD(6) + t303 * t214 + t308 * t215;
t229 = t303 * t259 + t308 * t260;
t212 = -t228 * mrSges(7,1) + t229 * mrSges(7,2);
t263 = qJD(6) + t265;
t219 = -t263 * mrSges(7,2) + t228 * mrSges(7,3);
t233 = qJDD(6) + t238;
t181 = m(7) * t186 + t233 * mrSges(7,1) - t197 * mrSges(7,3) - t229 * t212 + t263 * t219;
t187 = t303 * t188 + t308 * t189;
t196 = -t229 * qJD(6) + t308 * t214 - t303 * t215;
t220 = t263 * mrSges(7,1) - t229 * mrSges(7,3);
t182 = m(7) * t187 - t233 * mrSges(7,2) + t196 * mrSges(7,3) + t228 * t212 - t263 * t220;
t173 = t308 * t181 + t303 * t182;
t230 = -t259 * mrSges(6,1) + t260 * mrSges(6,2);
t243 = -t265 * mrSges(6,2) + t259 * mrSges(6,3);
t171 = m(6) * t190 + t238 * mrSges(6,1) - t215 * mrSges(6,3) - t260 * t230 + t265 * t243 + t173;
t244 = t265 * mrSges(6,1) - t260 * mrSges(6,3);
t326 = -t303 * t181 + t308 * t182;
t172 = m(6) * t191 - t238 * mrSges(6,2) + t214 * mrSges(6,3) + t259 * t230 - t265 * t244 + t326;
t327 = -t304 * t171 + t309 * t172;
t164 = m(5) * t207 - t295 * mrSges(5,2) + t241 * mrSges(5,3) + t272 * t253 - t296 * t262 + t327;
t206 = t310 * t221 - t305 * t222;
t261 = -t296 * mrSges(5,2) + t272 * mrSges(5,3);
t199 = -t295 * pkin(4) - t294 * pkin(9) + t273 * t254 - t206;
t192 = -t214 * pkin(5) - t258 * pkin(10) + t260 * t245 + t199;
t322 = m(7) * t192 - t196 * mrSges(7,1) + t197 * mrSges(7,2) - t228 * t219 + t229 * t220;
t317 = -m(6) * t199 + t214 * mrSges(6,1) - t215 * mrSges(6,2) + t259 * t243 - t260 * t244 - t322;
t177 = m(5) * t206 + t295 * mrSges(5,1) - t242 * mrSges(5,3) - t273 * t253 + t296 * t261 + t317;
t154 = t305 * t164 + t310 * t177;
t270 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t306 + Ifges(4,2) * t311) * qJD(1);
t271 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t306 + Ifges(4,4) * t311) * qJD(1);
t208 = Ifges(7,5) * t229 + Ifges(7,6) * t228 + Ifges(7,3) * t263;
t210 = Ifges(7,1) * t229 + Ifges(7,4) * t228 + Ifges(7,5) * t263;
t174 = -mrSges(7,1) * t192 + mrSges(7,3) * t187 + Ifges(7,4) * t197 + Ifges(7,2) * t196 + Ifges(7,6) * t233 - t229 * t208 + t263 * t210;
t209 = Ifges(7,4) * t229 + Ifges(7,2) * t228 + Ifges(7,6) * t263;
t175 = mrSges(7,2) * t192 - mrSges(7,3) * t186 + Ifges(7,1) * t197 + Ifges(7,4) * t196 + Ifges(7,5) * t233 + t228 * t208 - t263 * t209;
t223 = Ifges(6,5) * t260 + Ifges(6,6) * t259 + Ifges(6,3) * t265;
t225 = Ifges(6,1) * t260 + Ifges(6,4) * t259 + Ifges(6,5) * t265;
t156 = -mrSges(6,1) * t199 + mrSges(6,3) * t191 + Ifges(6,4) * t215 + Ifges(6,2) * t214 + Ifges(6,6) * t238 - pkin(5) * t322 + pkin(10) * t326 + t308 * t174 + t303 * t175 - t260 * t223 + t265 * t225;
t224 = Ifges(6,4) * t260 + Ifges(6,2) * t259 + Ifges(6,6) * t265;
t158 = mrSges(6,2) * t199 - mrSges(6,3) * t190 + Ifges(6,1) * t215 + Ifges(6,4) * t214 + Ifges(6,5) * t238 - pkin(10) * t173 - t303 * t174 + t308 * t175 + t259 * t223 - t265 * t224;
t248 = Ifges(5,4) * t273 + Ifges(5,2) * t272 + Ifges(5,6) * t296;
t249 = Ifges(5,1) * t273 + Ifges(5,4) * t272 + Ifges(5,5) * t296;
t318 = -mrSges(5,1) * t206 + mrSges(5,2) * t207 - Ifges(5,5) * t242 - Ifges(5,6) * t241 - Ifges(5,3) * t295 - pkin(4) * t317 - pkin(9) * t327 - t309 * t156 - t304 * t158 - t273 * t248 + t272 * t249;
t335 = mrSges(4,1) * t239 - mrSges(4,2) * t240 + Ifges(4,5) * t281 + Ifges(4,6) * t282 + Ifges(4,3) * qJDD(3) + pkin(3) * t154 + (t306 * t270 - t311 * t271) * qJD(1) - t318;
t279 = (-mrSges(4,1) * t311 + mrSges(4,2) * t306) * qJD(1);
t285 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t333;
t152 = m(4) * t239 + qJDD(3) * mrSges(4,1) - t281 * mrSges(4,3) + qJD(3) * t285 - t279 * t334 + t154;
t284 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t334;
t328 = t310 * t164 - t305 * t177;
t153 = m(4) * t240 - qJDD(3) * mrSges(4,2) + t282 * mrSges(4,3) - qJD(3) * t284 + t279 * t333 + t328;
t329 = -t306 * t152 + t311 * t153;
t145 = m(3) * t257 - t313 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t329;
t251 = -t313 * pkin(7) + t324;
t166 = t309 * t171 + t304 * t172;
t320 = m(5) * t227 - t241 * mrSges(5,1) + t242 * mrSges(5,2) - t272 * t261 + t273 * t262 + t166;
t316 = -m(4) * t251 + t282 * mrSges(4,1) - t281 * mrSges(4,2) - t284 * t334 + t285 * t333 - t320;
t160 = m(3) * t256 + qJDD(1) * mrSges(3,1) - t313 * mrSges(3,2) + t316;
t141 = t301 * t145 + t302 * t160;
t147 = t311 * t152 + t306 * t153;
t330 = t302 * t145 - t301 * t160;
t247 = Ifges(5,5) * t273 + Ifges(5,6) * t272 + Ifges(5,3) * t296;
t142 = mrSges(5,2) * t227 - mrSges(5,3) * t206 + Ifges(5,1) * t242 + Ifges(5,4) * t241 + Ifges(5,5) * t295 - pkin(9) * t166 - t304 * t156 + t309 * t158 + t272 * t247 - t296 * t248;
t321 = -mrSges(7,1) * t186 + mrSges(7,2) * t187 - Ifges(7,5) * t197 - Ifges(7,6) * t196 - Ifges(7,3) * t233 - t229 * t209 + t228 * t210;
t314 = mrSges(6,1) * t190 - mrSges(6,2) * t191 + Ifges(6,5) * t215 + Ifges(6,6) * t214 + Ifges(6,3) * t238 + pkin(5) * t173 + t260 * t224 - t259 * t225 - t321;
t148 = -mrSges(5,1) * t227 + mrSges(5,3) * t207 + Ifges(5,4) * t242 + Ifges(5,2) * t241 + Ifges(5,6) * t295 - pkin(4) * t166 - t273 * t247 + t296 * t249 - t314;
t269 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t306 + Ifges(4,6) * t311) * qJD(1);
t134 = -mrSges(4,1) * t251 + mrSges(4,3) * t240 + Ifges(4,4) * t281 + Ifges(4,2) * t282 + Ifges(4,6) * qJDD(3) - pkin(3) * t320 + pkin(8) * t328 + qJD(3) * t271 + t305 * t142 + t310 * t148 - t269 * t334;
t137 = mrSges(4,2) * t251 - mrSges(4,3) * t239 + Ifges(4,1) * t281 + Ifges(4,4) * t282 + Ifges(4,5) * qJDD(3) - pkin(8) * t154 - qJD(3) * t270 + t310 * t142 - t305 * t148 + t269 * t333;
t323 = mrSges(3,1) * t256 - mrSges(3,2) * t257 + Ifges(3,3) * qJDD(1) + pkin(2) * t316 + pkin(7) * t329 + t311 * t134 + t306 * t137;
t319 = mrSges(2,1) * t287 - mrSges(2,2) * t288 + Ifges(2,3) * qJDD(1) + pkin(1) * t141 + t323;
t139 = m(2) * t288 - t313 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t330;
t138 = m(2) * t287 + qJDD(1) * mrSges(2,1) - t313 * mrSges(2,2) + t141;
t135 = -mrSges(3,1) * t300 + mrSges(3,3) * t257 + t313 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t147 - t335;
t132 = mrSges(3,2) * t300 - mrSges(3,3) * t256 + Ifges(3,5) * qJDD(1) - t313 * Ifges(3,6) - pkin(7) * t147 - t306 * t134 + t311 * t137;
t131 = -mrSges(2,2) * g(3) - mrSges(2,3) * t287 + Ifges(2,5) * qJDD(1) - t313 * Ifges(2,6) - qJ(2) * t141 + t302 * t132 - t301 * t135;
t130 = Ifges(2,6) * qJDD(1) + t313 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t288 + t301 * t132 + t302 * t135 - pkin(1) * (m(3) * t300 + t147) + qJ(2) * t330;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t312 * t131 - t307 * t130 - pkin(6) * (t312 * t138 + t307 * t139), t131, t132, t137, t142, t158, t175; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t307 * t131 + t312 * t130 + pkin(6) * (-t307 * t138 + t312 * t139), t130, t135, t134, t148, t156, t174; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t319, t319, t323, t335, -t318, t314, -t321;];
m_new  = t1;
