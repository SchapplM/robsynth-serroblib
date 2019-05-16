% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR7
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
% Datum: 2019-05-06 04:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:03:52
% EndTime: 2019-05-06 04:04:16
% DurationCPUTime: 14.63s
% Computational Cost: add. (273550->342), mult. (573589->420), div. (0->0), fcn. (406397->10), ass. (0->138)
t301 = sin(qJ(1));
t306 = cos(qJ(1));
t278 = -t306 * g(1) - t301 * g(2);
t322 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t278;
t335 = -pkin(1) - pkin(7);
t334 = mrSges(2,1) - mrSges(3,2);
t333 = Ifges(2,5) - Ifges(3,4);
t332 = (-Ifges(2,6) + Ifges(3,5));
t277 = t301 * g(1) - t306 * g(2);
t307 = qJD(1) ^ 2;
t321 = -t307 * qJ(2) + qJDD(2) - t277;
t252 = qJDD(1) * t335 + t321;
t300 = sin(qJ(3));
t305 = cos(qJ(3));
t243 = t300 * g(3) + t305 * t252;
t329 = qJD(1) * qJD(3);
t327 = t300 * t329;
t272 = t305 * qJDD(1) - t327;
t219 = (-t272 - t327) * pkin(8) + (-t300 * t305 * t307 + qJDD(3)) * pkin(3) + t243;
t244 = -t305 * g(3) + t300 * t252;
t271 = -t300 * qJDD(1) - t305 * t329;
t330 = qJD(1) * t305;
t276 = qJD(3) * pkin(3) - pkin(8) * t330;
t294 = t300 ^ 2;
t220 = -t294 * t307 * pkin(3) + t271 * pkin(8) - qJD(3) * t276 + t244;
t299 = sin(qJ(4));
t304 = cos(qJ(4));
t204 = t304 * t219 - t299 * t220;
t262 = (-t299 * t305 - t300 * t304) * qJD(1);
t230 = t262 * qJD(4) + t299 * t271 + t304 * t272;
t263 = (-t299 * t300 + t304 * t305) * qJD(1);
t287 = qJDD(3) + qJDD(4);
t288 = qJD(3) + qJD(4);
t185 = (t262 * t288 - t230) * pkin(9) + (t262 * t263 + t287) * pkin(4) + t204;
t205 = t299 * t219 + t304 * t220;
t229 = -t263 * qJD(4) + t304 * t271 - t299 * t272;
t250 = t288 * pkin(4) - t263 * pkin(9);
t258 = t262 ^ 2;
t187 = -t258 * pkin(4) + t229 * pkin(9) - t288 * t250 + t205;
t298 = sin(qJ(5));
t303 = cos(qJ(5));
t183 = t298 * t185 + t303 * t187;
t240 = t298 * t262 + t303 * t263;
t201 = -t240 * qJD(5) + t303 * t229 - t298 * t230;
t239 = t303 * t262 - t298 * t263;
t214 = -t239 * mrSges(6,1) + t240 * mrSges(6,2);
t283 = qJD(5) + t288;
t232 = t283 * mrSges(6,1) - t240 * mrSges(6,3);
t282 = qJDD(5) + t287;
t215 = -t239 * pkin(5) - t240 * pkin(10);
t281 = t283 ^ 2;
t179 = -t281 * pkin(5) + t282 * pkin(10) + t239 * t215 + t183;
t223 = -t271 * pkin(3) + t276 * t330 + (-pkin(8) * t294 + t335) * t307 + t322;
t192 = -t229 * pkin(4) - t258 * pkin(9) + t263 * t250 + t223;
t202 = t239 * qJD(5) + t298 * t229 + t303 * t230;
t180 = (-t239 * t283 - t202) * pkin(10) + (t240 * t283 - t201) * pkin(5) + t192;
t297 = sin(qJ(6));
t302 = cos(qJ(6));
t176 = -t297 * t179 + t302 * t180;
t224 = -t297 * t240 + t302 * t283;
t190 = t224 * qJD(6) + t302 * t202 + t297 * t282;
t199 = qJDD(6) - t201;
t225 = t302 * t240 + t297 * t283;
t207 = -t224 * mrSges(7,1) + t225 * mrSges(7,2);
t236 = qJD(6) - t239;
t208 = -t236 * mrSges(7,2) + t224 * mrSges(7,3);
t172 = m(7) * t176 + t199 * mrSges(7,1) - t190 * mrSges(7,3) - t225 * t207 + t236 * t208;
t177 = t302 * t179 + t297 * t180;
t189 = -t225 * qJD(6) - t297 * t202 + t302 * t282;
t209 = t236 * mrSges(7,1) - t225 * mrSges(7,3);
t173 = m(7) * t177 - t199 * mrSges(7,2) + t189 * mrSges(7,3) + t224 * t207 - t236 * t209;
t324 = -t297 * t172 + t302 * t173;
t159 = m(6) * t183 - t282 * mrSges(6,2) + t201 * mrSges(6,3) + t239 * t214 - t283 * t232 + t324;
t182 = t303 * t185 - t298 * t187;
t231 = -t283 * mrSges(6,2) + t239 * mrSges(6,3);
t178 = -t282 * pkin(5) - t281 * pkin(10) + t240 * t215 - t182;
t319 = -m(7) * t178 + t189 * mrSges(7,1) - t190 * mrSges(7,2) + t224 * t208 - t225 * t209;
t168 = m(6) * t182 + t282 * mrSges(6,1) - t202 * mrSges(6,3) - t240 * t214 + t283 * t231 + t319;
t153 = t298 * t159 + t303 * t168;
t241 = -t262 * mrSges(5,1) + t263 * mrSges(5,2);
t248 = -t288 * mrSges(5,2) + t262 * mrSges(5,3);
t150 = m(5) * t204 + t287 * mrSges(5,1) - t230 * mrSges(5,3) - t263 * t241 + t288 * t248 + t153;
t249 = t288 * mrSges(5,1) - t263 * mrSges(5,3);
t325 = t303 * t159 - t298 * t168;
t151 = m(5) * t205 - t287 * mrSges(5,2) + t229 * mrSges(5,3) + t262 * t241 - t288 * t249 + t325;
t144 = t304 * t150 + t299 * t151;
t161 = t302 * t172 + t297 * t173;
t331 = qJD(1) * t300;
t270 = (mrSges(4,1) * t300 + mrSges(4,2) * t305) * qJD(1);
t274 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t331;
t141 = m(4) * t243 + qJDD(3) * mrSges(4,1) - t272 * mrSges(4,3) + qJD(3) * t274 - t270 * t330 + t144;
t275 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t330;
t326 = -t299 * t150 + t304 * t151;
t142 = m(4) * t244 - qJDD(3) * mrSges(4,2) + t271 * mrSges(4,3) - qJD(3) * t275 - t270 * t331 + t326;
t138 = -t300 * t141 + t305 * t142;
t137 = t305 * t141 + t300 * t142;
t257 = -qJDD(1) * pkin(1) + t321;
t320 = -m(3) * t257 + (t307 * mrSges(3,3)) - t137;
t318 = m(6) * t192 - t201 * mrSges(6,1) + t202 * mrSges(6,2) - t239 * t231 + t240 * t232 + t161;
t193 = Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t236;
t195 = Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t236;
t165 = -mrSges(7,1) * t178 + mrSges(7,3) * t177 + Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t199 - t225 * t193 + t236 * t195;
t194 = Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t236;
t166 = mrSges(7,2) * t178 - mrSges(7,3) * t176 + Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t199 + t224 * t193 - t236 * t194;
t210 = Ifges(6,5) * t240 + Ifges(6,6) * t239 + Ifges(6,3) * t283;
t211 = Ifges(6,4) * t240 + Ifges(6,2) * t239 + Ifges(6,6) * t283;
t145 = mrSges(6,2) * t192 - mrSges(6,3) * t182 + Ifges(6,1) * t202 + Ifges(6,4) * t201 + Ifges(6,5) * t282 - pkin(10) * t161 - t297 * t165 + t302 * t166 + t239 * t210 - t283 * t211;
t212 = Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t283;
t313 = mrSges(7,1) * t176 - mrSges(7,2) * t177 + Ifges(7,5) * t190 + Ifges(7,6) * t189 + Ifges(7,3) * t199 + t225 * t194 - t224 * t195;
t146 = -mrSges(6,1) * t192 + mrSges(6,3) * t183 + Ifges(6,4) * t202 + Ifges(6,2) * t201 + Ifges(6,6) * t282 - pkin(5) * t161 - t240 * t210 + t283 * t212 - t313;
t233 = Ifges(5,5) * t263 + Ifges(5,6) * t262 + Ifges(5,3) * t288;
t235 = Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * t288;
t133 = -mrSges(5,1) * t223 + mrSges(5,3) * t205 + Ifges(5,4) * t230 + Ifges(5,2) * t229 + Ifges(5,6) * t287 - pkin(4) * t318 + pkin(9) * t325 + t298 * t145 + t303 * t146 - t263 * t233 + t288 * t235;
t234 = Ifges(5,4) * t263 + Ifges(5,2) * t262 + Ifges(5,6) * t288;
t139 = mrSges(5,2) * t223 - mrSges(5,3) * t204 + Ifges(5,1) * t230 + Ifges(5,4) * t229 + Ifges(5,5) * t287 - pkin(9) * t153 + t303 * t145 - t298 * t146 + t262 * t233 - t288 * t234;
t251 = t307 * t335 + t322;
t259 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t305 - Ifges(4,6) * t300) * qJD(1);
t261 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t305 - Ifges(4,4) * t300) * qJD(1);
t312 = m(5) * t223 - t229 * mrSges(5,1) + t230 * mrSges(5,2) - t262 * t248 + t263 * t249 + t318;
t130 = -mrSges(4,1) * t251 + mrSges(4,3) * t244 + Ifges(4,4) * t272 + Ifges(4,2) * t271 + Ifges(4,6) * qJDD(3) - pkin(3) * t312 + pkin(8) * t326 + qJD(3) * t261 + t304 * t133 + t299 * t139 - t259 * t330;
t260 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t305 - Ifges(4,2) * t300) * qJD(1);
t132 = mrSges(4,2) * t251 - mrSges(4,3) * t243 + Ifges(4,1) * t272 + Ifges(4,4) * t271 + Ifges(4,5) * qJDD(3) - pkin(8) * t144 - qJD(3) * t260 - t299 * t133 + t304 * t139 - t259 * t331;
t255 = t307 * pkin(1) - t322;
t317 = mrSges(3,2) * t257 - mrSges(3,3) * t255 + Ifges(3,1) * qJDD(1) - pkin(7) * t137 - t300 * t130 + t305 * t132;
t156 = -m(4) * t251 + t271 * mrSges(4,1) - t272 * mrSges(4,2) - t274 * t331 - t275 * t330 - t312;
t316 = -mrSges(3,1) * t255 - pkin(2) * t156 - pkin(7) * t138 - t305 * t130 - t300 * t132;
t315 = -mrSges(6,1) * t182 + mrSges(6,2) * t183 - Ifges(6,5) * t202 - Ifges(6,6) * t201 - Ifges(6,3) * t282 - pkin(5) * t319 - pkin(10) * t324 - t302 * t165 - t297 * t166 - t240 * t211 + t239 * t212;
t310 = -m(3) * t255 + t307 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t156;
t314 = -mrSges(2,2) * t278 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t320) + qJ(2) * t310 + mrSges(2,1) * t277 + Ifges(2,3) * qJDD(1) + t317;
t311 = -mrSges(5,1) * t204 + mrSges(5,2) * t205 - Ifges(5,5) * t230 - Ifges(5,6) * t229 - Ifges(5,3) * t287 - pkin(4) * t153 - t263 * t234 + t262 * t235 + t315;
t309 = mrSges(4,1) * t243 - mrSges(4,2) * t244 + Ifges(4,5) * t272 + Ifges(4,6) * t271 + Ifges(4,3) * qJDD(3) + pkin(3) * t144 + t260 * t330 + t261 * t331 - t311;
t308 = -mrSges(3,1) * t257 - pkin(2) * t137 - t309;
t154 = m(2) * t278 - t307 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t310;
t136 = -m(3) * g(3) + t138;
t134 = m(2) * t277 - t307 * mrSges(2,2) + qJDD(1) * t334 + t320;
t129 = -t308 + t333 * qJDD(1) + (t332 * t307) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t277 - qJ(2) * t136;
t128 = mrSges(2,3) * t278 - pkin(1) * t136 + g(3) * t334 - qJDD(1) * t332 + t307 * t333 + t316;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t306 * t129 - t301 * t128 - pkin(6) * (t306 * t134 + t301 * t154), t129, t317, t132, t139, t145, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t129 + t306 * t128 + pkin(6) * (-t301 * t134 + t306 * t154), t128, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t307 * Ifges(3,5)) + t308, t130, t133, t146, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t314, t314, mrSges(3,2) * g(3) + t307 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t316, t309, -t311, -t315, t313;];
m_new  = t1;
