% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR16_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:52
% EndTime: 2019-12-31 20:45:05
% DurationCPUTime: 6.84s
% Computational Cost: add. (103479->331), mult. (234267->415), div. (0->0), fcn. (163664->10), ass. (0->135)
t327 = -2 * qJD(3);
t282 = cos(pkin(5));
t275 = t282 * qJD(1) + qJD(2);
t285 = sin(qJ(2));
t281 = sin(pkin(5));
t312 = qJD(1) * t281;
t305 = t285 * t312;
t326 = (pkin(2) * t275 + t327) * t305;
t286 = sin(qJ(1));
t290 = cos(qJ(1));
t270 = t286 * g(1) - t290 * g(2);
t291 = qJD(1) ^ 2;
t324 = pkin(7) * t281;
t253 = qJDD(1) * pkin(1) + t291 * t324 + t270;
t271 = -t290 * g(1) - t286 * g(2);
t309 = qJDD(1) * t281;
t254 = -t291 * pkin(1) + pkin(7) * t309 + t271;
t289 = cos(qJ(2));
t317 = t282 * t285;
t319 = t281 * t285;
t216 = -g(3) * t319 + t253 * t317 + t289 * t254;
t255 = (-pkin(2) * t289 - qJ(3) * t285) * t312;
t273 = t275 ^ 2;
t274 = t282 * qJDD(1) + qJDD(2);
t311 = qJD(1) * t289;
t306 = t281 * t311;
t193 = t273 * pkin(2) - t274 * qJ(3) - t255 * t306 + t275 * t327 - t216;
t325 = -pkin(2) - pkin(8);
t323 = t282 * g(3);
t322 = mrSges(3,1) - mrSges(4,2);
t321 = Ifges(3,4) + Ifges(4,6);
t320 = t281 ^ 2 * t291;
t318 = t281 * t289;
t316 = t282 * t289;
t258 = pkin(3) * t305 - t275 * pkin(8);
t259 = (qJD(2) * t311 + qJDD(1) * t285) * t281;
t260 = -qJD(2) * t305 + t289 * t309;
t308 = t289 ^ 2 * t320;
t186 = -pkin(3) * t308 - t323 - t259 * qJ(3) + t325 * t260 + (-t253 + (-qJ(3) * t275 * t289 - t258 * t285) * qJD(1)) * t281 + t326;
t313 = g(3) * t318 + t285 * t254;
t302 = -t273 * qJ(3) + t255 * t305 + qJDD(3) + t313;
t188 = t259 * pkin(3) + t325 * t274 + (-pkin(3) * t275 * t312 - pkin(8) * t285 * t320 - t253 * t282) * t289 + t302;
t284 = sin(qJ(4));
t288 = cos(qJ(4));
t182 = t288 * t186 + t284 * t188;
t239 = -t284 * t275 - t288 * t306;
t240 = t288 * t275 - t284 * t306;
t218 = -t239 * pkin(4) - t240 * pkin(9);
t248 = qJDD(4) + t259;
t265 = qJD(4) + t305;
t263 = t265 ^ 2;
t178 = -t263 * pkin(4) + t248 * pkin(9) + t239 * t218 + t182;
t185 = t260 * pkin(3) - pkin(8) * t308 + t275 * t258 - t193;
t213 = -t240 * qJD(4) - t288 * t260 - t284 * t274;
t214 = t239 * qJD(4) - t284 * t260 + t288 * t274;
t179 = (-t239 * t265 - t214) * pkin(9) + (t240 * t265 - t213) * pkin(4) + t185;
t283 = sin(qJ(5));
t287 = cos(qJ(5));
t175 = -t283 * t178 + t287 * t179;
t219 = -t283 * t240 + t287 * t265;
t191 = t219 * qJD(5) + t287 * t214 + t283 * t248;
t220 = t287 * t240 + t283 * t265;
t201 = -t219 * mrSges(6,1) + t220 * mrSges(6,2);
t238 = qJD(5) - t239;
t203 = -t238 * mrSges(6,2) + t219 * mrSges(6,3);
t211 = qJDD(5) - t213;
t171 = m(6) * t175 + t211 * mrSges(6,1) - t191 * mrSges(6,3) - t220 * t201 + t238 * t203;
t176 = t287 * t178 + t283 * t179;
t190 = -t220 * qJD(5) - t283 * t214 + t287 * t248;
t204 = t238 * mrSges(6,1) - t220 * mrSges(6,3);
t172 = m(6) * t176 - t211 * mrSges(6,2) + t190 * mrSges(6,3) + t219 * t201 - t238 * t204;
t161 = t287 * t171 + t283 * t172;
t229 = Ifges(4,1) * t275 + (-Ifges(4,4) * t285 - Ifges(4,5) * t289) * t312;
t315 = Ifges(3,3) * t275 + (Ifges(3,5) * t285 + Ifges(3,6) * t289) * t312 + t229;
t227 = Ifges(4,5) * t275 + (-Ifges(4,6) * t285 - Ifges(4,3) * t289) * t312;
t314 = -Ifges(3,6) * t275 - (Ifges(3,4) * t285 + Ifges(3,2) * t289) * t312 + t227;
t307 = t253 * t316;
t215 = t307 - t313;
t250 = -t275 * mrSges(3,2) + mrSges(3,3) * t306;
t251 = -mrSges(4,1) * t306 - t275 * mrSges(4,3);
t256 = (mrSges(4,2) * t289 - mrSges(4,3) * t285) * t312;
t257 = (-mrSges(3,1) * t289 + mrSges(3,2) * t285) * t312;
t217 = -t239 * mrSges(5,1) + t240 * mrSges(5,2);
t222 = t265 * mrSges(5,1) - t240 * mrSges(5,3);
t304 = -t283 * t171 + t287 * t172;
t158 = m(5) * t182 - t248 * mrSges(5,2) + t213 * mrSges(5,3) + t239 * t217 - t265 * t222 + t304;
t181 = -t284 * t186 + t288 * t188;
t221 = -t265 * mrSges(5,2) + t239 * mrSges(5,3);
t177 = -t248 * pkin(4) - t263 * pkin(9) + t240 * t218 - t181;
t300 = -m(6) * t177 + t190 * mrSges(6,1) - t191 * mrSges(6,2) + t219 * t203 - t220 * t204;
t167 = m(5) * t181 + t248 * mrSges(5,1) - t214 * mrSges(5,3) - t240 * t217 + t265 * t221 + t300;
t152 = t284 * t158 + t288 * t167;
t200 = -t274 * pkin(2) + t302 - t307;
t301 = -m(4) * t200 - t259 * mrSges(4,1) - t152;
t150 = m(3) * t215 - t259 * mrSges(3,3) + (t250 - t251) * t275 + t322 * t274 + (-t256 - t257) * t305 + t301;
t249 = t275 * mrSges(3,1) - mrSges(3,3) * t305;
t159 = -m(5) * t185 + t213 * mrSges(5,1) - t214 * mrSges(5,2) + t239 * t221 - t240 * t222 - t161;
t252 = mrSges(4,1) * t305 + t275 * mrSges(4,2);
t294 = -m(4) * t193 + t274 * mrSges(4,3) + t275 * t252 + t256 * t306 - t159;
t156 = t294 + t257 * t306 + (mrSges(3,3) + mrSges(4,1)) * t260 - t275 * t249 - t274 * mrSges(3,2) + m(3) * t216;
t144 = -t285 * t150 + t289 * t156;
t153 = t288 * t158 - t284 * t167;
t230 = -t281 * t253 - t323;
t194 = -t260 * pkin(2) + (-t275 * t306 - t259) * qJ(3) + t230 + t326;
t303 = m(4) * t194 - t259 * mrSges(4,3) + t251 * t306 + t153;
t149 = m(3) * t230 + t259 * mrSges(3,2) - t322 * t260 + (-t250 * t289 + (t249 - t252) * t285) * t312 + t303;
t141 = -t281 * t149 + t150 * t316 + t156 * t317;
t226 = Ifges(3,5) * t275 + (Ifges(3,1) * t285 + Ifges(3,4) * t289) * t312;
t195 = Ifges(6,5) * t220 + Ifges(6,6) * t219 + Ifges(6,3) * t238;
t197 = Ifges(6,1) * t220 + Ifges(6,4) * t219 + Ifges(6,5) * t238;
t165 = -mrSges(6,1) * t177 + mrSges(6,3) * t176 + Ifges(6,4) * t191 + Ifges(6,2) * t190 + Ifges(6,6) * t211 - t220 * t195 + t238 * t197;
t196 = Ifges(6,4) * t220 + Ifges(6,2) * t219 + Ifges(6,6) * t238;
t166 = mrSges(6,2) * t177 - mrSges(6,3) * t175 + Ifges(6,1) * t191 + Ifges(6,4) * t190 + Ifges(6,5) * t211 + t219 * t195 - t238 * t196;
t205 = Ifges(5,5) * t240 + Ifges(5,6) * t239 + Ifges(5,3) * t265;
t206 = Ifges(5,4) * t240 + Ifges(5,2) * t239 + Ifges(5,6) * t265;
t146 = mrSges(5,2) * t185 - mrSges(5,3) * t181 + Ifges(5,1) * t214 + Ifges(5,4) * t213 + Ifges(5,5) * t248 - pkin(9) * t161 - t283 * t165 + t287 * t166 + t239 * t205 - t265 * t206;
t207 = Ifges(5,1) * t240 + Ifges(5,4) * t239 + Ifges(5,5) * t265;
t293 = mrSges(6,1) * t175 - mrSges(6,2) * t176 + Ifges(6,5) * t191 + Ifges(6,6) * t190 + Ifges(6,3) * t211 + t220 * t196 - t219 * t197;
t147 = -mrSges(5,1) * t185 + mrSges(5,3) * t182 + Ifges(5,4) * t214 + Ifges(5,2) * t213 + Ifges(5,6) * t248 - pkin(4) * t161 - t240 * t205 + t265 * t207 - t293;
t228 = Ifges(4,4) * t275 + (-Ifges(4,2) * t285 - Ifges(4,6) * t289) * t312;
t297 = mrSges(4,2) * t200 - mrSges(4,3) * t193 + Ifges(4,1) * t274 - Ifges(4,4) * t259 - Ifges(4,5) * t260 - pkin(8) * t152 + t288 * t146 - t284 * t147 + t228 * t306;
t133 = t297 + pkin(2) * (-t274 * mrSges(4,2) - t275 * t251 + t301) + Ifges(3,3) * t274 + Ifges(3,6) * t260 + Ifges(3,5) * t259 + mrSges(3,1) * t215 - mrSges(3,2) * t216 + (-t226 * t289 + (-pkin(2) * t256 - t314) * t285) * t312 + qJ(3) * (t260 * mrSges(4,1) + t294);
t151 = t260 * mrSges(4,2) - t252 * t305 + t303;
t295 = -mrSges(4,1) * t193 + mrSges(4,2) * t194 - pkin(3) * t159 - pkin(8) * t153 - t284 * t146 - t288 * t147;
t135 = -mrSges(3,1) * t230 + mrSges(3,3) * t216 - pkin(2) * t151 + (t226 - t228) * t275 + (Ifges(3,6) - Ifges(4,5)) * t274 + (Ifges(3,2) + Ifges(4,3)) * t260 + t321 * t259 - t315 * t305 + t295;
t296 = -mrSges(5,1) * t181 + mrSges(5,2) * t182 - Ifges(5,5) * t214 - Ifges(5,6) * t213 - Ifges(5,3) * t248 - pkin(4) * t300 - pkin(9) * t304 - t287 * t165 - t283 * t166 - t240 * t206 + t239 * t207;
t292 = -mrSges(4,1) * t200 + mrSges(4,3) * t194 - pkin(3) * t152 + t296;
t137 = t315 * t306 - t292 + t314 * t275 + (Ifges(3,5) - Ifges(4,4)) * t274 + t321 * t260 + (Ifges(3,1) + Ifges(4,2)) * t259 + mrSges(3,2) * t230 - mrSges(3,3) * t215 - qJ(3) * t151;
t299 = mrSges(2,1) * t270 - mrSges(2,2) * t271 + Ifges(2,3) * qJDD(1) + pkin(1) * t141 + t282 * t133 + t135 * t318 + t137 * t319 + t144 * t324;
t142 = m(2) * t271 - t291 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t144;
t140 = t282 * t149 + (t150 * t289 + t156 * t285) * t281;
t138 = m(2) * t270 + qJDD(1) * mrSges(2,1) - t291 * mrSges(2,2) + t141;
t131 = -mrSges(2,2) * g(3) - mrSges(2,3) * t270 + Ifges(2,5) * qJDD(1) - t291 * Ifges(2,6) - t285 * t135 + t289 * t137 + (-t140 * t281 - t141 * t282) * pkin(7);
t130 = mrSges(2,1) * g(3) + mrSges(2,3) * t271 + t291 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t140 - t281 * t133 + (pkin(7) * t144 + t135 * t289 + t137 * t285) * t282;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t290 * t131 - t286 * t130 - pkin(6) * (t290 * t138 + t286 * t142), t131, t137, -t227 * t305 + t297, t146, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t286 * t131 + t290 * t130 + pkin(6) * (-t286 * t138 + t290 * t142), t130, t135, Ifges(4,4) * t274 - Ifges(4,2) * t259 - Ifges(4,6) * t260 - t275 * t227 - t229 * t306 + t292, t147, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t299, t299, t133, Ifges(4,5) * t274 - Ifges(4,6) * t259 - Ifges(4,3) * t260 + t275 * t228 + t229 * t305 - t295, -t296, t293;];
m_new = t1;
