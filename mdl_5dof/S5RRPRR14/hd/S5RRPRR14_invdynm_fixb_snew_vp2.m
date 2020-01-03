% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:36:05
% EndTime: 2019-12-31 20:36:35
% DurationCPUTime: 17.94s
% Computational Cost: add. (302897->324), mult. (692040->427), div. (0->0), fcn. (541861->12), ass. (0->136)
t271 = sin(pkin(5));
t304 = pkin(7) * t271;
t273 = cos(pkin(5));
t303 = t273 * g(3);
t276 = sin(qJ(2));
t302 = t271 * t276;
t280 = cos(qJ(2));
t301 = t271 * t280;
t300 = t273 * t276;
t299 = t273 * t280;
t297 = qJD(1) * t271;
t254 = (-pkin(2) * t280 - qJ(3) * t276) * t297;
t266 = t273 * qJD(1) + qJD(2);
t264 = t266 ^ 2;
t265 = t273 * qJDD(1) + qJDD(2);
t296 = qJD(1) * t280;
t277 = sin(qJ(1));
t281 = cos(qJ(1));
t261 = t277 * g(1) - t281 * g(2);
t282 = qJD(1) ^ 2;
t252 = qJDD(1) * pkin(1) + t282 * t304 + t261;
t262 = -t281 * g(1) - t277 * g(2);
t295 = qJDD(1) * t271;
t253 = -t282 * pkin(1) + pkin(7) * t295 + t262;
t298 = t252 * t300 + t280 * t253;
t212 = -t264 * pkin(2) + t265 * qJ(3) + (-g(3) * t276 + t254 * t296) * t271 + t298;
t256 = (qJD(2) * t296 + qJDD(1) * t276) * t271;
t294 = t276 * t297;
t257 = -qJD(2) * t294 + t280 * t295;
t213 = -t257 * pkin(2) - t303 - t256 * qJ(3) + (-t252 + (pkin(2) * t276 - qJ(3) * t280) * t266 * qJD(1)) * t271;
t270 = sin(pkin(10));
t272 = cos(pkin(10));
t245 = t270 * t266 + t272 * t294;
t182 = -0.2e1 * qJD(3) * t245 - t270 * t212 + t272 * t213;
t233 = t272 * t256 + t270 * t265;
t244 = t272 * t266 - t270 * t294;
t293 = t271 * t296;
t179 = (-t244 * t293 - t233) * pkin(8) + (t244 * t245 - t257) * pkin(3) + t182;
t183 = 0.2e1 * qJD(3) * t244 + t272 * t212 + t270 * t213;
t232 = -t270 * t256 + t272 * t265;
t234 = -pkin(3) * t293 - t245 * pkin(8);
t243 = t244 ^ 2;
t181 = -t243 * pkin(3) + t232 * pkin(8) + t234 * t293 + t183;
t275 = sin(qJ(4));
t279 = cos(qJ(4));
t176 = t275 * t179 + t279 * t181;
t227 = t275 * t244 + t279 * t245;
t198 = -t227 * qJD(4) + t279 * t232 - t275 * t233;
t226 = t279 * t244 - t275 * t245;
t206 = -t226 * mrSges(5,1) + t227 * mrSges(5,2);
t260 = qJD(4) - t293;
t217 = t260 * mrSges(5,1) - t227 * mrSges(5,3);
t249 = qJDD(4) - t257;
t207 = -t226 * pkin(4) - t227 * pkin(9);
t259 = t260 ^ 2;
t173 = -t259 * pkin(4) + t249 * pkin(9) + t226 * t207 + t176;
t223 = -g(3) * t301 + t252 * t299 - t276 * t253;
t211 = -t265 * pkin(2) - t264 * qJ(3) + t254 * t294 + qJDD(3) - t223;
t187 = -t232 * pkin(3) - t243 * pkin(8) + t245 * t234 + t211;
t199 = t226 * qJD(4) + t275 * t232 + t279 * t233;
t177 = (-t226 * t260 - t199) * pkin(9) + (t227 * t260 - t198) * pkin(4) + t187;
t274 = sin(qJ(5));
t278 = cos(qJ(5));
t170 = -t274 * t173 + t278 * t177;
t214 = -t274 * t227 + t278 * t260;
t186 = t214 * qJD(5) + t278 * t199 + t274 * t249;
t215 = t278 * t227 + t274 * t260;
t193 = -t214 * mrSges(6,1) + t215 * mrSges(6,2);
t197 = qJDD(5) - t198;
t225 = qJD(5) - t226;
t200 = -t225 * mrSges(6,2) + t214 * mrSges(6,3);
t166 = m(6) * t170 + t197 * mrSges(6,1) - t186 * mrSges(6,3) - t215 * t193 + t225 * t200;
t171 = t278 * t173 + t274 * t177;
t185 = -t215 * qJD(5) - t274 * t199 + t278 * t249;
t201 = t225 * mrSges(6,1) - t215 * mrSges(6,3);
t167 = m(6) * t171 - t197 * mrSges(6,2) + t185 * mrSges(6,3) + t214 * t193 - t225 * t201;
t290 = -t274 * t166 + t278 * t167;
t153 = m(5) * t176 - t249 * mrSges(5,2) + t198 * mrSges(5,3) + t226 * t206 - t260 * t217 + t290;
t175 = t279 * t179 - t275 * t181;
t216 = -t260 * mrSges(5,2) + t226 * mrSges(5,3);
t172 = -t249 * pkin(4) - t259 * pkin(9) + t227 * t207 - t175;
t289 = -m(6) * t172 + t185 * mrSges(6,1) - t186 * mrSges(6,2) + t214 * t200 - t215 * t201;
t162 = m(5) * t175 + t249 * mrSges(5,1) - t199 * mrSges(5,3) - t227 * t206 + t260 * t216 + t289;
t148 = t275 * t153 + t279 * t162;
t228 = -t244 * mrSges(4,1) + t245 * mrSges(4,2);
t230 = mrSges(4,2) * t293 + t244 * mrSges(4,3);
t146 = m(4) * t182 - t257 * mrSges(4,1) - t233 * mrSges(4,3) - t245 * t228 - t230 * t293 + t148;
t231 = -mrSges(4,1) * t293 - t245 * mrSges(4,3);
t291 = t279 * t153 - t275 * t162;
t147 = m(4) * t183 + t257 * mrSges(4,2) + t232 * mrSges(4,3) + t244 * t228 + t231 * t293 + t291;
t140 = t272 * t146 + t270 * t147;
t155 = t278 * t166 + t274 * t167;
t224 = -g(3) * t302 + t298;
t250 = t266 * mrSges(3,1) - mrSges(3,3) * t294;
t255 = (-mrSges(3,1) * t280 + mrSges(3,2) * t276) * t297;
t292 = -t270 * t146 + t272 * t147;
t138 = m(3) * t224 - t265 * mrSges(3,2) + t257 * mrSges(3,3) - t266 * t250 + t255 * t293 + t292;
t251 = -t266 * mrSges(3,2) + mrSges(3,3) * t293;
t287 = m(5) * t187 - t198 * mrSges(5,1) + t199 * mrSges(5,2) - t226 * t216 + t227 * t217 + t155;
t284 = -m(4) * t211 + t232 * mrSges(4,1) - t233 * mrSges(4,2) + t244 * t230 - t245 * t231 - t287;
t150 = m(3) * t223 + t265 * mrSges(3,1) - t256 * mrSges(3,3) + t266 * t251 - t255 * t294 + t284;
t135 = t280 * t138 - t276 * t150;
t238 = -t271 * t252 - t303;
t139 = m(3) * t238 - t257 * mrSges(3,1) + t256 * mrSges(3,2) + (t250 * t276 - t251 * t280) * t297 + t140;
t130 = t138 * t300 - t271 * t139 + t150 * t299;
t188 = Ifges(6,5) * t215 + Ifges(6,6) * t214 + Ifges(6,3) * t225;
t190 = Ifges(6,1) * t215 + Ifges(6,4) * t214 + Ifges(6,5) * t225;
t159 = -mrSges(6,1) * t172 + mrSges(6,3) * t171 + Ifges(6,4) * t186 + Ifges(6,2) * t185 + Ifges(6,6) * t197 - t215 * t188 + t225 * t190;
t189 = Ifges(6,4) * t215 + Ifges(6,2) * t214 + Ifges(6,6) * t225;
t160 = mrSges(6,2) * t172 - mrSges(6,3) * t170 + Ifges(6,1) * t186 + Ifges(6,4) * t185 + Ifges(6,5) * t197 + t214 * t188 - t225 * t189;
t202 = Ifges(5,5) * t227 + Ifges(5,6) * t226 + Ifges(5,3) * t260;
t203 = Ifges(5,4) * t227 + Ifges(5,2) * t226 + Ifges(5,6) * t260;
t141 = mrSges(5,2) * t187 - mrSges(5,3) * t175 + Ifges(5,1) * t199 + Ifges(5,4) * t198 + Ifges(5,5) * t249 - pkin(9) * t155 - t274 * t159 + t278 * t160 + t226 * t202 - t260 * t203;
t204 = Ifges(5,1) * t227 + Ifges(5,4) * t226 + Ifges(5,5) * t260;
t285 = mrSges(6,1) * t170 - mrSges(6,2) * t171 + Ifges(6,5) * t186 + Ifges(6,6) * t185 + Ifges(6,3) * t197 + t215 * t189 - t214 * t190;
t142 = -mrSges(5,1) * t187 + mrSges(5,3) * t176 + Ifges(5,4) * t199 + Ifges(5,2) * t198 + Ifges(5,6) * t249 - pkin(4) * t155 - t227 * t202 + t260 * t204 - t285;
t218 = Ifges(4,5) * t245 + Ifges(4,6) * t244 - Ifges(4,3) * t293;
t220 = Ifges(4,1) * t245 + Ifges(4,4) * t244 - Ifges(4,5) * t293;
t131 = -mrSges(4,1) * t211 + mrSges(4,3) * t183 + Ifges(4,4) * t233 + Ifges(4,2) * t232 - Ifges(4,6) * t257 - pkin(3) * t287 + pkin(8) * t291 + t275 * t141 + t279 * t142 - t245 * t218 - t220 * t293;
t219 = Ifges(4,4) * t245 + Ifges(4,2) * t244 - Ifges(4,6) * t293;
t132 = mrSges(4,2) * t211 - mrSges(4,3) * t182 + Ifges(4,1) * t233 + Ifges(4,4) * t232 - Ifges(4,5) * t257 - pkin(8) * t148 + t279 * t141 - t275 * t142 + t244 * t218 + t219 * t293;
t236 = Ifges(3,6) * t266 + (Ifges(3,4) * t276 + Ifges(3,2) * t280) * t297;
t237 = Ifges(3,5) * t266 + (Ifges(3,1) * t276 + Ifges(3,4) * t280) * t297;
t122 = Ifges(3,5) * t256 + Ifges(3,6) * t257 + Ifges(3,3) * t265 + mrSges(3,1) * t223 - mrSges(3,2) * t224 + t270 * t132 + t272 * t131 + pkin(2) * t284 + qJ(3) * t292 + (t236 * t276 - t237 * t280) * t297;
t235 = Ifges(3,3) * t266 + (Ifges(3,5) * t276 + Ifges(3,6) * t280) * t297;
t124 = mrSges(3,2) * t238 - mrSges(3,3) * t223 + Ifges(3,1) * t256 + Ifges(3,4) * t257 + Ifges(3,5) * t265 - qJ(3) * t140 - t270 * t131 + t272 * t132 + t235 * t293 - t266 * t236;
t286 = -mrSges(5,1) * t175 + mrSges(5,2) * t176 - Ifges(5,5) * t199 - Ifges(5,6) * t198 - Ifges(5,3) * t249 - pkin(4) * t289 - pkin(9) * t290 - t278 * t159 - t274 * t160 - t227 * t203 + t226 * t204;
t283 = -mrSges(4,1) * t182 + mrSges(4,2) * t183 - Ifges(4,5) * t233 - Ifges(4,6) * t232 - pkin(3) * t148 - t245 * t219 + t244 * t220 + t286;
t126 = -pkin(2) * t140 + t283 + (Ifges(3,2) + Ifges(4,3)) * t257 - t235 * t294 + mrSges(3,3) * t224 - mrSges(3,1) * t238 + Ifges(3,4) * t256 + Ifges(3,6) * t265 + t266 * t237;
t288 = mrSges(2,1) * t261 - mrSges(2,2) * t262 + Ifges(2,3) * qJDD(1) + pkin(1) * t130 + t273 * t122 + t124 * t302 + t126 * t301 + t135 * t304;
t133 = m(2) * t262 - t282 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t135;
t129 = t273 * t139 + (t138 * t276 + t150 * t280) * t271;
t127 = m(2) * t261 + qJDD(1) * mrSges(2,1) - t282 * mrSges(2,2) + t130;
t120 = -mrSges(2,2) * g(3) - mrSges(2,3) * t261 + Ifges(2,5) * qJDD(1) - t282 * Ifges(2,6) + t280 * t124 - t276 * t126 + (-t129 * t271 - t130 * t273) * pkin(7);
t119 = mrSges(2,1) * g(3) + mrSges(2,3) * t262 + t282 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t129 - t271 * t122 + (pkin(7) * t135 + t124 * t276 + t126 * t280) * t273;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t281 * t120 - t277 * t119 - pkin(6) * (t281 * t127 + t277 * t133), t120, t124, t132, t141, t160; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t277 * t120 + t281 * t119 + pkin(6) * (-t277 * t127 + t281 * t133), t119, t126, t131, t142, t159; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t288, t288, t122, -Ifges(4,3) * t257 - t283, -t286, t285;];
m_new = t1;
