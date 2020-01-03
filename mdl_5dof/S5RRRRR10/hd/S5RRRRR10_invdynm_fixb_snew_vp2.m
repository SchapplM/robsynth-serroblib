% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:19
% EndTime: 2019-12-31 22:32:56
% DurationCPUTime: 19.41s
% Computational Cost: add. (344632->324), mult. (739808->425), div. (0->0), fcn. (584081->12), ass. (0->138)
t270 = sin(pkin(5));
t304 = pkin(7) * t270;
t271 = cos(pkin(5));
t303 = t271 * g(3);
t275 = sin(qJ(2));
t302 = t270 * t275;
t280 = cos(qJ(2));
t301 = t270 * t280;
t300 = t271 * t275;
t299 = t271 * t280;
t297 = qJD(1) * t270;
t254 = (-pkin(2) * t280 - pkin(8) * t275) * t297;
t266 = t271 * qJD(1) + qJD(2);
t264 = t266 ^ 2;
t265 = t271 * qJDD(1) + qJDD(2);
t296 = qJD(1) * t280;
t276 = sin(qJ(1));
t281 = cos(qJ(1));
t261 = t276 * g(1) - t281 * g(2);
t282 = qJD(1) ^ 2;
t251 = qJDD(1) * pkin(1) + t282 * t304 + t261;
t262 = -t281 * g(1) - t276 * g(2);
t295 = qJDD(1) * t270;
t252 = -t282 * pkin(1) + pkin(7) * t295 + t262;
t298 = t251 * t300 + t280 * t252;
t212 = -t264 * pkin(2) + t265 * pkin(8) + (-g(3) * t275 + t254 * t296) * t270 + t298;
t255 = (qJD(2) * t296 + qJDD(1) * t275) * t270;
t294 = t275 * t297;
t256 = -qJD(2) * t294 + t280 * t295;
t213 = -t256 * pkin(2) - t255 * pkin(8) - t303 + (-t251 + (pkin(2) * t275 - pkin(8) * t280) * t266 * qJD(1)) * t270;
t274 = sin(qJ(3));
t279 = cos(qJ(3));
t186 = -t274 * t212 + t279 * t213;
t243 = t279 * t266 - t274 * t294;
t225 = t243 * qJD(3) + t279 * t255 + t274 * t265;
t244 = t274 * t266 + t279 * t294;
t248 = qJDD(3) - t256;
t293 = t270 * t296;
t260 = qJD(3) - t293;
t179 = (t243 * t260 - t225) * pkin(9) + (t243 * t244 + t248) * pkin(3) + t186;
t187 = t279 * t212 + t274 * t213;
t224 = -t244 * qJD(3) - t274 * t255 + t279 * t265;
t234 = t260 * pkin(3) - t244 * pkin(9);
t242 = t243 ^ 2;
t181 = -t242 * pkin(3) + t224 * pkin(9) - t260 * t234 + t187;
t273 = sin(qJ(4));
t278 = cos(qJ(4));
t177 = t273 * t179 + t278 * t181;
t230 = t273 * t243 + t278 * t244;
t196 = -t230 * qJD(4) + t278 * t224 - t273 * t225;
t229 = t278 * t243 - t273 * t244;
t206 = -t229 * mrSges(5,1) + t230 * mrSges(5,2);
t258 = qJD(4) + t260;
t217 = t258 * mrSges(5,1) - t230 * mrSges(5,3);
t247 = qJDD(4) + t248;
t207 = -t229 * pkin(4) - t230 * pkin(10);
t257 = t258 ^ 2;
t173 = -t257 * pkin(4) + t247 * pkin(10) + t229 * t207 + t177;
t226 = -g(3) * t301 + t251 * t299 - t275 * t252;
t211 = -t265 * pkin(2) - t264 * pkin(8) + t254 * t294 - t226;
t185 = -t224 * pkin(3) - t242 * pkin(9) + t244 * t234 + t211;
t197 = t229 * qJD(4) + t273 * t224 + t278 * t225;
t174 = (-t229 * t258 - t197) * pkin(10) + (t230 * t258 - t196) * pkin(4) + t185;
t272 = sin(qJ(5));
t277 = cos(qJ(5));
t170 = -t272 * t173 + t277 * t174;
t214 = -t272 * t230 + t277 * t258;
t184 = t214 * qJD(5) + t277 * t197 + t272 * t247;
t195 = qJDD(5) - t196;
t215 = t277 * t230 + t272 * t258;
t199 = -t214 * mrSges(6,1) + t215 * mrSges(6,2);
t228 = qJD(5) - t229;
t200 = -t228 * mrSges(6,2) + t214 * mrSges(6,3);
t166 = m(6) * t170 + t195 * mrSges(6,1) - t184 * mrSges(6,3) - t215 * t199 + t228 * t200;
t171 = t277 * t173 + t272 * t174;
t183 = -t215 * qJD(5) - t272 * t197 + t277 * t247;
t201 = t228 * mrSges(6,1) - t215 * mrSges(6,3);
t167 = m(6) * t171 - t195 * mrSges(6,2) + t183 * mrSges(6,3) + t214 * t199 - t228 * t201;
t290 = -t272 * t166 + t277 * t167;
t153 = m(5) * t177 - t247 * mrSges(5,2) + t196 * mrSges(5,3) + t229 * t206 - t258 * t217 + t290;
t176 = t278 * t179 - t273 * t181;
t216 = -t258 * mrSges(5,2) + t229 * mrSges(5,3);
t172 = -t247 * pkin(4) - t257 * pkin(10) + t230 * t207 - t176;
t289 = -m(6) * t172 + t183 * mrSges(6,1) - t184 * mrSges(6,2) + t214 * t200 - t215 * t201;
t162 = m(5) * t176 + t247 * mrSges(5,1) - t197 * mrSges(5,3) - t230 * t206 + t258 * t216 + t289;
t148 = t273 * t153 + t278 * t162;
t231 = -t243 * mrSges(4,1) + t244 * mrSges(4,2);
t232 = -t260 * mrSges(4,2) + t243 * mrSges(4,3);
t146 = m(4) * t186 + t248 * mrSges(4,1) - t225 * mrSges(4,3) - t244 * t231 + t260 * t232 + t148;
t233 = t260 * mrSges(4,1) - t244 * mrSges(4,3);
t291 = t278 * t153 - t273 * t162;
t147 = m(4) * t187 - t248 * mrSges(4,2) + t224 * mrSges(4,3) + t243 * t231 - t260 * t233 + t291;
t140 = t279 * t146 + t274 * t147;
t155 = t277 * t166 + t272 * t167;
t227 = -g(3) * t302 + t298;
t249 = t266 * mrSges(3,1) - mrSges(3,3) * t294;
t253 = (-mrSges(3,1) * t280 + mrSges(3,2) * t275) * t297;
t292 = -t274 * t146 + t279 * t147;
t138 = m(3) * t227 - t265 * mrSges(3,2) + t256 * mrSges(3,3) - t266 * t249 + t253 * t293 + t292;
t250 = -t266 * mrSges(3,2) + mrSges(3,3) * t293;
t287 = m(5) * t185 - t196 * mrSges(5,1) + t197 * mrSges(5,2) - t229 * t216 + t230 * t217 + t155;
t284 = -m(4) * t211 + t224 * mrSges(4,1) - t225 * mrSges(4,2) + t243 * t232 - t244 * t233 - t287;
t150 = m(3) * t226 + t265 * mrSges(3,1) - t255 * mrSges(3,3) + t266 * t250 - t253 * t294 + t284;
t135 = t280 * t138 - t275 * t150;
t238 = -t270 * t251 - t303;
t139 = m(3) * t238 - t256 * mrSges(3,1) + t255 * mrSges(3,2) + (t249 * t275 - t250 * t280) * t297 + t140;
t130 = t138 * t300 - t270 * t139 + t150 * t299;
t188 = Ifges(6,5) * t215 + Ifges(6,6) * t214 + Ifges(6,3) * t228;
t190 = Ifges(6,1) * t215 + Ifges(6,4) * t214 + Ifges(6,5) * t228;
t159 = -mrSges(6,1) * t172 + mrSges(6,3) * t171 + Ifges(6,4) * t184 + Ifges(6,2) * t183 + Ifges(6,6) * t195 - t215 * t188 + t228 * t190;
t189 = Ifges(6,4) * t215 + Ifges(6,2) * t214 + Ifges(6,6) * t228;
t160 = mrSges(6,2) * t172 - mrSges(6,3) * t170 + Ifges(6,1) * t184 + Ifges(6,4) * t183 + Ifges(6,5) * t195 + t214 * t188 - t228 * t189;
t202 = Ifges(5,5) * t230 + Ifges(5,6) * t229 + Ifges(5,3) * t258;
t203 = Ifges(5,4) * t230 + Ifges(5,2) * t229 + Ifges(5,6) * t258;
t141 = mrSges(5,2) * t185 - mrSges(5,3) * t176 + Ifges(5,1) * t197 + Ifges(5,4) * t196 + Ifges(5,5) * t247 - pkin(10) * t155 - t272 * t159 + t277 * t160 + t229 * t202 - t258 * t203;
t204 = Ifges(5,1) * t230 + Ifges(5,4) * t229 + Ifges(5,5) * t258;
t285 = mrSges(6,1) * t170 - mrSges(6,2) * t171 + Ifges(6,5) * t184 + Ifges(6,6) * t183 + Ifges(6,3) * t195 + t215 * t189 - t214 * t190;
t142 = -mrSges(5,1) * t185 + mrSges(5,3) * t177 + Ifges(5,4) * t197 + Ifges(5,2) * t196 + Ifges(5,6) * t247 - pkin(4) * t155 - t230 * t202 + t258 * t204 - t285;
t218 = Ifges(4,5) * t244 + Ifges(4,6) * t243 + Ifges(4,3) * t260;
t220 = Ifges(4,1) * t244 + Ifges(4,4) * t243 + Ifges(4,5) * t260;
t131 = -mrSges(4,1) * t211 + mrSges(4,3) * t187 + Ifges(4,4) * t225 + Ifges(4,2) * t224 + Ifges(4,6) * t248 - pkin(3) * t287 + pkin(9) * t291 + t273 * t141 + t278 * t142 - t244 * t218 + t260 * t220;
t219 = Ifges(4,4) * t244 + Ifges(4,2) * t243 + Ifges(4,6) * t260;
t132 = mrSges(4,2) * t211 - mrSges(4,3) * t186 + Ifges(4,1) * t225 + Ifges(4,4) * t224 + Ifges(4,5) * t248 - pkin(9) * t148 + t278 * t141 - t273 * t142 + t243 * t218 - t260 * t219;
t236 = Ifges(3,6) * t266 + (Ifges(3,4) * t275 + Ifges(3,2) * t280) * t297;
t237 = Ifges(3,5) * t266 + (Ifges(3,1) * t275 + Ifges(3,4) * t280) * t297;
t122 = Ifges(3,5) * t255 + Ifges(3,6) * t256 + Ifges(3,3) * t265 + mrSges(3,1) * t226 - mrSges(3,2) * t227 + t274 * t132 + t279 * t131 + pkin(2) * t284 + pkin(8) * t292 + (t236 * t275 - t237 * t280) * t297;
t235 = Ifges(3,3) * t266 + (Ifges(3,5) * t275 + Ifges(3,6) * t280) * t297;
t124 = mrSges(3,2) * t238 - mrSges(3,3) * t226 + Ifges(3,1) * t255 + Ifges(3,4) * t256 + Ifges(3,5) * t265 - pkin(8) * t140 - t274 * t131 + t279 * t132 + t235 * t293 - t266 * t236;
t286 = -mrSges(5,1) * t176 + mrSges(5,2) * t177 - Ifges(5,5) * t197 - Ifges(5,6) * t196 - Ifges(5,3) * t247 - pkin(4) * t289 - pkin(10) * t290 - t277 * t159 - t272 * t160 - t230 * t203 + t229 * t204;
t283 = mrSges(4,1) * t186 - mrSges(4,2) * t187 + Ifges(4,5) * t225 + Ifges(4,6) * t224 + Ifges(4,3) * t248 + pkin(3) * t148 + t244 * t219 - t243 * t220 - t286;
t126 = -mrSges(3,1) * t238 + mrSges(3,3) * t227 + Ifges(3,4) * t255 + Ifges(3,2) * t256 + Ifges(3,6) * t265 - pkin(2) * t140 - t235 * t294 + t266 * t237 - t283;
t288 = mrSges(2,1) * t261 - mrSges(2,2) * t262 + Ifges(2,3) * qJDD(1) + pkin(1) * t130 + t271 * t122 + t124 * t302 + t126 * t301 + t135 * t304;
t133 = m(2) * t262 - t282 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t135;
t129 = t271 * t139 + (t138 * t275 + t150 * t280) * t270;
t127 = m(2) * t261 + qJDD(1) * mrSges(2,1) - t282 * mrSges(2,2) + t130;
t120 = -mrSges(2,2) * g(3) - mrSges(2,3) * t261 + Ifges(2,5) * qJDD(1) - t282 * Ifges(2,6) + t280 * t124 - t275 * t126 + (-t129 * t270 - t130 * t271) * pkin(7);
t119 = mrSges(2,1) * g(3) + mrSges(2,3) * t262 + t282 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t129 - t270 * t122 + (pkin(7) * t135 + t124 * t275 + t126 * t280) * t271;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t281 * t120 - t276 * t119 - pkin(6) * (t281 * t127 + t276 * t133), t120, t124, t132, t141, t160; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t276 * t120 + t281 * t119 + pkin(6) * (-t276 * t127 + t281 * t133), t119, t126, t131, t142, t159; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t288, t288, t122, t283, -t286, t285;];
m_new = t1;
