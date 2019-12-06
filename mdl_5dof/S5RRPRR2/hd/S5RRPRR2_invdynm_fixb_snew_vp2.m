% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:22
% EndTime: 2019-12-05 18:27:46
% DurationCPUTime: 13.54s
% Computational Cost: add. (212504->312), mult. (506882->398), div. (0->0), fcn. (361789->10), ass. (0->124)
t270 = sin(qJ(2));
t274 = cos(qJ(2));
t292 = qJD(1) * qJD(2);
t248 = t270 * qJDD(1) + t274 * t292;
t271 = sin(qJ(1));
t275 = cos(qJ(1));
t255 = -t275 * g(1) - t271 * g(2);
t276 = qJD(1) ^ 2;
t243 = -t276 * pkin(1) + qJDD(1) * pkin(6) + t255;
t295 = t270 * t243;
t296 = pkin(2) * t276;
t208 = qJDD(2) * pkin(2) - t248 * qJ(3) - t295 + (qJ(3) * t292 + t270 * t296 - g(3)) * t274;
t229 = -t270 * g(3) + t274 * t243;
t249 = t274 * qJDD(1) - t270 * t292;
t294 = qJD(1) * t270;
t251 = qJD(2) * pkin(2) - qJ(3) * t294;
t265 = t274 ^ 2;
t209 = t249 * qJ(3) - qJD(2) * t251 - t265 * t296 + t229;
t266 = sin(pkin(9));
t267 = cos(pkin(9));
t238 = (t266 * t274 + t267 * t270) * qJD(1);
t185 = -0.2e1 * qJD(3) * t238 + t267 * t208 - t266 * t209;
t227 = t267 * t248 + t266 * t249;
t237 = (-t266 * t270 + t267 * t274) * qJD(1);
t173 = (qJD(2) * t237 - t227) * pkin(7) + (t237 * t238 + qJDD(2)) * pkin(3) + t185;
t186 = 0.2e1 * qJD(3) * t237 + t266 * t208 + t267 * t209;
t226 = -t266 * t248 + t267 * t249;
t232 = qJD(2) * pkin(3) - t238 * pkin(7);
t236 = t237 ^ 2;
t175 = -t236 * pkin(3) + t226 * pkin(7) - qJD(2) * t232 + t186;
t269 = sin(qJ(4));
t273 = cos(qJ(4));
t160 = t273 * t173 - t269 * t175;
t218 = t273 * t237 - t269 * t238;
t193 = t218 * qJD(4) + t269 * t226 + t273 * t227;
t219 = t269 * t237 + t273 * t238;
t261 = qJDD(2) + qJDD(4);
t262 = qJD(2) + qJD(4);
t157 = (t218 * t262 - t193) * pkin(8) + (t218 * t219 + t261) * pkin(4) + t160;
t161 = t269 * t173 + t273 * t175;
t192 = -t219 * qJD(4) + t273 * t226 - t269 * t227;
t213 = t262 * pkin(4) - t219 * pkin(8);
t214 = t218 ^ 2;
t158 = -t214 * pkin(4) + t192 * pkin(8) - t262 * t213 + t161;
t268 = sin(qJ(5));
t272 = cos(qJ(5));
t155 = t272 * t157 - t268 * t158;
t201 = t272 * t218 - t268 * t219;
t169 = t201 * qJD(5) + t268 * t192 + t272 * t193;
t202 = t268 * t218 + t272 * t219;
t181 = -t201 * mrSges(6,1) + t202 * mrSges(6,2);
t259 = qJD(5) + t262;
t194 = -t259 * mrSges(6,2) + t201 * mrSges(6,3);
t258 = qJDD(5) + t261;
t152 = m(6) * t155 + t258 * mrSges(6,1) - t169 * mrSges(6,3) - t202 * t181 + t259 * t194;
t156 = t268 * t157 + t272 * t158;
t168 = -t202 * qJD(5) + t272 * t192 - t268 * t193;
t195 = t259 * mrSges(6,1) - t202 * mrSges(6,3);
t153 = m(6) * t156 - t258 * mrSges(6,2) + t168 * mrSges(6,3) + t201 * t181 - t259 * t195;
t144 = t272 * t152 + t268 * t153;
t203 = -t218 * mrSges(5,1) + t219 * mrSges(5,2);
t211 = -t262 * mrSges(5,2) + t218 * mrSges(5,3);
t141 = m(5) * t160 + t261 * mrSges(5,1) - t193 * mrSges(5,3) - t219 * t203 + t262 * t211 + t144;
t212 = t262 * mrSges(5,1) - t219 * mrSges(5,3);
t288 = -t268 * t152 + t272 * t153;
t142 = m(5) * t161 - t261 * mrSges(5,2) + t192 * mrSges(5,3) + t218 * t203 - t262 * t212 + t288;
t137 = t273 * t141 + t269 * t142;
t222 = -t237 * mrSges(4,1) + t238 * mrSges(4,2);
t230 = -qJD(2) * mrSges(4,2) + t237 * mrSges(4,3);
t134 = m(4) * t185 + qJDD(2) * mrSges(4,1) - t227 * mrSges(4,3) + qJD(2) * t230 - t238 * t222 + t137;
t231 = qJD(2) * mrSges(4,1) - t238 * mrSges(4,3);
t289 = -t269 * t141 + t273 * t142;
t135 = m(4) * t186 - qJDD(2) * mrSges(4,2) + t226 * mrSges(4,3) - qJD(2) * t231 + t237 * t222 + t289;
t128 = t267 * t134 + t266 * t135;
t228 = -t274 * g(3) - t295;
t240 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t270 + Ifges(3,2) * t274) * qJD(1);
t241 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t270 + Ifges(3,4) * t274) * qJD(1);
t216 = Ifges(4,4) * t238 + Ifges(4,2) * t237 + Ifges(4,6) * qJD(2);
t217 = Ifges(4,1) * t238 + Ifges(4,4) * t237 + Ifges(4,5) * qJD(2);
t197 = Ifges(5,4) * t219 + Ifges(5,2) * t218 + Ifges(5,6) * t262;
t198 = Ifges(5,1) * t219 + Ifges(5,4) * t218 + Ifges(5,5) * t262;
t177 = Ifges(6,4) * t202 + Ifges(6,2) * t201 + Ifges(6,6) * t259;
t178 = Ifges(6,1) * t202 + Ifges(6,4) * t201 + Ifges(6,5) * t259;
t282 = -mrSges(6,1) * t155 + mrSges(6,2) * t156 - Ifges(6,5) * t169 - Ifges(6,6) * t168 - Ifges(6,3) * t258 - t202 * t177 + t201 * t178;
t280 = -mrSges(5,1) * t160 + mrSges(5,2) * t161 - Ifges(5,5) * t193 - Ifges(5,6) * t192 - Ifges(5,3) * t261 - pkin(4) * t144 - t219 * t197 + t218 * t198 + t282;
t278 = -mrSges(4,1) * t185 + mrSges(4,2) * t186 - Ifges(4,5) * t227 - Ifges(4,6) * t226 - Ifges(4,3) * qJDD(2) - pkin(3) * t137 - t238 * t216 + t237 * t217 + t280;
t297 = mrSges(3,1) * t228 - mrSges(3,2) * t229 + Ifges(3,5) * t248 + Ifges(3,6) * t249 + Ifges(3,3) * qJDD(2) + pkin(2) * t128 + (t270 * t240 - t274 * t241) * qJD(1) - t278;
t293 = qJD(1) * t274;
t254 = t271 * g(1) - t275 * g(2);
t247 = (-mrSges(3,1) * t274 + mrSges(3,2) * t270) * qJD(1);
t253 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t293;
t126 = m(3) * t228 + qJDD(2) * mrSges(3,1) - t248 * mrSges(3,3) + qJD(2) * t253 - t247 * t294 + t128;
t252 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t294;
t290 = -t266 * t134 + t267 * t135;
t127 = m(3) * t229 - qJDD(2) * mrSges(3,2) + t249 * mrSges(3,3) - qJD(2) * t252 + t247 * t293 + t290;
t291 = -t270 * t126 + t274 * t127;
t285 = -qJDD(1) * pkin(1) - t254;
t210 = -t249 * pkin(2) + qJDD(3) + t251 * t294 + (-qJ(3) * t265 - pkin(6)) * t276 + t285;
t183 = -t226 * pkin(3) - t236 * pkin(7) + t238 * t232 + t210;
t163 = -t192 * pkin(4) - t214 * pkin(8) + t219 * t213 + t183;
t287 = m(6) * t163 - t168 * mrSges(6,1) + t169 * mrSges(6,2) - t201 * t194 + t202 * t195;
t176 = Ifges(6,5) * t202 + Ifges(6,6) * t201 + Ifges(6,3) * t259;
t145 = -mrSges(6,1) * t163 + mrSges(6,3) * t156 + Ifges(6,4) * t169 + Ifges(6,2) * t168 + Ifges(6,6) * t258 - t202 * t176 + t259 * t178;
t146 = mrSges(6,2) * t163 - mrSges(6,3) * t155 + Ifges(6,1) * t169 + Ifges(6,4) * t168 + Ifges(6,5) * t258 + t201 * t176 - t259 * t177;
t196 = Ifges(5,5) * t219 + Ifges(5,6) * t218 + Ifges(5,3) * t262;
t129 = -mrSges(5,1) * t183 + mrSges(5,3) * t161 + Ifges(5,4) * t193 + Ifges(5,2) * t192 + Ifges(5,6) * t261 - pkin(4) * t287 + pkin(8) * t288 + t272 * t145 + t268 * t146 - t219 * t196 + t262 * t198;
t130 = mrSges(5,2) * t183 - mrSges(5,3) * t160 + Ifges(5,1) * t193 + Ifges(5,4) * t192 + Ifges(5,5) * t261 - pkin(8) * t144 - t268 * t145 + t272 * t146 + t218 * t196 - t262 * t197;
t215 = Ifges(4,5) * t238 + Ifges(4,6) * t237 + Ifges(4,3) * qJD(2);
t283 = m(5) * t183 - t192 * mrSges(5,1) + t193 * mrSges(5,2) - t218 * t211 + t219 * t212 + t287;
t123 = -mrSges(4,1) * t210 + mrSges(4,3) * t186 + Ifges(4,4) * t227 + Ifges(4,2) * t226 + Ifges(4,6) * qJDD(2) - pkin(3) * t283 + pkin(7) * t289 + qJD(2) * t217 + t273 * t129 + t269 * t130 - t238 * t215;
t124 = mrSges(4,2) * t210 - mrSges(4,3) * t185 + Ifges(4,1) * t227 + Ifges(4,4) * t226 + Ifges(4,5) * qJDD(2) - pkin(7) * t137 - qJD(2) * t216 - t269 * t129 + t273 * t130 + t237 * t215;
t239 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t270 + Ifges(3,6) * t274) * qJD(1);
t242 = -t276 * pkin(6) + t285;
t281 = m(4) * t210 - t226 * mrSges(4,1) + t227 * mrSges(4,2) - t237 * t230 + t238 * t231 + t283;
t116 = -mrSges(3,1) * t242 + mrSges(3,3) * t229 + Ifges(3,4) * t248 + Ifges(3,2) * t249 + Ifges(3,6) * qJDD(2) - pkin(2) * t281 + qJ(3) * t290 + qJD(2) * t241 + t267 * t123 + t266 * t124 - t239 * t294;
t118 = mrSges(3,2) * t242 - mrSges(3,3) * t228 + Ifges(3,1) * t248 + Ifges(3,4) * t249 + Ifges(3,5) * qJDD(2) - qJ(3) * t128 - qJD(2) * t240 - t266 * t123 + t267 * t124 + t239 * t293;
t279 = -m(3) * t242 + t249 * mrSges(3,1) - t248 * mrSges(3,2) - t252 * t294 + t253 * t293 - t281;
t284 = mrSges(2,1) * t254 - mrSges(2,2) * t255 + Ifges(2,3) * qJDD(1) + pkin(1) * t279 + pkin(6) * t291 + t274 * t116 + t270 * t118;
t147 = m(2) * t254 + qJDD(1) * mrSges(2,1) - t276 * mrSges(2,2) + t279;
t122 = t274 * t126 + t270 * t127;
t120 = m(2) * t255 - t276 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t291;
t119 = mrSges(2,1) * g(3) + mrSges(2,3) * t255 + t276 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t122 - t297;
t114 = -mrSges(2,2) * g(3) - mrSges(2,3) * t254 + Ifges(2,5) * qJDD(1) - t276 * Ifges(2,6) - pkin(6) * t122 - t270 * t116 + t274 * t118;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t275 * t114 - t271 * t119 - pkin(5) * (t271 * t120 + t275 * t147), t114, t118, t124, t130, t146; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t271 * t114 + t275 * t119 + pkin(5) * (t275 * t120 - t271 * t147), t119, t116, t123, t129, t145; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t284, t284, t297, -t278, -t280, -t282;];
m_new = t1;
