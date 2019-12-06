% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:01
% EndTime: 2019-12-05 18:38:19
% DurationCPUTime: 13.55s
% Computational Cost: add. (226298->312), mult. (520453->398), div. (0->0), fcn. (371609->10), ass. (0->124)
t270 = sin(qJ(2));
t274 = cos(qJ(2));
t292 = qJD(1) * qJD(2);
t248 = qJDD(1) * t270 + t274 * t292;
t271 = sin(qJ(1));
t275 = cos(qJ(1));
t255 = -g(1) * t275 - g(2) * t271;
t276 = qJD(1) ^ 2;
t243 = -pkin(1) * t276 + qJDD(1) * pkin(6) + t255;
t295 = t243 * t270;
t296 = pkin(2) * t276;
t208 = qJDD(2) * pkin(2) - pkin(7) * t248 - t295 + (pkin(7) * t292 + t270 * t296 - g(3)) * t274;
t231 = -g(3) * t270 + t274 * t243;
t249 = qJDD(1) * t274 - t270 * t292;
t294 = qJD(1) * t270;
t253 = qJD(2) * pkin(2) - pkin(7) * t294;
t265 = t274 ^ 2;
t209 = pkin(7) * t249 - qJD(2) * t253 - t265 * t296 + t231;
t269 = sin(qJ(3));
t273 = cos(qJ(3));
t186 = t273 * t208 - t209 * t269;
t240 = (-t269 * t270 + t273 * t274) * qJD(1);
t214 = qJD(3) * t240 + t248 * t273 + t249 * t269;
t241 = (t269 * t274 + t270 * t273) * qJD(1);
t262 = qJDD(2) + qJDD(3);
t263 = qJD(2) + qJD(3);
t173 = (t240 * t263 - t214) * qJ(4) + (t240 * t241 + t262) * pkin(3) + t186;
t187 = t269 * t208 + t273 * t209;
t213 = -qJD(3) * t241 - t248 * t269 + t249 * t273;
t233 = pkin(3) * t263 - qJ(4) * t241;
t236 = t240 ^ 2;
t175 = -pkin(3) * t236 + qJ(4) * t213 - t233 * t263 + t187;
t266 = sin(pkin(9));
t267 = cos(pkin(9));
t228 = t240 * t266 + t241 * t267;
t160 = -0.2e1 * qJD(4) * t228 + t267 * t173 - t175 * t266;
t193 = t213 * t266 + t214 * t267;
t227 = t240 * t267 - t241 * t266;
t157 = (t227 * t263 - t193) * pkin(8) + (t227 * t228 + t262) * pkin(4) + t160;
t161 = 0.2e1 * qJD(4) * t227 + t266 * t173 + t267 * t175;
t192 = t213 * t267 - t214 * t266;
t218 = pkin(4) * t263 - pkin(8) * t228;
t224 = t227 ^ 2;
t158 = -pkin(4) * t224 + pkin(8) * t192 - t218 * t263 + t161;
t268 = sin(qJ(5));
t272 = cos(qJ(5));
t155 = t157 * t272 - t158 * t268;
t201 = t227 * t272 - t228 * t268;
t169 = qJD(5) * t201 + t192 * t268 + t193 * t272;
t202 = t227 * t268 + t228 * t272;
t181 = -mrSges(6,1) * t201 + mrSges(6,2) * t202;
t260 = qJD(5) + t263;
t194 = -mrSges(6,2) * t260 + mrSges(6,3) * t201;
t259 = qJDD(5) + t262;
t152 = m(6) * t155 + mrSges(6,1) * t259 - mrSges(6,3) * t169 - t181 * t202 + t194 * t260;
t156 = t157 * t268 + t158 * t272;
t168 = -qJD(5) * t202 + t192 * t272 - t193 * t268;
t195 = mrSges(6,1) * t260 - mrSges(6,3) * t202;
t153 = m(6) * t156 - mrSges(6,2) * t259 + mrSges(6,3) * t168 + t181 * t201 - t195 * t260;
t144 = t272 * t152 + t268 * t153;
t203 = -mrSges(5,1) * t227 + mrSges(5,2) * t228;
t216 = -mrSges(5,2) * t263 + mrSges(5,3) * t227;
t141 = m(5) * t160 + mrSges(5,1) * t262 - mrSges(5,3) * t193 - t203 * t228 + t216 * t263 + t144;
t217 = mrSges(5,1) * t263 - mrSges(5,3) * t228;
t288 = -t152 * t268 + t272 * t153;
t142 = m(5) * t161 - mrSges(5,2) * t262 + mrSges(5,3) * t192 + t203 * t227 - t217 * t263 + t288;
t137 = t267 * t141 + t266 * t142;
t229 = -mrSges(4,1) * t240 + mrSges(4,2) * t241;
t232 = -mrSges(4,2) * t263 + mrSges(4,3) * t240;
t134 = m(4) * t186 + mrSges(4,1) * t262 - mrSges(4,3) * t214 - t229 * t241 + t232 * t263 + t137;
t234 = mrSges(4,1) * t263 - mrSges(4,3) * t241;
t289 = -t141 * t266 + t267 * t142;
t135 = m(4) * t187 - mrSges(4,2) * t262 + mrSges(4,3) * t213 + t229 * t240 - t234 * t263 + t289;
t128 = t273 * t134 + t269 * t135;
t230 = -g(3) * t274 - t295;
t238 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t270 + Ifges(3,2) * t274) * qJD(1);
t239 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t270 + Ifges(3,4) * t274) * qJD(1);
t222 = Ifges(4,4) * t241 + Ifges(4,2) * t240 + Ifges(4,6) * t263;
t223 = Ifges(4,1) * t241 + Ifges(4,4) * t240 + Ifges(4,5) * t263;
t197 = Ifges(5,4) * t228 + Ifges(5,2) * t227 + Ifges(5,6) * t263;
t198 = Ifges(5,1) * t228 + Ifges(5,4) * t227 + Ifges(5,5) * t263;
t177 = Ifges(6,4) * t202 + Ifges(6,2) * t201 + Ifges(6,6) * t260;
t178 = Ifges(6,1) * t202 + Ifges(6,4) * t201 + Ifges(6,5) * t260;
t282 = -mrSges(6,1) * t155 + mrSges(6,2) * t156 - Ifges(6,5) * t169 - Ifges(6,6) * t168 - Ifges(6,3) * t259 - t202 * t177 + t201 * t178;
t280 = -mrSges(5,1) * t160 + mrSges(5,2) * t161 - Ifges(5,5) * t193 - Ifges(5,6) * t192 - Ifges(5,3) * t262 - pkin(4) * t144 - t228 * t197 + t227 * t198 + t282;
t278 = -mrSges(4,1) * t186 + mrSges(4,2) * t187 - Ifges(4,5) * t214 - Ifges(4,6) * t213 - Ifges(4,3) * t262 - pkin(3) * t137 - t241 * t222 + t240 * t223 + t280;
t297 = mrSges(3,1) * t230 - mrSges(3,2) * t231 + Ifges(3,5) * t248 + Ifges(3,6) * t249 + Ifges(3,3) * qJDD(2) + pkin(2) * t128 + (t270 * t238 - t274 * t239) * qJD(1) - t278;
t293 = qJD(1) * t274;
t254 = g(1) * t271 - t275 * g(2);
t247 = (-mrSges(3,1) * t274 + mrSges(3,2) * t270) * qJD(1);
t252 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t293;
t126 = m(3) * t230 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t248 + qJD(2) * t252 - t247 * t294 + t128;
t251 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t294;
t290 = -t269 * t134 + t273 * t135;
t127 = m(3) * t231 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t249 - qJD(2) * t251 + t247 * t293 + t290;
t291 = -t126 * t270 + t274 * t127;
t285 = -qJDD(1) * pkin(1) - t254;
t215 = -pkin(2) * t249 + t253 * t294 + (-pkin(7) * t265 - pkin(6)) * t276 + t285;
t183 = -pkin(3) * t213 - qJ(4) * t236 + t241 * t233 + qJDD(4) + t215;
t163 = -pkin(4) * t192 - pkin(8) * t224 + t218 * t228 + t183;
t287 = m(6) * t163 - t168 * mrSges(6,1) + t169 * mrSges(6,2) - t201 * t194 + t202 * t195;
t176 = Ifges(6,5) * t202 + Ifges(6,6) * t201 + Ifges(6,3) * t260;
t145 = -mrSges(6,1) * t163 + mrSges(6,3) * t156 + Ifges(6,4) * t169 + Ifges(6,2) * t168 + Ifges(6,6) * t259 - t176 * t202 + t178 * t260;
t146 = mrSges(6,2) * t163 - mrSges(6,3) * t155 + Ifges(6,1) * t169 + Ifges(6,4) * t168 + Ifges(6,5) * t259 + t176 * t201 - t177 * t260;
t196 = Ifges(5,5) * t228 + Ifges(5,6) * t227 + Ifges(5,3) * t263;
t129 = -mrSges(5,1) * t183 + mrSges(5,3) * t161 + Ifges(5,4) * t193 + Ifges(5,2) * t192 + Ifges(5,6) * t262 - pkin(4) * t287 + pkin(8) * t288 + t272 * t145 + t268 * t146 - t228 * t196 + t263 * t198;
t130 = mrSges(5,2) * t183 - mrSges(5,3) * t160 + Ifges(5,1) * t193 + Ifges(5,4) * t192 + Ifges(5,5) * t262 - pkin(8) * t144 - t145 * t268 + t146 * t272 + t196 * t227 - t197 * t263;
t221 = Ifges(4,5) * t241 + Ifges(4,6) * t240 + Ifges(4,3) * t263;
t283 = m(5) * t183 - t192 * mrSges(5,1) + t193 * mrSges(5,2) - t227 * t216 + t228 * t217 + t287;
t123 = -mrSges(4,1) * t215 + mrSges(4,3) * t187 + Ifges(4,4) * t214 + Ifges(4,2) * t213 + Ifges(4,6) * t262 - pkin(3) * t283 + qJ(4) * t289 + t267 * t129 + t266 * t130 - t241 * t221 + t263 * t223;
t124 = mrSges(4,2) * t215 - mrSges(4,3) * t186 + Ifges(4,1) * t214 + Ifges(4,4) * t213 + Ifges(4,5) * t262 - qJ(4) * t137 - t129 * t266 + t130 * t267 + t221 * t240 - t222 * t263;
t237 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t270 + Ifges(3,6) * t274) * qJD(1);
t242 = -pkin(6) * t276 + t285;
t281 = m(4) * t215 - t213 * mrSges(4,1) + t214 * mrSges(4,2) - t240 * t232 + t241 * t234 + t283;
t116 = -mrSges(3,1) * t242 + mrSges(3,3) * t231 + Ifges(3,4) * t248 + Ifges(3,2) * t249 + Ifges(3,6) * qJDD(2) - pkin(2) * t281 + pkin(7) * t290 + qJD(2) * t239 + t273 * t123 + t269 * t124 - t237 * t294;
t118 = mrSges(3,2) * t242 - mrSges(3,3) * t230 + Ifges(3,1) * t248 + Ifges(3,4) * t249 + Ifges(3,5) * qJDD(2) - pkin(7) * t128 - qJD(2) * t238 - t123 * t269 + t124 * t273 + t237 * t293;
t279 = -m(3) * t242 + t249 * mrSges(3,1) - t248 * mrSges(3,2) - t251 * t294 + t252 * t293 - t281;
t284 = mrSges(2,1) * t254 - mrSges(2,2) * t255 + Ifges(2,3) * qJDD(1) + pkin(1) * t279 + pkin(6) * t291 + t274 * t116 + t270 * t118;
t147 = m(2) * t254 + qJDD(1) * mrSges(2,1) - t276 * mrSges(2,2) + t279;
t122 = t126 * t274 + t127 * t270;
t120 = m(2) * t255 - mrSges(2,1) * t276 - qJDD(1) * mrSges(2,2) + t291;
t119 = mrSges(2,1) * g(3) + mrSges(2,3) * t255 + t276 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t122 - t297;
t114 = -mrSges(2,2) * g(3) - mrSges(2,3) * t254 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t276 - pkin(6) * t122 - t116 * t270 + t118 * t274;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t275 * t114 - t271 * t119 - pkin(5) * (t120 * t271 + t147 * t275), t114, t118, t124, t130, t146; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t271 * t114 + t275 * t119 + pkin(5) * (t120 * t275 - t147 * t271), t119, t116, t123, t129, t145; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t284, t284, t297, -t278, -t280, -t282;];
m_new = t1;
