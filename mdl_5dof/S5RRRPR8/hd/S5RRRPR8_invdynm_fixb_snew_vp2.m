% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:13
% EndTime: 2019-12-31 21:19:23
% DurationCPUTime: 5.20s
% Computational Cost: add. (58272->313), mult. (121167->376), div. (0->0), fcn. (78765->8), ass. (0->119)
t263 = sin(qJ(2));
t266 = cos(qJ(2));
t287 = qJD(1) * qJD(2);
t243 = t263 * qJDD(1) + t266 * t287;
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t250 = -t267 * g(1) - t264 * g(2);
t268 = qJD(1) ^ 2;
t237 = -t268 * pkin(1) + qJDD(1) * pkin(6) + t250;
t292 = t263 * t237;
t295 = pkin(2) * t268;
t185 = qJDD(2) * pkin(2) - t243 * pkin(7) - t292 + (pkin(7) * t287 + t263 * t295 - g(3)) * t266;
t220 = -t263 * g(3) + t266 * t237;
t244 = t266 * qJDD(1) - t263 * t287;
t289 = qJD(1) * t263;
t248 = qJD(2) * pkin(2) - pkin(7) * t289;
t260 = t266 ^ 2;
t186 = t244 * pkin(7) - qJD(2) * t248 - t260 * t295 + t220;
t262 = sin(qJ(3));
t296 = cos(qJ(3));
t167 = t296 * t185 - t262 * t186;
t168 = t262 * t185 + t296 * t186;
t235 = (t262 * t266 + t296 * t263) * qJD(1);
t199 = t235 * qJD(3) + t262 * t243 - t296 * t244;
t288 = qJD(1) * t266;
t234 = t262 * t289 - t296 * t288;
t200 = -t234 * qJD(3) + t296 * t243 + t262 * t244;
t258 = qJD(2) + qJD(3);
t207 = Ifges(4,4) * t235 - Ifges(4,2) * t234 + Ifges(4,6) * t258;
t215 = -t234 * mrSges(5,2) - t235 * mrSges(5,3);
t223 = t234 * mrSges(5,1) - t258 * mrSges(5,3);
t257 = qJDD(2) + qJDD(3);
t225 = t235 * pkin(4) - t258 * pkin(8);
t230 = t234 ^ 2;
t249 = t264 * g(1) - t267 * g(2);
t282 = -qJDD(1) * pkin(1) - t249;
t201 = -t244 * pkin(2) + t248 * t289 + (-pkin(7) * t260 - pkin(6)) * t268 + t282;
t293 = t234 * t258;
t297 = -2 * qJD(4);
t272 = (-t200 + t293) * qJ(4) + t201 + (t258 * pkin(3) + t297) * t235;
t157 = -t230 * pkin(4) - t235 * t225 + (pkin(3) + pkin(8)) * t199 + t272;
t213 = t234 * pkin(3) - t235 * qJ(4);
t256 = t258 ^ 2;
t165 = -t257 * pkin(3) - t256 * qJ(4) + t235 * t213 + qJDD(4) - t167;
t158 = (t234 * t235 - t257) * pkin(8) + (t200 + t293) * pkin(4) + t165;
t261 = sin(qJ(5));
t265 = cos(qJ(5));
t154 = -t261 * t157 + t265 * t158;
t217 = t265 * t234 - t261 * t258;
t174 = t217 * qJD(5) + t261 * t199 + t265 * t257;
t218 = t261 * t234 + t265 * t258;
t182 = -t217 * mrSges(6,1) + t218 * mrSges(6,2);
t198 = qJDD(5) + t200;
t229 = qJD(5) + t235;
t202 = -t229 * mrSges(6,2) + t217 * mrSges(6,3);
t151 = m(6) * t154 + t198 * mrSges(6,1) - t174 * mrSges(6,3) - t218 * t182 + t229 * t202;
t155 = t265 * t157 + t261 * t158;
t173 = -t218 * qJD(5) + t265 * t199 - t261 * t257;
t203 = t229 * mrSges(6,1) - t218 * mrSges(6,3);
t152 = m(6) * t155 - t198 * mrSges(6,2) + t173 * mrSges(6,3) + t217 * t182 - t229 * t203;
t139 = t265 * t151 + t261 * t152;
t279 = -t256 * pkin(3) + t257 * qJ(4) - t234 * t213 + t168;
t160 = -t199 * pkin(4) - t230 * pkin(8) + ((2 * qJD(4)) + t225) * t258 + t279;
t175 = Ifges(6,5) * t218 + Ifges(6,6) * t217 + Ifges(6,3) * t229;
t177 = Ifges(6,1) * t218 + Ifges(6,4) * t217 + Ifges(6,5) * t229;
t142 = -mrSges(6,1) * t160 + mrSges(6,3) * t155 + Ifges(6,4) * t174 + Ifges(6,2) * t173 + Ifges(6,6) * t198 - t218 * t175 + t229 * t177;
t176 = Ifges(6,4) * t218 + Ifges(6,2) * t217 + Ifges(6,6) * t229;
t143 = mrSges(6,2) * t160 - mrSges(6,3) * t154 + Ifges(6,1) * t174 + Ifges(6,4) * t173 + Ifges(6,5) * t198 + t217 * t175 - t229 * t176;
t163 = t258 * t297 - t279;
t204 = Ifges(5,5) * t258 - Ifges(5,6) * t235 + Ifges(5,3) * t234;
t276 = -mrSges(5,2) * t165 + mrSges(5,3) * t163 - Ifges(5,1) * t257 + Ifges(5,4) * t200 - Ifges(5,5) * t199 + pkin(8) * t139 + t261 * t142 - t265 * t143 + t235 * t204;
t156 = -m(6) * t160 + t173 * mrSges(6,1) - t174 * mrSges(6,2) + t217 * t202 - t218 * t203;
t224 = t235 * mrSges(5,1) + t258 * mrSges(5,2);
t277 = -m(5) * t163 + t257 * mrSges(5,3) + t258 * t224 - t156;
t280 = -m(5) * t165 - t200 * mrSges(5,1) - t235 * t215 - t139;
t206 = Ifges(5,4) * t258 - Ifges(5,2) * t235 + Ifges(5,6) * t234;
t290 = Ifges(4,1) * t235 - Ifges(4,4) * t234 + Ifges(4,5) * t258 - t206;
t300 = -mrSges(4,2) * t168 + pkin(3) * (-t257 * mrSges(5,2) - t258 * t223 + t280) + qJ(4) * (-t199 * mrSges(5,1) - t234 * t215 + t277) + mrSges(4,1) * t167 + t235 * t207 - Ifges(4,6) * t199 + Ifges(4,5) * t200 + Ifges(4,3) * t257 - t276 + t290 * t234;
t214 = t234 * mrSges(4,1) + t235 * mrSges(4,2);
t221 = -t258 * mrSges(4,2) - t234 * mrSges(4,3);
t135 = m(4) * t167 - t200 * mrSges(4,3) - t235 * t214 + (t221 - t223) * t258 + (mrSges(4,1) - mrSges(5,2)) * t257 + t280;
t222 = t258 * mrSges(4,1) - t235 * mrSges(4,3);
t146 = m(4) * t168 - t257 * mrSges(4,2) - t258 * t222 + (-t214 - t215) * t234 + (-mrSges(4,3) - mrSges(5,1)) * t199 + t277;
t131 = t296 * t135 + t262 * t146;
t219 = -t266 * g(3) - t292;
t232 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t263 + Ifges(3,2) * t266) * qJD(1);
t233 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t263 + Ifges(3,4) * t266) * qJD(1);
t299 = mrSges(3,1) * t219 - mrSges(3,2) * t220 + Ifges(3,5) * t243 + Ifges(3,6) * t244 + Ifges(3,3) * qJDD(2) + pkin(2) * t131 + (t263 * t232 - t266 * t233) * qJD(1) + t300;
t294 = Ifges(4,4) + Ifges(5,6);
t140 = -t261 * t151 + t265 * t152;
t208 = Ifges(5,1) * t258 - Ifges(5,4) * t235 + Ifges(5,5) * t234;
t291 = -Ifges(4,5) * t235 + Ifges(4,6) * t234 - Ifges(4,3) * t258 - t208;
t242 = (-mrSges(3,1) * t266 + mrSges(3,2) * t263) * qJD(1);
t247 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t288;
t129 = m(3) * t219 + qJDD(2) * mrSges(3,1) - t243 * mrSges(3,3) + qJD(2) * t247 - t242 * t289 + t131;
t246 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t289;
t284 = -t262 * t135 + t296 * t146;
t130 = m(3) * t220 - qJDD(2) * mrSges(3,2) + t244 * mrSges(3,3) - qJD(2) * t246 + t242 * t288 + t284;
t285 = -t263 * t129 + t266 * t130;
t162 = t199 * pkin(3) + t272;
t136 = m(5) * t162 - t199 * mrSges(5,2) - t200 * mrSges(5,3) - t234 * t223 - t235 * t224 + t140;
t275 = -mrSges(5,1) * t163 + mrSges(5,2) * t162 - pkin(4) * t156 - pkin(8) * t140 - t265 * t142 - t261 * t143;
t123 = -mrSges(4,1) * t201 + mrSges(4,3) * t168 - pkin(3) * t136 + t290 * t258 + (Ifges(4,6) - Ifges(5,5)) * t257 + t291 * t235 + t294 * t200 + (-Ifges(4,2) - Ifges(5,3)) * t199 + t275;
t278 = mrSges(6,1) * t154 - mrSges(6,2) * t155 + Ifges(6,5) * t174 + Ifges(6,6) * t173 + Ifges(6,3) * t198 + t218 * t176 - t217 * t177;
t273 = mrSges(5,1) * t165 - mrSges(5,3) * t162 + pkin(4) * t139 + t278;
t127 = t273 + (-t207 + t204) * t258 + (Ifges(4,5) - Ifges(5,4)) * t257 + t291 * t234 + (Ifges(4,1) + Ifges(5,2)) * t200 - t294 * t199 + mrSges(4,2) * t201 - mrSges(4,3) * t167 - qJ(4) * t136;
t231 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t263 + Ifges(3,6) * t266) * qJD(1);
t236 = -t268 * pkin(6) + t282;
t274 = m(4) * t201 + t199 * mrSges(4,1) + t200 * mrSges(4,2) + t234 * t221 + t235 * t222 + t136;
t119 = -mrSges(3,1) * t236 + mrSges(3,3) * t220 + Ifges(3,4) * t243 + Ifges(3,2) * t244 + Ifges(3,6) * qJDD(2) - pkin(2) * t274 + pkin(7) * t284 + qJD(2) * t233 + t296 * t123 + t262 * t127 - t231 * t289;
t122 = mrSges(3,2) * t236 - mrSges(3,3) * t219 + Ifges(3,1) * t243 + Ifges(3,4) * t244 + Ifges(3,5) * qJDD(2) - pkin(7) * t131 - qJD(2) * t232 - t262 * t123 + t296 * t127 + t231 * t288;
t270 = -m(3) * t236 + t244 * mrSges(3,1) - t243 * mrSges(3,2) - t246 * t289 + t247 * t288 - t274;
t281 = mrSges(2,1) * t249 - mrSges(2,2) * t250 + Ifges(2,3) * qJDD(1) + pkin(1) * t270 + pkin(6) * t285 + t266 * t119 + t263 * t122;
t132 = m(2) * t249 + qJDD(1) * mrSges(2,1) - t268 * mrSges(2,2) + t270;
t126 = t266 * t129 + t263 * t130;
t124 = m(2) * t250 - t268 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t285;
t120 = mrSges(2,1) * g(3) + mrSges(2,3) * t250 + t268 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t126 - t299;
t117 = -mrSges(2,2) * g(3) - mrSges(2,3) * t249 + Ifges(2,5) * qJDD(1) - t268 * Ifges(2,6) - pkin(6) * t126 - t263 * t119 + t266 * t122;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t267 * t117 - t264 * t120 - pkin(5) * (t264 * t124 + t267 * t132), t117, t122, t127, -t234 * t206 - t276, t143; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t117 + t267 * t120 + pkin(5) * (t267 * t124 - t264 * t132), t120, t119, t123, Ifges(5,4) * t257 - Ifges(5,2) * t200 + Ifges(5,6) * t199 - t258 * t204 + t234 * t208 - t273, t142; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t281, t281, t299, t300, Ifges(5,5) * t257 - Ifges(5,6) * t200 + Ifges(5,3) * t199 + t258 * t206 + t235 * t208 - t275, t278;];
m_new = t1;
