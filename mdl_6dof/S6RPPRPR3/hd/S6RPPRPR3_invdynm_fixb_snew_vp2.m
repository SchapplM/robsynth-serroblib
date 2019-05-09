% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 14:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:09:00
% EndTime: 2019-05-05 14:09:12
% DurationCPUTime: 7.15s
% Computational Cost: add. (100308->300), mult. (203919->368), div. (0->0), fcn. (123185->10), ass. (0->124)
t268 = sin(qJ(1));
t271 = cos(qJ(1));
t244 = t268 * g(1) - t271 * g(2);
t235 = qJDD(1) * pkin(1) + t244;
t245 = -t271 * g(1) - t268 * g(2);
t273 = qJD(1) ^ 2;
t237 = -t273 * pkin(1) + t245;
t263 = sin(pkin(9));
t265 = cos(pkin(9));
t210 = t263 * t235 + t265 * t237;
t289 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t210;
t262 = sin(pkin(10));
t264 = cos(pkin(10));
t267 = sin(qJ(4));
t270 = cos(qJ(4));
t221 = (t262 * t270 + t264 * t267) * qJD(1);
t303 = 2 * qJD(5);
t302 = -pkin(2) - pkin(7);
t301 = pkin(4) * t273;
t300 = mrSges(3,1) - mrSges(4,2);
t299 = Ifges(3,5) - Ifges(4,4);
t298 = -Ifges(3,6) + Ifges(4,5);
t209 = t265 * t235 - t263 * t237;
t286 = -t273 * qJ(3) + qJDD(3) - t209;
t194 = t302 * qJDD(1) + t286;
t189 = t270 * t194;
t295 = qJD(1) * qJD(4);
t239 = t270 * qJDD(1) - t267 * t295;
t259 = -g(3) + qJDD(2);
t172 = qJDD(4) * pkin(4) - t239 * qJ(5) + t189 + (-qJ(5) * t295 - t270 * t301 - t259) * t267;
t186 = t267 * t194 + t270 * t259;
t238 = -t267 * qJDD(1) - t270 * t295;
t296 = qJD(1) * t270;
t242 = qJD(4) * pkin(4) - qJ(5) * t296;
t258 = t267 ^ 2;
t173 = t238 * qJ(5) - qJD(4) * t242 - t258 * t301 + t186;
t168 = t262 * t172 + t264 * t173 - t221 * t303;
t297 = qJD(1) * t267;
t222 = -t262 * t297 + t264 * t296;
t202 = mrSges(6,1) * t221 + mrSges(6,2) * t222;
t211 = t264 * t238 - t262 * t239;
t216 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t222;
t203 = pkin(5) * t221 - pkin(8) * t222;
t272 = qJD(4) ^ 2;
t165 = -t272 * pkin(5) + qJDD(4) * pkin(8) - t203 * t221 + t168;
t175 = -t238 * pkin(4) + qJDD(5) + t242 * t296 + (-qJ(5) * t258 + t302) * t273 + t289;
t212 = t262 * t238 + t264 * t239;
t169 = (qJD(4) * t221 - t212) * pkin(8) + (qJD(4) * t222 - t211) * pkin(5) + t175;
t266 = sin(qJ(6));
t269 = cos(qJ(6));
t162 = -t266 * t165 + t269 * t169;
t213 = t269 * qJD(4) - t266 * t222;
t182 = t213 * qJD(6) + t266 * qJDD(4) + t269 * t212;
t214 = t266 * qJD(4) + t269 * t222;
t187 = -mrSges(7,1) * t213 + mrSges(7,2) * t214;
t220 = qJD(6) + t221;
t191 = -mrSges(7,2) * t220 + mrSges(7,3) * t213;
t208 = qJDD(6) - t211;
t158 = m(7) * t162 + mrSges(7,1) * t208 - mrSges(7,3) * t182 - t187 * t214 + t191 * t220;
t163 = t269 * t165 + t266 * t169;
t181 = -t214 * qJD(6) + t269 * qJDD(4) - t266 * t212;
t192 = mrSges(7,1) * t220 - mrSges(7,3) * t214;
t159 = m(7) * t163 - mrSges(7,2) * t208 + mrSges(7,3) * t181 + t187 * t213 - t192 * t220;
t290 = -t266 * t158 + t269 * t159;
t145 = m(6) * t168 - qJDD(4) * mrSges(6,2) + t211 * mrSges(6,3) - qJD(4) * t216 - t221 * t202 + t290;
t288 = -t264 * t172 + t262 * t173;
t167 = -0.2e1 * qJD(5) * t222 - t288;
t215 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t221;
t164 = -qJDD(4) * pkin(5) - t272 * pkin(8) + (t303 + t203) * t222 + t288;
t284 = -m(7) * t164 + t181 * mrSges(7,1) - mrSges(7,2) * t182 + t213 * t191 - t192 * t214;
t154 = m(6) * t167 + qJDD(4) * mrSges(6,1) - mrSges(6,3) * t212 + qJD(4) * t215 - t202 * t222 + t284;
t137 = t262 * t145 + t264 * t154;
t185 = -t267 * t259 + t189;
t236 = (mrSges(5,1) * t267 + mrSges(5,2) * t270) * qJD(1);
t241 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t297;
t134 = m(5) * t185 + qJDD(4) * mrSges(5,1) - t239 * mrSges(5,3) + qJD(4) * t241 - t236 * t296 + t137;
t243 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t296;
t291 = t264 * t145 - t262 * t154;
t135 = m(5) * t186 - qJDD(4) * mrSges(5,2) + t238 * mrSges(5,3) - qJD(4) * t243 - t236 * t297 + t291;
t129 = t270 * t134 + t267 * t135;
t197 = -qJDD(1) * pkin(2) + t286;
t285 = -m(4) * t197 + t273 * mrSges(4,3) - t129;
t126 = m(3) * t209 - t273 * mrSges(3,2) + t300 * qJDD(1) + t285;
t193 = t302 * t273 + t289;
t147 = t269 * t158 + t266 * t159;
t283 = m(6) * t175 - mrSges(6,1) * t211 + t212 * mrSges(6,2) + t215 * t221 + t222 * t216 + t147;
t142 = -m(5) * t193 + mrSges(5,1) * t238 - t239 * mrSges(5,2) - t241 * t297 - t243 * t296 - t283;
t195 = t273 * pkin(2) - t289;
t277 = -m(4) * t195 + t273 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t142;
t140 = m(3) * t210 - t273 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t277;
t124 = t265 * t126 + t263 * t140;
t292 = -t263 * t126 + t265 * t140;
t130 = -t267 * t134 + t270 * t135;
t128 = m(4) * t259 + t130;
t176 = Ifges(7,5) * t214 + Ifges(7,6) * t213 + Ifges(7,3) * t220;
t178 = Ifges(7,1) * t214 + Ifges(7,4) * t213 + Ifges(7,5) * t220;
t151 = -mrSges(7,1) * t164 + mrSges(7,3) * t163 + Ifges(7,4) * t182 + Ifges(7,2) * t181 + Ifges(7,6) * t208 - t176 * t214 + t178 * t220;
t177 = Ifges(7,4) * t214 + Ifges(7,2) * t213 + Ifges(7,6) * t220;
t152 = mrSges(7,2) * t164 - mrSges(7,3) * t162 + Ifges(7,1) * t182 + Ifges(7,4) * t181 + Ifges(7,5) * t208 + t176 * t213 - t177 * t220;
t198 = Ifges(6,5) * t222 - Ifges(6,6) * t221 + (Ifges(6,3) * qJD(4));
t199 = Ifges(6,4) * t222 - Ifges(6,2) * t221 + Ifges(6,6) * qJD(4);
t131 = mrSges(6,2) * t175 - mrSges(6,3) * t167 + Ifges(6,1) * t212 + Ifges(6,4) * t211 + Ifges(6,5) * qJDD(4) - pkin(8) * t147 - qJD(4) * t199 - t266 * t151 + t269 * t152 - t221 * t198;
t200 = Ifges(6,1) * t222 - Ifges(6,4) * t221 + Ifges(6,5) * qJD(4);
t278 = mrSges(7,1) * t162 - mrSges(7,2) * t163 + Ifges(7,5) * t182 + Ifges(7,6) * t181 + Ifges(7,3) * t208 + t177 * t214 - t178 * t213;
t132 = -mrSges(6,1) * t175 + mrSges(6,3) * t168 + Ifges(6,4) * t212 + Ifges(6,2) * t211 + Ifges(6,6) * qJDD(4) - pkin(5) * t147 + qJD(4) * t200 - t198 * t222 - t278;
t226 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t270 - Ifges(5,6) * t267) * qJD(1);
t228 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t270 - Ifges(5,4) * t267) * qJD(1);
t118 = -mrSges(5,1) * t193 + mrSges(5,3) * t186 + Ifges(5,4) * t239 + Ifges(5,2) * t238 + Ifges(5,6) * qJDD(4) - pkin(4) * t283 + qJ(5) * t291 + qJD(4) * t228 + t262 * t131 + t264 * t132 - t226 * t296;
t227 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t270 - Ifges(5,2) * t267) * qJD(1);
t120 = mrSges(5,2) * t193 - mrSges(5,3) * t185 + Ifges(5,1) * t239 + Ifges(5,4) * t238 + Ifges(5,5) * qJDD(4) - qJ(5) * t137 - qJD(4) * t227 + t264 * t131 - t262 * t132 - t226 * t297;
t282 = mrSges(4,2) * t197 - mrSges(4,3) * t195 + Ifges(4,1) * qJDD(1) - pkin(7) * t129 - t267 * t118 + t270 * t120;
t281 = -mrSges(4,1) * t195 - pkin(3) * t142 - pkin(7) * t130 - t270 * t118 - t267 * t120;
t280 = -mrSges(6,1) * t167 + mrSges(6,2) * t168 - Ifges(6,5) * t212 - Ifges(6,6) * t211 - Ifges(6,3) * qJDD(4) - pkin(5) * t284 - pkin(8) * t290 - t269 * t151 - t266 * t152 - t222 * t199 - t221 * t200;
t279 = -mrSges(3,2) * t210 + pkin(2) * (-qJDD(1) * mrSges(4,2) + t285) + qJ(3) * t277 + mrSges(3,1) * t209 + Ifges(3,3) * qJDD(1) + t282;
t276 = mrSges(2,1) * t244 - mrSges(2,2) * t245 + Ifges(2,3) * qJDD(1) + pkin(1) * t124 + t279;
t275 = -mrSges(5,1) * t185 + mrSges(5,2) * t186 - Ifges(5,5) * t239 - Ifges(5,6) * t238 - Ifges(5,3) * qJDD(4) - pkin(4) * t137 - t227 * t296 - t228 * t297 + t280;
t274 = -mrSges(4,1) * t197 - pkin(3) * t129 + t275;
t122 = m(2) * t245 - t273 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t292;
t121 = m(2) * t244 + qJDD(1) * mrSges(2,1) - t273 * mrSges(2,2) + t124;
t117 = -t274 + t298 * t273 + (mrSges(3,2) - mrSges(4,3)) * t259 + t299 * qJDD(1) - mrSges(3,3) * t209 - qJ(3) * t128;
t116 = mrSges(3,3) * t210 - pkin(2) * t128 - t298 * qJDD(1) - t300 * t259 + t299 * t273 + t281;
t115 = -mrSges(2,2) * g(3) - mrSges(2,3) * t244 + Ifges(2,5) * qJDD(1) - t273 * Ifges(2,6) - qJ(2) * t124 - t263 * t116 + t265 * t117;
t114 = Ifges(2,6) * qJDD(1) + t273 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t245 + t263 * t117 + t265 * t116 - pkin(1) * (m(3) * t259 + t128) + qJ(2) * t292;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t271 * t115 - t268 * t114 - pkin(6) * (t271 * t121 + t268 * t122), t115, t117, t282, t120, t131, t152; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t268 * t115 + t271 * t114 + pkin(6) * (-t268 * t121 + t271 * t122), t114, t116, mrSges(4,3) * t259 + Ifges(4,4) * qJDD(1) - t273 * Ifges(4,5) + t274, t118, t132, t151; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t276, t276, t279, -mrSges(4,2) * t259 + t273 * Ifges(4,4) + Ifges(4,5) * qJDD(1) - t281, -t275, -t280, t278;];
m_new  = t1;
