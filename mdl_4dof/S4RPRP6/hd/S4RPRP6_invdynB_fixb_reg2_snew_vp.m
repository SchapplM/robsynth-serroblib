% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRP6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynB_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:15
% EndTime: 2019-12-31 16:46:17
% DurationCPUTime: 1.18s
% Computational Cost: add. (1786->212), mult. (3771->260), div. (0->0), fcn. (1793->4), ass. (0->162)
t258 = (qJD(1) * qJD(4));
t288 = -2 * t258;
t284 = -pkin(5) - pkin(1);
t226 = qJD(3) ^ 2;
t222 = sin(qJ(3));
t219 = t222 ^ 2;
t227 = qJD(1) ^ 2;
t267 = t219 * t227;
t206 = -t226 - t267;
t224 = cos(qJ(3));
t264 = t222 * t227;
t252 = t224 * t264;
t202 = qJDD(3) - t252;
t268 = t202 * t224;
t159 = t206 * t222 + t268;
t269 = t202 * t222;
t163 = t206 * t224 - t269;
t287 = pkin(2) * t159 - qJ(2) * t163;
t220 = t224 ^ 2;
t266 = t220 * t227;
t208 = -t226 - t266;
t201 = qJDD(3) + t252;
t271 = t201 * t222;
t161 = t208 * t224 - t271;
t270 = t201 * t224;
t166 = -t208 * t222 - t270;
t286 = pkin(2) * t161 - qJ(2) * t166;
t259 = qJD(1) * qJD(3);
t248 = t224 * t259;
t256 = t222 * qJDD(1);
t192 = -t248 - t256;
t260 = qJD(1) * t224;
t200 = qJD(3) * pkin(3) - qJ(4) * t260;
t223 = sin(qJ(1));
t225 = cos(qJ(1));
t203 = t223 * g(1) - t225 * g(2);
t239 = qJDD(2) - t203;
t231 = -t227 * qJ(2) + t239;
t175 = t284 * qJDD(1) + t231;
t262 = t224 * g(3) - t222 * t175;
t232 = t192 * qJ(4) - qJD(3) * t200 + t222 * t288 - t262;
t285 = -t200 * t260 + (qJ(4) * t219 - t284) * t227 - qJDD(4);
t261 = t219 + t220;
t198 = t261 * t227;
t283 = pkin(2) * t198;
t191 = 0.2e1 * t248 + t256;
t130 = -t159 * t225 + t191 * t223;
t282 = pkin(4) * t130;
t249 = t222 * t259;
t254 = t224 * qJDD(1);
t194 = -0.2e1 * t249 + t254;
t131 = -t161 * t225 + t194 * t223;
t281 = pkin(4) * t131;
t195 = t261 * qJDD(1);
t272 = t195 * t225;
t148 = -t198 * t223 + t272;
t280 = pkin(4) * t148;
t279 = t192 * pkin(3);
t276 = qJDD(1) * pkin(1);
t193 = -t249 + t254;
t238 = -t193 - t249;
t146 = t222 * g(3) + t224 * t175;
t251 = -qJDD(3) * pkin(3) - t146;
t229 = t238 * qJ(4) + t224 * t288 - t251;
t119 = -pkin(3) * t252 + t229;
t275 = t119 * t222;
t274 = t119 * t224;
t273 = t195 * t223;
t204 = t225 * g(1) + t223 * g(2);
t218 = qJDD(1) * qJ(2);
t237 = t204 - t218;
t257 = qJD(2) * qJD(1);
t230 = t237 - 0.2e1 * t257;
t170 = -t284 * t227 + t230;
t265 = t222 * t170;
t263 = t224 * t170;
t255 = t223 * qJDD(1);
t253 = t225 * qJDD(1);
t214 = 0.2e1 * t257;
t233 = t214 - t237;
t176 = -t227 * pkin(1) + t233;
t178 = -t231 + t276;
t132 = t225 * t176 - t178 * t223;
t155 = -t203 * t223 - t225 * t204;
t247 = t223 * t252;
t246 = t225 * t252;
t244 = qJ(2) * t191 + t284 * t159;
t243 = qJ(2) * t194 + t284 * t161;
t242 = -qJ(2) * t198 - t284 * t195;
t196 = -t223 * t227 + t253;
t241 = pkin(4) * t196 + t223 * g(3);
t197 = t225 * t227 + t255;
t240 = -pkin(4) * t197 + t225 * g(3);
t122 = t146 * t224 - t222 * t262;
t123 = -t146 * t222 - t224 * t262;
t129 = t176 * t223 + t178 * t225;
t154 = t203 * t225 - t204 * t223;
t235 = pkin(2) * t191 + t284 * t163;
t234 = pkin(2) * t194 + t284 * t166;
t228 = t230 + t285;
t207 = t226 - t266;
t205 = -t226 + t267;
t199 = (-t219 + t220) * t227;
t187 = t261 * t259;
t177 = -pkin(2) * t195 - pkin(3) * t254;
t174 = t225 * qJDD(3) - t223 * t187;
t173 = t223 * qJDD(3) + t225 * t187;
t172 = -t222 * t193 - t220 * t259;
t171 = -t224 * t192 - t219 * t259;
t165 = -t207 * t222 + t268;
t164 = (t193 - t249) * t224;
t162 = t205 * t224 - t271;
t160 = -t207 * t224 - t269;
t158 = -t205 * t222 - t270;
t157 = (-t192 + t248) * t222;
t156 = -pkin(3) * t194 - qJ(4) * t201;
t149 = -t198 * t225 - t273;
t145 = pkin(4) * t149;
t144 = -t191 * t224 - t194 * t222;
t143 = t191 * t222 - t194 * t224;
t142 = -t171 * t223 - t246;
t141 = -t172 * t223 + t246;
t140 = t171 * t225 - t247;
t139 = t172 * t225 + t247;
t138 = -t160 * t223 + t224 * t253;
t137 = -t158 * t223 - t222 * t253;
t136 = t160 * t225 + t223 * t254;
t135 = t158 * t225 - t222 * t255;
t134 = t161 * t223 + t194 * t225;
t133 = t159 * t223 + t191 * t225;
t128 = pkin(4) * t134;
t127 = pkin(4) * t133;
t126 = -t143 * t223 + t199 * t225;
t125 = t143 * t225 + t199 * t223;
t124 = t228 + t279;
t121 = -pkin(3) * t267 + t232;
t120 = -qJ(4) * t208 + t233 - t279 - t285;
t118 = (pkin(3) * t264 + (2 * t258)) * t224 + (-t238 + t254) * qJ(4) + t251;
t117 = -t123 - t283;
t116 = qJ(4) * t206 + (-t191 + t192) * pkin(3) + t228;
t115 = t262 + t286;
t114 = t146 + t287;
t113 = -qJ(4) * t256 + (t198 - t267) * pkin(3) + t232;
t112 = t235 - t263;
t111 = t234 + t265;
t110 = t122 * t223 - t170 * t225;
t109 = -t122 * t225 - t170 * t223;
t108 = pkin(3) * t124 + qJ(4) * t121;
t107 = (t208 + t267) * pkin(3) - t232 + t286;
t106 = (t202 - t252) * pkin(3) + t229 + t287;
t105 = pkin(2) * t122 - qJ(2) * t123;
t104 = t121 * t224 - t275;
t103 = t121 * t222 + t274;
t102 = -pkin(2) * t170 + t284 * t123;
t101 = qJ(4) * t269 - t224 * t116 + t235;
t100 = -t222 * t120 - t224 * t156 + t234;
t99 = -t113 * t224 - t118 * t222 - t283;
t98 = t103 * t223 - t124 * t225;
t97 = -t103 * t225 - t124 * t223;
t96 = pkin(2) * t103 + pkin(3) * t119 - qJ(2) * t104;
t95 = -pkin(2) * t124 + qJ(4) * t275 + t284 * t104 - t224 * t108;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t197, -t196, 0, t155, 0, 0, 0, 0, 0, 0, 0, t197, t196, t132, 0, 0, 0, 0, 0, 0, t133, t134, t149, t110, 0, 0, 0, 0, 0, 0, t133, t134, t149, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t196, -t197, 0, t154, 0, 0, 0, 0, 0, 0, 0, -t196, t197, t129, 0, 0, 0, 0, 0, 0, t130, t131, t148, t109, 0, 0, 0, 0, 0, 0, t130, t131, t148, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t163, t166, 0, t123, 0, 0, 0, 0, 0, 0, t163, t166, 0, t104; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t196, 0, -t197, 0, -t241, -t240, -t154, -pkin(4) * t154, 0, -t196, t197, 0, 0, 0, -t129, t241, t240, -pkin(4) * t129 + (-pkin(1) * t223 + qJ(2) * t225) * g(3), t141, t126, t138, t142, t137, t174, -t112 * t223 + t114 * t225 - t282, -t111 * t223 + t115 * t225 - t281, -pkin(2) * t272 - t117 * t223 - t280, -pkin(4) * t109 - t102 * t223 + t105 * t225, t141, t126, t138, t142, t137, t174, -t101 * t223 + t106 * t225 - t282, -t100 * t223 + t107 * t225 - t281, t177 * t225 - t223 * t99 - t280, -pkin(4) * t97 - t223 * t95 + t225 * t96; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t197, 0, t196, 0, t240, -t241, t155, pkin(4) * t155, 0, -t197, -t196, 0, 0, 0, t132, -t240, t241, pkin(4) * t132 + (pkin(1) * t225 + qJ(2) * t223) * g(3), t139, t125, t136, t140, t135, t173, t112 * t225 + t114 * t223 + t127, t111 * t225 + t115 * t223 + t128, -pkin(2) * t273 + t117 * t225 + t145, pkin(4) * t110 + t102 * t225 + t105 * t223, t139, t125, t136, t140, t135, t173, t101 * t225 + t106 * t223 + t127, t100 * t225 + t107 * t223 + t128, t177 * t223 + t225 * t99 + t145, pkin(4) * t98 + t223 * t96 + t225 * t95; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t203, t204, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t239 - 0.2e1 * t276, -t204 + t214 + 0.2e1 * t218, pkin(1) * t178 + qJ(2) * t176, t164, t144, t165, t157, t162, 0, t244 - t265, t243 - t263, -t122 + t242, -qJ(2) * t170 + t284 * t122, t164, t144, t165, t157, t162, 0, -qJ(4) * t268 - t116 * t222 + t244, t120 * t224 - t156 * t222 + t243, -t113 * t222 + t118 * t224 + t242, -qJ(2) * t124 - qJ(4) * t274 + t284 * t103 - t222 * t108;];
tauB_reg = t1;
