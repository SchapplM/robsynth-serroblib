% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPR5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:34
% EndTime: 2019-12-31 17:03:37
% DurationCPUTime: 1.82s
% Computational Cost: add. (3441->220), mult. (4776->300), div. (0->0), fcn. (2684->6), ass. (0->150)
t226 = (qJD(1) + qJD(2));
t224 = t226 ^ 2;
t233 = cos(qJ(2));
t225 = qJDD(1) + qJDD(2);
t230 = sin(qJ(2));
t256 = t225 * t230;
t196 = t224 * t233 + t256;
t254 = t225 * t233;
t198 = t224 * t230 - t254;
t231 = sin(qJ(1));
t234 = cos(qJ(1));
t156 = t196 * t234 - t198 * t231;
t247 = pkin(5) * t196 - t233 * g(3);
t285 = -pkin(5) * t198 + g(3) * t230;
t287 = pkin(4) * t156 + t231 * t285 + t234 * t247;
t244 = t196 * t231 + t234 * t198;
t277 = pkin(4) * t244 + t231 * t247 - t234 * t285;
t213 = g(1) * t231 - t234 * g(2);
t203 = qJDD(1) * pkin(1) + t213;
t214 = g(1) * t234 + g(2) * t231;
t237 = qJD(1) ^ 2;
t204 = -pkin(1) * t237 - t214;
t165 = -t233 * t203 + t204 * t230;
t166 = t203 * t230 + t204 * t233;
t245 = t165 * t230 + t233 * t166;
t130 = t165 * t233 - t166 * t230;
t272 = t130 * t234;
t290 = -t231 * t245 + t272;
t273 = t130 * t231;
t107 = t234 * t245 + t273;
t289 = pkin(1) * t196;
t288 = pkin(1) * t198;
t238 = (2 * qJD(3) * t226) + t166;
t274 = qJ(3) * t225;
t143 = -t224 * pkin(2) + t238 + t274;
t217 = t225 * pkin(2);
t239 = qJDD(3) + t165 - t217;
t148 = -t224 * qJ(3) + t239;
t122 = t143 * t230 - t148 * t233;
t246 = t233 * t143 + t148 * t230;
t101 = -t122 * t231 + t234 * t246;
t100 = t122 * t234 + t231 * t246;
t276 = pkin(2) + pkin(6);
t136 = -t224 * pkin(6) + t143;
t229 = sin(qJ(4));
t271 = t136 * t229;
t232 = cos(qJ(4));
t270 = t136 * t232;
t227 = t229 ^ 2;
t228 = t232 ^ 2;
t253 = t227 + t228;
t195 = t253 * t225;
t265 = t195 * t230;
t264 = t195 * t233;
t251 = t224 * t229 * t232;
t205 = qJDD(4) + t251;
t263 = t205 * t229;
t262 = t205 * t232;
t206 = qJDD(4) - t251;
t261 = t206 * t229;
t260 = t206 * t232;
t259 = t224 * t227;
t258 = t224 * t228;
t257 = t225 * t229;
t255 = t225 * t232;
t252 = qJD(4) * t226;
t250 = t229 * t252;
t249 = t232 * t252;
t178 = -t213 * t231 - t234 * t214;
t242 = t230 * t251;
t241 = t233 * t251;
t208 = qJDD(1) * t234 - t231 * t237;
t240 = -pkin(4) * t208 - g(3) * t231;
t142 = -pkin(6) * t225 + t148;
t135 = -g(3) * t232 + t142 * t229;
t134 = g(3) * t229 + t142 * t232;
t111 = t134 * t232 + t135 * t229;
t112 = -t134 * t229 + t135 * t232;
t177 = t213 * t234 - t214 * t231;
t236 = qJD(4) ^ 2;
t235 = pkin(1) * g(3);
t212 = -t236 - t258;
t211 = t236 - t258;
t210 = -t236 - t259;
t209 = -t236 + t259;
t207 = qJDD(1) * t231 + t234 * t237;
t201 = (-t227 + t228) * t224;
t200 = t253 * t224;
t193 = -0.2e1 * t250 + t255;
t192 = -t250 + t255;
t191 = -t249 - t257;
t190 = 0.2e1 * t249 + t257;
t188 = -pkin(4) * t207 + g(3) * t234;
t187 = t253 * t252;
t176 = qJDD(4) * t233 - t187 * t230;
t175 = qJDD(4) * t230 + t187 * t233;
t174 = -t192 * t229 - t228 * t252;
t173 = -t191 * t232 - t227 * t252;
t172 = -t212 * t229 - t262;
t171 = t210 * t232 - t261;
t170 = t212 * t232 - t263;
t169 = -t211 * t232 - t261;
t168 = t210 * t229 + t260;
t167 = -t209 * t229 - t262;
t159 = -t200 * t233 - t265;
t154 = -t200 * t230 + t264;
t153 = t190 * t229 - t193 * t232;
t152 = -t169 * t230 + t232 * t254;
t151 = -t167 * t230 - t229 * t254;
t150 = t169 * t233 + t230 * t255;
t149 = t167 * t233 - t229 * t256;
t147 = -t173 * t230 - t241;
t146 = -t174 * t230 + t241;
t145 = t173 * t233 - t242;
t144 = t174 * t233 + t242;
t140 = t170 * t230 + t193 * t233;
t139 = t168 * t230 + t190 * t233;
t138 = -t170 * t233 + t193 * t230;
t137 = -t168 * t233 + t190 * t230;
t133 = -t153 * t230 + t201 * t233;
t132 = t153 * t233 + t201 * t230;
t127 = pkin(5) * t245 + t235;
t126 = -t154 * t231 + t159 * t234;
t125 = t154 * t234 + t159 * t231;
t120 = -t138 * t231 + t140 * t234;
t119 = -t137 * t231 + t139 * t234;
t118 = t138 * t234 + t140 * t231;
t117 = t137 * t234 + t139 * t231;
t116 = -pkin(5) * t122 + (-pkin(2) * t230 + qJ(3) * t233) * g(3);
t115 = pkin(5) * t246 + t235 + (pkin(2) * t233 + qJ(3) * t230) * g(3);
t114 = pkin(3) * t170 - qJ(3) * t172 - t135;
t113 = pkin(3) * t168 - qJ(3) * t171 + t134;
t110 = pkin(3) * t190 - t276 * t171 + t270;
t109 = pkin(3) * t193 - t276 * t172 - t271;
t108 = -pkin(3) * t200 - t112;
t105 = t111 * t230 + t136 * t233;
t104 = -t111 * t233 + t136 * t230;
t103 = -pkin(3) * t264 - pkin(5) * t154 - t108 * t230;
t102 = -pkin(3) * t265 + pkin(5) * t159 + t108 * t233;
t99 = pkin(3) * t111 - qJ(3) * t112;
t98 = pkin(3) * t136 - t276 * t112;
t97 = -pkin(5) * t138 - t109 * t230 + t114 * t233;
t96 = -pkin(5) * t137 - t110 * t230 + t113 * t233;
t95 = -pkin(1) * t172 + pkin(5) * t140 + t109 * t233 + t114 * t230;
t94 = -pkin(1) * t171 + pkin(5) * t139 + t110 * t233 + t113 * t230;
t93 = -t104 * t231 + t105 * t234;
t92 = t104 * t234 + t105 * t231;
t91 = -pkin(5) * t104 - t230 * t98 + t233 * t99;
t90 = -pkin(1) * t112 + pkin(5) * t105 + t230 * t99 + t233 * t98;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t207, -t208, 0, t178, 0, 0, 0, 0, 0, 0, -t156, t244, 0, t107, 0, 0, 0, 0, 0, 0, 0, t156, -t244, t101, 0, 0, 0, 0, 0, 0, t119, t120, t126, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t208, -t207, 0, t177, 0, 0, 0, 0, 0, 0, -t244, -t156, 0, -t290, 0, 0, 0, 0, 0, 0, 0, t244, t156, t100, 0, 0, 0, 0, 0, 0, t117, t118, t125, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t171, t172, 0, t112; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t208, 0, -t207, 0, t240, -t188, -t177, -pkin(4) * t177, 0, 0, -t244, 0, -t156, 0, t277, t287, t290, pkin(4) * t290 + pkin(5) * t272 - t127 * t231, 0, t244, t156, 0, 0, 0, -t100, -t277, -t287, -pkin(4) * t100 - t115 * t231 + t116 * t234, -t144 * t231 + t146 * t234, -t132 * t231 + t133 * t234, -t150 * t231 + t152 * t234, -t145 * t231 + t147 * t234, -t149 * t231 + t151 * t234, -t175 * t231 + t176 * t234, -pkin(4) * t117 - t231 * t94 + t234 * t96, -pkin(4) * t118 - t231 * t95 + t234 * t97, -pkin(4) * t125 - t102 * t231 + t103 * t234, -pkin(4) * t92 - t231 * t90 + t234 * t91; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t207, 0, t208, 0, t188, t240, t178, pkin(4) * t178, 0, 0, t156, 0, -t244, 0, -t287, t277, t107, pkin(4) * t107 + pkin(5) * t273 + t127 * t234, 0, -t156, t244, 0, 0, 0, t101, t287, -t277, pkin(4) * t101 + t115 * t234 + t116 * t231, t144 * t234 + t146 * t231, t132 * t234 + t133 * t231, t150 * t234 + t152 * t231, t145 * t234 + t147 * t231, t149 * t234 + t151 * t231, t175 * t234 + t176 * t231, pkin(4) * t119 + t231 * t96 + t234 * t94, pkin(4) * t120 + t231 * t97 + t234 * t95, pkin(4) * t126 + t102 * t234 + t103 * t231, pkin(4) * t93 + t231 * t91 + t234 * t90; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t213, t214, 0, 0, 0, 0, 0, 0, 0, t225, -t165 - t288, -t166 - t289, 0, -pkin(1) * t130, t225, 0, 0, 0, 0, 0, 0, -t217 + t239 + t288, t238 + 0.2e1 * t274 + t289, pkin(1) * t122 - pkin(2) * t148 + qJ(3) * t143, (t192 - t250) * t232, -t190 * t232 - t193 * t229, -t211 * t229 + t260, (-t191 + t249) * t229, t209 * t232 - t263, 0, pkin(1) * t137 + qJ(3) * t190 - t276 * t168 + t271, pkin(1) * t138 + qJ(3) * t193 - t276 * t170 + t270, pkin(1) * t154 - qJ(3) * t200 + t276 * t195 - t111, pkin(1) * t104 + qJ(3) * t136 - t276 * t111;];
tauB_reg = t1;
