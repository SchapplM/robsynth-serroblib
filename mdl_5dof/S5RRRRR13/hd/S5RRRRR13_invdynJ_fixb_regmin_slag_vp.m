% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR13
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
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR13_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:16
% EndTime: 2024-09-27 17:32:19
% DurationCPUTime: 2.03s
% Computational Cost: add. (4413->330), mult. (6773->456), div. (0->0), fcn. (4484->22), ass. (0->234)
t166 = sin(pkin(5));
t173 = cos(qJ(5));
t174 = cos(qJ(4));
t268 = t173 * t174;
t238 = t166 * t268;
t168 = sin(qJ(5));
t169 = sin(qJ(4));
t274 = t168 * t169;
t316 = -t166 * t274 + t238;
t196 = t168 * t174 + t169 * t173;
t315 = t196 * t166;
t159 = qJDD(1) + qJDD(2);
t149 = qJDD(3) + t159;
t161 = qJD(1) + qJD(2);
t150 = qJD(3) + t161;
t302 = qJD(4) + qJD(5);
t182 = t302 * t196;
t21 = -t149 * t238 + (t149 * t274 + t150 * t182) * t166;
t314 = t316 * t150;
t176 = cos(qJ(2));
t290 = pkin(1) * qJD(1);
t117 = pkin(2) * t161 + t176 * t290;
t170 = sin(qJ(3));
t171 = sin(qJ(2));
t258 = qJD(3) * t171;
t230 = qJD(1) * t258;
t217 = pkin(1) * t230;
t120 = t170 * t217;
t175 = cos(qJ(3));
t261 = qJD(2) * t176;
t231 = qJD(1) * t261;
t252 = qJDD(1) * t171;
t187 = (t231 + t252) * pkin(1);
t300 = pkin(1) * t176;
t147 = qJDD(1) * t300;
t247 = t171 * t290;
t88 = pkin(2) * t159 - qJD(2) * t247 + t147;
t312 = -t170 * t88 - (qJD(3) * t117 + t187) * t175 + t120;
t259 = qJD(3) * t170;
t246 = pkin(2) * t259;
t269 = t171 * t175;
t194 = -t170 * t176 - t269;
t94 = t194 * t290;
t213 = t94 + t246;
t167 = cos(pkin(5));
t275 = t167 * t174;
t270 = t170 * t171;
t193 = t175 * t176 - t270;
t95 = t193 * t290;
t308 = t169 * t95 - t275 * t94 + (-t169 * t175 - t170 * t275) * qJD(3) * pkin(2);
t255 = qJD(4) * t174;
t234 = t167 * t255;
t276 = t167 * t169;
t78 = t175 * t117 - t170 * t247;
t79 = -t117 * t170 - t175 * t247;
t307 = -pkin(3) * t234 + t174 * t78 + t276 * t79;
t145 = pkin(2) * t175 + pkin(3);
t257 = qJD(3) * t175;
t306 = -t145 * t234 + t276 * t94 + (-pkin(2) * t257 + t95) * t174;
t256 = qJD(4) * t169;
t305 = t166 * (pkin(4) * t256 + t79);
t115 = t145 * t276;
t156 = t166 * pkin(9);
t121 = pkin(2) * t170 + t156;
t266 = t174 * t121 + t115;
t165 = qJ(1) + qJ(2);
t158 = qJ(3) + t165;
t143 = sin(t158);
t144 = cos(t158);
t304 = g(1) * t144 + g(2) * t143;
t303 = g(1) * t143 - g(2) * t144;
t281 = t150 * t167;
t119 = qJD(4) + t281;
t253 = qJD(4) - t119;
t282 = t150 * t166;
t57 = pkin(9) * t282 - t79;
t301 = g(3) * t166 + t253 * t57;
t214 = pkin(10) * t282 + t57;
t297 = pkin(3) * t150;
t65 = t78 + t297;
t245 = t65 * t276;
t27 = t174 * t214 + t245;
t299 = pkin(2) * t149;
t298 = pkin(3) * t149;
t296 = pkin(4) * t174;
t295 = pkin(10) * t166;
t69 = t315 * t150;
t292 = t69 * t314;
t146 = pkin(2) + t300;
t123 = t175 * t146;
t92 = -pkin(1) * t270 + pkin(3) + t123;
t76 = t92 * t276;
t264 = pkin(1) * t269 + t170 * t146;
t86 = t156 + t264;
t291 = t174 * t86 + t76;
t289 = t150 * t79;
t288 = t173 * t27;
t237 = t171 * t257;
t56 = -t146 * t259 + (qJD(2) * t194 - t237) * pkin(1);
t287 = t56 * t150;
t160 = t166 ^ 2;
t286 = t150 ^ 2 * t160;
t285 = t149 * t166;
t284 = t149 * t167;
t283 = t149 * t174;
t280 = t150 * t174;
t279 = t160 * t174;
t278 = t166 * t169;
t277 = t166 * t174;
t271 = t169 * t174;
t164 = qJ(4) + qJ(5);
t136 = pkin(3) * t276;
t265 = pkin(9) * t277 + t136;
t153 = sin(t165);
t155 = cos(t165);
t263 = g(1) * t155 + g(2) * t153;
t162 = t169 ^ 2;
t262 = -t174 ^ 2 + t162;
t254 = qJD(5) * t168;
t151 = pkin(5) + t164;
t130 = cos(t151) / 0.2e1;
t229 = pkin(5) - t164;
t139 = cos(t229);
t251 = t139 / 0.2e1 + t130;
t32 = pkin(9) * t285 - t312;
t81 = t175 * t88;
t205 = -t117 * t259 + t81;
t178 = (-t170 * t252 + (-t170 * t261 - t237) * qJD(1)) * pkin(1) + t205;
t33 = t178 + t298;
t250 = t174 * t32 + t65 * t234 + t33 * t276;
t55 = t146 * t257 + (qJD(2) * t193 - t170 * t258) * pkin(1);
t249 = t174 * t55 + t92 * t234 + t56 * t276;
t248 = sin(t151) / 0.2e1;
t244 = t166 * (-pkin(9) - pkin(10));
t242 = t150 * t278;
t236 = t150 * t256;
t235 = t166 * t256;
t233 = -t65 - t297;
t232 = -t86 - t295;
t228 = -t121 - t295;
t54 = t65 * t275;
t26 = -t169 * t214 + t54;
t19 = pkin(4) * t119 + t26;
t191 = t236 - t283;
t8 = -t256 * t57 + t250;
t6 = -t191 * t295 + t8;
t227 = qJD(5) * t19 + t6;
t125 = pkin(4) * t235;
t226 = t213 * t166 + t125;
t225 = -t169 * t55 + t56 * t275;
t118 = qJDD(4) + t284;
t223 = t118 + t284;
t222 = t119 + t281;
t221 = qJD(1) * (-qJD(2) + t161);
t220 = qJD(2) * (-qJD(1) - t161);
t219 = t167 * t246;
t218 = t169 * t244;
t211 = g(1) * t153 - g(2) * t155 + t147;
t37 = -t169 * t78 + t275 * t79;
t134 = pkin(10) * t277;
t80 = t134 + t265;
t210 = qJD(5) * t80 + t37 - (t174 * t244 - t136) * qJD(4);
t137 = pkin(3) * t275;
t157 = t167 * pkin(4);
t71 = t137 + t157 + t218;
t209 = -qJD(4) * t218 - qJD(5) * t71 + t307;
t116 = t145 * t275;
t53 = t169 * t228 + t116 + t157;
t208 = -qJD(5) * t53 - (qJD(4) * t228 - t219) * t169 + t306;
t58 = t134 + t266;
t207 = qJD(5) * t58 - (t174 * t228 - t115) * qJD(4) - t308;
t206 = sin(t229);
t204 = t232 * t169;
t203 = -t289 + t298;
t202 = t149 * t92 + t287;
t201 = -t168 * t19 - t288;
t77 = t92 * t275;
t36 = t157 + t77 + t204;
t41 = t134 + t291;
t200 = -t168 * t41 + t173 * t36;
t199 = t168 * t36 + t173 * t41;
t198 = qJD(4) * (-t150 * t92 - t65);
t102 = t248 - t206 / 0.2e1;
t154 = cos(t164);
t197 = -t102 * t144 - t143 * t154;
t64 = t102 * t143 - t144 * t154;
t30 = t33 * t275;
t83 = -t143 * t174 - t144 * t276;
t85 = -t143 * t276 + t144 * t174;
t192 = -g(1) * t83 - g(2) * t85 + t33 * t279 + (-t169 * t32 + t30 + (-t174 * t57 - t245) * qJD(4)) * t167;
t82 = t143 * t169 - t144 * t275;
t84 = -t143 * t275 - t144 * t169;
t190 = -g(1) * t82 - g(2) * t84 - t8 * t167;
t18 = (pkin(4) * t191 - t33) * t166;
t5 = pkin(4) * t118 + t30 + (-pkin(10) * t285 - t32) * t169 - t27 * qJD(4);
t184 = qJD(5) * t201 - t168 * t6 + t173 * t5;
t43 = (-pkin(4) * t280 - t65) * t166;
t50 = t182 * t166;
t188 = -g(1) * t197 + g(2) * t64 + t184 * t167 - t18 * t316 + t43 * t50;
t22 = t27 * t254;
t49 = t302 * (t268 - t274) * t166;
t152 = sin(t164);
t59 = t143 * t152 - t144 * t251;
t62 = -t143 * t251 - t144 * t152;
t186 = -g(1) * t59 - g(2) * t62 - (t168 * t5 + t227 * t173 - t22) * t167 + t18 * t315 + t43 * t49;
t20 = t149 * t315 + t314 * t302;
t183 = (-pkin(2) * t150 - t117) * qJD(3) - t187;
t114 = qJD(5) + t119;
t181 = -g(1) * t64 - g(2) * t197 - g(3) * (t130 - t139 / 0.2e1) - t43 * t314 + t22 + (-t27 * t114 - t5) * t168;
t180 = t304 + t312;
t179 = -g(1) * t62 + g(2) * t59 - g(3) * (t248 + t206 / 0.2e1) - t43 * t69 + t184;
t177 = cos(qJ(1));
t172 = sin(qJ(1));
t113 = (-pkin(3) - t296) * t166;
t111 = qJDD(5) + t118;
t99 = t118 * t167;
t96 = t111 * t167;
t93 = (-t145 - t296) * t166;
t70 = (t149 * t162 + 0.2e1 * t174 * t236) * t160;
t66 = (-t92 - t296) * t166;
t46 = 0.2e1 * (-qJD(4) * t150 * t262 + t149 * t271) * t160;
t42 = -t166 * t56 + t125;
t35 = (t169 * t223 + t222 * t255) * t166;
t34 = (t174 * t223 - t222 * t256) * t166;
t25 = -t314 ^ 2 + t69 ^ 2;
t15 = t114 * t69 - t21;
t14 = -t114 * t314 + t20;
t13 = (t174 * t232 - t76) * qJD(4) + t225;
t12 = qJD(4) * t204 + t249;
t11 = t20 * t315 + t49 * t69;
t10 = t111 * t316 - t114 * t50 - t167 * t21;
t9 = t111 * t315 + t114 * t49 + t167 * t20;
t3 = t20 * t316 - t21 * t315 + t314 * t49 - t50 * t69;
t1 = [qJDD(1), g(1) * t172 - g(2) * t177, g(1) * t177 + g(2) * t172, t159, (t159 * t176 + t171 * t220) * pkin(1) + t211, ((-qJDD(1) - t159) * t171 + t176 * t220) * pkin(1) + t263, t149, t123 * t149 + t287 + (-t175 * t230 + (-t231 + (-qJDD(1) - t149) * t171) * t170) * pkin(1) + t205 + t303, -t149 * t264 - t55 * t150 + t180, t70, t46, t35, t34, t99, (-qJD(4) * t291 + t225) * t119 + (-t169 * t86 + t77) * t118 + (t169 * t198 + t174 * t202) * t160 + t192, -(-t256 * t86 + t249) * t119 - t291 * t118 + (t174 * t198 + (-t202 - t33) * t169) * t160 + t190, t11, t3, t9, t10, t96, (-qJD(5) * t199 - t12 * t168 + t13 * t173) * t114 + t200 * t111 - t42 * t314 + t66 * t21 + t188, -(qJD(5) * t200 + t12 * t173 + t13 * t168) * t114 - t199 * t111 + t42 * t69 + t66 * t20 + t186; 0, 0, 0, t159, pkin(1) * t171 * t221 + t211, (t176 * t221 - t252) * pkin(1) + t263, t149, -t150 * t94 + t81 + (-t217 + t299) * t175 + t183 * t170 + t303, t150 * t95 + t120 + (-t88 - t299) * t170 + t183 * t175 + t304, t70, t46, t35, t34, t99, (-t121 * t169 + t116) * t118 + (-t266 * qJD(4) + t308) * t119 + (-t65 * t256 + t145 * t283 + (-t145 * t256 - t174 * t213) * t150) * t160 + t192, -t266 * t118 + ((qJD(4) * t121 + t219) * t169 + t306) * t119 + ((-t145 * t150 - t65) * t255 + (-t145 * t149 + t150 * t213 - t33) * t169) * t160 + t190, t11, t3, t9, t10, t96, (-t168 * t58 + t173 * t53) * t111 + t93 * t21 - t226 * t314 + (t168 * t208 - t173 * t207) * t114 + t188, -(t168 * t53 + t173 * t58) * t111 + t93 * t20 + t226 * t69 + (t168 * t207 + t173 * t208) * t114 + t186; 0, 0, 0, 0, 0, 0, t149, t178 - t289 + t303, t150 * t78 + t180, t70, t46, t35, t34, t99, (-pkin(9) * t278 + t137) * t118 - t37 * t119 + t203 * t279 + (t160 * t169 * t233 - t119 * t265) * qJD(4) + t192, -t265 * t118 + (pkin(9) * t235 + t307) * t119 + (t233 * t255 + (-t203 - t33) * t169) * t160 + t190, t11, t3, t9, t10, t96, (-t168 * t80 + t173 * t71) * t111 + t113 * t21 - t314 * t305 + (t168 * t209 - t173 * t210) * t114 + t188, -(t168 * t71 + t173 * t80) * t111 + t113 * t20 + t69 * t305 + (t168 * t210 + t173 * t209) * t114 + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271 * t286, t262 * t286, (t149 * t169 + t253 * t280) * t166, (-t150 * t169 * t253 + t283) * t166, t118, -g(1) * t84 + g(2) * t82 + t30 - t301 * t174 + (-t32 + (t150 * t160 - t167 * t253) * t65) * t169, t65 * t150 * t279 + g(1) * t85 - g(2) * t83 + t54 * t119 + t301 * t169 - t250, -t292, t25, t14, t15, t111, -(-t168 * t26 - t288) * t114 + (t173 * t111 - t114 * t254 + t242 * t314) * pkin(4) + t179, (t26 * t114 - t227) * t173 + (-qJD(5) * t173 * t114 - t168 * t111 - t242 * t69) * pkin(4) + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t292, t25, t14, t15, t111, -t114 * t201 + t179, (-t6 + (-qJD(5) + t114) * t19) * t173 + t181;];
tau_reg = t1;
