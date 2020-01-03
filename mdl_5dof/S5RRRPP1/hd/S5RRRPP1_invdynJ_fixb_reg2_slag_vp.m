% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:50:03
% EndTime: 2019-12-31 20:50:06
% DurationCPUTime: 2.19s
% Computational Cost: add. (3353->341), mult. (5117->390), div. (0->0), fcn. (3270->12), ass. (0->202)
t156 = sin(pkin(8));
t158 = sin(qJ(3));
t161 = cos(qJ(3));
t157 = -qJ(4) - pkin(7);
t216 = qJD(3) * t157;
t181 = -t158 * qJD(4) + t161 * t216;
t251 = cos(pkin(8));
t213 = t251 * t161;
t243 = t156 * t158;
t185 = t213 - t243;
t162 = cos(qJ(2));
t231 = qJD(1) * t162;
t223 = pkin(1) * t231;
t143 = t161 * qJD(4);
t93 = t158 * t216 + t143;
t255 = t156 * t181 - t185 * t223 + t251 * t93;
t153 = t158 ^ 2;
t154 = t161 ^ 2;
t232 = t153 + t154;
t149 = qJDD(1) + qJDD(2);
t159 = sin(qJ(2));
t227 = qJDD(1) * t159;
t229 = qJD(2) * t162;
t89 = t149 * pkin(7) + (qJD(1) * t229 + t227) * pkin(1);
t215 = t232 * t89;
t155 = qJ(1) + qJ(2);
t145 = cos(t155);
t132 = g(2) * t145;
t254 = pkin(1) * qJD(1);
t224 = t159 * t254;
t258 = t162 * pkin(1);
t234 = -qJD(2) * t224 + qJDD(1) * t258;
t262 = t149 * pkin(2);
t88 = -t234 - t262;
t282 = t88 + t132;
t150 = qJD(1) + qJD(2);
t272 = t232 * t162;
t281 = t150 * t272;
t103 = t150 * pkin(7) + t224;
t261 = t150 * pkin(2);
t104 = -t223 - t261;
t280 = t103 * t272 + t104 * t159;
t214 = t251 * t158;
t100 = t156 * t161 + t214;
t273 = t100 * t150;
t259 = t161 * pkin(3);
t136 = pkin(2) + t259;
t83 = -t136 * t150 + qJD(4) - t223;
t198 = t150 * t213;
t84 = t150 * t243 - t198;
t36 = t84 * pkin(4) - qJ(5) * t273 + t83;
t242 = t158 * t149;
t193 = -t149 * t213 + t156 * t242;
t94 = t100 * qJD(3);
t49 = t150 * t94 + t193;
t228 = qJD(3) * t158;
t218 = t156 * t228;
t173 = t100 * t149 - t150 * t218;
t196 = qJD(3) * t213;
t50 = t150 * t196 + t173;
t219 = t150 * t228;
t54 = pkin(3) * t219 - t136 * t149 + qJDD(4) - t234;
t167 = t49 * pkin(4) - t50 * qJ(5) + t54;
t9 = -qJD(5) * t273 + t167;
t279 = -t185 * t9 + t36 * t94;
t95 = t196 - t218;
t278 = -t9 * t100 - t36 * t95;
t277 = -t185 * t54 + t83 * t94;
t276 = t54 * t100 + t83 * t95;
t266 = t273 ^ 2;
t79 = t84 ^ 2;
t275 = -t79 - t266;
t274 = -t79 + t266;
t144 = sin(t155);
t271 = g(1) * t145 + g(2) * t144;
t133 = g(1) * t144;
t270 = t132 - t133;
t186 = qJ(4) * t149 + qJD(4) * t150 + t89;
t208 = qJ(4) * t150 + t103;
t188 = qJD(3) * t208;
t30 = qJDD(3) * pkin(3) - t186 * t158 - t161 * t188;
t35 = -t158 * t188 + t186 * t161;
t13 = t156 * t30 + t251 * t35;
t152 = qJDD(3) * qJ(5);
t10 = qJD(3) * qJD(5) + t13 + t152;
t12 = -t156 * t35 + t251 * t30;
t250 = qJDD(3) * pkin(4);
t11 = qJDD(5) - t250 - t12;
t73 = t208 * t161;
t253 = t156 * t73;
t72 = t208 * t158;
t71 = qJD(3) * pkin(3) - t72;
t41 = t251 * t71 - t253;
t39 = -qJD(3) * pkin(4) + qJD(5) - t41;
t67 = t251 * t73;
t42 = t156 * t71 + t67;
t40 = qJD(3) * qJ(5) + t42;
t269 = t10 * t185 + t11 * t100 + t39 * t95 - t40 * t94;
t268 = -t12 * t100 + t13 * t185 - t41 * t95 - t42 * t94;
t151 = qJ(3) + pkin(8);
t141 = sin(t151);
t142 = cos(t151);
t267 = g(3) * t141 + t142 * t271 - t13;
t265 = pkin(3) * t158;
t160 = sin(qJ(1));
t264 = g(1) * t160;
t263 = g(3) * t161;
t260 = t160 * pkin(1);
t257 = t273 * t84;
t256 = -t100 * t223 + t156 * t93 - t251 * t181;
t252 = t104 * t228 + t161 * t133;
t248 = t141 * t144;
t247 = t141 * t145;
t246 = t142 * t145;
t245 = t145 * t157;
t244 = t150 * t158;
t241 = t161 * t149;
t44 = -t156 * t72 + t67;
t240 = t44 * qJD(3);
t135 = t159 * pkin(1) + pkin(7);
t239 = -qJ(4) - t135;
t45 = -t251 * t72 - t253;
t238 = qJD(5) - t45;
t237 = g(1) * t248 - g(2) * t247;
t236 = t145 * pkin(2) + t144 * pkin(7);
t233 = t153 - t154;
t230 = qJD(2) * t159;
t226 = t104 * qJD(3) * t161 + t282 * t158;
t225 = pkin(1) * t229;
t139 = pkin(3) * t228;
t148 = t150 ^ 2;
t222 = t158 * t148 * t161;
t113 = t145 * t136;
t221 = pkin(4) * t246 + qJ(5) * t247 + t113;
t220 = t150 * t230;
t212 = t239 * t158;
t210 = t232 * t149;
t207 = -t144 * t157 + t113;
t206 = qJD(3) * t239;
t205 = -t271 + t215;
t204 = t150 * t224;
t203 = -t234 + t270;
t201 = t161 * t219;
t43 = t94 * pkin(4) - t95 * qJ(5) - t100 * qJD(5) + t139;
t200 = -t43 + t224;
t199 = -g(2) * t246 + t142 * t133;
t197 = g(1) * (-t144 * pkin(2) + t145 * pkin(7));
t163 = cos(qJ(1));
t194 = -g(2) * t163 + t264;
t192 = -t185 * t49 + t84 * t94;
t146 = t161 * qJ(4);
t119 = t161 * pkin(7) + t146;
t70 = t251 * t119 + t157 * t243;
t191 = t70 * qJDD(3) + t237;
t190 = -t142 * pkin(4) - t141 * qJ(5);
t189 = -t144 * t136 - t245;
t62 = qJD(3) * t94 - qJDD(3) * t185;
t187 = g(1) * t247 + g(2) * t248 - g(3) * t142 + t12;
t69 = t156 * t119 - t157 * t214;
t184 = -t69 * qJDD(3) + t199;
t169 = (-qJD(4) - t225) * t158 + t161 * t206;
t64 = t158 * t206 + t161 * t225 + t143;
t32 = t156 * t169 + t251 * t64;
t98 = t161 * t135 + t146;
t58 = t156 * t212 + t251 * t98;
t182 = -t32 * qJD(3) - t58 * qJDD(3) - t237;
t59 = -pkin(4) * t185 - t100 * qJ(5) - t136;
t180 = -t104 * t150 + t271 - t89;
t31 = t156 * t64 - t251 * t169;
t57 = t156 * t98 - t251 * t212;
t178 = t273 * t31 - t32 * t84 - t58 * t49 + t57 * t50 - t271;
t164 = qJD(3) ^ 2;
t177 = pkin(7) * t164 - t204 - t262;
t1 = t100 * t49 - t185 * t50 + t273 * t94 + t95 * t84;
t137 = -pkin(2) - t258;
t176 = pkin(1) * t220 + t135 * t164 + t137 * t149;
t175 = -t31 * qJD(3) - t57 * qJDD(3) + t199;
t174 = -t273 * t36 - qJDD(5) + t187;
t172 = -pkin(7) * qJDD(3) + (t223 - t261) * qJD(3);
t171 = -qJDD(3) * t135 + (t137 * t150 - t225) * qJD(3);
t170 = (-g(1) * (-t136 + t190) + g(2) * t157) * t144;
t168 = -t255 * t84 + t256 * t273 - t70 * t49 + t69 * t50 - t271;
t166 = 0.2e1 * t273 * qJD(3) + t193;
t147 = t163 * pkin(1);
t140 = pkin(1) * t230;
t127 = -t251 * pkin(3) - pkin(4);
t125 = t156 * pkin(3) + qJ(5);
t117 = -t136 - t258;
t115 = qJDD(3) * t161 - t164 * t158;
t114 = qJDD(3) * t158 + t164 * t161;
t102 = t140 + t139;
t91 = t154 * t149 - 0.2e1 * t201;
t90 = t153 * t149 + 0.2e1 * t201;
t66 = -0.2e1 * t233 * t150 * qJD(3) + 0.2e1 * t158 * t241;
t61 = t95 * qJD(3) + t100 * qJDD(3);
t53 = t59 - t258;
t46 = pkin(3) * t244 + pkin(4) * t273 + t84 * qJ(5);
t37 = t140 + t43;
t34 = (t198 + t84) * qJD(3) + t173;
t33 = (t198 - t84) * qJD(3) + t173;
t16 = t50 * t100 + t273 * t95;
t2 = [0, 0, 0, 0, 0, qJDD(1), t194, g(1) * t163 + g(2) * t160, 0, 0, 0, 0, 0, 0, 0, t149, (t149 * t162 - t220) * pkin(1) - t203, ((-qJDD(1) - t149) * t159 + (-qJD(1) - t150) * t229) * pkin(1) + t271, 0, (t194 + (t159 ^ 2 + t162 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t90, t66, t114, t91, t115, 0, t171 * t158 + (-t176 - t282) * t161 + t252, t171 * t161 + (t176 - t133) * t158 + t226, pkin(1) * qJD(2) * t281 + t135 * t210 + t205, t88 * t137 - t197 - g(2) * (t147 + t236) + t135 * t215 + (t280 * qJD(2) + t264) * pkin(1), t16, -t1, t61, t192, -t62, 0, t102 * t84 + t117 * t49 + t175 + t277, t102 * t273 + t117 * t50 + t182 + t276, t178 + t268, t13 * t58 + t42 * t32 - t12 * t57 - t41 * t31 + t54 * t117 + t83 * t102 - g(1) * (t189 - t260) - g(2) * (t147 + t207), t16, t61, t1, 0, t62, t192, t37 * t84 + t53 * t49 + t175 + t279, t178 + t269, -t273 * t37 - t53 * t50 - t182 + t278, t10 * t58 + t40 * t32 + t9 * t53 + t36 * t37 + t11 * t57 + t39 * t31 - g(1) * (-t245 - t260) - g(2) * (t147 + t221) + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -t203 + t204, (-t227 + (-qJD(2) + t150) * t231) * pkin(1) + t271, 0, 0, t90, t66, t114, t91, t115, 0, t172 * t158 + (-t177 - t282) * t161 + t252, t172 * t161 + (t177 - t133) * t158 + t226, pkin(7) * t210 - t254 * t281 + t205, -t88 * pkin(2) + pkin(7) * t215 - g(2) * t236 - t280 * t254 - t197, t16, -t1, t61, t192, -t62, 0, -t84 * t224 - t136 * t49 + (t84 * t265 - t256) * qJD(3) + t184 + t277, -t273 * t224 - t136 * t50 + (t265 * t273 - t255) * qJD(3) - t191 + t276, t168 + t268, t13 * t70 - t12 * t69 - t54 * t136 - g(1) * t189 - g(2) * t207 + (-t224 + t139) * t83 + t255 * t42 - t256 * t41, t16, t61, t1, 0, t62, t192, -t256 * qJD(3) - t200 * t84 + t59 * t49 + t184 + t279, t168 + t269, t255 * qJD(3) + t200 * t273 - t59 * t50 + t191 + t278, g(1) * t245 - g(2) * t221 + t10 * t70 + t11 * t69 - t200 * t36 + t255 * t40 + t256 * t39 + t9 * t59 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, t233 * t148, t242, t222, t241, qJDD(3), t158 * t180 - t263, g(3) * t158 + t161 * t180, 0, 0, t257, t274, t34, -t257, -t193, qJDD(3), t240 - t83 * t273 + (t251 * qJDD(3) - t84 * t244) * pkin(3) + t187, t45 * qJD(3) + t83 * t84 + (-qJDD(3) * t156 - t244 * t273) * pkin(3) + t267, (t42 - t44) * t273 + (-t41 + t45) * t84 + (-t156 * t49 - t251 * t50) * pkin(3), t41 * t44 - t42 * t45 + (t251 * t12 - t263 + t13 * t156 + (-t150 * t83 + t271) * t158) * pkin(3), t257, t34, -t274, qJDD(3), t193, -t257, t240 - t46 * t84 + (pkin(4) - t127) * qJDD(3) + t174, -t125 * t49 + t127 * t50 + (t40 - t44) * t273 + (t39 - t238) * t84, t125 * qJDD(3) - t36 * t84 + t46 * t273 + t152 + (0.2e1 * qJD(5) - t45) * qJD(3) - t267, t10 * t125 + t11 * t127 - t36 * t46 - t39 * t44 - g(3) * (-t190 + t259) + t238 * t40 + t271 * (pkin(4) * t141 - qJ(5) * t142 + t265); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t33, t275, t273 * t41 + t42 * t84 + t270 + t54, 0, 0, 0, 0, 0, 0, t166, t275, -t33, t40 * t84 + (-qJD(5) - t39) * t273 + t167 + t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t257, t34, -t266 - t164, -t40 * qJD(3) - t174 - t250;];
tau_reg = t2;
