% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:20
% EndTime: 2019-03-09 00:04:29
% DurationCPUTime: 3.11s
% Computational Cost: add. (5868->369), mult. (14390->503), div. (0->0), fcn. (10826->10), ass. (0->212)
t153 = sin(qJ(4));
t154 = sin(qJ(3));
t157 = cos(qJ(3));
t273 = cos(qJ(4));
t176 = -t153 * t154 + t273 * t157;
t158 = cos(qJ(2));
t150 = sin(pkin(6));
t230 = qJD(1) * t150;
t211 = t158 * t230;
t105 = t176 * t211;
t277 = -pkin(9) - pkin(8);
t216 = qJD(3) * t277;
t128 = t154 * t216;
t129 = t157 * t216;
t133 = t277 * t154;
t134 = t277 * t157;
t280 = t273 * t133 + t153 * t134;
t62 = t280 * qJD(4) + t273 * t128 + t153 * t129;
t292 = -t62 + t105;
t235 = t153 * t157;
t126 = t273 * t154 + t235;
t221 = qJD(3) + qJD(4);
t98 = t221 * t126;
t291 = t98 * qJD(2);
t156 = cos(qJ(5));
t223 = qJD(5) * t156;
t207 = qJD(2) * t273;
t227 = qJD(2) * t154;
t121 = t153 * t227 - t157 * t207;
t241 = t121 * t156;
t290 = t223 + t241;
t152 = sin(qJ(5));
t224 = qJD(5) * t152;
t242 = t121 * t152;
t289 = t224 + t242;
t117 = qJD(5) + t121;
t155 = sin(qJ(2));
t212 = t155 * t230;
t202 = -t277 * qJD(2) + t212;
t151 = cos(pkin(6));
t229 = qJD(1) * t151;
t103 = t154 * t229 + t202 * t157;
t206 = qJD(4) * t273;
t225 = qJD(4) * t153;
t102 = -t202 * t154 + t157 * t229;
t228 = qJD(2) * t150;
t205 = qJD(1) * t228;
t195 = t158 * t205;
t75 = t102 * qJD(3) + t157 * t195;
t76 = -t103 * qJD(3) - t154 * t195;
t263 = qJD(3) * pkin(3);
t96 = t102 + t263;
t16 = -t103 * t225 + t153 * t76 + t96 * t206 + t273 * t75;
t222 = qJD(2) * qJD(3);
t204 = t154 * t222;
t118 = pkin(3) * t204 + t155 * t205;
t97 = t221 * t176;
t162 = t97 * qJD(2);
t42 = pkin(4) * t291 - pkin(10) * t162 + t118;
t95 = t273 * t103;
t57 = t153 * t96 + t95;
t50 = t221 * pkin(10) + t57;
t146 = -pkin(3) * t157 - pkin(2);
t115 = t146 * qJD(2) - t211;
t123 = -qJD(2) * t235 - t154 * t207;
t71 = pkin(4) * t121 + pkin(10) * t123 + t115;
t178 = t152 * t42 + t156 * t16 + t71 * t223 - t50 * t224;
t265 = qJ(6) * t291;
t2 = qJD(6) * t117 + t178 + t265;
t203 = t152 * t16 - t156 * t42 + t50 * t223 + t71 * t224;
t276 = pkin(5) * t291;
t4 = t203 - t276;
t288 = t4 * t152 + t2 * t156;
t256 = t156 * t291;
t287 = (t117 * t224 - t256) * pkin(10);
t111 = t153 * t133 - t273 * t134;
t91 = -pkin(4) * t176 - pkin(10) * t126 + t146;
t252 = t156 * t111 + t152 * t91;
t219 = t154 * t263;
t52 = pkin(4) * t98 - pkin(10) * t97 + t219;
t286 = t111 * t224 - t91 * t223 + t292 * t156 + (t212 - t52) * t152;
t28 = t152 * t71 + t156 * t50;
t20 = qJ(6) * t117 + t28;
t200 = t156 * t221;
t106 = -t123 * t152 - t200;
t172 = t156 * t123 - t152 * t221;
t94 = t153 * t103;
t56 = t273 * t96 - t94;
t49 = -t221 * pkin(4) - t56;
t31 = t106 * pkin(5) + qJ(6) * t172 + t49;
t285 = t20 * t123 - t31 * t241;
t17 = t103 * t206 + t153 * t75 + t96 * t225 - t273 * t76;
t53 = -qJD(5) * t200 - t123 * t224 - t156 * t162;
t54 = -qJD(5) * t172 + t152 * t162;
t5 = pkin(5) * t54 + qJ(6) * t53 + qJD(6) * t172 + t17;
t284 = -t5 * t156 + t31 * t224;
t283 = t17 * t156 - t49 * t224;
t58 = t153 * t102 + t95;
t194 = pkin(3) * t225 - t58;
t249 = t111 * qJD(4) - t126 * t211 + t153 * t128 - t273 * t129;
t282 = t289 * pkin(5) - t290 * qJ(6) - qJD(6) * t152;
t281 = t212 - t219;
t279 = t172 ^ 2;
t278 = t117 ^ 2;
t275 = qJ(6) * t98 - qJD(6) * t176 - t286;
t78 = t105 * t152 - t156 * t212;
t274 = pkin(5) * t98 - t252 * qJD(5) - t152 * t62 + t156 * t52 + t78;
t272 = pkin(5) * t123;
t271 = t5 * t152;
t190 = pkin(5) * t152 - qJ(6) * t156;
t191 = t156 * pkin(5) + t152 * qJ(6);
t269 = -t190 * t97 - (t191 * qJD(5) - qJD(6) * t156) * t126 - t249;
t90 = -pkin(4) * t123 + pkin(10) * t121;
t268 = t152 * t90 + t156 * t56;
t59 = t273 * t102 - t94;
t80 = pkin(3) * t227 + t90;
t267 = t152 * t80 + t156 * t59;
t266 = pkin(3) * qJD(4);
t264 = qJD(2) * pkin(2);
t262 = t172 * t31;
t144 = pkin(3) * t153 + pkin(10);
t261 = t144 * t291;
t260 = t152 * t53;
t259 = t152 * t291;
t258 = t152 * t97;
t257 = t156 * t54;
t255 = t156 * t97;
t251 = t282 + t194;
t250 = -t57 + t282;
t248 = t106 * t117;
t247 = t106 * t152;
t246 = t172 * t106;
t245 = t172 * t117;
t244 = t172 * t156;
t243 = t117 * t123;
t240 = t123 * t121;
t239 = t126 * t156;
t238 = t150 * t155;
t237 = t150 * t158;
t160 = qJD(2) ^ 2;
t236 = t150 * t160;
t159 = qJD(3) ^ 2;
t234 = t159 * t154;
t233 = t159 * t157;
t27 = -t152 * t50 + t156 * t71;
t232 = qJD(6) - t27;
t231 = t154 ^ 2 - t157 ^ 2;
t226 = qJD(2) * t155;
t220 = t273 * pkin(3);
t217 = t155 * t236;
t215 = t152 * t273;
t214 = t156 * t273;
t209 = t150 * t226;
t208 = t158 * t228;
t201 = t117 * t156;
t199 = pkin(3) * t206;
t198 = t154 * t208;
t197 = t157 * t208;
t196 = -t28 * t123 + t17 * t152 + t49 * t223;
t189 = t49 * t121 - t261;
t19 = -pkin(5) * t117 + t232;
t188 = -t152 * t20 + t156 * t19;
t187 = -t152 * t56 + t156 * t90;
t132 = -pkin(4) - t191;
t185 = -t19 * t123 + t284;
t184 = t27 * t123 - t283;
t183 = t31 * t223 + t271;
t182 = t117 * t28 - t203;
t119 = t151 * t157 - t154 * t238;
t120 = t151 * t154 + t157 * t238;
t82 = t153 * t119 + t273 * t120;
t69 = t152 * t82 + t156 * t237;
t70 = -t152 * t237 + t156 * t82;
t181 = t115 * t123 - t17;
t179 = -t126 * t224 + t255;
t177 = t273 * t119 - t153 * t120;
t174 = t264 * qJD(2);
t100 = -qJD(3) * t120 - t198;
t99 = qJD(3) * t119 + t197;
t37 = t177 * qJD(4) + t153 * t100 + t273 * t99;
t11 = qJD(5) * t70 + t152 * t37 - t156 * t209;
t38 = t82 * qJD(4) - t273 * t100 + t153 * t99;
t173 = t38 * t106 - t11 * t117 - t177 * t54 - t291 * t69;
t171 = (-t117 * t223 - t259) * pkin(10);
t170 = -0.2e1 * qJD(3) * t264;
t169 = -t144 * t224 + t156 * t199;
t10 = -qJD(5) * t69 + t152 * t209 + t156 * t37;
t168 = t10 * t117 + t172 * t38 - t177 * t53 + t291 * t70;
t167 = t290 * t19 - t289 * t20 + t288;
t166 = qJD(5) * t188 + t288;
t165 = -t260 - t257 + (-t244 + t247) * qJD(5);
t164 = t115 * t121 - t16;
t145 = -t220 - pkin(4);
t124 = -t220 + t132;
t114 = t123 * qJ(6);
t72 = -t121 ^ 2 + t123 ^ 2;
t68 = -t123 * t221 - t291;
t67 = t121 * t221 + t162;
t64 = -pkin(5) * t172 + qJ(6) * t106;
t61 = t190 * t126 - t280;
t44 = pkin(5) * t176 + t111 * t152 - t156 * t91;
t43 = -qJ(6) * t176 + t252;
t34 = -t53 + t248;
t33 = -t187 + t272;
t32 = -t114 + t268;
t30 = t152 * t59 - t156 * t80 + t272;
t29 = -t114 + t267;
t24 = t117 * t201 - t123 * t172 + t259;
t23 = -t106 * t123 - t278 * t152 + t256;
t21 = -t172 * t201 - t260;
t6 = (-t53 - t248) * t156 + (-t54 + t245) * t152;
t1 = [0, 0, -t217, -t158 * t236, 0, 0, 0, 0, 0, -t157 * t217 + (t100 - t198) * qJD(3), t154 * t217 + (-t99 - t197) * qJD(3), 0, 0, 0, 0, 0, -t38 * t221 + (t121 * t226 - t158 * t291) * t150, -t37 * t221 + (-t155 * t123 - t158 * t97) * t228, 0, 0, 0, 0, 0, t173, -t168, t173, -t10 * t106 - t11 * t172 - t53 * t69 - t54 * t70, t168, t10 * t20 + t11 * t19 - t177 * t5 + t2 * t70 + t31 * t38 + t4 * t69; 0, 0, 0, 0, 0.2e1 * t157 * t204, -0.2e1 * t231 * t222, t233, -t234, 0, -pkin(8) * t233 + t154 * t170, pkin(8) * t234 + t157 * t170, -t123 * t97 + t126 * t162, -t97 * t121 + t123 * t98 - t126 * t291 + t162 * t176, t97 * t221, -t98 * t221, 0, t115 * t98 - t118 * t176 - t281 * t121 + t146 * t291 - t249 * t221, t115 * t97 + t118 * t126 + t281 * t123 + t146 * t162 + t292 * t221, -t172 * t179 - t239 * t53 (-t106 * t156 + t152 * t172) * t97 + (t260 - t257 + (t244 + t247) * qJD(5)) * t126, t117 * t179 - t172 * t98 + t176 * t53 + t239 * t291, -t126 * t259 - t106 * t98 + t176 * t54 + (-t126 * t223 - t258) * t117, t117 * t98 - t176 * t291, t203 * t176 + t27 * t98 - t280 * t54 + t78 * t117 + t249 * t106 + ((-qJD(5) * t111 + t52) * t117 + t91 * t291 + t49 * qJD(5) * t126) * t156 + ((-qJD(5) * t91 - t62) * t117 - t111 * t291 + t17 * t126 + t49 * t97) * t152, t286 * t117 + t283 * t126 - t172 * t249 + t176 * t178 - t252 * t291 + t49 * t255 - t28 * t98 + t280 * t53, -t269 * t106 + t274 * t117 + t183 * t126 + t176 * t4 - t19 * t98 + t31 * t258 - t291 * t44 + t54 * t61, -t43 * t54 - t44 * t53 + t188 * t97 + t274 * t172 - t275 * t106 + (-t152 * t2 + t156 * t4 + (-t152 * t19 - t156 * t20) * qJD(5)) * t126, t275 * t117 + t284 * t126 - t172 * t269 - t176 * t2 + t20 * t98 - t31 * t255 + t291 * t43 + t53 * t61, -t274 * t19 + t2 * t43 + t275 * t20 - t269 * t31 + t4 * t44 + t5 * t61; 0, 0, 0, 0, -t154 * t160 * t157, t231 * t160, 0, 0, 0, t154 * t174, t157 * t174, -t240, t72, t67, t68, 0, t58 * t221 + (-t121 * t227 - t221 * t225) * pkin(3) + t181, t59 * t221 + (t123 * t227 - t206 * t221) * pkin(3) + t164, t21, t6, t24, t23, t243, t145 * t54 + t189 * t152 + t194 * t106 + ((-qJD(5) * t144 - t80) * t156 + (-t199 + t59) * t152) * t117 + t184, -t145 * t53 + t189 * t156 - t194 * t172 + (-t169 + t267) * t117 + t196, t124 * t54 + (t121 * t31 - t261) * t152 + t251 * t106 + (-t144 * t223 - t152 * t199 + t30) * t117 + t185, t29 * t106 + t30 * t172 + (-t106 * t214 - t172 * t215) * t266 + t165 * t144 + t167, t124 * t53 - t271 + (-qJD(5) * t31 + t261) * t156 + t251 * t172 + (-t29 + t169) * t117 + t285, t5 * t124 - t19 * t30 - t20 * t29 + t251 * t31 + (t19 * t215 + t20 * t214) * t266 + t166 * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, t72, t67, t68, 0, t221 * t57 + t181, t221 * t56 + t164, t21, t6, t24, t23, t243, -pkin(4) * t54 - t57 * t106 - t117 * t187 + t242 * t49 + t171 + t184, pkin(4) * t53 + t268 * t117 + t172 * t57 + t49 * t241 + t196 + t287, t106 * t250 + t117 * t33 + t132 * t54 + t242 * t31 + t171 + t185, pkin(10) * t165 + t106 * t32 + t172 * t33 + t167, -t117 * t32 + t132 * t53 + t172 * t250 - t183 + t285 - t287, pkin(10) * t166 + t132 * t5 - t19 * t33 - t20 * t32 + t250 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, -t106 ^ 2 + t279, t34, -t54 - t245, t291, t172 * t49 + t182, t106 * t49 + t117 * t27 - t178, -t106 * t64 + t182 + t262 + 0.2e1 * t276, pkin(5) * t53 - qJ(6) * t54 - (t20 - t28) * t172 + (t19 - t232) * t106, 0.2e1 * t265 - t106 * t31 - t172 * t64 + (0.2e1 * qJD(6) - t27) * t117 + t178, -pkin(5) * t4 + qJ(6) * t2 - t19 * t28 + t20 * t232 - t31 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291 - t246, t34, -t278 - t279, -t117 * t20 - t262 + t4;];
tauc_reg  = t1;
