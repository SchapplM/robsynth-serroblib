% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:36
% EndTime: 2019-03-09 12:48:45
% DurationCPUTime: 3.27s
% Computational Cost: add. (4517->380), mult. (10120->518), div. (0->0), fcn. (6425->6), ass. (0->212)
t176 = cos(qJ(2));
t172 = sin(qJ(4));
t173 = sin(qJ(2));
t253 = t172 * t173;
t195 = pkin(4) * t176 - pkin(9) * t253;
t243 = qJD(1) * t176;
t158 = pkin(7) * t243;
t129 = pkin(3) * t243 + t158;
t175 = cos(qJ(4));
t244 = qJD(1) * t173;
t162 = pkin(2) * t244;
t202 = pkin(8) * t173 - qJ(3) * t176;
t93 = qJD(1) * t202 + t162;
t210 = t175 * t129 - t172 * t93;
t239 = qJD(4) * t172;
t177 = -pkin(2) - pkin(8);
t276 = pkin(9) - t177;
t298 = -qJD(1) * t195 + t276 * t239 - t210;
t134 = t276 * t175;
t226 = t175 * t244;
t261 = t172 * t129 + t175 * t93;
t294 = pkin(9) * t226 + qJD(4) * t134 + t261;
t234 = t172 * qJD(2);
t120 = -t175 * t243 - t234;
t222 = t172 * t243;
t241 = qJD(2) * t175;
t121 = -t222 + t241;
t171 = sin(qJ(5));
t174 = cos(qJ(5));
t200 = t120 * t171 + t174 * t121;
t279 = t200 ^ 2;
t62 = -t174 * t120 + t121 * t171;
t60 = t62 ^ 2;
t297 = -t60 + t279;
t278 = pkin(3) + pkin(7);
t122 = t171 * t175 + t172 * t174;
t187 = t122 * t173;
t284 = qJD(4) + qJD(5);
t68 = t284 * t122;
t270 = qJD(1) * t187 + t68;
t236 = qJD(5) * t171;
t251 = t174 * t175;
t255 = t171 * t172;
t269 = -t171 * t239 - t172 * t236 + t174 * t226 - t244 * t255 + t251 * t284;
t296 = qJ(6) * t62;
t295 = t200 * t62;
t157 = pkin(7) * t244;
t293 = qJD(3) + t157;
t151 = qJD(4) + t244;
t143 = qJD(5) + t151;
t238 = qJD(4) * t175;
t155 = qJD(2) * t238;
t237 = qJD(4) * t176;
t225 = t172 * t237;
t186 = t173 * t241 + t225;
t182 = qJD(1) * t186 - t155;
t235 = qJD(5) * t174;
t231 = qJD(1) * qJD(2);
t220 = t173 * t231;
t77 = qJD(4) * t120 + t172 * t220;
t29 = -t120 * t235 + t121 * t236 - t171 * t182 - t174 * t77;
t292 = t143 * t62 - t29;
t154 = t176 * t231;
t149 = pkin(7) * t154;
t112 = pkin(3) * t154 + t149;
t150 = pkin(2) * t220;
t233 = t173 * qJD(3);
t185 = qJD(2) * t202 - t233;
t72 = qJD(1) * t185 + t150;
t212 = t175 * t112 - t172 * t72;
t252 = t173 * qJ(3);
t219 = -pkin(1) - t252;
t286 = t176 * t177;
t113 = t219 + t286;
t83 = t113 * qJD(1);
t232 = pkin(3) * t244 + t293;
t88 = t177 * qJD(2) + t232;
t48 = t172 * t88 + t175 * t83;
t183 = -qJD(4) * t48 + t212;
t15 = pkin(4) * t154 - t77 * pkin(9) + t183;
t230 = -t172 * t112 - t175 * t72 - t88 * t238;
t193 = -t239 * t83 - t230;
t19 = pkin(9) * t182 + t193;
t47 = -t172 * t83 + t175 * t88;
t39 = -pkin(9) * t121 + t47;
t33 = pkin(4) * t151 + t39;
t40 = pkin(9) * t120 + t48;
t211 = -t171 * t15 - t174 * t19 - t33 * t235 + t40 * t236;
t168 = qJD(2) * qJ(3);
t106 = t168 + t129;
t70 = -pkin(4) * t120 + t106;
t291 = t62 * t70 + t211;
t290 = -0.2e1 * t231;
t35 = pkin(5) * t62 + qJD(6) + t70;
t289 = t200 * t35;
t288 = t298 * t174;
t287 = qJ(6) * t200;
t198 = -t251 + t255;
t94 = t198 * t176;
t133 = t276 * t172;
t246 = -t174 * t133 - t171 * t134;
t169 = t173 ^ 2;
t170 = t176 ^ 2;
t245 = t169 - t170;
t285 = -t133 * t236 + t134 * t235 - t298 * t171 + t294 * t174;
t227 = pkin(4) * t175 + pkin(3);
t242 = qJD(2) * t173;
t71 = (-pkin(7) - t227) * t242 - pkin(4) * t225;
t38 = t174 * t40;
t11 = t171 * t33 + t38;
t217 = t174 * t15 - t171 * t19;
t184 = -qJD(5) * t11 + t217;
t283 = -t70 * t200 + t184;
t30 = qJD(5) * t200 + t171 * t77 - t174 * t182;
t282 = t143 * t200 - t30;
t281 = t198 * t29 - t200 * t270;
t262 = pkin(4) * t238 + t227 * t244 + t293;
t1 = pkin(5) * t154 + t29 * qJ(6) - qJD(6) * t200 + t184;
t2 = -qJ(6) * t30 - qJD(6) * t62 - t211;
t36 = t171 * t40;
t10 = t174 * t33 - t36;
t6 = t10 - t287;
t5 = pkin(5) * t143 + t6;
t7 = t11 - t296;
t280 = -t1 * t198 + t2 * t122 + t269 * t7 - t270 * t5;
t277 = t5 - t6;
t275 = -qJ(6) * t269 - qJD(6) * t122 - t285;
t274 = -pkin(5) * t243 + qJ(6) * t270 - t246 * qJD(5) + t198 * qJD(6) + t171 * t294 + t288;
t273 = t174 * t39 - t36;
t138 = t278 * t173;
t125 = t175 * t138;
t218 = pkin(9) * t176 - t113;
t56 = t173 * pkin(4) + t172 * t218 + t125;
t250 = t175 * t176;
t124 = t172 * t138;
t260 = t175 * t113 + t124;
t59 = -pkin(9) * t250 + t260;
t271 = t171 * t56 + t174 * t59;
t268 = qJD(2) * pkin(2);
t266 = t77 * t172;
t265 = t77 * t175;
t128 = t278 * t242;
t167 = qJD(2) * qJD(3);
t90 = -qJD(1) * t128 + t167;
t264 = t90 * t172;
t263 = t90 * t175;
t259 = t106 * t173;
t258 = t121 * t176;
t257 = t151 * t173;
t256 = t151 * t177;
t254 = t172 * t120;
t179 = qJD(1) ^ 2;
t249 = t176 * t179;
t178 = qJD(2) ^ 2;
t248 = t178 * t173;
t247 = t178 * t176;
t153 = t172 * pkin(4) + qJ(3);
t139 = t278 * t176;
t135 = -pkin(2) * t176 + t219;
t107 = qJD(1) * t135;
t240 = qJD(2) * t176;
t229 = t173 * t249;
t228 = t175 * t257;
t104 = pkin(4) * t250 + t139;
t224 = t151 * t238;
t223 = t175 * t237;
t221 = t269 * t143;
t130 = t278 * t240;
t161 = pkin(2) * t242;
t79 = t161 + t185;
t209 = t175 * t130 - t172 * t79;
t26 = t195 * qJD(2) + (t175 * t218 - t124) * qJD(4) + t209;
t189 = -t113 * t239 + t172 * t130 + t138 * t238 + t175 * t79;
t28 = pkin(9) * t186 + t189;
t216 = -t171 * t28 + t174 * t26;
t215 = -t171 * t39 - t38;
t214 = -t171 * t59 + t174 * t56;
t208 = pkin(1) * t290;
t207 = qJD(3) - t268;
t206 = t133 * t171 - t174 * t134;
t204 = -t143 * t270 - t198 * t154;
t199 = t121 * t175 + t254;
t197 = -0.2e1 * qJD(2) * t107;
t196 = t151 * t172;
t190 = -qJ(3) * t240 - t233;
t81 = qJD(1) * t190 + t150;
t98 = t161 + t190;
t194 = pkin(7) * t178 + qJD(1) * t98 + t81;
t192 = t171 * t26 + t174 * t28 + t56 * t235 - t236 * t59;
t131 = pkin(7) * t220 - t167;
t132 = t157 + t207;
t137 = -t158 - t168;
t180 = -t131 * t176 + (t132 * t176 + (t137 + t158) * t173) * qJD(2);
t53 = t155 * pkin(4) + qJD(1) * t71 + t167;
t12 = t30 * pkin(5) + t53;
t156 = pkin(4) * t174 + pkin(5);
t142 = t175 * t154;
t141 = t173 * t154;
t126 = -qJ(3) * t243 + t162;
t95 = t122 * t176;
t87 = t107 * t244;
t55 = -qJ(6) * t122 + t246;
t54 = qJ(6) * t198 + t206;
t42 = t176 * t68 - t198 * t242;
t41 = qJD(2) * t187 + t284 * t94;
t23 = qJ(6) * t94 + t271;
t22 = pkin(5) * t173 + qJ(6) * t95 + t214;
t9 = t273 - t287;
t8 = t215 + t296;
t4 = qJ(6) * t42 + qJD(6) * t94 + t192;
t3 = pkin(5) * t240 - t41 * qJ(6) - qJD(5) * t271 + t95 * qJD(6) + t216;
t13 = [0, 0, 0, 0.2e1 * t141, t245 * t290, t247, -t248, 0, -pkin(7) * t247 + t173 * t208, pkin(7) * t248 + t176 * t208, t180, t173 * t197 + t176 * t194, -t173 * t194 + t176 * t197, pkin(7) * t180 + t107 * t98 + t81 * t135, -t176 * t266 + (t173 * t234 - t223) * t121, t199 * t242 + (-t172 * (t175 * t220 - t155) - t265 + (-t175 * t120 + (t121 - t222) * t172) * qJD(4)) * t176, -t151 * t223 + t77 * t173 + (t258 + (-qJD(1) * t170 + t257) * t172) * qJD(2), t151 * t225 + (qJD(4) * t222 - t155) * t173 + (t120 * t176 + (qJD(1) * t245 + t257) * t175) * qJD(2), t151 * t240 + t141, t209 * t151 + t128 * t120 + t139 * t155 + ((-qJD(1) * t139 - t106) * t241 + t212) * t173 + (-t151 * t260 - t173 * t48) * qJD(4) + (-t106 * t239 + t47 * qJD(2) + t263 + ((-t113 * t172 + t125) * qJD(2) - t139 * t239) * qJD(1)) * t176, -t189 * t151 - t128 * t121 + t139 * t77 + ((qJD(2) * t106 + qJD(4) * t83) * t172 + t230) * t173 + (-t106 * t238 - t264 + (-qJD(1) * t260 - t48) * qJD(2)) * t176, t200 * t41 + t29 * t95, t200 * t42 - t29 * t94 + t30 * t95 - t41 * t62, t143 * t41 - t173 * t29 + (-qJD(1) * t95 + t200) * t240, t143 * t42 - t173 * t30 + (qJD(1) * t94 - t62) * t240, t143 * t240 + t141, t216 * t143 + t217 * t173 + t71 * t62 + t104 * t30 - t53 * t94 - t70 * t42 + (-t11 * t173 - t143 * t271) * qJD(5) + (qJD(1) * t214 + t10) * t240, -t192 * t143 + t211 * t173 + t71 * t200 - t104 * t29 - t53 * t95 + t70 * t41 + (-qJD(1) * t271 - t11) * t240, t1 * t95 + t2 * t94 - t200 * t3 + t22 * t29 - t23 * t30 - t4 * t62 - t41 * t5 + t42 * t7, t2 * t23 + t7 * t4 + t1 * t22 + t5 * t3 + t12 * (-pkin(5) * t94 + t104) + (-t42 * pkin(5) + t71) * t35; 0, 0, 0, -t229, t245 * t179, 0, 0, 0, t179 * pkin(1) * t173, pkin(1) * t249 ((-t137 - t168) * t173 + (-t132 + t207) * t176) * qJD(1), -t126 * t243 + t87, 0.2e1 * t167 + (t107 * t176 + t126 * t173) * qJD(1), -t131 * qJ(3) - t137 * qJD(3) - t107 * t126 + (-t137 * t173 + (-t132 - t268) * t176) * qJD(1) * pkin(7), -t121 * t196 + t265, -t175 * t155 - t266 - t199 * qJD(4) + (t172 * t223 + (-t254 + (-t121 + t241) * t175) * t173) * qJD(1), -t151 * t239 + t142 + (-t151 * t253 - t258) * qJD(1), -t224 + (-t228 + (-t120 - t234) * t176) * qJD(1), -t151 * t243, qJ(3) * t155 + t264 - t210 * t151 - t232 * t120 + (t106 * t175 - t172 * t256) * qJD(4) + ((-qJ(3) * t239 - t47) * t176 + (t259 + (-t252 + t286) * qJD(2)) * t175) * qJD(1), qJ(3) * t77 + t263 + t261 * t151 + t232 * t121 + (-t106 * t172 - t175 * t256) * qJD(4) + (t48 * t176 + (-t177 * t240 - t259) * t172) * qJD(1), t281, t29 * t122 + t198 * t30 - t200 * t269 + t270 * t62, -t200 * t243 + t204, -t221 + (-qJD(2) * t122 + t62) * t243, -t143 * t243, t53 * t122 + t153 * t30 + t269 * t70 + t262 * t62 + (t133 * t235 + (qJD(5) * t134 + t294) * t171 + t288) * t143 + (qJD(2) * t206 - t10) * t243, -t53 * t198 - t153 * t29 - t270 * t70 + t262 * t200 + t285 * t143 + (-qJD(2) * t246 + t11) * t243, -t200 * t274 - t275 * t62 + t54 * t29 - t55 * t30 - t280, t2 * t55 + t1 * t54 + t12 * (pkin(5) * t122 + t153) + t275 * t7 + t274 * t5 + (pkin(5) * t269 + t262) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, -t169 * t179 - t178, qJD(2) * t137 + t149 + t87, 0, 0, 0, 0, 0, qJD(2) * t120 - t151 * t196 + t142, -t224 - qJD(2) * t121 + (-t176 * t234 - t228) * qJD(1), 0, 0, 0, 0, 0, -qJD(2) * t62 + t204, -t221 + (-t122 * t243 - t200) * qJD(2), -t122 * t30 - t269 * t62 - t281, -t35 * qJD(2) + t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 * t120, -t120 ^ 2 + t121 ^ 2, -t120 * t151 + t77, t121 * t151 + t182, t154, -t106 * t121 + t48 * t151 + t183, -t106 * t120 + t151 * t47 - t193, t295, t297, t292, t282, t154, -t215 * t143 + (-t121 * t62 - t143 * t236 + t154 * t174) * pkin(4) + t283, t273 * t143 + (-t121 * t200 - t143 * t235 - t154 * t171) * pkin(4) + t291, t156 * t29 - t5 * t62 + t7 * t200 + t9 * t62 + t8 * t200 + (-t171 * t30 + (t171 * t200 - t174 * t62) * qJD(5)) * pkin(4), -pkin(5) * t289 + t1 * t156 - t5 * t8 - t7 * t9 + (-t35 * t121 + t2 * t171 + (-t171 * t5 + t174 * t7) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, t297, t292, t282, t154, t11 * t143 + t283, t10 * t143 + t291, pkin(5) * t29 - t277 * t62, t277 * t7 + (t1 - t289) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 - t279, t200 * t5 + t7 * t62 + t12;];
tauc_reg  = t13;
