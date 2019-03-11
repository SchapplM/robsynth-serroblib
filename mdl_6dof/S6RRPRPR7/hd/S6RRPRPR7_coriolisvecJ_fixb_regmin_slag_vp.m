% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:43
% EndTime: 2019-03-09 10:47:53
% DurationCPUTime: 3.61s
% Computational Cost: add. (4654->332), mult. (11108->451), div. (0->0), fcn. (7845->8), ass. (0->199)
t152 = qJD(2) - qJD(4);
t162 = cos(qJ(6));
t157 = sin(pkin(10));
t244 = cos(pkin(10));
t160 = sin(qJ(4));
t161 = sin(qJ(2));
t163 = cos(qJ(4));
t164 = cos(qJ(2));
t111 = t161 * t160 + t164 * t163;
t94 = t111 * qJD(1);
t233 = qJD(1) * t161;
t220 = t163 * t233;
t232 = qJD(1) * t164;
t96 = -t160 * t232 + t220;
t187 = -t157 * t94 + t244 * t96;
t226 = qJD(6) * t162;
t159 = sin(qJ(6));
t227 = qJD(6) * t159;
t229 = qJD(4) * t160;
t219 = t164 * t229;
t231 = qJD(2) * t161;
t182 = t163 * t231 + t219;
t225 = qJD(1) * qJD(2);
t216 = t164 * t225;
t237 = qJD(4) * t220 + t160 * t216;
t174 = qJD(1) * t182 - t237;
t183 = t111 * qJD(4);
t217 = t161 * t225;
t61 = -qJD(1) * t183 + t160 * t217 + t163 * t216;
t34 = t157 * t174 + t244 * t61;
t14 = -t152 * t226 + t162 * t34 - t187 * t227;
t196 = t159 * t152 - t162 * t187;
t15 = -qJD(6) * t196 + t159 * t34;
t58 = t157 * t96 + t244 * t94;
t290 = qJD(6) + t58;
t279 = t196 * t290;
t44 = t162 * t152 + t159 * t187;
t283 = t290 * t44;
t300 = (t14 - t283) * t162 + (-t15 + t279) * t159;
t264 = t196 * t187;
t33 = t157 * t61 - t244 * t174;
t252 = t159 * t33;
t293 = t162 * t290;
t277 = -t290 * t293 - t252;
t297 = t264 - t277;
t144 = pkin(7) * t233;
t296 = -pkin(8) * t233 + qJD(3) + t144;
t255 = t14 * t159;
t295 = t196 * t293 - t255;
t145 = pkin(7) * t232;
t117 = -pkin(8) * t232 + t145;
t165 = -pkin(2) - pkin(3);
t120 = t163 * qJ(3) + t160 * t165;
t292 = qJD(4) * t120 + t163 * t117 + t296 * t160;
t195 = -t160 * qJ(3) + t163 * t165;
t291 = qJD(4) * t195 - t160 * t117 + t296 * t163;
t289 = -t96 * t152 + t174;
t288 = qJD(6) - t290;
t287 = pkin(5) * t187 + t58 * pkin(9);
t221 = t165 * qJD(2);
t79 = t221 + t296;
t154 = qJD(2) * qJ(3);
t97 = t117 + t154;
t212 = -t160 * t97 + t163 * t79;
t247 = t96 * qJ(5);
t41 = t212 - t247;
t38 = -t152 * pkin(4) + t41;
t198 = t160 * t79 + t163 * t97;
t249 = t94 * qJ(5);
t42 = t198 - t249;
t39 = t244 * t42;
t17 = t157 * t38 + t39;
t13 = -t152 * pkin(9) + t17;
t98 = -qJD(1) * pkin(1) - pkin(2) * t232 - qJ(3) * t233;
t78 = pkin(3) * t232 - t98;
t54 = t94 * pkin(4) + qJD(5) + t78;
t20 = pkin(5) * t58 - pkin(9) * t187 + t54;
t4 = t162 * t13 + t159 * t20;
t285 = t4 * t187;
t284 = -t247 + t291;
t263 = t187 * t44;
t282 = -t249 + t292;
t281 = -t94 * t152 + t61;
t199 = t159 * t13 - t162 * t20;
t280 = t187 * t199;
t278 = t290 * t187;
t135 = pkin(7) * t216;
t105 = -pkin(8) * t216 + t135;
t268 = pkin(7) - pkin(8);
t116 = t268 * t231;
t153 = qJD(2) * qJD(3);
t84 = -qJD(1) * t116 + t153;
t172 = -t198 * qJD(4) + t163 * t105 - t160 * t84;
t276 = t78 * t96 - t172;
t275 = -0.2e1 * t225;
t228 = qJD(4) * t163;
t188 = t160 * t105 + t163 * t84 + t79 * t228 - t97 * t229;
t11 = qJ(5) * t174 - t94 * qJD(5) + t188;
t169 = -t61 * qJ(5) - t96 * qJD(5) + t172;
t1 = t157 * t11 - t244 * t169;
t138 = t157 * pkin(4) + pkin(9);
t267 = t96 * pkin(4);
t271 = (qJD(6) * t138 + t267 + t287) * t290 + t1;
t142 = qJ(3) * t232;
t86 = t165 * t233 + t142;
t180 = t86 - t267;
t114 = -pkin(4) + t195;
t68 = t157 * t114 + t244 * t120;
t66 = -pkin(9) + t68;
t270 = (qJD(6) * t66 + t180 - t287) * t290 - t1;
t254 = t157 * t42;
t16 = t244 * t38 - t254;
t12 = t152 * pkin(5) - t16;
t2 = t244 * t11 + t157 * t169;
t112 = -t164 * t160 + t161 * t163;
t124 = t268 * t161;
t125 = t268 * t164;
t181 = -t112 * qJ(5) + t163 * t124 - t160 * t125;
t194 = -t160 * t124 - t163 * t125;
t52 = -t111 * qJ(5) - t194;
t29 = t157 * t181 + t244 * t52;
t64 = -t157 * t111 + t244 * t112;
t201 = t1 * t64 - t29 * t33;
t122 = -t164 * pkin(2) - t161 * qJ(3) - pkin(1);
t107 = t164 * pkin(3) - t122;
t189 = t111 * pkin(4) + t107;
t63 = t244 * t111 + t157 * t112;
t27 = t63 * pkin(5) - t64 * pkin(9) + t189;
t230 = qJD(2) * t164;
t69 = t160 * t230 + t161 * t228 - t182;
t70 = qJD(2) * t111 - t183;
t36 = -t157 * t69 + t244 * t70;
t118 = qJD(2) * t125;
t175 = qJD(4) * t194 + t160 * t116 + t163 * t118;
t170 = -t70 * qJ(5) - t112 * qJD(5) + t175;
t184 = -t163 * t116 + t160 * t118 + t124 * t228 - t125 * t229;
t23 = -t69 * qJ(5) - t111 * qJD(5) + t184;
t8 = t157 * t170 + t244 * t23;
t269 = t12 * t36 - (qJD(6) * t27 + t8) * t290 - (qJD(6) * t20 + t2) * t63 + t201;
t266 = t12 * t64;
t265 = t27 * t33;
t261 = t96 * t94;
t260 = t157 * t284 + t244 * t282;
t259 = -t157 * t282 + t244 * t284;
t110 = t157 * t163 + t244 * t160;
t258 = t152 * t110;
t186 = -t157 * t160 + t244 * t163;
t257 = t152 * t186;
t256 = qJD(2) * pkin(2);
t250 = t159 * t290;
t30 = t162 * t33;
t167 = qJD(1) ^ 2;
t243 = t164 * t167;
t166 = qJD(2) ^ 2;
t242 = t166 * t161;
t241 = t166 * t164;
t148 = t161 * qJD(3);
t236 = qJ(3) * t216 + qJD(1) * t148;
t235 = qJ(3) * t230 + t148;
t155 = t161 ^ 2;
t234 = -t164 ^ 2 + t155;
t224 = t94 ^ 2 - t96 ^ 2;
t223 = t64 * t227;
t222 = t161 * t243;
t208 = pkin(1) * t275;
t207 = qJD(3) - t256;
t206 = qJD(1) * t122 + t98;
t203 = t152 ^ 2;
t202 = t161 * t221;
t200 = t290 * t36 + t33 * t64;
t193 = -t227 * t290 - t250 * t58 + t30;
t192 = qJD(6) * t110 + t233;
t77 = pkin(2) * t217 - t236;
t87 = pkin(2) * t231 - t235;
t191 = -pkin(7) * t166 - qJD(1) * t87 - t77;
t67 = t244 * t114 - t157 * t120;
t73 = t202 + t235;
t179 = t78 * t94 - t188;
t19 = t244 * t41 - t254;
t178 = -t138 * t33 + (t12 + t19) * t290;
t177 = t69 * pkin(4) + t73;
t176 = -t66 * t33 + (-t12 - t259) * t290;
t119 = -pkin(7) * t217 + t153;
t121 = t144 + t207;
t123 = t145 + t154;
t171 = t119 * t164 + (t121 * t164 + (-t123 + t145) * t161) * qJD(2);
t168 = t237 * pkin(4) + (-pkin(4) * t219 + (-pkin(4) * t163 + t165) * t231) * qJD(1) + t236;
t139 = -t244 * pkin(4) - pkin(5);
t113 = pkin(2) * t233 - t142;
t71 = qJD(1) * t202 + t236;
t65 = pkin(5) - t67;
t35 = t157 * t70 + t244 * t69;
t28 = t157 * t52 - t244 * t181;
t18 = t157 * t41 + t39;
t9 = t35 * pkin(5) - t36 * pkin(9) + t177;
t7 = t157 * t23 - t244 * t170;
t6 = t33 * pkin(5) - t34 * pkin(9) + t168;
t5 = t162 * t6;
t3 = [0, 0, 0, 0.2e1 * t161 * t216, t234 * t275, t241, -t242, 0, -pkin(7) * t241 + t161 * t208, pkin(7) * t242 + t164 * t208, t191 * t164 + t206 * t231, t171, t191 * t161 - t206 * t230, pkin(7) * t171 + t77 * t122 + t98 * t87, t61 * t112 + t96 * t70, -t61 * t111 + t112 * t174 - t96 * t69 - t70 * t94, -t70 * t152, t69 * t152, 0, -t107 * t174 + t71 * t111 - t175 * t152 + t78 * t69 + t73 * t94, t107 * t61 + t71 * t112 + t184 * t152 + t78 * t70 + t73 * t96, -t16 * t36 - t17 * t35 + t187 * t7 - t2 * t63 + t28 * t34 - t58 * t8 + t201, t1 * t28 - t16 * t7 + t168 * t189 + t17 * t8 + t177 * t54 + t2 * t29, t196 * t223 + (t14 * t64 - t196 * t36) * t162 (t159 * t196 - t162 * t44) * t36 + (-t255 - t15 * t162 + (t159 * t44 + t162 * t196) * qJD(6)) * t64, t14 * t63 + t162 * t200 - t196 * t35 - t223 * t290, -t226 * t290 * t64 - t15 * t63 - t159 * t200 - t44 * t35, t290 * t35 + t33 * t63, t28 * t15 - t199 * t35 + t7 * t44 + t5 * t63 + (t265 + t9 * t290 + (-t13 * t63 - t29 * t290 + t266) * qJD(6)) * t162 + t269 * t159, t28 * t14 - t4 * t35 - t7 * t196 + (-(-qJD(6) * t29 + t9) * t290 - t265 - (-qJD(6) * t13 + t6) * t63 - qJD(6) * t266) * t159 + t269 * t162; 0, 0, 0, -t222, t234 * t167, 0, 0, 0, t167 * pkin(1) * t161, pkin(1) * t243 (t113 * t164 - t161 * t98) * qJD(1) ((t123 - t154) * t161 + (-t121 + t207) * t164) * qJD(1), 0.2e1 * t153 + (t113 * t161 + t164 * t98) * qJD(1), t119 * qJ(3) + t123 * qJD(3) - t98 * t113 + (t123 * t161 + (-t121 - t256) * t164) * qJD(1) * pkin(7), -t261, t224, -t281, -t289, 0, t152 * t292 - t86 * t94 + t276, t152 * t291 - t86 * t96 - t179, -t68 * t33 - t67 * t34 + (-t17 + t260) * t187 + (t16 - t259) * t58, -t1 * t67 - t260 * t16 + t259 * t17 - t54 * t180 + t2 * t68, t295, -t300, -t297, t250 * t290 - t263 - t30, t278, t65 * t15 + t176 * t159 - t270 * t162 + t260 * t44 - t280, t65 * t14 + t270 * t159 + t176 * t162 - t196 * t260 - t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, 0, -t155 * t167 - t166, -t123 * qJD(2) + t98 * t233 + t135, 0, 0, 0, 0, 0, -t160 * t203 - t94 * t233, -t163 * t203 - t96 * t233, -t110 * t33 - t186 * t34 - t187 * t258 + t257 * t58, -t1 * t186 + t2 * t110 + t258 * t16 - t257 * t17 - t54 * t233, 0, 0, 0, 0, 0, -t110 * t252 - t186 * t15 - t258 * t44 + (t257 * t159 - t192 * t162) * t290, -t110 * t30 - t186 * t14 + t258 * t196 + (t192 * t159 + t257 * t162) * t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, -t224, t281, t289, 0, -t198 * t152 - t276, -t212 * t152 + t179 (-t157 * t33 - t244 * t34) * pkin(4) - (t16 - t19) * t58 + (t17 - t18) * t187, t16 * t18 - t17 * t19 + (-t244 * t1 + t157 * t2 - t54 * t96) * pkin(4), -t295, t300, t297, t193 + t263, -t278, t139 * t15 + t178 * t159 - t271 * t162 - t18 * t44 + t280, t139 * t14 + t271 * t159 + t178 * t162 + t18 * t196 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187 ^ 2 - t58 ^ 2, t16 * t187 + t17 * t58 + t168, 0, 0, 0, 0, 0, t193 - t263, t264 + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196 * t44, t196 ^ 2 - t44 ^ 2, t14 + t283, -t15 - t279, t33, t12 * t196 - t159 * t2 - t288 * t4 + t5, t12 * t44 - t159 * t6 - t162 * t2 + t199 * t288;];
tauc_reg  = t3;
