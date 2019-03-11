% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:01
% EndTime: 2019-03-09 03:13:10
% DurationCPUTime: 2.69s
% Computational Cost: add. (3231->355), mult. (7011->443), div. (0->0), fcn. (3912->6), ass. (0->209)
t117 = sin(pkin(9)) * pkin(1) + pkin(7);
t246 = pkin(4) + t117;
t134 = sin(qJ(5));
t136 = cos(qJ(5));
t205 = t136 * qJD(3);
t137 = cos(qJ(3));
t215 = qJD(1) * t137;
t96 = -t134 * t215 + t205;
t265 = qJD(5) * t96;
t103 = t117 * qJD(1);
t124 = t137 * qJD(2);
t135 = sin(qJ(3));
t75 = t135 * t103 - t124;
t264 = -qJD(4) - t75;
t208 = t135 * qJD(1);
t116 = qJD(5) + t208;
t138 = -pkin(3) - pkin(8);
t185 = pkin(4) * qJD(1) + t103;
t255 = t124 - qJD(4);
t219 = t185 * t135 - t255;
t38 = t138 * qJD(3) + t219;
t118 = -cos(pkin(9)) * pkin(1) - pkin(2);
t157 = -t135 * qJ(4) + t118;
t72 = t138 * t137 + t157;
t50 = t72 * qJD(1);
t15 = -t134 * t50 + t136 * t38;
t211 = qJD(5) * t136;
t212 = qJD(5) * t134;
t204 = qJD(1) * qJD(3);
t188 = t135 * t204;
t115 = pkin(3) * t188;
t172 = pkin(8) * t135 - qJ(4) * t137;
t206 = t135 * qJD(4);
t146 = t172 * qJD(3) - t206;
t43 = qJD(1) * t146 + t115;
t207 = t135 * qJD(2);
t44 = (t185 * t137 + t207) * qJD(3);
t152 = -t134 * t44 - t136 * t43 - t38 * t211 + t50 * t212;
t263 = -t15 * t116 - t152;
t16 = t134 * t38 + t136 * t50;
t186 = t134 * t43 - t136 * t44 + t50 * t211 + t38 * t212;
t262 = t16 * t116 - t186;
t163 = t96 * t116;
t209 = t134 * qJD(3);
t94 = t136 * t215 + t209;
t164 = t94 * t116;
t54 = qJD(5) * t94 - t134 * t188;
t187 = t136 * t204;
t55 = -t135 * t187 + t265;
t261 = (t54 + t164) * t134 - (t55 + t163) * t136;
t12 = t116 * qJ(6) + t16;
t119 = t137 * t204;
t182 = pkin(5) * t119;
t2 = -t182 + t186;
t260 = -t12 * t116 + t2;
t196 = t136 * t208;
t258 = (t196 + t211) * t116;
t76 = t137 * t103 + t207;
t61 = -qJD(3) * qJ(4) - t76;
t58 = -qJD(3) * pkin(3) - t264;
t87 = t246 * t135;
t242 = t134 * t87 + t136 * t72;
t214 = qJD(3) * t135;
t122 = pkin(3) * t214;
t59 = t122 + t146;
t213 = qJD(3) * t137;
t83 = t246 * t213;
t9 = -qJD(5) * t242 - t134 * t59 + t136 * t83;
t168 = t134 * t15 - t136 * t16;
t254 = -qJD(5) * t168 - t134 * t152 - t136 * t186;
t178 = qJ(6) * t119;
t1 = t116 * qJD(6) - t152 + t178;
t218 = qJD(6) - t15;
t11 = -t116 * pkin(5) + t218;
t170 = t11 * t134 + t12 * t136;
t253 = qJD(5) * t170 + t1 * t134 - t2 * t136;
t252 = t96 ^ 2;
t251 = t116 ^ 2;
t121 = pkin(4) * t215;
t48 = t121 - t61;
t18 = t94 * pkin(5) - t96 * qJ(6) + t48;
t250 = t18 * t96;
t128 = qJD(3) * qJD(4);
t64 = qJD(3) * t124 - t103 * t214;
t53 = -t128 - t64;
t32 = -pkin(4) * t188 - t53;
t6 = t55 * pkin(5) + t54 * qJ(6) - t96 * qJD(6) + t32;
t249 = t6 * t134;
t248 = t6 * t136;
t247 = t96 * t94;
t174 = pkin(5) * t136 + qJ(6) * t134;
t156 = -pkin(4) - t174;
t245 = (qJD(1) * t156 - t103) * t135 - t174 * qJD(5) + t136 * qJD(6) + t255;
t223 = t134 * t138;
t79 = t94 * t211;
t244 = -t138 * t79 - t55 * t223;
t45 = t135 * t55;
t243 = t94 * t213 + t45;
t57 = t121 + t76;
t123 = pkin(3) * t208;
t77 = t172 * qJD(1) + t123;
t23 = t134 * t57 + t136 * t77;
t46 = t135 * t54;
t241 = t96 * t213 - t46;
t194 = t135 * t205;
t210 = qJD(5) * t137;
t192 = t134 * t210;
t93 = t116 * t192;
t240 = t116 * t194 + t93;
t238 = t134 * t96;
t237 = t136 * t54;
t236 = t136 * t94;
t235 = t137 * t54;
t234 = t137 * t96;
t233 = t138 * t54;
t232 = t138 * t96;
t229 = t32 * t134;
t228 = t32 * t136;
t47 = t55 * t134;
t65 = t76 * qJD(3);
t227 = t65 * t135;
t226 = t65 * t137;
t86 = -t137 * pkin(3) + t157;
t62 = qJD(1) * t86;
t224 = t116 * t135;
t222 = t135 * t136;
t139 = qJD(3) ^ 2;
t125 = t139 * t135;
t126 = t139 * t137;
t221 = t18 * qJD(3);
t220 = t48 * qJD(3);
t88 = t246 * t137;
t130 = t135 ^ 2;
t131 = t137 ^ 2;
t217 = t130 - t131;
t104 = qJD(1) * t118;
t216 = qJD(1) * t131;
t203 = t94 * t196 + t47 + t79;
t202 = t94 ^ 2 - t252;
t201 = t94 * t214;
t200 = t96 * t214;
t199 = t116 * t223;
t198 = t116 * t136 * t138;
t197 = t134 * t216;
t195 = t135 * t209;
t193 = t138 * t213;
t191 = t136 * t210;
t190 = t116 * t215;
t108 = t134 * t119;
t183 = -qJD(3) * t96 - t108;
t181 = t96 * t194;
t180 = t116 * t191;
t179 = t135 * t119;
t110 = t136 * t119;
t177 = t75 * qJD(3) + t64;
t173 = -t134 * pkin(5) + t136 * qJ(6);
t171 = -t11 * t136 + t12 * t134;
t169 = t134 * t16 + t136 * t15;
t22 = -t134 * t77 + t136 * t57;
t28 = -t134 * t72 + t136 * t87;
t165 = -0.2e1 * qJD(3) * t62;
t162 = t216 - t224;
t161 = 0.2e1 * qJD(3) * t104;
t160 = t116 * t134;
t69 = t94 * t195;
t159 = -t69 - t181;
t158 = t116 * t195 - t180;
t150 = -qJ(4) * t213 - t206;
t63 = qJD(1) * t150 + t115;
t84 = t122 + t150;
t154 = qJD(1) * t84 + t117 * t139 + t63;
t8 = t134 * t83 + t136 * t59 + t87 * t211 - t72 * t212;
t151 = t94 * t215 - t108 - t258;
t149 = -t131 * t187 + t240;
t148 = t55 - t163;
t147 = t47 + (t236 - t238) * qJD(5);
t145 = -qJD(3) * t94 - t116 * t160 + t110;
t144 = t134 * t163 - t203 + t237;
t143 = t227 - t53 * t137 + (t135 * t61 + t137 * t58) * qJD(3);
t142 = t227 + t64 * t137 + (-t135 * t76 + t137 * t75) * qJD(3);
t141 = -t94 * t192 + (t137 * t55 - t201) * t136;
t140 = qJD(1) ^ 2;
t114 = t135 * t140 * t137;
t106 = -0.2e1 * t179;
t105 = 0.2e1 * t179;
t102 = t217 * t140;
t101 = qJ(4) - t173;
t100 = t138 * t110;
t98 = -qJ(4) * t215 + t123;
t85 = -0.2e1 * t217 * t204;
t82 = t246 * t214;
t67 = (t116 + t208) * t213;
t51 = t62 * t208;
t42 = t174 * t137 + t88;
t41 = pkin(5) * t96 + qJ(6) * t94;
t39 = t136 * t235;
t26 = t164 - t54;
t25 = -t135 * pkin(5) - t28;
t24 = t135 * qJ(6) + t242;
t21 = -t116 * t212 + t110 + (-t134 * t224 - t234) * qJD(1);
t20 = -pkin(5) * t215 - t22;
t19 = qJ(6) * t215 + t23;
t17 = (t173 * qJD(5) + qJD(6) * t134) * t137 + (-t117 + t156) * t214;
t14 = -t160 * t96 - t237;
t13 = -t96 * t191 + (t200 + t235) * t134;
t10 = -qJD(3) * t197 + t158 + t241;
t7 = -pkin(5) * t213 - t9;
t5 = qJ(6) * t213 + t135 * qJD(6) + t8;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t85, t126, t106, -t125, 0, -t117 * t126 + t135 * t161, t117 * t125 + t137 * t161, t142, t142 * t117, 0, -t126, t125, t105, t85, t106, t143, t135 * t165 + t137 * t154, -t135 * t154 + t137 * t165, t117 * t143 + t62 * t84 + t63 * t86, t13, t181 + t39 - t69 + (t47 + (t236 + t238) * qJD(5)) * t137, t10, t141, -t45 + (-t136 * t216 - t137 * t94) * qJD(3) + t240, t67, t9 * t116 + t88 * t55 - t82 * t94 + (-t205 * t48 - t186) * t135 + (-t48 * t212 + t228 + (qJD(1) * t28 + t15) * qJD(3)) * t137, -t8 * t116 - t88 * t54 - t82 * t96 + (t209 * t48 + t152) * t135 + (-t48 * t211 - t229 + (-qJD(1) * t242 - t16) * qJD(3)) * t137, t28 * t54 - t242 * t55 - t8 * t94 - t9 * t96 - t168 * t214 + (qJD(5) * t169 - t134 * t186 + t136 * t152) * t137, t15 * t9 - t152 * t242 + t16 * t8 - t186 * t28 + t32 * t88 - t48 * t82, t13, t10, -t39 + (-t210 * t94 - t200) * t136 + (t201 + (-t55 - t265) * t137) * t134, t67, -t149 + t243, t141, -t7 * t116 + t17 * t94 + t42 * t55 + (-t18 * t205 - t2) * t135 + (-t18 * t212 + t248 + (-qJD(1) * t25 - t11) * qJD(3)) * t137, -t24 * t55 - t25 * t54 - t5 * t94 + t7 * t96 + t170 * t214 + (qJD(5) * t171 - t1 * t136 - t134 * t2) * t137, t5 * t116 - t17 * t96 + t42 * t54 + (-t18 * t209 + t1) * t135 + (t18 * t211 + t249 + (qJD(1) * t24 + t12) * qJD(3)) * t137, t1 * t24 + t11 * t7 + t12 * t5 + t17 * t18 + t2 * t25 + t42 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t126, 0, t64 * t135 - t226 + (t135 * t75 + t137 * t76) * qJD(3), 0, 0, 0, 0, 0, 0, 0, t125, t126, -t53 * t135 - t226 + (t135 * t58 - t137 * t61) * qJD(3), 0, 0, 0, 0, 0, 0, t149 + t243, t162 * t209 + t180 + t241 (t147 - t237) * t137 + t159 (qJD(3) * t169 + t32) * t135 + (t220 - t254) * t137, 0, 0, 0, 0, 0, 0, -t162 * t205 + t243 + t93, t137 * t147 + t159 - t39, t46 + (-t197 - t234) * qJD(3) + t158 (qJD(3) * t171 + t6) * t135 + (t221 - t253) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t102, 0, t114, 0, 0, -t104 * t208, -t104 * t215 - t177, 0, 0, 0, 0, 0, -t114, t102, t114, 0, -t98 * t215 + t51, 0.2e1 * t128 + (t135 * t98 + t137 * t62) * qJD(1) + t177, -pkin(3) * t65 - qJ(4) * t53 + t264 * t61 - t58 * t76 - t62 * t98, t14, t261, t21, t203, t151, -t190, qJ(4) * t55 - t22 * t116 + t229 + t100 + t219 * t94 + (t136 * t48 - t199) * qJD(5) + (-t137 * t15 + t48 * t222) * qJD(1), -qJ(4) * t54 + t23 * t116 + t228 + t219 * t96 + (-t134 * t48 - t198) * qJD(5) + (t137 * t16 + (-t135 * t48 - t193) * t134) * qJD(1), t22 * t96 + t23 * t94 + (t233 - t262) * t136 + (t15 * t208 + t152 + (t15 + t232) * qJD(5)) * t134 + t244, t32 * qJ(4) + t254 * t138 - t15 * t22 - t16 * t23 + t219 * t48, t14, t21, -t261, -t190, -t151, t136 * t164 + t47, t101 * t55 + t20 * t116 + t249 + t100 - t245 * t94 + (t136 * t18 - t199) * qJD(5) + (t11 * t137 + t18 * t222) * qJD(1), t19 * t94 - t20 * t96 + (t233 + t260) * t136 + (-t11 * t208 - t1 + (-t11 + t232) * qJD(5)) * t134 + t244, t101 * t54 - t19 * t116 - t248 + t245 * t96 + (t134 * t18 + t198) * qJD(5) + (-t12 * t137 + (t135 * t18 + t193) * t134) * qJD(1), t6 * t101 - t11 * t20 - t12 * t19 + t253 * t138 - t245 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t130 * t140 - t139, t51 + (t61 + t76) * qJD(3), 0, 0, 0, 0, 0, 0, t145, -t251 * t136 + t183, t144, t134 * t263 + t136 * t262 - t220, 0, 0, 0, 0, 0, 0, t145, t144, -t183 + t258, -t221 - t260 * t136 + (t11 * t116 + t1) * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, -t202, t26, -t247, -t148, t119, -t48 * t96 + t262, t48 * t94 - t263, 0, 0, t247, t26, t202, t119, t148, -t247, -t41 * t94 + 0.2e1 * t182 - t250 + t262, pkin(5) * t54 - qJ(6) * t55 + (t12 - t16) * t96 + (t11 - t218) * t94, 0.2e1 * t178 - t18 * t94 + t41 * t96 + (0.2e1 * qJD(6) - t15) * t116 - t152, -pkin(5) * t2 + qJ(6) * t1 - t11 * t16 + t12 * t218 - t18 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119 + t247, t26, -t251 - t252, t250 + t260;];
tauc_reg  = t3;
