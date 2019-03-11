% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:39:03
% EndTime: 2019-03-09 03:39:08
% DurationCPUTime: 2.62s
% Computational Cost: add. (3737->312), mult. (9086->435), div. (0->0), fcn. (6778->10), ass. (0->175)
t145 = sin(pkin(11));
t147 = cos(pkin(11));
t151 = sin(qJ(3));
t199 = t151 * qJD(1);
t154 = cos(qJ(3));
t204 = qJD(1) * t154;
t115 = -t145 * t199 + t147 * t204;
t238 = qJD(5) + qJD(6);
t248 = t115 - t238;
t149 = sin(qJ(6));
t150 = sin(qJ(5));
t152 = cos(qJ(6));
t153 = cos(qJ(5));
t127 = t149 * t153 + t152 * t150;
t228 = t248 * t127;
t125 = t145 * t154 + t147 * t151;
t117 = t125 * qJD(1);
t197 = t153 * qJD(3);
t96 = t150 * t117 - t197;
t98 = t150 * qJD(3) + t153 * t117;
t172 = t149 * t96 - t152 * t98;
t40 = t149 * t98 + t152 * t96;
t247 = t172 * t40;
t136 = sin(pkin(10)) * pkin(1) + pkin(7);
t206 = qJ(4) + t136;
t203 = qJD(5) * t150;
t211 = t150 * t115;
t246 = t203 - t211;
t245 = t172 ^ 2 - t40 ^ 2;
t181 = t206 * qJD(1);
t100 = t151 * qJD(2) + t181 * t154;
t214 = t147 * t100;
t227 = qJD(3) * pkin(3);
t99 = t154 * qJD(2) - t181 * t151;
t94 = t99 + t227;
t37 = t145 * t94 + t214;
t34 = qJD(3) * pkin(8) + t37;
t138 = -cos(pkin(10)) * pkin(1) - pkin(2);
t170 = -t154 * pkin(3) + t138;
t162 = t170 * qJD(1);
t114 = qJD(4) + t162;
t55 = -t115 * pkin(4) - t117 * pkin(8) + t114;
t22 = t150 * t55 + t153 * t34;
t14 = -t96 * pkin(9) + t22;
t201 = qJD(6) * t149;
t12 = t14 * t201;
t91 = t145 * t100;
t36 = t147 * t94 - t91;
t33 = -qJD(3) * pkin(4) - t36;
t26 = t96 * pkin(5) + t33;
t244 = t26 * t40 + t12;
t200 = qJD(6) * t152;
t196 = qJD(1) * qJD(3);
t187 = t154 * t196;
t188 = t151 * t196;
t109 = -t145 * t188 + t147 * t187;
t58 = qJD(5) * t197 + t153 * t109 - t117 * t203;
t59 = t98 * qJD(5) + t150 * t109;
t10 = -t149 * t59 + t152 * t58 - t96 * t200 - t98 * t201;
t111 = qJD(5) - t115;
t107 = qJD(6) + t111;
t243 = t40 * t107 + t10;
t116 = t125 * qJD(3);
t108 = qJD(1) * t116;
t195 = qJD(1) * qJD(4);
t240 = -t100 * qJD(3) - t151 * t195;
t82 = t99 * qJD(3) + t154 * t195;
t29 = t240 * t145 + t147 * t82;
t133 = pkin(3) * t188;
t61 = t108 * pkin(4) - t109 * pkin(8) + t133;
t57 = t153 * t61;
t159 = -t22 * qJD(5) - t150 * t29 + t57;
t2 = t108 * pkin(5) - t58 * pkin(9) + t159;
t202 = qJD(5) * t153;
t167 = t150 * t61 + t153 * t29 + t55 * t202 - t34 * t203;
t5 = -t59 * pkin(9) + t167;
t192 = -t149 * t5 + t152 * t2;
t224 = t152 * t14;
t21 = -t150 * t34 + t153 * t55;
t13 = -t98 * pkin(9) + t21;
t9 = t111 * pkin(5) + t13;
t4 = t149 * t9 + t224;
t242 = -t4 * qJD(6) + t26 * t172 + t192;
t158 = t172 * qJD(6) - t149 * t58 - t152 * t59;
t241 = -t107 * t172 + t158;
t76 = t127 * t125;
t124 = t145 * t151 - t147 * t154;
t119 = t124 * qJD(3);
t209 = t153 * t119;
t239 = -t125 * t203 - t209;
t126 = t149 * t150 - t152 * t153;
t229 = t248 * t126;
t237 = -t229 * t107 - t127 * t108;
t103 = t153 * t108;
t236 = t125 * t103 + t111 * t239;
t235 = t10 * t124 - t116 * t172;
t135 = t145 * pkin(3) + pkin(8);
t234 = pkin(9) + t135;
t210 = t150 * t119;
t216 = t125 * t153;
t217 = t125 * t150;
t19 = -t201 * t217 + (t238 * t216 - t210) * t152 + t239 * t149;
t233 = -t19 * t107 - t76 * t108;
t232 = t98 * t116 + t58 * t124;
t45 = t147 * t99 - t91;
t71 = pkin(3) * t199 + t117 * pkin(4) - t115 * pkin(8);
t231 = t150 * t71 + t153 * t45;
t120 = t206 * t151;
t121 = t206 * t154;
t80 = -t145 * t120 + t147 * t121;
t74 = t153 * t80;
t75 = t124 * pkin(4) - t125 * pkin(8) + t170;
t230 = t150 * t75 + t74;
t226 = t117 * t40;
t225 = t117 * t96;
t28 = t145 * t82 - t147 * t240;
t223 = t28 * t153;
t222 = t172 * t117;
t221 = t58 * t150;
t220 = t96 * t111;
t219 = t98 * t111;
t218 = t98 * t117;
t213 = t150 * t108;
t155 = qJD(3) ^ 2;
t208 = t155 * t151;
t207 = t155 * t154;
t205 = t151 ^ 2 - t154 ^ 2;
t130 = qJD(1) * t138;
t194 = t151 * t227;
t137 = -t147 * pkin(3) - pkin(4);
t190 = t125 * t202;
t189 = qJD(6) * t9 + t5;
t44 = t145 * t99 + t214;
t184 = qJD(5) * t234;
t183 = qJD(3) * t206;
t101 = t154 * qJD(4) - t151 * t183;
t102 = -t151 * qJD(4) - t154 * t183;
t53 = t145 * t101 - t147 * t102;
t79 = t147 * t120 + t145 * t121;
t182 = t111 * t153;
t180 = t228 * t107 - t126 * t108;
t179 = t246 * pkin(5) - t44;
t122 = t234 * t150;
t178 = -pkin(9) * t211 + qJD(6) * t122 + t150 * t184 + t231;
t123 = t234 * t153;
t64 = t153 * t71;
t177 = t117 * pkin(5) + qJD(6) * t123 - t150 * t45 + t64 + (-pkin(9) * t115 + t184) * t153;
t18 = t126 * t119 - t238 * t76;
t77 = t126 * t125;
t176 = -t18 * t107 + t77 * t108;
t175 = -t80 * t108 + t28 * t125;
t174 = -t116 * t40 + t124 * t158;
t173 = -t116 * t96 - t124 * t59;
t169 = 0.2e1 * qJD(3) * t130;
t168 = -t246 * t111 + t103;
t54 = t147 * t101 + t145 * t102;
t72 = t116 * pkin(4) + t119 * pkin(8) + t194;
t166 = t150 * t72 + t153 * t54 + t75 * t202 - t80 * t203;
t165 = t190 - t210;
t163 = -t135 * t108 + t111 * t33;
t157 = -t111 * t165 - t125 * t213;
t156 = qJD(1) ^ 2;
t128 = -t153 * pkin(5) + t137;
t83 = t108 * t124;
t70 = t153 * t75;
t65 = t153 * t72;
t49 = pkin(5) * t217 + t79;
t25 = pkin(5) * t165 + t53;
t24 = -pkin(9) * t217 + t230;
t23 = t124 * pkin(5) - pkin(9) * t216 - t150 * t80 + t70;
t17 = t59 * pkin(5) + t28;
t7 = -pkin(9) * t165 + t166;
t6 = pkin(9) * t209 + t116 * pkin(5) - t150 * t54 + t65 + (-t74 + (pkin(9) * t125 - t75) * t150) * qJD(5);
t3 = -t149 * t14 + t152 * t9;
t1 = [0, 0, 0, 0, 0.2e1 * t151 * t187, -0.2e1 * t205 * t196, t207, -t208, 0, -t136 * t207 + t151 * t169, t136 * t208 + t154 * t169, t79 * t109 + t54 * t115 - t37 * t116 + t53 * t117 + t36 * t119 - t29 * t124 + t175, t28 * t79 + t29 * t80 - t36 * t53 + t37 * t54 + (t114 + t162) * t194, -t98 * t209 + (t58 * t153 - t98 * t203) * t125 -(-t150 * t98 - t153 * t96) * t119 + (-t221 - t153 * t59 + (t150 * t96 - t153 * t98) * qJD(5)) * t125, t232 + t236, t157 + t173, t111 * t116 + t83 (-t202 * t80 + t65) * t111 + t70 * t108 + (-t202 * t34 + t57) * t124 + t21 * t116 + t53 * t96 + t79 * t59 + t33 * t190 + ((-qJD(5) * t75 - t54) * t111 + (-qJD(5) * t55 - t29) * t124 - t33 * t119 + t175) * t150, -t166 * t111 - t230 * t108 - t167 * t124 - t22 * t116 + t53 * t98 + t79 * t58 - t33 * t209 + (-t203 * t33 + t223) * t125, -t10 * t77 - t172 * t18, -t10 * t76 - t158 * t77 + t172 * t19 - t18 * t40, -t176 + t235, t174 + t233, t107 * t116 + t83 (-t149 * t7 + t152 * t6) * t107 + (-t149 * t24 + t152 * t23) * t108 + t192 * t124 + t3 * t116 + t25 * t40 - t49 * t158 + t17 * t76 + t26 * t19 + ((-t149 * t23 - t152 * t24) * t107 - t4 * t124) * qJD(6), t49 * t10 - t4 * t116 + t12 * t124 - t17 * t77 + t26 * t18 - t25 * t172 + (-(-qJD(6) * t24 + t6) * t107 - t23 * t108 - t2 * t124) * t149 + (-(qJD(6) * t23 + t7) * t107 - t24 * t108 - t189 * t124) * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, -t207, -t125 * t108 + t124 * t109 - t119 * t115 + t116 * t117, -t36 * t116 - t37 * t119 + t28 * t124 + t29 * t125, 0, 0, 0, 0, 0, t157 - t173, t232 - t236, 0, 0, 0, 0, 0, -t174 + t233, t176 + t235; 0, 0, 0, 0, -t151 * t156 * t154, t205 * t156, 0, 0, 0, -t130 * t199, -t130 * t204 (t37 - t44) * t117 + (t36 - t45) * t115 + (-t108 * t145 - t109 * t147) * pkin(3), t36 * t44 - t37 * t45 + (-t114 * t199 + t145 * t29 - t147 * t28) * pkin(3), t182 * t98 + t221 (t58 - t220) * t153 + (-t59 - t219) * t150, t111 * t182 + t213 - t218, t168 + t225, -t111 * t117, -t21 * t117 + t137 * t59 - t223 - t44 * t96 + (-t135 * t202 - t64) * t111 + (t45 * t111 + t163) * t150, t22 * t117 + t137 * t58 + t28 * t150 - t44 * t98 + (t135 * t203 + t231) * t111 + t163 * t153, t10 * t127 - t172 * t229, -t10 * t126 + t127 * t158 - t172 * t228 - t229 * t40, t222 - t237, t180 + t226, -t107 * t117 (-t152 * t122 - t149 * t123) * t108 - t128 * t158 + t17 * t126 - t3 * t117 + t179 * t40 - t228 * t26 + (t149 * t178 - t152 * t177) * t107 -(-t149 * t122 + t152 * t123) * t108 + t128 * t10 + t17 * t127 + t4 * t117 - t179 * t172 + t229 * t26 + (t149 * t177 + t152 * t178) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 ^ 2 - t117 ^ 2, -t37 * t115 + t36 * t117 + t133, 0, 0, 0, 0, 0, t168 - t225, -t111 ^ 2 * t153 - t213 - t218, 0, 0, 0, 0, 0, t180 - t226, t222 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t96, -t96 ^ 2 + t98 ^ 2, t58 + t220, t219 - t59, t108, t22 * t111 - t33 * t98 + t159, t21 * t111 + t33 * t96 - t167, -t247, t245, t243, t241, t108 -(-t149 * t13 - t224) * t107 + (-t107 * t201 + t152 * t108 - t98 * t40) * pkin(5) + t242 (-t14 * t107 - t2) * t149 + (t13 * t107 - t189) * t152 + (-t107 * t200 - t149 * t108 + t172 * t98) * pkin(5) + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t247, t245, t243, t241, t108, t4 * t107 + t242, t3 * t107 - t149 * t2 - t152 * t189 + t244;];
tauc_reg  = t1;
