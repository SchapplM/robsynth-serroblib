% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:24
% EndTime: 2019-03-09 02:16:31
% DurationCPUTime: 3.07s
% Computational Cost: add. (5322->358), mult. (11640->437), div. (0->0), fcn. (8111->6), ass. (0->181)
t124 = sin(qJ(5));
t125 = cos(qJ(5));
t187 = qJD(5) * t125;
t188 = qJD(5) * t124;
t121 = sin(pkin(9));
t122 = cos(pkin(9));
t126 = cos(qJ(4));
t221 = sin(qJ(4));
t96 = t121 * t126 + t122 * t221;
t134 = t96 * qJD(3);
t123 = -pkin(1) - qJ(3);
t235 = t123 * qJD(1);
t106 = qJD(2) + t235;
t172 = -pkin(7) * qJD(1) + t106;
t83 = t172 * t121;
t84 = t172 * t122;
t236 = t126 * t84 - t221 * t83;
t35 = -qJD(1) * t134 + qJD(4) * t236;
t57 = t126 * t83 + t221 * t84;
t53 = qJD(4) * pkin(8) + t57;
t120 = qJD(1) * qJ(2);
t114 = qJD(3) + t120;
t116 = t121 * pkin(3);
t100 = qJD(1) * t116 + t114;
t133 = qJD(1) * t96;
t195 = t122 * t126;
t178 = qJD(1) * t195;
t179 = t221 * t121;
t91 = -qJD(1) * t179 + t178;
t54 = pkin(4) * t133 - pkin(8) * t91 + t100;
t119 = qJD(1) * qJD(2);
t176 = qJD(4) * t221;
t164 = qJD(1) * t176;
t189 = qJD(4) * t126;
t175 = qJD(1) * t189;
t191 = t121 * t175 + t122 * t164;
t103 = t121 * t164;
t82 = t122 * t175 - t103;
t55 = pkin(4) * t82 + pkin(8) * t191 + t119;
t173 = t124 * t35 - t125 * t55 + t187 * t53 + t188 * t54;
t20 = t124 * t54 + t125 * t53;
t242 = qJD(5) + t133;
t220 = t20 * t242;
t232 = -t173 + t220;
t97 = -t179 + t195;
t174 = t242 ^ 2;
t16 = qJ(6) * t242 + t20;
t246 = t16 * t242;
t186 = t125 * qJD(4);
t73 = t124 * t91 - t186;
t245 = t242 * t73;
t19 = -t124 * t53 + t125 * t54;
t3 = t124 * t55 + t125 * t35 + t187 * t54 - t188 * t53;
t244 = -t19 * t242 + t3;
t243 = t97 * qJD(3);
t45 = -qJD(5) * t186 + t125 * t191 + t188 * t91;
t75 = qJD(4) * t124 + t125 * t91;
t93 = t96 * qJD(4);
t157 = -t97 * t45 - t93 * t75;
t168 = qJD(5) * t96 + qJD(1);
t94 = -t121 * t176 + t122 * t189;
t201 = t125 * t94;
t80 = t125 * t82;
t241 = (t124 * t168 - t201) * t242 - t96 * t80 - t157;
t46 = qJD(5) * t75 - t124 * t191;
t198 = t46 * t125;
t199 = t45 * t124;
t202 = t125 * t75;
t240 = (qJD(5) * (t124 * t73 - t202) - t198 + t199) * t97 + (t124 * t75 + t125 * t73) * t93;
t239 = t46 * t97 - t73 * t93;
t207 = -pkin(7) + t123;
t98 = t207 * t121;
t99 = t207 * t122;
t69 = -t126 * t99 + t221 * t98;
t237 = -t242 * t75 + t46;
t190 = t121 ^ 2 + t122 ^ 2;
t234 = t190 * qJD(3);
t224 = pkin(5) * t82;
t2 = t173 - t224;
t233 = -t2 + t246;
t231 = -t191 * t97 - t91 * t93;
t52 = -qJD(4) * pkin(4) - t236;
t21 = pkin(5) * t73 - qJ(6) * t75 + t52;
t222 = pkin(8) * t82;
t230 = t21 * t242 - t222;
t208 = t97 * t82;
t155 = -t242 * t93 + t208;
t180 = t97 * t187;
t229 = t124 * t155 + t180 * t242 + t96 * t46 + t94 * t73;
t110 = qJ(2) + t116;
t66 = pkin(4) * t96 - pkin(8) * t97 + t110;
t70 = t126 * t98 + t221 * t99;
t204 = t124 * t66 + t125 * t70;
t47 = -qJD(4) * t69 - t134;
t64 = pkin(4) * t94 + pkin(8) * t93 + qJD(2);
t10 = -qJD(5) * t204 - t124 * t47 + t125 * t64;
t226 = t75 ^ 2;
t225 = t91 ^ 2;
t184 = 0.2e1 * t119;
t223 = pkin(8) * t75;
t141 = t243 * qJD(1);
t36 = t57 * qJD(4) + t141;
t219 = t36 * t69;
t218 = t36 * t97;
t217 = t73 * t133;
t216 = t75 * t73;
t214 = t75 * t91;
t213 = t82 * t96;
t212 = t242 * t91;
t211 = t91 * t73;
t210 = t91 * t133;
t152 = pkin(5) * t124 - qJ(6) * t125;
t206 = t124 * qJD(6) - t152 * t242 + t57;
t205 = -t124 * t46 - t187 * t73;
t67 = pkin(4) * t91 + pkin(8) * t133;
t26 = t124 * t67 + t125 * t236;
t78 = t124 * t82;
t203 = t124 * t133;
t197 = t82 * qJ(6);
t196 = qJD(4) * t133;
t194 = t94 * qJD(4);
t192 = qJD(6) - t19;
t185 = pkin(8) * qJD(5) * t242;
t183 = t73 ^ 2 - t226;
t181 = t97 * t188;
t71 = t75 * t188;
t171 = t124 * t242;
t169 = qJD(1) * t190;
t7 = t46 * pkin(5) + t45 * qJ(6) - t75 * qJD(6) + t36;
t167 = -t7 - t185;
t166 = t36 + t185;
t165 = t203 * t75 + t71;
t1 = qJD(6) * t242 + t197 + t3;
t163 = t1 * t96 + t16 * t94;
t15 = -pkin(5) * t242 + t192;
t162 = t15 * t94 + t2 * t96;
t161 = t173 * t96 - t19 * t94;
t160 = t20 * t94 + t3 * t96;
t159 = t21 * t93 - t7 * t97;
t158 = -t52 * t93 + t218;
t156 = -t45 * t96 + t75 * t94;
t154 = t133 * t94 + t213;
t153 = (qJD(5) * t73 - t45) * pkin(8);
t151 = t124 * t16 - t125 * t15;
t150 = t124 * t20 + t125 * t19;
t25 = -t124 * t236 + t125 * t67;
t28 = -t124 * t70 + t125 * t66;
t144 = t78 + (t125 * t133 + t187) * t242;
t143 = t80 + (-t188 - t203) * t242;
t142 = t21 * t75 + t173;
t139 = t242 * t52 - t222;
t9 = t124 * t64 + t125 * t47 + t187 * t66 - t188 * t70;
t138 = (t45 - t217) * t125 + t205;
t137 = t171 * t73 - t198;
t136 = -t143 - t211;
t131 = t236 * t93 - t35 * t96 - t57 * t94 + t218;
t130 = t124 * t239 + t180 * t73;
t48 = qJD(4) * t70 + t243;
t129 = -t96 * t78 + (-t124 * t94 - t125 * t168) * t242 - t239;
t128 = -t96 * t198 - t73 * t201 + t168 * t202 + (t168 * t73 + t156) * t124;
t127 = qJD(1) ^ 2;
t101 = -pkin(5) * t125 - qJ(6) * t124 - pkin(4);
t88 = t133 ^ 2;
t85 = t93 * qJD(4);
t39 = pkin(5) * t75 + qJ(6) * t73;
t38 = pkin(8) * t198;
t34 = t152 * t97 + t69;
t32 = t242 * t94 + t213;
t24 = -pkin(5) * t96 - t28;
t23 = qJ(6) * t96 + t204;
t22 = -t45 + t245;
t18 = -pkin(5) * t91 - t25;
t17 = qJ(6) * t91 + t26;
t14 = t144 - t214;
t13 = t202 * t242 - t199;
t12 = (qJ(6) * qJD(5) * t97 - pkin(5) * t93) * t124 + (qJ(6) * t93 + (pkin(5) * qJD(5) - qJD(6)) * t97) * t125 + t48;
t11 = t125 * t157 - t71 * t97;
t8 = -t94 * pkin(5) - t10;
t6 = qJ(6) * t94 + qJD(6) * t96 + t9;
t5 = t125 * t155 - t181 * t242 + t156;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, qJ(2) * t184, 0, 0, 0, 0, 0, 0, t121 * t184, t122 * t184, 0.2e1 * qJD(3) * t169 (t114 + t120) * qJD(2) + (-t106 - t235) * t234, t231, t133 * t93 + t191 * t96 - t91 * t94 - t208, -t85, t154, -t194, 0, 0.2e1 * qJD(2) * t133 - t48 * qJD(4) + t100 * t94 + t110 * t82, -t110 * t191 - t100 * t93 - t47 * qJD(4) + (qJD(1) * t97 + t91) * qJD(2), -t133 * t47 - t191 * t69 + t48 * t91 - t70 * t82 + t131, t35 * t70 + t219 + t57 * t47 - t236 * t48 + (qJD(1) * t110 + t100) * qJD(2), t11, t240, t5, t130, -t229, t32, t10 * t242 + t124 * t158 + t180 * t52 + t28 * t82 + t69 * t46 + t48 * t73 - t161, t125 * t158 - t181 * t52 - t204 * t82 - t242 * t9 - t69 * t45 + t48 * t75 - t160, -t10 * t75 + t28 * t45 - t204 * t46 - t9 * t73 + t150 * t93 + (-t3 * t124 + t173 * t125 + (t124 * t19 - t125 * t20) * qJD(5)) * t97, t10 * t19 - t173 * t28 + t20 * t9 + t204 * t3 + t48 * t52 + t219, t11, t5, -t240, t32, t229, t130, t12 * t73 - t124 * t159 + t180 * t21 - t24 * t82 - t242 * t8 + t34 * t46 - t162, -t23 * t46 - t24 * t45 - t6 * t73 + t8 * t75 + t151 * t93 + (-t1 * t124 + t2 * t125 + (-t124 * t15 - t125 * t16) * qJD(5)) * t97, -t12 * t75 + t125 * t159 + t181 * t21 + t23 * t82 + t242 * t6 + t34 * t45 + t163, t1 * t23 + t12 * t21 + t15 * t8 + t16 * t6 + t2 * t24 + t34 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t127 * qJ(2), 0, 0, 0, 0, 0, 0, -t127 * t121, -t127 * t122, 0 (-t114 - t234) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t133 - t85, -qJD(1) * t91 - t194, -t154 - t231, -qJD(1) * t100 - t131, 0, 0, 0, 0, 0, 0, t129, t241, t128 (-t168 * t19 + t160) * t125 + (-t168 * t20 + t161) * t124 - t158, 0, 0, 0, 0, 0, 0, t129, t128, -t241 (t15 * t168 + t163) * t125 + (-t16 * t168 + t162) * t124 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190 * t127, t106 * t169 + t119, 0, 0, 0, 0, 0, 0, -t103 + (t91 + t178) * qJD(4), -t191 - t196, -t88 - t225, t133 * t57 + t236 * t91 + t119, 0, 0, 0, 0, 0, 0, t143 - t211, -t125 * t174 - t214 - t78, t171 * t75 + t138, t124 * t244 + t125 * t232 - t52 * t91, 0, 0, 0, 0, 0, 0, -t124 * t174 - t211 + t80, t138 + t165, t144 + t214, -t21 * t91 + t233 * t125 + (t15 * t242 + t1) * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, -t88 + t225, -t191 + t196, -t210, t103 + (t91 - t178) * qJD(4), 0, -t100 * t91 - t141 (qJD(3) + t100) * t133, 0, 0, t13 (-t45 - t217) * t125 - t165 + t205, t14, t137, -t136, -t212, -pkin(4) * t46 + t124 * t139 - t125 * t166 - t19 * t91 - t242 * t25 - t57 * t73, pkin(4) * t45 + t124 * t166 + t125 * t139 + t20 * t91 + t242 * t26 - t57 * t75, t25 * t75 + t26 * t73 - t38 + (-t19 * t133 + t3 + (-t19 + t223) * qJD(5)) * t125 + (t153 - t232) * t124, -t36 * pkin(4) - t19 * t25 - t20 * t26 - t52 * t57 + (-qJD(5) * t150 + t124 * t173 + t3 * t125) * pkin(8), t13, t14, t71 + (t133 * t75 + t46) * t124 + (t45 + t245) * t125, -t212, t136, t137, t101 * t46 + t124 * t230 + t125 * t167 + t15 * t91 + t18 * t242 - t206 * t73, t17 * t73 - t18 * t75 - t38 + (t15 * t133 + t1 + (t15 + t223) * qJD(5)) * t125 + (t153 - t233) * t124, t101 * t45 + t124 * t167 - t125 * t230 - t16 * t91 - t17 * t242 + t206 * t75, t7 * t101 - t15 * t18 - t16 * t17 - t206 * t21 + (-qJD(5) * t151 + t1 * t125 + t2 * t124) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, -t183, t22, -t216, -t237, t82, -t52 * t75 + t232, t52 * t73 - t244, 0, 0, t216, t22, t183, t82, t237, -t216, -t39 * t73 - t142 + t220 + 0.2e1 * t224, pkin(5) * t45 - t46 * qJ(6) + (t16 - t20) * t75 + (t15 - t192) * t73, 0.2e1 * t197 - t21 * t73 + t39 * t75 + (0.2e1 * qJD(6) - t19) * t242 + t3, -t2 * pkin(5) + t1 * qJ(6) - t15 * t20 + t16 * t192 - t21 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82 + t216, t22, -t226 - t174, t142 - t224 - t246;];
tauc_reg  = t4;
