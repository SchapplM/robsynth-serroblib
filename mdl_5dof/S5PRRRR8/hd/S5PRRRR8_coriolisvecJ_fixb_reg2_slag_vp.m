% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:41
% EndTime: 2019-12-05 17:16:51
% DurationCPUTime: 2.89s
% Computational Cost: add. (4387->347), mult. (11215->490), div. (0->0), fcn. (8432->10), ass. (0->197)
t121 = sin(qJ(2));
t116 = sin(pkin(5));
t192 = qJD(1) * t116;
t167 = t121 * t192;
t120 = sin(qJ(3));
t186 = qJD(3) * t120;
t136 = pkin(3) * t186 - t167;
t123 = cos(qJ(3));
t232 = cos(qJ(4));
t170 = t232 * t123;
t119 = sin(qJ(4));
t197 = t119 * t120;
t137 = t170 - t197;
t235 = -pkin(8) - pkin(7);
t100 = t235 * t123;
t99 = t235 * t120;
t138 = t119 * t100 + t232 * t99;
t124 = cos(qJ(2));
t166 = t124 * t192;
t171 = qJD(3) * t235;
t196 = t119 * t123;
t95 = t120 * t171;
t220 = t138 * qJD(4) - t137 * t166 + t171 * t196 + t232 * t95;
t152 = qJD(3) * t170;
t163 = qJD(4) * t232;
t181 = qJD(3) + qJD(4);
t67 = -t123 * t163 + t181 * t197 - t152;
t93 = t232 * t120 + t196;
t68 = t181 * t93;
t238 = -t68 * pkin(4) - t67 * pkin(9) - t136;
t79 = -t232 * t100 + t119 * t99;
t219 = t79 * qJD(4) + t119 * t95 - t235 * t152 - t93 * t166;
t118 = sin(qJ(5));
t122 = cos(qJ(5));
t97 = qJD(2) * pkin(7) + t167;
t160 = pkin(8) * qJD(2) + t97;
t117 = cos(pkin(5));
t191 = qJD(1) * t120;
t165 = t117 * t191;
t72 = t160 * t123 + t165;
t176 = t232 * t72;
t198 = t117 * t123;
t104 = qJD(1) * t198;
t145 = t160 * t120;
t71 = t104 - t145;
t66 = qJD(3) * pkin(3) + t71;
t36 = t119 * t66 + t176;
t31 = t181 * pkin(9) + t36;
t113 = -t123 * pkin(3) - pkin(2);
t83 = t113 * qJD(2) - t166;
t189 = qJD(2) * t120;
t88 = -qJD(2) * t170 + t119 * t189;
t90 = t93 * qJD(2);
t47 = t88 * pkin(4) - t90 * pkin(9) + t83;
t142 = t118 * t31 - t122 * t47;
t19 = t118 * t47 + t122 * t31;
t237 = t118 * t142 + t122 * t19;
t236 = t181 * qJD(2);
t59 = t137 * t236;
t77 = t118 * t181 + t122 * t90;
t34 = t77 * qJD(5) + t118 * t59;
t190 = qJD(2) * t116;
t162 = qJD(1) * t190;
t151 = t124 * t162;
t141 = t120 * t151;
t185 = qJD(4) * t119;
t202 = qJD(3) * t104 + t123 * t151;
t50 = -qJD(3) * t145 + t202;
t130 = -t119 * t141 - t72 * t185 + t232 * t50;
t132 = t72 * qJD(3);
t10 = -t119 * t132 + t66 * t163 + t130;
t60 = t68 * qJD(2);
t182 = qJD(2) * qJD(3);
t161 = t120 * t182;
t85 = pkin(3) * t161 + t121 * t162;
t27 = t60 * pkin(4) - t59 * pkin(9) + t85;
t2 = -t142 * qJD(5) + t122 * t10 + t118 * t27;
t62 = -pkin(4) * t137 - t93 * pkin(9) + t113;
t37 = -t118 * t79 + t122 * t62;
t234 = t37 * qJD(5) - t238 * t118 + t220 * t122;
t38 = t118 * t62 + t122 * t79;
t233 = t38 * qJD(5) + t220 * t118 + t238 * t122;
t49 = t232 * (-t132 - t141);
t159 = t119 * t50 - t49;
t11 = t36 * qJD(4) + t159;
t201 = t116 * t121;
t86 = -t120 * t201 + t198;
t87 = t117 * t120 + t123 * t201;
t140 = -t119 * t87 + t232 * t86;
t231 = t11 * t140;
t230 = t11 * t138;
t229 = t11 * t93;
t1 = t2 * t122;
t228 = t60 * t137;
t157 = t122 * t181;
t75 = t118 * t90 - t157;
t84 = qJD(5) + t88;
t227 = t75 * t84;
t226 = t77 * t75;
t225 = t77 * t84;
t224 = t83 * t90;
t223 = t84 * t90;
t222 = t90 * t88;
t221 = t93 * t60;
t218 = pkin(3) * qJD(4);
t217 = qJD(2) * pkin(2);
t214 = t118 * t60;
t213 = t118 * t75;
t212 = t118 * t88;
t211 = t119 * t72;
t210 = t120 * t97;
t208 = t122 * t60;
t207 = t122 * t77;
t206 = t122 * t88;
t205 = t123 * t97;
t184 = qJD(5) * t118;
t33 = -qJD(5) * t157 - t122 * t59 + t90 * t184;
t204 = t33 * t118;
t203 = t34 * t122;
t200 = t116 * t124;
t126 = qJD(2) ^ 2;
t199 = t116 * t126;
t125 = qJD(3) ^ 2;
t195 = t125 * t120;
t194 = t125 * t123;
t114 = t120 ^ 2;
t115 = t123 ^ 2;
t193 = t114 - t115;
t188 = qJD(2) * t121;
t187 = qJD(2) * t123;
t183 = qJD(5) * t122;
t180 = t142 * t206 - t19 * t212 + t1;
t178 = pkin(3) * t189;
t177 = t232 * t66;
t175 = t93 * t184;
t174 = t93 * t183;
t173 = t121 * t199;
t172 = t120 * t126 * t123;
t35 = t177 - t211;
t30 = -t181 * pkin(4) - t35;
t28 = t30 * t184;
t29 = t30 * t183;
t169 = t116 * t188;
t168 = t124 * t190;
t164 = t142 * t90 + t28;
t158 = t122 * t84;
t156 = t11 * t118 + t19 * t90 + t29;
t155 = pkin(3) * t163;
t154 = t120 * t168;
t153 = t123 * t168;
t150 = t123 * t161;
t39 = t119 * t71 + t176;
t149 = pkin(3) * t185 - t39;
t61 = t90 * pkin(4) + t88 * pkin(9);
t98 = -t166 - t217;
t148 = -t98 - t166;
t147 = -t30 * t67 + t229;
t146 = -t67 * t84 + t221;
t143 = t118 * t19 - t122 * t142;
t55 = t119 * t86 + t232 * t87;
t45 = -t118 * t55 - t122 * t200;
t139 = t118 * t200 - t122 * t55;
t81 = t165 + t205;
t135 = qJD(3) * (-t148 - t217);
t3 = -qJD(5) * t19 - t118 * t10 + t122 * t27;
t134 = -t143 * qJD(5) - t3 * t118;
t111 = t119 * pkin(3) + pkin(9);
t133 = -t111 * t60 - t84 * t155 + t30 * t88;
t131 = t119 * (-pkin(8) * t187 - t81);
t129 = t134 + t1;
t57 = -t97 * t186 + t202;
t58 = -qJD(3) * t205 + (-qJD(3) * t117 - t168) * t191;
t80 = t104 - t210;
t128 = -t58 * t120 + t57 * t123 + (-t120 * t81 - t123 * t80) * qJD(3);
t127 = t83 * t88 - t130;
t112 = -t232 * pkin(3) - pkin(4);
t70 = -t87 * qJD(3) - t154;
t69 = t86 * qJD(3) + t153;
t53 = t61 + t178;
t48 = -t88 ^ 2 + t90 ^ 2;
t44 = t90 * t181 - t93 * t236;
t43 = (qJD(2) * t137 + t88) * t181;
t40 = t232 * t71 - t211;
t25 = t55 * qJD(4) + t119 * t69 - t232 * t70;
t24 = t140 * qJD(4) + t119 * t70 + t232 * t69;
t23 = t118 * t61 + t122 * t35;
t22 = -t118 * t35 + t122 * t61;
t21 = t118 * t53 + t122 * t40;
t20 = -t118 * t40 + t122 * t53;
t17 = t84 * t158 - t77 * t90 + t214;
t16 = -t84 ^ 2 * t118 + t75 * t90 + t208;
t15 = t84 * t213 - t203;
t14 = t77 * t158 - t204;
t8 = t139 * qJD(5) - t118 * t24 + t122 * t169;
t7 = t45 * qJD(5) + t118 * t169 + t122 * t24;
t4 = (-t33 - t227) * t122 + (-t34 - t225) * t118;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, -t124 * t199, 0, 0, 0, 0, 0, 0, 0, 0, -t123 * t173 + (t70 - t154) * qJD(3), t120 * t173 + (-t69 - t153) * qJD(3), (-t120 * t70 + t123 * t69 + (-t120 * t87 - t123 * t86) * qJD(3)) * qJD(2), t57 * t87 + t58 * t86 + t81 * t69 + t80 * t70 + (t98 - t166) * t169, 0, 0, 0, 0, 0, 0, -t25 * t181 + (-t124 * t60 + t88 * t188) * t116, -t24 * t181 + (-t124 * t59 + t188 * t90) * t116, -t140 * t59 - t24 * t88 + t25 * t90 - t55 * t60, t10 * t55 - t231 + t36 * t24 - t35 * t25 + (-t124 * t85 + t188 * t83) * t116, 0, 0, 0, 0, 0, 0, -t140 * t34 + t25 * t75 + t45 * t60 + t8 * t84, t139 * t60 + t140 * t33 + t25 * t77 - t7 * t84, t139 * t34 + t45 * t33 - t7 * t75 - t8 * t77, -t139 * t2 - t142 * t8 + t19 * t7 + t30 * t25 + t3 * t45 - t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t150, -0.2e1 * t193 * t182, t194, -0.2e1 * t150, -t195, 0, -pkin(7) * t194 + t120 * t135, pkin(7) * t195 + t123 * t135, (-t114 - t115) * t151 + t128, ((t120 * t80 - t123 * t81) * t124 + (-t98 - t217) * t121) * t192 + t128 * pkin(7), t59 * t93 - t90 * t67, t137 * t59 + t67 * t88 - t90 * t68 - t221, -t67 * t181, t88 * t68 - t228, -t68 * t181, 0, t113 * t60 + t136 * t88 - t137 * t85 - t219 * t181 + t83 * t68, t113 * t59 + t136 * t90 - t220 * t181 - t83 * t67 + t85 * t93, t10 * t137 - t138 * t59 + t219 * t90 - t220 * t88 + t35 * t67 - t36 * t68 - t79 * t60 + t229, t10 * t79 + t85 * t113 + t136 * t83 - t219 * t35 + t220 * t36 - t230, -t77 * t175 + (-t33 * t93 - t67 * t77) * t122, (t118 * t77 + t122 * t75) * t67 + (t204 - t203 + (-t207 + t213) * qJD(5)) * t93, t122 * t146 + t137 * t33 - t175 * t84 + t77 * t68, t75 * t174 + (t34 * t93 - t67 * t75) * t118, -t118 * t146 + t137 * t34 - t174 * t84 - t75 * t68, t84 * t68 - t228, t147 * t118 - t137 * t3 - t138 * t34 - t142 * t68 + t219 * t75 - t233 * t84 + t93 * t29 + t37 * t60, t147 * t122 + t137 * t2 + t138 * t33 - t19 * t68 + t219 * t77 - t234 * t84 - t93 * t28 - t38 * t60, t37 * t33 - t38 * t34 + t233 * t77 - t234 * t75 + t143 * t67 + (-t237 * qJD(5) - t118 * t2 - t122 * t3) * t93, t142 * t233 + t234 * t19 + t2 * t38 + t219 * t30 + t3 * t37 - t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, t193 * t126, 0, t172, 0, 0, t148 * t189, -t98 * t187 + (t80 + t210) * qJD(3) - t202, 0, 0, t222, t48, t43, -t222, t44, 0, -t72 * t163 + t49 + t39 * t181 - t88 * t178 - t224 + (-qJD(4) * t66 - t181 * t218 - t50) * t119, (-t177 + t40) * qJD(4) + (-t131 + t40) * qJD(3) + (-t163 * t181 - t189 * t90) * pkin(3) + t127, (t36 - t39) * t90 + (-t35 + t40) * t88 + (-t232 * t59 - t119 * t60 + (t119 * t90 - t232 * t88) * qJD(4)) * pkin(3), t35 * t39 - t36 * t40 + (-t83 * t189 - t232 * t11 + t10 * t119 + (-t119 * t35 + t232 * t36) * qJD(4)) * pkin(3), t14, t4, t17, t15, t16, -t223, t112 * t34 - t20 * t84 + t149 * t75 + (-qJD(5) * t111 * t84 - t11) * t122 + t133 * t118 + t164, -t112 * t33 + (t111 * t184 + t21) * t84 + t149 * t77 + t133 * t122 + t156, t20 * t77 + t21 * t75 + (-t75 * t155 - t111 * t34 + (t111 * t77 + t142) * qJD(5)) * t122 + (t77 * t155 - t111 * t33 - t3 + (t111 * t75 - t19) * qJD(5)) * t118 + t180, t11 * t112 + t142 * t20 - t19 * t21 - t30 * t39 + (t119 * t30 + t237 * t232) * t218 + t129 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t48, t43, -t222, t44, 0, t36 * qJD(3) - t159 - t224, (-t177 + t35) * qJD(4) + (-t131 + t35) * qJD(3) + t127, 0, 0, t14, t4, t17, t15, t16, -t223, t30 * t212 - pkin(4) * t34 - t11 * t122 - t22 * t84 - t36 * t75 + (-t183 * t84 - t214) * pkin(9) + t164, t30 * t206 + pkin(4) * t33 + t23 * t84 - t36 * t77 + (t184 * t84 - t208) * pkin(9) + t156, t22 * t77 + t23 * t75 + (-t204 - t203 + (t207 + t213) * qJD(5)) * pkin(9) + t134 + t180, -t11 * pkin(4) + pkin(9) * t129 + t142 * t22 - t19 * t23 - t30 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, -t75 ^ 2 + t77 ^ 2, -t33 + t227, -t226, t225 - t34, t60, t19 * t84 - t30 * t77 + t3, -t142 * t84 + t30 * t75 - t2, 0, 0;];
tauc_reg = t5;
