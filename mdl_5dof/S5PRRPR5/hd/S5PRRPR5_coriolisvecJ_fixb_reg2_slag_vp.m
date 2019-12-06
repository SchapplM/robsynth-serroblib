% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:54
% EndTime: 2019-12-05 16:28:07
% DurationCPUTime: 2.47s
% Computational Cost: add. (2978->307), mult. (7921->446), div. (0->0), fcn. (6008->10), ass. (0->171)
t109 = sin(qJ(2));
t105 = sin(pkin(5));
t168 = qJD(1) * t105;
t147 = t109 * t168;
t108 = sin(qJ(3));
t163 = qJD(3) * t108;
t211 = pkin(3) * t163 - t147;
t104 = sin(pkin(10));
t111 = cos(qJ(3));
t189 = -qJ(4) - pkin(7);
t141 = qJD(3) * t189;
t119 = -t108 * qJD(4) + t111 * t141;
t177 = cos(pkin(10));
t139 = t177 * t111;
t176 = t104 * t108;
t120 = t139 - t176;
t112 = cos(qJ(2));
t146 = t112 * t168;
t158 = t111 * qJD(4);
t72 = t108 * t141 + t158;
t186 = t104 * t119 - t120 * t146 + t177 * t72;
t140 = t177 * t108;
t84 = t104 * t111 + t140;
t77 = t84 * qJD(3);
t80 = t120 * qJD(3);
t210 = -t77 * pkin(4) + t80 * pkin(8) - t211;
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t166 = qJD(2) * t105;
t143 = qJD(1) * t166;
t133 = t112 * t143;
t106 = cos(pkin(5));
t172 = t106 * t111;
t94 = qJD(1) * t172;
t185 = qJD(3) * t94 + t111 * t133;
t87 = qJD(2) * pkin(7) + t147;
t41 = -t163 * t87 + t185;
t30 = (-qJ(4) * t163 + t158) * qJD(2) + t41;
t144 = t177 * t30;
t126 = -qJD(4) - t146;
t162 = qJD(3) * t111;
t148 = qJ(4) * t162;
t167 = qJD(1) * t108;
t145 = t106 * t167;
t63 = t111 * t87 + t145;
t206 = -t63 * qJD(3) + (t108 * t126 - t148) * qJD(2);
t12 = t104 * t206 + t144;
t136 = qJ(4) * qJD(2) + t87;
t54 = t111 * t136 + t145;
t47 = t177 * t54;
t53 = -t108 * t136 + t94;
t49 = qJD(3) * pkin(3) + t53;
t18 = t104 * t49 + t47;
t14 = qJD(3) * pkin(8) + t18;
t100 = -t111 * pkin(3) - pkin(2);
t70 = qJD(2) * t100 + qJD(4) - t146;
t165 = qJD(2) * t108;
t93 = qJD(2) * t139;
t75 = t104 * t165 - t93;
t78 = t84 * qJD(2);
t26 = t75 * pkin(4) - t78 * pkin(8) + t70;
t124 = t107 * t14 - t110 * t26;
t68 = qJD(2) * t77;
t157 = qJD(2) * qJD(3);
t142 = t108 * t157;
t90 = t104 * t142;
t69 = qJD(3) * t93 - t90;
t73 = pkin(3) * t142 + t109 * t143;
t24 = t68 * pkin(4) - t69 * pkin(8) + t73;
t1 = -qJD(5) * t124 + t107 * t24 + t110 * t12;
t71 = qJD(5) + t75;
t209 = t124 * t71 + t1;
t8 = t107 * t26 + t110 * t14;
t2 = -qJD(5) * t8 - t107 * t12 + t110 * t24;
t208 = t8 * t71 + t2;
t138 = t107 * t71;
t61 = t107 * qJD(3) + t110 * t78;
t207 = t61 * t138;
t205 = t78 ^ 2;
t45 = -pkin(4) * t120 - t84 * pkin(8) + t100;
t89 = t189 * t111;
t58 = t189 * t176 - t177 * t89;
t15 = -t107 * t58 + t110 * t45;
t202 = qJD(5) * t15 - t210 * t107 + t186 * t110;
t16 = t107 * t45 + t110 * t58;
t201 = qJD(5) * t16 + t186 * t107 + t210 * t110;
t200 = pkin(3) * t108;
t11 = t104 * t30 - t177 * t206;
t175 = t105 * t109;
t81 = -t108 * t175 + t172;
t82 = t106 * t108 + t111 * t175;
t39 = t104 * t82 - t177 * t81;
t199 = t11 * t39;
t57 = -t104 * t89 - t189 * t140;
t198 = t11 * t57;
t197 = t11 * t84;
t159 = t110 * qJD(3);
t59 = t107 * t78 - t159;
t196 = t59 * t75;
t195 = t61 * t59;
t194 = t61 * t78;
t193 = t68 * t120;
t192 = t78 * t59;
t191 = t78 * t75;
t190 = t84 * t68;
t160 = qJD(5) * t110;
t181 = t107 * t69;
t32 = qJD(5) * t61 + t181;
t188 = -t107 * t32 - t59 * t160;
t187 = t104 * t72 - t119 * t177 - t84 * t146;
t184 = qJD(2) * pkin(2);
t183 = t104 * t54;
t182 = t107 * t68;
t180 = t108 * t87;
t161 = qJD(5) * t107;
t31 = -qJD(5) * t159 - t110 * t69 + t161 * t78;
t179 = t31 * t107;
t178 = t32 * t110;
t174 = t105 * t112;
t114 = qJD(2) ^ 2;
t173 = t105 * t114;
t113 = qJD(3) ^ 2;
t171 = t113 * t108;
t170 = t113 * t111;
t102 = t108 ^ 2;
t103 = t111 ^ 2;
t169 = t102 - t103;
t164 = qJD(2) * t109;
t155 = pkin(3) * t165;
t154 = t84 * t161;
t153 = t84 * t160;
t152 = t109 * t173;
t151 = t108 * t114 * t111;
t150 = t105 * t164;
t149 = t112 * t166;
t137 = t110 * t71;
t135 = t108 * t149;
t134 = t111 * t149;
t132 = t111 * t142;
t97 = t104 * pkin(3) + pkin(8);
t131 = qJD(5) * t71 * t97 + t11;
t88 = -t146 - t184;
t130 = -t88 - t146;
t129 = -t107 * t8 + t110 * t124;
t17 = t177 * t49 - t183;
t13 = -qJD(3) * pkin(4) - t17;
t128 = t13 * t80 + t197;
t127 = t71 * t80 + t190;
t123 = t110 * t68 - t75 * t138 - t161 * t71;
t40 = t104 * t81 + t177 * t82;
t28 = -t107 * t40 - t110 * t174;
t122 = t107 * t174 - t110 * t40;
t121 = t13 * t71 - t68 * t97;
t117 = qJD(3) * (-t130 - t184);
t42 = -t87 * t162 + (-qJD(3) * t106 - t149) * t167;
t62 = t94 - t180;
t116 = -t42 * t108 + t41 * t111 + (-t108 * t63 - t111 * t62) * qJD(3);
t98 = -pkin(3) * t177 - pkin(4);
t74 = t75 ^ 2;
t52 = -qJD(3) * t82 - t135;
t51 = qJD(3) * t81 + t134;
t35 = t78 * pkin(4) + t75 * pkin(8) + t155;
t22 = t177 * t53 - t183;
t21 = t104 * t52 + t177 * t51;
t20 = t104 * t53 + t47;
t19 = t104 * t51 - t177 * t52;
t10 = t107 * t35 + t110 * t22;
t9 = -t107 * t22 + t110 * t35;
t6 = qJD(5) * t122 - t107 * t21 + t110 * t150;
t5 = qJD(5) * t28 + t107 * t150 + t110 * t21;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t112 * t173, 0, 0, 0, 0, 0, 0, 0, 0, -t111 * t152 + (t52 - t135) * qJD(3), t108 * t152 + (-t51 - t134) * qJD(3), (-t108 * t52 + t111 * t51 + (-t108 * t82 - t111 * t81) * qJD(3)) * qJD(2), t41 * t82 + t42 * t81 + t63 * t51 + t62 * t52 + (t88 - t146) * t150, 0, 0, 0, 0, 0, 0, -t19 * qJD(3) + (-t112 * t68 + t164 * t75) * t105, -t21 * qJD(3) + (-t112 * t69 + t164 * t78) * t105, t19 * t78 - t21 * t75 + t39 * t69 - t40 * t68, t199 + t12 * t40 - t17 * t19 + t18 * t21 + (-t112 * t73 + t164 * t70) * t105, 0, 0, 0, 0, 0, 0, t19 * t59 + t28 * t68 + t39 * t32 + t6 * t71, t122 * t68 + t19 * t61 - t39 * t31 - t5 * t71, t122 * t32 + t28 * t31 - t5 * t59 - t6 * t61, -t1 * t122 - t124 * t6 + t13 * t19 + t2 * t28 + t5 * t8 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132, -0.2e1 * t169 * t157, t170, -0.2e1 * t132, -t171, 0, -pkin(7) * t170 + t108 * t117, pkin(7) * t171 + t111 * t117, (-t102 - t103) * t133 + t116, ((t108 * t62 - t111 * t63) * t112 + (-t88 - t184) * t109) * t168 + t116 * pkin(7), t69 * t84 + t78 * t80, t120 * t69 - t80 * t75 - t78 * t77 - t190, t80 * qJD(3), t75 * t77 - t193, -t77 * qJD(3), 0, -t75 * t147 + t100 * t68 + t70 * t77 - t73 * t120 + (t75 * t200 - t187) * qJD(3), -t78 * t147 + t100 * t69 + t70 * t80 + t73 * t84 + (t78 * t200 - t186) * qJD(3), t12 * t120 - t17 * t80 - t18 * t77 - t186 * t75 + t187 * t78 + t57 * t69 - t58 * t68 + t197, t73 * t100 + t12 * t58 - t187 * t17 + t186 * t18 + t211 * t70 + t198, -t61 * t154 + (-t31 * t84 + t61 * t80) * t110, (-t107 * t61 - t110 * t59) * t80 + (t179 - t178 + (t107 * t59 - t110 * t61) * qJD(5)) * t84, t110 * t127 + t120 * t31 - t154 * t71 + t61 * t77, t59 * t153 + (t32 * t84 + t59 * t80) * t107, -t107 * t127 + t120 * t32 - t153 * t71 - t59 * t77, t71 * t77 - t193, t128 * t107 - t120 * t2 - t124 * t77 + t13 * t153 + t15 * t68 + t187 * t59 - t201 * t71 + t57 * t32, t1 * t120 + t128 * t110 - t13 * t154 - t16 * t68 + t187 * t61 - t202 * t71 - t57 * t31 - t8 * t77, t15 * t31 - t16 * t32 + t129 * t80 + t201 * t61 - t202 * t59 + (-t1 * t107 - t2 * t110 + (-t107 * t124 - t110 * t8) * qJD(5)) * t84, t1 * t16 + t124 * t201 + t187 * t13 + t15 * t2 + t202 * t8 + t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t169 * t114, 0, t151, 0, 0, t130 * t165, -t88 * t111 * qJD(2) + (t62 + t180) * qJD(3) - t185, 0, 0, t191, -t74 + t205, -t90 + (t93 + t75) * qJD(3), -t191, 0, 0, t20 * qJD(3) - t155 * t75 - t70 * t78 - t11, -t144 + t70 * t75 + (t104 * t63 + t22) * qJD(3) + (t104 * t148 + (-pkin(3) * t78 - t104 * t126) * t108) * qJD(2), (t18 - t20) * t78 + (-t17 + t22) * t75 + (-t104 * t68 - t177 * t69) * pkin(3), t17 * t20 - t18 * t22 + (t104 * t12 - t11 * t177 - t165 * t70) * pkin(3), t137 * t61 - t179, (-t31 - t196) * t110 - t207 + t188, t137 * t71 + t182 - t194, t138 * t59 - t178, t123 + t192, -t71 * t78, t107 * t121 - t110 * t131 + t124 * t78 - t20 * t59 + t98 * t32 - t9 * t71, t10 * t71 + t107 * t131 + t110 * t121 - t20 * t61 - t98 * t31 + t8 * t78, t10 * t59 + t9 * t61 + (-t32 * t97 + t124 * t75 + t1 + (t61 * t97 + t124) * qJD(5)) * t110 + (-t31 * t97 - t75 * t8 - t2 + (t59 * t97 - t8) * qJD(5)) * t107, -t8 * t10 + t11 * t98 - t13 * t20 + t124 * t9 + (qJD(5) * t129 + t1 * t110 - t2 * t107) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78 * qJD(3), -t90 + (t93 - t75) * qJD(3), -t74 - t205, t17 * t78 + t18 * t75 + t73, 0, 0, 0, 0, 0, 0, t123 - t192, -t110 * t71 ^ 2 - t182 - t194, (t31 - t196) * t110 + t207 + t188, t107 * t209 + t208 * t110 - t13 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, -t59 ^ 2 + t61 ^ 2, t59 * t71 - t31, -t195, -t181 + (-qJD(5) + t71) * t61, t68, -t13 * t61 + t208, t13 * t59 - t209, 0, 0;];
tauc_reg = t3;
