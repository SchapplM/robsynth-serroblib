% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:06
% EndTime: 2020-01-03 12:02:10
% DurationCPUTime: 1.10s
% Computational Cost: add. (2687->204), mult. (5054->276), div. (0->0), fcn. (3173->8), ass. (0->160)
t114 = cos(pkin(9));
t120 = cos(qJ(2));
t183 = pkin(1) * qJD(1);
t157 = t120 * t183;
t113 = sin(pkin(9));
t117 = sin(qJ(2));
t159 = t117 * t183;
t95 = t113 * t159;
t77 = t114 * t157 - t95;
t118 = cos(qJ(5));
t119 = cos(qJ(4));
t162 = qJD(4) * t119;
t167 = t118 * t119;
t200 = -qJD(5) * t167 - t118 * t162;
t110 = qJD(1) + qJD(2);
t115 = sin(qJ(5));
t116 = sin(qJ(4));
t85 = t115 * t119 + t118 * t116;
t72 = t85 * t110;
t163 = qJD(4) * t116;
t105 = t119 * qJD(3);
t68 = t77 * qJD(2);
t172 = qJD(4) * t105 + t119 * t68;
t88 = t110 * pkin(2) + t157;
t62 = t113 * t88 + t114 * t159;
t58 = t110 * pkin(7) + t62;
t26 = -t58 * t163 + t172;
t178 = t116 * t68;
t161 = t116 * qJD(3);
t46 = t119 * t58 + t161;
t27 = -t46 * qJD(4) - t178;
t179 = t116 * t58;
t45 = t105 - t179;
t199 = -t27 * t116 + t26 * t119 + (-t116 * t46 - t119 * t45) * qJD(4);
t148 = pkin(8) * t110 + t58;
t135 = t148 * t116;
t19 = -qJD(4) * t135 + t172;
t39 = t105 - t135;
t36 = qJD(4) * pkin(4) + t39;
t198 = (qJD(5) * t36 + t19) * t118;
t109 = qJD(4) + qJD(5);
t40 = t148 * t119 + t161;
t104 = t120 * pkin(1) + pkin(2);
t169 = t114 * t117;
t184 = pkin(1) * t169 + t113 * t104;
t74 = pkin(7) + t184;
t197 = -pkin(8) - t74;
t196 = t114 * pkin(2);
t195 = t119 * pkin(4);
t153 = -pkin(3) - t195;
t61 = t114 * t88 - t95;
t47 = t153 * t110 - t61;
t194 = t47 * t72;
t168 = t115 * t116;
t155 = t110 * t168;
t70 = -t110 * t167 + t155;
t193 = t72 * t70;
t100 = t113 * pkin(2) + pkin(7);
t192 = -pkin(8) - t100;
t82 = t192 * t116;
t107 = t119 * pkin(8);
t83 = t119 * t100 + t107;
t48 = -t115 * t83 + t118 * t82;
t144 = qJD(4) * t192;
t79 = t116 * t144;
t80 = t119 * t144;
t84 = -t167 + t168;
t191 = t48 * qJD(5) + t115 * t80 + t118 * t79 + t84 * t77;
t49 = t115 * t82 + t118 * t83;
t190 = -t49 * qJD(5) - t115 * t79 + t118 * t80 + t85 * t77;
t55 = t109 * t85;
t152 = t110 * t163;
t128 = t113 * t120 + t169;
t126 = pkin(1) * t128;
t76 = qJD(2) * t126;
t67 = qJD(1) * t76;
t56 = pkin(4) * t152 + t67;
t189 = t47 * t55 + t56 * t84;
t42 = t55 * t110;
t134 = t109 * t168;
t54 = t134 + t200;
t188 = -t85 * t42 + t54 * t70;
t187 = -t47 * t54 + t56 * t85;
t57 = -t110 * pkin(3) - t61;
t186 = t67 * t116 + t57 * t162;
t185 = t200 * t110;
t182 = pkin(1) * qJD(2);
t181 = t110 * t57;
t180 = t115 * t40;
t177 = t118 * t40;
t50 = t54 * t109;
t75 = qJD(1) * t126;
t176 = t75 * t110;
t175 = t76 * t110;
t174 = t77 * t110;
t170 = t113 * t117;
t78 = (t114 * t120 - t170) * t182;
t173 = t78 * t110;
t171 = t110 * t116;
t121 = qJD(4) ^ 2;
t166 = t121 * t116;
t111 = t116 ^ 2;
t112 = t119 ^ 2;
t165 = t111 - t112;
t164 = t111 + t112;
t160 = pkin(4) * t171;
t158 = pkin(4) * t163;
t12 = t118 * t36 - t180;
t13 = t115 * t36 + t177;
t20 = -t40 * qJD(4) - t178;
t145 = -qJD(5) * t180 + t115 * t20;
t3 = t145 + t198;
t146 = -t115 * t19 + t118 * t20;
t4 = -t13 * qJD(5) + t146;
t156 = t12 * t54 - t13 * t55 - t3 * t84 - t4 * t85;
t108 = t110 ^ 2;
t154 = t116 * t108 * t119;
t149 = -pkin(4) * t109 - t36;
t147 = qJD(4) * t197;
t142 = -pkin(1) * t170 + t114 * t104;
t140 = t119 * t152;
t73 = -pkin(3) - t142;
t139 = -t75 + t158;
t138 = (-qJD(2) + t110) * t183;
t137 = (-qJD(1) - t110) * t182;
t41 = t110 * t134 + t185;
t136 = -t84 * t41 + t72 * t55;
t133 = t121 * t74 + t175;
t59 = t197 * t116;
t60 = t119 * t74 + t107;
t29 = -t115 * t60 + t118 * t59;
t30 = t115 * t59 + t118 * t60;
t132 = t116 * t45 - t119 * t46;
t131 = qJD(4) * (t110 * t73 - t78);
t130 = t100 * t121 - t176;
t101 = -pkin(3) - t196;
t129 = qJD(4) * (t101 * t110 + t77);
t127 = t47 * t70 - t145;
t125 = t128 * qJD(1) * t182;
t106 = t121 * t119;
t91 = t153 - t196;
t87 = -0.2e1 * t140;
t86 = 0.2e1 * t140;
t69 = -0.2e1 * t165 * t110 * qJD(4);
t66 = t73 - t195;
t65 = t76 + t158;
t52 = t57 * t163;
t51 = t55 * t109;
t34 = -t116 * t78 + t119 * t147;
t33 = t116 * t147 + t119 * t78;
t28 = -t70 ^ 2 + t72 ^ 2;
t22 = -t185 + (-t155 + t70) * t109;
t15 = t118 * t39 - t180;
t14 = -t115 * t39 - t177;
t11 = t42 * t84 + t70 * t55;
t10 = -t41 * t85 - t72 * t54;
t7 = -t30 * qJD(5) - t115 * t33 + t118 * t34;
t6 = t29 * qJD(5) + t115 * t34 + t118 * t33;
t5 = -t136 + t188;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t137, t120 * t137, 0, 0, 0, 0, 0, 0, 0, 0, -t125 - t175, -t68 - t173, 0, -t67 * t142 + t68 * t184 - t61 * t76 + t62 * t78, t86, t69, t106, t87, -t166, 0, t52 + t116 * t131 + (-t133 - t67) * t119, t133 * t116 + t119 * t131 + t186, t164 * t173 + t199, -t132 * t78 + t199 * t74 + t57 * t76 + t67 * t73, t10, t5, -t50, t11, -t51, 0, t7 * t109 + t66 * t42 + t65 * t70 + t189, -t6 * t109 - t66 * t41 + t65 * t72 + t187, t29 * t41 - t30 * t42 - t6 * t70 - t7 * t72 + t156, t12 * t7 + t13 * t6 + t4 * t29 + t3 * t30 + t47 * t65 + t56 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t138, t120 * t138, 0, 0, 0, 0, 0, 0, 0, 0, -t125 + t176, -t68 + t174, 0, t61 * t75 - t62 * t77 + (t113 * t68 - t114 * t67) * pkin(2), t86, t69, t106, t87, -t166, 0, t52 + t116 * t129 + (-t130 - t67) * t119, t130 * t116 + t119 * t129 + t186, -t164 * t174 + t199, t100 * t199 + t67 * t101 + t132 * t77 - t57 * t75, t10, t5, -t50, t11, -t51, 0, t190 * t109 + t139 * t70 + t91 * t42 + t189, -t191 * t109 + t139 * t72 - t91 * t41 + t187, -t190 * t72 - t191 * t70 + t48 * t41 - t49 * t42 + t156, t190 * t12 + t191 * t13 + t139 * t47 + t3 * t49 + t4 * t48 + t56 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t106, 0, -t132 * qJD(4) + t26 * t116 + t27 * t119, 0, 0, 0, 0, 0, 0, -t51, t50, t136 + t188, -t12 * t55 - t13 * t54 + t3 * t85 - t4 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t165 * t108, 0, t154, 0, 0, (-t68 - t181) * t116, -t119 * t181 + (t45 + t179) * qJD(4) - t172, 0, 0, t193, t28, t22, -t193, 0, 0, -t70 * t160 - t14 * t109 - t194 + (t149 * t115 - t177) * qJD(5) + t146, -t72 * t160 + t15 * t109 + (t149 * qJD(5) - t19) * t118 + t127, (t13 + t14) * t72 + (-t12 + t15) * t70 + (-t115 * t42 + t118 * t41 + (t115 * t72 - t118 * t70) * qJD(5)) * pkin(4), -t12 * t14 - t13 * t15 + (-t47 * t171 + t115 * t3 + t118 * t4 + (-t115 * t12 + t118 * t13) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t28, t22, -t193, 0, 0, t13 * t109 - t194 + t4, t12 * t109 + t127 - t198, 0, 0;];
tauc_reg = t1;
