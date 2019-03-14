% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:51
% EndTime: 2019-03-09 04:29:59
% DurationCPUTime: 2.87s
% Computational Cost: add. (2661->263), mult. (6417->449), div. (0->0), fcn. (5440->8), ass. (0->145)
t111 = cos(qJ(4));
t109 = sin(qJ(4));
t177 = sin(pkin(10));
t142 = t177 * t109;
t178 = cos(pkin(10));
t76 = -t178 * t111 + t142;
t110 = sin(qJ(3));
t112 = cos(qJ(3));
t165 = t112 * qJD(3);
t168 = qJD(4) * t111;
t202 = t109 * t165 + t110 * t168;
t140 = qJD(3) * t177;
t133 = t112 * t140;
t141 = qJD(3) * t178;
t135 = t112 * t141;
t144 = t178 * t109;
t77 = t177 * t111 + t144;
t67 = t77 * qJD(4);
t38 = t109 * t133 + t110 * t67 - t111 * t135;
t62 = t76 * t110;
t184 = t38 * t76 + t62 * t67;
t68 = t76 * qJD(4);
t37 = -t109 * t135 + t68 * t110 - t111 * t133;
t61 = t77 * t110;
t198 = -t77 * t37 - t68 * t61;
t201 = t198 - t184;
t200 = t198 + t184;
t199 = -0.4e1 * t110;
t104 = t110 * qJD(3);
t197 = t61 * t104 + t37 * t112;
t185 = -t62 * t37 + t38 * t61;
t105 = t109 ^ 2;
t107 = t111 ^ 2;
t139 = qJD(4) * (t105 - t107);
t106 = t110 ^ 2;
t170 = -t112 ^ 2 + t106;
t138 = t170 * qJD(3);
t196 = 2 * qJD(6);
t174 = t110 * t111;
t98 = sin(pkin(9)) * pkin(1) + pkin(7);
t176 = qJD(4) * t98;
t194 = pkin(8) * t112;
t195 = pkin(3) * t110;
t129 = -t194 + t195;
t121 = t129 * qJD(3);
t100 = -cos(pkin(9)) * pkin(1) - pkin(2);
t191 = t112 * pkin(3);
t118 = t100 - t191;
t192 = t110 * pkin(8);
t114 = t118 - t192;
t64 = t111 * t114;
t182 = -qJD(4) * t64 - t109 * t121;
t15 = (-qJ(5) * qJD(4) - qJD(3) * t98) * t174 + (-t110 * qJD(5) + (-qJ(5) * qJD(3) - t176) * t112) * t109 - t182;
t166 = t111 * qJD(5);
t173 = t111 * t112;
t151 = t109 * t104;
t181 = t111 * t121 + t98 * t151;
t187 = -qJ(5) - pkin(8);
t193 = t110 * pkin(4);
t79 = t98 * t173;
t9 = -t110 * t166 + (-qJ(5) * t173 + t193) * qJD(3) + (-t79 + (-t187 * t110 - t118) * t109) * qJD(4) + t181;
t4 = t178 * t15 + t177 * t9;
t189 = t61 * t37;
t180 = t109 * t98;
t39 = -qJ(5) * t174 + t64 + (-pkin(4) - t180) * t112;
t175 = t109 * t110;
t51 = t109 * t114 + t79;
t45 = -qJ(5) * t175 + t51;
t14 = t177 * t39 + t178 * t45;
t66 = pkin(4) * t175 + t110 * t98;
t172 = t62 * qJD(6);
t169 = qJD(4) * t109;
t167 = qJD(4) * t112;
t164 = t112 * qJD(6);
t163 = -0.2e1 * t189;
t162 = 0.2e1 * t76 * t67;
t161 = -0.2e1 * pkin(3) * qJD(4);
t160 = qJ(6) * t104 + t4;
t80 = t98 * t165;
t52 = t202 * pkin(4) + t80;
t159 = 0.2e1 * qJD(3) * t100;
t103 = pkin(4) * t169;
t158 = t112 * t180;
t156 = t106 * t176;
t102 = -t111 * pkin(4) - pkin(3);
t155 = t109 * t167;
t153 = t111 * t167;
t152 = t105 * t165;
t149 = t109 * t168;
t148 = t110 * t165;
t147 = t111 * t165;
t145 = qJD(4) * t187;
t117 = -t109 * qJD(5) + t111 * t145;
t65 = t109 * t145 + t166;
t40 = -t178 * t117 + t177 * t65;
t41 = t177 * t117 + t178 * t65;
t87 = t187 * t111;
t53 = -t187 * t144 - t177 * t87;
t54 = t187 * t142 - t178 * t87;
t146 = t53 * t40 + t54 * t41;
t137 = t109 * t147;
t136 = t106 * t149;
t131 = t177 * t15 - t178 * t9;
t130 = -t191 - t192;
t128 = -t37 * t76 + t61 * t67;
t126 = t77 * t67 - t68 * t76;
t50 = t64 - t158;
t125 = -t109 * t51 - t111 * t50;
t124 = t109 * t50 - t111 * t51;
t123 = -t53 * t104 + t40 * t112;
t122 = t54 * t104 - t41 * t112;
t43 = t76 * t104 - t112 * t67;
t44 = t77 * t104 + t68 * t112;
t120 = -t37 * t53 - t38 * t54 + t61 * t40 - t62 * t41;
t119 = t54 * t37 - t53 * t38 - t40 * t62 - t41 * t61;
t70 = t111 * t104 + t155;
t13 = -t177 * t45 + t178 * t39;
t116 = 0.2e1 * t40 * t77 - 0.2e1 * t41 * t76 - 0.2e1 * t53 * t68 - 0.2e1 * t54 * t67;
t27 = t62 * t38;
t115 = -0.2e1 * t148 + 0.2e1 * t27 - 0.2e1 * t189;
t20 = t70 * t98 + t182;
t21 = -t51 * qJD(4) + t181;
t113 = t125 * qJD(4) - t21 * t109 - t20 * t111;
t99 = -t178 * pkin(4) - pkin(5);
t96 = t177 * pkin(4) + qJ(6);
t94 = t107 * t165;
t92 = -0.2e1 * t148;
t86 = t107 * t148;
t85 = t105 * t148;
t72 = t151 - t153;
t69 = t110 * t169 - t147;
t55 = t110 * t139 - t137;
t49 = -0.2e1 * t77 * t68;
t46 = t76 * pkin(5) - t77 * qJ(6) + t102;
t31 = t67 * pkin(5) + t68 * qJ(6) - t77 * qJD(6) + t103;
t30 = t61 * pkin(5) + t62 * qJ(6) + t66;
t22 = 0.2e1 * t27;
t16 = -0.2e1 * t62 * t104 + 0.2e1 * t112 * t38;
t11 = -t38 * t77 + t62 * t68;
t10 = t112 * pkin(5) - t13;
t8 = -t112 * qJ(6) + t14;
t5 = -t37 * pkin(5) + t38 * qJ(6) + t172 + t52;
t2 = -pkin(5) * t104 + t131;
t1 = t160 - t164;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t148, -0.2e1 * t138, 0, t92, 0, 0, t110 * t159, t112 * t159, 0, 0, 0.2e1 * t86 - 0.2e1 * t136, 0.2e1 * t106 * t139 + t137 * t199, 0.2e1 * t110 * t155 + 0.2e1 * t111 * t138, 0.2e1 * t85 + 0.2e1 * t136, -0.2e1 * t109 * t138 + 0.2e1 * t110 * t153, t92, 0.2e1 * t111 * t156 - 0.2e1 * t21 * t112 + 0.2e1 * (t50 + 0.2e1 * t158) * t104, -0.2e1 * t109 * t156 - 0.2e1 * t20 * t112 + 0.2e1 * (-t51 + 0.2e1 * t79) * t104, 0.2e1 * t125 * t165 + 0.2e1 * (qJD(4) * t124 + t109 * t20 - t111 * t21) * t110, 0.2e1 * t148 * t98 ^ 2 - 0.2e1 * t51 * t20 + 0.2e1 * t50 * t21, t22, 0.2e1 * t185, t16, t163, -0.2e1 * t197, t92, 0.2e1 * t13 * t104 + 0.2e1 * t112 * t131 - 0.2e1 * t66 * t37 + 0.2e1 * t52 * t61, -0.2e1 * t14 * t104 + 0.2e1 * t4 * t112 - 0.2e1 * t66 * t38 - 0.2e1 * t52 * t62, 0.2e1 * t13 * t38 - 0.2e1 * t131 * t62 + 0.2e1 * t14 * t37 - 0.2e1 * t4 * t61, -0.2e1 * t13 * t131 + 0.2e1 * t14 * t4 + 0.2e1 * t66 * t52, t22, t16, -0.2e1 * t185, t92, 0.2e1 * t197, t163, -0.2e1 * t10 * t104 + 0.2e1 * t2 * t112 - 0.2e1 * t30 * t37 + 0.2e1 * t5 * t61, -0.2e1 * t1 * t61 - 0.2e1 * t10 * t38 - 0.2e1 * t2 * t62 + 0.2e1 * t8 * t37, -0.2e1 * t1 * t112 + 0.2e1 * t8 * t104 + 0.2e1 * t30 * t38 + 0.2e1 * t5 * t62, 0.2e1 * t8 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t30 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * t110 + (-t124 * t112 + t170 * t98) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t104 - t52 * t112 + t13 * t37 + t131 * t61 - t14 * t38 - t4 * t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t62 - t10 * t37 + t30 * t104 - t5 * t112 + t2 * t61 - t8 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t85 + 0.2e1 * t86 - 0.2e1 * t148, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, 0, -t104, 0, -t80, t98 * t104, 0, 0, -t55, t149 * t199 - t152 + t94, t72, t55, t70, 0 (pkin(8) * t173 + (-pkin(3) * t111 + t180) * t110) * qJD(4) + (t109 * t130 - t79) * qJD(3) (t109 * t129 + t98 * t174) * qJD(4) + (t111 * t130 + t158) * qJD(3), t113, -pkin(3) * t80 + pkin(8) * t113, t11, -t201, t44, t128, -t43, 0, -t102 * t37 + t61 * t103 + t52 * t76 + t66 * t67 + t123, -t102 * t38 - t62 * t103 + t52 * t77 - t66 * t68 - t122, t13 * t68 + t131 * t77 - t14 * t67 - t4 * t76 + t119, t52 * t102 + t66 * t103 - t13 * t40 + t131 * t53 + t14 * t41 + t4 * t54, t11, t44, t201, 0, t43, t128, t30 * t67 + t31 * t61 - t46 * t37 + t5 * t76 + t123, -t1 * t76 - t10 * t68 + t2 * t77 - t8 * t67 + t119, t30 * t68 + t31 * t62 + t46 * t38 - t5 * t77 + t122, t1 * t54 + t10 * t40 + t2 * t53 + t30 * t31 + t8 * t41 + t5 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t165, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t72, t94 + t152 (-t195 + (t105 + t107) * t194) * qJD(3), 0, 0, 0, 0, 0, 0, t43, t44, t200, -pkin(4) * t155 + t102 * t104 + t120, 0, 0, 0, 0, 0, 0, t43, t200, -t44, t46 * t104 - t112 * t31 + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t149, -0.2e1 * t139, 0, -0.2e1 * t149, 0, 0, t109 * t161, t111 * t161, 0, 0, t49, -0.2e1 * t126, 0, t162, 0, 0, 0.2e1 * t102 * t67 + 0.2e1 * t76 * t103, -0.2e1 * t102 * t68 + 0.2e1 * t77 * t103, t116, 0.2e1 * t102 * t103 + 0.2e1 * t146, t49, 0, 0.2e1 * t126, 0, 0, t162, 0.2e1 * t31 * t76 + 0.2e1 * t46 * t67, t116, -0.2e1 * t31 * t77 + 0.2e1 * t46 * t68, 0.2e1 * t46 * t31 + 0.2e1 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, -t202, t104, t21, t20, 0, 0, 0, 0, -t38, 0, t37, t104, t141 * t193 - t131, -t140 * t193 - t4 (t177 * t37 + t178 * t38) * pkin(4) (-t131 * t178 + t177 * t4) * pkin(4), 0, -t38, 0, t104, -t37, 0 (pkin(5) - t99) * t104 - t131, -qJD(6) * t61 + t96 * t37 - t99 * t38, t96 * t104 + t160 - 0.2e1 * t164, t8 * qJD(6) + t1 * t96 + t2 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t69, 0, 0, 0, 0, 0, 0, 0, 0, t37, t38, 0 (-t177 * t38 + t178 * t37) * pkin(4), 0, 0, 0, 0, 0, 0, t37, 0, -t38, -t37 * t99 - t38 * t96 - t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, 0, -t169, 0, -pkin(8) * t168, pkin(8) * t169, 0, 0, 0, 0, -t68, 0, -t67, 0, -t40, -t41 (-t177 * t67 + t178 * t68) * pkin(4) (t177 * t41 - t178 * t40) * pkin(4), 0, -t68, 0, 0, t67, 0, -t40, -qJD(6) * t76 - t96 * t67 - t99 * t68, t41, t54 * qJD(6) + t40 * t99 + t41 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t96 * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, t52, 0, 0, 0, 0, 0, 0, -t37, 0, t38, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t68, 0, t103, 0, 0, 0, 0, 0, 0, t67, 0, t68, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t38, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;