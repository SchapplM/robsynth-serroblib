% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:04
% EndTime: 2019-12-31 18:52:12
% DurationCPUTime: 2.53s
% Computational Cost: add. (3660->289), mult. (9720->382), div. (0->0), fcn. (7118->6), ass. (0->150)
t119 = cos(pkin(8));
t112 = -t119 * pkin(2) - pkin(1);
t104 = t112 * qJD(1) + qJD(2);
t205 = qJD(2) + t104;
t118 = sin(pkin(8));
t121 = sin(qJ(3));
t191 = cos(qJ(3));
t130 = -t121 * t118 + t191 * t119;
t204 = t130 * qJD(1);
t125 = qJD(3) * t204;
t203 = qJD(4) * qJD(3) + t125;
t120 = sin(qJ(4));
t122 = cos(qJ(4));
t155 = qJD(4) * t122;
t156 = qJD(4) * t120;
t126 = t130 * qJD(2);
t178 = pkin(6) + qJ(2);
t105 = t178 * t118;
t102 = qJD(1) * t105;
t106 = t178 * t119;
t103 = qJD(1) * t106;
t198 = -t191 * t102 - t121 * t103;
t43 = qJD(1) * t126 + qJD(3) * t198;
t101 = t191 * t118 + t121 * t119;
t95 = t101 * qJD(1);
t52 = -pkin(3) * t204 - t95 * pkin(7) + t104;
t98 = t101 * qJD(3);
t87 = qJD(1) * t98;
t58 = t87 * pkin(3) - pkin(7) * t125;
t75 = -t121 * t102 + t191 * t103;
t70 = qJD(3) * pkin(7) + t75;
t145 = -t120 * t58 - t122 * t43 - t52 * t155 + t70 * t156;
t31 = -t120 * t70 + t122 * t52;
t88 = qJD(4) - t204;
t139 = -t31 * t88 - t145;
t32 = t120 * t52 + t122 * t70;
t81 = -t122 * qJD(3) + t120 * t95;
t24 = -t81 * qJ(5) + t32;
t190 = t24 * t88;
t11 = -qJD(4) * t32 - t120 * t43 + t122 * t58;
t134 = t203 * t122 - t95 * t156;
t124 = -qJ(5) * t134 + t11;
t192 = t87 * pkin(4);
t83 = t120 * qJD(3) + t122 * t95;
t2 = -t83 * qJD(5) + t124 + t192;
t202 = t2 + t190;
t144 = t120 * t88;
t201 = t83 * t144;
t200 = t32 * t88 + t11;
t199 = t122 * t87 - t88 * t156;
t197 = -t191 * t105 - t121 * t106;
t195 = t83 ^ 2;
t194 = t95 ^ 2;
t193 = t81 * pkin(4);
t129 = t101 * qJD(2);
t44 = qJD(1) * t129 + t75 * qJD(3);
t187 = t44 * t197;
t186 = t81 * t88;
t185 = t81 * t204;
t184 = t83 * t81;
t183 = t83 * t88;
t182 = t83 * t95;
t181 = t88 * t95;
t180 = t95 * t81;
t179 = t95 * t204;
t177 = -qJ(5) - pkin(7);
t23 = -t83 * qJ(5) + t31;
t16 = t88 * pkin(4) + t23;
t176 = t16 - t23;
t146 = qJD(4) * t177;
t158 = t122 * qJ(5);
t71 = t95 * pkin(3) - pkin(7) * t204;
t35 = -t120 * t198 + t122 * t71;
t175 = t95 * pkin(4) + t120 * qJD(5) - t122 * t146 - t158 * t204 + t35;
t161 = t120 * qJ(5);
t36 = t120 * t71 + t122 * t198;
t174 = -t122 * qJD(5) - t120 * t146 - t161 * t204 + t36;
t49 = t203 * t120 + t95 * t155;
t173 = -t120 * t49 - t81 * t155;
t73 = -pkin(3) * t130 - t101 * pkin(7) + t112;
t79 = -t121 * t105 + t191 * t106;
t76 = t122 * t79;
t39 = t120 * t73 + t76;
t172 = t120 * t83;
t171 = t120 * t87;
t170 = t120 * t204;
t97 = t130 * qJD(3);
t169 = t120 * t97;
t148 = qJD(5) + t193;
t69 = -qJD(3) * pkin(3) - t198;
t41 = t148 + t69;
t168 = t122 * t41;
t167 = t122 * t97;
t26 = t49 * pkin(4) + t44;
t166 = t26 * t120;
t165 = t26 * t122;
t164 = t134 * t120;
t163 = t49 * t122;
t162 = t87 * t130;
t157 = t118 ^ 2 + t119 ^ 2;
t53 = t197 * qJD(3) + t126;
t72 = t98 * pkin(3) - t97 * pkin(7);
t153 = t120 * t72 + t122 * t53 + t73 * t155;
t149 = t101 * t155;
t147 = -t120 * t53 + t122 * t72;
t38 = -t120 * t79 + t122 * t73;
t143 = t122 * t88;
t142 = t157 * qJD(1) ^ 2;
t141 = pkin(7) * qJD(4) * t88 + t44;
t133 = t49 * qJ(5) + t145;
t3 = -t81 * qJD(5) - t133;
t140 = -t88 * t16 + t3;
t138 = t120 * t32 + t122 * t31;
t137 = -qJ(5) * t97 - qJD(5) * t101;
t136 = t88 * t170 + t199;
t135 = 0.2e1 * t157 * qJD(2) * qJD(1);
t132 = -pkin(7) * t87 + t88 * t69;
t54 = t79 * qJD(3) + t129;
t114 = -t122 * pkin(4) - pkin(3);
t108 = t177 * t122;
t107 = t177 * t120;
t91 = t204 ^ 2;
t80 = t81 ^ 2;
t55 = t120 * t101 * pkin(4) - t197;
t45 = pkin(4) * t170 + t75;
t40 = t88 * t98 - t162;
t37 = -t80 + t195;
t34 = (t149 + t169) * pkin(4) + t54;
t33 = -t101 * t161 + t39;
t30 = -t49 + t183;
t29 = t134 + t186;
t27 = -pkin(4) * t130 - t101 * t158 + t38;
t22 = -t88 ^ 2 * t122 - t171 - t182;
t21 = t88 * t143 + t171 - t182;
t20 = t136 + t180;
t19 = t136 - t180;
t18 = t81 * t144 - t163;
t17 = t83 * t143 + t164;
t15 = -t173 * t101 + t81 * t169;
t14 = t83 * t167 + (t122 * t134 - t83 * t156) * t101;
t13 = -t39 * qJD(4) + t147;
t12 = -t79 * t156 + t153;
t9 = -t88 * t169 + t49 * t130 - t81 * t98 + (-t88 * t155 - t171) * t101;
t8 = t199 * t101 - t130 * t134 + t88 * t167 + t83 * t98;
t7 = -qJ(5) * t149 + (-qJD(4) * t79 + t137) * t120 + t153;
t6 = t98 * pkin(4) + t137 * t122 + (-t76 + (qJ(5) * t101 - t73) * t120) * qJD(4) + t147;
t5 = (t134 + t185) * t122 - t201 + t173;
t4 = (-t134 + t185) * t122 + t201 + t173;
t1 = -(t122 * t81 + t172) * t97 + (-t164 - t163 + (t120 * t81 - t122 * t83) * qJD(4)) * t101;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, qJ(2) * t135, t101 * t125 + t95 * t97, -t101 * t87 + t125 * t130 + t204 * t97 - t95 * t98, t97 * qJD(3), -t204 * t98 - t162, -t98 * qJD(3), 0, -t54 * qJD(3) + t104 * t98 + t112 * t87, -t53 * qJD(3) + t104 * t97 + t112 * t125, t44 * t101 - t125 * t197 + t130 * t43 - t198 * t97 + t204 * t53 + t54 * t95 - t75 * t98 - t79 * t87, -t198 * t54 + t43 * t79 + t75 * t53 - t187, t14, t1, t8, t15, t9, t40, t69 * t169 - t11 * t130 + t13 * t88 + t31 * t98 + t38 * t87 - t197 * t49 + t54 * t81 + (t44 * t120 + t155 * t69) * t101, t69 * t167 - t145 * t130 - t12 * t88 - t32 * t98 - t39 * t87 - t197 * t134 + t54 * t83 + (t44 * t122 - t156 * t69) * t101, -t12 * t81 - t13 * t83 - t38 * t134 - t39 * t49 - t138 * t97 + (t145 * t120 - t11 * t122 + (t120 * t31 - t122 * t32) * qJD(4)) * t101, t11 * t38 + t32 * t12 + t31 * t13 - t145 * t39 + t69 * t54 - t187, t14, t1, t8, t15, t9, t40, t41 * t169 - t2 * t130 + t16 * t98 + t27 * t87 + t34 * t81 + t55 * t49 + t6 * t88 + (t155 * t41 + t166) * t101, t41 * t167 + t3 * t130 - t24 * t98 - t33 * t87 + t34 * t83 + t55 * t134 - t7 * t88 + (-t156 * t41 + t165) * t101, -t27 * t134 - t33 * t49 - t6 * t83 - t7 * t81 - (t120 * t24 + t122 * t16) * t97 + (-t3 * t120 - t2 * t122 + (t120 * t16 - t122 * t24) * qJD(4)) * t101, t16 * t6 + t2 * t27 + t24 * t7 + t26 * t55 + t3 * t33 + t41 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -qJ(2) * t142, 0, 0, 0, 0, 0, 0, 0.2e1 * t95 * qJD(3), 0.2e1 * t125, -t91 - t194, t198 * t95 - t204 * t75, 0, 0, 0, 0, 0, 0, t19, t22, t4, t139 * t120 + t200 * t122 - t69 * t95, 0, 0, 0, 0, 0, 0, t19, t22, t4, t140 * t120 + t202 * t122 - t41 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, -t91 + t194, 0, t179, 0, 0, -t205 * t95, -t205 * t204, 0, 0, t17, t5, t21, t18, t20, -t181, -pkin(3) * t49 + t120 * t132 - t122 * t141 - t31 * t95 - t35 * t88 - t75 * t81, -pkin(3) * t134 + t120 * t141 + t122 * t132 + t32 * t95 + t36 * t88 - t75 * t83, t35 * t83 + t36 * t81 + ((qJD(4) * t83 - t49) * pkin(7) + t139) * t122 + ((qJD(4) * t81 + t134) * pkin(7) - t200) * t120, -t44 * pkin(3) - t31 * t35 - t32 * t36 - t69 * t75 + (-qJD(4) * t138 - t11 * t120 - t122 * t145) * pkin(7), t17, t5, t21, t18, t20, -t181, t107 * t87 + t114 * t49 - t165 - t16 * t95 - t45 * t81 - t175 * t88 + (-t41 * t204 + (t41 + t193) * qJD(4)) * t120, -t204 * t168 + t108 * t87 + t114 * t134 + t166 + t24 * t95 - t45 * t83 + t174 * t88 + (pkin(4) * t172 + t168) * qJD(4), -t107 * t134 + t108 * t49 - t202 * t120 + t140 * t122 + t174 * t81 + t175 * t83, t2 * t107 - t3 * t108 + t26 * t114 + (pkin(4) * t156 - t45) * t41 - t174 * t24 - t175 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t37, t29, -t184, t30, t87, -t69 * t83 + t200, t69 * t81 - t139, 0, 0, t184, t37, t29, -t184, t30, t87, 0.2e1 * t192 + t190 + (-t148 - t41) * t83 + t124, -t195 * pkin(4) + t23 * t88 + (qJD(5) + t41) * t81 + t133, -pkin(4) * t134 - t176 * t81, t176 * t24 + (-t41 * t83 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 + t183, t134 - t186, -t80 - t195, t16 * t83 + t24 * t81 + t26;];
tauc_reg = t10;
