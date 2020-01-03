% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:12
% EndTime: 2019-12-31 18:22:17
% DurationCPUTime: 1.84s
% Computational Cost: add. (2652->267), mult. (6421->404), div. (0->0), fcn. (4237->8), ass. (0->157)
t124 = sin(qJ(5));
t126 = cos(qJ(5));
t120 = sin(pkin(9));
t122 = cos(pkin(9));
t173 = qJD(3) * t122;
t125 = sin(qJ(3));
t176 = qJD(1) * t125;
t90 = -t120 * t176 + t173;
t165 = t122 * t176;
t174 = qJD(3) * t120;
t91 = t165 + t174;
t42 = t124 * t91 - t126 * t90;
t203 = t42 ^ 2;
t45 = t124 * t90 + t126 * t91;
t202 = t45 ^ 2;
t127 = cos(qJ(3));
t175 = qJD(1) * t127;
t110 = -qJD(5) + t175;
t201 = t110 * t42;
t170 = qJD(5) * t126;
t186 = t120 * t124;
t200 = -qJD(5) * t186 + t122 * t170;
t168 = qJD(1) * qJD(3);
t158 = t127 * t168;
t151 = t120 * t158;
t171 = qJD(3) * t127;
t162 = t122 * t171;
t152 = t126 * t162;
t18 = (qJD(5) * t91 + t151) * t124 - qJD(1) * t152 - t90 * t170;
t112 = sin(pkin(8)) * pkin(1) + pkin(6);
t104 = t112 * qJD(1);
t169 = t125 * qJD(2);
t78 = t104 * t127 + t169;
t62 = qJD(3) * qJ(4) + t78;
t113 = -cos(pkin(8)) * pkin(1) - pkin(2);
t88 = -pkin(3) * t127 - qJ(4) * t125 + t113;
t65 = t88 * qJD(1);
t23 = -t120 * t62 + t122 * t65;
t14 = -pkin(4) * t175 - pkin(7) * t91 + t23;
t24 = t120 * t65 + t122 * t62;
t16 = pkin(7) * t90 + t24;
t145 = t124 * t16 - t126 * t14;
t182 = t122 * t127;
t140 = pkin(4) * t125 - pkin(7) * t182;
t135 = t140 * qJD(3);
t115 = qJD(2) * t171;
t95 = t125 * t104;
t58 = t115 + (qJD(4) - t95) * qJD(3);
t148 = pkin(3) * t125 - qJ(4) * t127;
t82 = t148 * qJD(3) - t125 * qJD(4);
t66 = t82 * qJD(1);
t20 = -t120 * t58 + t122 * t66;
t15 = qJD(1) * t135 + t20;
t21 = t120 * t66 + t122 * t58;
t17 = -pkin(7) * t151 + t21;
t1 = -t145 * qJD(5) + t124 * t15 + t126 * t17;
t77 = qJD(2) * t127 - t95;
t97 = t148 * qJD(1);
t36 = -t120 * t77 + t122 * t97;
t22 = t140 * qJD(1) + t36;
t166 = t120 * t175;
t37 = t120 * t97 + t122 * t77;
t27 = -pkin(7) * t166 + t37;
t196 = pkin(7) + qJ(4);
t101 = t196 * t120;
t102 = t196 * t122;
t50 = -t101 * t126 - t102 * t124;
t93 = -t126 * t122 + t186;
t199 = -t93 * qJD(4) + t50 * qJD(5) - t124 * t22 - t126 * t27;
t51 = -t101 * t124 + t102 * t126;
t94 = t120 * t126 + t122 * t124;
t198 = -t94 * qJD(4) - t51 * qJD(5) + t124 * t27 - t126 * t22;
t197 = t45 * t42;
t137 = t94 * t127;
t133 = qJD(3) * t137;
t19 = qJD(1) * t133 + t45 * qJD(5);
t163 = t120 * t171;
t85 = t94 * qJD(5);
t32 = t124 * t163 + t125 * t85 - t152;
t74 = t93 * t125;
t195 = t74 * t19 + t32 * t42;
t157 = t125 * t168;
t33 = t200 * t125 + t133;
t73 = t94 * t125;
t194 = t33 * t110 - t73 * t157;
t193 = qJD(1) * t137 - t85;
t192 = -t93 * t175 - t200;
t172 = qJD(3) * t125;
t164 = t112 * t172;
t39 = t120 * t164 + t122 * t82;
t96 = t112 * t182;
t47 = t120 * t88 + t96;
t190 = t125 * t90;
t68 = t78 * qJD(3);
t189 = t68 * t125;
t188 = t68 * t127;
t119 = t127 ^ 2;
t129 = qJD(1) ^ 2;
t187 = t119 * t129;
t185 = t120 * t125;
t184 = t120 * t127;
t183 = t122 * t125;
t128 = qJD(3) ^ 2;
t181 = t128 * t125;
t180 = t128 * t127;
t118 = t125 ^ 2;
t179 = t118 - 0.2e1 * t119;
t178 = t118 - t119;
t105 = qJD(1) * t113;
t177 = qJD(1) * t120;
t167 = t112 * t184;
t159 = t122 * t168;
t156 = pkin(4) * t120 + t112;
t155 = t127 * t18 + t45 * t172;
t154 = -qJD(3) * pkin(3) + qJD(4);
t153 = -t90 + t173;
t150 = t127 * t157;
t149 = -t18 * t73 + t33 * t45;
t147 = -t120 * t20 + t122 * t21;
t146 = -t120 * t23 + t122 * t24;
t6 = t124 * t14 + t126 * t16;
t76 = t122 * t88;
t31 = -pkin(7) * t183 + t76 + (-t112 * t120 - pkin(4)) * t127;
t35 = -pkin(7) * t185 + t47;
t9 = -t124 * t35 + t126 * t31;
t10 = t124 * t31 + t126 * t35;
t143 = qJD(1) * (-t91 + t174);
t142 = qJD(1) * t153;
t141 = 0.2e1 * qJD(3) * t105;
t139 = t127 * t19 - t42 * t172;
t138 = t127 * t143;
t136 = -t110 * t32 + t74 * t157;
t60 = t154 - t77;
t48 = pkin(4) * t151 + t68;
t132 = -qJ(4) * t172 + (t154 - t60) * t127;
t2 = -t6 * qJD(5) - t124 * t17 + t126 * t15;
t67 = -t104 * t172 + t115;
t130 = t189 + t67 * t127 + (-t125 * t78 - t127 * t77) * qJD(3);
t117 = t122 ^ 2;
t116 = t120 ^ 2;
t114 = -pkin(4) * t122 - pkin(3);
t109 = t125 * t129 * t127;
t106 = -0.2e1 * t150;
t81 = t156 * t125;
t80 = t91 * t172;
t72 = t156 * t171;
t70 = t120 * t82;
t69 = t90 * t162;
t52 = t169 + (pkin(4) * t177 + t104) * t127;
t46 = t76 - t167;
t40 = -t122 * t164 + t70;
t38 = -pkin(4) * t90 + t60;
t30 = t70 + (-pkin(7) * t184 - t112 * t183) * qJD(3);
t26 = t135 + t39;
t4 = -t10 * qJD(5) - t124 * t30 + t126 * t26;
t3 = t9 * qJD(5) + t124 * t26 + t126 * t30;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t150, -0.2e1 * t178 * t168, t180, t106, -t181, 0, -t112 * t180 + t125 * t141, t112 * t181 + t127 * t141, t130, t130 * t112, (t117 * t176 + t122 * t91) * t171, t69 + (-t91 - 0.2e1 * t165) * t163, t179 * t159 + t80, (t116 * t176 - t120 * t90) * t171, (-t179 * t177 + t190) * qJD(3), t106, t68 * t185 + (-qJD(1) * t39 - t20) * t127 + ((-t112 * t90 + t120 * t60) * t127 + (t23 + (t46 + t167) * qJD(1)) * t125) * qJD(3), t68 * t183 + (qJD(1) * t40 + t21) * t127 + ((t112 * t91 + t122 * t60) * t127 + (-t24 + (-t47 + t96) * qJD(1)) * t125) * qJD(3), -t39 * t91 + t40 * t90 + (-t120 * t21 - t122 * t20) * t125 + (-t120 * t24 - t122 * t23 + (-t120 * t47 - t122 * t46) * qJD(1)) * t171, t20 * t46 + t21 * t47 + t23 * t39 + t24 * t40 + (t60 * t171 + t189) * t112, t18 * t74 - t32 * t45, -t149 + t195, -t136 + t155, t19 * t73 + t33 * t42, t139 + t194, (-t110 - t175) * t172, -t110 * t4 - t127 * t2 + t19 * t81 + t33 * t38 + t42 * t72 + t48 * t73 + (qJD(1) * t9 - t145) * t172, t1 * t127 + t110 * t3 - t18 * t81 - t32 * t38 + t45 * t72 - t48 * t74 + (-qJD(1) * t10 - t6) * t172, -t1 * t73 - t10 * t19 - t145 * t32 + t18 * t9 + t2 * t74 - t3 * t42 - t33 * t6 - t4 * t45, t1 * t10 - t145 * t4 + t2 * t9 + t3 * t6 + t38 * t72 + t48 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -t180, 0, t67 * t125 - t188 + (-t125 * t77 + t127 * t78) * qJD(3), 0, 0, 0, 0, 0, 0, (-t118 * t177 - t190) * qJD(3), -t118 * t159 + t80, t91 * t163 + t69, -t188 + t147 * t125 + (t125 * t60 + t146 * t127) * qJD(3), 0, 0, 0, 0, 0, 0, -t139 + t194, t136 + t155, t149 + t195, -t1 * t74 - t127 * t48 + t145 * t33 + t172 * t38 - t2 * t73 - t32 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t178 * t129, 0, t109, 0, 0, -t105 * t176, -t105 * t175 - t115 + (t77 + t95) * qJD(3), 0, 0, t122 * t138, (t120 * t91 - t122 * t90 + (-t116 + t117) * qJD(3)) * t175, t122 * t187 + t125 * t143, -t153 * t166, -t120 * t187 + t125 * t142, t109, -t122 * t68 + t78 * t90 + (t132 * t120 - t125 * t23 + t127 * t36) * qJD(1), t120 * t68 - t78 * t91 + (t132 * t122 + t125 * t24 - t127 * t37) * qJD(1), t36 * t91 - t37 * t90 + (qJD(4) * t90 + t23 * t175 + t21) * t122 + (qJD(4) * t91 + t24 * t175 - t20) * t120, -pkin(3) * t68 + t147 * qJ(4) + t146 * qJD(4) - t23 * t36 - t24 * t37 - t60 * t78, -t18 * t94 - t192 * t45, t18 * t93 - t19 * t94 + t192 * t42 + t193 * t45, t192 * t110 + (qJD(3) * t94 - t45) * t176, t19 * t93 - t193 * t42, -t193 * t110 + (-qJD(3) * t93 + t42) * t176, t110 * t176, t114 * t19 - t42 * t52 + t48 * t93 - t193 * t38 - t198 * t110 + (qJD(3) * t50 + t145) * t176, -t114 * t18 - t45 * t52 + t48 * t94 - t192 * t38 + t199 * t110 + (-qJD(3) * t51 + t6) * t176, -t1 * t93 - t145 * t192 + t18 * t50 - t19 * t51 + t193 * t6 - t198 * t45 - t199 * t42 - t2 * t94, t1 * t51 + t114 * t48 - t145 * t198 + t199 * t6 + t2 * t50 - t38 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, t127 * t142, -t90 ^ 2 - t91 ^ 2, t23 * t91 - t24 * t90 + t68, 0, 0, 0, 0, 0, 0, -t110 * t45 + t19, -t18 + t201, -t202 - t203, -t145 * t45 + t42 * t6 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t202 - t203, -t18 - t201, -t197, -t158 * t94 + (-qJD(5) - t110) * t45, t157, -t110 * t6 - t38 * t45 + t2, t110 * t145 + t38 * t42 - t1, 0, 0;];
tauc_reg = t5;
