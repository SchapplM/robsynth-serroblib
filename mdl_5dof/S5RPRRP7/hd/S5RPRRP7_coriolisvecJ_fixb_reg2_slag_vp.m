% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:36
% EndTime: 2019-12-31 18:45:42
% DurationCPUTime: 1.91s
% Computational Cost: add. (2105->286), mult. (4995->374), div. (0->0), fcn. (2949->6), ass. (0->159)
t90 = sin(qJ(4));
t159 = qJD(3) * t90;
t91 = sin(qJ(3));
t162 = qJD(1) * t91;
t92 = cos(qJ(4));
t69 = t92 * t162 + t159;
t93 = cos(qJ(3));
t161 = qJD(1) * t93;
t79 = -qJD(4) + t161;
t178 = t69 * t79;
t153 = t92 * qJD(3);
t67 = t90 * t162 - t153;
t180 = t67 * t79;
t151 = qJD(1) * qJD(3);
t133 = t92 * t151;
t156 = qJD(4) * t90;
t139 = t91 * t156;
t150 = qJD(3) * qJD(4);
t39 = qJD(1) * t139 - t93 * t133 - t92 * t150;
t155 = qJD(4) * t92;
t138 = t91 * t155;
t157 = qJD(3) * t93;
t142 = t90 * t157;
t102 = t138 + t142;
t40 = t102 * qJD(1) + t90 * t150;
t204 = (t39 - t180) * t92 + (t40 - t178) * t90;
t154 = t91 * qJD(2);
t81 = sin(pkin(8)) * pkin(1) + pkin(6);
t75 = t81 * qJD(1);
t54 = t93 * t75 + t154;
t44 = qJD(3) * pkin(7) + t54;
t197 = qJD(2) * t93 - t91 * t75;
t45 = t197 * qJD(3);
t82 = -cos(pkin(8)) * pkin(1) - pkin(2);
t64 = -pkin(3) * t93 - pkin(7) * t91 + t82;
t48 = t64 * qJD(1);
t123 = pkin(3) * t91 - pkin(7) * t93;
t73 = t123 * qJD(3);
t63 = qJD(1) * t73;
t130 = t44 * t155 + t48 * t156 + t90 * t45 - t92 * t63;
t19 = t44 * t92 + t48 * t90;
t108 = -t19 * t79 - t130;
t140 = t79 * t156;
t175 = t90 * t93;
t201 = ((t67 + t153) * t91 - t79 * t175) * qJD(1) + t140;
t14 = -qJ(5) * t79 + t19;
t83 = t91 * t151;
t126 = pkin(4) * t83;
t2 = -t126 + t130;
t200 = t14 * t79 + t2;
t171 = t92 * t93;
t199 = t81 * t171 + t90 * t64;
t198 = t40 + t178;
t195 = qJD(4) * t199 - t92 * t73;
t194 = t69 ^ 2;
t193 = pkin(7) * t69;
t192 = pkin(7) * t79;
t46 = t54 * qJD(3);
t5 = t40 * pkin(4) + t39 * qJ(5) - t69 * qJD(5) + t46;
t191 = t5 * t90;
t190 = t5 * t92;
t43 = -qJD(3) * pkin(3) - t197;
t17 = pkin(4) * t67 - qJ(5) * t69 + t43;
t188 = t17 * t69;
t186 = t40 * t92;
t185 = t43 * t90;
t184 = t43 * t92;
t183 = t46 * t90;
t182 = t46 * t91;
t181 = t46 * t92;
t179 = t69 * t67;
t177 = t81 * t90;
t176 = t90 * t91;
t174 = t91 * t92;
t173 = t92 * t64;
t170 = t93 * t39;
t169 = t93 * t40;
t94 = qJD(3) ^ 2;
t168 = t94 * t91;
t167 = t94 * t93;
t121 = pkin(4) * t90 - qJ(5) * t92;
t166 = t154 + (t121 * qJD(1) + t75) * t93 - t121 * qJD(4) + t90 * qJD(5);
t141 = t93 * t153;
t30 = t40 * t174;
t165 = -t67 * t141 - t30;
t72 = t123 * qJD(1);
t24 = t197 * t92 + t90 * t72;
t164 = t64 * t155 + t90 * t73;
t86 = t91 ^ 2;
t163 = -t93 ^ 2 + t86;
t76 = qJD(1) * t82;
t158 = qJD(3) * t91;
t18 = -t44 * t90 + t48 * t92;
t152 = qJD(5) - t18;
t149 = t90 * t192;
t148 = t92 * t192;
t95 = qJD(1) ^ 2;
t146 = t91 * t95 * t93;
t51 = t69 * t138;
t145 = t69 * t142 - t39 * t176 + t51;
t144 = pkin(7) * t158;
t143 = pkin(7) * t153;
t137 = t79 * t155;
t136 = t79 * t162;
t135 = t67 ^ 2 - t194;
t134 = pkin(4) + t177;
t77 = t86 * t133;
t131 = t77 - t170;
t128 = qJD(4) * t67 - t39;
t127 = t91 * t137;
t125 = t93 * t83;
t124 = qJ(5) * t83;
t122 = pkin(4) * t92 + qJ(5) * t90;
t13 = pkin(4) * t79 + t152;
t120 = t13 * t92 - t14 * t90;
t119 = t13 * t90 + t14 * t92;
t118 = -t18 * t92 - t19 * t90;
t117 = t18 * t90 - t19 * t92;
t23 = -t197 * t90 + t72 * t92;
t116 = t128 * pkin(7);
t115 = qJD(1) * t86 - t79 * t93;
t113 = 0.2e1 * qJD(3) * t76;
t109 = t121 + t81;
t107 = t67 * t139 + t165;
t106 = -t48 * t155 + t44 * t156 - t92 * t45 - t90 * t63;
t105 = t115 * t90;
t104 = t127 + t169;
t103 = -t139 + t141;
t101 = -t180 * t90 - t186;
t1 = -qJD(5) * t79 - t106 + t124;
t100 = t120 * qJD(4) + t1 * t92 + t2 * t90;
t99 = t118 * qJD(4) - t106 * t92 + t130 * t90;
t98 = t45 * t93 + t182 + (-t197 * t93 - t54 * t91) * qJD(3);
t97 = t102 * t67 + t40 * t176;
t59 = t67 * t158;
t96 = -qJD(3) * t105 + t127 - t169 + t59;
t74 = -pkin(3) - t122;
t65 = t79 * t141;
t58 = t69 * t158;
t47 = (-t79 - t161) * t158;
t41 = t109 * t91;
t33 = pkin(7) * t186;
t32 = pkin(4) * t69 + qJ(5) * t67;
t28 = -t81 * t175 + t173;
t27 = t134 * t93 - t173;
t26 = -qJ(5) * t93 + t199;
t22 = -pkin(4) * t162 - t23;
t21 = qJ(5) * t162 + t24;
t20 = -t39 - t180;
t16 = (t122 * qJD(4) - qJD(5) * t92) * t91 + t109 * t157;
t15 = -t137 + (t79 * t171 + (-t69 + t159) * t91) * qJD(1);
t12 = t158 * t177 - t195;
t11 = (-t91 * t153 - t93 * t156) * t81 + t164;
t10 = -t178 * t92 - t39 * t90;
t9 = -t134 * t158 + t195;
t8 = t103 * t69 - t39 * t174;
t7 = (-t81 * t156 - qJD(5)) * t93 + (-t81 * t92 + qJ(5)) * t158 + t164;
t6 = t79 * t139 + t170 + t58 - t65 + t77;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t125, -0.2e1 * t163 * t151, t167, -0.2e1 * t125, -t168, 0, t91 * t113 - t81 * t167, t93 * t113 + t81 * t168, t98, t98 * t81, t8, t107 - t145, t6, t97, (-t67 * t91 - t105) * qJD(3) + t104, t47, -t12 * t79 + (t130 + (t67 * t81 + t185) * qJD(3)) * t93 + (t43 * t155 + t40 * t81 + t183 + (qJD(1) * t28 + t18) * qJD(3)) * t91, t11 * t79 + (-t106 + (t69 * t81 + t184) * qJD(3)) * t93 + (-t43 * t156 - t39 * t81 + t181 + (-qJD(1) * t199 - t19) * qJD(3)) * t91, -t11 * t67 - t12 * t69 + t28 * t39 - t199 * t40 + t118 * t157 + (qJD(4) * t117 + t106 * t90 + t130 * t92) * t91, t19 * t11 + t18 * t12 - t130 * t28 - t106 * t199 + (t43 * t157 + t182) * t81, t8, t6, t103 * t67 + t145 + t30, t47, t115 * t159 - t104 + t59, t97, t16 * t67 + t41 * t40 + t9 * t79 + (t17 * t159 + t2) * t93 + (t17 * t155 + t191 + (-qJD(1) * t27 - t13) * qJD(3)) * t91, -t26 * t40 - t27 * t39 - t67 * t7 + t69 * t9 + t120 * t157 + (-qJD(4) * t119 - t1 * t90 + t2 * t92) * t91, -t16 * t69 + t41 * t39 - t7 * t79 + (-t17 * t153 - t1) * t93 + (t17 * t156 - t190 + (qJD(1) * t26 + t14) * qJD(3)) * t91, t1 * t26 + t13 * t9 + t14 * t7 + t16 * t17 + t2 * t27 + t41 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, -t167, 0, t45 * t91 - t46 * t93 + (-t197 * t91 + t54 * t93) * qJD(3), 0, 0, 0, 0, 0, 0, t96, t103 * t79 - t131 + t58, t51 + (t128 * t91 + t69 * t157) * t90 + t165, (-qJD(3) * t117 - t46) * t93 + (qJD(3) * t43 + t99) * t91, 0, 0, 0, 0, 0, 0, t96, t107 + t145, -t65 + (-qJD(3) * t69 + t140) * t91 + t131, (qJD(3) * t119 - t5) * t93 + (qJD(3) * t17 + t100) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, t163 * t95, 0, t146, 0, 0, -t76 * t162, -t76 * t161, 0, 0, t10, -t204, t15, t101, t201, t136, -pkin(3) * t40 + t23 * t79 - t181 - t54 * t67 + (t148 + t185) * qJD(4) + (-t18 * t91 + (-t43 * t93 - t144) * t90) * qJD(1), pkin(3) * t39 - t24 * t79 + t183 - t54 * t69 + (-t149 + t184) * qJD(4) + (-t43 * t171 + (t19 - t143) * t91) * qJD(1), t23 * t69 + t24 * t67 - t33 + (t18 * t161 - t106 + (-t18 + t193) * qJD(4)) * t92 + (t116 - t108) * t90, -t46 * pkin(3) + pkin(7) * t99 - t18 * t23 - t19 * t24 - t43 * t54, t10, t15, t204, t136, -t201, t101, -t22 * t79 + t40 * t74 - t190 - t166 * t67 + (t17 * t90 + t148) * qJD(4) + (t13 * t91 + (-t17 * t93 - t144) * t90) * qJD(1), t21 * t67 - t22 * t69 - t33 + (-t13 * t161 + t1 + (t13 + t193) * qJD(4)) * t92 + (t116 + t200) * t90, t21 * t79 + t39 * t74 - t191 + t166 * t69 + (-t17 * t92 + t149) * qJD(4) + (t17 * t171 + (-t14 + t143) * t91) * qJD(1), pkin(7) * t100 - t13 * t22 - t14 * t21 - t166 * t17 + t5 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t135, t20, -t179, -t198, t83, -t43 * t69 + t108, -t18 * t79 + t43 * t67 + t106, 0, 0, t179, t20, t135, t83, t198, -t179, -t32 * t67 + t108 + 0.2e1 * t126 - t188, pkin(4) * t39 - t40 * qJ(5) + (t14 - t19) * t69 + (t13 - t152) * t67, 0.2e1 * t124 - t17 * t67 + t32 * t69 + (-0.2e1 * qJD(5) + t18) * t79 - t106, -t2 * pkin(4) + t1 * qJ(5) - t13 * t19 + t14 * t152 - t17 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83 + t179, t20, -t79 ^ 2 - t194, t188 + t200;];
tauc_reg = t3;
