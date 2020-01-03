% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:59
% EndTime: 2019-12-31 19:48:04
% DurationCPUTime: 1.61s
% Computational Cost: add. (1379->248), mult. (3355->365), div. (0->0), fcn. (2089->6), ass. (0->153)
t118 = sin(qJ(2));
t163 = t118 * qJD(1);
t101 = qJD(5) + t163;
t196 = qJD(5) - t101;
t190 = pkin(3) + pkin(6);
t114 = sin(pkin(8));
t115 = cos(pkin(8));
t117 = sin(qJ(5));
t119 = cos(qJ(5));
t76 = t119 * t114 + t117 * t115;
t127 = t76 * t118;
t64 = t76 * qJD(5);
t188 = -qJD(1) * t127 - t64;
t120 = cos(qJ(2));
t170 = qJD(1) * t120;
t154 = t115 * t170;
t165 = t114 * qJD(2);
t72 = t154 + t165;
t155 = t114 * t170;
t164 = t115 * qJD(2);
t74 = -t155 + t164;
t136 = t117 * t72 - t119 * t74;
t195 = t101 * t136;
t161 = qJD(1) * qJD(2);
t194 = -0.2e1 * t161;
t103 = pkin(6) * t163;
t193 = qJD(3) + t103;
t111 = qJD(2) * qJ(3);
t192 = qJD(4) + t111;
t116 = -pkin(2) - qJ(4);
t168 = qJD(2) * t120;
t104 = pkin(6) * t170;
t105 = pkin(3) * t170;
t83 = t104 + t105;
t62 = t83 + t192;
t191 = t118 * (-t62 + t192) - t116 * t168;
t152 = t118 * t161;
t144 = t119 * t152;
t145 = t117 * t152;
t12 = -t136 * qJD(5) + t114 * t145 - t115 * t144;
t189 = -pkin(7) + t116;
t100 = pkin(2) * t152;
t134 = -qJ(3) * t120 + qJ(4) * t118;
t162 = t118 * qJD(3);
t124 = t134 * qJD(2) - t120 * qJD(4) - t162;
t27 = t124 * qJD(1) + t100;
t151 = t120 * t161;
t99 = pkin(6) * t151;
t57 = t99 + (-qJD(4) + t105) * qJD(2);
t10 = t114 * t57 + t115 * t27;
t169 = qJD(2) * t118;
t106 = pkin(2) * t169;
t36 = t106 + t124;
t84 = t190 * t168;
t18 = t114 * t84 + t115 * t36;
t150 = -t118 * qJ(3) - pkin(1);
t71 = t116 * t120 + t150;
t46 = t71 * qJD(1);
t172 = pkin(3) * t163 + t193;
t51 = t116 * qJD(2) + t172;
t16 = t114 * t51 + t115 * t46;
t107 = pkin(2) * t163;
t58 = t134 * qJD(1) + t107;
t26 = t114 * t83 + t115 * t58;
t95 = t190 * t118;
t32 = t114 * t95 + t115 * t71;
t157 = t115 * t163;
t166 = qJD(5) * t119;
t167 = qJD(5) * t117;
t177 = t117 * t114;
t187 = -t114 * t167 + t115 * t166 + t119 * t157 - t163 * t177;
t186 = qJD(2) * pkin(2);
t15 = -t114 * t46 + t115 * t51;
t185 = t120 * t15;
t184 = t120 * t16;
t28 = t117 * t74 + t119 * t72;
t183 = t28 * t101;
t110 = qJD(2) * qJD(3);
t82 = t190 * t169;
t56 = -qJD(1) * t82 + t110;
t182 = t56 * t114;
t92 = -t120 * pkin(2) + t150;
t70 = qJD(1) * t92;
t112 = t118 ^ 2;
t122 = qJD(1) ^ 2;
t181 = t112 * t122;
t180 = t114 * t118;
t179 = t115 * t118;
t178 = t115 * t120;
t176 = t120 * t122;
t121 = qJD(2) ^ 2;
t175 = t121 * t118;
t174 = t121 * t120;
t159 = -pkin(4) * t115 - pkin(3);
t173 = -t159 * t163 + t193;
t96 = t190 * t120;
t171 = -t120 ^ 2 + t112;
t160 = t118 * t176;
t132 = pkin(4) * t120 - pkin(7) * t180;
t125 = t132 * qJD(2);
t9 = -t114 * t27 + t115 * t57;
t4 = qJD(1) * t125 + t9;
t147 = pkin(7) * t157;
t5 = qJD(2) * t147 + t10;
t158 = -t117 * t5 + t119 * t4;
t153 = t187 * t101;
t17 = -t114 * t36 + t115 * t84;
t25 = -t114 * t58 + t115 * t83;
t149 = pkin(1) * t194;
t148 = qJD(3) - t186;
t133 = -t119 * t115 + t177;
t146 = t188 * t101 - t133 * t151;
t143 = t117 * t4 + t119 * t5;
t6 = pkin(4) * t163 - t74 * pkin(7) + t15;
t7 = -t72 * pkin(7) + t16;
t142 = t117 * t7 - t119 * t6;
t2 = t117 * t6 + t119 * t7;
t141 = t10 * t114 + t9 * t115;
t139 = -t114 * t15 + t115 * t16;
t79 = t115 * t95;
t20 = t118 * pkin(4) + t79 + (pkin(7) * t120 - t71) * t114;
t23 = -pkin(7) * t178 + t32;
t138 = -t117 * t23 + t119 * t20;
t137 = t117 * t20 + t119 * t23;
t135 = -0.2e1 * qJD(2) * t70;
t128 = -qJ(3) * t168 - t162;
t47 = t128 * qJD(1) + t100;
t60 = t106 + t128;
t131 = pkin(6) * t121 + qJD(1) * t60 + t47;
t86 = t189 * t115;
t130 = qJD(4) * t114 - qJD(5) * t86 + t147 + t26;
t85 = t189 * t114;
t129 = t132 * qJD(1) + qJD(4) * t115 + qJD(5) * t85 + t25;
t11 = t114 * t144 + t115 * t145 - t72 * t166 - t74 * t167;
t49 = (-pkin(6) + t159) * t169;
t53 = t133 * t120;
t90 = pkin(6) * t152 - t110;
t91 = t103 + t148;
t94 = -t104 - t111;
t123 = -t90 * t120 + (t120 * t91 + (t94 + t104) * t118) * qJD(2);
t102 = t114 * pkin(4) + qJ(3);
t80 = -qJ(3) * t170 + t107;
t63 = pkin(4) * t178 + t96;
t54 = t76 * t120;
t50 = t70 * t163;
t35 = qJD(1) * t49 + t110;
t33 = t72 * pkin(4) + t62;
t31 = -t114 * t71 + t79;
t22 = t120 * t64 - t133 * t169;
t21 = qJD(2) * t127 + qJD(5) * t53;
t13 = t118 * pkin(7) * t164 + t18;
t8 = t125 + t17;
t1 = [0, 0, 0, 0.2e1 * t118 * t151, t171 * t194, t174, -t175, 0, -pkin(6) * t174 + t118 * t149, pkin(6) * t175 + t120 * t149, t123, t118 * t135 + t131 * t120, -t131 * t118 + t120 * t135, t123 * pkin(6) + t47 * t92 + t70 * t60, t56 * t178 - t82 * t72 + (qJD(1) * t17 + t9) * t118 + (-t62 * t179 + t185 + (t120 * t31 - t96 * t179) * qJD(1)) * qJD(2), -t120 * t182 - t82 * t74 + (-qJD(1) * t18 - t10) * t118 + (t62 * t180 - t184 + (-t120 * t32 + t96 * t180) * qJD(1)) * qJD(2), -t17 * t74 - t18 * t72 + (-t10 * t115 + t114 * t9) * t120 + ((-t114 * t31 + t115 * t32) * qJD(1) + t139) * t169, t10 * t32 + t15 * t17 + t16 * t18 + t9 * t31 + t56 * t96 - t62 * t82, -t11 * t54 - t136 * t21, t11 * t53 + t54 * t12 - t136 * t22 - t21 * t28, t21 * t101 + t11 * t118 + (-qJD(1) * t54 - t136) * t168, t22 * t101 - t12 * t118 + (qJD(1) * t53 - t28) * t168, (t101 + t163) * t168, (-t117 * t13 + t119 * t8) * t101 + t158 * t118 + t49 * t28 + t63 * t12 - t35 * t53 - t33 * t22 + (-t101 * t137 - t118 * t2) * qJD(5) + (qJD(1) * t138 - t142) * t168, -(t117 * t8 + t119 * t13) * t101 - t143 * t118 - t49 * t136 + t63 * t11 - t35 * t54 + t33 * t21 + (-t101 * t138 + t118 * t142) * qJD(5) + (-qJD(1) * t137 - t2) * t168; 0, 0, 0, -t160, t171 * t122, 0, 0, 0, t122 * pkin(1) * t118, pkin(1) * t176, ((-t94 - t111) * t118 + (t148 - t91) * t120) * qJD(1), -t170 * t80 + t50, 0.2e1 * t110 + (t118 * t80 + t120 * t70) * qJD(1), -t90 * qJ(3) - t94 * qJD(3) - t70 * t80 + (-t118 * t94 + (-t91 - t186) * t120) * qJD(1) * pkin(6), t182 + t172 * t72 + (-t115 * t191 - t118 * t25 - t185) * qJD(1), t56 * t115 + t172 * t74 + (t114 * t191 + t118 * t26 + t184) * qJD(1), t25 * t74 + t26 * t72 + (qJD(4) * t74 - t16 * t163 - t9) * t115 + (qJD(4) * t72 + t15 * t163 - t10) * t114, t56 * qJ(3) - t15 * t25 - t16 * t26 + t172 * t62 + t141 * t116 + (-t114 * t16 - t115 * t15) * qJD(4), -t11 * t133 - t136 * t188, -t11 * t76 + t12 * t133 + t136 * t187 - t188 * t28, t136 * t170 + t146, -t153 + (-qJD(2) * t76 + t28) * t170, -t101 * t170, t102 * t12 + t35 * t76 + t187 * t33 + t173 * t28 + (t117 * t130 - t119 * t129) * t101 + ((-t117 * t85 + t119 * t86) * qJD(2) + t142) * t170, t102 * t11 - t35 * t133 + t188 * t33 - t173 * t136 + (t117 * t129 + t119 * t130) * t101 + (-(t117 * t86 + t119 * t85) * qJD(2) + t2) * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, -t121 - t181, qJD(2) * t94 + t50 + t99, -t114 * t181 + (-t72 + t154) * qJD(2), -t115 * t181 + (-t74 - t155) * qJD(2), (t114 * t74 - t115 * t72) * t163, -t62 * qJD(2) + t139 * t163 + t141, 0, 0, 0, 0, 0, -qJD(2) * t28 + t146, -t153 + (-t170 * t76 + t136) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t74 - t164) * t163, (-t72 + t165) * t163, -t72 ^ 2 - t74 ^ 2, t15 * t74 + t16 * t72 + t56, 0, 0, 0, 0, 0, t12 - t195, t11 - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136 * t28, t136 ^ 2 - t28 ^ 2, t11 + t183, -t12 - t195, t151, t33 * t136 - t196 * t2 + t158, t142 * t196 + t33 * t28 - t143;];
tauc_reg = t1;
