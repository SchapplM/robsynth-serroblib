% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:48
% EndTime: 2021-01-15 22:24:55
% DurationCPUTime: 1.38s
% Computational Cost: add. (3265->241), mult. (8683->310), div. (0->0), fcn. (6035->6), ass. (0->152)
t117 = qJD(2) + qJD(3);
t120 = sin(pkin(8));
t121 = cos(pkin(8));
t124 = cos(qJ(3));
t125 = cos(qJ(2));
t164 = qJD(1) * t125;
t156 = t124 * t164;
t122 = sin(qJ(3));
t123 = sin(qJ(2));
t165 = qJD(1) * t123;
t157 = t122 * t165;
t80 = -t156 + t157;
t82 = -t122 * t164 - t124 * t165;
t52 = t120 * t82 - t121 * t80;
t196 = t117 * t52;
t143 = -t120 * t80 - t121 * t82;
t195 = t143 ^ 2;
t161 = qJD(1) * qJD(2);
t194 = -0.2e1 * t161;
t191 = pkin(6) + pkin(7);
t100 = t191 * t125;
t95 = qJD(1) * t100;
t172 = t124 * t95;
t99 = t191 * t123;
t93 = qJD(1) * t99;
t141 = t122 * t93 - t172;
t175 = qJ(4) * t80;
t134 = t141 + t175;
t160 = pkin(2) * t120 * t122;
t162 = qJD(3) * t124;
t83 = t122 * t95;
t176 = -t124 * t93 - t83;
t77 = t82 * qJ(4);
t44 = t77 + t176;
t177 = qJD(3) * t160 + t120 * t134 + (-pkin(2) * t162 + t44) * t121;
t113 = -pkin(2) * t125 - pkin(1);
t98 = t113 * qJD(1);
t63 = pkin(3) * t80 + qJD(4) + t98;
t22 = -pkin(4) * t52 - qJ(5) * t143 + t63;
t185 = t22 * t143;
t174 = qJD(2) * pkin(2);
t87 = -t93 + t174;
t153 = t124 * t87 - t83;
t42 = t153 + t77;
t92 = t122 * t125 + t123 * t124;
t193 = qJD(1) * t92;
t159 = qJD(2) * t191;
t147 = qJD(1) * t159;
t88 = t123 * t147;
t192 = (qJD(3) * t87 - t88) * t124;
t190 = pkin(3) * t82;
t163 = qJD(3) * t122;
t89 = t125 * t147;
t151 = -t122 * t89 - t95 * t163;
t62 = t117 * t92;
t57 = t62 * qJD(1);
t13 = -qJ(4) * t57 - qJD(4) * t80 + t151 + t192;
t142 = -t122 * t87 - t172;
t152 = t122 * t88 - t124 * t89;
t130 = t142 * qJD(3) + t152;
t154 = t125 * t161;
t56 = qJD(3) * t156 - t117 * t157 + t124 * t154;
t14 = -qJ(4) * t56 + qJD(4) * t82 + t130;
t2 = t120 * t13 - t121 * t14;
t132 = -qJ(4) * t92 - t100 * t122 - t124 * t99;
t140 = -t100 * t124 + t122 * t99;
t91 = t122 * t123 - t124 * t125;
t48 = -qJ(4) * t91 - t140;
t29 = t120 * t48 - t121 * t132;
t189 = t2 * t29;
t94 = t123 * t159;
t96 = t125 * t159;
t129 = t140 * qJD(3) + t122 * t94 - t124 * t96;
t61 = t117 * t91;
t128 = qJ(4) * t61 - qJD(4) * t92 + t129;
t136 = -t100 * t163 - t122 * t96 - t124 * t94 - t99 * t162;
t23 = -qJ(4) * t62 - qJD(4) * t91 + t136;
t6 = t120 * t23 - t121 * t128;
t188 = t117 * t6;
t7 = t120 * t128 + t121 * t23;
t187 = t117 * t7;
t43 = -t142 - t175;
t39 = t121 * t43;
t19 = t120 * t42 + t39;
t186 = t19 * t143;
t184 = t22 * t52;
t183 = t143 * t63;
t182 = t52 * t63;
t181 = t82 * t80;
t180 = t98 * t82;
t3 = t120 * t14 + t121 * t13;
t171 = t121 * t122;
t179 = t120 * t44 - t121 * t134 - (t120 * t124 + t171) * qJD(3) * pkin(2);
t178 = -qJD(5) + t177;
t36 = pkin(3) * t117 + t42;
t18 = t120 * t36 + t39;
t173 = t120 * t43;
t112 = pkin(2) * t124 + pkin(3);
t76 = pkin(2) * t171 + t120 * t112;
t127 = qJD(1) ^ 2;
t170 = t125 * t127;
t126 = qJD(2) ^ 2;
t169 = t126 * t123;
t168 = t126 * t125;
t20 = t121 * t42 - t173;
t167 = qJD(5) - t20;
t166 = t123 ^ 2 - t125 ^ 2;
t116 = t117 * qJD(5);
t1 = t116 + t3;
t115 = t123 * t174;
t114 = pkin(2) * t165;
t158 = t179 * t143;
t155 = -pkin(2) * t117 - t87;
t46 = pkin(3) * t57 + qJD(2) * t114;
t55 = pkin(3) * t62 + t115;
t31 = t120 * t56 + t121 * t57;
t149 = pkin(1) * t194;
t148 = t117 * t20 - t3;
t17 = t121 * t36 - t173;
t15 = -pkin(4) * t117 + qJD(5) - t17;
t16 = qJ(5) * t117 + t18;
t146 = t143 * t16 - t15 * t52;
t145 = t143 * t18 + t17 * t52;
t144 = -t52 ^ 2 - t195;
t32 = -t120 * t57 + t121 * t56;
t65 = pkin(3) * t91 + t113;
t139 = t98 * t80 - t151;
t138 = t117 * t19 - t2;
t137 = t117 * t143 + t31;
t75 = t112 * t121 - t160;
t27 = pkin(4) * t143 - qJ(5) * t52 - t190;
t135 = t179 * t117 - t2;
t133 = t32 + t196;
t30 = t120 * t132 + t121 * t48;
t60 = -t120 * t91 + t121 * t92;
t131 = t143 * t6 + t2 * t60 + t29 * t32 - t30 * t31 + t52 * t7;
t5 = pkin(4) * t31 - qJ(5) * t32 - qJD(5) * t143 + t46;
t110 = -pkin(3) * t121 - pkin(4);
t109 = pkin(3) * t120 + qJ(5);
t70 = -pkin(4) - t75;
t69 = qJ(5) + t76;
t64 = t114 - t190;
t59 = t120 * t92 + t121 * t91;
t45 = -t80 ^ 2 + t82 ^ 2;
t38 = (-t193 - t82) * t117;
t37 = t117 * t80 + t56;
t34 = -t120 * t62 - t121 * t61;
t33 = -t120 * t61 + t121 * t62;
t28 = pkin(4) * t59 - qJ(5) * t60 + t65;
t26 = t114 + t27;
t8 = pkin(4) * t33 - qJ(5) * t34 - qJD(5) * t60 + t55;
t4 = [0, 0, 0, 0.2e1 * t123 * t154, t166 * t194, t168, -t169, 0, -pkin(6) * t168 + t123 * t149, pkin(6) * t169 + t125 * t149, t56 * t92 + t61 * t82, -t56 * t91 - t57 * t92 + t61 * t80 + t62 * t82, -t61 * t117, -t62 * t117, 0, t113 * t57 + t98 * t62 + t129 * t117 + (qJD(1) * t91 + t80) * t115, t113 * t56 - t98 * t61 - t136 * t117 + (-t82 + t193) * t115, t31 * t65 + t33 * t63 + t46 * t59 - t52 * t55 - t188, t143 * t55 + t32 * t65 + t34 * t63 + t46 * t60 - t187, -t17 * t34 - t18 * t33 - t3 * t59 + t131, -t17 * t6 + t18 * t7 + t3 * t30 + t46 * t65 + t55 * t63 + t189, t22 * t33 + t28 * t31 + t5 * t59 - t52 * t8 - t188, -t1 * t59 + t15 * t34 - t16 * t33 + t131, -t143 * t8 - t22 * t34 - t28 * t32 - t5 * t60 + t187, t1 * t30 + t15 * t6 + t16 * t7 + t22 * t8 + t28 * t5 + t189; 0, 0, 0, -t123 * t170, t166 * t127, 0, 0, 0, t127 * pkin(1) * t123, pkin(1) * t170, -t181, t45, t37, t38, 0, -t80 * t114 + t180 - t141 * t117 + (t155 * t122 - t172) * qJD(3) + t152, t82 * t114 + t176 * t117 + (t155 * qJD(3) + t88) * t124 + t139, t52 * t64 + t135 - t183, t177 * t117 - t143 * t64 - t182 - t3, -t177 * t52 - t31 * t76 - t32 * t75 + t145 - t158, t179 * t17 - t177 * t18 - t2 * t75 + t3 * t76 - t63 * t64, t26 * t52 + t135 - t185, -t178 * t52 - t31 * t69 + t32 * t70 + t146 - t158, -t178 * t117 + t143 * t26 + t1 + t184, t1 * t69 - t179 * t15 - t178 * t16 + t2 * t70 - t22 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, t45, t37, t38, 0, -t142 * t117 + t130 + t180, t153 * t117 + t139 - t192, -t190 * t52 + t138 - t183, t143 * t190 + t148 - t182, -t186 - t20 * t52 + (-t120 * t31 - t121 * t32) * pkin(3) + t145, t17 * t19 - t18 * t20 + (t120 * t3 - t121 * t2 + t63 * t82) * pkin(3), t27 * t52 + t138 - t185, -t109 * t31 + t110 * t32 + t167 * t52 + t146 - t186, t143 * t27 + 0.2e1 * t116 - t148 + t184, t1 * t109 + t110 * t2 - t15 * t19 + t167 * t16 - t22 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t133, t144, t143 * t17 - t18 * t52 + t46, t137, t144, -t133, -t143 * t15 - t16 * t52 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143 * t52, t32 - t196, -t117 ^ 2 - t195, -t117 * t16 + t185 + t2;];
tauc_reg = t4;
