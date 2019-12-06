% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:41
% EndTime: 2019-12-05 17:16:49
% DurationCPUTime: 1.79s
% Computational Cost: add. (1739->210), mult. (4474->333), div. (0->0), fcn. (3430->10), ass. (0->139)
t178 = cos(qJ(4));
t131 = qJD(2) * t178;
t94 = sin(qJ(3));
t151 = qJD(2) * t94;
t93 = sin(qJ(4));
t97 = cos(qJ(3));
t182 = t97 * t131 - t93 * t151;
t87 = qJD(3) + qJD(4);
t38 = t182 * t87;
t58 = qJD(5) - t182;
t181 = qJD(5) - t58;
t179 = pkin(7) + pkin(8);
t90 = sin(pkin(5));
t154 = qJD(1) * t90;
t95 = sin(qJ(2));
t140 = t95 * t154;
t129 = t179 * qJD(2) + t140;
t91 = cos(pkin(5));
t153 = qJD(1) * t91;
t48 = -t129 * t94 + t97 * t153;
t156 = qJD(3) * pkin(3);
t104 = t94 * t156 - t140;
t106 = t178 * t97 - t93 * t94;
t49 = t129 * t97 + t94 * t153;
t165 = t93 * t49;
t43 = t48 + t156;
t21 = t178 * t43 - t165;
t16 = -t87 * pkin(4) - t21;
t71 = t179 * t94;
t72 = t179 * t97;
t107 = -t178 * t71 - t93 * t72;
t98 = cos(qJ(2));
t139 = t98 * t154;
t135 = qJD(3) * t179;
t67 = t94 * t135;
t68 = t97 * t135;
t160 = -t107 * qJD(4) + t106 * t139 + t178 * t67 + t93 * t68;
t86 = -t97 * pkin(3) - pkin(2);
t57 = t86 * qJD(2) - t139;
t164 = t93 * t97;
t64 = -qJD(2) * t164 - t94 * t131;
t29 = -pkin(4) * t182 + t64 * pkin(9) + t57;
t130 = qJD(4) * t178;
t149 = qJD(4) * t93;
t152 = qJD(2) * t90;
t134 = qJD(1) * t152;
t122 = t98 * t134;
t32 = t48 * qJD(3) + t97 * t122;
t33 = -t49 * qJD(3) - t94 * t122;
t3 = t43 * t130 - t49 * t149 + t178 * t32 + t93 * t33;
t66 = t178 * t94 + t164;
t45 = t87 * t66;
t39 = t45 * qJD(2);
t132 = -t178 * t33 + t93 * t32;
t136 = t178 * t49;
t22 = t93 * t43 + t136;
t4 = t22 * qJD(4) + t132;
t41 = -pkin(4) * t106 - t66 * pkin(9) + t86;
t44 = t87 * t106;
t56 = t178 * t72 - t93 * t71;
t180 = (qJD(5) * t29 + t3) * t106 + t16 * t44 + t4 * t66 + (-qJD(5) * t41 + t160) * t58 - t56 * t39;
t92 = sin(qJ(5));
t96 = cos(qJ(5));
t116 = t96 * t64 - t92 * t87;
t20 = -t116 * qJD(5) + t92 * t38;
t177 = t16 * t182;
t176 = t16 * t66;
t147 = qJD(5) * t96;
t148 = qJD(5) * t92;
t19 = t87 * t147 + t64 * t148 + t96 * t38;
t175 = t19 * t92;
t174 = t41 * t39;
t52 = -t92 * t64 - t96 * t87;
t173 = t52 * t58;
t172 = t116 * t58;
t171 = t58 * t64;
t170 = t64 * t182;
t169 = t90 * t95;
t168 = t90 * t98;
t166 = t92 * t39;
t163 = t96 * t39;
t99 = qJD(3) ^ 2;
t162 = t99 * t94;
t161 = t99 * t97;
t159 = t56 * qJD(4) - t66 * t139 + t178 * t68 - t93 * t67;
t146 = qJD(2) * qJD(3);
t133 = t94 * t146;
t59 = pkin(3) * t133 + t95 * t134;
t158 = t94 ^ 2 - t97 ^ 2;
t157 = qJD(2) * pkin(2);
t100 = qJD(2) ^ 2;
t155 = t100 * t90;
t150 = qJD(2) * t95;
t144 = t95 * t155;
t143 = pkin(3) * t151;
t142 = t90 * t150;
t141 = t98 * t152;
t127 = t58 * t96;
t40 = -t64 * pkin(4) - pkin(9) * t182;
t84 = t93 * pkin(3) + pkin(9);
t126 = qJD(5) * t84 + t143 + t40;
t125 = t94 * t141;
t124 = t97 * t141;
t17 = t87 * pkin(9) + t22;
t9 = t96 * t17 + t92 * t29;
t123 = t16 * t147 + t4 * t92 - t9 * t64;
t23 = t93 * t48 + t136;
t121 = pkin(3) * t149 - t23;
t120 = t45 * pkin(4) - t44 * pkin(9) + t104;
t118 = -t84 * t39 - t177;
t117 = t92 * t17 - t96 * t29;
t114 = -t117 * t64 + t16 * t148 - t4 * t96;
t24 = t178 * t48 - t165;
t113 = -pkin(3) * t130 + t24;
t60 = -t94 * t169 + t91 * t97;
t61 = t97 * t169 + t91 * t94;
t36 = t178 * t61 + t93 * t60;
t112 = -t96 * t168 - t92 * t36;
t111 = t92 * t168 - t96 * t36;
t110 = t57 * t64 - t132;
t109 = -t66 * t148 + t96 * t44;
t108 = t178 * t60 - t93 * t61;
t105 = t157 * qJD(2);
t103 = -0.2e1 * qJD(3) * t157;
t101 = -t182 * t57 - t3;
t85 = -t178 * pkin(3) - pkin(4);
t47 = -t61 * qJD(3) - t125;
t46 = t60 * qJD(3) + t124;
t30 = -t182 ^ 2 + t64 ^ 2;
t28 = (-qJD(2) * t66 - t64) * t87;
t13 = t39 * pkin(4) - t38 * pkin(9) + t59;
t12 = t96 * t13;
t11 = t36 * qJD(4) - t178 * t47 + t93 * t46;
t10 = t108 * qJD(4) + t178 * t46 + t93 * t47;
t7 = -t116 * t64 + t58 * t127 + t166;
t6 = -t58 ^ 2 * t92 - t52 * t64 + t163;
t5 = -t116 * t127 + t175;
t1 = (t19 - t173) * t96 + (-t20 + t172) * t92;
t2 = [0, 0, -t144, -t98 * t155, 0, 0, 0, 0, 0, -t97 * t144 + (t47 - t125) * qJD(3), t94 * t144 + (-t46 - t124) * qJD(3), 0, 0, 0, 0, 0, -t11 * t87 + (-t150 * t182 - t39 * t98) * t90, -t10 * t87 + (-t64 * t150 - t38 * t98) * t90, 0, 0, 0, 0, 0, (qJD(5) * t111 - t92 * t10 + t142 * t96) * t58 + t112 * t39 + t11 * t52 - t108 * t20, -(qJD(5) * t112 + t96 * t10 + t142 * t92) * t58 + t111 * t39 - t11 * t116 - t108 * t19; 0, 0, 0, 0, 0.2e1 * t97 * t133, -0.2e1 * t158 * t146, t161, -t162, 0, -pkin(7) * t161 + t94 * t103, pkin(7) * t162 + t103 * t97, t38 * t66 - t64 * t44, t106 * t38 + t182 * t44 - t66 * t39 + t64 * t45, t44 * t87, -t45 * t87, 0, -t104 * t182 - t106 * t59 - t159 * t87 + t86 * t39 + t57 * t45, -t104 * t64 + t160 * t87 + t86 * t38 + t57 * t44 + t59 * t66, t19 * t96 * t66 - t109 * t116, (t116 * t92 - t52 * t96) * t44 + (-t175 - t20 * t96 + (t116 * t96 + t52 * t92) * qJD(5)) * t66, -t106 * t19 + t109 * t58 - t116 * t45 + t66 * t163, -t66 * t166 + t20 * t106 - t52 * t45 + (-t66 * t147 - t92 * t44) * t58, -t106 * t39 + t58 * t45, -t12 * t106 - t107 * t20 - t117 * t45 + t159 * t52 + (t174 + t120 * t58 + (t106 * t17 - t56 * t58 + t176) * qJD(5)) * t96 + t180 * t92, -t107 * t19 - t9 * t45 - t159 * t116 + (-t174 + (-qJD(5) * t17 + t13) * t106 - qJD(5) * t176 + (qJD(5) * t56 - t120) * t58) * t92 + t180 * t96; 0, 0, 0, 0, -t94 * t100 * t97, t158 * t100, 0, 0, 0, t94 * t105, t97 * t105, t170, t30, 0, t28, 0, t182 * t143 + t23 * t87 + (-t136 + (-pkin(3) * t87 - t43) * t93) * qJD(4) + t110, t24 * t87 + (-t130 * t87 + t64 * t151) * pkin(3) + t101, t5, t1, t7, t6, t171, t85 * t20 + t118 * t92 + t121 * t52 + (t113 * t92 - t126 * t96) * t58 + t114, t85 * t19 + t118 * t96 - t121 * t116 + (t113 * t96 + t126 * t92) * t58 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t30, 0, t28, 0, t110 + (-qJD(4) + t87) * t22, t21 * t87 + t101, t5, t1, t7, t6, t171, -pkin(4) * t20 - (-t92 * t21 + t96 * t40) * t58 - t22 * t52 - t92 * t177 + (-t58 * t147 - t166) * pkin(9) + t114, -pkin(4) * t19 + (t96 * t21 + t92 * t40) * t58 + t22 * t116 - t96 * t177 + (t58 * t148 - t163) * pkin(9) + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 * t52, t116 ^ 2 - t52 ^ 2, t19 + t173, -t20 - t172, t39, t16 * t116 - t181 * t9 - t92 * t3 + t12, t181 * t117 - t92 * t13 + t16 * t52 - t96 * t3;];
tauc_reg = t2;
