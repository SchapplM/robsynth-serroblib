% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:22
% DurationCPUTime: 1.19s
% Computational Cost: add. (3838->212), mult. (9847->331), div. (0->0), fcn. (7883->12), ass. (0->139)
t100 = cos(qJ(3));
t99 = cos(qJ(4));
t140 = qJD(4) * t99;
t95 = sin(qJ(4));
t141 = qJD(4) * t95;
t96 = sin(qJ(3));
t163 = t95 * t96;
t101 = cos(qJ(2));
t97 = sin(qJ(2));
t146 = t100 * t97;
t111 = t101 * t96 + t146;
t91 = sin(pkin(5));
t144 = qJD(1) * t91;
t60 = t111 * t144;
t63 = (t100 * t101 - t96 * t97) * t91;
t61 = qJD(1) * t63;
t82 = t100 * pkin(2) + pkin(3);
t152 = t95 * t60 - t99 * t61 + t82 * t140 + (-t96 * t141 + (t100 * t99 - t163) * qJD(3)) * pkin(2);
t162 = t96 * t99;
t151 = t99 * t60 + t95 * t61 - t82 * t141 + (-t96 * t140 + (-t100 * t95 - t162) * qJD(3)) * pkin(2);
t87 = qJD(2) + qJD(3);
t84 = qJD(4) + t87;
t92 = cos(pkin(6));
t167 = t92 * t84;
t76 = qJD(5) + t167;
t176 = (qJD(5) - t76) * t84;
t93 = cos(pkin(5));
t143 = qJD(1) * t93;
t130 = t97 * t144;
t126 = t101 * t144;
t75 = qJD(2) * pkin(2) + t126;
t57 = t100 * t130 + t96 * t75;
t164 = t95 * t57;
t118 = t96 * t130;
t56 = t100 * t75 - t118;
t53 = t87 * pkin(3) + t56;
t30 = t99 * t53 - t164;
t29 = t84 * pkin(4) + t30;
t90 = sin(pkin(6));
t109 = t90 * t143 + t29 * t92;
t171 = t84 * t90;
t159 = t99 * t57;
t31 = t95 * t53 + t159;
t28 = pkin(10) * t171 + t31;
t94 = sin(qJ(5));
t98 = cos(qJ(5));
t12 = t109 * t94 + t98 * t28;
t64 = t111 * t91;
t44 = t99 * t63 - t95 * t64;
t113 = t44 * t92 + t90 * t93;
t45 = t95 * t63 + t99 * t64;
t24 = t113 * t94 + t98 * t45;
t108 = t111 * qJD(2);
t103 = (-qJD(3) * t146 - t108) * t144;
t142 = qJD(3) * t75;
t37 = -t96 * t142 + t103;
t122 = -t57 * t141 + t95 * t37;
t116 = qJD(2) * t126;
t150 = t87 * t118;
t36 = t100 * (t116 + t142) - t150;
t9 = (qJD(4) * t53 + t36) * t99 + t122;
t123 = -t95 * t36 + t99 * t37;
t10 = -qJD(4) * t31 + t123;
t86 = t90 ^ 2;
t175 = t10 * t86;
t174 = t10 * t94;
t27 = t92 * t143 - t90 * t29;
t173 = t27 * t90;
t172 = t84 * t86;
t170 = t86 * t98;
t169 = t90 * t94;
t168 = t90 * t98;
t166 = t92 * t94;
t165 = t92 * t98;
t71 = pkin(4) * t166 + pkin(10) * t168;
t135 = qJD(5) * t71;
t158 = t31 * t165 + t94 * t30 - t135;
t14 = -t31 * t166 + t98 * t30;
t70 = pkin(4) * t165 - pkin(10) * t169;
t157 = -qJD(5) * t70 + t14;
t149 = pkin(2) * t162 + t95 * t82;
t85 = t90 * pkin(10);
t62 = t85 + t149;
t117 = -pkin(2) * t163 + t99 * t82;
t67 = pkin(4) + t117;
t40 = t67 * t165 - t94 * t62;
t139 = qJD(5) * t40;
t156 = t151 * t166 + t152 * t98 + t139;
t41 = t67 * t166 + t98 * t62;
t138 = qJD(5) * t41;
t155 = t151 * t165 - t152 * t94 - t138;
t77 = t95 * pkin(3) + t85;
t81 = t99 * pkin(3) + pkin(4);
t59 = t81 * t166 + t98 * t77;
t136 = qJD(5) * t59;
t147 = pkin(3) * qJD(4);
t32 = -t95 * t56 - t159;
t33 = t99 * t56 - t164;
t154 = -t32 * t165 + t94 * t33 - t136 + (-t95 * t165 - t94 * t99) * t147;
t58 = t81 * t165 - t94 * t77;
t137 = qJD(5) * t58;
t153 = t32 * t166 + t98 * t33 - t137 - (-t95 * t166 + t98 * t99) * t147;
t148 = t94 ^ 2 - t98 ^ 2;
t145 = qJD(2) ^ 2 * t91;
t134 = qJD(5) * t94;
t133 = qJD(5) * t98;
t128 = t90 * t134;
t8 = t10 * t165;
t4 = -t12 * qJD(5) - t94 * t9 + t8;
t131 = t10 * t170 + t27 * t128 + t4 * t92;
t129 = qJD(5) * t172;
t127 = t90 * t133;
t125 = -pkin(3) * t84 - t53;
t11 = t109 * t98 - t94 * t28;
t3 = qJD(5) * t11 + t10 * t166 + t98 * t9;
t124 = t27 * t127 - t3 * t92;
t121 = t76 + t167;
t83 = t84 ^ 2;
t120 = t94 * t83 * t170;
t115 = pkin(3) * t141 + t32;
t114 = t94 * t98 * t129;
t112 = (-pkin(2) * t87 - t75) * qJD(3);
t23 = t113 * t98 - t94 * t45;
t105 = -t4 * t94 + (-t11 * t98 - t12 * t94) * qJD(5);
t104 = -qJD(5) * t109 - t27 * t171 - t9;
t69 = -0.2e1 * t114;
t68 = 0.2e1 * t114;
t55 = -0.2e1 * t148 * t129;
t51 = t121 * t127;
t50 = t121 * t128;
t47 = (-qJD(3) * t111 - t108) * t91;
t46 = t87 * t63;
t34 = -t90 * t44 + t93 * t92;
t16 = -qJD(4) * t45 - t95 * t46 + t99 * t47;
t15 = qJD(4) * t44 + t99 * t46 + t95 * t47;
t6 = -t24 * qJD(5) - t94 * t15 + t16 * t165;
t5 = qJD(5) * t23 + t98 * t15 + t16 * t166;
t1 = t3 * t168;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * t145, -t101 * t145, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t87, -t46 * t87, 0, t36 * t64 + t37 * t63 + t57 * t46 + t56 * t47, 0, 0, 0, 0, 0, 0, t16 * t84, -t15 * t84, 0, t10 * t44 + t31 * t15 + t30 * t16 + t9 * t45, 0, 0, 0, 0, 0, 0, t6 * t76 + (t34 * t128 + t16 * t170) * t84, -t5 * t76 + (-t16 * t86 * t94 + t34 * t127) * t84, (t5 * t98 - t6 * t94 + (-t23 * t98 - t24 * t94) * qJD(5)) * t171, t11 * t6 + t12 * t5 + t4 * t23 + t3 * t24 + (-t10 * t34 - t16 * t27) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t96 + t60 * t87 + t103, t61 * t87 + (t112 - t116) * t100 + t150, 0, t56 * t60 - t57 * t61 + (t100 * t37 + t36 * t96 + (t100 * t57 - t56 * t96) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t151 * t84 + t10, -t152 * t84 - t9, 0, t10 * t117 + t9 * t149 + t151 * t30 + t152 * t31, t68, t55, t51, t69, -t50, 0, t155 * t76 + (-t67 * t134 + t151 * t98) * t172 + t131, -t156 * t76 + (-t174 + (-t67 * t133 - t151 * t94) * t84) * t86 + t124, t1 + (((-t139 + t156) * t98 + (-t138 - t155) * t94) * t84 + t105) * t90, t155 * t11 + t156 * t12 - t151 * t173 + t67 * t175 + t3 * t41 + t4 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t87 + t37, t56 * t87 - t36, 0, 0, 0, 0, 0, 0, 0, 0, -t32 * t84 + (t125 * t95 - t159) * qJD(4) + t123, t33 * t84 + (t125 * qJD(4) - t36) * t99 - t122, 0, -t30 * t32 - t31 * t33 + (t10 * t99 + t9 * t95 + (-t30 * t95 + t31 * t99) * qJD(4)) * pkin(3), t68, t55, t51, t69, -t50, 0, t154 * t76 + (-t115 * t98 - t81 * t134) * t172 + t131, t153 * t76 + (-t174 + (t115 * t94 - t81 * t133) * t84) * t86 + t124, t1 + (((-t137 - t153) * t98 + (-t136 - t154) * t94) * t84 + t105) * t90, t154 * t11 + t115 * t173 - t153 * t12 + t81 * t175 + t3 * t59 + t4 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t84 + t10, t30 * t84 - t9, 0, 0, t68, t55, t51, t69, -t50, 0, (-pkin(4) * t134 + t31 * t98) * t172 + t158 * t76 + t131, t157 * t76 + (-t174 + (-pkin(4) * t133 - t31 * t94) * t84) * t86 + t124, t1 + ((-t14 * t98 + (-t135 - t158) * t94) * t84 + t105) * t90, pkin(4) * t175 + t158 * t11 - t157 * t12 - t31 * t173 + t3 * t71 + t4 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, t148 * t86 * t83, t168 * t176, t120, -t169 * t176, 0, t104 * t94 + t12 * t76 - t28 * t133 + t8, t11 * t76 + (qJD(5) * t28 - t10 * t92) * t94 + t104 * t98, 0, 0;];
tauc_reg = t2;
