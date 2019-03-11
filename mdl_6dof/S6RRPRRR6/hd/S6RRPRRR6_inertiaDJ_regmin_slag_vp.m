% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:25
% EndTime: 2019-03-09 13:53:31
% DurationCPUTime: 1.62s
% Computational Cost: add. (2564->208), mult. (5500->321), div. (0->0), fcn. (5264->8), ass. (0->134)
t164 = sin(qJ(4));
t166 = cos(qJ(4));
t95 = sin(qJ(2));
t97 = cos(qJ(2));
t114 = t164 * t95 + t166 * t97;
t101 = (qJD(2) - qJD(4)) * t114;
t93 = sin(qJ(6));
t91 = t93 ^ 2;
t96 = cos(qJ(6));
t92 = t96 ^ 2;
t174 = (t91 - t92) * qJD(6);
t175 = -pkin(2) * t97 - qJ(3) * t95;
t168 = pkin(7) - pkin(8);
t72 = t168 * t95;
t172 = qJD(4) + qJD(5);
t98 = -pkin(2) - pkin(3);
t67 = qJ(3) * t166 + t164 * t98;
t171 = qJD(3) * t164 + qJD(4) * t67;
t170 = -0.2e1 * t174;
t169 = 0.2e1 * qJD(3);
t165 = cos(qJ(5));
t73 = t168 * t97;
t113 = -t164 * t73 + t166 * t72;
t139 = t97 * t164;
t140 = t95 * t166;
t64 = t140 - t139;
t34 = -pkin(9) * t64 + t113;
t112 = t164 * t72 + t166 * t73;
t35 = -pkin(9) * t114 + t112;
t94 = sin(qJ(5));
t22 = -t165 * t34 + t35 * t94;
t86 = t97 * qJD(2);
t84 = pkin(7) * t86;
t127 = pkin(8) * t86 - t84;
t115 = t164 * t127;
t149 = t95 * qJD(2);
t131 = t164 * qJD(2);
t134 = qJD(4) * t166;
t155 = -t131 * t97 - t134 * t95;
t100 = -t115 + t155 * pkin(9) + (pkin(9) * t139 + t113) * qJD(4) + (pkin(9) - t168) * t149 * t166;
t23 = t165 * t35 + t34 * t94;
t102 = t127 * t166 - t131 * t72;
t133 = qJD(4) * t164;
t99 = pkin(9) * t101 + t133 * t72 + t134 * t73 + t102;
t5 = qJD(5) * t23 + t100 * t94 + t165 * t99;
t87 = qJD(6) * t96;
t167 = t22 * t87 + t5 * t93;
t110 = t94 * t114;
t44 = t165 * t64 - t110;
t163 = t44 * t93;
t162 = t44 * t96;
t116 = t133 * t97 + t155;
t128 = qJD(2) * t140;
t105 = -t128 - t116;
t106 = t165 * t114;
t152 = qJD(5) * t94;
t14 = qJD(5) * t106 - t101 * t165 + t105 * t94 + t152 * t64;
t161 = t92 * t14;
t160 = t94 * t67;
t159 = t96 * t14;
t132 = qJD(5) * t165;
t15 = -qJD(5) * t110 + t101 * t94 + t105 * t165 + t132 * t64;
t158 = t96 * t15;
t138 = t165 * t67;
t111 = qJ(3) * t164 - t166 * t98;
t66 = -pkin(4) - t111;
t117 = t66 * t94 + t138;
t54 = -qJD(3) * t166 + qJD(4) * t111;
t135 = t165 * t171 - t94 * t54;
t26 = qJD(5) * t117 + t135;
t45 = -t165 * t66 + pkin(5) + t160;
t157 = t26 * t93 + t45 * t87;
t144 = pkin(4) * t152;
t142 = t165 * pkin(4);
t83 = -t142 - pkin(5);
t156 = t144 * t93 + t83 * t87;
t154 = qJ(3) * t86 + qJD(3) * t95;
t150 = qJD(6) * t93;
t148 = -0.2e1 * pkin(1) * qJD(2);
t70 = -pkin(1) + t175;
t147 = -t132 * t66 + t165 * t54 + t171 * t94;
t146 = pkin(5) * t150;
t145 = pkin(5) * t87;
t143 = pkin(7) * t149;
t141 = t93 * t87;
t137 = 0.4e1 * t93 * t162;
t37 = t45 * t150;
t136 = t26 * t96 - t37;
t60 = pkin(3) * t97 - t70;
t129 = pkin(4) * t132;
t126 = t165 * t166;
t43 = t64 * t94 + t106;
t49 = pkin(4) * t114 + t60;
t16 = pkin(5) * t43 - pkin(10) * t44 + t49;
t125 = t16 * t96 - t23 * t93;
t124 = t16 * t93 + t23 * t96;
t46 = -pkin(10) + t117;
t123 = t43 * t46 - t44 * t45;
t63 = t164 * t94 - t126;
t65 = t164 * t165 + t166 * t94;
t122 = t43 * t65 - t44 * t63;
t82 = pkin(4) * t94 + pkin(10);
t121 = t43 * t82 - t44 * t83;
t68 = t83 * t150;
t120 = -t144 * t96 + t68;
t119 = -t14 * t93 + t44 * t87;
t118 = t150 * t44 + t159;
t9 = t150 * t43 - t158;
t25 = t152 * t67 + t147;
t109 = -t14 * t45 - t15 * t46 + t25 * t43 + t26 * t44;
t47 = (qJD(5) * t164 + t133) * t94 - t172 * t126;
t48 = t172 * t65;
t108 = -t14 * t63 - t15 * t65 + t43 * t47 + t44 * t48;
t107 = qJD(2) * t175 + qJD(3) * t97;
t103 = -t14 * t83 - t15 * t82 + (-t165 * t43 + t44 * t94) * qJD(5) * pkin(4);
t29 = -t116 * pkin(4) + (-pkin(4) * t166 + t98) * t149 + t154;
t76 = -0.2e1 * t141;
t75 = 0.2e1 * t141;
t57 = pkin(2) * t149 - t154;
t50 = t98 * t149 + t154;
t42 = t44 ^ 2;
t31 = -t150 * t63 + t48 * t96;
t30 = t48 * t93 + t63 * t87;
t28 = qJD(4) * t112 + t102;
t27 = -t113 * qJD(4) + t128 * t168 + t115;
t20 = t22 * t150;
t10 = t15 * t93 + t43 * t87;
t8 = t93 * t159 + t174 * t44;
t7 = qJD(6) * t137 - t14 * t91 + t161;
t6 = t15 * pkin(5) + t14 * pkin(10) + t29;
t4 = -t100 * t165 - t132 * t34 + t152 * t35 + t94 * t99;
t2 = -qJD(6) * t124 + t4 * t93 + t6 * t96;
t1 = -qJD(6) * t125 + t4 * t96 - t6 * t93;
t3 = [0, 0, 0, 0.2e1 * t95 * t86, 0.2e1 * (-t95 ^ 2 + t97 ^ 2) * qJD(2), 0, 0, 0, t95 * t148, t97 * t148, 0.2e1 * t149 * t70 - 0.2e1 * t57 * t97, 0, -0.2e1 * t57 * t95 - 0.2e1 * t70 * t86, 0.2e1 * t70 * t57, 0.2e1 * t64 * t101, -0.2e1 * t101 * t114 - 0.2e1 * t105 * t64, 0, 0, 0, 0.2e1 * t105 * t60 + 0.2e1 * t114 * t50, 0.2e1 * t101 * t60 + 0.2e1 * t50 * t64, -0.2e1 * t44 * t14, 0.2e1 * t14 * t43 - 0.2e1 * t15 * t44, 0, 0, 0, 0.2e1 * t15 * t49 + 0.2e1 * t29 * t43, -0.2e1 * t14 * t49 + 0.2e1 * t29 * t44, -0.2e1 * t141 * t42 - 0.2e1 * t161 * t44, t137 * t14 + 0.2e1 * t174 * t42, -0.2e1 * t118 * t43 + 0.2e1 * t158 * t44, -0.2e1 * t119 * t43 - 0.2e1 * t15 * t163, 0.2e1 * t43 * t15, 0.2e1 * t119 * t22 + 0.2e1 * t125 * t15 + 0.2e1 * t163 * t5 + 0.2e1 * t2 * t43, 0.2e1 * t1 * t43 - 0.2e1 * t118 * t22 - 0.2e1 * t124 * t15 + 0.2e1 * t162 * t5; 0, 0, 0, 0, 0, t86, -t149, 0, -t84, t143, -t84, t107, -t143, t107 * pkin(7), 0, 0, -t101, t105, 0, t28, -t27, 0, 0, t14, t15, 0, t5, -t4, t8, t7, -t10, t9, 0, -t20 + (-qJD(6) * t123 + t5) * t96 + t109 * t93, t109 * t96 + t123 * t150 - t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, qJ(3) * t169, 0, 0, 0, 0, 0, 0.2e1 * t171, -0.2e1 * t54, 0, 0, 0, 0, 0, 0.2e1 * t26, -0.2e1 * t25, t75, t170, 0, 0, 0, 0.2e1 * t136, -0.2e1 * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t93 - t122 * t87, t108 * t96 + t122 * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t134, 0, 0, 0, 0, 0, t48, -t47, 0, 0, 0, 0, 0, t31, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t105, 0, -t28, t27, 0, 0, -t14, -t15, 0, -t5, t4, -t8, -t7, t10, -t9, 0, t20 + (-qJD(6) * t121 - t5) * t96 + t103 * t93, t103 * t96 + t121 * t150 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, t54, 0, 0, 0, 0, 0 (-t138 + (pkin(4) - t66) * t94) * qJD(5) - t135 (t142 + t160) * qJD(5) + t147, t76, -t170, 0, 0, 0, t37 - t68 + (-t26 + t144) * t96, -t156 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t134, 0, 0, 0, 0, 0, -t48, t47, 0, 0, 0, 0, 0, -t31, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t144, -0.2e1 * t129, t75, t170, 0, 0, 0, 0.2e1 * t120, 0.2e1 * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, -t5, t4, -t8, -t7, t10, -t9, 0, t20 + (pkin(5) * t14 - pkin(10) * t15) * t93 + (-t5 + (-pkin(5) * t44 - pkin(10) * t43) * qJD(6)) * t96, pkin(5) * t118 + pkin(10) * t9 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, t76, -t170, 0, 0, 0, -t136 + t146, t145 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0, 0, 0, 0, -t31, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, -t129, t75, t170, 0, 0, 0, t120 - t146, -t145 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t170, 0, 0, 0, -0.2e1 * t146, -0.2e1 * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, -t119, t15, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t150, 0, t25 * t93 - t46 * t87, t150 * t46 + t25 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t93 - t65 * t87, t150 * t65 + t47 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t150, 0, -t129 * t93 - t82 * t87, -t129 * t96 + t150 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t150, 0, -pkin(10) * t87, pkin(10) * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
