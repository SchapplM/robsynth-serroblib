% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:06
% EndTime: 2019-03-09 18:04:11
% DurationCPUTime: 1.60s
% Computational Cost: add. (5458->198), mult. (12024->329), div. (0->0), fcn. (12300->10), ass. (0->135)
t159 = qJD(2) + qJD(3);
t95 = sin(qJ(3));
t96 = sin(qJ(2));
t98 = cos(qJ(3));
t99 = cos(qJ(2));
t74 = t95 * t96 - t98 * t99;
t54 = t159 * t74;
t75 = t95 * t99 + t98 * t96;
t161 = t54 * qJ(4) - t75 * qJD(4);
t158 = pkin(7) + pkin(8);
t80 = t158 * t96;
t81 = t158 * t99;
t112 = -t98 * t80 - t95 * t81;
t160 = t112 * qJD(3);
t97 = cos(qJ(6));
t90 = t97 ^ 2;
t93 = sin(qJ(6));
t145 = t93 ^ 2 - t90;
t121 = t145 * qJD(6);
t128 = qJD(2) * t158;
t76 = t96 * t128;
t77 = t99 * t128;
t113 = t98 * t76 + t95 * t77;
t55 = t159 * t75;
t101 = -t55 * qJ(4) - t74 * qJD(4) - t113;
t111 = t95 * t80 - t98 * t81;
t123 = t95 * t76 - t98 * t77;
t91 = sin(pkin(11));
t92 = cos(pkin(11));
t20 = -t91 * t101 + t92 * (t123 + t161) + (t92 * t111 - t91 * t112) * qJD(3);
t38 = -t92 * t54 - t91 * t55;
t100 = t38 * pkin(9) - t20;
t115 = -t91 * t54 + t92 * t55;
t40 = t111 * qJD(3) + t123;
t21 = t92 * (t101 + t160) + t91 * (t40 + t161);
t104 = -t115 * pkin(9) + t21;
t155 = cos(qJ(5));
t49 = -t75 * qJ(4) + t112;
t50 = -t74 * qJ(4) - t111;
t30 = t92 * t49 - t91 * t50;
t53 = -t91 * t74 + t92 * t75;
t24 = -t53 * pkin(9) + t30;
t114 = t92 * t74 + t91 * t75;
t31 = t91 * t49 + t92 * t50;
t25 = -t114 * pkin(9) + t31;
t94 = sin(qJ(5));
t17 = t155 * t25 + t94 * t24;
t5 = t17 * qJD(5) + t155 * t100 + t94 * t104;
t3 = t5 * t93;
t157 = t91 * pkin(3);
t16 = -t155 * t24 + t94 * t25;
t88 = qJD(6) * t97;
t156 = t16 * t88 + t3;
t35 = t155 * t114 + t94 * t53;
t18 = -t35 * qJD(5) - t94 * t115 + t155 * t38;
t36 = -t94 * t114 + t155 * t53;
t154 = t36 * t18;
t153 = t36 * t97;
t152 = t91 * t95;
t151 = t92 * t95;
t19 = t36 * qJD(5) + t155 * t115 + t94 * t38;
t150 = t93 * t19;
t149 = t97 * t18;
t148 = t97 * t19;
t85 = t98 * pkin(2) + pkin(3);
t71 = pkin(2) * t151 + t91 * t85;
t129 = t155 * t71;
t70 = -pkin(2) * t152 + t92 * t85;
t65 = pkin(4) + t70;
t106 = t94 * t65 + t129;
t144 = pkin(2) * qJD(3);
t68 = (-t91 * t98 - t151) * t144;
t69 = (t92 * t98 - t152) * t144;
t124 = t155 * t68 - t94 * t69;
t33 = t106 * qJD(5) - t124;
t44 = -t155 * t65 + t94 * t71 - pkin(5);
t147 = t33 * t93 + t44 * t88;
t120 = t155 * t157;
t83 = t92 * pkin(3) + pkin(4);
t105 = t94 * t83 + t120;
t64 = t105 * qJD(5);
t137 = t94 * t157;
t66 = -t155 * t83 - pkin(5) + t137;
t146 = t64 * t93 + t66 * t88;
t142 = qJD(2) * t96;
t141 = qJD(2) * t99;
t140 = qJD(5) * t94;
t139 = qJD(6) * t93;
t136 = -0.2e1 * pkin(1) * qJD(2);
t122 = qJD(5) * t155;
t135 = -t65 * t122 - t155 * t69 - t94 * t68;
t134 = pkin(5) * t139;
t133 = pkin(5) * t88;
t87 = pkin(2) * t142;
t132 = t95 * t144;
t131 = t98 * t144;
t130 = t93 * t88;
t86 = -t99 * pkin(2) - pkin(1);
t51 = t55 * pkin(3) + t87;
t127 = -0.4e1 * t93 * t153;
t42 = t44 * t139;
t126 = -t33 * t97 + t42;
t57 = t66 * t139;
t125 = -t64 * t97 + t57;
t110 = t74 * pkin(3) + t86;
t41 = t114 * pkin(4) + t110;
t22 = t35 * pkin(5) - t36 * pkin(10) + t41;
t119 = t97 * t17 + t93 * t22;
t118 = t93 * t17 - t97 * t22;
t45 = pkin(10) + t106;
t117 = t35 * t45 - t36 * t44;
t67 = pkin(10) + t105;
t116 = t35 * t67 - t36 * t66;
t109 = t93 * t18 + t36 * t88;
t108 = t36 * t139 - t149;
t10 = t35 * t88 + t150;
t107 = t35 * t139 - t148;
t32 = t71 * t140 + t135;
t103 = t18 * t44 - t19 * t45 + t32 * t35 + t33 * t36;
t78 = t83 * t122;
t63 = qJD(5) * t137 - t78;
t102 = t18 * t66 - t19 * t67 + t35 * t63 + t36 * t64;
t28 = t115 * pkin(4) + t51;
t82 = 0.2e1 * t130;
t73 = -0.2e1 * t121;
t39 = t113 - t160;
t34 = t36 ^ 2;
t14 = t16 * t139;
t8 = -t36 * t121 + t93 * t149;
t7 = qJD(6) * t127 - t145 * t18;
t6 = t19 * pkin(5) - t18 * pkin(10) + t28;
t4 = t94 * t100 - t155 * t104 - t24 * t122 + t25 * t140;
t2 = -t119 * qJD(6) + t93 * t4 + t97 * t6;
t1 = t118 * qJD(6) + t97 * t4 - t93 * t6;
t9 = [0, 0, 0, 0.2e1 * t96 * t141, 0.2e1 * (-t96 ^ 2 + t99 ^ 2) * qJD(2), 0, 0, 0, t96 * t136, t99 * t136, -0.2e1 * t75 * t54, 0.2e1 * t54 * t74 - 0.2e1 * t75 * t55, 0, 0, 0, 0.2e1 * t86 * t55 + 0.2e1 * t74 * t87, -0.2e1 * t86 * t54 + 0.2e1 * t75 * t87, -0.2e1 * t21 * t114 - 0.2e1 * t31 * t115 - 0.2e1 * t20 * t53 - 0.2e1 * t30 * t38, 0.2e1 * t110 * t51 + 0.2e1 * t30 * t20 + 0.2e1 * t31 * t21, 0.2e1 * t154, -0.2e1 * t18 * t35 - 0.2e1 * t36 * t19, 0, 0, 0, 0.2e1 * t41 * t19 + 0.2e1 * t28 * t35, 0.2e1 * t41 * t18 + 0.2e1 * t28 * t36, -0.2e1 * t34 * t130 + 0.2e1 * t90 * t154, 0.2e1 * t34 * t121 + t18 * t127, -0.2e1 * t108 * t35 + 0.2e1 * t36 * t148, -0.2e1 * t109 * t35 - 0.2e1 * t36 * t150, 0.2e1 * t35 * t19, 0.2e1 * t109 * t16 - 0.2e1 * t118 * t19 + 0.2e1 * t2 * t35 + 0.2e1 * t36 * t3, 0.2e1 * t1 * t35 - 0.2e1 * t108 * t16 - 0.2e1 * t119 * t19 + 0.2e1 * t5 * t153; 0, 0, 0, 0, 0, t141, -t142, 0, -pkin(7) * t141, pkin(7) * t142, 0, 0, -t54, -t55, 0, t40, t39, -t69 * t114 - t71 * t115 - t70 * t38 - t68 * t53, t20 * t70 + t21 * t71 + t30 * t68 + t31 * t69, 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t107, 0, t14 + (-t117 * qJD(6) - t5) * t97 + t103 * t93, t103 * t97 + t117 * t139 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t132, -0.2e1 * t131, 0, 0.2e1 * t70 * t68 + 0.2e1 * t71 * t69, 0, 0, 0, 0, 0, -0.2e1 * t33, 0.2e1 * t32, t82, t73, 0, 0, 0, 0.2e1 * t126, 0.2e1 * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55, 0, t40, t39 (-t91 * t115 - t92 * t38) * pkin(3) (t20 * t92 + t21 * t91) * pkin(3), 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t107, 0, t14 + (-t116 * qJD(6) - t5) * t97 + t102 * t93, t102 * t97 + t116 * t139 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t131, 0 (t68 * t92 + t69 * t91) * pkin(3), 0, 0, 0, 0, 0 (-t120 - t129 + (-t65 - t83) * t94) * qJD(5) + t124, -t78 + (t71 + t157) * t140 + t135, t82, t73, 0, 0, 0, t42 + t57 + (-t33 - t64) * t97, t146 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t64, 0.2e1 * t63, t82, t73, 0, 0, 0, 0.2e1 * t125, 0.2e1 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, t19, t18, 0, 0, 0, 0, 0, -t107, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t107, 0, t14 + (-pkin(5) * t18 - pkin(10) * t19) * t93 + (-t5 + (-pkin(5) * t36 - pkin(10) * t35) * qJD(6)) * t97, t108 * pkin(5) + t107 * pkin(10) + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, t82, t73, 0, 0, 0, t126 - t134, -t133 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t63, t82, t73, 0, 0, 0, t125 - t134, -t133 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t73, 0, 0, 0, -0.2e1 * t134, -0.2e1 * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t109, t19, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t139, 0, t32 * t93 - t45 * t88, t45 * t139 + t32 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t139, 0, t63 * t93 - t67 * t88, t67 * t139 + t63 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t139, 0, -pkin(10) * t88, pkin(10) * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
