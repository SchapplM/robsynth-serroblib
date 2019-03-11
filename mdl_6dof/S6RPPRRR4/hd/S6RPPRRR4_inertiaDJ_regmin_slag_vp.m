% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:57
% EndTime: 2019-03-09 02:27:00
% DurationCPUTime: 1.21s
% Computational Cost: add. (1032->168), mult. (2300->313), div. (0->0), fcn. (2037->8), ass. (0->120)
t66 = sin(qJ(4));
t116 = t66 * qJD(4);
t62 = sin(pkin(10));
t63 = cos(pkin(10));
t70 = -pkin(1) - pkin(2);
t126 = t63 * qJ(2) + t62 * t70;
t41 = -pkin(7) + t126;
t117 = t63 * qJD(2);
t69 = cos(qJ(4));
t91 = t69 * t117;
t143 = t41 * t116 - t91;
t87 = -t62 * qJ(2) + t63 * t70;
t40 = pkin(3) - t87;
t83 = t69 * pkin(4) + t66 * pkin(8);
t26 = t40 + t83;
t65 = sin(qJ(5));
t23 = t65 * t26;
t68 = cos(qJ(5));
t127 = t68 * t69;
t30 = t41 * t127;
t142 = -t30 - t23;
t121 = qJD(5) * t66;
t115 = t69 * qJD(4);
t92 = t68 * t115;
t141 = -t65 * t121 + t92;
t140 = qJD(5) + qJD(6);
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t43 = t64 * t68 + t67 * t65;
t19 = t140 * t43;
t128 = t67 * t68;
t42 = t64 * t65 - t128;
t10 = t42 * t115 + t19 * t66;
t60 = t68 ^ 2;
t125 = t65 ^ 2 - t60;
t86 = qJD(5) * t125;
t59 = t66 ^ 2;
t85 = (-t69 ^ 2 + t59) * qJD(4);
t120 = qJD(5) * t68;
t118 = t62 * qJD(2);
t82 = -pkin(4) * t66 + pkin(8) * t69;
t33 = t82 * qJD(4) + t118;
t119 = qJD(5) * t69;
t102 = t65 * t119;
t93 = t68 * t116;
t71 = t93 + t102;
t8 = -t26 * t120 - t65 * t33 + t71 * t41 - t68 * t91;
t139 = 0.2e1 * qJD(2);
t138 = pkin(8) + pkin(9);
t98 = t65 * t115;
t38 = t66 * t120 + t98;
t5 = t38 * pkin(9) - t8;
t137 = t67 * t5;
t135 = t41 * t65;
t134 = t62 * t66;
t133 = t62 * t69;
t131 = t65 * t66;
t16 = pkin(9) * t131 - t142;
t132 = t64 * t16;
t130 = t66 * t68;
t129 = t67 * t16;
t123 = pkin(5) * qJD(6);
t122 = qJD(5) * t65;
t114 = -0.2e1 * pkin(4) * qJD(5);
t113 = t64 * t131;
t112 = t65 * t133;
t111 = t62 * t127;
t109 = pkin(5) * t122;
t108 = pkin(5) * t116;
t107 = pkin(8) * t120;
t106 = t64 * t123;
t105 = t67 * t123;
t101 = t68 * t119;
t100 = t59 * t120;
t99 = t68 * t117;
t97 = t65 * t120;
t96 = t66 * t115;
t94 = t62 * t116;
t73 = t143 * t65 + t68 * t33;
t4 = (-pkin(5) * t66 + pkin(9) * t127) * qJD(4) + (-t30 + (-pkin(9) * t66 - t26) * t65) * qJD(5) + t73;
t90 = t67 * t4 - t64 * t5;
t24 = t68 * t26;
t14 = pkin(9) * t130 + t24 + (pkin(5) - t135) * t69;
t89 = -t69 * pkin(5) - t14;
t88 = qJD(5) * t138;
t84 = t65 * t92;
t81 = t67 * t14 - t132;
t80 = t64 * t14 + t129;
t35 = -t65 * t63 + t111;
t75 = t68 * t63 + t112;
t79 = -t64 * t35 - t67 * t75;
t78 = t67 * t35 - t64 * t75;
t49 = t138 * t65;
t50 = t138 * t68;
t77 = -t67 * t49 - t64 * t50;
t76 = -t64 * t49 + t67 * t50;
t74 = t42 * t116 - t69 * t19;
t72 = t41 * t115 + t66 * t117;
t55 = -t68 * pkin(5) - pkin(4);
t52 = -0.2e1 * t96;
t45 = t68 * t88;
t44 = t65 * t88;
t39 = -t65 * t116 + t101;
t32 = t66 * t128 - t113;
t31 = t43 * t66;
t25 = (-pkin(5) * t65 + t41) * t66;
t21 = -t35 * qJD(5) + t65 * t94;
t20 = t75 * qJD(5) + t62 * t93;
t18 = t140 * t42;
t17 = -t38 * pkin(5) + t72;
t15 = -t43 * t116 - t18 * t69;
t13 = -t76 * qJD(6) + t64 * t44 - t67 * t45;
t12 = -t77 * qJD(6) + t67 * t44 + t64 * t45;
t11 = -qJD(6) * t113 + (t140 * t130 + t98) * t67 + t141 * t64;
t9 = t142 * qJD(5) + t73;
t7 = -t78 * qJD(6) + t64 * t20 + t67 * t21;
t6 = -t79 * qJD(6) + t67 * t20 - t64 * t21;
t2 = -t80 * qJD(6) + t90;
t1 = -t81 * qJD(6) - t64 * t4 - t137;
t3 = [0, 0, 0, 0, t139, qJ(2) * t139, 0.2e1 * t118, 0.2e1 * t117 (t126 * t63 - t87 * t62) * t139, 0.2e1 * t96, -0.2e1 * t85, 0, 0, 0, -0.2e1 * t40 * t116 + 0.2e1 * t69 * t118, -0.2e1 * t40 * t115 - 0.2e1 * t66 * t118, -0.2e1 * t59 * t97 + 0.2e1 * t60 * t96, 0.2e1 * t59 * t86 - 0.4e1 * t66 * t84, 0.2e1 * t66 * t102 + 0.2e1 * t68 * t85, 0.2e1 * t66 * t101 - 0.2e1 * t65 * t85, t52, -0.2e1 * t59 * t65 * t117 - 0.2e1 * t24 * t116 + 0.2e1 * t9 * t69 + 0.2e1 * (-t65 * t96 - t100) * t41, 0.2e1 * t8 * t69 + 0.2e1 * (t41 * t122 - t99) * t59 + 0.2e1 * (-t30 + t23) * t116, -0.2e1 * t32 * t10, 0.2e1 * t10 * t31 - 0.2e1 * t32 * t11, 0.2e1 * t10 * t69 + 0.2e1 * t32 * t116, 0.2e1 * t11 * t69 - 0.2e1 * t31 * t116, t52, -0.2e1 * t25 * t11 - 0.2e1 * t81 * t116 - 0.2e1 * t17 * t31 + 0.2e1 * t2 * t69, 0.2e1 * t1 * t69 + 0.2e1 * t25 * t10 + 0.2e1 * t80 * t116 - 0.2e1 * t17 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t116, t63 * t115, 0, 0, 0, 0, 0, -t62 * t100 + t21 * t69 + (t75 - 0.2e1 * t112) * t116, t59 * t62 * t122 + t20 * t69 + (t35 - 0.2e1 * t111) * t116, 0, 0, 0, 0, 0, -t11 * t134 + t7 * t69 + (-t31 * t133 - t79 * t66) * qJD(4), t10 * t134 + t6 * t69 + (-t32 * t133 + t66 * t78) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, t116, 0, -t72, t143, t66 * t86 - t84, t125 * t115 + 0.4e1 * t66 * t97, t39, -t71, 0 (-t107 + (pkin(4) * t65 - t41 * t68) * qJD(4)) * t69 + (pkin(8) * qJD(4) * t65 - t99 + (pkin(4) * t68 + t135) * qJD(5)) * t66 (t83 * qJD(4) + t41 * t121) * t68 + (t82 * qJD(5) + t72) * t65, t10 * t43 + t32 * t18, -t10 * t42 + t43 * t11 - t18 * t31 + t32 * t19, t15, t74, 0, -t31 * t109 - t55 * t11 - t77 * t116 + t13 * t69 + t17 * t42 + t25 * t19, t55 * t10 - t109 * t32 + t76 * t116 + t12 * t69 + t17 * t43 - t25 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t115, t94, 0, 0, 0, 0, 0, -t141 * t62, t38 * t62, 0, 0, 0, 0, 0, t62 * t10 (t43 * t115 - t18 * t66) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t115, 0, 0, 0, 0, 0, -t71, -t39, 0, 0, 0, 0, 0, t74, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t97, -0.2e1 * t86, 0, 0, 0, t65 * t114, t68 * t114, -0.2e1 * t43 * t18, 0.2e1 * t18 * t42 - 0.2e1 * t43 * t19, 0, 0, 0, 0.2e1 * t42 * t109 + 0.2e1 * t55 * t19, 0.2e1 * t109 * t43 - 0.2e1 * t55 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t38, -t116, t9, t8, 0, 0, t10, t11, -t116, -t67 * t108 + (t89 * t64 - t129) * qJD(6) + t90, -t137 + (-t4 + t108) * t64 + (t67 * t89 + t132) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20, 0, 0, 0, 0, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t141, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t122, 0, -t107, pkin(8) * t122, 0, 0, -t18, -t19, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t106, -0.2e1 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, -t116, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
