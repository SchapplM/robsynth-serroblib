% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:03
% EndTime: 2019-03-09 03:03:09
% DurationCPUTime: 2.00s
% Computational Cost: add. (2408->199), mult. (5082->338), div. (0->0), fcn. (4653->8), ass. (0->123)
t127 = cos(pkin(10));
t76 = sin(pkin(10));
t78 = sin(qJ(3));
t80 = cos(qJ(3));
t102 = sin(pkin(9)) * pkin(1) + pkin(7);
t93 = qJ(4) + t102;
t90 = qJD(3) * t93;
t86 = t80 * qJD(4) - t78 * t90;
t87 = -t78 * qJD(4) - t80 * t90;
t82 = t127 * t86 + t76 * t87;
t108 = t127 * t80;
t60 = t76 * t78 - t108;
t61 = t127 * t78 + t76 * t80;
t71 = -cos(pkin(9)) * pkin(1) - pkin(2);
t64 = -pkin(3) * t80 + t71;
t88 = pkin(4) * t60 - pkin(8) * t61 + t64;
t150 = -qJD(5) * t88 - t82;
t124 = t78 * qJD(3);
t55 = qJD(3) * t108 - t124 * t76;
t77 = sin(qJ(5));
t135 = t77 * t55;
t79 = cos(qJ(5));
t73 = qJD(5) * t79;
t33 = t61 * t73 + t135;
t149 = t33 * pkin(5);
t54 = t61 * qJD(3);
t148 = -t54 * t61 - t55 * t60;
t74 = t77 ^ 2;
t75 = t79 ^ 2;
t130 = t74 - t75;
t118 = pkin(3) * t124;
t89 = pkin(4) * t54 - pkin(8) * t55 + t118;
t119 = t150 * t79 - t77 * t89;
t125 = qJD(5) * t77;
t58 = t93 * t80;
t91 = t93 * t78;
t35 = t127 * t58 - t76 * t91;
t4 = t125 * t35 + t119;
t5 = t150 * t77 - t35 * t73 + t79 * t89;
t12 = -t35 * t77 + t79 * t88;
t13 = t79 * t35 + t77 * t88;
t99 = t12 * t77 - t13 * t79;
t147 = qJD(5) * t99 + t4 * t77 - t5 * t79;
t143 = t54 * pkin(5);
t126 = qJD(5) * t61;
t111 = qJ(6) * t126;
t123 = t79 * qJD(6);
t129 = t79 * qJ(6);
t81 = t111 * t77 - t123 * t61 - t129 * t55 + t5;
t1 = t81 + t143;
t6 = pkin(5) * t60 - t129 * t61 + t12;
t138 = t61 * t77;
t7 = -qJ(6) * t138 + t13;
t103 = t6 * t77 - t7 * t79;
t2 = t79 * t111 + (qJ(6) * t55 + qJD(5) * t35 + qJD(6) * t61) * t77 + t119;
t146 = qJD(5) * t103 - t1 * t79 + t2 * t77;
t145 = 0.2e1 * qJD(3);
t144 = 0.2e1 * qJD(5);
t142 = t79 * pkin(5);
t21 = -t127 * t87 + t76 * t86;
t34 = t127 * t91 + t58 * t76;
t141 = t34 * t21;
t140 = t60 * t54;
t139 = t61 * t55;
t137 = t61 * t79;
t47 = t74 * t55;
t48 = t75 * t55;
t136 = t77 * t54;
t134 = t79 * t54;
t133 = t79 * t55;
t132 = t133 * t60 + t134 * t61;
t69 = pkin(3) * t76 + pkin(8);
t128 = qJ(6) + t69;
t122 = t80 * qJD(3);
t40 = 0.2e1 * t140;
t121 = t71 * t145;
t70 = -pkin(3) * t127 - pkin(4);
t120 = t70 * t144;
t117 = pkin(5) * t125;
t116 = t60 * t125;
t115 = t61 * t125;
t113 = t77 * t73;
t112 = t78 * t122;
t63 = t70 - t142;
t110 = -t63 + t142;
t109 = -0.4e1 * t77 * t137;
t59 = t61 ^ 2;
t107 = t59 * t113;
t104 = -t6 * t79 - t7 * t77;
t101 = pkin(5) * t74 + t63 * t79;
t100 = -t12 * t79 - t13 * t77;
t98 = t21 * t60 + t34 * t54;
t96 = -t54 * t69 + t55 * t70;
t56 = t128 * t77;
t57 = t128 * t79;
t95 = t56 * t77 + t57 * t79;
t94 = t60 * t69 - t61 * t70;
t31 = t60 * t73 + t136;
t30 = t115 - t133;
t92 = t102 * qJD(3);
t85 = qJD(5) * t100 - t4 * t79 - t5 * t77;
t44 = t125 * t128 - t123;
t45 = -t77 * qJD(6) - t128 * t73;
t84 = -t44 * t79 - t45 * t77 + (t56 * t79 - t57 * t77) * qJD(5);
t66 = -0.2e1 * t113;
t65 = 0.2e1 * t113;
t62 = -0.2e1 * t130 * qJD(5);
t37 = t61 * t48;
t36 = t61 * t47;
t29 = t116 - t134;
t28 = -t47 - t48;
t20 = pkin(5) * t138 + t34;
t18 = 0.2e1 * t37 - 0.2e1 * t107;
t17 = 0.2e1 * t36 + 0.2e1 * t107;
t16 = t126 * t130 - t77 * t133;
t15 = qJD(5) * t109 - t47 + t48;
t14 = t130 * t144 * t59 + t109 * t55;
t11 = t21 + t149;
t10 = 0.2e1 * t36 + 0.2e1 * t37 + 0.2e1 * t140;
t9 = -0.2e1 * t136 * t61 - 0.2e1 * t33 * t60;
t8 = -0.2e1 * t115 * t60 + 0.2e1 * t132;
t3 = t148 * t79 + t132;
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t112 (-t78 ^ 2 + t80 ^ 2) * t145, 0, -0.2e1 * t112, 0, 0, t78 * t121, t80 * t121, 0, 0, 0.2e1 * t139, 0.2e1 * t148, 0, t40, 0, 0, 0.2e1 * t118 * t60 + 0.2e1 * t54 * t64, 0.2e1 * t118 * t61 + 0.2e1 * t55 * t64, 0.2e1 * t21 * t61 + 0.2e1 * t34 * t55 - 0.2e1 * t35 * t54 - 0.2e1 * t60 * t82, 0.2e1 * t118 * t64 + 0.2e1 * t35 * t82 + 0.2e1 * t141, t18, t14, t8, t17, t9, t40, 0.2e1 * t12 * t54 + 0.2e1 * t138 * t21 + 0.2e1 * t33 * t34 + 0.2e1 * t5 * t60, -0.2e1 * t13 * t54 + 0.2e1 * t137 * t21 - 0.2e1 * t30 * t34 + 0.2e1 * t4 * t60, 0.2e1 * t100 * t55 + 0.2e1 * t147 * t61, 0.2e1 * t12 * t5 - 0.2e1 * t13 * t4 + 0.2e1 * t141, t18, t14, t8, t17, t9, t40, 0.2e1 * t1 * t60 + 0.2e1 * t11 * t138 + 0.2e1 * t20 * t33 + 0.2e1 * t6 * t54, 0.2e1 * t11 * t137 + 0.2e1 * t2 * t60 - 0.2e1 * t20 * t30 - 0.2e1 * t7 * t54, 0.2e1 * t104 * t55 + 0.2e1 * t146 * t61, 0.2e1 * t1 * t6 + 0.2e1 * t11 * t20 - 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t55 + t61 * t82 + t98, 0, 0, 0, 0, 0, 0, 0, t3, 0, -t55 * t99 + t61 * t85 + t98, 0, 0, 0, 0, 0, 0, 0, t3, 0, t11 * t60 + t20 * t54 - t103 * t55 + (qJD(5) * t104 - t1 * t77 - t2 * t79) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t139 + 0.2e1 * t140, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, -t124, 0, -t80 * t92, t78 * t92, 0, 0, 0, 0, t55, 0, -t54, 0, -t21, -t82 (-t127 * t55 - t54 * t76) * pkin(3) (-t127 * t21 + t76 * t82) * pkin(3), -t16, t15, t31, t16, -t29, 0, -t21 * t79 + t96 * t77 + (t34 * t77 - t79 * t94) * qJD(5), t21 * t77 + t96 * t79 + (t34 * t79 + t77 * t94) * qJD(5), t85, t21 * t70 + t69 * t85, -t16, t15, t31, t16, -t29, 0, t63 * t135 - t11 * t79 + t45 * t60 - t56 * t54 + (t101 * t61 + t20 * t77) * qJD(5), t63 * t133 + t11 * t77 + t44 * t60 - t57 * t54 + (t110 * t138 + t20 * t79) * qJD(5) (-t45 * t61 + t55 * t56 - t2 + (-t57 * t61 - t6) * qJD(5)) * t79 + (t44 * t61 - t55 * t57 - t1 + (-t56 * t61 - t7) * qJD(5)) * t77, -t1 * t56 + t11 * t63 + t117 * t20 - t2 * t57 - t44 * t7 + t45 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t122, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55, 0 (-t127 * t54 + t55 * t76) * pkin(3), 0, 0, 0, 0, 0, 0, t29, t31, -t28, t54 * t70 + (t74 + t75) * t69 * t55, 0, 0, 0, 0, 0, 0, t29, t31, -t28, pkin(5) * t116 + t54 * t63 + t55 * t95 + t61 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t62, 0, t66, 0, 0, t77 * t120, t79 * t120, 0, 0, t65, t62, 0, t66, 0, 0, -0.2e1 * t110 * t125, t101 * t144, 0.2e1 * t84, 0.2e1 * t117 * t63 - 0.2e1 * t44 * t57 - 0.2e1 * t45 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t55, 0, t118, 0, 0, 0, 0, 0, 0, -t29, -t31, t28, -t147, 0, 0, 0, 0, 0, 0, -t29, -t31, t28, -t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t95 - t44 * t77 + t45 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, -t33, t54, t5, t4, 0, 0, 0, 0, -t30, 0, -t33, t54, t81 + 0.2e1 * t143, t2, t30 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t30, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t30, 0, -t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, -t125, 0, -t69 * t73, t69 * t125, 0, 0, 0, 0, t73, 0, -t125, 0, t45, t44, -pkin(5) * t73, t45 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t73, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t73, 0, -t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t30, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t73, 0, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t19;