% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:18
% EndTime: 2019-03-09 08:08:21
% DurationCPUTime: 1.13s
% Computational Cost: add. (1555->159), mult. (3663->307), div. (0->0), fcn. (3383->8), ass. (0->104)
t77 = sin(pkin(10));
t75 = t77 ^ 2;
t79 = cos(pkin(10));
t117 = t79 ^ 2 + t75;
t130 = t117 * qJD(4);
t82 = cos(qJ(6));
t112 = qJD(6) * t82;
t80 = sin(qJ(6));
t113 = qJD(6) * t80;
t53 = t77 * t112 - t79 * t113;
t61 = t80 * t77 + t82 * t79;
t63 = t82 * t77 - t80 * t79;
t129 = pkin(4) + pkin(5);
t115 = cos(pkin(9));
t78 = sin(pkin(9));
t81 = sin(qJ(2));
t83 = cos(qJ(2));
t62 = t115 * t81 + t78 * t83;
t128 = pkin(8) * t62;
t71 = t78 * pkin(2) + qJ(4);
t127 = -pkin(8) + t71;
t119 = -qJ(3) - pkin(7);
t102 = qJD(2) * t119;
t49 = t83 * qJD(3) + t81 * t102;
t84 = -t81 * qJD(3) + t83 * t102;
t28 = -t115 * t84 + t78 * t49;
t67 = t119 * t81;
t68 = t119 * t83;
t37 = -t115 * t67 - t78 * t68;
t126 = t37 * t28;
t50 = t62 * qJD(2);
t125 = t77 * t50;
t103 = t115 * t83;
t109 = t81 * qJD(2);
t51 = qJD(2) * t103 - t78 * t109;
t45 = t77 * t51;
t124 = t79 * t50;
t74 = pkin(2) * t109;
t23 = t50 * pkin(3) - t51 * qJ(4) - t62 * qJD(4) + t74;
t29 = t115 * t49 + t78 * t84;
t10 = t77 * t23 + t79 * t29;
t104 = -t83 * pkin(2) - pkin(1);
t60 = t78 * t81 - t103;
t34 = t60 * pkin(3) - t62 * qJ(4) + t104;
t38 = -t115 * t68 + t78 * t67;
t19 = t77 * t34 + t79 * t38;
t118 = t71 * t130;
t116 = qJ(5) * t79;
t111 = t75 * qJD(5);
t110 = t77 * qJD(5);
t108 = t83 * qJD(2);
t107 = -0.2e1 * pkin(1) * qJD(2);
t15 = t60 * qJ(5) + t19;
t26 = t77 * t29;
t9 = t79 * t23 - t26;
t35 = t77 * t38;
t18 = t79 * t34 - t35;
t5 = t50 * qJ(5) + t60 * qJD(5) + t10;
t6 = -t50 * pkin(4) - t9;
t101 = t5 * t77 - t6 * t79;
t73 = -t115 * pkin(2) - pkin(3);
t100 = t10 * t79 - t9 * t77;
t99 = t10 * t77 + t9 * t79;
t11 = t77 * t128 + t15;
t8 = t35 + (-t34 - t128) * t79 - t129 * t60;
t98 = t82 * t11 + t80 * t8;
t97 = t80 * t11 - t82 * t8;
t96 = pkin(4) * t77 - t116;
t88 = -t79 * t62 * qJD(5) + t28;
t12 = t96 * t51 + t88;
t20 = t96 * t62 + t37;
t95 = t12 * t62 + t20 * t51;
t94 = t28 * t62 + t37 * t51;
t52 = t61 * qJD(6);
t93 = t63 * t50 - t52 * t60;
t55 = t127 * t77;
t56 = t127 * t79;
t92 = t82 * t55 - t80 * t56;
t91 = t80 * t55 + t82 * t56;
t90 = -qJD(4) * t60 - t50 * t71;
t89 = -t129 * t77 + t116;
t87 = -t77 * qJ(5) + t73;
t54 = -t79 * pkin(4) + t87;
t86 = t51 * t54 + t90;
t85 = t51 * t73 + t90;
t59 = 0.2e1 * t130;
t46 = t79 * t51;
t44 = t129 * t79 - t87;
t32 = t61 * t62;
t31 = t63 * t62;
t30 = t117 * t51;
t25 = t63 * qJD(4) - t91 * qJD(6);
t24 = -t61 * qJD(4) - t92 * qJD(6);
t22 = t61 * t50 + t53 * t60;
t17 = t89 * t62 - t37;
t16 = -t60 * pkin(4) - t18;
t14 = -t63 * t51 + t62 * t52;
t13 = -t61 * t51 - t53 * t62;
t7 = t89 * t51 - t88;
t4 = pkin(8) * t45 + t5;
t3 = t26 + (-pkin(8) * t51 - t23) * t79 - t129 * t50;
t2 = -t98 * qJD(6) + t82 * t3 - t80 * t4;
t1 = t97 * qJD(6) - t80 * t3 - t82 * t4;
t21 = [0, 0, 0, 0.2e1 * t81 * t108, 0.2e1 * (-t81 ^ 2 + t83 ^ 2) * qJD(2), 0, 0, 0, t81 * t107, t83 * t107, -0.2e1 * t29 * t60 - 0.2e1 * t38 * t50 + 0.2e1 * t94, 0.2e1 * t104 * t74 + 0.2e1 * t38 * t29 + 0.2e1 * t126, 0.2e1 * t18 * t50 + 0.2e1 * t9 * t60 + 0.2e1 * t94 * t77, -0.2e1 * t10 * t60 - 0.2e1 * t19 * t50 + 0.2e1 * t94 * t79, -0.2e1 * t99 * t62 + 0.2e1 * (-t18 * t79 - t19 * t77) * t51, 0.2e1 * t19 * t10 + 0.2e1 * t18 * t9 + 0.2e1 * t126, -0.2e1 * t16 * t50 - 0.2e1 * t6 * t60 + 0.2e1 * t95 * t77, -0.2e1 * t101 * t62 + 0.2e1 * (-t15 * t77 + t16 * t79) * t51, 0.2e1 * t15 * t50 + 0.2e1 * t5 * t60 - 0.2e1 * t95 * t79, 0.2e1 * t20 * t12 + 0.2e1 * t15 * t5 + 0.2e1 * t16 * t6, -0.2e1 * t32 * t13, -0.2e1 * t13 * t31 - 0.2e1 * t32 * t14, 0.2e1 * t13 * t60 - 0.2e1 * t32 * t50, 0.2e1 * t14 * t60 - 0.2e1 * t31 * t50, 0.2e1 * t60 * t50, 0.2e1 * t17 * t14 - 0.2e1 * t2 * t60 - 0.2e1 * t7 * t31 + 0.2e1 * t97 * t50, -0.2e1 * t1 * t60 - 0.2e1 * t17 * t13 + 0.2e1 * t7 * t32 + 0.2e1 * t98 * t50; 0, 0, 0, 0, 0, t108, -t109, 0, -pkin(7) * t108, pkin(7) * t109 (-t115 * t51 - t50 * t78) * pkin(2) (-t115 * t28 + t29 * t78) * pkin(2), -t28 * t79 + t85 * t77, t28 * t77 + t85 * t79, t100, t28 * t73 + t100 * t71 + (-t18 * t77 + t19 * t79) * qJD(4), -t62 * t111 - t12 * t79 + t86 * t77, t5 * t79 + t6 * t77, -t12 * t77 + (t62 * t110 - t86) * t79, t12 * t54 + (qJD(4) * t15 + t5 * t71) * t79 + (qJD(4) * t16 - qJD(5) * t20 + t6 * t71) * t77, -t13 * t63 - t32 * t52, t13 * t61 - t63 * t14 - t52 * t31 - t32 * t53, -t93, t22, 0, -t31 * t110 + t44 * t14 + t17 * t53 - t25 * t60 - t92 * t50 + t7 * t61, t32 * t110 - t44 * t13 - t17 * t52 - t24 * t60 + t91 * t50 + t7 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0.2e1 * t118, 0.2e1 * t79 * t110, t59, 0.2e1 * t111, -0.2e1 * t54 * t110 + 0.2e1 * t118, -0.2e1 * t63 * t52, 0.2e1 * t52 * t61 - 0.2e1 * t63 * t53, 0, 0, 0, 0.2e1 * t61 * t110 + 0.2e1 * t44 * t53, 0.2e1 * t63 * t110 - 0.2e1 * t44 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t124, -t125, -t30, t99, t124, -t30, t125, t101, 0, 0, 0, 0, 0, t22, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t46, 0, t28, t45, 0, -t46, t12, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, 0, 0, 0, 0, -t53, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t46, 0, t6, 0, 0, 0, 0, 0, t60 * t113 - t82 * t50, t60 * t112 + t80 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * qJD(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -t50, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53, 0, t25, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t21;
