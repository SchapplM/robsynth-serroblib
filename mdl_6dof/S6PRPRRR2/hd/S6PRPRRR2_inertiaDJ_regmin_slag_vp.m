% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:55
% EndTime: 2019-03-08 20:29:59
% DurationCPUTime: 1.20s
% Computational Cost: add. (1048->171), mult. (2959->312), div. (0->0), fcn. (2863->12), ass. (0->122)
t74 = sin(qJ(4));
t137 = -0.4e1 * t74;
t70 = cos(pkin(12));
t60 = -t70 * pkin(2) - pkin(3);
t78 = cos(qJ(4));
t92 = -t78 * pkin(4) - t74 * pkin(9);
t45 = t60 + t92;
t73 = sin(qJ(5));
t36 = t73 * t45;
t77 = cos(qJ(5));
t126 = t77 * t78;
t68 = sin(pkin(12));
t59 = t68 * pkin(2) + pkin(8);
t49 = t59 * t126;
t136 = -t49 - t36;
t114 = t78 * qJD(4);
t100 = t77 * t114;
t118 = qJD(5) * t73;
t135 = -t74 * t118 + t100;
t66 = t77 ^ 2;
t123 = t73 ^ 2 - t66;
t95 = t123 * qJD(5);
t134 = qJD(5) + qJD(6);
t117 = qJD(5) * t77;
t116 = qJD(5) * t78;
t106 = t73 * t116;
t63 = t74 * qJD(4);
t40 = t63 * t77 + t106;
t91 = pkin(4) * t74 - pkin(9) * t78;
t52 = t91 * qJD(4);
t10 = -t45 * t117 + t40 * t59 - t73 * t52;
t133 = pkin(9) + pkin(10);
t76 = cos(qJ(6));
t104 = t73 * t114;
t80 = t117 * t74 + t104;
t9 = -t80 * pkin(10) - t10;
t132 = t76 * t9;
t131 = t59 * t73;
t129 = t73 * t74;
t22 = -pkin(10) * t129 - t136;
t72 = sin(qJ(6));
t130 = t72 * t22;
t128 = t74 * t77;
t127 = t76 * t22;
t101 = t59 * t63;
t124 = t73 * t101 + t77 * t52;
t65 = t74 ^ 2;
t122 = -t78 ^ 2 + t65;
t69 = sin(pkin(6));
t121 = qJD(2) * t69;
t75 = sin(qJ(2));
t79 = cos(qJ(2));
t32 = (t68 * t79 + t70 * t75) * t69;
t71 = cos(pkin(6));
t25 = t32 * t74 - t71 * t78;
t120 = qJD(4) * t25;
t119 = qJD(4) * t77;
t115 = qJD(6) * t72;
t113 = -0.2e1 * pkin(4) * qJD(5);
t112 = 0.2e1 * qJD(4) * t60;
t111 = pkin(5) * t118;
t110 = pkin(5) * t63;
t109 = pkin(5) * t115;
t108 = qJD(6) * t76 * pkin(5);
t105 = t77 * t116;
t103 = t73 * t117;
t102 = t74 * t114;
t99 = t59 * t114;
t8 = (pkin(5) * t74 - pkin(10) * t126) * qJD(4) + (-t49 + (pkin(10) * t74 - t45) * t73) * qJD(5) + t124;
t98 = -t72 * t9 + t76 * t8;
t37 = t77 * t45;
t21 = -pkin(10) * t128 + t37 + (-pkin(5) - t131) * t78;
t97 = t78 * pkin(5) - t21;
t96 = qJD(5) * t133;
t94 = t122 * qJD(4);
t93 = t73 * t100;
t26 = t32 * t78 + t71 * t74;
t84 = t68 * t75 - t70 * t79;
t31 = t84 * t69;
t18 = -t26 * t73 + t31 * t77;
t19 = t26 * t77 + t31 * t73;
t90 = t76 * t18 - t72 * t19;
t89 = t72 * t18 + t76 * t19;
t88 = t76 * t21 - t130;
t87 = t72 * t21 + t127;
t56 = t133 * t73;
t57 = t133 * t77;
t86 = -t76 * t56 - t72 * t57;
t85 = -t72 * t56 + t76 * t57;
t48 = t72 * t77 + t76 * t73;
t47 = t72 * t73 - t76 * t77;
t24 = t134 * t48;
t83 = -t78 * t24 + t47 * t63;
t30 = t84 * t121;
t16 = t26 * qJD(4) - t30 * t74;
t82 = t117 * t25 + t16 * t73;
t81 = t118 * t25 - t16 * t77;
t62 = -t77 * pkin(5) - pkin(4);
t58 = -0.2e1 * t102;
t51 = t77 * t96;
t50 = t73 * t96;
t42 = t63 * t73 - t105;
t38 = (pkin(5) * t73 + t59) * t74;
t35 = t47 * t74;
t34 = t48 * t74;
t29 = qJD(2) * t32;
t27 = t80 * pkin(5) + t99;
t23 = t134 * t47;
t20 = t23 * t78 + t48 * t63;
t17 = -t30 * t78 - t120;
t15 = -t85 * qJD(6) + t72 * t50 - t76 * t51;
t14 = -t86 * qJD(6) + t76 * t50 + t72 * t51;
t13 = -t115 * t129 + (t134 * t128 + t104) * t76 + t135 * t72;
t12 = -t76 * t100 + t72 * t104 + t24 * t74;
t11 = t136 * qJD(5) + t124;
t6 = t18 * qJD(5) + t17 * t77 + t29 * t73;
t5 = -t19 * qJD(5) - t17 * t73 + t29 * t77;
t4 = -t87 * qJD(6) + t98;
t3 = -t88 * qJD(6) - t72 * t8 - t132;
t2 = -t89 * qJD(6) + t76 * t5 - t72 * t6;
t1 = -t90 * qJD(6) - t72 * t5 - t76 * t6;
t7 = [0, 0, 0, 0, 0.2e1 * t31 * t29 - 0.2e1 * t32 * t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t75 * t121, -t79 * t121 (-t29 * t70 - t30 * t68) * pkin(2), 0, 0, 0, 0, 0, -t29 * t78 + t31 * t63, t31 * t114 + t29 * t74, 0, 0, 0, 0, 0 (t73 * t120 - t5) * t78 + (qJD(4) * t18 + t82) * t74 (t25 * t119 + t6) * t78 + (-qJD(4) * t19 - t81) * t74, 0, 0, 0, 0, 0, t25 * t13 + t16 * t34 - t2 * t78 + t63 * t90, -t1 * t78 - t25 * t12 - t16 * t35 - t63 * t89; 0, 0, 0, 0, 0, 0.2e1 * t102, -0.2e1 * t94, 0, 0, 0, t74 * t112, t78 * t112, 0.2e1 * t66 * t102 - 0.2e1 * t65 * t103, t93 * t137 + 0.2e1 * t65 * t95, 0.2e1 * t74 * t106 + 0.2e1 * t122 * t119, 0.2e1 * t74 * t105 - 0.2e1 * t73 * t94, t58, 0.2e1 * t37 * t63 - 0.2e1 * t11 * t78 + 0.2e1 * (t73 * t102 + t117 * t65) * t59, -0.2e1 * t65 * t59 * t118 - 0.2e1 * t10 * t78 + 0.2e1 * (-t36 + t49) * t63, 0.2e1 * t35 * t12, 0.2e1 * t12 * t34 + 0.2e1 * t35 * t13, 0.2e1 * t12 * t78 - 0.2e1 * t35 * t63, 0.2e1 * t78 * t13 - 0.2e1 * t34 * t63, t58, 0.2e1 * t38 * t13 + 0.2e1 * t27 * t34 - 0.2e1 * t4 * t78 + 0.2e1 * t63 * t88, -0.2e1 * t38 * t12 - 0.2e1 * t27 * t35 - 0.2e1 * t3 * t78 - 0.2e1 * t63 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, t81, t82, 0, 0, 0, 0, 0, t16 * t47 + t25 * t24, t16 * t48 - t25 * t23; 0, 0, 0, 0, 0, 0, 0, t114, -t63, 0, -t99, t101, -t74 * t95 + t93, t103 * t137 - t123 * t114, t42, t40, 0 (pkin(9) * t126 + (-pkin(4) * t77 + t131) * t74) * qJD(5) + (t92 * t73 - t49) * qJD(4) (t59 * t128 + t91 * t73) * qJD(5) + (t78 * t131 + t92 * t77) * qJD(4), -t12 * t48 + t35 * t23, t12 * t47 - t48 * t13 + t23 * t34 + t35 * t24, t20, -t83, 0, t34 * t111 + t62 * t13 - t15 * t78 + t38 * t24 + t27 * t47 + t63 * t86, -t35 * t111 - t62 * t12 - t14 * t78 - t38 * t23 + t27 * t48 - t63 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t114, 0, 0, 0, 0, 0, -t40, t42, 0, 0, 0, 0, 0, t83, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t103, -0.2e1 * t95, 0, 0, 0, t73 * t113, t77 * t113, -0.2e1 * t48 * t23, 0.2e1 * t23 * t47 - 0.2e1 * t48 * t24, 0, 0, 0, 0.2e1 * t47 * t111 + 0.2e1 * t62 * t24, 0.2e1 * t48 * t111 - 0.2e1 * t62 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, -t80, t63, t11, t10, 0, 0, -t12, -t13, t63, t76 * t110 + (t97 * t72 - t127) * qJD(6) + t98, -t132 + (-t8 - t110) * t72 + (t97 * t76 + t130) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t135, 0, 0, 0, 0, 0, -t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, -t118, 0, -pkin(9) * t117, pkin(9) * t118, 0, 0, -t23, -t24, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t109, -0.2e1 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, t63, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
