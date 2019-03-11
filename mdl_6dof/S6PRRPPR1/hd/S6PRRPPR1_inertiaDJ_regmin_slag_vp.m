% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:12
% EndTime: 2019-03-08 21:02:16
% DurationCPUTime: 1.29s
% Computational Cost: add. (1513->165), mult. (3973->332), div. (0->0), fcn. (3918->12), ass. (0->109)
t72 = sin(pkin(12));
t75 = cos(pkin(12));
t77 = sin(qJ(6));
t80 = cos(qJ(6));
t131 = -t77 * t72 + t80 * t75;
t58 = t80 * t72 + t77 * t75;
t51 = t58 * qJD(6);
t130 = t131 * qJD(6);
t129 = 0.2e1 * t130;
t73 = sin(pkin(11));
t64 = t73 * pkin(3) + qJ(5);
t128 = pkin(9) + t64;
t115 = cos(pkin(11));
t82 = cos(qJ(2));
t113 = qJD(2) * t82;
t74 = sin(pkin(6));
t105 = t74 * t113;
t81 = cos(qJ(3));
t79 = sin(qJ(2));
t122 = t74 * t79;
t76 = cos(pkin(6));
t78 = sin(qJ(3));
t86 = -t78 * t122 + t76 * t81;
t41 = t86 * qJD(3) + t81 * t105;
t52 = t81 * t122 + t76 * t78;
t83 = -t52 * qJD(3) - t78 * t105;
t17 = -t115 * t83 + t73 * t41;
t34 = -t115 * t86 + t73 * t52;
t127 = t34 * t17;
t117 = -qJ(4) - pkin(8);
t102 = qJD(3) * t117;
t47 = t81 * qJD(4) + t78 * t102;
t84 = -t78 * qJD(4) + t81 * t102;
t30 = -t115 * t84 + t73 * t47;
t104 = t115 * t78;
t61 = t117 * t81;
t42 = -t117 * t104 - t73 * t61;
t126 = t42 * t30;
t57 = t73 * t81 + t104;
t125 = t57 * t72;
t103 = t115 * t81;
t112 = t78 * qJD(3);
t49 = qJD(3) * t103 - t73 * t112;
t124 = t72 * t49;
t123 = t73 * t78;
t121 = t74 * t82;
t120 = t75 * t49;
t48 = t57 * qJD(3);
t69 = pkin(3) * t112;
t23 = t48 * pkin(4) - t49 * qJ(5) - t57 * qJD(5) + t69;
t31 = t115 * t47 + t73 * t84;
t8 = t72 * t23 + t75 * t31;
t55 = -t103 + t123;
t68 = -t81 * pkin(3) - pkin(2);
t38 = t55 * pkin(4) - t57 * qJ(5) + t68;
t43 = -t115 * t61 + t117 * t123;
t16 = t72 * t38 + t75 * t43;
t116 = t72 ^ 2 + t75 ^ 2;
t114 = qJD(2) * t79;
t111 = t81 * qJD(3);
t110 = -0.2e1 * pkin(2) * qJD(3);
t109 = t82 * t112;
t106 = t74 * t114;
t7 = t75 * t23 - t72 * t31;
t15 = t75 * t38 - t72 * t43;
t101 = 0.2e1 * t116 * qJD(5);
t100 = t7 * t75 + t8 * t72;
t99 = -t7 * t72 + t8 * t75;
t67 = -t115 * pkin(3) - pkin(4);
t10 = -pkin(9) * t125 + t16;
t9 = -t75 * t57 * pkin(9) + t55 * pkin(5) + t15;
t98 = t80 * t10 + t77 * t9;
t97 = t77 * t10 - t80 * t9;
t18 = t115 * t41 + t73 * t83;
t13 = t75 * t106 - t72 * t18;
t14 = t72 * t106 + t75 * t18;
t96 = t13 * t75 + t14 * t72;
t95 = -t13 * t72 + t14 * t75;
t94 = t17 * t42 + t34 * t30;
t93 = t17 * t57 + t34 * t49;
t35 = t115 * t52 + t73 * t86;
t26 = -t75 * t121 - t72 * t35;
t27 = -t72 * t121 + t75 * t35;
t92 = t80 * t26 - t77 * t27;
t91 = t77 * t26 + t80 * t27;
t90 = t30 * t57 + t42 * t49;
t89 = t130 * t55 + t58 * t48;
t53 = t128 * t72;
t54 = t128 * t75;
t88 = -t80 * t53 - t77 * t54;
t87 = -t77 * t53 + t80 * t54;
t85 = -qJD(5) * t55 - t48 * t64 + t49 * t67;
t60 = -t75 * pkin(5) + t67;
t33 = t131 * t57;
t32 = t58 * t57;
t29 = pkin(5) * t125 + t42;
t25 = -t58 * qJD(5) - t87 * qJD(6);
t24 = -qJD(5) * t131 - t88 * qJD(6);
t22 = t131 * t48 - t51 * t55;
t19 = pkin(5) * t124 + t30;
t12 = t130 * t57 + t58 * t49;
t11 = t131 * t49 - t57 * t51;
t6 = -pkin(9) * t124 + t8;
t5 = t48 * pkin(5) - pkin(9) * t120 + t7;
t4 = -t91 * qJD(6) + t80 * t13 - t77 * t14;
t3 = -t92 * qJD(6) - t77 * t13 - t80 * t14;
t2 = -t98 * qJD(6) + t80 * t5 - t77 * t6;
t1 = t97 * qJD(6) - t77 * t5 - t80 * t6;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t74 ^ 2 * t79 * t113 + 0.2e1 * t35 * t18 + 0.2e1 * t127, 0, 0, 0, 0.2e1 * t26 * t13 + 0.2e1 * t27 * t14 + 0.2e1 * t127, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t106, -t105, 0, 0, 0, 0, 0 (-t81 * t114 - t109) * t74 (-t82 * t111 + t78 * t114) * t74, -t18 * t55 - t35 * t48 + t93, t18 * t43 + t35 * t31 + (-pkin(3) * t109 + t68 * t114) * t74 + t94, t13 * t55 + t26 * t48 + t93 * t72, -t14 * t55 - t27 * t48 + t93 * t75, -t96 * t57 + (-t26 * t75 - t27 * t72) * t49, t13 * t15 + t14 * t16 + t26 * t7 + t27 * t8 + t94, 0, 0, 0, 0, 0, t34 * t12 + t17 * t32 + t4 * t55 + t92 * t48, t34 * t11 + t17 * t33 + t3 * t55 - t91 * t48; 0, 0, 0, 0, 0.2e1 * t78 * t111, 0.2e1 * (-t78 ^ 2 + t81 ^ 2) * qJD(3), 0, 0, 0, t78 * t110, t81 * t110, -0.2e1 * t31 * t55 - 0.2e1 * t43 * t48 + 0.2e1 * t90, 0.2e1 * t43 * t31 + 0.2e1 * t68 * t69 + 0.2e1 * t126, 0.2e1 * t15 * t48 + 0.2e1 * t7 * t55 + 0.2e1 * t90 * t72, -0.2e1 * t16 * t48 - 0.2e1 * t8 * t55 + 0.2e1 * t90 * t75, -0.2e1 * t100 * t57 + 0.2e1 * (-t15 * t75 - t16 * t72) * t49, 0.2e1 * t15 * t7 + 0.2e1 * t16 * t8 + 0.2e1 * t126, 0.2e1 * t33 * t11, -0.2e1 * t11 * t32 - 0.2e1 * t33 * t12, 0.2e1 * t11 * t55 + 0.2e1 * t33 * t48, -0.2e1 * t12 * t55 - 0.2e1 * t32 * t48, 0.2e1 * t55 * t48, 0.2e1 * t29 * t12 + 0.2e1 * t19 * t32 + 0.2e1 * t2 * t55 - 0.2e1 * t97 * t48, 0.2e1 * t1 * t55 + 0.2e1 * t29 * t11 + 0.2e1 * t19 * t33 - 0.2e1 * t98 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t41, 0 (-t115 * t17 + t18 * t73) * pkin(3), -t17 * t75, t17 * t72, t95, t17 * t67 + t95 * t64 + (-t26 * t72 + t27 * t75) * qJD(5), 0, 0, 0, 0, 0, -t131 * t17 + t34 * t51, t130 * t34 + t17 * t58; 0, 0, 0, 0, 0, 0, t111, -t112, 0, -pkin(8) * t111, pkin(8) * t112 (-t115 * t49 - t48 * t73) * pkin(3) (-t115 * t30 + t31 * t73) * pkin(3), -t30 * t75 + t85 * t72, t30 * t72 + t85 * t75, t99, t30 * t67 + t99 * t64 + (-t15 * t72 + t16 * t75) * qJD(5), t11 * t58 + t130 * t33, t11 * t131 - t58 * t12 - t130 * t32 - t33 * t51, t89, t22, 0, t60 * t12 - t131 * t19 + t25 * t55 + t29 * t51 + t88 * t48, t60 * t11 + t130 * t29 + t19 * t58 + t24 * t55 - t87 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t64 * t101, t58 * t129, 0.2e1 * t130 * t131 - 0.2e1 * t58 * t51, 0, 0, 0, 0.2e1 * t60 * t51, t60 * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, 0, 0, t96, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t75 * t48, -t72 * t48, -t116 * t49, t100, 0, 0, 0, 0, 0, t22, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t120, 0, t30, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12, t48, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t51, 0, t25, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t20;
