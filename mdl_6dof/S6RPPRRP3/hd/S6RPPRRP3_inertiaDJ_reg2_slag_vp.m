% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:48
% EndTime: 2019-03-09 02:03:54
% DurationCPUTime: 1.70s
% Computational Cost: add. (1101->168), mult. (2199->270), div. (0->0), fcn. (1600->6), ass. (0->113)
t57 = sin(qJ(5));
t52 = t57 ^ 2;
t59 = cos(qJ(5));
t54 = t59 ^ 2;
t116 = t52 + t54;
t60 = cos(qJ(4));
t49 = t60 * qJD(4);
t28 = t116 * t49;
t58 = sin(qJ(4));
t93 = t58 * t49;
t43 = 0.2e1 * t93;
t117 = t52 - t54;
t38 = t117 * qJD(5);
t108 = qJ(6) * qJD(4);
t135 = (pkin(5) * qJD(5) - qJD(6)) * t60 + t58 * t108;
t110 = t58 * qJD(4);
t76 = t59 * pkin(5) + t57 * qJ(6);
t132 = t76 * qJD(5) - t59 * qJD(6);
t47 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t75 = pkin(5) * t57 - qJ(6) * t59;
t72 = -t47 + t75;
t5 = -t72 * t110 + t132 * t60;
t36 = t59 * t58 * t47;
t126 = t58 * pkin(4);
t124 = t60 * pkin(8);
t48 = sin(pkin(9)) * pkin(1) + qJ(3);
t73 = t48 - t124;
t69 = t73 + t126;
t119 = t57 * t69 + t36;
t7 = t58 * qJ(6) + t119;
t125 = t59 * pkin(4);
t122 = t57 * t47;
t89 = -pkin(5) + t122;
t8 = -t59 * t73 + (t89 - t125) * t58;
t81 = t57 * t8 + t59 * t7;
t134 = t81 * qJD(4) - t5;
t24 = t75 * qJD(5) - t57 * qJD(6);
t120 = t60 * t24;
t19 = t72 * t60;
t37 = -pkin(4) - t76;
t133 = (t37 * t58 + t124) * qJD(4) - qJD(5) * t19 - t120;
t127 = pkin(8) * t58;
t128 = pkin(4) * t60;
t82 = t127 + t128;
t68 = t82 * qJD(4) + qJD(3);
t61 = -qJD(5) * t119 + t59 * t68;
t131 = (-0.1e1 + t116) * t43;
t130 = 0.2e1 * qJD(3);
t129 = 0.2e1 * qJD(6);
t123 = t37 * t60;
t121 = t59 * t60;
t118 = pkin(8) * t28;
t53 = t58 ^ 2;
t55 = t60 ^ 2;
t115 = t53 - t55;
t114 = t53 + t55;
t113 = qJD(4) * t19;
t112 = qJD(5) * t57;
t50 = qJD(5) * t59;
t111 = qJD(5) * t60;
t107 = -0.2e1 * pkin(4) * qJD(5);
t106 = t58 * t122;
t105 = t48 * t130;
t67 = t59 * t69;
t92 = t47 * t49;
t104 = qJD(5) * t67 + t57 * t68 + t59 * t92;
t102 = pkin(8) * t112;
t101 = pkin(8) * t50;
t100 = t47 * t112;
t99 = t58 * t112;
t98 = t57 * t111;
t97 = t59 * t111;
t96 = t57 * t50;
t95 = t47 * t110;
t94 = t59 * t110;
t90 = t60 * t108;
t88 = t116 * t58;
t87 = qJD(4) * t115;
t85 = t57 * t94;
t84 = t55 * t96;
t83 = pkin(8) * t88;
t80 = t57 * t7 - t59 * t8;
t9 = t67 - t106;
t78 = t119 * t59 - t57 * t9;
t77 = t119 * t57 + t59 * t9;
t71 = pkin(5) * t110 - qJ(6) * t111;
t70 = -t5 + (t123 - t127) * qJD(5);
t1 = t90 + (qJD(6) - t100) * t58 + t104;
t2 = t89 * t49 - t61;
t65 = -t80 * qJD(5) + t1 * t59 + t2 * t57;
t3 = t47 * t99 - t104;
t4 = -t57 * t92 + t61;
t64 = -t77 * qJD(5) - t3 * t59 - t4 * t57;
t62 = t65 + t113;
t51 = qJD(4) * t53;
t42 = -0.2e1 * t96;
t41 = 0.2e1 * t96;
t34 = -t57 * t110 + t97;
t33 = t57 * t49 + t58 * t50;
t32 = t114 * t50;
t31 = t94 + t98;
t30 = -t59 * t49 + t99;
t29 = t114 * t112;
t27 = qJD(4) * t88;
t21 = -0.2e1 * t54 * t93 - 0.2e1 * t84;
t20 = -0.2e1 * t52 * t93 + 0.2e1 * t84;
t18 = t117 * t111 + t85;
t17 = -t57 * t87 + t58 * t97;
t16 = -t117 * t110 + 0.4e1 * t60 * t96;
t15 = -0.2e1 * t58 * t98 - 0.2e1 * t59 * t87;
t14 = t55 * t38 + 0.2e1 * t60 * t85;
t6 = t51 + (-t116 * t115 - t55) * qJD(4);
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t105, -0.2e1 * t93, -0.2e1 * t55 * qJD(4) + 0.2e1 * t51, 0, t43, 0, 0, 0.2e1 * qJD(3) * t58 + 0.2e1 * t48 * t49, 0.2e1 * qJD(3) * t60 - 0.2e1 * t48 * t110, 0, t105, t21, 0.2e1 * t14, t15, t20, -0.2e1 * t17, t43, -0.2e1 * t55 * t47 * t50 + 0.2e1 * t4 * t58 + 0.2e1 * (t9 + 0.2e1 * t106) * t49, 0.2e1 * t55 * t100 + 0.2e1 * t3 * t58 + 0.2e1 * (-t119 + 0.2e1 * t36) * t49, 0.2e1 * t77 * t110 + 0.2e1 * (-t78 * qJD(5) + t3 * t57 - t4 * t59) * t60, -0.2e1 * t47 ^ 2 * t93 - 0.2e1 * t119 * t3 + 0.2e1 * t9 * t4, t21, t15, -0.2e1 * t14, t43, 0.2e1 * t17, t20, 0.2e1 * (-t57 * t113 - t2) * t58 + 0.2e1 * (-qJD(4) * t8 + t19 * t50 + t5 * t57) * t60, 0.2e1 * t80 * t110 + 0.2e1 * (-t81 * qJD(5) - t1 * t57 + t2 * t59) * t60, 0.2e1 * (t59 * t113 + t1) * t58 + 0.2e1 * (qJD(4) * t7 + t19 * t112 - t5 * t59) * t60, 0.2e1 * t7 * t1 + 0.2e1 * t19 * t5 + 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t60 + (t115 * t47 - t78 * t58) * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t58 + t62 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t29, 0, t78 * t49 + (t64 - 0.2e1 * t92) * t58, 0, 0, 0, 0, 0, 0, -t32, 0, -t29, t134 * t60 + t62 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, -t49, 0, -t95, -t92, 0, 0, -t18, -t16, t33, t18, -t30, 0 (-t60 * t122 - t82 * t59) * qJD(5) + (-t36 + (-t124 + t126) * t57) * qJD(4) (-t47 * t121 + t82 * t57) * qJD(5) + (-pkin(8) * t121 + (t122 + t125) * t58) * qJD(4), t64, -pkin(4) * t95 + t64 * pkin(8), -t18, t33, t16, 0, t30, t18, -t133 * t57 + t70 * t59, t65, t133 * t59 + t70 * t57, t65 * pkin(8) + t19 * t24 + t5 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t110, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, -t27 (-t83 - t128) * qJD(4), 0, 0, 0, 0, 0, 0, t30, -t27, -t33, t58 * t24 + (-t83 + t123) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t49, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t34, t28, -pkin(4) * t110 + t118, 0, 0, 0, 0, 0, 0, -t31, t28, t34, t37 * t110 + t118 - t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -0.2e1 * t38, 0, t42, 0, 0, t57 * t107, t59 * t107, 0, 0, t41, 0, 0.2e1 * t38, 0, 0, t42, 0.2e1 * t37 * t112 - 0.2e1 * t24 * t59, 0, -0.2e1 * t24 * t57 - 0.2e1 * t37 * t50, 0.2e1 * t37 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t34, t49, t4, t3, 0, 0, 0, -t31, 0, t49, t34, 0 (0.2e1 * pkin(5) - t122) * t49 + t61, t135 * t57 + t71 * t59, 0.2e1 * t90 + (t129 - t100) * t58 + t104, -t2 * pkin(5) + t1 * qJ(6) + t7 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t31, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, -t31, -t135 * t59 + t71 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t30, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, -t30, -t132 * t58 - t75 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t112, 0, -t101, t102, 0, 0, 0, t50, 0, 0, t112, 0, -t101, -t132, -t102, -t132 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, qJ(6) * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t31, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
