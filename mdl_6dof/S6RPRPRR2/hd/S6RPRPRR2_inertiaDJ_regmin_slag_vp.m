% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:38:57
% EndTime: 2019-03-09 03:39:01
% DurationCPUTime: 1.11s
% Computational Cost: add. (1717->170), mult. (3867->299), div. (0->0), fcn. (3775->10), ass. (0->114)
t121 = cos(pkin(11));
t71 = sin(pkin(10)) * pkin(1) + pkin(7);
t122 = qJ(4) + t71;
t78 = sin(pkin(11));
t82 = sin(qJ(3));
t133 = t78 * t82;
t85 = cos(qJ(3));
t57 = t122 * t85;
t40 = t121 * t57 - t122 * t133;
t84 = cos(qJ(5));
t31 = t84 * t40;
t99 = t121 * t85;
t61 = -t99 + t133;
t100 = t121 * t82;
t62 = t78 * t85 + t100;
t73 = -cos(pkin(10)) * pkin(1) - pkin(2);
t88 = -t85 * pkin(3) + t73;
t32 = t61 * pkin(4) - t62 * pkin(8) + t88;
t81 = sin(qJ(5));
t125 = t81 * t32 + t31;
t120 = qJD(5) * t81;
t108 = t62 * t120;
t116 = t82 * qJD(3);
t56 = qJD(3) * t99 - t78 * t116;
t126 = t84 * t56;
t143 = t126 - t108;
t114 = qJD(5) + qJD(6);
t142 = 0.2e1 * qJD(5);
t55 = t62 * qJD(3);
t141 = t55 * pkin(5);
t140 = t61 * pkin(5);
t119 = qJD(5) * t84;
t97 = qJD(3) * t122;
t47 = t85 * qJD(4) - t82 * t97;
t86 = -t82 * qJD(4) - t85 * t97;
t24 = t121 * t47 + t78 * t86;
t75 = pkin(3) * t116;
t30 = t55 * pkin(4) - t56 * pkin(8) + t75;
t6 = -t32 * t119 + t40 * t120 - t84 * t24 - t81 * t30;
t129 = t81 * t56;
t87 = t62 * t119 + t129;
t5 = -t87 * pkin(9) - t6;
t83 = cos(qJ(6));
t139 = t83 * t5;
t70 = t78 * pkin(3) + pkin(8);
t138 = pkin(9) + t70;
t137 = t61 * t55;
t136 = t62 * t56;
t135 = t62 * t81;
t134 = t62 * t84;
t13 = -pkin(9) * t135 + t125;
t80 = sin(qJ(6));
t132 = t80 * t13;
t131 = t80 * t81;
t130 = t81 * t55;
t128 = t83 * t13;
t127 = t84 * t55;
t124 = t61 * t126 + t62 * t127;
t77 = t84 ^ 2;
t123 = t81 ^ 2 - t77;
t118 = qJD(6) * t80;
t117 = qJD(6) * t83;
t115 = t85 * qJD(3);
t72 = -t121 * pkin(3) - pkin(4);
t113 = t72 * t142;
t112 = 0.2e1 * t115;
t111 = pkin(5) * t120;
t110 = pkin(5) * t118;
t109 = pkin(5) * t117;
t107 = t81 * t119;
t103 = -t81 * t24 + t84 * t30;
t4 = -pkin(9) * t126 + t141 + (-t31 + (pkin(9) * t62 - t32) * t81) * qJD(5) + t103;
t106 = t83 * t4 - t80 * t5;
t102 = t84 * t32 - t81 * t40;
t12 = -pkin(9) * t134 + t102 + t140;
t105 = -t12 - t140;
t104 = -0.4e1 * t81 * t134;
t23 = -t121 * t86 + t78 * t47;
t39 = t122 * t100 + t78 * t57;
t101 = qJD(5) * t138;
t98 = t123 * qJD(5);
t96 = t83 * t12 - t132;
t95 = t80 * t12 + t128;
t44 = t114 * t131 - t84 * t117 - t83 * t119;
t64 = t80 * t84 + t83 * t81;
t94 = t44 * t61 - t64 * t55;
t45 = t114 * t64;
t63 = -t83 * t84 + t131;
t93 = t61 * t45 + t55 * t63;
t92 = -t55 * t70 + t56 * t72;
t58 = t138 * t81;
t59 = t138 * t84;
t91 = -t83 * t58 - t80 * t59;
t90 = -t80 * t58 + t83 * t59;
t89 = t61 * t70 - t62 * t72;
t37 = t61 * t119 + t130;
t65 = -t84 * pkin(5) + t72;
t60 = t62 ^ 2;
t54 = t84 * t101;
t53 = t81 * t101;
t43 = 0.2e1 * t137;
t35 = t61 * t120 - t127;
t34 = t63 * t62;
t33 = t64 * t62;
t22 = pkin(5) * t135 + t39;
t16 = -t90 * qJD(6) + t80 * t53 - t83 * t54;
t15 = -t91 * qJD(6) + t83 * t53 + t80 * t54;
t14 = t87 * pkin(5) + t23;
t11 = -t118 * t135 + (t114 * t134 + t129) * t83 + t143 * t80;
t10 = -t83 * t126 + t80 * t129 + t45 * t62;
t7 = -t125 * qJD(5) + t103;
t2 = -t95 * qJD(6) + t106;
t1 = -t96 * qJD(6) - t80 * t4 - t139;
t3 = [0, 0, 0, 0, t82 * t112, 0.2e1 * (-t82 ^ 2 + t85 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t73 * t116, t73 * t112, 0.2e1 * t23 * t62 - 0.2e1 * t24 * t61 + 0.2e1 * t39 * t56 - 0.2e1 * t40 * t55, 0.2e1 * t39 * t23 + 0.2e1 * t40 * t24 + 0.2e1 * t88 * t75, -0.2e1 * t60 * t107 + 0.2e1 * t77 * t136, t123 * t60 * t142 + t56 * t104, -0.2e1 * t61 * t108 + 0.2e1 * t124, -0.2e1 * t62 * t130 - 0.2e1 * t87 * t61, t43, 0.2e1 * t102 * t55 + 0.2e1 * t23 * t135 + 0.2e1 * t87 * t39 + 0.2e1 * t7 * t61, -0.2e1 * t125 * t55 + 0.2e1 * t23 * t134 + 0.2e1 * t143 * t39 + 0.2e1 * t6 * t61, 0.2e1 * t34 * t10, 0.2e1 * t10 * t33 + 0.2e1 * t34 * t11, -0.2e1 * t10 * t61 - 0.2e1 * t34 * t55, -0.2e1 * t61 * t11 - 0.2e1 * t55 * t33, t43, 0.2e1 * t22 * t11 + 0.2e1 * t14 * t33 + 0.2e1 * t2 * t61 + 0.2e1 * t96 * t55, 0.2e1 * t1 * t61 - 0.2e1 * t22 * t10 - 0.2e1 * t14 * t34 - 0.2e1 * t95 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t61 + t24 * t62 + t39 * t55 + t40 * t56, 0, 0, 0, 0, 0, 0 (-t62 * t55 - t56 * t61) * t84 + t124, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t136 + 0.2e1 * t137, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t115, -t116, 0, -t71 * t115, t71 * t116 (-t121 * t56 - t55 * t78) * pkin(3) (-t121 * t23 + t24 * t78) * pkin(3), t81 * t126 - t62 * t98, qJD(5) * t104 - t123 * t56, t37, -t35, 0, -t23 * t84 + t92 * t81 + (t39 * t81 - t89 * t84) * qJD(5), t23 * t81 + t92 * t84 + (t39 * t84 + t89 * t81) * qJD(5), -t10 * t64 + t34 * t44, t10 * t63 - t64 * t11 + t44 * t33 + t34 * t45, -t94, -t93, 0, t65 * t11 + t33 * t111 + t14 * t63 + t16 * t61 + t22 * t45 + t91 * t55, -t65 * t10 - t34 * t111 + t14 * t64 + t15 * t61 - t22 * t44 - t90 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t115, 0 (-t121 * t55 + t56 * t78) * pkin(3), 0, 0, 0, 0, 0, t35, t37, 0, 0, 0, 0, 0, t93, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t107, -0.2e1 * t98, 0, 0, 0, t81 * t113, t84 * t113, -0.2e1 * t64 * t44, 0.2e1 * t44 * t63 - 0.2e1 * t64 * t45, 0, 0, 0, 0.2e1 * t63 * t111 + 0.2e1 * t65 * t45, 0.2e1 * t64 * t111 - 0.2e1 * t65 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, -t35, -t37, 0, 0, 0, 0, 0, -t93, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, -t87, t55, t7, t6, 0, 0, -t10, -t11, t55, t83 * t141 + (t105 * t80 - t128) * qJD(6) + t106, -t139 + (-t4 - t141) * t80 + (t105 * t83 + t132) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t143, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t120, 0, -t70 * t119, t70 * t120, 0, 0, -t44, -t45, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t110, -0.2e1 * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, t55, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
