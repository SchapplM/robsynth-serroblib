% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:25
% EndTime: 2022-01-20 11:17:27
% DurationCPUTime: 0.65s
% Computational Cost: add. (770->113), mult. (1948->180), div. (0->0), fcn. (1587->8), ass. (0->108)
t136 = qJD(4) + qJD(5);
t83 = sin(qJ(5));
t84 = sin(qJ(4));
t86 = cos(qJ(5));
t87 = cos(qJ(4));
t90 = t83 * t87 + t84 * t86;
t38 = t136 * t90;
t81 = sin(pkin(9));
t23 = t38 * t81;
t137 = 0.2e1 * t23;
t135 = 2 * qJD(4);
t134 = pkin(8) * t81;
t88 = cos(qJ(2));
t133 = t88 * pkin(1);
t82 = cos(pkin(9));
t132 = t38 * t82;
t117 = pkin(1) * qJD(2);
t105 = t88 * t117;
t69 = qJD(3) + t105;
t79 = t81 ^ 2;
t56 = t79 * t69;
t131 = t81 * t84;
t130 = t82 * t84;
t129 = t82 * t87;
t128 = t86 * t87;
t112 = qJD(5) * t83;
t114 = qJD(4) * t84;
t99 = t81 * t114;
t24 = t136 * t81 * t128 - t112 * t131 - t83 * t99;
t113 = qJD(4) * t87;
t98 = t81 * t113;
t66 = pkin(4) * t98;
t40 = t81 * t69 + t66;
t43 = t90 * t81;
t75 = pkin(4) * t131;
t85 = sin(qJ(2));
t76 = t85 * pkin(1) + qJ(3);
t45 = t81 * t76 + t75;
t127 = t45 * t24 + t40 * t43;
t48 = t81 * qJD(3) + t66;
t53 = t81 * qJ(3) + t75;
t126 = t53 * t24 + t48 * t43;
t115 = qJD(3) * t82;
t59 = -t82 * pkin(3) - t81 * pkin(7) - pkin(2);
t125 = t59 * t113 + t87 * t115;
t100 = t79 * t113;
t124 = t76 * t100 + t84 * t56;
t106 = t85 * t117;
t47 = t59 - t133;
t107 = t84 * t106 + t47 * t113 + t69 * t129;
t73 = t82 * t114;
t15 = t76 * t73 - t107;
t123 = -t15 * t82 + t87 * t56;
t80 = t82 ^ 2;
t122 = t80 * t69 + t56;
t77 = t79 * qJD(3);
t121 = qJ(3) * t100 + t84 * t77;
t97 = qJ(3) * t114;
t35 = t82 * t97 - t125;
t120 = -t35 * t82 + t87 * t77;
t119 = t80 * qJD(3) + t77;
t118 = t79 + t80;
t116 = qJ(3) * t84;
t111 = t87 * t134;
t108 = t76 * t129;
t16 = -t69 * t130 + t87 * t106 + (-t84 * t47 - t108) * qJD(4);
t65 = pkin(8) * t99;
t10 = t65 + t16;
t96 = t47 - t134;
t22 = t96 * t87 + (-t76 * t84 - pkin(4)) * t82;
t29 = t96 * t84 + t108;
t9 = (-t76 * t130 - t111) * qJD(4) + t107;
t2 = -t83 * t10 - t86 * t9 + (-t22 * t86 + t29 * t83) * qJD(5);
t89 = t83 * t84 - t128;
t44 = t89 * t81;
t110 = -t2 * t82 - t45 * t23 - t40 * t44;
t25 = (-t82 * t116 - t111) * qJD(4) + t125;
t102 = qJ(3) * t129;
t36 = -t84 * t115 + (-t84 * t59 - t102) * qJD(4);
t26 = t65 + t36;
t95 = t59 - t134;
t31 = t95 * t87 + (-pkin(4) - t116) * t82;
t39 = t95 * t84 + t102;
t5 = -t86 * t25 - t83 * t26 + (-t31 * t86 + t39 * t83) * qJD(5);
t109 = -t53 * t23 - t48 * t44 - t5 * t82;
t104 = pkin(4) * t112;
t103 = qJD(5) * t86 * pkin(4);
t101 = t79 * t114;
t94 = t118 * t69;
t93 = t118 * qJD(3);
t92 = t81 * t82 * t135;
t91 = t82 * t106;
t3 = -t83 * t9 + t86 * t10 + (-t22 * t83 - t29 * t86) * qJD(5);
t6 = -t83 * t25 + t86 * t26 + (-t31 * t83 - t39 * t86) * qJD(5);
t74 = t82 * t113;
t64 = t82 * t103;
t63 = t82 * t104;
t58 = -0.2e1 * t84 * t100;
t55 = t87 * t92;
t54 = t84 * t92;
t42 = (t84 ^ 2 - t87 ^ 2) * t79 * t135;
t37 = t136 * t89;
t32 = t37 * t82;
t20 = 0.2e1 * t24 * t82;
t19 = t82 * t137;
t11 = t44 * t137;
t7 = 0.2e1 * t23 * t43 + 0.2e1 * t44 * t24;
t1 = [0, 0, 0, 0, -0.2e1 * t106, -0.2e1 * t105, -0.2e1 * t91, 0.2e1 * t122, 0.2e1 * (-pkin(2) - t133) * t106 + 0.2e1 * t76 * t94, t58, t42, t54, t55, 0, -0.2e1 * t16 * t82 + 0.2e1 * t124, -0.2e1 * t76 * t101 + 0.2e1 * t123, t11, t7, t19, t20, 0, -0.2e1 * t3 * t82 + 0.2e1 * t127, 0.2e1 * t110; 0, 0, 0, 0, -t106, -t105, -t91, t119 + t122, -pkin(2) * t106 + qJ(3) * t94 + t76 * t93, t58, t42, t54, t55, 0, (-t16 - t36) * t82 + t121 + t124, (-qJ(3) - t76) * t101 + t120 + t123, t11, t7, t19, t20, 0, (-t3 - t6) * t82 + t126 + t127, t109 + t110; 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t119, 0.2e1 * qJ(3) * t93, t58, t42, t54, t55, 0, -0.2e1 * t36 * t82 + 0.2e1 * t121, -0.2e1 * t79 * t97 + 0.2e1 * t120, t11, t7, t19, t20, 0, -0.2e1 * t6 * t82 + 0.2e1 * t126, 0.2e1 * t109; 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, 0, 0, 0, 0, t73, t74, 0, 0, 0, 0, 0, t132, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t74, 0, 0, 0, 0, 0, t132, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t98, 0, t16, t15, 0, 0, -t23, -t24, 0, t63 + t3, t64 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t98, 0, t36, t35, 0, 0, -t23, -t24, 0, t63 + t6, t64 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t113, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t104, -0.2e1 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
