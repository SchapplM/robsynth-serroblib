% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:43:42
% EndTime: 2019-05-05 21:43:47
% DurationCPUTime: 1.68s
% Computational Cost: add. (981->145), mult. (1901->262), div. (0->0), fcn. (2048->6), ass. (0->92)
t60 = sin(pkin(9));
t61 = cos(pkin(9));
t62 = sin(qJ(4));
t64 = cos(qJ(4));
t38 = t60 * t64 + t61 * t62;
t65 = cos(qJ(3));
t93 = t65 * t38;
t117 = t93 ^ 2;
t36 = t60 * t62 - t61 * t64;
t116 = t36 ^ 2;
t52 = -t64 * pkin(4) - pkin(3);
t115 = 0.2e1 * t52;
t114 = 0.2e1 * t64;
t113 = 2 * qJ(2);
t112 = t60 * pkin(4);
t111 = t61 * pkin(4);
t63 = sin(qJ(3));
t110 = t63 * pkin(4);
t109 = t65 * pkin(3);
t66 = -pkin(1) - pkin(7);
t100 = t62 * t66;
t42 = t63 * pkin(3) - t65 * pkin(8) + qJ(2);
t35 = t64 * t42;
t88 = qJ(5) * t65;
t13 = -t64 * t88 + t35 + (pkin(4) - t100) * t63;
t95 = t64 * t66;
t85 = t63 * t95;
t17 = t85 + (t42 - t88) * t62;
t4 = t60 * t13 + t61 * t17;
t90 = -qJ(5) - pkin(8);
t43 = t90 * t64;
t81 = t90 * t62;
t19 = -t60 * t43 - t61 * t81;
t108 = t19 * t63;
t21 = -t61 * t43 + t60 * t81;
t107 = t21 * t63;
t106 = t93 * t36;
t101 = t62 * t65;
t51 = t64 * t65;
t32 = -t60 * t101 + t61 * t51;
t105 = t32 * t93;
t104 = t38 * t36;
t33 = t38 * t63;
t103 = t62 * t63;
t102 = t62 * t64;
t99 = t63 * t93;
t98 = t63 * t36;
t97 = t63 * t66;
t96 = t64 * t63;
t94 = t65 * t36;
t92 = t65 * t63;
t91 = t65 * t66;
t55 = t62 ^ 2;
t57 = t64 ^ 2;
t89 = t55 + t57;
t56 = t63 ^ 2;
t58 = t65 ^ 2;
t45 = t56 + t58;
t87 = -0.2e1 * t92;
t86 = t19 ^ 2 + t21 ^ 2;
t1 = t63 * qJ(6) + t4;
t84 = t62 * t51;
t83 = t19 * t33 - t21 * t98;
t82 = t19 * t32 - t21 * t93;
t80 = t89 * t63;
t79 = -t61 * t13 + t60 * t17;
t78 = t32 * t33 + t93 * t98;
t77 = t33 * t38 + t36 * t98;
t39 = pkin(4) * t101 - t91;
t76 = t33 ^ 2 + t98 ^ 2 + t58;
t75 = -pkin(8) * t63 - t109;
t22 = -t62 * t97 + t35;
t23 = t62 * t42 + t85;
t74 = -t22 * t62 + t23 * t64;
t73 = -t33 * t63 - t65 * t93;
t72 = -t32 * t36 - t38 * t93;
t71 = t65 * t32 - t63 * t98;
t70 = 0.2e1 * t19 * t38 - 0.2e1 * t21 * t36;
t67 = qJ(2) ^ 2;
t59 = t66 ^ 2;
t53 = t58 * t59;
t49 = pkin(5) + t111;
t47 = qJ(6) + t112;
t41 = t45 * t66;
t34 = t38 ^ 2;
t26 = t32 ^ 2;
t24 = 0.2e1 * t32 * t63;
t16 = t32 * t38;
t12 = t36 * pkin(5) - t38 * qJ(6) + t52;
t5 = pkin(5) * t93 - t32 * qJ(6) + t39;
t2 = -t63 * pkin(5) + t79;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t113, pkin(1) ^ 2 + t67, t58, t87, 0, t56, 0, 0, t63 * t113, t65 * t113, -0.2e1 * t41, t56 * t59 + t53 + t67, t57 * t58, -0.2e1 * t58 * t102, t92 * t114, t55 * t58, t62 * t87, t56, -0.2e1 * t58 * t100 + 0.2e1 * t22 * t63, -0.2e1 * t23 * t63 - 0.2e1 * t58 * t95, 0.2e1 * (-t22 * t64 - t23 * t62) * t65, t22 ^ 2 + t23 ^ 2 + t53, t26, -0.2e1 * t105, t24, t117, -0.2e1 * t99, t56, 0.2e1 * t39 * t93 - 0.2e1 * t63 * t79, 0.2e1 * t39 * t32 - 0.2e1 * t4 * t63, 0.2e1 * t32 * t79 - 0.2e1 * t4 * t93, t39 ^ 2 + t4 ^ 2 + t79 ^ 2, t26, t24, 0.2e1 * t105, t56, 0.2e1 * t99, t117, -0.2e1 * t2 * t63 + 0.2e1 * t5 * t93, -0.2e1 * t1 * t93 + 0.2e1 * t2 * t32, 0.2e1 * t1 * t63 - 0.2e1 * t5 * t32, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t45, t41, 0, 0, 0, 0, 0, 0, -t45 * t62, -t45 * t64, 0, t58 * t66 + t74 * t63, 0, 0, 0, 0, 0, 0, t73, -t71, t78, t33 * t79 - t39 * t65 - t4 * t98, 0, 0, 0, 0, 0, 0, t73, t78, t71, -t1 * t98 + t2 * t33 - t5 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t56 + t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t63, 0, t91, -t97, 0, 0, t84 (-t55 + t57) * t65, t103, -t84, t96, 0, t75 * t62 + t64 * t91, -t62 * t91 + t75 * t64, t74, pkin(3) * t91 + t74 * pkin(8), t16, t72, t33, t106, -t98, 0, t39 * t36 + t52 * t93 - t108, t52 * t32 + t39 * t38 - t107, -t4 * t36 + t38 * t79 + t82, t19 * t79 + t4 * t21 + t39 * t52, t16, t33, -t72, 0, t98, t106, t12 * t93 + t5 * t36 - t108, -t1 * t36 + t2 * t38 + t82, -t12 * t32 - t5 * t38 + t107, t1 * t21 + t5 * t12 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t63, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t101, t80, pkin(8) * t80 + t109, 0, 0, 0, 0, 0, 0, -t94, -t93, t77, -t65 * t52 + t83, 0, 0, 0, 0, 0, 0, -t94, t77, t93, -t65 * t12 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t55, 0.2e1 * t102, 0, t57, 0, 0, pkin(3) * t114, -0.2e1 * pkin(3) * t62, 0.2e1 * t89 * pkin(8), t89 * pkin(8) ^ 2 + pkin(3) ^ 2, t34, -0.2e1 * t104, 0, t116, 0, 0, t36 * t115, t38 * t115, t70, t52 ^ 2 + t86, t34, 0, 0.2e1 * t104, 0, 0, t116, 0.2e1 * t12 * t36, t70, -0.2e1 * t12 * t38, t12 ^ 2 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, -t101, t63, t22, -t23, 0, 0, 0, 0, t32, 0, -t93, t63, t61 * t110 - t79, -t60 * t110 - t4 (-t32 * t61 - t60 * t93) * pkin(4) (t4 * t60 - t61 * t79) * pkin(4), 0, t32, 0, t63, t93, 0 (pkin(5) + t49) * t63 - t79, -t49 * t32 - t47 * t93, t47 * t63 + t1, t1 * t47 - t2 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t96, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t98, 0 (-t33 * t61 - t60 * t98) * pkin(4), 0, 0, 0, 0, 0, 0, -t33, 0, -t98, -t33 * t49 - t47 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t64, 0, -t62 * pkin(8), -t64 * pkin(8), 0, 0, 0, 0, t38, 0, -t36, 0, -t19, -t21 (-t36 * t60 - t38 * t61) * pkin(4) (-t19 * t61 + t21 * t60) * pkin(4), 0, t38, 0, 0, t36, 0, -t19, -t47 * t36 - t49 * t38, t21, -t19 * t49 + t21 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t111, -0.2e1 * t112, 0 (t60 ^ 2 + t61 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t49, 0, 0.2e1 * t47, t47 ^ 2 + t49 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t32, 0, t39, 0, 0, 0, 0, 0, 0, t93, 0, -t32, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t38, 0, t52, 0, 0, 0, 0, 0, 0, t36, 0, -t38, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t32, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
