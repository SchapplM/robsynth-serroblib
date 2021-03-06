% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR13_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:33
% EndTime: 2019-12-31 21:46:38
% DurationCPUTime: 1.31s
% Computational Cost: add. (1022->150), mult. (2393->280), div. (0->0), fcn. (2565->8), ass. (0->101)
t59 = sin(qJ(3));
t52 = t59 ^ 2;
t62 = cos(qJ(3));
t54 = t62 ^ 2;
t113 = t52 + t54;
t57 = cos(pkin(5));
t56 = sin(pkin(5));
t60 = sin(qJ(2));
t99 = t56 * t60;
t27 = -t57 * t62 + t59 * t99;
t112 = t27 ^ 2;
t29 = t57 * t59 + t62 * t99;
t26 = t29 ^ 2;
t111 = -0.2e1 * t29;
t110 = 0.2e1 * t56;
t109 = -0.2e1 * t59;
t108 = 0.2e1 * t62;
t107 = 2 * qJ(4);
t106 = pkin(3) + pkin(9);
t105 = pkin(1) * t60;
t63 = cos(qJ(2));
t104 = pkin(1) * t63;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t98 = t56 * t63;
t14 = t27 * t61 + t58 * t98;
t103 = t14 * t58;
t102 = t27 * t62;
t101 = t29 * t58;
t24 = t29 * t59;
t50 = t56 ^ 2;
t100 = t50 * t63;
t97 = t58 * t59;
t96 = t58 * t62;
t95 = t58 * t106;
t94 = t59 * t62;
t15 = -t27 * t58 + t61 * t98;
t93 = t61 * t15;
t92 = t61 * t58;
t91 = t61 * t62;
t90 = t61 * t106;
t84 = pkin(7) * t98;
t21 = t84 + (pkin(8) + t105) * t57;
t22 = (-pkin(2) * t63 - pkin(8) * t60 - pkin(1)) * t56;
t11 = t62 * t21 + t59 * t22;
t89 = t113 * pkin(8) ^ 2;
t51 = t58 ^ 2;
t53 = t61 ^ 2;
t39 = t51 + t53;
t88 = qJ(4) * t62;
t87 = t27 * t111;
t86 = t57 * t110;
t85 = -0.2e1 * t94;
t83 = t29 * t98;
t82 = t27 * t98;
t81 = t59 * t98;
t80 = t62 * t98;
t79 = t58 * t91;
t78 = qJ(4) * t98;
t77 = -t59 * qJ(4) - pkin(2);
t10 = -t59 * t21 + t62 * t22;
t76 = pkin(8) * t81;
t75 = pkin(8) * t80;
t41 = pkin(3) * t98;
t9 = -t10 + t41;
t3 = t29 * pkin(4) + pkin(9) * t98 + t9;
t40 = pkin(7) * t99;
t20 = t40 + (-pkin(2) - t104) * t57;
t67 = -t29 * qJ(4) + t20;
t4 = t106 * t27 + t67;
t1 = t61 * t3 - t58 * t4;
t2 = t58 * t3 + t61 * t4;
t74 = t1 * t61 + t2 * t58;
t8 = t78 - t11;
t73 = t9 * t59 - t8 * t62;
t72 = -pkin(3) * t59 + t88;
t71 = -t10 * t59 + t11 * t62;
t33 = -t106 * t62 + t77;
t47 = t59 * pkin(8);
t37 = t59 * pkin(4) + t47;
t12 = -t58 * t33 + t61 * t37;
t13 = t61 * t33 + t58 * t37;
t6 = t12 * t61 + t13 * t58;
t70 = -t59 * t27 + t29 * t62;
t69 = -t106 * t59 + t88;
t68 = (t24 - t102) * pkin(8);
t65 = qJ(4) ^ 2;
t49 = t62 * pkin(8);
t45 = t61 * t59;
t43 = t50 * t63 ^ 2;
t42 = 0.2e1 * t94;
t38 = t62 * pkin(4) + t49;
t36 = -t62 * pkin(3) + t77;
t35 = 0.2e1 * t113 * pkin(8);
t34 = t39 * t106;
t32 = t57 * t105 + t84;
t31 = t57 * t104 - t40;
t25 = t29 * t61;
t7 = t27 * pkin(3) + t67;
t5 = -t27 * pkin(4) - t8;
t16 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t50 * t60 ^ 2, 0.2e1 * t60 * t100, t60 * t86, t43, t63 * t86, t57 ^ 2, 0.2e1 * pkin(1) * t100 + 0.2e1 * t31 * t57, -0.2e1 * t50 * t105 - 0.2e1 * t32 * t57, (-t31 * t60 + t32 * t63) * t110, t50 * pkin(1) ^ 2 + t31 ^ 2 + t32 ^ 2, t26, t87, -0.2e1 * t83, t112, 0.2e1 * t82, t43, -0.2e1 * t10 * t98 + 0.2e1 * t20 * t27, 0.2e1 * t11 * t98 + 0.2e1 * t20 * t29, -0.2e1 * t10 * t29 - 0.2e1 * t11 * t27, t10 ^ 2 + t11 ^ 2 + t20 ^ 2, t43, 0.2e1 * t83, -0.2e1 * t82, t26, t87, t112, 0.2e1 * t8 * t27 + 0.2e1 * t9 * t29, -0.2e1 * t7 * t27 - 0.2e1 * t9 * t98, -0.2e1 * t7 * t29 + 0.2e1 * t8 * t98, t7 ^ 2 + t8 ^ 2 + t9 ^ 2, t15 ^ 2, -0.2e1 * t15 * t14, t15 * t111, t14 ^ 2, 0.2e1 * t14 * t29, t26, 0.2e1 * t1 * t29 - 0.2e1 * t5 * t14, -0.2e1 * t5 * t15 - 0.2e1 * t2 * t29, 0.2e1 * t1 * t15 + 0.2e1 * t2 * t14, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, t98, t57, t31, -t32, 0, 0, t24, t70, -t81, -t102, -t80, 0, -pkin(2) * t27 - t20 * t62 + t76, -pkin(2) * t29 + t20 * t59 + t75, t68 + t71, -t20 * pkin(2) + t71 * pkin(8), 0, t81, t80, t24, t70, -t102, t68 + t73, -t36 * t27 + t7 * t62 - t76, -t36 * t29 - t7 * t59 - t75, t73 * pkin(8) + t7 * t36, t15 * t96, (t93 - t103) * t62, -t15 * t59 - t29 * t96, -t14 * t91, t14 * t59 - t29 * t91, t24, t1 * t59 + t12 * t29 - t38 * t14 + t5 * t91, -t13 * t29 - t38 * t15 - t2 * t59 - t5 * t96, t12 * t15 + t13 * t14 + (t1 * t58 - t2 * t61) * t62, t1 * t12 + t2 * t13 + t5 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t52, t42, 0, t54, 0, 0, pkin(2) * t108, pkin(2) * t109, t35, pkin(2) ^ 2 + t89, 0, 0, 0, t52, t42, t54, t35, t36 * t108, t36 * t109, t36 ^ 2 + t89, t51 * t54, 0.2e1 * t54 * t92, t58 * t85, t53 * t54, t61 * t85, t52, 0.2e1 * t12 * t59 + 0.2e1 * t38 * t91, -0.2e1 * t13 * t59 - 0.2e1 * t38 * t96, (t12 * t58 - t13 * t61) * t108, t12 ^ 2 + t13 ^ 2 + t38 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t27, -t98, t10, -t11, 0, 0, -t98, -t29, t27, 0, 0, 0, -t29 * pkin(3) - qJ(4) * t27, -t10 + 0.2e1 * t41, -0.2e1 * t78 + t11, -t9 * pkin(3) - t8 * qJ(4), -t93, t61 * t14 + t15 * t58, t25, -t103, -t101, 0, -qJ(4) * t14 - t29 * t90 + t5 * t58, -qJ(4) * t15 + t29 * t95 + t5 * t61, (-t106 * t15 - t1) * t61 + (-t106 * t14 - t2) * t58, t5 * qJ(4) - t106 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t62, 0, -t47, -t49, 0, 0, 0, -t59, -t62, 0, 0, 0, t72, t47, t49, t72 * pkin(8), -t79, (t51 - t53) * t62, t45, t79, -t97, 0, t38 * t58 + t69 * t61, t38 * t61 - t69 * t58, -t6, t38 * qJ(4) - t106 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3), t107, pkin(3) ^ 2 + t65, t53, -0.2e1 * t92, 0, t51, 0, 0, t58 * t107, t61 * t107, 0.2e1 * t34, t106 ^ 2 * t39 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t98, 0, t9, 0, 0, 0, 0, 0, 0, t25, -t101, t93 + t103, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, t47, 0, 0, 0, 0, 0, 0, t45, -t97, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, t14, t29, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, 0, -t91, t59, t12, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t58, 0, -t90, t95, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t16;
