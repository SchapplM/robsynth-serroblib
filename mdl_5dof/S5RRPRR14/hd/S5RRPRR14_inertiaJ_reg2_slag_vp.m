% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR14_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:37
% EndTime: 2019-12-31 20:38:42
% DurationCPUTime: 1.26s
% Computational Cost: add. (1732->144), mult. (4035->310), div. (0->0), fcn. (4600->10), ass. (0->95)
t56 = sin(pkin(10));
t58 = cos(pkin(10));
t61 = sin(qJ(4));
t98 = cos(qJ(4));
t44 = t98 * t56 + t61 * t58;
t108 = -0.2e1 * t44;
t59 = cos(pkin(5));
t57 = sin(pkin(5));
t62 = sin(qJ(2));
t90 = t57 * t62;
t34 = t56 * t90 - t59 * t58;
t36 = t59 * t56 + t58 * t90;
t19 = t98 * t34 + t61 * t36;
t107 = t19 ^ 2;
t82 = pkin(8) + qJ(3);
t45 = t82 * t58;
t73 = t82 * t56;
t24 = t61 * t45 + t98 * t73;
t106 = t24 ^ 2;
t42 = t61 * t56 - t98 * t58;
t105 = t42 ^ 2;
t104 = -0.2e1 * t19;
t50 = -t58 * pkin(3) - pkin(2);
t103 = 0.2e1 * t50;
t102 = 0.2e1 * t57;
t101 = 0.2e1 * t58;
t100 = pkin(1) * t62;
t64 = cos(qJ(2));
t99 = pkin(1) * t64;
t21 = -t61 * t34 + t98 * t36;
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t89 = t57 * t64;
t13 = t60 * t21 + t63 * t89;
t97 = t13 * t63;
t15 = t63 * t21 - t60 * t89;
t96 = t15 * t60;
t95 = t15 * t63;
t94 = t19 * t42;
t93 = t34 * t58;
t92 = t36 * t56;
t52 = t57 ^ 2;
t91 = t52 * t64;
t88 = t59 * t62;
t87 = t60 * t19;
t86 = t60 * t42;
t85 = t60 * t44;
t84 = t60 * t63;
t18 = t63 * t19;
t83 = t63 * t44;
t76 = pkin(7) * t89;
t30 = t76 + (qJ(3) + t100) * t59;
t31 = (-pkin(2) * t64 - qJ(3) * t62 - pkin(1)) * t57;
t17 = t58 * t30 + t56 * t31;
t51 = t56 ^ 2;
t53 = t58 ^ 2;
t81 = t51 + t53;
t54 = t60 ^ 2;
t55 = t63 ^ 2;
t80 = t54 + t55;
t79 = t42 * t108;
t78 = -0.2e1 * t89;
t77 = 0.2e1 * t89;
t75 = t60 * t83;
t74 = qJ(3) * t89;
t16 = -t56 * t30 + t58 * t31;
t72 = -pkin(4) * t44 - pkin(9) * t42;
t11 = -t34 * pkin(8) + t17;
t8 = -pkin(3) * t89 - t36 * pkin(8) + t16;
t6 = t98 * t11 + t61 * t8;
t4 = -pkin(9) * t89 + t6;
t47 = pkin(7) * t90;
t33 = t47 + (-pkin(2) - t99) * t59;
t22 = t34 * pkin(3) + t33;
t7 = t19 * pkin(4) - t21 * pkin(9) + t22;
t1 = -t60 * t4 + t63 * t7;
t2 = t63 * t4 + t60 * t7;
t71 = t1 * t63 + t2 * t60;
t70 = -t1 * t60 + t2 * t63;
t23 = t42 * pkin(4) - t44 * pkin(9) + t50;
t26 = t98 * t45 - t61 * t73;
t10 = t60 * t23 + t63 * t26;
t9 = t63 * t23 - t60 * t26;
t69 = t10 * t63 - t9 * t60;
t68 = t10 * t60 + t9 * t63;
t67 = -t16 * t56 + t17 * t58;
t5 = -t61 * t11 + t98 * t8;
t48 = t52 * t64 ^ 2;
t40 = t44 ^ 2;
t39 = pkin(1) * t88 + t76;
t38 = t59 * t99 - t47;
t37 = t63 * t42;
t12 = t60 * t13;
t3 = pkin(4) * t89 - t5;
t14 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t52 * t62 ^ 2, 0.2e1 * t62 * t91, t88 * t102, t48, t59 * t77, t59 ^ 2, 0.2e1 * pkin(1) * t91 + 0.2e1 * t38 * t59, -0.2e1 * t52 * t100 - 0.2e1 * t39 * t59, (-t38 * t62 + t39 * t64) * t102, t52 * pkin(1) ^ 2 + t38 ^ 2 + t39 ^ 2, t36 ^ 2, -0.2e1 * t36 * t34, t36 * t78, t34 ^ 2, t34 * t77, t48, -0.2e1 * t16 * t89 + 0.2e1 * t33 * t34, 0.2e1 * t17 * t89 + 0.2e1 * t33 * t36, -0.2e1 * t16 * t36 - 0.2e1 * t17 * t34, t16 ^ 2 + t17 ^ 2 + t33 ^ 2, t21 ^ 2, t21 * t104, t21 * t78, t107, t19 * t77, t48, 0.2e1 * t22 * t19 - 0.2e1 * t5 * t89, 0.2e1 * t22 * t21 + 0.2e1 * t6 * t89, -0.2e1 * t6 * t19 - 0.2e1 * t5 * t21, t22 ^ 2 + t5 ^ 2 + t6 ^ 2, t15 ^ 2, -0.2e1 * t15 * t13, 0.2e1 * t15 * t19, t13 ^ 2, t13 * t104, t107, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t13, 0.2e1 * t3 * t15 - 0.2e1 * t2 * t19, -0.2e1 * t1 * t15 - 0.2e1 * t2 * t13, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, 0, t89, t59, t38, -t39, 0, 0, t92, -t56 * t34 + t36 * t58, -t56 * t89, -t93, -t58 * t89, 0, -pkin(2) * t34 - t33 * t58 + t56 * t74, -pkin(2) * t36 + t33 * t56 + t58 * t74, (t92 - t93) * qJ(3) + t67, -t33 * pkin(2) + t67 * qJ(3), t21 * t44, -t44 * t19 - t21 * t42, -t44 * t89, t94, t42 * t89, 0, t50 * t19 + t22 * t42 + t24 * t89, t50 * t21 + t22 * t44 + t26 * t89, -t26 * t19 + t24 * t21 - t6 * t42 - t5 * t44, t22 * t50 - t5 * t24 + t6 * t26, t15 * t83, (-t96 - t97) * t44, t15 * t42 + t19 * t83, t13 * t85, -t13 * t42 - t19 * t85, t94, t1 * t42 + t24 * t13 + t9 * t19 + t3 * t85, -t10 * t19 + t24 * t15 - t2 * t42 + t3 * t83, -t10 * t13 - t9 * t15 - t71 * t44, t1 * t9 + t2 * t10 + t3 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t51, t56 * t101, 0, t53, 0, 0, pkin(2) * t101, -0.2e1 * pkin(2) * t56, 0.2e1 * t81 * qJ(3), t81 * qJ(3) ^ 2 + pkin(2) ^ 2, t40, t79, 0, t105, 0, 0, t42 * t103, t44 * t103, 0.2e1 * t24 * t44 - 0.2e1 * t26 * t42, t26 ^ 2 + t50 ^ 2 + t106, t55 * t40, -0.2e1 * t40 * t84, 0.2e1 * t42 * t83, t54 * t40, t60 * t79, t105, 0.2e1 * t24 * t85 + 0.2e1 * t9 * t42, -0.2e1 * t10 * t42 + 0.2e1 * t24 * t83, t68 * t108, t10 ^ 2 + t9 ^ 2 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t36, 0, t33, 0, 0, 0, 0, 0, 0, t19, t21, 0, t22, 0, 0, 0, 0, 0, 0, t18, -t87, -t12 - t95, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t56, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t42, t44, 0, t50, 0, 0, 0, 0, 0, 0, t37, -t86, -t80 * t44, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, -t89, t5, -t6, 0, 0, t96, -t12 + t95, t87, -t97, t18, 0, -pkin(4) * t13 - pkin(9) * t87 - t3 * t63, -pkin(4) * t15 - pkin(9) * t18 + t3 * t60, (t96 - t97) * pkin(9) + t70, -t3 * pkin(4) + t70 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t42, 0, -t24, -t26, 0, 0, t75, (-t54 + t55) * t44, t86, -t75, t37, 0, -t24 * t63 + t72 * t60, t24 * t60 + t72 * t63, t69, -t24 * pkin(4) + t69 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t54, 0.2e1 * t84, 0, t55, 0, 0, 0.2e1 * pkin(4) * t63, -0.2e1 * pkin(4) * t60, 0.2e1 * t80 * pkin(9), t80 * pkin(9) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t13, t19, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, -t85, t42, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t63, 0, -t60 * pkin(9), -t63 * pkin(9), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t14;
