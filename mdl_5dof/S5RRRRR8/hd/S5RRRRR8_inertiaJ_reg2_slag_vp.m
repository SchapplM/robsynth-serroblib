% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:47
% EndTime: 2019-12-31 22:25:52
% DurationCPUTime: 1.25s
% Computational Cost: add. (1284->136), mult. (2515->259), div. (0->0), fcn. (2848->8), ass. (0->97)
t79 = sin(qJ(3));
t111 = t79 * pkin(2);
t66 = pkin(8) + t111;
t78 = sin(qJ(4));
t73 = t78 ^ 2;
t82 = cos(qJ(4));
t75 = t82 ^ 2;
t94 = t73 + t75;
t96 = t94 * t66;
t115 = -pkin(7) - pkin(6);
t84 = cos(qJ(2));
t59 = t115 * t84;
t83 = cos(qJ(3));
t80 = sin(qJ(2));
t91 = t115 * t80;
t35 = -t79 * t59 - t83 * t91;
t121 = t35 ^ 2;
t51 = t79 * t80 - t83 * t84;
t46 = t51 ^ 2;
t120 = 0.2e1 * t51;
t107 = t83 * pkin(2);
t68 = -t82 * pkin(4) - pkin(3);
t56 = t68 - t107;
t119 = 0.2e1 * t56;
t118 = 0.2e1 * t68;
t69 = -t84 * pkin(2) - pkin(1);
t117 = 0.2e1 * t69;
t116 = 0.2e1 * t84;
t77 = sin(qJ(5));
t54 = t79 * t84 + t83 * t80;
t28 = t51 * pkin(3) - t54 * pkin(8) + t69;
t38 = -t83 * t59 + t79 * t91;
t10 = t82 * t28 - t78 * t38;
t101 = t82 * t54;
t113 = t51 * pkin(4);
t8 = -pkin(9) * t101 + t10 + t113;
t81 = cos(qJ(5));
t102 = t82 * t38;
t9 = t102 + (-pkin(9) * t54 + t28) * t78;
t3 = -t77 * t9 + t81 * t8;
t109 = t81 * t9;
t4 = t77 * t8 + t109;
t49 = t77 * t78 - t81 * t82;
t53 = t77 * t82 + t81 * t78;
t114 = -t3 * t53 - t4 * t49;
t112 = t77 * pkin(4);
t110 = t81 * pkin(4);
t108 = t82 * pkin(8);
t67 = -pkin(3) - t107;
t106 = pkin(3) - t67;
t105 = t35 * t82;
t104 = t78 * t54;
t103 = t78 * t82;
t100 = t82 * t66;
t43 = (-pkin(9) - t66) * t78;
t72 = t82 * pkin(9);
t44 = t72 + t100;
t25 = t81 * t43 - t77 * t44;
t26 = t77 * t43 + t81 * t44;
t99 = -t25 * t53 - t26 * t49;
t57 = (-pkin(9) - pkin(8)) * t78;
t58 = t72 + t108;
t34 = t81 * t57 - t77 * t58;
t37 = t77 * t57 + t81 * t58;
t98 = -t34 * t53 - t37 * t49;
t97 = t56 + t68;
t95 = t94 * pkin(8);
t74 = t80 ^ 2;
t76 = t84 ^ 2;
t93 = t74 + t76;
t92 = -0.2e1 * t54 * t51;
t90 = -pkin(3) * t54 - pkin(8) * t51;
t11 = t78 * t28 + t102;
t5 = -t10 * t78 + t11 * t82;
t89 = -t51 * t66 + t54 * t67;
t62 = 0.2e1 * t103;
t48 = t54 ^ 2;
t47 = t53 ^ 2;
t45 = t49 ^ 2;
t42 = t82 * t51;
t41 = t78 * t51;
t40 = t78 * t101;
t33 = t53 * t51;
t32 = t49 * t51;
t31 = -0.2e1 * t53 * t49;
t30 = t35 * t78;
t29 = (-t73 + t75) * t54;
t27 = (-t49 * t77 - t53 * t81) * pkin(4);
t23 = t81 * t101 - t77 * t104;
t21 = t53 * t54;
t18 = pkin(4) * t104 + t35;
t15 = t23 * t53;
t14 = t21 * t49;
t13 = t18 * t53;
t12 = t18 * t49;
t6 = -t53 * t21 - t23 * t49;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t74, t80 * t116, 0, t76, 0, 0, pkin(1) * t116, -0.2e1 * pkin(1) * t80, 0.2e1 * t93 * pkin(6), t93 * pkin(6) ^ 2 + pkin(1) ^ 2, t48, t92, 0, t46, 0, 0, t51 * t117, t54 * t117, 0.2e1 * t35 * t54 - 0.2e1 * t38 * t51, t38 ^ 2 + t69 ^ 2 + t121, t75 * t48, -0.2e1 * t48 * t103, t101 * t120, t73 * t48, t78 * t92, t46, 0.2e1 * t10 * t51 + 0.2e1 * t35 * t104, 0.2e1 * t35 * t101 - 0.2e1 * t11 * t51, 0.2e1 * (-t10 * t82 - t11 * t78) * t54, t10 ^ 2 + t11 ^ 2 + t121, t23 ^ 2, -0.2e1 * t23 * t21, t23 * t120, t21 ^ 2, -t21 * t120, t46, 0.2e1 * t18 * t21 + 0.2e1 * t3 * t51, 0.2e1 * t18 * t23 - 0.2e1 * t4 * t51, -0.2e1 * t4 * t21 - 0.2e1 * t3 * t23, t18 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, t84, 0, -t80 * pkin(6), -t84 * pkin(6), 0, 0, 0, 0, t54, 0, -t51, 0, -t35, -t38, (-t51 * t79 - t54 * t83) * pkin(2), (-t35 * t83 + t38 * t79) * pkin(2), t40, t29, t41, -t40, t42, 0, t78 * t89 - t105, t82 * t89 + t30, t5, t35 * t67 + t5 * t66, t15, t6, t33, t14, -t32, 0, t56 * t21 + t25 * t51 + t12, t56 * t23 - t26 * t51 + t13, -t26 * t21 - t25 * t23 + t114, t18 * t56 + t3 * t25 + t4 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t107, -0.2e1 * t111, 0, (t79 ^ 2 + t83 ^ 2) * pkin(2) ^ 2, t73, t62, 0, t75, 0, 0, -0.2e1 * t67 * t82, 0.2e1 * t67 * t78, 0.2e1 * t96, t94 * t66 ^ 2 + t67 ^ 2, t47, t31, 0, t45, 0, 0, t49 * t119, t53 * t119, 0.2e1 * t99, t25 ^ 2 + t26 ^ 2 + t56 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t51, 0, -t35, -t38, 0, 0, t40, t29, t41, -t40, t42, 0, t78 * t90 - t105, t82 * t90 + t30, t5, -t35 * pkin(3) + pkin(8) * t5, t15, t6, t33, t14, -t32, 0, t68 * t21 + t34 * t51 + t12, t68 * t23 - t37 * t51 + t13, -t37 * t21 - t34 * t23 + t114, t18 * t68 + t3 * t34 + t4 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t107, -t111, 0, 0, t73, t62, 0, t75, 0, 0, t106 * t82, -t106 * t78, t95 + t96, -t67 * pkin(3) + pkin(8) * t96, t47, t31, 0, t45, 0, 0, t97 * t49, t97 * t53, t98 + t99, t25 * t34 + t26 * t37 + t56 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t73, t62, 0, t75, 0, 0, 0.2e1 * pkin(3) * t82, -0.2e1 * pkin(3) * t78, 0.2e1 * t95, t94 * pkin(8) ^ 2 + pkin(3) ^ 2, t47, t31, 0, t45, 0, 0, t49 * t118, t53 * t118, 0.2e1 * t98, t34 ^ 2 + t37 ^ 2 + t68 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, -t104, t51, t10, -t11, 0, 0, 0, 0, t23, 0, -t21, t51, t51 * t110 + t3, -t109 + (-t8 - t113) * t77, (-t21 * t77 - t23 * t81) * pkin(4), (t3 * t81 + t4 * t77) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, t82, 0, -t78 * t66, -t100, 0, 0, 0, 0, t53, 0, -t49, 0, t25, -t26, t27, (t25 * t81 + t26 * t77) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, t82, 0, -t78 * pkin(8), -t108, 0, 0, 0, 0, t53, 0, -t49, 0, t34, -t37, t27, (t34 * t81 + t37 * t77) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t110, -0.2e1 * t112, 0, (t77 ^ 2 + t81 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t21, t51, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t49, 0, t25, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t49, 0, t34, -t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t110, -t112, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
