% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:17
% EndTime: 2019-12-31 21:17:21
% DurationCPUTime: 1.09s
% Computational Cost: add. (1118->119), mult. (2199->230), div. (0->0), fcn. (2497->8), ass. (0->89)
t77 = sin(pkin(9));
t73 = t77 ^ 2;
t78 = cos(pkin(9));
t74 = t78 ^ 2;
t93 = t73 + t74;
t94 = t93 * qJ(4);
t104 = cos(qJ(5));
t79 = sin(qJ(5));
t116 = t104 * t78 - t79 * t77;
t109 = -pkin(7) - pkin(6);
t83 = cos(qJ(2));
t61 = t109 * t83;
t80 = sin(qJ(3));
t82 = cos(qJ(3));
t81 = sin(qJ(2));
t90 = t109 * t81;
t36 = -t80 * t61 - t82 * t90;
t115 = t36 ^ 2;
t52 = t80 * t81 - t82 * t83;
t50 = t52 ^ 2;
t114 = 0.2e1 * t52;
t106 = t82 * pkin(2);
t66 = -t78 * pkin(4) - pkin(3);
t56 = t66 - t106;
t113 = 0.2e1 * t56;
t112 = 0.2e1 * t66;
t69 = -t83 * pkin(2) - pkin(1);
t111 = 0.2e1 * t69;
t110 = 0.2e1 * t83;
t54 = t80 * t83 + t82 * t81;
t100 = t78 * t54;
t27 = t52 * pkin(3) - t54 * qJ(4) + t69;
t38 = -t82 * t61 + t80 * t90;
t9 = t78 * t27 - t77 * t38;
t7 = t52 * pkin(4) - pkin(8) * t100 + t9;
t10 = t77 * t27 + t78 * t38;
t102 = t77 * t54;
t8 = -pkin(8) * t102 + t10;
t3 = t104 * t7 - t79 * t8;
t4 = t104 * t8 + t79 * t7;
t49 = t104 * t77 + t79 * t78;
t108 = t116 * t4 - t3 * t49;
t107 = t80 * pkin(2);
t68 = -pkin(3) - t106;
t105 = pkin(3) - t68;
t103 = t36 * t78;
t101 = t77 * t78;
t65 = qJ(4) + t107;
t43 = (-pkin(8) - t65) * t77;
t72 = t78 * pkin(8);
t44 = t78 * t65 + t72;
t25 = t104 * t43 - t79 * t44;
t26 = t104 * t44 + t79 * t43;
t98 = t116 * t26 - t25 * t49;
t57 = (-pkin(8) - qJ(4)) * t77;
t58 = t78 * qJ(4) + t72;
t32 = t104 * t57 - t79 * t58;
t33 = t104 * t58 + t79 * t57;
t97 = t116 * t33 - t32 * t49;
t96 = t56 + t66;
t95 = t93 * t65;
t75 = t81 ^ 2;
t76 = t83 ^ 2;
t92 = t75 + t76;
t91 = -0.2e1 * t54 * t52;
t5 = t10 * t78 - t9 * t77;
t88 = -pkin(3) * t54 - qJ(4) * t52;
t87 = -t52 * t65 + t54 * t68;
t62 = 0.2e1 * t101;
t51 = t54 ^ 2;
t46 = t49 ^ 2;
t45 = t116 ^ 2;
t42 = t78 * t52;
t41 = t77 * t52;
t39 = t77 * t100;
t35 = t49 * t52;
t34 = t116 * t52;
t31 = 0.2e1 * t49 * t116;
t30 = t36 * t77;
t28 = (-t73 + t74) * t54;
t22 = t116 * t54;
t20 = t49 * t54;
t17 = pkin(4) * t102 + t36;
t14 = t22 * t49;
t13 = t20 * t116;
t12 = t17 * t49;
t11 = t17 * t116;
t6 = t116 * t22 - t49 * t20;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t75, t81 * t110, 0, t76, 0, 0, pkin(1) * t110, -0.2e1 * pkin(1) * t81, 0.2e1 * t92 * pkin(6), t92 * pkin(6) ^ 2 + pkin(1) ^ 2, t51, t91, 0, t50, 0, 0, t52 * t111, t54 * t111, 0.2e1 * t36 * t54 - 0.2e1 * t38 * t52, t38 ^ 2 + t69 ^ 2 + t115, t74 * t51, -0.2e1 * t51 * t101, t100 * t114, t73 * t51, t77 * t91, t50, 0.2e1 * t36 * t102 + 0.2e1 * t9 * t52, -0.2e1 * t10 * t52 + 0.2e1 * t36 * t100, 0.2e1 * (-t10 * t77 - t78 * t9) * t54, t10 ^ 2 + t9 ^ 2 + t115, t22 ^ 2, -0.2e1 * t22 * t20, t22 * t114, t20 ^ 2, -t20 * t114, t50, 0.2e1 * t17 * t20 + 0.2e1 * t3 * t52, 0.2e1 * t17 * t22 - 0.2e1 * t4 * t52, -0.2e1 * t4 * t20 - 0.2e1 * t3 * t22, t17 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, t83, 0, -t81 * pkin(6), -t83 * pkin(6), 0, 0, 0, 0, t54, 0, -t52, 0, -t36, -t38, (-t52 * t80 - t54 * t82) * pkin(2), (-t36 * t82 + t38 * t80) * pkin(2), t39, t28, t41, -t39, t42, 0, t87 * t77 - t103, t87 * t78 + t30, t5, t36 * t68 + t5 * t65, t14, t6, t35, -t13, t34, 0, t56 * t20 + t25 * t52 - t11, t56 * t22 - t26 * t52 + t12, -t26 * t20 - t25 * t22 + t108, t17 * t56 + t3 * t25 + t4 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t106, -0.2e1 * t107, 0, (t80 ^ 2 + t82 ^ 2) * pkin(2) ^ 2, t73, t62, 0, t74, 0, 0, -0.2e1 * t68 * t78, 0.2e1 * t68 * t77, 0.2e1 * t95, t93 * t65 ^ 2 + t68 ^ 2, t46, t31, 0, t45, 0, 0, -t116 * t113, t49 * t113, 0.2e1 * t98, t25 ^ 2 + t26 ^ 2 + t56 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t52, 0, -t36, -t38, 0, 0, t39, t28, t41, -t39, t42, 0, t88 * t77 - t103, t88 * t78 + t30, t5, -t36 * pkin(3) + t5 * qJ(4), t14, t6, t35, -t13, t34, 0, t66 * t20 + t32 * t52 - t11, t66 * t22 - t33 * t52 + t12, -t33 * t20 - t32 * t22 + t108, t17 * t66 + t3 * t32 + t4 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t106, -t107, 0, 0, t73, t62, 0, t74, 0, 0, t105 * t78, -t105 * t77, t94 + t95, -t68 * pkin(3) + t65 * t94, t46, t31, 0, t45, 0, 0, -t96 * t116, t96 * t49, t97 + t98, t25 * t32 + t26 * t33 + t56 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t73, t62, 0, t74, 0, 0, 0.2e1 * pkin(3) * t78, -0.2e1 * pkin(3) * t77, 0.2e1 * t94, t93 * qJ(4) ^ 2 + pkin(3) ^ 2, t46, t31, 0, t45, 0, 0, -t116 * t112, t49 * t112, 0.2e1 * t97, t32 ^ 2 + t33 ^ 2 + t66 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t100, 0, t36, 0, 0, 0, 0, 0, 0, t20, t22, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t77, 0, t68, 0, 0, 0, 0, 0, 0, -t116, t49, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t77, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t116, t49, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t20, t52, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t116, 0, t25, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t116, 0, t32, -t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
