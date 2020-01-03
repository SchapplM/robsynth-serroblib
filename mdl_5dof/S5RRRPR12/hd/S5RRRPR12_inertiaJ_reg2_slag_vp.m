% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t112 = cos(qJ(5));
t69 = sin(pkin(10));
t71 = cos(pkin(10));
t73 = sin(qJ(5));
t122 = t112 * t71 - t73 * t69;
t70 = sin(pkin(5));
t75 = sin(qJ(2));
t101 = t70 * t75;
t72 = cos(pkin(5));
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t39 = t74 * t101 - t72 * t76;
t38 = t39 ^ 2;
t121 = -0.2e1 * t39;
t120 = 0.2e1 * t39;
t41 = t76 * t101 + t72 * t74;
t119 = -0.2e1 * t41;
t61 = -t71 * pkin(4) - pkin(3);
t118 = 0.2e1 * t61;
t117 = 0.2e1 * t70;
t116 = -0.2e1 * t76;
t115 = pkin(1) * t75;
t77 = cos(qJ(2));
t114 = pkin(1) * t77;
t113 = pkin(8) * t69;
t63 = t74 * pkin(8);
t100 = t70 * t77;
t88 = pkin(7) * t100;
t33 = t88 + (pkin(8) + t115) * t72;
t34 = (-pkin(2) * t77 - pkin(8) * t75 - pkin(1)) * t70;
t19 = -t74 * t33 + t76 * t34;
t18 = pkin(3) * t100 - t19;
t111 = t18 * t69;
t110 = t18 * t71;
t23 = t71 * t100 + t41 * t69;
t109 = t23 * t71;
t25 = -t69 * t100 + t41 * t71;
t108 = t25 * t69;
t107 = t39 * t76;
t106 = t41 * t74;
t65 = t70 ^ 2;
t105 = t65 * t77;
t104 = t69 * t71;
t103 = t69 * t74;
t102 = t69 * t76;
t99 = t71 * t74;
t98 = t71 * t76;
t97 = t72 * t75;
t95 = t74 * t39;
t94 = t74 * t76;
t93 = pkin(9) + qJ(4);
t55 = pkin(7) * t101;
t32 = t55 + (-pkin(2) - t114) * t72;
t16 = t39 * pkin(3) - t41 * qJ(4) + t32;
t20 = t76 * t33 + t74 * t34;
t17 = -qJ(4) * t100 + t20;
t6 = t69 * t16 + t71 * t17;
t50 = -t76 * pkin(3) - t74 * qJ(4) - pkin(2);
t31 = pkin(8) * t98 + t69 * t50;
t64 = t69 ^ 2;
t66 = t71 ^ 2;
t92 = t64 + t66;
t91 = qJ(4) * t39;
t90 = 0.2e1 * t100;
t89 = 0.2e1 * t94;
t87 = t74 * t100;
t86 = t76 * t100;
t85 = t69 * t99;
t5 = t71 * t16 - t69 * t17;
t83 = -t5 * t69 + t6 * t71;
t82 = -pkin(3) * t74 + qJ(4) * t76;
t81 = -t19 * t74 + t20 * t76;
t45 = t71 * t50;
t30 = -pkin(8) * t102 + t45;
t80 = -t30 * t69 + t31 * t71;
t48 = t112 * t69 + t73 * t71;
t79 = pkin(8) ^ 2;
t68 = t76 ^ 2;
t67 = t74 ^ 2;
t62 = t67 * t79;
t58 = t65 * t77 ^ 2;
t52 = t93 * t71;
t51 = t93 * t69;
t49 = pkin(4) * t103 + t63;
t43 = pkin(1) * t97 + t88;
t42 = t72 * t114 - t55;
t37 = t122 * t74;
t35 = t48 * t74;
t28 = t112 * t52 - t73 * t51;
t27 = -t112 * t51 - t73 * t52;
t26 = -pkin(9) * t103 + t31;
t21 = -pkin(9) * t99 + t45 + (-pkin(4) - t113) * t76;
t12 = t112 * t25 - t73 * t23;
t10 = t112 * t23 + t73 * t25;
t9 = t112 * t26 + t73 * t21;
t8 = t112 * t21 - t73 * t26;
t7 = t23 * pkin(4) + t18;
t4 = -t23 * pkin(9) + t6;
t3 = t39 * pkin(4) - t25 * pkin(9) + t5;
t2 = t112 * t4 + t73 * t3;
t1 = t112 * t3 - t73 * t4;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t65 * t75 ^ 2, 0.2e1 * t75 * t105, t97 * t117, t58, t72 * t90, t72 ^ 2, 0.2e1 * pkin(1) * t105 + 0.2e1 * t42 * t72, -0.2e1 * t65 * t115 - 0.2e1 * t43 * t72, (-t42 * t75 + t43 * t77) * t117, t65 * pkin(1) ^ 2 + t42 ^ 2 + t43 ^ 2, t41 ^ 2, t39 * t119, t100 * t119, t38, t39 * t90, t58, -0.2e1 * t19 * t100 + 0.2e1 * t32 * t39, 0.2e1 * t20 * t100 + 0.2e1 * t32 * t41, -0.2e1 * t19 * t41 - 0.2e1 * t20 * t39, t19 ^ 2 + t20 ^ 2 + t32 ^ 2, t25 ^ 2, -0.2e1 * t25 * t23, t25 * t120, t23 ^ 2, t23 * t121, t38, 0.2e1 * t18 * t23 + 0.2e1 * t5 * t39, 0.2e1 * t18 * t25 - 0.2e1 * t6 * t39, -0.2e1 * t6 * t23 - 0.2e1 * t5 * t25, t18 ^ 2 + t5 ^ 2 + t6 ^ 2, t12 ^ 2, -0.2e1 * t12 * t10, t12 * t120, t10 ^ 2, t10 * t121, t38, 0.2e1 * t1 * t39 + 0.2e1 * t7 * t10, 0.2e1 * t7 * t12 - 0.2e1 * t2 * t39, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, t100, t72, t42, -t43, 0, 0, t106, t41 * t76 - t95, -t87, -t107, -t86, 0, -pkin(2) * t39 + pkin(8) * t87 - t32 * t76, -pkin(2) * t41 + pkin(8) * t86 + t32 * t74, (t106 - t107) * pkin(8) + t81, -t32 * pkin(2) + t81 * pkin(8), t25 * t99, (-t108 - t109) * t74, -t25 * t76 + t71 * t95, t23 * t103, t23 * t76 - t69 * t95, -t107, t30 * t39 - t5 * t76 + (pkin(8) * t23 + t111) * t74, -t31 * t39 + t6 * t76 + (pkin(8) * t25 + t110) * t74, -t31 * t23 - t30 * t25 + (-t5 * t71 - t6 * t69) * t74, t18 * t63 + t5 * t30 + t6 * t31, t12 * t37, -t37 * t10 - t12 * t35, -t12 * t76 + t37 * t39, t10 * t35, t10 * t76 - t35 * t39, -t107, -t1 * t76 + t49 * t10 + t7 * t35 + t8 * t39, t49 * t12 + t2 * t76 + t7 * t37 - t9 * t39, -t1 * t37 - t9 * t10 - t8 * t12 - t2 * t35, t1 * t8 + t2 * t9 + t7 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t67, t89, 0, t68, 0, 0, 0.2e1 * pkin(2) * t76, -0.2e1 * pkin(2) * t74, 0.2e1 * (t67 + t68) * pkin(8), pkin(2) ^ 2 + t68 * t79 + t62, t66 * t67, -0.2e1 * t67 * t104, -0.2e1 * t71 * t94, t64 * t67, t69 * t89, t68, 0.2e1 * t67 * t113 - 0.2e1 * t30 * t76, 0.2e1 * t67 * pkin(8) * t71 + 0.2e1 * t31 * t76, 0.2e1 * (-t30 * t71 - t31 * t69) * t74, t30 ^ 2 + t31 ^ 2 + t62, t37 ^ 2, -0.2e1 * t37 * t35, t37 * t116, t35 ^ 2, -t35 * t116, t68, 0.2e1 * t49 * t35 - 0.2e1 * t8 * t76, 0.2e1 * t49 * t37 + 0.2e1 * t9 * t76, -0.2e1 * t9 * t35 - 0.2e1 * t8 * t37, t49 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t39, -t100, t19, -t20, 0, 0, t108, -t69 * t23 + t25 * t71, t69 * t39, -t109, t71 * t39, 0, -pkin(3) * t23 - t69 * t91 - t110, -pkin(3) * t25 - t71 * t91 + t111, (t108 - t109) * qJ(4) + t83, -t18 * pkin(3) + qJ(4) * t83, t12 * t48, -t48 * t10 + t12 * t122, t48 * t39, -t10 * t122, t122 * t39, 0, t61 * t10 - t122 * t7 + t27 * t39, t61 * t12 - t28 * t39 + t7 * t48, -t1 * t48 - t28 * t10 - t27 * t12 + t122 * t2, t1 * t27 + t2 * t28 + t7 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, t76, 0, -t63, -t76 * pkin(8), 0, 0, t85, (-t64 + t66) * t74, -t102, -t85, -t98, 0, -pkin(8) * t99 + t82 * t69, pkin(8) * t103 + t71 * t82, t80, -pkin(3) * t63 + qJ(4) * t80, t37 * t48, t122 * t37 - t48 * t35, -t48 * t76, -t35 * t122, -t122 * t76, 0, -t122 * t49 - t27 * t76 + t61 * t35, t28 * t76 + t61 * t37 + t49 * t48, t122 * t9 - t27 * t37 - t28 * t35 - t8 * t48, t8 * t27 + t9 * t28 + t49 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t64, 0.2e1 * t104, 0, t66, 0, 0, 0.2e1 * pkin(3) * t71, -0.2e1 * pkin(3) * t69, 0.2e1 * t92 * qJ(4), qJ(4) ^ 2 * t92 + pkin(3) ^ 2, t48 ^ 2, 0.2e1 * t48 * t122, 0, t122 ^ 2, 0, 0, -t122 * t118, t48 * t118, 0.2e1 * t122 * t28 - 0.2e1 * t27 * t48, t27 ^ 2 + t28 ^ 2 + t61 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t25, 0, t18, 0, 0, 0, 0, 0, 0, t10, t12, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t99, 0, t63, 0, 0, 0, 0, 0, 0, t35, t37, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t69, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t122, t48, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, t39, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t35, -t76, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, t122, 0, t27, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t11;
