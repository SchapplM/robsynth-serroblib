% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:42
% EndTime: 2019-12-31 22:43:48
% DurationCPUTime: 1.62s
% Computational Cost: add. (2074->209), mult. (4885->421), div. (0->0), fcn. (5437->10), ass. (0->109)
t66 = sin(pkin(5));
t71 = sin(qJ(2));
t100 = t66 * t71;
t67 = cos(pkin(5));
t70 = sin(qJ(3));
t74 = cos(qJ(3));
t39 = t70 * t100 - t67 * t74;
t38 = t39 ^ 2;
t125 = -0.2e1 * t39;
t124 = 0.2e1 * t39;
t41 = t74 * t100 + t67 * t70;
t123 = -0.2e1 * t41;
t73 = cos(qJ(4));
t59 = -pkin(4) * t73 - pkin(3);
t122 = 0.2e1 * t59;
t121 = 0.2e1 * t66;
t120 = -0.2e1 * t70;
t119 = -0.2e1 * t74;
t118 = 0.2e1 * t74;
t117 = -pkin(10) - pkin(9);
t116 = pkin(1) * t71;
t75 = cos(qJ(2));
t115 = pkin(1) * t75;
t114 = pkin(3) * t73;
t113 = pkin(4) * t39;
t69 = sin(qJ(4));
t112 = pkin(8) * t69;
t68 = sin(qJ(5));
t111 = t68 * pkin(4);
t110 = t70 * pkin(8);
t72 = cos(qJ(5));
t109 = t72 * pkin(4);
t99 = t66 * t75;
t23 = t41 * t69 + t73 * t99;
t54 = pkin(7) * t100;
t32 = t54 + (-pkin(2) - t115) * t67;
t15 = pkin(3) * t39 - pkin(9) * t41 + t32;
t85 = pkin(7) * t99;
t33 = t85 + (pkin(8) + t116) * t67;
t34 = (-pkin(2) * t75 - pkin(8) * t71 - pkin(1)) * t66;
t19 = t33 * t74 + t34 * t70;
t17 = -pkin(9) * t99 + t19;
t7 = t15 * t69 + t17 * t73;
t5 = -pkin(10) * t23 + t7;
t108 = t72 * t5;
t18 = -t70 * t33 + t34 * t74;
t16 = pkin(3) * t99 - t18;
t107 = t16 * t69;
t106 = t16 * t73;
t105 = t23 * t73;
t25 = t41 * t73 - t69 * t99;
t104 = t25 * t69;
t103 = t39 * t74;
t102 = t41 * t70;
t61 = t66 ^ 2;
t101 = t61 * t75;
t98 = t67 * t71;
t97 = t69 * t39;
t96 = t69 * t70;
t95 = t69 * t73;
t94 = t69 * t74;
t49 = -pkin(3) * t74 - pkin(9) * t70 - pkin(2);
t90 = t73 * t74;
t86 = pkin(8) * t90;
t26 = t86 + (-pkin(10) * t70 + t49) * t69;
t93 = t72 * t26;
t92 = t73 * t39;
t91 = t73 * t70;
t62 = t69 ^ 2;
t64 = t73 ^ 2;
t89 = t62 + t64;
t88 = 0.2e1 * t99;
t87 = t70 * t118;
t84 = t70 * t99;
t83 = t74 * t99;
t82 = t69 * t91;
t6 = t73 * t15 - t17 * t69;
t4 = -pkin(10) * t25 + t113 + t6;
t1 = t72 * t4 - t5 * t68;
t44 = t73 * t49;
t22 = -pkin(10) * t91 + t44 + (-pkin(4) - t112) * t74;
t9 = t72 * t22 - t26 * t68;
t81 = -t6 * t69 + t7 * t73;
t80 = -t18 * t70 + t19 * t74;
t30 = -pkin(8) * t94 + t44;
t31 = t49 * t69 + t86;
t79 = -t30 * t69 + t31 * t73;
t47 = t68 * t73 + t69 * t72;
t77 = pkin(8) ^ 2;
t65 = t74 ^ 2;
t63 = t70 ^ 2;
t60 = t63 * t77;
t56 = t61 * t75 ^ 2;
t51 = t117 * t73;
t50 = t117 * t69;
t48 = (pkin(4) * t69 + pkin(8)) * t70;
t45 = t68 * t69 - t72 * t73;
t43 = pkin(1) * t98 + t85;
t42 = t67 * t115 - t54;
t37 = -t68 * t96 + t72 * t91;
t35 = t47 * t70;
t28 = t50 * t68 - t51 * t72;
t27 = t50 * t72 + t51 * t68;
t13 = -t23 * t68 + t25 * t72;
t11 = t72 * t23 + t25 * t68;
t10 = t22 * t68 + t93;
t8 = pkin(4) * t23 + t16;
t2 = t4 * t68 + t108;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t61 * t71 ^ 2, 0.2e1 * t71 * t101, t98 * t121, t56, t67 * t88, t67 ^ 2, 0.2e1 * pkin(1) * t101 + 0.2e1 * t42 * t67, -0.2e1 * t61 * t116 - 0.2e1 * t43 * t67, (-t42 * t71 + t43 * t75) * t121, pkin(1) ^ 2 * t61 + t42 ^ 2 + t43 ^ 2, t41 ^ 2, t39 * t123, t99 * t123, t38, t39 * t88, t56, -0.2e1 * t18 * t99 + 0.2e1 * t32 * t39, 0.2e1 * t19 * t99 + 0.2e1 * t32 * t41, -0.2e1 * t18 * t41 - 0.2e1 * t19 * t39, t18 ^ 2 + t19 ^ 2 + t32 ^ 2, t25 ^ 2, -0.2e1 * t25 * t23, t25 * t124, t23 ^ 2, t23 * t125, t38, 0.2e1 * t16 * t23 + 0.2e1 * t39 * t6, 0.2e1 * t16 * t25 - 0.2e1 * t39 * t7, -0.2e1 * t23 * t7 - 0.2e1 * t25 * t6, t16 ^ 2 + t6 ^ 2 + t7 ^ 2, t13 ^ 2, -0.2e1 * t13 * t11, t13 * t124, t11 ^ 2, t11 * t125, t38, 0.2e1 * t1 * t39 + 0.2e1 * t11 * t8, 0.2e1 * t13 * t8 - 0.2e1 * t2 * t39, -0.2e1 * t1 * t13 - 0.2e1 * t11 * t2, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, t99, t67, t42, -t43, 0, 0, t102, -t39 * t70 + t41 * t74, -t84, -t103, -t83, 0, -pkin(2) * t39 + pkin(8) * t84 - t32 * t74, -pkin(2) * t41 + pkin(8) * t83 + t32 * t70, (t102 - t103) * pkin(8) + t80, -t32 * pkin(2) + pkin(8) * t80, t25 * t91, (-t104 - t105) * t70, -t25 * t74 + t39 * t91, t23 * t96, t23 * t74 - t39 * t96, -t103, t30 * t39 - t6 * t74 + (pkin(8) * t23 + t107) * t70, -t31 * t39 + t7 * t74 + (pkin(8) * t25 + t106) * t70, -t23 * t31 - t25 * t30 + (-t6 * t73 - t69 * t7) * t70, t110 * t16 + t30 * t6 + t31 * t7, t13 * t37, -t11 * t37 - t13 * t35, -t13 * t74 + t37 * t39, t11 * t35, t11 * t74 - t35 * t39, -t103, -t1 * t74 + t11 * t48 + t35 * t8 + t39 * t9, -t10 * t39 + t13 * t48 + t2 * t74 + t37 * t8, -t1 * t37 - t10 * t11 - t13 * t9 - t2 * t35, t1 * t9 + t10 * t2 + t48 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t63, t87, 0, t65, 0, 0, pkin(2) * t118, pkin(2) * t120, 0.2e1 * (t63 + t65) * pkin(8), pkin(2) ^ 2 + t65 * t77 + t60, t64 * t63, -0.2e1 * t63 * t95, t90 * t120, t62 * t63, t69 * t87, t65, 0.2e1 * t112 * t63 - 0.2e1 * t30 * t74, 0.2e1 * pkin(8) * t63 * t73 + 0.2e1 * t31 * t74, 0.2e1 * (-t30 * t73 - t31 * t69) * t70, t30 ^ 2 + t31 ^ 2 + t60, t37 ^ 2, -0.2e1 * t37 * t35, t37 * t119, t35 ^ 2, -t35 * t119, t65, 0.2e1 * t35 * t48 - 0.2e1 * t74 * t9, 0.2e1 * t10 * t74 + 0.2e1 * t37 * t48, -0.2e1 * t10 * t35 - 0.2e1 * t37 * t9, t10 ^ 2 + t48 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t39, -t99, t18, -t19, 0, 0, t104, -t23 * t69 + t25 * t73, t97, -t105, t92, 0, -pkin(3) * t23 - pkin(9) * t97 - t106, -pkin(3) * t25 - pkin(9) * t92 + t107, (t104 - t105) * pkin(9) + t81, -t16 * pkin(3) + pkin(9) * t81, t13 * t47, -t11 * t47 - t13 * t45, t47 * t39, t11 * t45, -t45 * t39, 0, t11 * t59 + t27 * t39 + t45 * t8, t13 * t59 - t28 * t39 + t47 * t8, -t1 * t47 - t11 * t28 - t13 * t27 - t2 * t45, t1 * t27 + t2 * t28 + t59 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, t74, 0, -t110, -t74 * pkin(8), 0, 0, t82, (-t62 + t64) * t70, -t94, -t82, -t90, 0, -pkin(8) * t91 + (-pkin(3) * t70 + pkin(9) * t74) * t69, pkin(9) * t90 + (t112 - t114) * t70, t79, -pkin(3) * t110 + pkin(9) * t79, t37 * t47, -t35 * t47 - t37 * t45, -t47 * t74, t35 * t45, t45 * t74, 0, -t27 * t74 + t35 * t59 + t45 * t48, t28 * t74 + t37 * t59 + t47 * t48, -t10 * t45 - t27 * t37 - t28 * t35 - t47 * t9, t10 * t28 + t27 * t9 + t48 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t62, 0.2e1 * t95, 0, t64, 0, 0, 0.2e1 * t114, -0.2e1 * pkin(3) * t69, 0.2e1 * t89 * pkin(9), pkin(9) ^ 2 * t89 + pkin(3) ^ 2, t47 ^ 2, -0.2e1 * t47 * t45, 0, t45 ^ 2, 0, 0, t45 * t122, t47 * t122, -0.2e1 * t27 * t47 - 0.2e1 * t28 * t45, t27 ^ 2 + t28 ^ 2 + t59 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, t39, t6, -t7, 0, 0, 0, 0, t13, 0, -t11, t39, t109 * t39 + t1, -t108 + (-t4 - t113) * t68, (-t11 * t68 - t13 * t72) * pkin(4), (t1 * t72 + t2 * t68) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, -t96, -t74, t30, -t31, 0, 0, 0, 0, t37, 0, -t35, -t74, -t109 * t74 + t9, -t93 + (pkin(4) * t74 - t22) * t68, (-t35 * t68 - t37 * t72) * pkin(4), (t10 * t68 + t72 * t9) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t73, 0, -t69 * pkin(9), -t73 * pkin(9), 0, 0, 0, 0, t47, 0, -t45, 0, t27, -t28, (-t45 * t68 - t47 * t72) * pkin(4), (t27 * t72 + t28 * t68) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t109, -0.2e1 * t111, 0, (t68 ^ 2 + t72 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t11, t39, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t35, -t74, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t45, 0, t27, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t109, -t111, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
