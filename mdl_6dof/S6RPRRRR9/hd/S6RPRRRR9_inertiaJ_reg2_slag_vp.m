% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t71 = cos(qJ(4));
t72 = cos(qJ(3));
t54 = t71 * t72;
t66 = sin(qJ(5));
t70 = cos(qJ(5));
t67 = sin(qJ(4));
t94 = t67 * t72;
t34 = t70 * t54 - t66 * t94;
t65 = sin(qJ(6));
t69 = cos(qJ(6));
t42 = t66 * t71 + t67 * t70;
t88 = t72 * t42;
t11 = t34 * t65 + t69 * t88;
t109 = -0.2e1 * t11;
t40 = t66 * t67 - t70 * t71;
t56 = -pkin(4) * t71 - pkin(3);
t29 = pkin(5) * t40 + t56;
t108 = 0.2e1 * t29;
t107 = 0.2e1 * t56;
t68 = sin(qJ(3));
t106 = 0.2e1 * t68;
t105 = 0.2e1 * t71;
t104 = 2 * qJ(2);
t103 = -pkin(9) - pkin(8);
t102 = pkin(5) * t68;
t101 = t65 * pkin(5);
t100 = t66 * pkin(4);
t58 = t69 * pkin(5);
t45 = pkin(3) * t68 - pkin(8) * t72 + qJ(2);
t38 = t71 * t45;
t73 = -pkin(1) - pkin(7);
t93 = t67 * t73;
t18 = -pkin(9) * t54 + t38 + (pkin(4) - t93) * t68;
t89 = t71 * t73;
t82 = t68 * t89;
t23 = t82 + (-pkin(9) * t72 + t45) * t67;
t91 = t70 * t23;
t9 = t18 * t66 + t91;
t5 = -pkin(10) * t88 + t9;
t99 = t69 * t5;
t59 = t70 * pkin(4);
t98 = t72 * pkin(3);
t97 = t42 * t68;
t96 = t67 * t68;
t95 = t67 * t71;
t92 = t68 * t73;
t90 = t71 * t68;
t87 = t72 * t68;
t86 = t72 * t73;
t60 = t67 ^ 2;
t62 = t71 ^ 2;
t85 = t60 + t62;
t61 = t68 ^ 2;
t63 = t72 ^ 2;
t51 = t61 + t63;
t84 = -0.2e1 * t87;
t83 = t69 * t100;
t81 = t67 * t54;
t8 = t70 * t18 - t23 * t66;
t4 = -pkin(10) * t34 + t102 + t8;
t1 = t69 * t4 - t5 * t65;
t80 = t85 * t68;
t46 = t103 * t67;
t47 = t103 * t71;
t24 = t70 * t46 + t47 * t66;
t39 = pkin(4) * t94 - t86;
t55 = t59 + pkin(5);
t35 = -t65 * t100 + t69 * t55;
t79 = -pkin(8) * t68 - t98;
t2 = t4 * t65 + t99;
t26 = -t67 * t92 + t38;
t27 = t45 * t67 + t82;
t78 = -t26 * t67 + t27 * t71;
t25 = t46 * t66 - t47 * t70;
t74 = qJ(2) ^ 2;
t64 = t73 ^ 2;
t57 = t63 * t64;
t44 = t51 * t73;
t36 = t55 * t65 + t83;
t33 = -t66 * t96 + t70 * t90;
t22 = -t40 * t65 + t42 * t69;
t20 = t69 * t40 + t42 * t65;
t19 = pkin(5) * t88 + t39;
t16 = -pkin(10) * t40 + t25;
t15 = -pkin(10) * t42 + t24;
t14 = t34 * t69 - t65 * t88;
t13 = t33 * t69 - t65 * t97;
t10 = -t33 * t65 - t69 * t97;
t7 = t15 * t65 + t16 * t69;
t6 = t15 * t69 - t16 * t65;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t104, pkin(1) ^ 2 + t74, t63, t84, 0, t61, 0, 0, t68 * t104, t72 * t104, -0.2e1 * t44, t61 * t64 + t57 + t74, t62 * t63, -0.2e1 * t63 * t95, t87 * t105, t60 * t63, t67 * t84, t61, 0.2e1 * t26 * t68 - 0.2e1 * t63 * t93, -0.2e1 * t27 * t68 - 0.2e1 * t63 * t89, 0.2e1 * (-t26 * t71 - t27 * t67) * t72, t26 ^ 2 + t27 ^ 2 + t57, t34 ^ 2, -0.2e1 * t34 * t88, t34 * t106, t88 ^ 2, -t88 * t106, t61, 0.2e1 * t39 * t88 + 0.2e1 * t68 * t8, 0.2e1 * t34 * t39 - 0.2e1 * t68 * t9, -0.2e1 * t34 * t8 - 0.2e1 * t88 * t9, t39 ^ 2 + t8 ^ 2 + t9 ^ 2, t14 ^ 2, t14 * t109, t14 * t106, t11 ^ 2, t68 * t109, t61, 0.2e1 * t1 * t68 + 0.2e1 * t11 * t19, 0.2e1 * t14 * t19 - 0.2e1 * t2 * t68, -0.2e1 * t1 * t14 - 0.2e1 * t11 * t2, t1 ^ 2 + t19 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t51, t44, 0, 0, 0, 0, 0, 0, -t51 * t67, -t51 * t71, 0, t63 * t73 + t68 * t78, 0, 0, 0, 0, 0, 0, -t68 * t97 - t72 * t88, -t33 * t68 - t34 * t72, -t33 * t88 + t34 * t97, t33 * t9 - t39 * t72 - t8 * t97, 0, 0, 0, 0, 0, 0, t10 * t68 - t11 * t72, -t13 * t68 - t14 * t72, -t10 * t14 - t11 * t13, t1 * t10 + t13 * t2 - t19 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t85 + t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 ^ 2 + t97 ^ 2 + t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t13 ^ 2 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, -t68, 0, t86, -t92, 0, 0, t81 (-t60 + t62) * t72, t96, -t81, t90, 0, t79 * t67 + t71 * t86, -t67 * t86 + t79 * t71, t78, pkin(3) * t86 + pkin(8) * t78, t34 * t42, -t34 * t40 - t42 * t88, t97, t88 * t40, -t40 * t68, 0, t24 * t68 + t39 * t40 + t56 * t88, -t25 * t68 + t34 * t56 + t39 * t42, -t24 * t34 - t25 * t88 - t40 * t9 - t42 * t8, t24 * t8 + t25 * t9 + t39 * t56, t14 * t22, -t11 * t22 - t14 * t20, t22 * t68, t11 * t20, -t20 * t68, 0, t11 * t29 + t19 * t20 + t6 * t68, t14 * t29 + t19 * t22 - t68 * t7, -t1 * t22 - t11 * t7 - t14 * t6 - t2 * t20, t1 * t6 + t19 * t29 + t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t68, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t94, t80, pkin(8) * t80 + t98, 0, 0, 0, 0, 0, 0, -t72 * t40, -t88, -t33 * t40 + t42 * t97, -t24 * t97 + t25 * t33 - t56 * t72, 0, 0, 0, 0, 0, 0, -t72 * t20, -t72 * t22, -t10 * t22 - t13 * t20, t10 * t6 + t13 * t7 - t29 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t60, 0.2e1 * t95, 0, t62, 0, 0, pkin(3) * t105, -0.2e1 * pkin(3) * t67, 0.2e1 * t85 * pkin(8), pkin(8) ^ 2 * t85 + pkin(3) ^ 2, t42 ^ 2, -0.2e1 * t42 * t40, 0, t40 ^ 2, 0, 0, t40 * t107, t42 * t107, -0.2e1 * t24 * t42 - 0.2e1 * t25 * t40, t24 ^ 2 + t25 ^ 2 + t56 ^ 2, t22 ^ 2, -0.2e1 * t22 * t20, 0, t20 ^ 2, 0, 0, t20 * t108, t22 * t108, -0.2e1 * t20 * t7 - 0.2e1 * t22 * t6, t29 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t94, t68, t26, -t27, 0, 0, 0, 0, t34, 0, -t88, t68, t68 * t59 + t8, -t91 + (-pkin(4) * t68 - t18) * t66 (-t34 * t70 - t66 * t88) * pkin(4) (t66 * t9 + t70 * t8) * pkin(4), 0, 0, t14, 0, -t11, t68, t35 * t68 + t1, -t36 * t68 - t2, -t11 * t36 - t14 * t35, t1 * t35 + t2 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t90, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t33, 0 (t33 * t66 - t70 * t97) * pkin(4), 0, 0, 0, 0, 0, 0, t10, -t13, 0, t10 * t35 + t13 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, t71, 0, -t67 * pkin(8), -t71 * pkin(8), 0, 0, 0, 0, t42, 0, -t40, 0, t24, -t25 (-t40 * t66 - t42 * t70) * pkin(4) (t24 * t70 + t25 * t66) * pkin(4), 0, 0, t22, 0, -t20, 0, t6, -t7, -t20 * t36 - t22 * t35, t35 * t6 + t36 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t100, 0 (t66 ^ 2 + t70 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t36, 0, t35 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t88, t68, t8, -t9, 0, 0, 0, 0, t14, 0, -t11, t68, t68 * t58 + t1, -t99 + (-t4 - t102) * t65 (-t11 * t65 - t14 * t69) * pkin(5) (t1 * t69 + t2 * t65) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t33, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t13, 0 (t10 * t69 + t13 * t65) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, -t40, 0, t24, -t25, 0, 0, 0, 0, t22, 0, -t20, 0, t6, -t7 (-t20 * t65 - t22 * t69) * pkin(5) (t6 * t69 + t65 * t7) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t100, 0, 0, 0, 0, 0, 0, 0, 1, t35 + t58, -t83 + (-pkin(5) - t55) * t65, 0 (t35 * t69 + t36 * t65) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t101, 0 (t65 ^ 2 + t69 ^ 2) * pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t11, t68, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t20, 0, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t101, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t3;
