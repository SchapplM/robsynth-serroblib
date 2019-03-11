% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t55 = sin(qJ(4));
t48 = t55 ^ 2;
t57 = cos(qJ(4));
t50 = t57 ^ 2;
t78 = t48 + t50;
t43 = qJ(5) * t57;
t101 = -pkin(4) * t55 + t43;
t53 = pkin(4) + qJ(6);
t77 = t55 * qJ(5);
t100 = -t53 * t57 - t77;
t99 = 0.2e1 * t53;
t98 = -0.2e1 * t55;
t97 = 0.2e1 * t57;
t58 = cos(qJ(3));
t96 = 0.2e1 * t58;
t95 = 2 * qJ(2);
t56 = sin(qJ(3));
t93 = pkin(8) * t56;
t92 = t56 * pkin(4);
t91 = t58 * pkin(3);
t51 = t58 ^ 2;
t59 = -pkin(1) - pkin(7);
t90 = t51 * t59;
t36 = t55 * t56;
t88 = t55 * t57;
t37 = t55 * t58;
t87 = t56 * t59;
t40 = t57 * t58;
t86 = t57 * t59;
t10 = -pkin(3) + t100;
t85 = t58 * t10;
t67 = -t57 * pkin(4) - t77;
t20 = -pkin(3) + t67;
t84 = t58 * t20;
t83 = t58 * t56;
t82 = t58 * t59;
t19 = t56 * pkin(3) - t58 * pkin(8) + qJ(2);
t81 = -t57 * t19 + t55 * t87;
t7 = t55 * t19 + t56 * t86;
t80 = t78 * t93;
t79 = t78 * pkin(8) ^ 2;
t49 = t56 ^ 2;
t29 = t49 + t51;
t76 = t56 * qJ(5);
t75 = t55 * t83;
t74 = t51 * t88;
t73 = t57 * t83;
t72 = -t59 - t43;
t71 = -t91 - t93;
t4 = -t76 - t7;
t5 = t81 - t92;
t70 = -t4 * t57 + t5 * t55;
t69 = t55 * t81 + t7 * t57;
t68 = -t84 + t93;
t45 = t55 * pkin(8);
t24 = t55 * pkin(5) + t45;
t47 = t57 * pkin(8);
t25 = t57 * pkin(5) + t47;
t66 = t24 * t55 + t25 * t57;
t65 = -pkin(5) * t37 + t7;
t64 = -pkin(5) * t40 - t81;
t62 = qJ(2) ^ 2;
t61 = qJ(5) ^ 2;
t60 = 0.2e1 * qJ(5);
t52 = t59 ^ 2;
t42 = 0.2e1 * t76;
t41 = t51 * t52;
t39 = t57 * t56;
t38 = t50 * t51;
t35 = t48 * t51;
t34 = 0.2e1 * t88;
t32 = pkin(4) * t37;
t30 = t57 * t76;
t26 = t55 * t40;
t23 = 0.2e1 * t73;
t22 = -0.2e1 * t74;
t21 = 0.2e1 * t75;
t18 = 0.2e1 * t78 * pkin(8);
t17 = t29 * t59;
t16 = t29 * t57;
t15 = t29 * t55;
t14 = (-t48 + t50) * t58;
t13 = t78 * t56;
t9 = t78 * t49 + t51;
t8 = t72 * t58 + t32;
t3 = t32 + (qJ(6) * t55 + t72) * t58;
t2 = t65 + t76;
t1 = -t53 * t56 - t64;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t95, pkin(1) ^ 2 + t62, t51, -0.2e1 * t83, 0, t49, 0, 0, t56 * t95, t58 * t95, -0.2e1 * t17, t49 * t52 + t41 + t62, t38, t22, t23, t35, -0.2e1 * t75, t49, -0.2e1 * t55 * t90 - 0.2e1 * t56 * t81, -0.2e1 * t51 * t86 - 0.2e1 * t7 * t56 (-t55 * t7 + t57 * t81) * t96, t7 ^ 2 + t81 ^ 2 + t41, t49, -0.2e1 * t73, t21, t38, t22, t35 (t4 * t55 + t5 * t57) * t96, -0.2e1 * t8 * t37 + 0.2e1 * t5 * t56, -0.2e1 * t4 * t56 - 0.2e1 * t8 * t40, t4 ^ 2 + t5 ^ 2 + t8 ^ 2, t49, t21, t23, t35, 0.2e1 * t74, t38 (t1 * t57 - t2 * t55) * t96, 0.2e1 * t2 * t56 - 0.2e1 * t3 * t40, -0.2e1 * t1 * t56 + 0.2e1 * t3 * t37, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t29, t17, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, t69 * t56 + t90, 0, 0, 0, 0, 0, 0, 0, t15, t16, t70 * t56 - t8 * t58, 0, 0, 0, 0, 0, 0, 0, t16, -t15, -t3 * t58 + (t1 * t55 + t2 * t57) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t56, 0, t82, -t87, 0, 0, t26, t14, t36, -t26, t39, 0, t71 * t55 + t57 * t82, -t55 * t82 + t71 * t57, t69, pkin(3) * t82 + t69 * pkin(8), 0, -t36, -t39, t26, t14, -t26, t70, t68 * t55 + t8 * t57, -t8 * t55 + t68 * t57, t70 * pkin(8) + t8 * t20, 0, -t39, t36, -t26, -t14, t26 (t24 * t58 + t2) * t57 + (-t25 * t58 + t1) * t55, t25 * t56 - t3 * t55 - t57 * t85, -t24 * t56 - t3 * t57 + t55 * t85, t1 * t24 + t3 * t10 + t2 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t56, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t37, t13, t80 + t91, 0, 0, 0, 0, 0, 0, t13, -t40, t37, t80 - t84, 0, 0, 0, 0, 0, 0, t13, t37, t40, t66 * t56 - t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t48, t34, 0, t50, 0, 0, pkin(3) * t97, pkin(3) * t98, t18, pkin(3) ^ 2 + t79, 0, 0, 0, t48, t34, t50, t18, t20 * t97, t20 * t98, t20 ^ 2 + t79, 0, 0, 0, t50, -0.2e1 * t88, t48, 0.2e1 * t66, t10 * t98, -0.2e1 * t10 * t57, t10 ^ 2 + t24 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t37, t56, -t81, -t7, 0, 0, t56, -t40, t37, 0, 0, 0, t67 * t58, t81 - 0.2e1 * t92, t42 + t7, -t5 * pkin(4) - t4 * qJ(5), t56, t37, t40, 0, 0, 0, t100 * t58, t42 + t65, t56 * t99 + t64, t2 * qJ(5) - t1 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t39, -pkin(4) * t36 + t30, 0, 0, 0, 0, 0, 0, 0, t39, -t36, -t53 * t36 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t57, 0, -t45, -t47, 0, 0, 0, -t55, -t57, 0, 0, 0, t101, t45, t47, t101 * pkin(8), 0, -t57, t55, 0, 0, 0, -t53 * t55 + t43, t25, -t24, t25 * qJ(5) - t24 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(4), t60, pkin(4) ^ 2 + t61, 1, 0, 0, 0, 0, 0, 0, t60, t99, t53 ^ 2 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t56, 0, t5, 0, 0, 0, 0, 0, 0, t40, 0, -t56, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, t45, 0, 0, 0, 0, 0, 0, t55, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0, 0, -1, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t56, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
