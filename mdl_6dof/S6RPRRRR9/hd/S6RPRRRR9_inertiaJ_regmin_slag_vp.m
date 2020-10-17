% Calculate minimal parameter regressor of joint inertia matrix for
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
% MM_reg [((6+1)*6/2)x34]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:28:52
% EndTime: 2019-05-06 04:28:55
% DurationCPUTime: 0.92s
% Computational Cost: add. (734->125), mult. (1512->202), div. (0->0), fcn. (1775->8), ass. (0->82)
t57 = sin(qJ(5));
t58 = sin(qJ(4));
t61 = cos(qJ(5));
t62 = cos(qJ(4));
t37 = t57 * t58 - t61 * t62;
t49 = -t62 * pkin(4) - pkin(3);
t27 = t37 * pkin(5) + t49;
t92 = 0.2e1 * t27;
t91 = 0.2e1 * t49;
t59 = sin(qJ(3));
t90 = -0.2e1 * t59;
t89 = 0.2e1 * t59;
t88 = 0.2e1 * t62;
t87 = 2 * qJ(2);
t86 = pkin(8) + pkin(9);
t56 = sin(qJ(6));
t85 = t56 * pkin(5);
t84 = t57 * pkin(4);
t83 = t59 * pkin(5);
t60 = cos(qJ(6));
t50 = t60 * pkin(5);
t38 = t57 * t62 + t61 * t58;
t63 = cos(qJ(3));
t29 = t63 * t38;
t40 = t59 * pkin(3) - t63 * pkin(8) + qJ(2);
t35 = t62 * t40;
t47 = t62 * t63;
t64 = -pkin(1) - pkin(7);
t77 = t58 * t64;
t17 = -pkin(9) * t47 + t35 + (pkin(4) - t77) * t59;
t73 = t62 * t64;
t66 = t59 * t73;
t21 = t66 + (-pkin(9) * t63 + t40) * t58;
t75 = t61 * t21;
t9 = t57 * t17 + t75;
t5 = -t29 * pkin(10) + t9;
t82 = t60 * t5;
t51 = t61 * pkin(4);
t81 = t38 * t59;
t80 = t58 * t59;
t79 = t58 * t62;
t78 = t58 * t63;
t76 = t59 * t64;
t74 = t62 * t59;
t72 = t63 * t37;
t71 = t63 * t59;
t70 = t63 * t64;
t53 = t59 ^ 2;
t55 = t63 ^ 2;
t69 = -t53 - t55;
t68 = -0.2e1 * t71;
t67 = t60 * t84;
t8 = t61 * t17 - t57 * t21;
t4 = pkin(10) * t72 + t8 + t83;
t1 = t60 * t4 - t56 * t5;
t41 = t86 * t58;
t42 = t86 * t62;
t22 = -t61 * t41 - t57 * t42;
t36 = pkin(4) * t78 - t70;
t48 = t51 + pkin(5);
t32 = t60 * t48 - t56 * t84;
t65 = -pkin(3) * t63 - pkin(8) * t59;
t2 = t56 * t4 + t82;
t23 = -t57 * t41 + t61 * t42;
t54 = t62 ^ 2;
t52 = t58 ^ 2;
t33 = t56 * t48 + t67;
t30 = -t57 * t80 + t61 * t74;
t25 = t58 * t40 + t66;
t24 = -t58 * t76 + t35;
t20 = -t56 * t37 + t60 * t38;
t19 = t60 * t37 + t56 * t38;
t18 = t29 * pkin(5) + t36;
t15 = -t37 * pkin(10) + t23;
t14 = -t38 * pkin(10) + t22;
t13 = -t56 * t29 - t60 * t72;
t12 = t60 * t30 - t56 * t81;
t11 = t60 * t29 - t56 * t72;
t10 = -t56 * t30 - t60 * t81;
t7 = t56 * t14 + t60 * t15;
t6 = t60 * t14 - t56 * t15;
t3 = [1, 0, 0, -2 * pkin(1), t87, pkin(1) ^ 2 + qJ(2) ^ 2, t55, t68, 0, 0, 0, t59 * t87, t63 * t87, t54 * t55, -0.2e1 * t55 * t79, t71 * t88, t58 * t68, t53, 0.2e1 * t24 * t59 - 0.2e1 * t55 * t77, -0.2e1 * t25 * t59 - 0.2e1 * t55 * t73, t72 ^ 2, 0.2e1 * t72 * t29, -t72 * t89, t29 * t90, t53, 0.2e1 * t36 * t29 + 0.2e1 * t8 * t59, -0.2e1 * t36 * t72 - 0.2e1 * t9 * t59, t13 ^ 2, -0.2e1 * t13 * t11, t13 * t89, t11 * t90, t53, 0.2e1 * t1 * t59 + 0.2e1 * t18 * t11, 0.2e1 * t18 * t13 - 0.2e1 * t2 * t59; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t58, t69 * t62, 0, 0, 0, 0, 0, -t63 * t29 - t59 * t81, -t30 * t59 + t63 * t72, 0, 0, 0, 0, 0, t10 * t59 - t63 * t11, -t12 * t59 - t63 * t13; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t63, -t59, 0, t70, -t76, t58 * t47 (-t52 + t54) * t63, t80, t74, 0, t65 * t58 + t62 * t70, -t58 * t70 + t65 * t62, -t72 * t38, -t38 * t29 + t37 * t72, t81, -t37 * t59, 0, t22 * t59 + t49 * t29 + t36 * t37, -t23 * t59 + t36 * t38 - t49 * t72, t13 * t20, -t20 * t11 - t13 * t19, t20 * t59, -t19 * t59, 0, t27 * t11 + t18 * t19 + t6 * t59, t27 * t13 + t18 * t20 - t7 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t59, 0, 0, 0, 0, 0, t47, -t78, 0, 0, 0, 0, 0, -t72, -t29, 0, 0, 0, 0, 0, -t63 * t19, -t63 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t52, 0.2e1 * t79, 0, 0, 0, pkin(3) * t88, -0.2e1 * pkin(3) * t58, t38 ^ 2, -0.2e1 * t38 * t37, 0, 0, 0, t37 * t91, t38 * t91, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t92, t20 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t78, t59, t24, -t25, 0, 0, -t72, -t29, t59, t59 * t51 + t8, -t75 + (-t59 * pkin(4) - t17) * t57, 0, 0, t13, -t11, t59, t32 * t59 + t1, -t33 * t59 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t74, 0, 0, 0, 0, 0, -t81, -t30, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t62, 0, -t58 * pkin(8), -t62 * pkin(8), 0, 0, t38, -t37, 0, t22, -t23, 0, 0, t20, -t19, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t84, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t29, t59, t8, -t9, 0, 0, t13, -t11, t59, t59 * t50 + t1, -t82 + (-t4 - t83) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t30, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, 0, t22, -t23, 0, 0, t20, -t19, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t51, -t84, 0, 0, 0, 0, 1, t32 + t50, -t67 + (-pkin(5) - t48) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t11, t59, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
