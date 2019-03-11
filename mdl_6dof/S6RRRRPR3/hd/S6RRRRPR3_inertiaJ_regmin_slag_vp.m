% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t63 = cos(qJ(3));
t54 = t63 * pkin(2);
t49 = t54 + pkin(3);
t58 = sin(qJ(4));
t62 = cos(qJ(4));
t59 = sin(qJ(3));
t86 = t59 * pkin(2);
t74 = t62 * t86;
t94 = -(pkin(3) + t49) * t58 - t74;
t60 = sin(qJ(2));
t64 = cos(qJ(2));
t34 = t59 * t60 - t63 * t64;
t35 = t59 * t64 + t63 * t60;
t20 = t62 * t34 + t58 * t35;
t21 = -t58 * t34 + t62 * t35;
t50 = -t64 * pkin(2) - pkin(1);
t26 = t34 * pkin(3) + t50;
t68 = -t21 * qJ(5) + t26;
t8 = t20 * pkin(4) + t68;
t93 = -0.2e1 * t8;
t92 = 0.2e1 * t26;
t91 = 0.2e1 * t50;
t90 = 0.2e1 * t64;
t89 = pkin(4) + pkin(10);
t88 = pkin(7) + pkin(8);
t87 = t58 * pkin(3);
t53 = t62 * pkin(3);
t48 = -t53 - pkin(4);
t84 = t21 * t20;
t69 = t58 * t49 + t74;
t28 = qJ(5) + t69;
t83 = t28 * t20;
t45 = qJ(5) + t87;
t82 = t45 * t20;
t57 = sin(qJ(6));
t81 = t57 * t20;
t80 = t57 * t21;
t61 = cos(qJ(6));
t79 = t61 * t20;
t78 = t61 * t57;
t77 = -t62 * t49 + t58 * t86;
t76 = qJ(5) * t20;
t75 = 0.2e1 * t84;
t30 = -pkin(4) + t77;
t67 = -0.2e1 * pkin(4);
t73 = t67 + t77;
t39 = t88 * t60;
t40 = t88 * t64;
t22 = -t63 * t39 - t59 * t40;
t14 = -t35 * pkin(9) + t22;
t23 = t59 * t39 - t63 * t40;
t15 = -t34 * pkin(9) - t23;
t9 = -t62 * t14 + t58 * t15;
t10 = t58 * t14 + t62 * t15;
t27 = -pkin(10) + t30;
t72 = -t21 * t27 + t83;
t43 = -pkin(10) + t48;
t71 = -t21 * t43 + t82;
t70 = t21 * t89 + t76;
t66 = 0.2e1 * qJ(5);
t56 = t61 ^ 2;
t55 = t57 ^ 2;
t52 = qJ(5) * t61;
t51 = qJ(5) * t57;
t44 = -0.2e1 * t78;
t38 = t45 * t61;
t37 = t45 * t57;
t25 = t28 * t61;
t24 = t28 * t57;
t19 = t21 ^ 2;
t18 = t20 ^ 2;
t17 = t61 * t21;
t16 = t20 * t78;
t12 = (-t55 + t56) * t20;
t7 = -t20 * pkin(5) + t10;
t6 = t21 * pkin(5) + t9;
t5 = t20 * t89 + t68;
t4 = t7 * t61;
t3 = t7 * t57;
t2 = t61 * t5 + t57 * t6;
t1 = -t57 * t5 + t61 * t6;
t11 = [1, 0, 0, t60 ^ 2, t60 * t90, 0, 0, 0, pkin(1) * t90, -0.2e1 * pkin(1) * t60, t35 ^ 2, -0.2e1 * t35 * t34, 0, 0, 0, t34 * t91, t35 * t91, t19, -0.2e1 * t84, 0, 0, 0, t20 * t92, t21 * t92, -0.2e1 * t10 * t20 + 0.2e1 * t9 * t21, t20 * t93, t21 * t93, t10 ^ 2 + t8 ^ 2 + t9 ^ 2, t55 * t18, 0.2e1 * t18 * t78, t57 * t75, t61 * t75, t19, 0.2e1 * t1 * t21 - 0.2e1 * t7 * t79, -0.2e1 * t2 * t21 + 0.2e1 * t7 * t81; 0, 0, 0, 0, 0, t60, t64, 0, -t60 * pkin(7), -t64 * pkin(7), 0, 0, t35, -t34, 0, t22, t23, 0, 0, t21, -t20, 0, -t9, -t10, t30 * t21 - t83, t9, t10, t10 * t28 + t9 * t30, t16, t12, t17, -t80, 0, -t61 * t72 + t3, t57 * t72 + t4; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t54, -0.2e1 * t86, 0, 0, 0, 0, 1, -0.2e1 * t77, -0.2e1 * t69, 0, 0.2e1 * t30, 0.2e1 * t28, t28 ^ 2 + t30 ^ 2, t56, t44, 0, 0, 0, 0.2e1 * t24, 0.2e1 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, 0, t22, t23, 0, 0, t21, -t20, 0, -t9, -t10, t48 * t21 - t82, t9, t10, t10 * t45 + t9 * t48, t16, t12, t17, -t80, 0, -t61 * t71 + t3, t57 * t71 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t54, -t86, 0, 0, 0, 0, 1, t53 - t77, t94, 0, -t53 + t73, t66 - t94, t28 * t45 + t30 * t48, t56, t44, 0, 0, 0, t37 + t24, t38 + t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t53, -0.2e1 * t87, 0, 0.2e1 * t48, 0.2e1 * t45, t45 ^ 2 + t48 ^ 2, t56, t44, 0, 0, 0, 0.2e1 * t37, 0.2e1 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t9, -t10, -pkin(4) * t21 - t76, t9, t10, -t9 * pkin(4) + t10 * qJ(5), t16, t12, t17, -t80, 0, -t61 * t70 + t3, t57 * t70 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t77, -t69, 0, t73, t66 + t69, -t30 * pkin(4) + t28 * qJ(5), t56, t44, 0, 0, 0, t51 + t24, t52 + t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t53, -t87, 0, t67 - t53, t66 + t87, -t48 * pkin(4) + t45 * qJ(5), t56, t44, 0, 0, 0, t51 + t37, t52 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t67, t66, pkin(4) ^ 2 + qJ(5) ^ 2, t56, t44, 0, 0, 0, 0.2e1 * t51, 0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, t9, 0, 0, 0, 0, 0, t17, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t48, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t79, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t57, 0, t61 * t27, -t57 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t57, 0, t61 * t43, -t57 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t57, 0, -t61 * t89, t57 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;
