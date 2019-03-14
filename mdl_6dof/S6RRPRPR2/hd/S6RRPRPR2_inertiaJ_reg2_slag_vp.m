% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t55 = sin(pkin(10));
t56 = cos(pkin(10));
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t38 = t55 * t59 - t56 * t61;
t40 = t55 * t61 + t56 * t59;
t58 = sin(qJ(4));
t91 = cos(qJ(4));
t23 = t38 * t91 + t40 * t58;
t25 = -t38 * t58 + t40 * t91;
t49 = -t61 * pkin(2) - pkin(1);
t29 = t38 * pkin(3) + t49;
t68 = -t25 * qJ(5) + t29;
t9 = pkin(4) * t23 + t68;
t100 = -0.2e1 * t9;
t21 = t23 ^ 2;
t22 = t25 ^ 2;
t92 = t56 * pkin(2);
t48 = pkin(3) + t92;
t93 = t55 * pkin(2);
t36 = t48 * t58 + t91 * t93;
t32 = qJ(5) + t36;
t99 = t32 ^ 2;
t98 = 0.2e1 * t29;
t97 = 0.2e1 * t32;
t96 = 0.2e1 * t49;
t95 = 0.2e1 * t61;
t63 = 0.2e1 * qJ(5);
t94 = pkin(4) + pkin(9);
t90 = t25 * t23;
t89 = t32 * t23;
t57 = sin(qJ(6));
t50 = t57 ^ 2;
t88 = t50 * t23;
t87 = t57 * t23;
t86 = t57 * t25;
t60 = cos(qJ(6));
t85 = t60 * t23;
t20 = t60 * t25;
t84 = t60 * t57;
t83 = -qJ(3) - pkin(7);
t52 = t60 ^ 2;
t82 = t50 + t52;
t51 = t59 ^ 2;
t53 = t61 ^ 2;
t81 = t51 + t53;
t80 = qJ(5) * t23;
t79 = t32 * qJ(5);
t78 = qJ(5) + t32;
t77 = -0.2e1 * t90;
t76 = 0.2e1 * t90;
t42 = t83 * t59;
t43 = t83 * t61;
t28 = t42 * t55 - t43 * t56;
t17 = -pkin(8) * t38 + t28;
t27 = t42 * t56 + t43 * t55;
t69 = -pkin(8) * t40 + t27;
t10 = t17 * t58 - t69 * t91;
t12 = t17 * t91 + t58 * t69;
t75 = t10 ^ 2 + t12 ^ 2;
t41 = t82 * t94;
t6 = t23 * t94 + t68;
t7 = pkin(5) * t25 + t10;
t2 = -t57 * t6 + t60 * t7;
t3 = t57 * t7 + t6 * t60;
t1 = t2 * t60 + t3 * t57;
t74 = -t2 * t57 + t3 * t60;
t73 = -t48 * t91 + t58 * t93;
t34 = -pkin(4) + t73;
t31 = -pkin(9) + t34;
t72 = -t25 * t31 + t89;
t71 = t25 * t94 + t80;
t70 = 0.2e1 * t10 * t25 - 0.2e1 * t12 * t23;
t65 = -0.2e1 * pkin(4);
t64 = qJ(5) ^ 2;
t46 = -0.2e1 * t84;
t19 = t52 * t23;
t18 = t23 * t84;
t15 = t82 * t31;
t14 = t19 - t88;
t8 = -pkin(5) * t23 + t12;
t5 = t8 * t60;
t4 = t8 * t57;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t51, t59 * t95, 0, t53, 0, 0, pkin(1) * t95, -0.2e1 * pkin(1) * t59, 0.2e1 * t81 * pkin(7), pkin(7) ^ 2 * t81 + pkin(1) ^ 2, t40 ^ 2, -0.2e1 * t40 * t38, 0, t38 ^ 2, 0, 0, t38 * t96, t40 * t96, -0.2e1 * t27 * t40 - 0.2e1 * t28 * t38, t27 ^ 2 + t28 ^ 2 + t49 ^ 2, t22, t77, 0, t21, 0, 0, t23 * t98, t25 * t98, t70, t29 ^ 2 + t75, 0, 0, 0, t22, t77, t21, t70, t23 * t100, t25 * t100, t9 ^ 2 + t75, t50 * t21, 0.2e1 * t21 * t84, t57 * t76, t52 * t21, t60 * t76, t22, 0.2e1 * t2 * t25 - 0.2e1 * t8 * t85, -0.2e1 * t25 * t3 + 0.2e1 * t8 * t87, 0.2e1 * t74 * t23, t2 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t61, 0, -t59 * pkin(7), -t61 * pkin(7), 0, 0, 0, 0, t40, 0, -t38, 0, t27, -t28 (-t38 * t55 - t40 * t56) * pkin(2) (t27 * t56 + t28 * t55) * pkin(2), 0, 0, t25, 0, -t23, 0, -t10, -t12, -t23 * t36 + t25 * t73, t10 * t73 + t12 * t36, 0, -t25, t23, 0, 0, 0, t25 * t34 - t89, t10, t12, t10 * t34 + t12 * t32, t18, t14, t20, -t18, -t86, 0, -t60 * t72 + t4, t57 * t72 + t5, -t1, t1 * t31 + t8 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t92, -0.2e1 * t93, 0 (t55 ^ 2 + t56 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t73, -0.2e1 * t36, 0, t36 ^ 2 + t73 ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t34, t97, t34 ^ 2 + t99, t52, t46, 0, t50, 0, 0, t57 * t97, t60 * t97, -0.2e1 * t15, t31 ^ 2 * t82 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, 0, t49, 0, 0, 0, 0, 0, 0, t23, t25, 0, t29, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, t9, 0, 0, 0, 0, 0, 0, -t86, -t20, t19 + t88, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, 0, -t10, -t12, 0, 0, 0, -t25, t23, 0, 0, 0, -pkin(4) * t25 - t80, t10, t12, -pkin(4) * t10 + qJ(5) * t12, t18, t14, t20, -t18, -t86, 0, -t60 * t71 + t4, t57 * t71 + t5, -t1, t8 * qJ(5) - t1 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t73, -t36, 0, 0, 1, 0, 0, 0, 0, 0, 0, t65 + t73, t63 + t36, -pkin(4) * t34 + t79, t52, t46, 0, t50, 0, 0, t78 * t57, t78 * t60 (-t31 + t94) * t82, -t31 * t41 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, t65, t63, pkin(4) ^ 2 + t64, t52, t46, 0, t50, 0, 0, t57 * t63, t60 * t63, 0.2e1 * t41, t82 * t94 ^ 2 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, t10, 0, 0, 0, 0, 0, 0, t20, -t86, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, t85, t25, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, -t57, 0, t60 * t31, -t57 * t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, -t57, 0, -t60 * t94, t57 * t94, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t11;