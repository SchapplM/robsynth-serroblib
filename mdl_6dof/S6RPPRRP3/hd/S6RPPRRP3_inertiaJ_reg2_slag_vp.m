% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t51 = cos(qJ(5));
t45 = t51 ^ 2;
t52 = cos(qJ(4));
t38 = t45 * t52;
t49 = sin(qJ(5));
t43 = t49 ^ 2;
t75 = t43 * t52;
t17 = t38 + t75;
t65 = t43 + t45;
t89 = 2 * pkin(5);
t47 = sin(pkin(9));
t82 = t47 * pkin(1);
t31 = qJ(3) + t82;
t88 = t31 ^ 2;
t87 = 0.2e1 * t31;
t86 = -0.2e1 * t49;
t85 = 0.2e1 * t51;
t84 = 0.2e1 * t52;
t50 = sin(qJ(4));
t83 = pkin(8) * t50;
t48 = cos(pkin(9));
t81 = t48 * pkin(1);
t80 = t49 * pkin(8);
t79 = t50 * pkin(4);
t78 = t51 * pkin(8);
t77 = t52 * pkin(4);
t12 = -t52 * pkin(8) + t31 + t79;
t33 = -pkin(2) - t81;
t26 = -pkin(7) + t33;
t39 = t51 * t50;
t4 = t49 * t12 + t26 * t39;
t76 = t26 * t49;
t46 = t52 ^ 2;
t74 = t46 * t26;
t73 = t49 * t51;
t36 = t49 * t52;
t72 = t50 * t26;
t40 = t51 * t52;
t56 = -t51 * pkin(5) - t49 * qJ(6);
t20 = -pkin(4) + t56;
t71 = t52 * t20;
t70 = t52 * t26;
t69 = t52 * t50;
t68 = t65 * t83;
t67 = t17 * pkin(8);
t66 = t65 * pkin(8) ^ 2;
t44 = t50 ^ 2;
t64 = t44 + t46;
t63 = t50 * qJ(6);
t62 = t49 * t69;
t61 = t46 * t73;
t60 = -t77 - t83;
t1 = t63 + t4;
t9 = t51 * t12;
t2 = -t9 + (-pkin(5) + t76) * t50;
t59 = t1 * t51 + t2 * t49;
t3 = -t49 * t72 + t9;
t58 = -t3 * t49 + t4 * t51;
t57 = -t71 + t83;
t55 = -pkin(5) * t49 + t51 * qJ(6);
t37 = t45 * t46;
t35 = t49 * t50;
t34 = t43 * t46;
t24 = t26 ^ 2;
t23 = t49 * t40;
t22 = t69 * t85;
t21 = t46 * t24;
t19 = 0.2e1 * t65 * pkin(8);
t18 = t64 * t51;
t16 = t65 * t50;
t15 = t64 * t49;
t14 = -t38 + t75;
t11 = t37 + t34 + t44;
t10 = t44 * t65 + t46;
t7 = t64 * t26;
t6 = (-0.1e1 + t65) * t69;
t5 = (-t26 - t55) * t52;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t81, -0.2e1 * t82, 0 (t47 ^ 2 + t48 ^ 2) * pkin(1) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t33, t87, t33 ^ 2 + t88, t46, -0.2e1 * t69, 0, t44, 0, 0, t50 * t87, t31 * t84, -0.2e1 * t7, t44 * t24 + t21 + t88, t37, -0.2e1 * t61, t22, t34, -0.2e1 * t62, t44, 0.2e1 * t3 * t50 - 0.2e1 * t49 * t74, -0.2e1 * t4 * t50 - 0.2e1 * t51 * t74 (-t3 * t51 - t4 * t49) * t84, t3 ^ 2 + t4 ^ 2 + t21, t37, t22, 0.2e1 * t61, t44, 0.2e1 * t62, t34, -0.2e1 * t2 * t50 + 0.2e1 * t36 * t5 (-t1 * t49 + t2 * t51) * t84, 0.2e1 * t1 * t50 - 0.2e1 * t40 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t58 - t72) * t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t50 + t52 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t7, 0, 0, 0, 0, 0, 0, -t15, -t18, 0, t50 * t58 + t74, 0, 0, 0, 0, 0, 0, -t15, 0, t18, -t5 * t52 + t50 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, -t50, 0, t70, -t72, 0, 0, t23, -t14, t35, -t23, t39, 0, t49 * t60 + t51 * t70, -t49 * t70 + t51 * t60, t58, pkin(4) * t70 + pkin(8) * t58, t23, t35, t14, 0, -t39, -t23, -t49 * t57 - t5 * t51, t59, -t5 * t49 + t51 * t57, pkin(8) * t59 + t5 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t52, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t35, t17, t67 - t79, 0, 0, 0, 0, 0, 0, -t39, t17, -t35, t50 * t20 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t50, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t36, t16, t68 + t77, 0, 0, 0, 0, 0, 0, t40, t16, t36, t68 - t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, 0.2e1 * t73, 0, t45, 0, 0, pkin(4) * t85, pkin(4) * t86, t19, pkin(4) ^ 2 + t66, t43, 0, -0.2e1 * t73, 0, 0, t45, -0.2e1 * t20 * t51, t19, t20 * t86, t20 ^ 2 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t36, t50, t3, -t4, 0, 0, 0, t40, 0, t50, t36, 0, t9 + (t89 - t76) * t50, t56 * t52, 0.2e1 * t63 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t40, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, t40, t55 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t39, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, t39, t55 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t51, 0, -t80, -t78, 0, 0, 0, t49, 0, 0, -t51, 0, -t80, t55, t78, t55 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t89, 0, 0.2e1 * qJ(6) (pkin(5) ^ 2) + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t40, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t8;
