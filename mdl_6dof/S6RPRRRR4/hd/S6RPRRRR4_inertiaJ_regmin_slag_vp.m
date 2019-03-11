% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t50 = sin(pkin(11));
t51 = cos(pkin(11));
t55 = sin(qJ(3));
t77 = cos(qJ(3));
t31 = t55 * t50 - t77 * t51;
t39 = -t51 * pkin(2) - pkin(1);
t23 = t31 * pkin(3) + t39;
t32 = t77 * t50 + t55 * t51;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t67 = t57 * t31 + t54 * t32;
t18 = t67 * pkin(4) + t23;
t84 = 0.2e1 * t18;
t83 = 0.2e1 * t23;
t82 = 0.2e1 * t39;
t52 = sin(qJ(6));
t81 = pkin(5) * t52;
t53 = sin(qJ(5));
t70 = pkin(7) + qJ(2);
t33 = t70 * t50;
t34 = t70 * t51;
t66 = -t77 * t33 - t55 * t34;
t19 = -t32 * pkin(8) + t66;
t60 = t55 * t33 - t77 * t34;
t20 = -t31 * pkin(8) - t60;
t10 = t57 * t19 - t54 * t20;
t22 = -t54 * t31 + t57 * t32;
t59 = -t22 * pkin(9) + t10;
t76 = cos(qJ(5));
t11 = -t54 * t19 - t57 * t20;
t9 = -t67 * pkin(9) - t11;
t4 = t53 * t9 - t76 * t59;
t56 = cos(qJ(6));
t80 = t4 * t56;
t79 = t53 * pkin(4);
t78 = t54 * pkin(3);
t44 = t57 * pkin(3);
t42 = t44 + pkin(4);
t64 = -t76 * t42 + t53 * t78;
t25 = -pkin(5) + t64;
t75 = t25 * t56;
t43 = t76 * pkin(4);
t41 = -t43 - pkin(5);
t74 = t41 * t56;
t16 = t53 * t22 + t76 * t67;
t13 = t52 * t16;
t17 = t76 * t22 - t53 * t67;
t73 = t52 * t17;
t72 = t52 * t56;
t71 = t56 * t17;
t69 = t50 ^ 2 + t51 ^ 2;
t68 = -0.2e1 * t17 * t16;
t65 = t76 * t78;
t63 = -pkin(5) * t17 - pkin(10) * t16;
t28 = -t53 * t42 - t65;
t26 = pkin(10) - t28;
t62 = -t16 * t26 + t17 * t25;
t40 = pkin(10) + t79;
t61 = -t16 * t40 + t17 * t41;
t49 = t56 ^ 2;
t48 = t52 ^ 2;
t45 = pkin(5) * t56;
t37 = 0.2e1 * t72;
t36 = t41 * t52;
t24 = t25 * t52;
t15 = t17 ^ 2;
t14 = t56 * t16;
t12 = t52 * t71;
t7 = (-t48 + t49) * t17;
t6 = t16 * pkin(5) - t17 * pkin(10) + t18;
t5 = t53 * t59 + t76 * t9;
t3 = t4 * t52;
t2 = t56 * t5 + t52 * t6;
t1 = -t52 * t5 + t56 * t6;
t8 = [1, 0, 0, 0.2e1 * pkin(1) * t51, -0.2e1 * pkin(1) * t50, 0.2e1 * t69 * qJ(2), t69 * qJ(2) ^ 2 + pkin(1) ^ 2, t32 ^ 2, -0.2e1 * t32 * t31, 0, 0, 0, t31 * t82, t32 * t82, t22 ^ 2, -0.2e1 * t22 * t67, 0, 0, 0, t67 * t83, t22 * t83, t15, t68, 0, 0, 0, t16 * t84, t17 * t84, t49 * t15, -0.2e1 * t15 * t72, 0.2e1 * t16 * t71, t52 * t68, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t4 * t73, -0.2e1 * t2 * t16 + 0.2e1 * t4 * t71; 0, 0, 0, -t51, t50, 0, -pkin(1), 0, 0, 0, 0, 0, t31, t32, 0, 0, 0, 0, 0, t67, t22, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, t14, -t13; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, t66, t60, 0, 0, t22, -t67, 0, t10, t11, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t52 * t62 - t80, t56 * t62 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t78, 0, 0, 0, 0, 1, -0.2e1 * t64, 0.2e1 * t28, t48, t37, 0, 0, 0, -0.2e1 * t75, 0.2e1 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t67, 0, t10, t11, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t52 * t61 - t80, t56 * t61 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t44, -t78, 0, 0, 0, 0, 1, t43 - t64, -t65 + (-pkin(4) - t42) * t53, t48, t37, 0, 0, 0 (-t25 - t41) * t56, t36 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t79, t48, t37, 0, 0, 0, -0.2e1 * t74, 0.2e1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t52 * t63 - t80, t56 * t63 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t64, t28, t48, t37, 0, 0, 0, t45 - t75, t24 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t43, -t79, t48, t37, 0, 0, 0, t45 - t74, t36 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t48, t37, 0, 0, 0, 0.2e1 * t45, -0.2e1 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t73, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t56, 0, -t52 * t26, -t56 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t56, 0, -t52 * t40, -t56 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t56, 0, -t52 * pkin(10), -t56 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
