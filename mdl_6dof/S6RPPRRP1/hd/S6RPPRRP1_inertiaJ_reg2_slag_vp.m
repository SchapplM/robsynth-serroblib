% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t51 = cos(pkin(10));
t54 = sin(qJ(4));
t49 = sin(pkin(10));
t79 = cos(qJ(4));
t66 = t79 * t49;
t34 = t54 * t51 + t66;
t90 = -0.2e1 * t34;
t50 = sin(pkin(9));
t83 = t50 * pkin(1);
t41 = qJ(3) + t83;
t80 = pkin(7) + t41;
t28 = t80 * t51;
t10 = t54 * t28 + t80 * t66;
t89 = t10 ^ 2;
t74 = t54 * t49;
t32 = -t79 * t51 + t74;
t88 = t32 ^ 2;
t52 = cos(pkin(9));
t82 = t52 * pkin(1);
t43 = -pkin(2) - t82;
t35 = -t51 * pkin(3) + t43;
t87 = 0.2e1 * t35;
t86 = 0.2e1 * t49;
t85 = t32 * pkin(4);
t84 = t32 * pkin(5);
t55 = cos(qJ(5));
t81 = t55 * pkin(5);
t78 = t10 * t32;
t53 = sin(qJ(5));
t77 = t53 * t32;
t76 = t53 * t34;
t75 = t53 * t55;
t12 = t79 * t28 - t80 * t74;
t73 = t55 * t12;
t27 = t55 * t34;
t72 = -qJ(6) - pkin(8);
t45 = t49 ^ 2;
t46 = t51 ^ 2;
t71 = t45 + t46;
t47 = t53 ^ 2;
t48 = t55 ^ 2;
t38 = t47 + t48;
t70 = qJ(6) * t34;
t69 = t32 * t90;
t68 = pkin(5) * t76;
t67 = t38 * pkin(8);
t9 = -t34 * pkin(8) + t35 + t85;
t3 = -t53 * t12 + t55 * t9;
t65 = -pkin(4) * t34 - pkin(8) * t32;
t58 = -t55 * t70 + t3;
t1 = t58 + t84;
t2 = t73 + (t9 - t70) * t53;
t64 = t1 * t55 + t2 * t53;
t63 = -t1 * t53 + t2 * t55;
t4 = t53 * t9 + t73;
t62 = t3 * t55 + t4 * t53;
t61 = -t3 * t53 + t4 * t55;
t36 = t72 * t53;
t37 = t72 * t55;
t60 = t55 * t36 - t53 * t37;
t59 = -t36 * t53 - t37 * t55;
t44 = -pkin(4) - t81;
t40 = 0.2e1 * t75;
t31 = t34 ^ 2;
t26 = t55 * t32;
t25 = t48 * t34;
t24 = t48 * t31;
t21 = t47 * t34;
t20 = t47 * t31;
t18 = t53 * t27;
t17 = -0.2e1 * t31 * t75;
t16 = 0.2e1 * t32 * t27;
t15 = t53 * t69;
t14 = -t21 - t25;
t13 = -t21 + t25;
t7 = t24 + t20 + t88;
t6 = t10 + t68;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t82, -0.2e1 * t83, 0 (t50 ^ 2 + t52 ^ 2) * pkin(1) ^ 2, t45, t51 * t86, 0, t46, 0, 0, -0.2e1 * t43 * t51, t43 * t86, 0.2e1 * t71 * t41, t71 * t41 ^ 2 + t43 ^ 2, t31, t69, 0, t88, 0, 0, t32 * t87, t34 * t87, 0.2e1 * t10 * t34 - 0.2e1 * t12 * t32, t12 ^ 2 + t35 ^ 2 + t89, t24, t17, t16, t20, t15, t88, 0.2e1 * t10 * t76 + 0.2e1 * t3 * t32, 0.2e1 * t10 * t27 - 0.2e1 * t4 * t32, t62 * t90, t3 ^ 2 + t4 ^ 2 + t89, t24, t17, t16, t20, t15, t88, 0.2e1 * t1 * t32 + 0.2e1 * t6 * t76, -0.2e1 * t2 * t32 + 0.2e1 * t6 * t27, t64 * t90, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t34 + t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t34 + t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t32 + t63 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t49, 0, t43, 0, 0, 0, 0, 0, 0, t32, t34, 0, t35, 0, 0, 0, 0, 0, 0, t26, -t77, t14, t62, 0, 0, 0, 0, 0, 0, t26, -t77, t14, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t32, 0, -t10, -t12, 0, 0, t18, t13, t77, -t18, t26, 0, -t10 * t55 + t65 * t53, t10 * t53 + t65 * t55, t61, -t10 * pkin(4) + t61 * pkin(8), t18, t13, t77, -t18, t26, 0, t36 * t32 + t44 * t76 - t6 * t55, t44 * t27 + t37 * t32 + t6 * t53, -t60 * t34 + t63, t1 * t36 - t2 * t37 + t6 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t34, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t77, -t14, t34 * t67 - t85, 0, 0, 0, 0, 0, 0, -t26, t77, -t14, t32 * t44 + t59 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t47, t40, 0, t48, 0, 0, 0.2e1 * pkin(4) * t55, -0.2e1 * pkin(4) * t53, 0.2e1 * t67, t38 * pkin(8) ^ 2 + pkin(4) ^ 2, t47, t40, 0, t48, 0, 0, -0.2e1 * t44 * t55, 0.2e1 * t44 * t53, 0.2e1 * t59, t36 ^ 2 + t37 ^ 2 + t44 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t76, t32, t3, -t4, 0, 0, 0, 0, t27, 0, -t76, t32, t58 + 0.2e1 * t84, -t2, -pkin(5) * t27, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t27, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t27, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t53, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t53, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t55, 0, -t53 * pkin(8), -t55 * pkin(8), 0, 0, 0, 0, t53, 0, t55, 0, t36, t37, -t53 * pkin(5), t36 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(5), 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t27, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t53, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;