% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t44 = sin(qJ(4));
t45 = sin(qJ(3));
t48 = cos(qJ(4));
t49 = cos(qJ(3));
t27 = t44 * t45 - t48 * t49;
t28 = t44 * t49 + t48 * t45;
t38 = -t49 * pkin(3) - pkin(2);
t54 = -t28 * qJ(5) + t38;
t15 = t27 * pkin(4) + t54;
t77 = -0.2e1 * t15;
t71 = t44 * pkin(3);
t34 = qJ(5) + t71;
t76 = 0.2e1 * t34;
t75 = 0.2e1 * t38;
t74 = 0.2e1 * t49;
t52 = 0.2e1 * qJ(5);
t73 = pkin(4) + pkin(10);
t72 = -pkin(9) - pkin(8);
t70 = t48 * pkin(3);
t69 = t28 * t27;
t68 = t34 * t27;
t41 = sin(pkin(6));
t67 = t41 * sin(qJ(2));
t50 = cos(qJ(2));
t66 = t41 * t50;
t43 = sin(qJ(6));
t65 = t43 * t27;
t64 = t43 * t28;
t47 = cos(qJ(6));
t63 = t47 * t27;
t62 = t47 * t43;
t61 = qJ(5) * t27;
t60 = qJ(5) + t34;
t59 = 0.2e1 * t69;
t58 = t27 * t66;
t57 = t28 * t66;
t37 = -pkin(4) - t70;
t30 = t72 * t45;
t31 = t72 * t49;
t17 = -t48 * t30 - t44 * t31;
t32 = -pkin(10) + t37;
t56 = -t28 * t32 + t68;
t18 = t44 * t30 - t48 * t31;
t55 = t28 * t73 + t61;
t53 = -0.2e1 * pkin(4);
t42 = cos(pkin(6));
t40 = t47 ^ 2;
t39 = t43 ^ 2;
t33 = -0.2e1 * t62;
t26 = t28 ^ 2;
t25 = t27 ^ 2;
t24 = t47 * t28;
t23 = t42 * t45 + t49 * t67;
t22 = t42 * t49 - t45 * t67;
t21 = t27 * t62;
t16 = (-t39 + t40) * t27;
t13 = t44 * t22 + t48 * t23;
t12 = -t48 * t22 + t44 * t23;
t11 = -t27 * pkin(5) + t18;
t10 = t28 * pkin(5) + t17;
t9 = t13 * t47;
t8 = t13 * t43;
t7 = t11 * t47;
t6 = t11 * t43;
t5 = t73 * t27 + t54;
t4 = -t43 * t12 + t47 * t66;
t3 = t47 * t12 + t43 * t66;
t2 = t43 * t10 + t47 * t5;
t1 = t47 * t10 - t43 * t5;
t14 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 ^ 2 * t50 ^ 2 + t12 ^ 2 + t13 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t66, -t67, 0, 0, 0, 0, 0, t49 * t66, -t45 * t66, 0, 0, 0, 0, 0, -t58, -t57, t12 * t28 - t13 * t27, t58, t57, t12 * t17 + t13 * t18 - t15 * t66, 0, 0, 0, 0, 0, -t13 * t63 + t3 * t28, t13 * t65 + t4 * t28; 0, 1, 0, 0, t45 ^ 2, t45 * t74, 0, 0, 0, pkin(2) * t74, -0.2e1 * pkin(2) * t45, t26, -0.2e1 * t69, 0, 0, 0, t27 * t75, t28 * t75, 0.2e1 * t17 * t28 - 0.2e1 * t18 * t27, t27 * t77, t28 * t77, t15 ^ 2 + t17 ^ 2 + t18 ^ 2, t39 * t25, 0.2e1 * t25 * t62, t43 * t59, t47 * t59, t26, 0.2e1 * t1 * t28 - 0.2e1 * t11 * t63, 0.2e1 * t11 * t65 - 0.2e1 * t2 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, 0, 0, 0, 0, 0, -t12, -t13, 0, t12, t13, t12 * t37 + t13 * t34, 0, 0, 0, 0, 0, t8, t9; 0, 0, 0, 0, 0, 0, t45, t49, 0, -t45 * pkin(8), -t49 * pkin(8), 0, 0, t28, -t27, 0, -t17, -t18, t37 * t28 - t68, t17, t18, t17 * t37 + t18 * t34, t21, t16, t24, -t64, 0, -t56 * t47 + t6, t56 * t43 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t70, -0.2e1 * t71, 0, 0.2e1 * t37, t76, t34 ^ 2 + t37 ^ 2, t40, t33, 0, 0, 0, t43 * t76, t47 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, t12, t13, -t12 * pkin(4) + t13 * qJ(5), 0, 0, 0, 0, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, -t17, -t18, -pkin(4) * t28 - t61, t17, t18, -t17 * pkin(4) + t18 * qJ(5), t21, t16, t24, -t64, 0, -t55 * t47 + t6, t55 * t43 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t70, -t71, 0, t53 - t70, t52 + t71, -t37 * pkin(4) + t34 * qJ(5), t40, t33, 0, 0, 0, t60 * t43, t60 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t53, t52, pkin(4) ^ 2 + qJ(5) ^ 2, t40, t33, 0, 0, 0, t43 * t52, t47 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, t17, 0, 0, 0, 0, 0, t24, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t63, t28, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t43, 0, t47 * t32, -t43 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t43, 0, -t47 * t73, t43 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t14;
