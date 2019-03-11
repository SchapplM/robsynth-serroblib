% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t43 = sin(pkin(11));
t45 = cos(pkin(11));
t49 = sin(qJ(4));
t60 = cos(qJ(4));
t26 = t43 * t49 - t45 * t60;
t70 = -0.2e1 * t26;
t69 = 0.2e1 * t26;
t46 = cos(pkin(10));
t36 = -pkin(1) * t46 - pkin(2);
t30 = -pkin(3) * t45 + t36;
t68 = 0.2e1 * t30;
t51 = cos(qJ(5));
t38 = -pkin(5) * t51 - pkin(4);
t67 = 0.2e1 * t38;
t66 = pkin(8) + pkin(9);
t65 = t26 * pkin(5);
t47 = sin(qJ(6));
t64 = t47 * pkin(5);
t50 = cos(qJ(6));
t63 = t50 * pkin(5);
t27 = t43 * t60 + t45 * t49;
t10 = pkin(4) * t26 - pkin(8) * t27 + t30;
t48 = sin(qJ(5));
t44 = sin(pkin(10));
t34 = pkin(1) * t44 + qJ(3);
t61 = pkin(7) + t34;
t22 = t61 * t43;
t23 = t61 * t45;
t14 = -t22 * t49 + t23 * t60;
t57 = t51 * t14;
t5 = t57 + (-pkin(9) * t27 + t10) * t48;
t62 = t50 * t5;
t28 = t47 * t48 - t50 * t51;
t15 = t26 * t28;
t29 = t47 * t51 + t48 * t50;
t16 = t29 * t26;
t20 = t48 * t26;
t59 = t48 * t27;
t58 = t48 * t51;
t56 = t51 * t27;
t55 = t43 ^ 2 + t45 ^ 2;
t54 = t27 * t70;
t6 = t10 * t51 - t14 * t48;
t4 = -pkin(9) * t56 + t6 + t65;
t1 = t4 * t50 - t47 * t5;
t53 = -pkin(4) * t27 - pkin(8) * t26;
t13 = t22 * t60 + t23 * t49;
t42 = t51 ^ 2;
t41 = t48 ^ 2;
t32 = t66 * t51;
t31 = t66 * t48;
t25 = t27 ^ 2;
t24 = t26 ^ 2;
t21 = t51 * t26;
t18 = -t31 * t47 + t32 * t50;
t17 = -t31 * t50 - t32 * t47;
t12 = -t47 * t59 + t50 * t56;
t11 = t29 * t27;
t8 = pkin(5) * t59 + t13;
t7 = t10 * t48 + t57;
t2 = t4 * t47 + t62;
t3 = [1, 0, 0 (t44 ^ 2 + t46 ^ 2) * pkin(1) ^ 2, -0.2e1 * t36 * t45, 0.2e1 * t36 * t43, 0.2e1 * t55 * t34, t34 ^ 2 * t55 + t36 ^ 2, t25, t54, 0, 0, 0, t26 * t68, t27 * t68, t42 * t25, -0.2e1 * t25 * t58, t56 * t69, t48 * t54, t24, 0.2e1 * t13 * t59 + 0.2e1 * t26 * t6, 0.2e1 * t13 * t56 - 0.2e1 * t26 * t7, t12 ^ 2, -0.2e1 * t12 * t11, t12 * t69, t11 * t70, t24, 0.2e1 * t1 * t26 + 0.2e1 * t11 * t8, 0.2e1 * t12 * t8 - 0.2e1 * t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t45, t43, 0, t36, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, t21, -t20, 0, 0, 0, 0, 0, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, -t13, -t14, t48 * t56 (-t41 + t42) * t27, t20, t21, 0, -t13 * t51 + t48 * t53, t13 * t48 + t51 * t53, t12 * t29, -t11 * t29 - t12 * t28, t16, -t15, 0, t11 * t38 + t17 * t26 + t28 * t8, t12 * t38 - t18 * t26 + t29 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, 0, 0, 0, 0, -t21, t20, 0, 0, 0, 0, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t41, 0.2e1 * t58, 0, 0, 0, 0.2e1 * pkin(4) * t51, -0.2e1 * pkin(4) * t48, t29 ^ 2, -0.2e1 * t29 * t28, 0, 0, 0, t28 * t67, t29 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t59, t26, t6, -t7, 0, 0, t12, -t11, t26, t26 * t63 + t1, -t62 + (-t4 - t65) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t56, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, 0, 0, 0, 0, -t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t48 * pkin(8), -t51 * pkin(8), 0, 0, t29, -t28, 0, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t63, -0.2e1 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, t26, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t63, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
