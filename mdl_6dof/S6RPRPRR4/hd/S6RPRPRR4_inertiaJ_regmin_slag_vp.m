% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(5));
t30 = t41 * pkin(5) + qJ(4);
t70 = 0.2e1 * t30;
t42 = sin(qJ(3));
t69 = 0.2e1 * t42;
t68 = 0.2e1 * qJ(4);
t46 = -pkin(3) - pkin(8);
t40 = sin(qJ(6));
t67 = t40 * pkin(5);
t66 = t42 * pkin(5);
t43 = cos(qJ(6));
t65 = t43 * pkin(5);
t44 = cos(qJ(5));
t45 = cos(qJ(3));
t39 = cos(pkin(10));
t29 = -t39 * pkin(1) - pkin(2);
t54 = t42 * qJ(4);
t48 = t29 - t54;
t11 = t46 * t45 + t48;
t51 = pkin(9) * t45 - t11;
t38 = sin(pkin(10));
t28 = t38 * pkin(1) + pkin(7);
t24 = t42 * t28;
t18 = t42 * pkin(4) + t24;
t62 = t41 * t18;
t5 = -t51 * t44 + t62;
t64 = t43 * t5;
t63 = t45 * pkin(3);
t61 = t41 * t42;
t60 = t41 * t45;
t20 = t40 * t44 + t43 * t41;
t59 = t42 * t20;
t58 = t44 * t41;
t57 = t44 * t45;
t56 = t45 * t42;
t25 = t45 * t28;
t19 = t45 * pkin(4) + t25;
t35 = t42 ^ 2;
t37 = t45 ^ 2;
t55 = t35 + t37;
t53 = t45 * qJ(4);
t52 = -0.2e1 * t56;
t15 = t44 * t18;
t4 = t51 * t41 + t15 + t66;
t1 = t43 * t4 - t40 * t5;
t50 = -t42 * pkin(3) + t53;
t49 = t42 * t46 + t53;
t36 = t44 ^ 2;
t34 = t41 ^ 2;
t32 = t44 * t46;
t31 = t44 * t42;
t23 = -t44 * pkin(9) + t32;
t22 = (-pkin(9) + t46) * t41;
t21 = -t40 * t41 + t43 * t44;
t17 = t48 - t63;
t16 = t21 * t42;
t13 = t20 * t45;
t12 = t40 * t60 - t43 * t57;
t10 = pkin(5) * t57 + t19;
t9 = t43 * t22 + t40 * t23;
t8 = -t40 * t22 + t43 * t23;
t7 = t44 * t11 + t62;
t6 = -t41 * t11 + t15;
t2 = t40 * t4 + t64;
t3 = [1, 0, 0 (t38 ^ 2 + t39 ^ 2) * pkin(1) ^ 2, t35, 0.2e1 * t56, 0, 0, 0, -0.2e1 * t29 * t45, t29 * t69, 0.2e1 * t55 * t28, 0.2e1 * t17 * t45, -0.2e1 * t17 * t42, t55 * t28 ^ 2 + t17 ^ 2, t34 * t37, 0.2e1 * t37 * t58, t41 * t52, t44 * t52, t35, 0.2e1 * t19 * t57 + 0.2e1 * t6 * t42, -0.2e1 * t19 * t60 - 0.2e1 * t7 * t42, t13 ^ 2, -0.2e1 * t13 * t12, -t13 * t69, t12 * t69, t35, 0.2e1 * t1 * t42 - 0.2e1 * t10 * t12, -0.2e1 * t10 * t13 - 0.2e1 * t2 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t42, t45, 0, -t24, -t25, t50, t24, t25, t50 * t28, -t41 * t57 (t34 - t36) * t45, t31, -t61, 0, t19 * t41 + t49 * t44, t19 * t44 - t49 * t41, -t13 * t21, t21 * t12 + t13 * t20, t16, -t59, 0, t10 * t20 - t30 * t12 + t8 * t42, t10 * t21 - t30 * t13 - t9 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, -t45, t42, t54 + t63, 0, 0, 0, 0, 0, t61, t31, 0, 0, 0, 0, 0, t59, t16; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t68, pkin(3) ^ 2 + qJ(4) ^ 2, t36, -0.2e1 * t58, 0, 0, 0, t41 * t68, t44 * t68, t21 ^ 2, -0.2e1 * t21 * t20, 0, 0, 0, t20 * t70, t21 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, t24, 0, 0, 0, 0, 0, t31, -t61, 0, 0, 0, 0, 0, t16, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t57, t42, t6, -t7, 0, 0, -t13, t12, t42, t42 * t65 + t1, -t64 + (-t4 - t66) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t60, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, t32, -t41 * t46, 0, 0, t21, -t20, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, 0, 0, 0, 0, t21, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12, t42, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
