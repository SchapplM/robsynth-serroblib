% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t45 = sin(qJ(5));
t49 = cos(qJ(5));
t47 = sin(qJ(3));
t50 = cos(qJ(4));
t56 = t50 * t47;
t46 = sin(qJ(4));
t60 = t46 * t47;
t21 = -t45 * t60 + t49 * t56;
t71 = -0.2e1 * t21;
t36 = -t50 * pkin(4) - pkin(3);
t70 = 0.2e1 * t36;
t69 = -0.2e1 * t47;
t51 = cos(qJ(3));
t68 = 0.2e1 * t51;
t67 = pkin(9) + pkin(10);
t66 = pkin(3) * t50;
t65 = pkin(8) * t46;
t64 = t45 * pkin(4);
t63 = t51 * pkin(4);
t43 = sin(pkin(6));
t62 = t43 * sin(qJ(2));
t61 = t43 * cos(qJ(2));
t59 = t46 * t50;
t58 = t46 * t51;
t29 = -t51 * pkin(3) - t47 * pkin(9) - pkin(2);
t55 = t50 * t51;
t53 = pkin(8) * t55;
t14 = t53 + (-pkin(10) * t47 + t29) * t46;
t57 = t49 * t14;
t37 = t47 * pkin(8);
t28 = pkin(4) * t60 + t37;
t54 = t47 * t68;
t24 = t50 * t29;
t10 = -pkin(10) * t56 + t24 + (-pkin(4) - t65) * t51;
t3 = t49 * t10 - t45 * t14;
t30 = t67 * t46;
t31 = t67 * t50;
t15 = -t49 * t30 - t45 * t31;
t4 = t45 * t10 + t57;
t16 = -t45 * t30 + t49 * t31;
t26 = t45 * t50 + t49 * t46;
t44 = cos(pkin(6));
t42 = t51 ^ 2;
t41 = t50 ^ 2;
t40 = t47 ^ 2;
t39 = t46 ^ 2;
t38 = t49 * pkin(4);
t35 = t38 + pkin(5);
t25 = t45 * t46 - t49 * t50;
t23 = t44 * t47 + t51 * t62;
t22 = -t44 * t51 + t47 * t62;
t20 = t26 * t47;
t19 = t25 * pkin(5) + t36;
t18 = t46 * t29 + t53;
t17 = -pkin(8) * t58 + t24;
t13 = t23 * t50 - t46 * t61;
t12 = -t23 * t46 - t50 * t61;
t11 = t20 * pkin(5) + t28;
t8 = -t25 * qJ(6) + t16;
t7 = -t26 * qJ(6) + t15;
t6 = t45 * t12 + t49 * t13;
t5 = t49 * t12 - t45 * t13;
t2 = -t20 * qJ(6) + t4;
t1 = -t51 * pkin(5) - t21 * qJ(6) + t3;
t9 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, t61, -t62, 0, 0, 0, 0, 0, t51 * t61, -t47 * t61, 0, 0, 0, 0, 0, -t12 * t51 + t22 * t60, t13 * t51 + t22 * t56, 0, 0, 0, 0, 0, t22 * t20 - t5 * t51, t22 * t21 + t6 * t51, -t6 * t20 - t5 * t21, t5 * t1 + t22 * t11 + t6 * t2; 0, 1, 0, 0, t40, t54, 0, 0, 0, pkin(2) * t68, pkin(2) * t69, t41 * t40, -0.2e1 * t40 * t59, t55 * t69, t46 * t54, t42, -0.2e1 * t17 * t51 + 0.2e1 * t40 * t65, 0.2e1 * t40 * pkin(8) * t50 + 0.2e1 * t18 * t51, t21 ^ 2, t20 * t71, t51 * t71, t20 * t68, t42, 0.2e1 * t28 * t20 - 0.2e1 * t3 * t51, 0.2e1 * t28 * t21 + 0.2e1 * t4 * t51, -0.2e1 * t1 * t21 - 0.2e1 * t2 * t20, t1 ^ 2 + t11 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, 0, 0, 0, 0, -t22 * t50, t22 * t46, 0, 0, 0, 0, 0, t22 * t25, t22 * t26, -t6 * t25 - t5 * t26, t22 * t19 + t5 * t7 + t6 * t8; 0, 0, 0, 0, 0, 0, t47, t51, 0, -t37, -t51 * pkin(8), t46 * t56 (-t39 + t41) * t47, -t58, -t55, 0, -pkin(8) * t56 + (-pkin(3) * t47 + pkin(9) * t51) * t46, pkin(9) * t55 + (t65 - t66) * t47, t21 * t26, -t26 * t20 - t21 * t25, -t26 * t51, t25 * t51, 0, -t15 * t51 + t36 * t20 + t28 * t25, t16 * t51 + t36 * t21 + t28 * t26, -t1 * t26 - t2 * t25 - t8 * t20 - t7 * t21, t1 * t7 + t11 * t19 + t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t39, 0.2e1 * t59, 0, 0, 0, 0.2e1 * t66, -0.2e1 * pkin(3) * t46, t26 ^ 2, -0.2e1 * t26 * t25, 0, 0, 0, t25 * t70, t26 * t70, -0.2e1 * t8 * t25 - 0.2e1 * t7 * t26, t19 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, 0, 0, 0, 0, 0, t5, -t6, 0, t5 * t35 + t6 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t60, -t51, t17, -t18, 0, 0, t21, -t20, -t51, -t49 * t63 + t3, -t57 + (-t10 + t63) * t45, -t20 * t64 - t35 * t21, t1 * t35 + t2 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t50, 0, -t46 * pkin(9), -t50 * pkin(9), 0, 0, t26, -t25, 0, t15, -t16, -t25 * t64 - t35 * t26, t7 * t35 + t8 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t64, 0, t45 ^ 2 * pkin(4) ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, -t51, t3, -t4, -pkin(5) * t21, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t15, -t16, -pkin(5) * t26, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t64, 0, t35 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t9;
