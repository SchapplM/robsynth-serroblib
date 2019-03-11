% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t42 = sin(pkin(13));
t46 = cos(pkin(13));
t50 = sin(qJ(6));
t53 = cos(qJ(6));
t28 = t42 * t50 - t46 * t53;
t51 = sin(qJ(4));
t22 = t28 * t51;
t72 = 0.2e1 * t22;
t36 = -pkin(5) * t46 - pkin(4);
t71 = 0.2e1 * t36;
t54 = cos(qJ(4));
t70 = 0.2e1 * t54;
t69 = pkin(9) * t42;
t37 = t51 * pkin(9);
t68 = t54 * pkin(9);
t67 = t42 * t51;
t44 = sin(pkin(7));
t52 = sin(qJ(3));
t66 = t44 * t52;
t55 = cos(qJ(3));
t65 = t44 * t55;
t64 = t46 * t51;
t47 = cos(pkin(12));
t48 = cos(pkin(7));
t63 = t47 * t48;
t62 = pkin(10) + qJ(5);
t31 = -pkin(4) * t54 - qJ(5) * t51 - pkin(3);
t20 = t31 * t42 + t46 * t68;
t61 = t42 ^ 2 + t46 ^ 2;
t43 = sin(pkin(12));
t45 = sin(pkin(6));
t49 = cos(pkin(6));
t13 = t49 * t66 + (t43 * t55 + t52 * t63) * t45;
t23 = -t44 * t45 * t47 + t48 * t49;
t10 = t13 * t54 + t23 * t51;
t12 = -t49 * t65 + (t43 * t52 - t55 * t63) * t45;
t3 = -t10 * t42 + t12 * t46;
t4 = t10 * t46 + t12 * t42;
t60 = -t3 * t42 + t4 * t46;
t59 = -pkin(4) * t51 + qJ(5) * t54;
t25 = t48 * t51 + t54 * t66;
t14 = -t25 * t42 - t46 * t65;
t15 = t25 * t46 - t42 * t65;
t58 = -t14 * t42 + t15 * t46;
t27 = t46 * t31;
t19 = -t42 * t68 + t27;
t57 = -t19 * t42 + t20 * t46;
t29 = t42 * t53 + t46 * t50;
t41 = t51 ^ 2;
t33 = t62 * t46;
t32 = t62 * t42;
t30 = pkin(5) * t67 + t37;
t24 = -t48 * t54 + t51 * t66;
t21 = t29 * t51;
t18 = -t32 * t50 + t33 * t53;
t17 = -t32 * t53 - t33 * t50;
t16 = -pkin(10) * t67 + t20;
t11 = -pkin(10) * t64 + t27 + (-pkin(5) - t69) * t54;
t9 = t13 * t51 - t23 * t54;
t8 = t14 * t50 + t15 * t53;
t7 = t14 * t53 - t15 * t50;
t6 = t11 * t50 + t16 * t53;
t5 = t11 * t53 - t16 * t50;
t2 = t3 * t50 + t4 * t53;
t1 = t3 * t53 - t4 * t50;
t26 = [1, t49 ^ 2 + (t43 ^ 2 + t47 ^ 2) * t45 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t3 + t15 * t4 + t24 * t9, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t15 ^ 2 + t24 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t12 * t54, t12 * t51, -t3 * t54 + t67 * t9, t4 * t54 + t64 * t9 (-t3 * t46 - t4 * t42) * t51, t19 * t3 + t20 * t4 + t37 * t9, 0, 0, 0, 0, 0, -t1 * t54 + t21 * t9, t2 * t54 - t22 * t9; 0, 0, 0, t65, -t66, 0, 0, 0, 0, 0, t54 * t65, -t51 * t65, -t14 * t54 + t24 * t67, t15 * t54 + t24 * t64 (-t14 * t46 - t15 * t42) * t51, t14 * t19 + t15 * t20 + t24 * t37, 0, 0, 0, 0, 0, t21 * t24 - t54 * t7, -t22 * t24 + t54 * t8; 0, 0, 1, 0, 0, t41, t51 * t70, 0, 0, 0, pkin(3) * t70, -0.2e1 * pkin(3) * t51, -0.2e1 * t19 * t54 + 0.2e1 * t41 * t69, 0.2e1 * pkin(9) * t41 * t46 + 0.2e1 * t20 * t54, 0.2e1 * (-t19 * t46 - t20 * t42) * t51, pkin(9) ^ 2 * t41 + t19 ^ 2 + t20 ^ 2, t22 ^ 2, t21 * t72, t54 * t72, t21 * t70, t54 ^ 2, 0.2e1 * t21 * t30 - 0.2e1 * t5 * t54, -0.2e1 * t22 * t30 + 0.2e1 * t54 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -t9 * t46, t9 * t42, t60, -t9 * pkin(4) + qJ(5) * t60, 0, 0, 0, 0, 0, t9 * t28, t9 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, -t24 * t46, t24 * t42, t58, -t24 * pkin(4) + qJ(5) * t58, 0, 0, 0, 0, 0, t24 * t28, t24 * t29; 0, 0, 0, 0, 0, 0, 0, t51, t54, 0, -t37, -t68, -pkin(9) * t64 + t42 * t59, pkin(9) * t67 + t46 * t59, t57, -pkin(4) * t37 + qJ(5) * t57, -t22 * t29, -t21 * t29 + t22 * t28, -t29 * t54, t28 * t54, 0, -t17 * t54 + t21 * t36 + t28 * t30, t18 * t54 - t22 * t36 + t29 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t46, -0.2e1 * pkin(4) * t42, 0.2e1 * t61 * qJ(5), qJ(5) ^ 2 * t61 + pkin(4) ^ 2, t29 ^ 2, -0.2e1 * t29 * t28, 0, 0, 0, t28 * t71, t29 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t64, 0, t37, 0, 0, 0, 0, 0, t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t42, 0, -pkin(4), 0, 0, 0, 0, 0, t28, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t21, -t54, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t26;
