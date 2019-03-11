% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t42 = sin(qJ(5));
t43 = sin(qJ(4));
t45 = cos(qJ(5));
t46 = cos(qJ(4));
t23 = t42 * t46 + t45 * t43;
t47 = cos(qJ(3));
t17 = t47 * t23;
t78 = -0.2e1 * t17;
t77 = -0.2e1 * t23;
t35 = -t46 * pkin(4) - pkin(3);
t76 = 0.2e1 * t35;
t75 = 0.2e1 * t46;
t74 = 2 * qJ(2);
t73 = -pkin(9) - pkin(8);
t37 = t42 * pkin(4);
t44 = sin(qJ(3));
t72 = t44 * pkin(5);
t71 = t45 * pkin(4);
t25 = t44 * pkin(3) - t47 * pkin(8) + qJ(2);
t48 = -pkin(1) - pkin(7);
t62 = t46 * t48;
t56 = t44 * t62;
t11 = t56 + (-pkin(9) * t47 + t25) * t43;
t20 = t46 * t25;
t32 = t46 * t47;
t65 = t43 * t48;
t9 = -pkin(9) * t32 + t20 + (pkin(4) - t65) * t44;
t4 = t45 * t11 + t42 * t9;
t26 = t73 * t46;
t55 = t73 * t43;
t12 = -t42 * t26 - t45 * t55;
t70 = t12 * t44;
t13 = -t45 * t26 + t42 * t55;
t69 = t13 * t44;
t16 = t23 * t44;
t68 = t43 * t44;
t67 = t43 * t46;
t66 = t43 * t47;
t64 = t44 * t48;
t63 = t46 * t44;
t22 = t42 * t43 - t45 * t46;
t61 = t47 * t22;
t60 = t47 * t44;
t59 = t47 * t48;
t39 = t44 ^ 2;
t41 = t47 ^ 2;
t58 = -t39 - t41;
t57 = -0.2e1 * t60;
t36 = t44 * qJ(6);
t1 = t36 + t4;
t54 = t42 * t11 - t45 * t9;
t21 = pkin(4) * t66 - t59;
t53 = -pkin(3) * t47 - pkin(8) * t44;
t52 = -t16 * t44 - t47 * t17;
t18 = -t42 * t68 + t45 * t63;
t19 = t45 * t32 - t42 * t66;
t51 = t18 * t44 + t47 * t19;
t50 = 0.2e1 * pkin(5);
t49 = 0.2e1 * qJ(6);
t40 = t46 ^ 2;
t38 = t43 ^ 2;
t33 = pkin(5) + t71;
t30 = t37 + qJ(6);
t15 = t43 * t25 + t56;
t14 = -t43 * t64 + t20;
t8 = t22 * pkin(5) - t23 * qJ(6) + t35;
t5 = t17 * pkin(5) - t19 * qJ(6) + t21;
t2 = t54 - t72;
t3 = [1, 0, 0, -2 * pkin(1), t74, pkin(1) ^ 2 + qJ(2) ^ 2, t41, t57, 0, 0, 0, t44 * t74, t47 * t74, t40 * t41, -0.2e1 * t41 * t67, t60 * t75, t43 * t57, t39, 0.2e1 * t14 * t44 - 0.2e1 * t41 * t65, -0.2e1 * t15 * t44 - 0.2e1 * t41 * t62, t19 ^ 2, t19 * t78, 0.2e1 * t19 * t44, t44 * t78, t39, 0.2e1 * t21 * t17 - 0.2e1 * t44 * t54, 0.2e1 * t21 * t19 - 0.2e1 * t4 * t44, 0.2e1 * t5 * t17 - 0.2e1 * t2 * t44, -0.2e1 * t1 * t17 + 0.2e1 * t2 * t19, 0.2e1 * t1 * t44 - 0.2e1 * t5 * t19, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t43, t58 * t46, 0, 0, 0, 0, 0, t52, -t51, t52, t16 * t19 - t18 * t17, t51, t1 * t18 + t2 * t16 - t5 * t47; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 + t18 ^ 2 + t41; 0, 0, 0, 0, 0, 0, 0, 0, t47, -t44, 0, t59, -t64, t43 * t32 (-t38 + t40) * t47, t68, t63, 0, t43 * t53 + t46 * t59, -t43 * t59 + t46 * t53, t19 * t23, -t23 * t17 - t19 * t22, t16, -t22 * t44, 0, t35 * t17 + t21 * t22 - t70, t35 * t19 + t21 * t23 - t69, t8 * t17 + t5 * t22 - t70, -t1 * t22 + t12 * t19 - t13 * t17 + t2 * t23, -t8 * t19 - t5 * t23 + t69, t1 * t13 + t2 * t12 + t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t44, 0, 0, 0, 0, 0, t32, -t66, 0, 0, 0, 0, 0, -t61, -t17, -t61, t16 * t23 - t18 * t22, t17, t16 * t12 + t18 * t13 - t47 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, 0.2e1 * t67, 0, 0, 0, pkin(3) * t75, -0.2e1 * pkin(3) * t43, t23 ^ 2, t22 * t77, 0, 0, 0, t22 * t76, t23 * t76, 0.2e1 * t8 * t22, 0.2e1 * t12 * t23 - 0.2e1 * t13 * t22, t8 * t77, t12 ^ 2 + t13 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t66, t44, t14, -t15, 0, 0, t19, -t17, t44, t44 * t71 - t54, -t44 * t37 - t4 (pkin(5) + t33) * t44 - t54, -t30 * t17 - t33 * t19, t30 * t44 + t1, t1 * t30 - t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t63, 0, 0, 0, 0, 0, -t16, -t18, -t16, 0, t18, -t16 * t33 + t18 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t46, 0, -t43 * pkin(8), -t46 * pkin(8), 0, 0, t23, -t22, 0, -t12, -t13, -t12, -t30 * t22 - t33 * t23, t13, -t12 * t33 + t13 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t37, 0.2e1 * t33, 0, 0.2e1 * t30, t30 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, t44, -t54, -t4, -t54 + 0.2e1 * t72, -pkin(5) * t19 - t17 * qJ(6), 0.2e1 * t36 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t18, -t16, 0, t18, -t16 * pkin(5) + t18 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t12, -t13, -t12, -pkin(5) * t23 - t22 * qJ(6), t13, -t12 * pkin(5) + t13 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t71, -t37, t50 + t71, 0, t49 + t37, t33 * pkin(5) + t30 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t50, 0, t49, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
