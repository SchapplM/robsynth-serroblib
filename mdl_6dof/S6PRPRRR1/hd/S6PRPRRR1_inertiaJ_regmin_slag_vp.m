% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = cos(pkin(12));
t31 = -t39 * pkin(2) - pkin(3);
t46 = cos(qJ(4));
t27 = -t46 * pkin(4) + t31;
t67 = 0.2e1 * t27;
t43 = sin(qJ(4));
t66 = 0.2e1 * t43;
t42 = sin(qJ(5));
t65 = t42 * pkin(4);
t45 = cos(qJ(6));
t37 = sin(pkin(12));
t38 = sin(pkin(6));
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t18 = (t37 * t47 + t39 * t44) * t38;
t40 = cos(pkin(6));
t14 = -t18 * t43 + t40 * t46;
t15 = t18 * t46 + t40 * t43;
t61 = cos(qJ(5));
t6 = -t61 * t14 + t42 * t15;
t64 = t6 * t45;
t52 = t61 * pkin(4);
t34 = -t52 - pkin(5);
t63 = pkin(5) - t34;
t30 = t37 * pkin(2) + pkin(8);
t62 = pkin(9) + t30;
t23 = t62 * t46;
t51 = t61 * t43;
t10 = t42 * t23 + t62 * t51;
t60 = t10 * t45;
t59 = t38 * t44;
t58 = t38 * t47;
t26 = t42 * t46 + t51;
t41 = sin(qJ(6));
t57 = t41 * t26;
t56 = t41 * t45;
t55 = t42 * t43;
t54 = t45 * t26;
t25 = -t61 * t46 + t55;
t53 = -0.2e1 * t26 * t25;
t50 = -pkin(5) * t26 - pkin(10) * t25;
t33 = pkin(10) + t65;
t49 = -t25 * t33 + t26 * t34;
t36 = t45 ^ 2;
t35 = t41 ^ 2;
t29 = 0.2e1 * t56;
t24 = t26 ^ 2;
t22 = t45 * t25;
t21 = t41 * t25;
t19 = t41 * t54;
t16 = t37 * t59 - t39 * t58;
t13 = (-t35 + t36) * t26;
t11 = t61 * t23 - t62 * t55;
t9 = t25 * pkin(5) - t26 * pkin(10) + t27;
t8 = t10 * t41;
t7 = t42 * t14 + t61 * t15;
t5 = t6 * t41;
t4 = t45 * t11 + t41 * t9;
t3 = -t41 * t11 + t45 * t9;
t2 = t16 * t41 + t45 * t7;
t1 = t16 * t45 - t41 * t7;
t12 = [1, 0, 0, 0, t16 ^ 2 + t18 ^ 2 + t40 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t58, -t59 (-t16 * t39 + t18 * t37) * pkin(2), 0, 0, 0, 0, 0, -t16 * t46, t16 * t43, 0, 0, 0, 0, 0, t16 * t25, t16 * t26, 0, 0, 0, 0, 0, t1 * t25 + t6 * t57, -t2 * t25 + t6 * t54; 0, 1, 0, 0 (t37 ^ 2 + t39 ^ 2) * pkin(2) ^ 2, t43 ^ 2, t46 * t66, 0, 0, 0, -0.2e1 * t31 * t46, t31 * t66, t24, t53, 0, 0, 0, t25 * t67, t26 * t67, t36 * t24, -0.2e1 * t24 * t56, 0.2e1 * t25 * t54, t41 * t53, t25 ^ 2, 0.2e1 * t10 * t57 + 0.2e1 * t3 * t25, 0.2e1 * t10 * t54 - 0.2e1 * t4 * t25; 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t64, t5; 0, 0, 0, 0, 0, 0, 0, t43, t46, 0, -t43 * t30, -t46 * t30, 0, 0, t26, -t25, 0, -t10, -t11, t19, t13, t21, t22, 0, t49 * t41 - t60, t49 * t45 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43, 0, 0, 0, 0, 0, -t25, -t26, 0, 0, 0, 0, 0, -t22, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t65, t35, t29, 0, 0, 0, -0.2e1 * t34 * t45, 0.2e1 * t34 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t64, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, -t10, -t11, t19, t13, t21, t22, 0, t50 * t41 - t60, t50 * t45 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, 0, 0, 0, 0, -t22, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t52, -t65, t35, t29, 0, 0, 0, t63 * t45, -t63 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t35, t29, 0, 0, 0, 0.2e1 * pkin(5) * t45, -0.2e1 * pkin(5) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t57, t25, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t45, 0, -t41 * t33, -t45 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t45, 0, -t41 * pkin(10), -t45 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t12;
