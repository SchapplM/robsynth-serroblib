% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR1
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
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t43 = cos(pkin(11));
t33 = -t43 * pkin(1) - pkin(2);
t50 = cos(qJ(3));
t28 = -t50 * pkin(3) + t33;
t46 = sin(qJ(4));
t47 = sin(qJ(3));
t49 = cos(qJ(4));
t53 = t46 * t47 - t49 * t50;
t18 = t53 * pkin(4) + t28;
t73 = 0.2e1 * t18;
t72 = 0.2e1 * t28;
t71 = 0.2e1 * t47;
t44 = sin(qJ(6));
t70 = pkin(5) * t44;
t45 = sin(qJ(5));
t42 = sin(pkin(11));
t32 = t42 * pkin(1) + pkin(7);
t66 = pkin(8) + t32;
t25 = t66 * t47;
t26 = t66 * t50;
t11 = -t49 * t25 - t46 * t26;
t27 = t46 * t50 + t49 * t47;
t52 = -t27 * pkin(9) + t11;
t65 = cos(qJ(5));
t12 = t46 * t25 - t49 * t26;
t9 = -t53 * pkin(9) - t12;
t4 = t45 * t9 - t65 * t52;
t48 = cos(qJ(6));
t69 = t4 * t48;
t68 = t45 * pkin(4);
t67 = t46 * pkin(3);
t38 = t49 * pkin(3);
t36 = t38 + pkin(4);
t57 = -t65 * t36 + t45 * t67;
t20 = -pkin(5) + t57;
t64 = t20 * t48;
t37 = t65 * pkin(4);
t35 = -t37 - pkin(5);
t63 = t35 * t48;
t17 = t65 * t27 - t45 * t53;
t62 = t44 * t17;
t61 = t44 * t48;
t60 = t48 * t17;
t16 = t45 * t27 + t65 * t53;
t59 = -0.2e1 * t17 * t16;
t58 = t65 * t67;
t56 = -pkin(5) * t17 - pkin(10) * t16;
t23 = -t45 * t36 - t58;
t21 = pkin(10) - t23;
t55 = -t16 * t21 + t17 * t20;
t34 = pkin(10) + t68;
t54 = -t16 * t34 + t17 * t35;
t41 = t48 ^ 2;
t40 = t44 ^ 2;
t39 = pkin(5) * t48;
t31 = 0.2e1 * t61;
t30 = t35 * t44;
t19 = t20 * t44;
t15 = t17 ^ 2;
t14 = t48 * t16;
t13 = t44 * t16;
t10 = t44 * t60;
t7 = (-t40 + t41) * t17;
t6 = t16 * pkin(5) - t17 * pkin(10) + t18;
t5 = t45 * t52 + t65 * t9;
t3 = t4 * t44;
t2 = t44 * t6 + t48 * t5;
t1 = -t44 * t5 + t48 * t6;
t8 = [1, 0, 0 (t42 ^ 2 + t43 ^ 2) * pkin(1) ^ 2, t47 ^ 2, t50 * t71, 0, 0, 0, -0.2e1 * t33 * t50, t33 * t71, t27 ^ 2, -0.2e1 * t27 * t53, 0, 0, 0, t53 * t72, t27 * t72, t15, t59, 0, 0, 0, t16 * t73, t17 * t73, t41 * t15, -0.2e1 * t15 * t61, 0.2e1 * t16 * t60, t44 * t59, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t4 * t62, -0.2e1 * t2 * t16 + 0.2e1 * t4 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t47, t50, 0, -t47 * t32, -t50 * t32, 0, 0, t27, -t53, 0, t11, t12, 0, 0, t17, -t16, 0, -t4, -t5, t10, t7, t13, t14, 0, t55 * t44 - t69, t55 * t48 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47, 0, 0, 0, 0, 0, -t53, -t27, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t67, 0, 0, 0, 0, 1, -0.2e1 * t57, 0.2e1 * t23, t40, t31, 0, 0, 0, -0.2e1 * t64, 0.2e1 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t53, 0, t11, t12, 0, 0, t17, -t16, 0, -t4, -t5, t10, t7, t13, t14, 0, t54 * t44 - t69, t54 * t48 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t27, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t67, 0, 0, 0, 0, 1, t37 - t57, -t58 + (-pkin(4) - t36) * t45, t40, t31, 0, 0, 0 (-t20 - t35) * t48, t30 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t68, t40, t31, 0, 0, 0, -0.2e1 * t63, 0.2e1 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t4, -t5, t10, t7, t13, t14, 0, t56 * t44 - t69, t56 * t48 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t57, t23, t40, t31, 0, 0, 0, t39 - t64, t19 - t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t68, t40, t31, 0, 0, 0, t39 - t63, t30 - t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t40, t31, 0, 0, 0, 0.2e1 * t39, -0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t62, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t48, 0, -t44 * t21, -t48 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t48, 0, -t44 * t34, -t48 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t48, 0, -t44 * pkin(10), -t48 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
