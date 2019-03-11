% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t44 = sin(qJ(4));
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t62 = cos(qJ(4));
t26 = t44 * t47 + t62 * t45;
t39 = sin(pkin(11));
t41 = cos(pkin(11));
t59 = t44 * t45;
t51 = -t62 * t47 + t59;
t17 = t39 * t26 + t41 * t51;
t70 = t17 ^ 2;
t42 = cos(pkin(10));
t34 = -t42 * pkin(1) - pkin(2);
t27 = -t47 * pkin(3) + t34;
t69 = 0.2e1 * t27;
t43 = sin(qJ(6));
t68 = 0.2e1 * t43;
t67 = 0.2e1 * t45;
t46 = cos(qJ(6));
t66 = -0.2e1 * t46;
t40 = sin(pkin(10));
t32 = t40 * pkin(1) + pkin(7);
t63 = pkin(8) + t32;
t54 = t62 * t63;
t10 = (t62 * qJ(5) + t54) * t47 + (-qJ(5) - t63) * t59;
t14 = t26 * t63;
t50 = -t26 * qJ(5) - t14;
t4 = t39 * t10 - t41 * t50;
t65 = t4 * t46;
t64 = t44 * pkin(3);
t12 = t43 * t17;
t19 = t41 * t26 - t39 * t51;
t61 = t43 * t19;
t60 = t43 * t46;
t57 = t46 * t19;
t36 = t62 * pkin(3);
t35 = t36 + pkin(4);
t23 = t41 * t35 - t39 * t64;
t21 = -pkin(5) - t23;
t33 = -t41 * pkin(4) - pkin(5);
t56 = t21 + t33;
t24 = t39 * t35 + t41 * t64;
t22 = pkin(9) + t24;
t53 = -t17 * t22 + t19 * t21;
t31 = t39 * pkin(4) + pkin(9);
t52 = -t17 * t31 + t19 * t33;
t20 = t51 * pkin(4) + t27;
t38 = t46 ^ 2;
t37 = t43 ^ 2;
t30 = 0.2e1 * t60;
t16 = t19 ^ 2;
t15 = -t47 * t54 + t63 * t59;
t13 = t46 * t17;
t11 = t43 * t57;
t8 = (-t37 + t38) * t19;
t7 = t17 * pkin(5) - t19 * pkin(9) + t20;
t6 = t41 * t10 + t39 * t50;
t3 = t4 * t43;
t2 = t43 * t7 + t46 * t6;
t1 = -t43 * t6 + t46 * t7;
t5 = [1, 0, 0 (t40 ^ 2 + t42 ^ 2) * pkin(1) ^ 2, t45 ^ 2, t47 * t67, 0, 0, 0, -0.2e1 * t34 * t47, t34 * t67, t26 ^ 2, -0.2e1 * t26 * t51, 0, 0, 0, t51 * t69, t26 * t69, -0.2e1 * t6 * t17 + 0.2e1 * t4 * t19, t20 ^ 2 + t4 ^ 2 + t6 ^ 2, t38 * t16, -0.2e1 * t16 * t60, 0.2e1 * t17 * t57, -0.2e1 * t17 * t61, t70, 0.2e1 * t1 * t17 + 0.2e1 * t4 * t61, -0.2e1 * t2 * t17 + 0.2e1 * t4 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t17 + t6 * t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 + t70, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t45, t47, 0, -t45 * t32, -t47 * t32, 0, 0, t26, -t51, 0, -t14, t15, -t24 * t17 - t23 * t19, -t4 * t23 + t6 * t24, t11, t8, t12, t13, 0, t53 * t43 - t65, t53 * t46 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t45, 0, 0, 0, 0, 0, -t51, -t26, 0, -t17 * t23 + t19 * t24, 0, 0, 0, 0, 0, -t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t64, 0, t23 ^ 2 + t24 ^ 2, t37, t30, 0, 0, 0, t21 * t66, t21 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t51, 0, -t14, t15 (-t17 * t39 - t19 * t41) * pkin(4) (t39 * t6 - t4 * t41) * pkin(4), t11, t8, t12, t13, 0, t52 * t43 - t65, t52 * t46 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t26, 0 (-t17 * t41 + t19 * t39) * pkin(4), 0, 0, 0, 0, 0, -t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t64, 0 (t23 * t41 + t24 * t39) * pkin(4), t37, t30, 0, 0, 0, -t56 * t46, t56 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t39 ^ 2 + t41 ^ 2) * pkin(4) ^ 2, t37, t30, 0, 0, 0, t33 * t66, t33 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, t13, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t61, t17, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t46, 0, -t43 * t22, -t46 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t46, 0, -t43 * t31, -t46 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
