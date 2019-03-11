% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(4));
t67 = -0.2e1 * t41;
t38 = cos(pkin(11));
t30 = -t38 * pkin(5) - pkin(4);
t68 = 0.2e1 * t30;
t66 = 2 * qJ(3);
t44 = cos(qJ(4));
t65 = t44 * pkin(4);
t36 = sin(pkin(11));
t40 = sin(qJ(6));
t43 = cos(qJ(6));
t22 = t40 * t36 - t43 * t38;
t64 = t22 * t41;
t23 = t43 * t36 + t40 * t38;
t63 = t23 * t41;
t35 = t44 ^ 2;
t46 = -pkin(2) - pkin(8);
t62 = t35 * t46;
t61 = t36 * t44;
t60 = t36 * t46;
t37 = sin(pkin(6));
t42 = sin(qJ(2));
t59 = t37 * t42;
t45 = cos(qJ(2));
t58 = t37 * t45;
t28 = t38 * t44;
t57 = t41 * t46;
t56 = t44 * t22;
t14 = t44 * t23;
t55 = t44 * t46;
t54 = pkin(9) + qJ(5);
t24 = t41 * pkin(4) - t44 * qJ(5) + qJ(3);
t12 = t36 * t24 + t38 * t57;
t53 = t36 ^ 2 + t38 ^ 2;
t34 = t41 ^ 2;
t52 = -t34 - t35;
t51 = t53 * qJ(5);
t39 = cos(pkin(6));
t18 = t39 * t44 - t41 * t58;
t7 = -t18 * t36 + t38 * t59;
t8 = t18 * t38 + t36 * t59;
t50 = -t7 * t36 + t8 * t38;
t49 = -qJ(5) * t41 - t65;
t20 = t38 * t24;
t11 = -t36 * t57 + t20;
t48 = -t11 * t36 + t12 * t38;
t26 = t54 * t38;
t25 = t54 * t36;
t21 = (pkin(5) * t36 - t46) * t44;
t17 = t39 * t41 + t44 * t58;
t10 = -t40 * t25 + t43 * t26;
t9 = -t43 * t25 - t40 * t26;
t6 = -pkin(9) * t61 + t12;
t5 = -pkin(9) * t28 + t20 + (pkin(5) - t60) * t41;
t4 = t40 * t7 + t43 * t8;
t3 = -t40 * t8 + t43 * t7;
t2 = t40 * t5 + t43 * t6;
t1 = -t40 * t6 + t43 * t5;
t13 = [1, 0, 0, 0, 0, 0, t39 ^ 2 + (t42 ^ 2 + t45 ^ 2) * t37 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 + t7 ^ 2 + t8 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t58, -t59, -t58, t59 (pkin(2) * t45 + qJ(3) * t42) * t37, 0, 0, 0, 0, 0, t41 * t59, t44 * t59, t17 * t61 + t7 * t41, t17 * t28 - t8 * t41 (-t36 * t8 - t38 * t7) * t44, t7 * t11 + t8 * t12 - t17 * t55, 0, 0, 0, 0, 0, t17 * t14 + t3 * t41, -t17 * t56 - t4 * t41; 0, 1, 0, 0, -0.2e1 * pkin(2), t66, pkin(2) ^ 2 + (qJ(3) ^ 2) t35, t44 * t67, 0, 0, 0, t41 * t66, t44 * t66, 0.2e1 * t11 * t41 - 0.2e1 * t35 * t60, -0.2e1 * t12 * t41 - 0.2e1 * t38 * t62, 0.2e1 * (-t11 * t38 - t12 * t36) * t44, t35 * t46 ^ 2 + t11 ^ 2 + t12 ^ 2, t56 ^ 2, 0.2e1 * t56 * t14, t56 * t67, t14 * t67, t34, 0.2e1 * t1 * t41 + 0.2e1 * t21 * t14, -0.2e1 * t2 * t41 - 0.2e1 * t21 * t56; 0, 0, 0, 0, 0, 0, -t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t44 + t50 * t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, t52 * t36, t52 * t38, 0, t48 * t41 + t62, 0, 0, 0, 0, 0, -t44 * t14 - t41 * t63, t41 * t64 + t44 * t56; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t34 + t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, -t17 * t38, t17 * t36, t50, -t17 * pkin(4) + t50 * qJ(5), 0, 0, 0, 0, 0, t17 * t22, t17 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, t55, -t57, t49 * t36 + t38 * t55, -t36 * t55 + t49 * t38, t48, pkin(4) * t55 + t48 * qJ(5), -t56 * t23, -t23 * t14 + t22 * t56, t63, -t64, 0, t30 * t14 + t21 * t22 + t9 * t41, -t10 * t41 + t21 * t23 - t30 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, t28, -t61, t53 * t41, t41 * t51 + t65, 0, 0, 0, 0, 0, -t56, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t38, -0.2e1 * pkin(4) * t36, 0.2e1 * t51, t53 * qJ(5) ^ 2 + pkin(4) ^ 2, t23 ^ 2, -0.2e1 * t23 * t22, 0, 0, 0, t22 * t68, t23 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t28, 0, -t55, 0, 0, 0, 0, 0, t14, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t36, 0, -pkin(4), 0, 0, 0, 0, 0, t22, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t14, t41, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
