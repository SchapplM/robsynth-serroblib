% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t71 = 2 * pkin(5);
t35 = sin(qJ(5));
t70 = -0.2e1 * t35;
t36 = sin(qJ(4));
t69 = -0.2e1 * t36;
t38 = cos(qJ(4));
t68 = 0.2e1 * t38;
t37 = cos(qJ(5));
t67 = pkin(4) * t37;
t66 = pkin(5) * t35;
t65 = pkin(8) * t38;
t64 = t35 * pkin(8);
t63 = t37 * pkin(8);
t33 = sin(pkin(9));
t34 = cos(pkin(9));
t39 = -pkin(1) - pkin(2);
t17 = t34 * qJ(2) + t33 * t39;
t14 = -pkin(7) + t17;
t55 = t38 * t14;
t15 = t33 * qJ(2) - t34 * t39;
t13 = pkin(3) + t15;
t9 = t38 * pkin(4) + t36 * pkin(8) + t13;
t4 = t35 * t9 + t37 * t55;
t62 = t14 * t35;
t30 = t36 ^ 2;
t61 = t30 * t35;
t60 = t30 * t37;
t31 = t37 ^ 2;
t59 = t31 * t36;
t23 = t35 * t36;
t58 = t35 * t37;
t57 = t36 * t33;
t56 = t37 * t36;
t26 = t37 * t38;
t54 = t38 * t33;
t29 = t35 ^ 2;
t53 = t29 + t31;
t52 = t37 * qJ(6);
t51 = t38 * qJ(6);
t50 = t36 * t68;
t49 = t35 * t57;
t48 = t33 * t56;
t47 = t53 * pkin(8);
t1 = t51 + t4;
t8 = t37 * t9;
t2 = -t8 + (-pkin(5) + t62) * t38;
t46 = t1 * t37 + t2 * t35;
t44 = t37 * pkin(5) + t35 * qJ(6);
t18 = -pkin(4) - t44;
t45 = t18 * t36 + t65;
t43 = t52 - t66;
t11 = t37 * t34 + t35 * t54;
t12 = -t35 * t34 + t37 * t54;
t42 = t11 * t35 + t12 * t37;
t41 = -t11 * t38 - t33 * t61;
t32 = t38 ^ 2;
t28 = t33 ^ 2;
t25 = t31 * t30;
t24 = t35 * t38;
t22 = t29 * t36;
t19 = t36 * t52;
t6 = t12 * t38 + t33 * t60;
t5 = t19 + (t14 - t66) * t36;
t3 = -t35 * t55 + t8;
t7 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t15, 0.2e1 * t17, t15 ^ 2 + t17 ^ 2, t30, t50, 0, 0, 0, t13 * t68, t13 * t69, t25, -0.2e1 * t30 * t58, t26 * t69, t35 * t50, t32, -0.2e1 * t14 * t61 + 0.2e1 * t3 * t38, -0.2e1 * t14 * t60 - 0.2e1 * t4 * t38, -0.2e1 * t2 * t38 - 0.2e1 * t5 * t23, 0.2e1 * (t1 * t35 - t2 * t37) * t36, 0.2e1 * t1 * t38 + 0.2e1 * t5 * t56, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -1, 0, -pkin(1), -t34, t33, -t15 * t34 + t17 * t33, 0, 0, 0, 0, 0, -t34 * t38, t36 * t34, 0, 0, 0, 0, 0, t41, -t6, t41 (-t11 * t37 + t12 * t35) * t36, t6, t1 * t12 + t2 * t11 + t5 * t57; 0, 0, 0, 0, 0, 1, 0, 0, t34 ^ 2 + t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 + t12 ^ 2 + t30 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t36 - t5 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t42 - t54) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t30 + t25 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t38, 0, -t36 * t14, -t55, -t35 * t56, t22 - t59, t24, t26, 0, -t14 * t56 + (pkin(4) * t36 - t65) * t35, -pkin(8) * t26 + (t62 + t67) * t36, -t45 * t35 - t5 * t37, t46, -t5 * t35 + t45 * t37, t46 * pkin(8) + t5 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t54, 0, 0, 0, 0, 0, -t48, t49, -t48, t42, -t49, t42 * pkin(8) + t18 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0, 0, 0, 0, 0, t26, -t24, t26, t22 + t59, t24, -t38 * t18 + t36 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t29, 0.2e1 * t58, 0, 0, 0, 0.2e1 * t67, pkin(4) * t70, -0.2e1 * t18 * t37, 0.2e1 * t47, t18 * t70, t53 * pkin(8) ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t23, t38, t3, -t4, t8 + (t71 - t62) * t38, t44 * t36, 0.2e1 * t51 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -t11, 0, t12, -t11 * pkin(5) + t12 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t56, -t23, 0, t56, -pkin(5) * t23 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t37, 0, -t64, -t63, -t64, t43, t63, t43 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t71, 0, 0.2e1 * qJ(6) (pkin(5) ^ 2) + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t56, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t7;
