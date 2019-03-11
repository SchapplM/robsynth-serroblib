% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRR2
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t47 = cos(qJ(5));
t33 = -pkin(5) * t47 - pkin(4);
t72 = 0.2e1 * t33;
t48 = cos(qJ(4));
t71 = -0.2e1 * t48;
t70 = 0.2e1 * t48;
t69 = pkin(9) + pkin(10);
t68 = pkin(4) * t47;
t42 = sin(qJ(6));
t67 = t42 * pkin(5);
t46 = cos(qJ(6));
t66 = t46 * pkin(5);
t40 = cos(pkin(12));
t31 = -pkin(2) * t40 - pkin(3);
t44 = sin(qJ(4));
t23 = -pkin(4) * t48 - pkin(9) * t44 + t31;
t43 = sin(qJ(5));
t38 = sin(pkin(12));
t30 = pkin(2) * t38 + pkin(8);
t53 = t48 * t30;
t51 = t47 * t53;
t9 = t51 + (-pkin(10) * t44 + t23) * t43;
t65 = t46 * t9;
t64 = t48 * pkin(5);
t25 = t42 * t47 + t43 * t46;
t63 = t25 * t48;
t62 = t30 * t43;
t39 = sin(pkin(6));
t45 = sin(qJ(2));
t61 = t39 * t45;
t49 = cos(qJ(2));
t60 = t39 * t49;
t59 = t43 * t44;
t58 = t43 * t47;
t57 = t43 * t48;
t56 = t47 * t44;
t55 = t47 * t48;
t24 = t42 * t43 - t46 * t47;
t54 = t48 * t24;
t52 = t44 * t70;
t21 = t47 * t23;
t8 = -pkin(10) * t56 + t21 + (-pkin(5) - t62) * t48;
t3 = -t42 * t9 + t46 * t8;
t41 = cos(pkin(6));
t37 = t48 ^ 2;
t36 = t47 ^ 2;
t35 = t44 ^ 2;
t34 = t43 ^ 2;
t28 = t69 * t47;
t27 = t69 * t43;
t22 = (pkin(5) * t43 + t30) * t44;
t20 = -t42 * t59 + t46 * t56;
t19 = t25 * t44;
t18 = (t38 * t49 + t40 * t45) * t39;
t16 = t38 * t61 - t40 * t60;
t15 = -t27 * t42 + t28 * t46;
t14 = -t27 * t46 - t28 * t42;
t13 = t18 * t48 + t41 * t44;
t12 = t18 * t44 - t41 * t48;
t11 = t23 * t43 + t51;
t10 = -t43 * t53 + t21;
t6 = t13 * t47 + t16 * t43;
t5 = -t13 * t43 + t16 * t47;
t4 = t42 * t8 + t65;
t2 = t42 * t5 + t46 * t6;
t1 = -t42 * t6 + t46 * t5;
t7 = [1, 0, 0, 0, t16 ^ 2 + t18 ^ 2 + t41 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t60, -t61 (-t16 * t40 + t18 * t38) * pkin(2), 0, 0, 0, 0, 0, -t16 * t48, t16 * t44, 0, 0, 0, 0, 0, t12 * t59 - t48 * t5, t12 * t56 + t48 * t6, 0, 0, 0, 0, 0, -t1 * t48 + t12 * t19, t12 * t20 + t2 * t48; 0, 1, 0, 0 (t38 ^ 2 + t40 ^ 2) * pkin(2) ^ 2, t35, t52, 0, 0, 0, t31 * t71, 0.2e1 * t31 * t44, t36 * t35, -0.2e1 * t35 * t58, -0.2e1 * t44 * t55, t43 * t52, t37, -0.2e1 * t10 * t48 + 0.2e1 * t35 * t62, 0.2e1 * t30 * t35 * t47 + 0.2e1 * t11 * t48, t20 ^ 2, -0.2e1 * t20 * t19, t20 * t71, t19 * t70, t37, 0.2e1 * t19 * t22 - 0.2e1 * t3 * t48, 0.2e1 * t20 * t22 + 0.2e1 * t4 * t48; 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t12 * t47, t12 * t43, 0, 0, 0, 0, 0, t12 * t24, t12 * t25; 0, 0, 0, 0, 0, 0, 0, t44, t48, 0, -t44 * t30, -t53, t43 * t56 (-t34 + t36) * t44, -t57, -t55, 0, -t30 * t56 + (-pkin(4) * t44 + pkin(9) * t48) * t43, pkin(9) * t55 + (t62 - t68) * t44, t20 * t25, -t19 * t25 - t20 * t24, -t63, t54, 0, -t14 * t48 + t19 * t33 + t22 * t24, t15 * t48 + t20 * t33 + t22 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t44, 0, 0, 0, 0, 0, t55, -t57, 0, 0, 0, 0, 0, -t54, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t34, 0.2e1 * t58, 0, 0, 0, 0.2e1 * t68, -0.2e1 * pkin(4) * t43, t25 ^ 2, -0.2e1 * t25 * t24, 0, 0, 0, t24 * t72, t25 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t59, -t48, t10, -t11, 0, 0, t20, -t19, -t48, -t46 * t64 + t3, -t65 + (-t8 + t64) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t56, 0, 0, 0, 0, 0, -t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t47, 0, -t43 * pkin(9), -t47 * pkin(9), 0, 0, t25, -t24, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t66, -0.2e1 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, -t48, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t66, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t7;
