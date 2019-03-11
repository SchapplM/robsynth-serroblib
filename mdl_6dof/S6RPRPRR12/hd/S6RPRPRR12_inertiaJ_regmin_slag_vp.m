% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t46 = cos(qJ(3));
t41 = sin(qJ(6));
t42 = sin(qJ(5));
t44 = cos(qJ(6));
t45 = cos(qJ(5));
t50 = t41 * t42 - t44 * t45;
t67 = t50 * t46;
t29 = t42 * pkin(5) + qJ(4);
t66 = 0.2e1 * t29;
t65 = -0.2e1 * t46;
t64 = 2 * qJ(2);
t63 = 0.2e1 * qJ(4);
t62 = t41 * pkin(5);
t61 = t44 * pkin(5);
t43 = sin(qJ(3));
t21 = t43 * pkin(3) - t46 * qJ(4) + qJ(2);
t16 = t43 * pkin(8) + t21;
t51 = pkin(9) * t43 + t16;
t48 = -pkin(1) - pkin(7);
t25 = (pkin(4) - t48) * t46;
t58 = t42 * t25;
t5 = t51 * t45 + t58;
t60 = t44 * t5;
t59 = t46 * pkin(5);
t18 = t41 * t45 + t44 * t42;
t13 = t18 * t46;
t30 = t42 * t43;
t57 = t42 * t46;
t12 = t43 * t18;
t56 = t45 * t42;
t31 = t45 * t43;
t55 = t46 * t43;
t54 = t46 * t48;
t37 = t43 ^ 2;
t39 = t46 ^ 2;
t28 = t37 + t39;
t53 = t43 * qJ(4);
t52 = 0.2e1 * t55;
t17 = t45 * t25;
t4 = -t51 * t42 + t17 + t59;
t1 = t44 * t4 - t41 * t5;
t26 = t46 * pkin(3) + t53;
t47 = -pkin(3) - pkin(8);
t49 = -t46 * t47 + t53;
t38 = t45 ^ 2;
t36 = t42 ^ 2;
t34 = t45 * t47;
t33 = t43 * t48;
t32 = t45 * t46;
t24 = -t45 * pkin(9) + t34;
t23 = -t43 * pkin(4) + t33;
t22 = (-pkin(9) + t47) * t42;
t20 = t28 * t48;
t14 = t33 + (-pkin(5) * t45 - pkin(4)) * t43;
t10 = t41 * t30 - t44 * t31;
t9 = t44 * t22 + t41 * t24;
t8 = -t41 * t22 + t44 * t24;
t7 = t45 * t16 + t58;
t6 = -t42 * t16 + t17;
t2 = t41 * t4 + t60;
t3 = [1, 0, 0, -2 * pkin(1), t64, pkin(1) ^ 2 + qJ(2) ^ 2, t39, -0.2e1 * t55, 0, 0, 0, t43 * t64, t46 * t64, -0.2e1 * t20, -0.2e1 * t21 * t43, t21 * t65, t28 * t48 ^ 2 + t21 ^ 2, t36 * t37, 0.2e1 * t37 * t56, t42 * t52, t45 * t52, t39, -0.2e1 * t23 * t31 + 0.2e1 * t6 * t46, 0.2e1 * t23 * t30 - 0.2e1 * t7 * t46, t12 ^ 2, -0.2e1 * t12 * t10, 0.2e1 * t12 * t46, t10 * t65, t39, 0.2e1 * t1 * t46 + 0.2e1 * t14 * t10, 0.2e1 * t14 * t12 - 0.2e1 * t2 * t46; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t28, 0, 0, t20, 0, 0, 0, 0, 0, -t28 * t45, t28 * t42, 0, 0, 0, 0, 0, t43 * t10 + t46 * t67, t43 * t12 + t13 * t46; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43, 0, t54, -t33, -t26, -t54, t33, t26 * t48, t42 * t31 (-t36 + t38) * t43, t32, -t57, 0, t23 * t42 - t49 * t45, t23 * t45 + t49 * t42, -t12 * t50, t10 * t50 - t12 * t18, -t67, -t13, 0, t29 * t10 + t14 * t18 + t8 * t46, t29 * t12 - t14 * t50 - t9 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43, 0, -t46, t43, t26, 0, 0, 0, 0, 0, t30, t31, 0, 0, 0, 0, 0, t12, -t43 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t63, pkin(3) ^ 2 + qJ(4) ^ 2, t38, -0.2e1 * t56, 0, 0, 0, t42 * t63, t45 * t63, t50 ^ 2, 0.2e1 * t50 * t18, 0, 0, 0, t18 * t66, -t50 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, -t54, 0, 0, 0, 0, 0, t32, -t57, 0, 0, 0, 0, 0, -t67, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, t46, t6, -t7, 0, 0, t12, -t10, t46, t44 * t59 + t1, -t60 + (-t4 - t59) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t57, 0, 0, 0, 0, 0, t67, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, t34, -t42 * t47, 0, 0, -t50, -t18, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, 0, 0, 0, 0, -t50, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t61, -0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10, t46, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t18, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t61, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
