% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t48 = cos(qJ(4));
t36 = -t48 * pkin(4) - pkin(3);
t71 = 0.2e1 * t36;
t49 = cos(qJ(3));
t70 = -0.2e1 * t49;
t69 = pkin(8) + pkin(9);
t68 = pkin(3) * t48;
t44 = sin(qJ(5));
t45 = sin(qJ(4));
t47 = cos(qJ(5));
t25 = t44 * t48 + t47 * t45;
t46 = sin(qJ(3));
t18 = t25 * t46;
t67 = t18 * pkin(5);
t66 = t44 * pkin(4);
t43 = cos(pkin(10));
t33 = -t43 * pkin(1) - pkin(2);
t23 = -t49 * pkin(3) - t46 * pkin(8) + t33;
t42 = sin(pkin(10));
t32 = t42 * pkin(1) + pkin(7);
t54 = t49 * t32;
t51 = t48 * t54;
t9 = t51 + (-pkin(9) * t46 + t23) * t45;
t65 = t47 * t9;
t64 = t49 * pkin(4);
t63 = t25 * t18;
t62 = t25 * t49;
t61 = t32 * t45;
t60 = t45 * t46;
t59 = t45 * t48;
t58 = t45 * t49;
t57 = t48 * t46;
t56 = t48 * t49;
t24 = t44 * t45 - t47 * t48;
t55 = t49 * t24;
t53 = t49 * t46;
t27 = t46 * t32;
t22 = pkin(4) * t60 + t27;
t52 = 0.2e1 * t53;
t21 = t48 * t23;
t6 = -pkin(9) * t57 + t21 + (-pkin(4) - t61) * t49;
t3 = -t44 * t9 + t47 * t6;
t28 = t69 * t45;
t29 = t69 * t48;
t14 = -t47 * t28 - t44 * t29;
t4 = t44 * t6 + t65;
t15 = -t44 * t28 + t47 * t29;
t41 = t49 ^ 2;
t40 = t48 ^ 2;
t39 = t46 ^ 2;
t38 = t45 ^ 2;
t37 = t47 * pkin(4);
t35 = t37 + pkin(5);
t20 = -t44 * t60 + t47 * t57;
t17 = t20 ^ 2;
t16 = t24 * pkin(5) + t36;
t13 = t45 * t23 + t51;
t12 = -t45 * t54 + t21;
t11 = t20 * t24;
t10 = t22 + t67;
t8 = -t24 * qJ(6) + t15;
t7 = -t25 * qJ(6) + t14;
t2 = -t18 * qJ(6) + t4;
t1 = -t49 * pkin(5) - t20 * qJ(6) + t3;
t5 = [1, 0, 0 (t42 ^ 2 + t43 ^ 2) * pkin(1) ^ 2, t39, t52, 0, 0, 0, t33 * t70, 0.2e1 * t33 * t46, t40 * t39, -0.2e1 * t39 * t59, -0.2e1 * t48 * t53, t45 * t52, t41, -0.2e1 * t12 * t49 + 0.2e1 * t39 * t61, 0.2e1 * t39 * t32 * t48 + 0.2e1 * t13 * t49, t17, -0.2e1 * t20 * t18, t20 * t70, -t18 * t70, t41, 0.2e1 * t22 * t18 - 0.2e1 * t3 * t49, 0.2e1 * t22 * t20 + 0.2e1 * t4 * t49, -0.2e1 * t1 * t20 - 0.2e1 * t2 * t18, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t18 - t10 * t49 + t2 * t20; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 ^ 2 + t17 + t41; 0, 0, 0, 0, 0, 0, t46, t49, 0, -t27, -t54, t45 * t57 (-t38 + t40) * t46, -t58, -t56, 0, -t32 * t57 + (-pkin(3) * t46 + pkin(8) * t49) * t45, pkin(8) * t56 + (t61 - t68) * t46, t20 * t25, -t11 - t63, -t62, t55, 0, -t14 * t49 + t36 * t18 + t22 * t24, t15 * t49 + t36 * t20 + t22 * t25, -t1 * t25 - t8 * t18 - t2 * t24 - t7 * t20, t1 * t7 + t10 * t16 + t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46, 0, 0, 0, 0, 0, t56, -t58, 0, 0, 0, 0, 0, -t55, -t62, -t11 + t63, -t49 * t16 - t18 * t7 + t20 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, 0.2e1 * t59, 0, 0, 0, 0.2e1 * t68, -0.2e1 * pkin(3) * t45, t25 ^ 2, -0.2e1 * t25 * t24, 0, 0, 0, t24 * t71, t25 * t71, -0.2e1 * t8 * t24 - 0.2e1 * t7 * t25, t16 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t60, -t49, t12, -t13, 0, 0, t20, -t18, -t49, -t47 * t64 + t3, -t65 + (-t6 + t64) * t44, -t18 * t66 - t35 * t20, t1 * t35 + t2 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t57, 0, 0, 0, 0, 0, -t18, -t20, 0, -t18 * t35 + t20 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t48, 0, -t45 * pkin(8), -t48 * pkin(8), 0, 0, t25, -t24, 0, t14, -t15, -t24 * t66 - t35 * t25, t7 * t35 + t8 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t66, 0, t44 ^ 2 * pkin(4) ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, -t49, t3, -t4, -pkin(5) * t20, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t20, 0, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, t14, -t15, -pkin(5) * t25, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t66, 0, t35 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;
