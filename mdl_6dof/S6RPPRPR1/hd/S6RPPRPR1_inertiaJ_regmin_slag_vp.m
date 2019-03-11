% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t43 = sin(pkin(10));
t46 = cos(pkin(10));
t49 = sin(qJ(4));
t64 = cos(qJ(4));
t26 = t43 * t49 - t46 * t64;
t69 = -0.2e1 * t26;
t47 = cos(pkin(9));
t36 = -pkin(1) * t47 - pkin(2);
t29 = -pkin(3) * t46 + t36;
t68 = 0.2e1 * t29;
t45 = cos(pkin(11));
t37 = -pkin(5) * t45 - pkin(4);
t67 = 0.2e1 * t37;
t66 = t26 * pkin(4);
t44 = sin(pkin(9));
t33 = pkin(1) * t44 + qJ(3);
t65 = pkin(7) + t33;
t28 = t43 * t64 + t46 * t49;
t11 = -qJ(5) * t28 + t29 + t66;
t21 = t65 * t43;
t22 = t65 * t46;
t15 = -t21 * t49 + t22 * t64;
t42 = sin(pkin(11));
t6 = t11 * t42 + t15 * t45;
t48 = sin(qJ(6));
t50 = cos(qJ(6));
t25 = t42 * t48 - t45 * t50;
t19 = t26 * t25;
t63 = t26 * t42;
t27 = t42 * t50 + t45 * t48;
t20 = t27 * t26;
t62 = t42 * t28;
t61 = t45 * t26;
t60 = t45 * t28;
t59 = pkin(8) + qJ(5);
t58 = t42 ^ 2 + t45 ^ 2;
t57 = t43 ^ 2 + t46 ^ 2;
t5 = t11 * t45 - t15 * t42;
t56 = t58 * qJ(5);
t55 = t42 * t6 + t45 * t5;
t54 = -t42 * t5 + t45 * t6;
t53 = -pkin(4) * t28 - qJ(5) * t26;
t14 = t21 * t64 + t22 * t49;
t31 = t59 * t45;
t30 = t59 * t42;
t24 = t28 ^ 2;
t23 = t26 ^ 2;
t18 = -t30 * t48 + t31 * t50;
t17 = -t30 * t50 - t31 * t48;
t16 = t58 * t28;
t13 = t25 * t28;
t12 = t27 * t28;
t7 = pkin(5) * t62 + t14;
t4 = -pkin(8) * t62 + t6;
t3 = pkin(5) * t26 - pkin(8) * t60 + t5;
t2 = t3 * t48 + t4 * t50;
t1 = t3 * t50 - t4 * t48;
t8 = [1, 0, 0 (t44 ^ 2 + t47 ^ 2) * pkin(1) ^ 2, -0.2e1 * t36 * t46, 0.2e1 * t36 * t43, 0.2e1 * t57 * t33, t33 ^ 2 * t57 + t36 ^ 2, t24, t28 * t69, 0, 0, 0, t26 * t68, t28 * t68, 0.2e1 * t14 * t62 + 0.2e1 * t26 * t5, 0.2e1 * t14 * t60 - 0.2e1 * t26 * t6, -0.2e1 * t55 * t28, t14 ^ 2 + t5 ^ 2 + t6 ^ 2, t13 ^ 2, 0.2e1 * t13 * t12, t13 * t69, t12 * t69, t23, 0.2e1 * t1 * t26 + 0.2e1 * t12 * t7, -0.2e1 * t13 * t7 - 0.2e1 * t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t26 + t28 * t54, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t58 + t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t46, t43, 0, t36, 0, 0, 0, 0, 0, t26, t28, t61, -t63, -t16, t55, 0, 0, 0, 0, 0, -t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, 0, -t14, -t15, -t14 * t45 + t42 * t53, t14 * t42 + t45 * t53, t54, -t14 * pkin(4) + qJ(5) * t54, -t13 * t27, -t12 * t27 + t13 * t25, t20, -t19, 0, t12 * t37 + t17 * t26 + t25 * t7, -t13 * t37 - t18 * t26 + t27 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t28, -t61, t63, t16, t28 * t56 - t66, 0, 0, 0, 0, 0, t19, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t45, -0.2e1 * pkin(4) * t42, 0.2e1 * t56, qJ(5) ^ 2 * t58 + pkin(4) ^ 2, t27 ^ 2, -0.2e1 * t27 * t25, 0, 0, 0, t25 * t67, t27 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t60, 0, t14, 0, 0, 0, 0, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t42, 0, -pkin(4), 0, 0, 0, 0, 0, t25, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, t26, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
