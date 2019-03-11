% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t42 = cos(pkin(12));
t34 = -pkin(5) * t42 - pkin(4);
t73 = 0.2e1 * t34;
t46 = sin(qJ(4));
t72 = 0.2e1 * t46;
t49 = cos(qJ(4));
t71 = -0.2e1 * t49;
t70 = 0.2e1 * t49;
t69 = t49 * pkin(4);
t39 = sin(pkin(12));
t45 = sin(qJ(6));
t48 = cos(qJ(6));
t25 = t39 * t48 + t42 * t45;
t68 = t25 * t49;
t40 = sin(pkin(11));
t32 = pkin(2) * t40 + pkin(8);
t67 = t32 * t39;
t66 = t39 * t46;
t41 = sin(pkin(6));
t47 = sin(qJ(2));
t65 = t41 * t47;
t50 = cos(qJ(2));
t64 = t41 * t50;
t63 = t42 * t46;
t29 = t46 * t32;
t24 = t39 * t45 - t42 * t48;
t62 = t49 * t24;
t61 = t49 * t32;
t60 = t49 * t39;
t59 = t49 * t42;
t58 = pkin(9) + qJ(5);
t43 = cos(pkin(11));
t33 = -pkin(2) * t43 - pkin(3);
t23 = -qJ(5) * t46 + t33 - t69;
t10 = t23 * t39 + t32 * t59;
t57 = t39 ^ 2 + t42 ^ 2;
t56 = t57 * qJ(5);
t17 = (t40 * t50 + t43 * t47) * t41;
t44 = cos(pkin(6));
t12 = t17 * t49 + t44 * t46;
t15 = t40 * t65 - t43 * t64;
t5 = -t12 * t39 + t15 * t42;
t6 = t12 * t42 + t15 * t39;
t55 = -t39 * t5 + t42 * t6;
t21 = t42 * t23;
t9 = -t32 * t60 + t21;
t54 = t10 * t42 - t39 * t9;
t53 = -pkin(4) * t46 + qJ(5) * t49;
t38 = t49 ^ 2;
t37 = t46 ^ 2;
t28 = t58 * t42;
t27 = t58 * t39;
t22 = pkin(5) * t66 + t29;
t19 = t24 * t46;
t18 = t25 * t46;
t14 = -t27 * t45 + t28 * t48;
t13 = -t27 * t48 - t28 * t45;
t11 = t17 * t46 - t44 * t49;
t8 = -pkin(9) * t66 + t10;
t7 = -pkin(9) * t63 + t21 + (-pkin(5) - t67) * t49;
t4 = t45 * t7 + t48 * t8;
t3 = -t45 * t8 + t48 * t7;
t2 = t45 * t5 + t48 * t6;
t1 = -t45 * t6 + t48 * t5;
t16 = [1, 0, 0, 0, t15 ^ 2 + t17 ^ 2 + t44 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 + t5 ^ 2 + t6 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t64, -t65 (-t15 * t43 + t17 * t40) * pkin(2), 0, 0, 0, 0, 0, -t15 * t49, t15 * t46, t11 * t66 - t49 * t5, t11 * t63 + t49 * t6 (-t39 * t6 - t42 * t5) * t46, t10 * t6 + t11 * t29 + t5 * t9, 0, 0, 0, 0, 0, -t1 * t49 + t11 * t18, -t11 * t19 + t2 * t49; 0, 1, 0, 0 (t40 ^ 2 + t43 ^ 2) * pkin(2) ^ 2, t37, t46 * t70, 0, 0, 0, t33 * t71, t33 * t72, 0.2e1 * t37 * t67 - 0.2e1 * t49 * t9, 0.2e1 * t32 * t37 * t42 + 0.2e1 * t10 * t49 (-t10 * t39 - t42 * t9) * t72, t32 ^ 2 * t37 + t10 ^ 2 + t9 ^ 2, t19 ^ 2, 0.2e1 * t19 * t18, -t19 * t71, t18 * t70, t38, 0.2e1 * t18 * t22 - 0.2e1 * t3 * t49, -0.2e1 * t19 * t22 + 0.2e1 * t4 * t49; 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t49 + t46 * t55, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t54 - t61) * t46, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t57 + t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -t11 * t42, t11 * t39, t55, -t11 * pkin(4) + qJ(5) * t55, 0, 0, 0, 0, 0, t11 * t24, t11 * t25; 0, 0, 0, 0, 0, 0, 0, t46, t49, 0, -t29, -t61, -t29 * t42 + t39 * t53, t29 * t39 + t42 * t53, t54, -pkin(4) * t29 + qJ(5) * t54, -t19 * t25, -t18 * t25 + t19 * t24, -t68, t62, 0, -t13 * t49 + t18 * t34 + t22 * t24, t14 * t49 - t19 * t34 + t22 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46, t59, -t60, t57 * t46, t46 * t56 + t69, 0, 0, 0, 0, 0, -t62, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t42, -0.2e1 * pkin(4) * t39, 0.2e1 * t56, qJ(5) ^ 2 * t57 + pkin(4) ^ 2, t25 ^ 2, -0.2e1 * t25 * t24, 0, 0, 0, t24 * t73, t25 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t63, 0, t29, 0, 0, 0, 0, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t39, 0, -pkin(4), 0, 0, 0, 0, 0, t24, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t18, -t49, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t16;
