% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:55
% EndTime: 2020-01-03 11:22:59
% DurationCPUTime: 0.61s
% Computational Cost: add. (539->83), mult. (1192->185), div. (0->0), fcn. (1305->8), ass. (0->57)
t42 = sin(pkin(9));
t45 = cos(pkin(9));
t47 = cos(pkin(7));
t44 = sin(pkin(7));
t46 = cos(pkin(8));
t60 = t46 * t44;
t18 = t42 * t60 + t45 * t47;
t68 = t18 ^ 2;
t67 = -0.2e1 * t18;
t66 = -0.2e1 * t44;
t65 = 0.2e1 * t47;
t24 = -pkin(2) * t47 - qJ(3) * t44 - pkin(1);
t43 = sin(pkin(8));
t54 = qJ(2) * t47;
t16 = t43 * t24 + t46 * t54;
t12 = -qJ(4) * t47 + t16;
t32 = t44 * qJ(2);
t17 = t32 + (pkin(3) * t43 - qJ(4) * t46) * t44;
t7 = t45 * t12 + t42 * t17;
t64 = t42 * t43;
t63 = t43 * t45;
t62 = t44 * t43;
t61 = t44 * t45;
t59 = t46 * t47;
t48 = sin(qJ(5));
t58 = t48 * t42;
t49 = cos(qJ(5));
t57 = t49 * t42;
t35 = t43 ^ 2;
t38 = t46 ^ 2;
t56 = t35 + t38;
t55 = t48 ^ 2 + t49 ^ 2;
t36 = t44 ^ 2;
t53 = t36 * qJ(2);
t52 = t44 * t65;
t15 = t24 * t46 - t43 * t54;
t13 = t47 * pkin(3) - t15;
t6 = -t12 * t42 + t17 * t45;
t51 = t15 * t46 + t16 * t43;
t50 = qJ(2) ^ 2;
t39 = t47 ^ 2;
t37 = t45 ^ 2;
t34 = t42 ^ 2;
t31 = t36 * t50;
t29 = t35 * t36;
t28 = t34 * t35;
t22 = -t46 * t48 + t49 * t63;
t21 = -t46 * t49 - t48 * t63;
t20 = -t47 * t42 + t45 * t60;
t11 = t20 * t49 + t48 * t62;
t9 = t20 * t48 - t49 * t62;
t5 = pkin(6) * t62 + t7;
t4 = -pkin(4) * t62 - t6;
t3 = pkin(4) * t18 - pkin(6) * t20 + t13;
t2 = t3 * t48 + t49 * t5;
t1 = t3 * t49 - t48 * t5;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, t52, 0, t39, 0, 0, pkin(1) * t65, pkin(1) * t66, 0.2e1 * (t36 + t39) * qJ(2), pkin(1) ^ 2 + t39 * t50 + t31, t38 * t36, -0.2e1 * t46 * t36 * t43, t59 * t66, t29, t43 * t52, t39, -0.2e1 * t15 * t47 + 0.2e1 * t43 * t53, 0.2e1 * t16 * t47 + 0.2e1 * t46 * t53, t51 * t66, t15 ^ 2 + t16 ^ 2 + t31, t20 ^ 2, t20 * t67, 0.2e1 * t20 * t62, t68, t62 * t67, t29, 0.2e1 * t13 * t18 + 0.2e1 * t6 * t62, 0.2e1 * t13 * t20 - 0.2e1 * t62 * t7, -0.2e1 * t18 * t7 - 0.2e1 * t20 * t6, t13 ^ 2 + t6 ^ 2 + t7 ^ 2, t11 ^ 2, -0.2e1 * t11 * t9, 0.2e1 * t11 * t18, t9 ^ 2, t9 * t67, t68, 0.2e1 * t1 * t18 + 0.2e1 * t4 * t9, 0.2e1 * t11 * t4 - 0.2e1 * t18 * t2, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t2 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t44, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t59, t43 * t47, -t56 * t44, t51, 0, 0, 0, 0, 0, 0, -t35 * t42 * t44 - t18 * t46, -t20 * t46 - t35 * t61, (-t18 * t45 + t20 * t42) * t43, -t13 * t46 + (-t42 * t6 + t45 * t7) * t43, 0, 0, 0, 0, 0, 0, t18 * t21 + t64 * t9, t11 * t64 - t18 * t22, -t11 * t21 - t22 * t9, t1 * t21 + t2 * t22 + t4 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t37 + t28 + t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 ^ 2 + t22 ^ 2 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t60, 0, t32, 0, 0, 0, 0, 0, 0, t43 * t61, -t42 * t62, -t18 * t42 - t20 * t45, t42 * t7 + t45 * t6, 0, 0, 0, 0, 0, 0, -t18 * t58 - t45 * t9, -t11 * t45 - t18 * t57, (t11 * t48 - t49 * t9) * t42, -t4 * t45 + (-t1 * t48 + t2 * t49) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t21 * t48 + t22 * t49 - t63) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 + t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t55 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t20, 0, t13, 0, 0, 0, 0, 0, 0, t49 * t18, -t48 * t18, -t11 * t49 - t48 * t9, t1 * t49 + t2 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t49 + t22 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, t18, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
