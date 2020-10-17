% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:21
% EndTime: 2019-12-31 22:12:23
% DurationCPUTime: 0.57s
% Computational Cost: add. (657->113), mult. (1591->235), div. (0->0), fcn. (1735->8), ass. (0->74)
t41 = cos(pkin(5));
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t40 = sin(pkin(5));
t44 = sin(qJ(2));
t65 = t40 * t44;
t23 = t41 * t43 + t46 * t65;
t42 = sin(qJ(4));
t45 = cos(qJ(4));
t47 = cos(qJ(2));
t64 = t40 * t47;
t13 = t23 * t42 + t45 * t64;
t78 = -0.2e1 * t13;
t77 = -0.2e1 * t23;
t76 = -0.2e1 * t43;
t75 = 0.2e1 * t46;
t74 = pkin(1) * t44;
t73 = pkin(1) * t47;
t72 = pkin(3) * t45;
t71 = pkin(8) * t42;
t70 = t42 * pkin(4);
t50 = pkin(7) * t64;
t20 = t50 + (pkin(8) + t74) * t41;
t21 = (-pkin(2) * t47 - pkin(8) * t44 - pkin(1)) * t40;
t10 = -t43 * t20 + t46 * t21;
t8 = pkin(3) * t64 - t10;
t69 = t8 * t42;
t68 = t8 * t45;
t14 = t23 * t45 - t42 * t64;
t67 = t14 * t42;
t36 = t40 ^ 2;
t66 = t36 * t47;
t63 = t41 * t44;
t22 = -t41 * t46 + t43 * t65;
t62 = t42 * t22;
t61 = t42 * t45;
t60 = t42 * t46;
t59 = t43 * t22;
t58 = t45 * t22;
t57 = t45 * t43;
t56 = t45 * t46;
t55 = -qJ(5) - pkin(9);
t54 = qJ(5) * t43;
t53 = 0.2e1 * t64;
t52 = t43 * t75;
t51 = pkin(8) * t56;
t49 = t43 * t64;
t48 = t46 * t64;
t32 = pkin(7) * t65;
t19 = t32 + (-pkin(2) - t73) * t41;
t7 = t22 * pkin(3) - t23 * pkin(9) + t19;
t11 = t46 * t20 + t43 * t21;
t9 = -pkin(9) * t64 + t11;
t3 = -t42 * t9 + t45 * t7;
t4 = t42 * t7 + t45 * t9;
t39 = t45 ^ 2;
t38 = t43 ^ 2;
t37 = t42 ^ 2;
t35 = -t45 * pkin(4) - pkin(3);
t30 = t55 * t45;
t29 = t55 * t42;
t28 = -t46 * pkin(3) - t43 * pkin(9) - pkin(2);
t27 = (pkin(8) + t70) * t43;
t26 = t45 * t28;
t25 = pkin(1) * t63 + t50;
t24 = t41 * t73 - t32;
t18 = t42 * t28 + t51;
t17 = -pkin(8) * t60 + t26;
t15 = t51 + (t28 - t54) * t42;
t12 = -t45 * t54 + t26 + (-pkin(4) - t71) * t46;
t5 = t13 * pkin(4) + t8;
t2 = -t13 * qJ(5) + t4;
t1 = t22 * pkin(4) - t14 * qJ(5) + t3;
t6 = [1, 0, 0, t36 * t44 ^ 2, 0.2e1 * t44 * t66, 0.2e1 * t40 * t63, t41 * t53, t41 ^ 2, 0.2e1 * pkin(1) * t66 + 0.2e1 * t24 * t41, -0.2e1 * t25 * t41 - 0.2e1 * t36 * t74, t23 ^ 2, t22 * t77, t64 * t77, t22 * t53, t36 * t47 ^ 2, -0.2e1 * t10 * t64 + 0.2e1 * t19 * t22, 0.2e1 * t11 * t64 + 0.2e1 * t19 * t23, t14 ^ 2, t14 * t78, 0.2e1 * t14 * t22, t22 * t78, t22 ^ 2, 0.2e1 * t8 * t13 + 0.2e1 * t3 * t22, 0.2e1 * t8 * t14 - 0.2e1 * t4 * t22, -0.2e1 * t1 * t14 - 0.2e1 * t2 * t13, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t65, t64, t41, t24, -t25, t23 * t43, t23 * t46 - t59, -t49, -t48, 0, -pkin(2) * t22 + pkin(8) * t49 - t19 * t46, -pkin(2) * t23 + pkin(8) * t48 + t19 * t43, t14 * t57, (-t13 * t45 - t67) * t43, -t14 * t46 + t22 * t57, t13 * t46 - t42 * t59, -t22 * t46, t17 * t22 - t3 * t46 + (pkin(8) * t13 + t69) * t43, -t18 * t22 + t4 * t46 + (pkin(8) * t14 + t68) * t43, -t12 * t14 - t15 * t13 + (-t1 * t45 - t2 * t42) * t43, t1 * t12 + t2 * t15 + t5 * t27; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, t52, 0, 0, 0, pkin(2) * t75, pkin(2) * t76, t39 * t38, -0.2e1 * t38 * t61, t56 * t76, t42 * t52, t46 ^ 2, -0.2e1 * t17 * t46 + 0.2e1 * t38 * t71, 0.2e1 * t38 * pkin(8) * t45 + 0.2e1 * t18 * t46, 0.2e1 * (-t12 * t45 - t15 * t42) * t43, t12 ^ 2 + t15 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, -t64, t10, -t11, t67, -t42 * t13 + t14 * t45, t62, t58, 0, -pkin(3) * t13 - pkin(9) * t62 - t68, -pkin(3) * t14 - pkin(9) * t58 + t69, -t1 * t42 + t30 * t13 - t29 * t14 + t2 * t45, t1 * t29 - t2 * t30 + t5 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t46, 0, -t43 * pkin(8), -t46 * pkin(8), t42 * t57, (-t37 + t39) * t43, -t60, -t56, 0, -pkin(8) * t57 + (-pkin(3) * t43 + pkin(9) * t46) * t42, pkin(9) * t56 + (t71 - t72) * t43, (-t29 * t43 + t15) * t45 + (t30 * t43 - t12) * t42, t12 * t29 - t15 * t30 + t27 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t37, 0.2e1 * t61, 0, 0, 0, 0.2e1 * t72, -0.2e1 * pkin(3) * t42, -0.2e1 * t29 * t42 - 0.2e1 * t30 * t45, t29 ^ 2 + t30 ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t22, t3, -t4, -pkin(4) * t14, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t42 * t43, -t46, t17, -t18, -pkin(4) * t57, t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t45, 0, -t42 * pkin(9), -t45 * pkin(9), -t70, t29 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
