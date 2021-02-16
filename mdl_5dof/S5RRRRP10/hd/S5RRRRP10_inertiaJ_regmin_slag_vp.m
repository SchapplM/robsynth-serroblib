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
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-16 00:35:52
% EndTime: 2021-01-16 00:35:56
% DurationCPUTime: 0.72s
% Computational Cost: add. (899->140), mult. (2143->279), div. (0->0), fcn. (2341->8), ass. (0->78)
t42 = cos(pkin(5));
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t41 = sin(pkin(5));
t45 = sin(qJ(2));
t68 = t41 * t45;
t23 = t42 * t44 + t47 * t68;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t48 = cos(qJ(2));
t67 = t41 * t48;
t13 = t23 * t43 + t46 * t67;
t82 = -0.2e1 * t13;
t81 = -0.2e1 * t23;
t80 = -0.2e1 * t44;
t79 = 0.2e1 * t47;
t78 = pkin(1) * t45;
t77 = pkin(1) * t48;
t76 = pkin(3) * t46;
t75 = pkin(8) * t43;
t22 = -t42 * t47 + t44 * t68;
t74 = t22 * pkin(4);
t73 = t43 * pkin(4);
t53 = pkin(7) * t67;
t20 = t53 + (pkin(8) + t78) * t42;
t21 = (-pkin(2) * t48 - pkin(8) * t45 - pkin(1)) * t41;
t10 = -t44 * t20 + t47 * t21;
t8 = pkin(3) * t67 - t10;
t72 = t8 * t43;
t71 = t8 * t46;
t14 = t23 * t46 - t43 * t67;
t70 = t14 * t43;
t37 = t41 ^ 2;
t69 = t37 * t48;
t66 = t42 * t45;
t65 = t43 * t22;
t64 = t43 * t44;
t63 = t43 * t46;
t62 = t43 * t47;
t61 = t44 * t22;
t60 = t46 * t22;
t35 = t46 * t44;
t59 = t46 * t47;
t58 = qJ(5) + pkin(9);
t57 = qJ(5) * t44;
t56 = 0.2e1 * t67;
t55 = t44 * t79;
t54 = pkin(8) * t59;
t52 = t44 * t67;
t51 = t47 * t67;
t32 = pkin(7) * t68;
t19 = t32 + (-pkin(2) - t77) * t42;
t7 = t22 * pkin(3) - t23 * pkin(9) + t19;
t11 = t47 * t20 + t44 * t21;
t9 = -pkin(9) * t67 + t11;
t3 = -t43 * t9 + t46 * t7;
t28 = -t47 * pkin(3) - t44 * pkin(9) - pkin(2);
t26 = t46 * t28;
t50 = -t46 * t57 + t26;
t4 = t43 * t7 + t46 * t9;
t49 = -t14 * qJ(5) + t3;
t40 = t46 ^ 2;
t39 = t44 ^ 2;
t38 = t43 ^ 2;
t36 = -t46 * pkin(4) - pkin(3);
t30 = t58 * t46;
t29 = t58 * t43;
t27 = (pkin(8) + t73) * t44;
t25 = pkin(1) * t66 + t53;
t24 = t42 * t77 - t32;
t18 = t43 * t28 + t54;
t17 = -pkin(8) * t62 + t26;
t15 = t54 + (t28 - t57) * t43;
t12 = (-pkin(4) - t75) * t47 + t50;
t5 = t13 * pkin(4) + t8;
t2 = -t13 * qJ(5) + t4;
t1 = t49 + t74;
t6 = [1, 0, 0, t37 * t45 ^ 2, 0.2e1 * t45 * t69, 0.2e1 * t41 * t66, t42 * t56, t42 ^ 2, 0.2e1 * pkin(1) * t69 + 0.2e1 * t24 * t42, -0.2e1 * t25 * t42 - 0.2e1 * t37 * t78, t23 ^ 2, t22 * t81, t67 * t81, t22 * t56, t37 * t48 ^ 2, -0.2e1 * t10 * t67 + 0.2e1 * t19 * t22, 0.2e1 * t11 * t67 + 0.2e1 * t19 * t23, t14 ^ 2, t14 * t82, 0.2e1 * t14 * t22, t22 * t82, t22 ^ 2, 0.2e1 * t8 * t13 + 0.2e1 * t3 * t22, 0.2e1 * t8 * t14 - 0.2e1 * t4 * t22, 0.2e1 * t1 * t22 + 0.2e1 * t5 * t13, 0.2e1 * t5 * t14 - 0.2e1 * t2 * t22, -0.2e1 * t1 * t14 - 0.2e1 * t2 * t13, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t68, t67, t42, t24, -t25, t23 * t44, t23 * t47 - t61, -t52, -t51, 0, -pkin(2) * t22 + pkin(8) * t52 - t19 * t47, -pkin(2) * t23 + pkin(8) * t51 + t19 * t44, t14 * t35, (-t13 * t46 - t70) * t44, -t14 * t47 + t22 * t35, t13 * t47 - t43 * t61, -t22 * t47, t17 * t22 - t3 * t47 + (pkin(8) * t13 + t72) * t44, -t18 * t22 + t4 * t47 + (pkin(8) * t14 + t71) * t44, -t1 * t47 + t12 * t22 + t27 * t13 + t5 * t64, t27 * t14 - t15 * t22 + t2 * t47 + t5 * t35, -t12 * t14 - t15 * t13 + (-t1 * t46 - t2 * t43) * t44, t1 * t12 + t2 * t15 + t5 * t27; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t39, t55, 0, 0, 0, pkin(2) * t79, pkin(2) * t80, t40 * t39, -0.2e1 * t39 * t63, t59 * t80, t43 * t55, t47 ^ 2, -0.2e1 * t17 * t47 + 0.2e1 * t39 * t75, 0.2e1 * t39 * pkin(8) * t46 + 0.2e1 * t18 * t47, -0.2e1 * t12 * t47 + 0.2e1 * t27 * t64, 0.2e1 * t15 * t47 + 0.2e1 * t27 * t35, 0.2e1 * (-t12 * t46 - t15 * t43) * t44, t12 ^ 2 + t15 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, -t67, t10, -t11, t70, -t43 * t13 + t14 * t46, t65, t60, 0, -pkin(3) * t13 - pkin(9) * t65 - t71, -pkin(3) * t14 - pkin(9) * t60 + t72, t36 * t13 - t29 * t22 - t5 * t46, t36 * t14 - t30 * t22 + t5 * t43, -t1 * t43 - t30 * t13 + t29 * t14 + t2 * t46, -t1 * t29 + t2 * t30 + t5 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t47, 0, -t44 * pkin(8), -t47 * pkin(8), t43 * t35, (-t38 + t40) * t44, -t62, -t59, 0, -pkin(8) * t35 + (-pkin(3) * t44 + pkin(9) * t47) * t43, pkin(9) * t59 + (t75 - t76) * t44, -t27 * t46 + t29 * t47 + t36 * t64, t27 * t43 + t30 * t47 + t36 * t35, (t29 * t44 + t15) * t46 + (-t30 * t44 - t12) * t43, -t12 * t29 + t15 * t30 + t27 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, 0.2e1 * t63, 0, 0, 0, 0.2e1 * t76, -0.2e1 * pkin(3) * t43, -0.2e1 * t36 * t46, 0.2e1 * t36 * t43, 0.2e1 * t29 * t43 + 0.2e1 * t30 * t46, t29 ^ 2 + t30 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t22, t3, -t4, t49 + 0.2e1 * t74, -t2, -t14 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t64, -t47, t17, -t18, (-0.2e1 * pkin(4) - t75) * t47 + t50, -t15, -pkin(4) * t35, t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t46, 0, -t43 * pkin(9), -t46 * pkin(9), -t29, -t30, -t73, -t29 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t35, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t43, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
