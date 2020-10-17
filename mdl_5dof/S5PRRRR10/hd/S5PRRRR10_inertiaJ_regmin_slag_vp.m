% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:06
% EndTime: 2019-12-05 17:26:08
% DurationCPUTime: 0.50s
% Computational Cost: add. (398->99), mult. (1133->214), div. (0->0), fcn. (1360->12), ass. (0->75)
t38 = cos(pkin(6));
t41 = sin(qJ(4));
t45 = cos(qJ(4));
t36 = sin(pkin(6));
t42 = sin(qJ(3));
t66 = t36 * t42;
t24 = t41 * t38 + t45 * t66;
t40 = sin(qJ(5));
t44 = cos(qJ(5));
t46 = cos(qJ(3));
t65 = t36 * t46;
t14 = t40 * t24 + t44 * t65;
t78 = -0.2e1 * t14;
t77 = -0.2e1 * t24;
t76 = -0.2e1 * t41;
t75 = 0.2e1 * t45;
t74 = pkin(2) * t42;
t73 = pkin(2) * t46;
t72 = pkin(4) * t44;
t71 = pkin(9) * t40;
t50 = pkin(8) * t65;
t20 = t50 + (pkin(9) + t74) * t38;
t21 = (-pkin(3) * t46 - pkin(9) * t42 - pkin(2)) * t36;
t10 = -t41 * t20 + t45 * t21;
t8 = pkin(4) * t65 - t10;
t70 = t8 * t40;
t69 = t8 * t44;
t15 = t44 * t24 - t40 * t65;
t68 = t15 * t40;
t32 = t36 ^ 2;
t67 = t32 * t46;
t37 = sin(pkin(5));
t43 = sin(qJ(2));
t64 = t37 * t43;
t47 = cos(qJ(2));
t63 = t37 * t47;
t62 = t38 * t42;
t61 = t38 * t47;
t23 = -t45 * t38 + t41 * t66;
t60 = t40 * t23;
t59 = t40 * t41;
t58 = t40 * t44;
t57 = t40 * t45;
t56 = t41 * t23;
t55 = t44 * t23;
t54 = t44 * t41;
t53 = t44 * t45;
t52 = 0.2e1 * t65;
t51 = t41 * t75;
t49 = t41 * t65;
t48 = t45 * t65;
t11 = t45 * t20 + t41 * t21;
t30 = pkin(8) * t66;
t19 = t30 + (-pkin(3) - t73) * t38;
t39 = cos(pkin(5));
t35 = t44 ^ 2;
t34 = t41 ^ 2;
t33 = t40 ^ 2;
t28 = -t45 * pkin(4) - t41 * pkin(10) - pkin(3);
t26 = pkin(2) * t62 + t50;
t25 = t38 * t73 - t30;
t22 = -t36 * t63 + t39 * t38;
t17 = pkin(9) * t53 + t40 * t28;
t16 = -pkin(9) * t57 + t44 * t28;
t13 = t39 * t66 + (t42 * t61 + t43 * t46) * t37;
t12 = -t37 * t46 * t61 - t39 * t65 + t42 * t64;
t9 = -pkin(10) * t65 + t11;
t7 = t23 * pkin(4) - t24 * pkin(10) + t19;
t6 = t13 * t45 + t22 * t41;
t5 = t13 * t41 - t22 * t45;
t4 = t12 * t40 + t6 * t44;
t3 = t12 * t44 - t6 * t40;
t2 = t40 * t7 + t44 * t9;
t1 = -t40 * t9 + t44 * t7;
t18 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t63, -t64, 0, 0, 0, 0, 0, -t12 * t38 - t22 * t65, -t13 * t38 + t22 * t66, 0, 0, 0, 0, 0, t12 * t23 + t5 * t65, t12 * t24 + t6 * t65, 0, 0, 0, 0, 0, t5 * t14 + t3 * t23, t5 * t15 - t4 * t23; 0, 1, 0, 0, t32 * t42 ^ 2, 0.2e1 * t42 * t67, 0.2e1 * t36 * t62, t38 * t52, t38 ^ 2, 0.2e1 * pkin(2) * t67 + 0.2e1 * t25 * t38, -0.2e1 * t26 * t38 - 0.2e1 * t32 * t74, t24 ^ 2, t23 * t77, t65 * t77, t23 * t52, t32 * t46 ^ 2, -0.2e1 * t10 * t65 + 0.2e1 * t19 * t23, 0.2e1 * t11 * t65 + 0.2e1 * t19 * t24, t15 ^ 2, t15 * t78, 0.2e1 * t15 * t23, t23 * t78, t23 ^ 2, 0.2e1 * t1 * t23 + 0.2e1 * t8 * t14, 0.2e1 * t8 * t15 - 0.2e1 * t2 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t12 * t45, t12 * t41, 0, 0, 0, 0, 0, -t3 * t45 + t5 * t59, t4 * t45 + t5 * t54; 0, 0, 0, 0, 0, 0, t66, t65, t38, t25, -t26, t24 * t41, t24 * t45 - t56, -t49, -t48, 0, -pkin(3) * t23 + pkin(9) * t49 - t19 * t45, -pkin(3) * t24 + pkin(9) * t48 + t19 * t41, t15 * t54, (-t14 * t44 - t68) * t41, -t15 * t45 + t23 * t54, t14 * t45 - t40 * t56, -t23 * t45, -t1 * t45 + t16 * t23 + (pkin(9) * t14 + t70) * t41, -t17 * t23 + t2 * t45 + (pkin(9) * t15 + t69) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t34, t51, 0, 0, 0, pkin(3) * t75, pkin(3) * t76, t35 * t34, -0.2e1 * t34 * t58, t53 * t76, t40 * t51, t45 ^ 2, -0.2e1 * t16 * t45 + 0.2e1 * t34 * t71, 0.2e1 * t34 * pkin(9) * t44 + 0.2e1 * t17 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t44, t5 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, -t65, t10, -t11, t68, -t40 * t14 + t15 * t44, t60, t55, 0, -pkin(4) * t14 - pkin(10) * t60 - t69, -pkin(4) * t15 - pkin(10) * t55 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t45, 0, -t41 * pkin(9), -t45 * pkin(9), t40 * t54, (-t33 + t35) * t41, -t57, -t53, 0, -pkin(9) * t54 + (-pkin(4) * t41 + pkin(10) * t45) * t40, pkin(10) * t53 + (t71 - t72) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t33, 0.2e1 * t58, 0, 0, 0, 0.2e1 * t72, -0.2e1 * pkin(4) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t23, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t59, -t45, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t44, 0, -t40 * pkin(10), -t44 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t18;
