% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:41:33
% EndTime: 2019-05-05 15:41:34
% DurationCPUTime: 0.55s
% Computational Cost: add. (354->91), mult. (628->154), div. (0->0), fcn. (725->8), ass. (0->65)
t42 = sin(qJ(6));
t45 = cos(qJ(6));
t44 = sin(qJ(4));
t46 = cos(qJ(5));
t54 = t46 * t44;
t43 = sin(qJ(5));
t58 = t43 * t44;
t16 = -t42 * t58 + t45 * t54;
t72 = -0.2e1 * t16;
t34 = -pkin(5) * t46 - pkin(4);
t71 = 0.2e1 * t34;
t70 = -0.2e1 * t44;
t47 = cos(qJ(4));
t69 = 0.2e1 * t47;
t68 = pkin(8) + pkin(9);
t67 = pkin(4) * t46;
t66 = t42 * pkin(5);
t65 = t45 * pkin(5);
t40 = sin(pkin(10));
t41 = cos(pkin(10));
t48 = -pkin(1) - pkin(2);
t24 = qJ(2) * t40 - t41 * t48;
t20 = pkin(3) + t24;
t14 = pkin(4) * t47 + pkin(8) * t44 + t20;
t26 = qJ(2) * t41 + t40 * t48;
t21 = -pkin(7) + t26;
t53 = t47 * t21;
t49 = t46 * t53;
t5 = t49 + (pkin(9) * t44 + t14) * t43;
t64 = t45 * t5;
t63 = t47 * pkin(5);
t62 = t21 * t43;
t23 = t42 * t46 + t43 * t45;
t61 = t23 * t47;
t37 = t44 ^ 2;
t60 = t37 * t43;
t59 = t37 * t46;
t57 = t43 * t46;
t56 = t43 * t47;
t55 = t44 * t40;
t33 = t46 * t47;
t22 = t42 * t43 - t45 * t46;
t52 = t47 * t22;
t51 = t47 * t40;
t50 = t44 * t69;
t12 = t46 * t14;
t4 = pkin(9) * t54 + t12 + (pkin(5) - t62) * t47;
t1 = t4 * t45 - t42 * t5;
t39 = t47 ^ 2;
t38 = t46 ^ 2;
t36 = t43 ^ 2;
t28 = t68 * t46;
t27 = t68 * t43;
t19 = -t41 * t43 + t46 * t51;
t18 = -t41 * t46 - t43 * t51;
t15 = t23 * t44;
t13 = (-pkin(5) * t43 + t21) * t44;
t11 = -t27 * t42 + t28 * t45;
t10 = -t27 * t45 - t28 * t42;
t9 = t18 * t42 + t19 * t45;
t8 = t18 * t45 - t19 * t42;
t7 = t14 * t43 + t49;
t6 = -t43 * t53 + t12;
t2 = t4 * t42 + t64;
t3 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t24, 0.2e1 * t26, t24 ^ 2 + t26 ^ 2, t37, t50, 0, 0, 0, t20 * t69, t20 * t70, t38 * t37, -0.2e1 * t37 * t57, t33 * t70, t43 * t50, t39, -0.2e1 * t21 * t60 + 0.2e1 * t47 * t6, -0.2e1 * t21 * t59 - 0.2e1 * t47 * t7, t16 ^ 2, t15 * t72, t47 * t72, t15 * t69, t39, 0.2e1 * t1 * t47 - 0.2e1 * t13 * t15, -0.2e1 * t13 * t16 - 0.2e1 * t2 * t47; 0, 0, 0, -1, 0, -pkin(1), -t41, t40, -t24 * t41 + t26 * t40, 0, 0, 0, 0, 0, -t41 * t47, t44 * t41, 0, 0, 0, 0, 0, t18 * t47 - t40 * t60, -t19 * t47 - t40 * t59, 0, 0, 0, 0, 0, -t15 * t55 + t47 * t8, -t16 * t55 - t47 * t9; 0, 0, 0, 0, 0, 1, 0, 0, t40 ^ 2 + t41 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t47, 0, -t44 * t21, -t53, -t43 * t54 (t36 - t38) * t44, t56, t33, 0, -t21 * t54 + (pkin(4) * t44 - pkin(8) * t47) * t43, -pkin(8) * t33 + (t62 + t67) * t44, -t16 * t23, t15 * t23 + t16 * t22, t61, -t52, 0, t10 * t47 + t13 * t22 - t15 * t34, -t11 * t47 + t13 * t23 - t16 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t51, 0, 0, 0, 0, 0, -t40 * t54, t43 * t55, 0, 0, 0, 0, 0, t22 * t55, t23 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t44, 0, 0, 0, 0, 0, t33, -t56, 0, 0, 0, 0, 0, -t52, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t36, 0.2e1 * t57, 0, 0, 0, 0.2e1 * t67, -0.2e1 * pkin(4) * t43, t23 ^ 2, -0.2e1 * t23 * t22, 0, 0, 0, t22 * t71, t23 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t58, t47, t6, -t7, 0, 0, -t16, t15, t47, t45 * t63 + t1, -t64 + (-t4 - t63) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, 0, 0, 0, 0, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t54, 0, 0, 0, 0, 0, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t46, 0, -t43 * pkin(8), -t46 * pkin(8), 0, 0, t23, -t22, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, t47, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
