% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP6
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
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:36:47
% EndTime: 2019-05-06 01:36:49
% DurationCPUTime: 0.61s
% Computational Cost: add. (963->100), mult. (1894->183), div. (0->0), fcn. (2271->8), ass. (0->70)
t52 = sin(pkin(10));
t53 = cos(pkin(10));
t56 = sin(qJ(3));
t73 = cos(qJ(3));
t31 = t56 * t52 - t53 * t73;
t81 = -0.2e1 * t31;
t80 = 0.2e1 * t31;
t43 = -t53 * pkin(2) - pkin(1);
t79 = 0.2e1 * t43;
t57 = cos(qJ(4));
t46 = -t57 * pkin(4) - pkin(3);
t78 = 0.2e1 * t46;
t77 = pkin(8) + pkin(9);
t76 = t31 * pkin(4);
t54 = sin(qJ(5));
t55 = sin(qJ(4));
t72 = cos(qJ(5));
t60 = t72 * t57;
t34 = t54 * t55 - t60;
t75 = t34 * pkin(5);
t74 = t54 * pkin(4);
t32 = t52 * t73 + t56 * t53;
t68 = t55 * t32;
t16 = t32 * t60 - t54 * t68;
t71 = t16 * t34;
t36 = t54 * t57 + t55 * t72;
t70 = t36 * t31;
t69 = t55 * t31;
t67 = t55 * t57;
t64 = pkin(7) + qJ(2);
t38 = t64 * t52;
t39 = t64 * t53;
t20 = -t56 * t38 + t39 * t73;
t66 = t57 * t20;
t65 = t57 * t32;
t63 = t52 ^ 2 + t53 ^ 2;
t62 = t32 * t81;
t47 = t72 * pkin(4);
t18 = t31 * pkin(3) - t32 * pkin(8) + t43;
t7 = t66 + (-pkin(9) * t32 + t18) * t55;
t61 = t72 * t7;
t9 = t57 * t18 - t55 * t20;
t6 = -pkin(9) * t65 + t76 + t9;
t3 = -t54 * t7 + t72 * t6;
t40 = t77 * t55;
t41 = t77 * t57;
t22 = -t72 * t40 - t54 * t41;
t59 = -pkin(3) * t32 - pkin(8) * t31;
t19 = t38 * t73 + t56 * t39;
t12 = pkin(4) * t68 + t19;
t4 = t54 * t6 + t61;
t23 = -t54 * t40 + t41 * t72;
t51 = t57 ^ 2;
t50 = t55 ^ 2;
t45 = t47 + pkin(5);
t33 = t36 ^ 2;
t29 = t32 ^ 2;
t28 = t31 ^ 2;
t27 = t57 * t31;
t24 = t46 + t75;
t21 = t34 * t31;
t15 = t36 * t32;
t14 = -t34 * qJ(6) + t23;
t13 = -t36 * qJ(6) + t22;
t11 = t36 * t15;
t10 = t55 * t18 + t66;
t8 = t15 * pkin(5) + t12;
t2 = -t15 * qJ(6) + t4;
t1 = t31 * pkin(5) - t16 * qJ(6) + t3;
t5 = [1, 0, 0, 0.2e1 * pkin(1) * t53, -0.2e1 * pkin(1) * t52, 0.2e1 * t63 * qJ(2), qJ(2) ^ 2 * t63 + pkin(1) ^ 2, t29, t62, 0, 0, 0, t31 * t79, t32 * t79, t51 * t29, -0.2e1 * t29 * t67, t65 * t80, t55 * t62, t28, 0.2e1 * t19 * t68 + 0.2e1 * t9 * t31, -0.2e1 * t10 * t31 + 0.2e1 * t19 * t65, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t80, t15 * t81, t28, 0.2e1 * t12 * t15 + 0.2e1 * t3 * t31, 0.2e1 * t12 * t16 - 0.2e1 * t4 * t31, -0.2e1 * t1 * t16 - 0.2e1 * t2 * t15, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, -t53, t52, 0, -pkin(1), 0, 0, 0, 0, 0, t31, t32, 0, 0, 0, 0, 0, t27, -t69, 0, 0, 0, 0, 0, -t21, -t70, -t11 + t71, -t1 * t34 + t2 * t36; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 ^ 2 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, -t19, -t20, t55 * t65 (-t50 + t51) * t32, t69, t27, 0, -t19 * t57 + t55 * t59, t19 * t55 + t57 * t59, t16 * t36, -t11 - t71, t70, -t21, 0, t12 * t34 + t46 * t15 + t22 * t31, t12 * t36 + t46 * t16 - t23 * t31, -t1 * t36 - t13 * t16 - t14 * t15 - t2 * t34, t1 * t13 + t2 * t14 + t8 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t13 + t36 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t50, 0.2e1 * t67, 0, 0, 0, 0.2e1 * pkin(3) * t57, -0.2e1 * pkin(3) * t55, t33, -0.2e1 * t36 * t34, 0, 0, 0, t34 * t78, t36 * t78, -0.2e1 * t13 * t36 - 0.2e1 * t14 * t34, t13 ^ 2 + t14 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t68, t31, t9, -t10, 0, 0, t16, -t15, t31, t31 * t47 + t3, -t61 + (-t6 - t76) * t54, -t15 * t74 - t45 * t16, t1 * t45 + t2 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t55, 0, 0, 0, 0, 0, -t34, -t36, 0, -t34 * t45 + t36 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t57, 0, -t55 * pkin(8), -t57 * pkin(8), 0, 0, t36, -t34, 0, t22, -t23, -t34 * t74 - t45 * t36, t13 * t45 + t14 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t47, -0.2e1 * t74, 0, t54 ^ 2 * pkin(4) ^ 2 + t45 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t31, t3, -t4, -pkin(5) * t16, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t36, 0, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34, 0, t22, -t23, -pkin(5) * t36, t13 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t47, -t74, 0, t45 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;
