% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:19:16
% EndTime: 2019-05-06 17:19:18
% DurationCPUTime: 0.62s
% Computational Cost: add. (1084->94), mult. (2025->170), div. (0->0), fcn. (2441->8), ass. (0->75)
t51 = sin(pkin(10));
t52 = cos(pkin(10));
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t34 = t51 * t58 + t52 * t55;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t60 = t51 * t55 - t52 * t58;
t21 = t57 * t34 - t54 * t60;
t83 = -0.2e1 * t21;
t47 = -t58 * pkin(2) - pkin(1);
t28 = t60 * pkin(3) + t47;
t82 = 0.2e1 * t28;
t81 = 0.2e1 * t58;
t80 = t51 * pkin(2);
t53 = sin(qJ(5));
t79 = t53 * pkin(5);
t56 = cos(qJ(5));
t78 = t56 * pkin(5);
t77 = t56 * pkin(9);
t45 = t52 * pkin(2) + pkin(3);
t31 = t57 * t45 - t54 * t80;
t29 = -pkin(4) - t31;
t76 = pkin(4) - t29;
t69 = -qJ(3) - pkin(7);
t40 = t69 * t55;
t42 = t69 * t58;
t23 = t52 * t40 + t51 * t42;
t14 = -t34 * pkin(8) + t23;
t24 = t51 * t40 - t52 * t42;
t15 = -t60 * pkin(8) + t24;
t10 = -t57 * t14 + t54 * t15;
t75 = t10 * t56;
t20 = t54 * t34 + t57 * t60;
t17 = t53 * t20;
t74 = t53 * t21;
t73 = t53 * t56;
t11 = t54 * t14 + t57 * t15;
t72 = t56 * t11;
t71 = t56 * t21;
t32 = -t54 * t45 - t57 * t80;
t30 = pkin(9) - t32;
t70 = t56 * t30;
t49 = t53 ^ 2;
t50 = t56 ^ 2;
t68 = t49 + t50;
t48 = t56 * qJ(6);
t67 = t20 * t83;
t9 = t20 * pkin(4) - t21 * pkin(9) + t28;
t4 = -t53 * t11 + t56 * t9;
t1 = t20 * pkin(5) - t21 * t48 + t4;
t3 = t72 + (-qJ(6) * t21 + t9) * t53;
t66 = -t1 * t53 + t3 * t56;
t65 = -pkin(4) * t21 - pkin(9) * t20;
t64 = t1 * t56 + t3 * t53;
t63 = -t20 * t30 + t21 * t29;
t25 = (-qJ(6) - t30) * t53;
t26 = t48 + t70;
t62 = t25 * t56 + t26 * t53;
t39 = (-qJ(6) - pkin(9)) * t53;
t41 = t48 + t77;
t61 = t56 * t39 + t53 * t41;
t46 = -pkin(4) - t78;
t44 = 0.2e1 * t73;
t38 = t41 * t56;
t27 = t29 - t78;
t22 = t26 * t56;
t19 = t21 ^ 2;
t18 = t56 * t20;
t16 = t53 * t71;
t12 = (-t49 + t50) * t21;
t8 = t10 * t53;
t6 = pkin(5) * t74 + t10;
t5 = t53 * t9 + t72;
t2 = [1, 0, 0, t55 ^ 2, t55 * t81, 0, 0, 0, pkin(1) * t81, -0.2e1 * pkin(1) * t55, -0.2e1 * t23 * t34 - 0.2e1 * t24 * t60, t23 ^ 2 + t24 ^ 2 + t47 ^ 2, t19, t67, 0, 0, 0, t20 * t82, t21 * t82, t50 * t19, -0.2e1 * t19 * t73, 0.2e1 * t20 * t71, t53 * t67, t20 ^ 2, 0.2e1 * t10 * t74 + 0.2e1 * t4 * t20, 0.2e1 * t10 * t71 - 0.2e1 * t5 * t20, t64 * t83, t1 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t55, t58, 0, -t55 * pkin(7), -t58 * pkin(7) (-t52 * t34 - t51 * t60) * pkin(2) (t23 * t52 + t24 * t51) * pkin(2), 0, 0, t21, -t20, 0, -t10, -t11, t16, t12, t17, t18, 0, t63 * t53 - t75, t63 * t56 + t8, -t62 * t21 + t66, t1 * t25 + t3 * t26 + t6 * t27; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t51 ^ 2 + t52 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t31, 0.2e1 * t32, t49, t44, 0, 0, 0, -0.2e1 * t29 * t56, 0.2e1 * t29 * t53, -0.2e1 * t25 * t53 + 0.2e1 * t22, t25 ^ 2 + t26 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, t18, -t17, -t68 * t21, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t10, -t11, t16, t12, t17, t18, 0, t65 * t53 - t75, t65 * t56 + t8, -t61 * t21 + t66, t1 * t39 + t3 * t41 + t6 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, t32, t49, t44, 0, 0, 0, t76 * t56, -t76 * t53, t22 + t38 + (-t25 - t39) * t53, t25 * t39 + t26 * t41 + t27 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t49, t44, 0, 0, 0, 0.2e1 * pkin(4) * t56, -0.2e1 * pkin(4) * t53, -0.2e1 * t39 * t53 + 0.2e1 * t38, t39 ^ 2 + t41 ^ 2 + t46 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t74, t20, t4, -t5, -pkin(5) * t71, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t56, 0, -t53 * t30, -t70, -t79, t25 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t53, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t56, 0, -t53 * pkin(9), -t77, -t79, t39 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t2;
