% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:32:01
% EndTime: 2019-05-06 19:32:04
% DurationCPUTime: 0.74s
% Computational Cost: add. (1349->86), mult. (2584->152), div. (0->0), fcn. (3282->10), ass. (0->76)
t54 = sin(pkin(11));
t55 = cos(pkin(11));
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t36 = -t54 * t59 + t55 * t62;
t49 = -t62 * pkin(2) - pkin(1);
t31 = -t36 * pkin(3) + t49;
t37 = t54 * t62 + t55 * t59;
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t68 = t61 * t36 - t58 * t37;
t18 = -pkin(4) * t68 + t31;
t85 = 0.2e1 * t18;
t27 = t58 * t36 + t61 * t37;
t84 = 0.2e1 * t27;
t83 = 0.2e1 * t62;
t82 = pkin(2) * t54;
t56 = sin(qJ(6));
t81 = pkin(5) * t56;
t57 = sin(qJ(5));
t72 = -qJ(3) - pkin(7);
t41 = t72 * t59;
t42 = t72 * t62;
t28 = t55 * t41 + t54 * t42;
t24 = -t37 * pkin(8) + t28;
t29 = t54 * t41 - t55 * t42;
t25 = t36 * pkin(8) + t29;
t10 = t61 * t24 - t58 * t25;
t64 = -t27 * pkin(9) + t10;
t78 = cos(qJ(5));
t11 = -t58 * t24 - t61 * t25;
t9 = pkin(9) * t68 - t11;
t4 = t57 * t9 - t64 * t78;
t60 = cos(qJ(6));
t80 = t4 * t60;
t79 = t57 * pkin(4);
t46 = t55 * pkin(2) + pkin(3);
t33 = t61 * t46 - t58 * t82;
t32 = pkin(4) + t33;
t34 = t58 * t46 + t61 * t82;
t69 = -t78 * t32 + t57 * t34;
t20 = -pkin(5) + t69;
t77 = t20 * t60;
t50 = t78 * pkin(4);
t48 = -t50 - pkin(5);
t76 = t48 * t60;
t16 = t57 * t27 - t68 * t78;
t13 = t56 * t16;
t17 = t27 * t78 + t57 * t68;
t75 = t56 * t17;
t74 = t56 * t60;
t73 = t60 * t17;
t71 = -0.2e1 * t17 * t16;
t70 = t78 * t34;
t67 = -pkin(5) * t17 - pkin(10) * t16;
t23 = -t57 * t32 - t70;
t21 = pkin(10) - t23;
t66 = -t16 * t21 + t17 * t20;
t47 = pkin(10) + t79;
t65 = -t16 * t47 + t17 * t48;
t53 = t60 ^ 2;
t52 = t56 ^ 2;
t51 = pkin(5) * t60;
t45 = 0.2e1 * t74;
t44 = t48 * t56;
t19 = t20 * t56;
t15 = t17 ^ 2;
t14 = t60 * t16;
t12 = t56 * t73;
t7 = (-t52 + t53) * t17;
t6 = t16 * pkin(5) - t17 * pkin(10) + t18;
t5 = t57 * t64 + t78 * t9;
t3 = t4 * t56;
t2 = t60 * t5 + t56 * t6;
t1 = -t56 * t5 + t60 * t6;
t8 = [1, 0, 0, t59 ^ 2, t59 * t83, 0, 0, 0, pkin(1) * t83, -0.2e1 * pkin(1) * t59, -0.2e1 * t28 * t37 + 0.2e1 * t29 * t36, t28 ^ 2 + t29 ^ 2 + t49 ^ 2, t27 ^ 2, t68 * t84, 0, 0, 0, -0.2e1 * t31 * t68, t31 * t84, t15, t71, 0, 0, 0, t16 * t85, t17 * t85, t53 * t15, -0.2e1 * t15 * t74, 0.2e1 * t16 * t73, t56 * t71, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t4 * t75, -0.2e1 * t2 * t16 + 0.2e1 * t4 * t73; 0, 0, 0, 0, 0, t59, t62, 0, -t59 * pkin(7), -t62 * pkin(7) (t36 * t54 - t37 * t55) * pkin(2) (t28 * t55 + t29 * t54) * pkin(2), 0, 0, t27, t68, 0, t10, t11, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t56 * t66 - t80, t60 * t66 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t54 ^ 2 + t55 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t34, 0, 0, 0, 0, 1, -0.2e1 * t69, 0.2e1 * t23, t52, t45, 0, 0, 0, -0.2e1 * t77, 0.2e1 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, -t68, t27, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, t14, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t68, 0, t10, t11, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t56 * t65 - t80, t60 * t65 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t33, -t34, 0, 0, 0, 0, 1, t50 - t69, -t70 + (-pkin(4) - t32) * t57, t52, t45, 0, 0, 0 (-t20 - t48) * t60, t44 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t79, t52, t45, 0, 0, 0, -0.2e1 * t76, 0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t56 * t67 - t80, t60 * t67 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t69, t23, t52, t45, 0, 0, 0, t51 - t77, t19 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t79, t52, t45, 0, 0, 0, t51 - t76, t44 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t52, t45, 0, 0, 0, 0.2e1 * t51, -0.2e1 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t75, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t60, 0, -t56 * t21, -t60 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t60, 0, -t56 * t47, -t60 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t60, 0, -t56 * pkin(10), -t60 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
