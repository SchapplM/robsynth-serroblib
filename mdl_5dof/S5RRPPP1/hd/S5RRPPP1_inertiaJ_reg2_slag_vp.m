% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:30
% EndTime: 2019-12-31 19:24:34
% DurationCPUTime: 0.98s
% Computational Cost: add. (732->127), mult. (1845->247), div. (0->0), fcn. (1921->6), ass. (0->78)
t58 = sin(qJ(2));
t59 = cos(qJ(2));
t55 = sin(pkin(5));
t77 = qJ(3) * t55;
t34 = -t59 * pkin(2) - t58 * t77 - pkin(1);
t57 = cos(pkin(5));
t66 = qJ(3) * t57 + pkin(7);
t35 = t66 * t58;
t56 = cos(pkin(8));
t91 = (t34 * t55 - t35 * t57) * t56;
t54 = sin(pkin(8));
t81 = t57 * t59;
t26 = t58 * t54 - t56 * t81;
t28 = t54 * t81 + t58 * t56;
t61 = (t26 * t54 - t28 * t56) * t55;
t24 = t26 ^ 2;
t25 = t28 ^ 2;
t90 = 0.2e1 * t55;
t89 = 0.2e1 * t59;
t88 = t55 * pkin(2);
t87 = t26 * t28;
t50 = t55 ^ 2;
t86 = t50 * t56;
t85 = t50 * t59;
t84 = t54 * t57;
t47 = t55 * t54;
t48 = t55 * t56;
t83 = t55 * t59;
t82 = t56 * t57;
t80 = pkin(3) + qJ(5);
t36 = t66 * t59;
t32 = t54 * t36;
t79 = pkin(3) * t83 + t32;
t31 = pkin(2) * t84 + t56 * t77;
t52 = t58 ^ 2;
t53 = t59 ^ 2;
t78 = t52 + t53;
t76 = -0.2e1 * t87;
t75 = t26 * t48;
t17 = t28 * t47;
t74 = t28 * t83;
t73 = t54 * t86;
t72 = t57 * t47;
t71 = t55 * t82;
t70 = t26 * t83;
t69 = t55 * t81;
t8 = t34 * t47 - t35 * t84 + t56 * t36;
t68 = -pkin(2) * t56 - pkin(3);
t67 = -qJ(4) * t54 - pkin(2);
t9 = t57 * t34 + t55 * t35;
t18 = -t57 * qJ(4) - t31;
t12 = t57 * t26 + t56 * t85;
t63 = -t28 * t57 + t54 * t85;
t62 = -t28 * qJ(4) + t9;
t5 = qJ(4) * t83 - t8;
t51 = t57 ^ 2;
t49 = t50 * t53;
t46 = t50 * t56 ^ 2;
t45 = t50 * t54 ^ 2;
t41 = t54 * t77;
t39 = -0.2e1 * t71;
t38 = 0.2e1 * t72;
t37 = 0.2e1 * t73;
t30 = pkin(2) * t82 - t41;
t20 = (-pkin(3) * t56 + t67) * t55;
t19 = t68 * t57 + t41;
t16 = -0.2e1 * t74;
t15 = -0.2e1 * t70;
t14 = (-t80 * t56 + t67) * t55;
t13 = pkin(4) * t48 - t18;
t10 = pkin(4) * t47 + t41 + (-qJ(5) + t68) * t57;
t7 = -t32 + t91;
t6 = t79 - t91;
t4 = t26 * pkin(3) + t62;
t3 = -t26 * pkin(4) - t5;
t2 = t80 * t26 + t62;
t1 = t35 * t82 + t28 * pkin(4) + (qJ(5) * t59 - t34 * t56) * t55 + t79;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t52, t58 * t89, 0, t53, 0, 0, pkin(1) * t89, -0.2e1 * pkin(1) * t58, 0.2e1 * t78 * pkin(7), t78 * pkin(7) ^ 2 + pkin(1) ^ 2, t25, t76, t16, t24, 0.2e1 * t70, t49, 0.2e1 * t9 * t26 - 0.2e1 * t7 * t83, 0.2e1 * t9 * t28 + 0.2e1 * t8 * t83, -0.2e1 * t8 * t26 - 0.2e1 * t7 * t28, t7 ^ 2 + t8 ^ 2 + t9 ^ 2, t49, 0.2e1 * t74, t15, t25, t76, t24, 0.2e1 * t5 * t26 + 0.2e1 * t6 * t28, -0.2e1 * t4 * t26 - 0.2e1 * t6 * t83, -0.2e1 * t4 * t28 + 0.2e1 * t5 * t83, t4 ^ 2 + t5 ^ 2 + t6 ^ 2, t49, t15, t16, t24, 0.2e1 * t87, t25, 0.2e1 * t1 * t28 - 0.2e1 * t3 * t26, -0.2e1 * t2 * t28 - 0.2e1 * t3 * t83, 0.2e1 * t1 * t83 + 0.2e1 * t2 * t26, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t59, 0, -t58 * pkin(7), -t59 * pkin(7), 0, 0, t17, -t61, -t63, -t75, -t12, -t69, t7 * t57 + (-pkin(2) * t26 - t30 * t59 - t56 * t9) * t55, -t8 * t57 + (-pkin(2) * t28 + t31 * t59 + t54 * t9) * t55, -t31 * t26 - t30 * t28 + (-t54 * t7 + t56 * t8) * t55, t7 * t30 + t8 * t31 - t9 * t88, -t69, t63, t12, t17, -t61, -t75, t18 * t26 + t19 * t28 + (-t5 * t56 + t54 * t6) * t55, -t20 * t26 + t6 * t57 + (-t19 * t59 + t4 * t56) * t55, -t20 * t28 - t5 * t57 + (t18 * t59 - t4 * t54) * t55, t5 * t18 + t6 * t19 + t4 * t20, -t69, t12, -t63, -t75, t61, t17, t10 * t28 - t13 * t26 + (t1 * t54 + t3 * t56) * t55, -t14 * t28 + t3 * t57 + (-t13 * t59 - t2 * t54) * t55, -t1 * t57 + t14 * t26 + (t10 * t59 - t2 * t56) * t55, t1 * t10 + t3 * t13 + t2 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45, t37, t38, t46, 0.2e1 * t71, t51, 0.2e1 * pkin(2) * t86 + 0.2e1 * t30 * t57, -0.2e1 * t50 * pkin(2) * t54 - 0.2e1 * t31 * t57, (-t30 * t54 + t31 * t56) * t90, t50 * pkin(2) ^ 2 + t30 ^ 2 + t31 ^ 2, t51, -0.2e1 * t72, t39, t45, t37, t46, (-t18 * t56 + t19 * t54) * t90, 0.2e1 * t19 * t57 + 0.2e1 * t20 * t48, -0.2e1 * t18 * t57 - 0.2e1 * t20 * t47, t18 ^ 2 + t19 ^ 2 + t20 ^ 2, t51, t39, t38, t46, -0.2e1 * t73, t45, (t10 * t54 + t13 * t56) * t90, 0.2e1 * t13 * t57 - 0.2e1 * t14 * t47, -0.2e1 * t10 * t57 - 0.2e1 * t14 * t48, t10 ^ 2 + t13 ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, t9, 0, 0, 0, 0, 0, 0, 0, -t26, -t28, t4, 0, 0, 0, 0, 0, 0, 0, -t28, t26, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, -t88, 0, 0, 0, 0, 0, 0, 0, t48, -t47, t20, 0, 0, 0, 0, 0, 0, 0, -t47, -t48, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t83, 0, t6, 0, 0, 0, 0, 0, 0, t28, 0, t83, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t57, 0, t19, 0, 0, 0, 0, 0, 0, t47, 0, -t57, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t83, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t57, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t11;
