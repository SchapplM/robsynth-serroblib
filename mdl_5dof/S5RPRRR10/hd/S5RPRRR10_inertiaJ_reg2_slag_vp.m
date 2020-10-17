% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:44
% EndTime: 2019-12-31 19:10:48
% DurationCPUTime: 0.84s
% Computational Cost: add. (950->102), mult. (1886->201), div. (0->0), fcn. (2198->8), ass. (0->70)
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t49 = sin(qJ(3));
t71 = cos(qJ(3));
t28 = t71 * t45 + t49 * t46;
t82 = -0.2e1 * t28;
t63 = pkin(6) + qJ(2);
t34 = t63 * t46;
t58 = t63 * t45;
t15 = t49 * t34 + t71 * t58;
t81 = t15 ^ 2;
t26 = t49 * t45 - t71 * t46;
t23 = t26 ^ 2;
t80 = 0.2e1 * t26;
t38 = -t46 * pkin(2) - pkin(1);
t79 = 0.2e1 * t38;
t51 = cos(qJ(4));
t40 = -t51 * pkin(4) - pkin(3);
t78 = 0.2e1 * t40;
t77 = 0.2e1 * t46;
t76 = -pkin(8) - pkin(7);
t75 = t26 * pkin(4);
t47 = sin(qJ(5));
t74 = t47 * pkin(4);
t50 = cos(qJ(5));
t73 = t50 * pkin(4);
t14 = t26 * pkin(3) - t28 * pkin(7) + t38;
t48 = sin(qJ(4));
t17 = t71 * t34 - t49 * t58;
t65 = t51 * t17;
t5 = t65 + (-pkin(8) * t28 + t14) * t48;
t72 = t50 * t5;
t64 = t51 * t28;
t67 = t48 * t28;
t12 = -t47 * t67 + t50 * t64;
t31 = t47 * t48 - t50 * t51;
t70 = t12 * t31;
t33 = t47 * t51 + t50 * t48;
t69 = t33 * t26;
t68 = t48 * t26;
t66 = t48 * t51;
t41 = t45 ^ 2;
t42 = t46 ^ 2;
t62 = t41 + t42;
t43 = t48 ^ 2;
t44 = t51 ^ 2;
t61 = t43 + t44;
t60 = t26 * t82;
t59 = t48 * t64;
t6 = t51 * t14 - t48 * t17;
t4 = -pkin(8) * t64 + t6 + t75;
t1 = t50 * t4 - t47 * t5;
t57 = -pkin(3) * t28 - pkin(7) * t26;
t7 = t48 * t14 + t65;
t56 = t7 * t48 + t6 * t51;
t55 = -t6 * t48 + t7 * t51;
t36 = t76 * t51;
t35 = t76 * t48;
t30 = t33 ^ 2;
t29 = t31 ^ 2;
t24 = t28 ^ 2;
t22 = t51 * t26;
t20 = t47 * t35 - t50 * t36;
t19 = t50 * t35 + t47 * t36;
t18 = t31 * t26;
t10 = t33 * t28;
t9 = pkin(4) * t67 + t15;
t8 = t33 * t10;
t2 = t47 * t4 + t72;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t41, t45 * t77, 0, t42, 0, 0, pkin(1) * t77, -0.2e1 * pkin(1) * t45, 0.2e1 * t62 * qJ(2), t62 * qJ(2) ^ 2 + pkin(1) ^ 2, t24, t60, 0, t23, 0, 0, t26 * t79, t28 * t79, 0.2e1 * t15 * t28 - 0.2e1 * t17 * t26, t17 ^ 2 + t38 ^ 2 + t81, t44 * t24, -0.2e1 * t24 * t66, t64 * t80, t43 * t24, t48 * t60, t23, 0.2e1 * t15 * t67 + 0.2e1 * t6 * t26, 0.2e1 * t15 * t64 - 0.2e1 * t7 * t26, t56 * t82, t6 ^ 2 + t7 ^ 2 + t81, t12 ^ 2, -0.2e1 * t12 * t10, t12 * t80, t10 ^ 2, -t10 * t80, t23, 0.2e1 * t1 * t26 + 0.2e1 * t9 * t10, 0.2e1 * t9 * t12 - 0.2e1 * t2 * t26, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t45, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t26, t28, 0, t38, 0, 0, 0, 0, 0, 0, t22, -t68, -t61 * t28, t56, 0, 0, 0, 0, 0, 0, -t18, -t69, -t8 + t70, -t1 * t31 + t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, 0, -t15, -t17, 0, 0, t59, (-t43 + t44) * t28, t68, -t59, t22, 0, -t15 * t51 + t57 * t48, t15 * t48 + t57 * t51, t55, -t15 * pkin(3) + t55 * pkin(7), t12 * t33, -t8 - t70, t69, t10 * t31, -t18, 0, t40 * t10 + t19 * t26 + t9 * t31, t40 * t12 - t20 * t26 + t9 * t33, -t1 * t33 - t20 * t10 - t19 * t12 - t2 * t31, t1 * t19 + t2 * t20 + t9 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t19 + t33 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, 0.2e1 * t66, 0, t44, 0, 0, 0.2e1 * pkin(3) * t51, -0.2e1 * pkin(3) * t48, 0.2e1 * t61 * pkin(7), t61 * pkin(7) ^ 2 + pkin(3) ^ 2, t30, -0.2e1 * t33 * t31, 0, t29, 0, 0, t31 * t78, t33 * t78, -0.2e1 * t19 * t33 - 0.2e1 * t20 * t31, t19 ^ 2 + t20 ^ 2 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t67, t26, t6, -t7, 0, 0, 0, 0, t12, 0, -t10, t26, t26 * t73 + t1, -t72 + (-t4 - t75) * t47, (-t10 * t47 - t12 * t50) * pkin(4), (t1 * t50 + t2 * t47) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, 0, (-t31 * t50 + t33 * t47) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, t51, 0, -t48 * pkin(7), -t51 * pkin(7), 0, 0, 0, 0, t33, 0, -t31, 0, t19, -t20, (-t31 * t47 - t33 * t50) * pkin(4), (t19 * t50 + t20 * t47) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t73, -0.2e1 * t74, 0, (t47 ^ 2 + t50 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, t26, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t31, 0, t19, -t20, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t73, -t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
