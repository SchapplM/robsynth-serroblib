% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR15_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:09
% EndTime: 2019-12-31 20:43:12
% DurationCPUTime: 0.86s
% Computational Cost: add. (599->102), mult. (1087->179), div. (0->0), fcn. (1116->6), ass. (0->77)
t48 = sin(qJ(2));
t42 = t48 ^ 2;
t51 = cos(qJ(2));
t44 = t51 ^ 2;
t88 = t42 + t44;
t46 = sin(qJ(5));
t47 = sin(qJ(4));
t49 = cos(qJ(5));
t50 = cos(qJ(4));
t18 = t46 * t50 + t49 * t47;
t20 = -t46 * t47 + t49 * t50;
t87 = (t18 * t46 + t20 * t49) * pkin(4);
t17 = t20 ^ 2;
t85 = t18 ^ 2;
t86 = t17 + t85;
t32 = t47 * pkin(4) + qJ(3);
t84 = 0.2e1 * t32;
t83 = -0.2e1 * t48;
t82 = 0.2e1 * t48;
t81 = 0.2e1 * t51;
t80 = 0.2e1 * qJ(3);
t52 = -pkin(2) - pkin(7);
t79 = t46 * pkin(4);
t78 = t48 * pkin(4);
t77 = t49 * pkin(4);
t61 = -t48 * qJ(3) - pkin(1);
t16 = t52 * t51 + t61;
t62 = pkin(8) * t51 - t16;
t37 = t48 * pkin(6);
t27 = t48 * pkin(3) + t37;
t72 = t47 * t27;
t6 = -t62 * t50 + t72;
t76 = t49 * t6;
t67 = t50 * t51;
t70 = t47 * t51;
t11 = t46 * t70 - t49 * t67;
t75 = t18 * t11;
t74 = t18 * t48;
t12 = t18 * t51;
t73 = t20 * t12;
t71 = t47 * t48;
t69 = t48 * t51;
t68 = t50 * t47;
t66 = t88 * pkin(6) ^ 2;
t39 = t51 * pkin(6);
t28 = t51 * pkin(3) + t39;
t41 = t47 ^ 2;
t43 = t50 ^ 2;
t30 = t41 + t43;
t65 = t51 * qJ(3);
t64 = -0.2e1 * t69;
t63 = t47 * t67;
t21 = t50 * t27;
t5 = t62 * t47 + t21 + t78;
t1 = -t46 * t6 + t49 * t5;
t2 = t46 * t5 + t76;
t60 = t1 * t20 + t2 * t18;
t7 = -t47 * t16 + t21;
t8 = t50 * t16 + t72;
t3 = t8 * t47 + t7 * t50;
t24 = (-pkin(8) + t52) * t47;
t35 = t50 * t52;
t25 = -t50 * pkin(8) + t35;
t10 = t49 * t24 + t46 * t25;
t9 = -t46 * t24 + t49 * t25;
t59 = t10 * t18 + t9 * t20;
t58 = -t48 * pkin(2) + t65;
t56 = t48 * t52 + t65;
t53 = qJ(3) ^ 2;
t34 = t50 * t48;
t31 = 0.2e1 * t69;
t26 = -t51 * pkin(2) + t61;
t23 = 0.2e1 * t88 * pkin(6);
t22 = t30 * t52;
t15 = t20 * t48;
t14 = pkin(4) * t67 + t28;
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t42, t31, 0, t44, 0, 0, pkin(1) * t81, pkin(1) * t83, t23, pkin(1) ^ 2 + t66, 0, 0, 0, t42, t31, t44, t23, t26 * t81, t26 * t83, t26 ^ 2 + t66, t41 * t44, 0.2e1 * t44 * t68, t47 * t64, t43 * t44, t50 * t64, t42, 0.2e1 * t28 * t67 + 0.2e1 * t7 * t48, -0.2e1 * t28 * t70 - 0.2e1 * t8 * t48, (t47 * t7 - t50 * t8) * t81, t28 ^ 2 + t7 ^ 2 + t8 ^ 2, t12 ^ 2, -0.2e1 * t12 * t11, -t12 * t82, t11 ^ 2, t11 * t82, t42, 0.2e1 * t1 * t48 - 0.2e1 * t14 * t11, -0.2e1 * t14 * t12 - 0.2e1 * t2 * t48, 0.2e1 * t1 * t12 + 0.2e1 * t2 * t11, t1 ^ 2 + t14 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, t51, 0, -t37, -t39, 0, 0, 0, -t48, -t51, 0, 0, 0, t58, t37, t39, t58 * pkin(6), -t63, (t41 - t43) * t51, t34, t63, -t71, 0, t28 * t47 + t56 * t50, t28 * t50 - t56 * t47, -t3, t28 * qJ(3) + t3 * t52, -t73, t20 * t11 + t12 * t18, t15, -t75, -t74, 0, -t32 * t11 + t14 * t18 + t9 * t48, -t10 * t48 - t32 * t12 + t14 * t20, t10 * t11 + t9 * t12 - t60, t1 * t9 + t2 * t10 + t14 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t80, pkin(2) ^ 2 + t53, t43, -0.2e1 * t68, 0, t41, 0, 0, t47 * t80, t50 * t80, -0.2e1 * t22, t30 * t52 ^ 2 + t53, t17, -0.2e1 * t20 * t18, 0, t85, 0, 0, t18 * t84, t20 * t84, -0.2e1 * t59, t10 ^ 2 + t32 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, t37, 0, 0, 0, 0, 0, 0, t34, -t71, 0, t3, 0, 0, 0, 0, 0, 0, t15, -t74, t73 + t75, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t30, t22, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, -t67, t48, t7, -t8, 0, 0, 0, 0, -t12, 0, t11, t48, t48 * t77 + t1, -t76 + (-t5 - t78) * t46, (t11 * t46 + t12 * t49) * pkin(4), (t1 * t49 + t2 * t46) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t47, 0, t35, -t47 * t52, 0, 0, 0, 0, t20, 0, -t18, 0, t9, -t10, -t87, (t10 * t46 + t49 * t9) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t77, -0.2e1 * t79, 0, (t46 ^ 2 + t49 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, t11, t48, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t77, -t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t4;
