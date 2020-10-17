% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:13
% EndTime: 2019-12-31 20:18:16
% DurationCPUTime: 0.86s
% Computational Cost: add. (999->80), mult. (1922->169), div. (0->0), fcn. (2252->8), ass. (0->65)
t43 = cos(pkin(9));
t70 = t43 * pkin(2);
t34 = pkin(3) + t70;
t45 = sin(qJ(4));
t67 = cos(qJ(4));
t42 = sin(pkin(9));
t71 = t42 * pkin(2);
t28 = t45 * t34 + t67 * t71;
t26 = pkin(8) + t28;
t44 = sin(qJ(5));
t38 = t44 ^ 2;
t47 = cos(qJ(5));
t40 = t47 ^ 2;
t60 = t38 + t40;
t62 = t60 * t26;
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t30 = t42 * t46 - t43 * t48;
t32 = t42 * t48 + t43 * t46;
t18 = -t45 * t30 + t32 * t67;
t77 = -0.2e1 * t18;
t63 = -qJ(3) - pkin(6);
t56 = t63 * t46;
t57 = t63 * t48;
t20 = t42 * t56 - t43 * t57;
t11 = -t30 * pkin(7) + t20;
t19 = t42 * t57 + t43 * t56;
t52 = -t32 * pkin(7) + t19;
t6 = t45 * t11 - t52 * t67;
t76 = t6 ^ 2;
t16 = t30 * t67 + t45 * t32;
t75 = t16 ^ 2;
t35 = -t48 * pkin(2) - pkin(1);
t23 = t30 * pkin(3) + t35;
t74 = 0.2e1 * t23;
t73 = 0.2e1 * t35;
t72 = 0.2e1 * t48;
t69 = t6 * t47;
t27 = t67 * t34 - t45 * t71;
t25 = -pkin(4) - t27;
t68 = pkin(4) - t25;
t13 = t44 * t16;
t66 = t44 * t18;
t65 = t44 * t47;
t64 = t47 * t18;
t61 = t60 * pkin(8);
t39 = t46 ^ 2;
t41 = t48 ^ 2;
t59 = t39 + t41;
t58 = t16 * t77;
t55 = -pkin(4) * t18 - pkin(8) * t16;
t5 = t16 * pkin(4) - t18 * pkin(8) + t23;
t8 = t11 * t67 + t45 * t52;
t2 = -t44 * t8 + t47 * t5;
t3 = t44 * t5 + t47 * t8;
t54 = t2 * t47 + t3 * t44;
t1 = -t2 * t44 + t3 * t47;
t53 = -t16 * t26 + t18 * t25;
t33 = 0.2e1 * t65;
t15 = t18 ^ 2;
t14 = t47 * t16;
t12 = t44 * t64;
t9 = (-t38 + t40) * t18;
t4 = t6 * t44;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t39, t46 * t72, 0, t41, 0, 0, pkin(1) * t72, -0.2e1 * pkin(1) * t46, 0.2e1 * t59 * pkin(6), pkin(6) ^ 2 * t59 + pkin(1) ^ 2, t32 ^ 2, -0.2e1 * t32 * t30, 0, t30 ^ 2, 0, 0, t30 * t73, t32 * t73, -0.2e1 * t19 * t32 - 0.2e1 * t20 * t30, t19 ^ 2 + t20 ^ 2 + t35 ^ 2, t15, t58, 0, t75, 0, 0, t16 * t74, t18 * t74, -0.2e1 * t8 * t16 + 0.2e1 * t6 * t18, t23 ^ 2 + t8 ^ 2 + t76, t40 * t15, -0.2e1 * t15 * t65, 0.2e1 * t16 * t64, t38 * t15, t44 * t58, t75, 0.2e1 * t2 * t16 + 0.2e1 * t6 * t66, -0.2e1 * t3 * t16 + 0.2e1 * t6 * t64, t54 * t77, t2 ^ 2 + t3 ^ 2 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t48, 0, -t46 * pkin(6), -t48 * pkin(6), 0, 0, 0, 0, t32, 0, -t30, 0, t19, -t20, (-t30 * t42 - t32 * t43) * pkin(2), (t19 * t43 + t20 * t42) * pkin(2), 0, 0, t18, 0, -t16, 0, -t6, -t8, -t28 * t16 - t27 * t18, -t6 * t27 + t8 * t28, t12, t9, t13, -t12, t14, 0, t44 * t53 - t69, t47 * t53 + t4, t1, t1 * t26 + t6 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t70, -0.2e1 * t71, 0, (t42 ^ 2 + t43 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, -0.2e1 * t28, 0, t27 ^ 2 + t28 ^ 2, t38, t33, 0, t40, 0, 0, -0.2e1 * t25 * t47, 0.2e1 * t25 * t44, 0.2e1 * t62, t26 ^ 2 * t60 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t32, 0, t35, 0, 0, 0, 0, 0, 0, t16, t18, 0, t23, 0, 0, 0, 0, 0, 0, t14, -t13, -t60 * t18, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, -t6, -t8, 0, 0, t12, t9, t13, -t12, t14, 0, t44 * t55 - t69, t47 * t55 + t4, t1, -t6 * pkin(4) + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t28, 0, 0, t38, t33, 0, t40, 0, 0, t68 * t47, -t68 * t44, t61 + t62, -t25 * pkin(4) + pkin(8) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t38, t33, 0, t40, 0, 0, 0.2e1 * pkin(4) * t47, -0.2e1 * pkin(4) * t44, 0.2e1 * t61, pkin(8) ^ 2 * t60 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t66, t16, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t47, 0, -t44 * t26, -t47 * t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t44, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t47, 0, -t44 * pkin(8), -t47 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t7;
