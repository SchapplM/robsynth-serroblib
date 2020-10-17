% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:36
% EndTime: 2020-01-03 11:57:38
% DurationCPUTime: 0.58s
% Computational Cost: add. (420->74), mult. (772->118), div. (0->0), fcn. (696->8), ass. (0->67)
t45 = sin(pkin(9));
t41 = t45 ^ 2;
t47 = cos(pkin(9));
t42 = t47 ^ 2;
t77 = t41 + t42;
t76 = -0.2e1 * t45;
t75 = 0.2e1 * t45;
t74 = -0.2e1 * t47;
t52 = cos(qJ(2));
t40 = t52 * pkin(1);
t39 = t40 + pkin(2);
t46 = sin(pkin(8));
t48 = cos(pkin(8));
t50 = sin(qJ(2));
t70 = t50 * pkin(1);
t58 = t48 * t70;
t18 = t39 * t46 + t58;
t15 = qJ(4) + t18;
t49 = sin(qJ(5));
t35 = t49 * t47;
t51 = cos(qJ(5));
t55 = -pkin(4) * t47 - pkin(7) * t45 - pkin(3);
t60 = -t48 * t39 + t46 * t70;
t8 = t55 + t60;
t2 = -t15 * t35 + t51 * t8;
t71 = t48 * pkin(2);
t19 = t55 - t71;
t72 = t46 * pkin(2);
t33 = qJ(4) + t72;
t5 = t19 * t51 - t33 * t35;
t73 = -t2 - t5;
t64 = t51 * t47;
t3 = t15 * t64 + t49 * t8;
t66 = t41 * t51;
t69 = t15 * t66 + t3 * t47;
t6 = t19 * t49 + t33 * t64;
t68 = t33 * t66 + t6 * t47;
t24 = t41 * t33;
t67 = t41 * t49;
t25 = t42 * t33;
t65 = t49 * t45;
t63 = t77 * t15;
t16 = -pkin(3) + t60;
t38 = -pkin(3) - t71;
t62 = t16 + t38;
t61 = t25 + t24;
t43 = t49 ^ 2;
t44 = t51 ^ 2;
t59 = t43 + t44;
t30 = t47 * t75;
t57 = t2 * t51 + t3 * t49;
t56 = t6 * t49 + t5 * t51;
t37 = t51 * t45;
t36 = t44 * t41;
t34 = t43 * t41;
t32 = t33 ^ 2;
t28 = -0.2e1 * t49 * t66;
t27 = t64 * t76;
t26 = t49 * t30;
t23 = t41 * t32;
t21 = t33 * t67;
t20 = t59 * t45;
t14 = t15 ^ 2;
t11 = t41 * t14;
t9 = t15 * t67;
t7 = t15 * t24;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t70, 0, (t50 ^ 2 + t52 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t60, -0.2e1 * t18, 0, t18 ^ 2 + t60 ^ 2, t41, t30, 0, t42, 0, 0, t16 * t74, t16 * t75, 0.2e1 * t63, t14 * t42 + t16 ^ 2 + t11, t36, t28, t27, t34, t26, t42, -0.2e1 * t2 * t47 + 0.2e1 * t9, 0.2e1 * t69, t57 * t76, t2 ^ 2 + t3 ^ 2 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t40, -t70, 0, 0, 0, 0, 0, 0, 0, 1, -t60 + t71, -t58 + (-pkin(2) - t39) * t46, 0, (t18 * t46 - t48 * t60) * pkin(2), t41, t30, 0, t42, 0, 0, -t62 * t47, t62 * t45, t61 + t63, t15 * t25 + t16 * t38 + t7, t36, t28, t27, t34, t26, t42, t47 * t73 + t21 + t9, t68 + t69, (t73 * t51 + (-t3 - t6) * t49) * t45, t2 * t5 + t3 * t6 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t72, 0, (t46 ^ 2 + t48 ^ 2) * pkin(2) ^ 2, t41, t30, 0, t42, 0, 0, t38 * t74, t38 * t75, 0.2e1 * t61, t32 * t42 + t38 ^ 2 + t23, t36, t28, t27, t34, t26, t42, -0.2e1 * t47 * t5 + 0.2e1 * t21, 0.2e1 * t68, t56 * t76, t5 ^ 2 + t6 ^ 2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t15 * t47 - t2 * t49 + t3 * t51) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t33 * t47 - t49 * t5 + t51 * t6) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 + t34 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t45, 0, t16, 0, 0, 0, 0, 0, 0, -t64, t35, -t20, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t45, 0, t38, 0, 0, 0, 0, 0, 0, -t64, t35, -t20, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t65, -t47, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t65, -t47, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
