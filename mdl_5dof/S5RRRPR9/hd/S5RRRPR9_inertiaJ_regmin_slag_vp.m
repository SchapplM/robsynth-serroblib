% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:38
% EndTime: 2019-12-31 21:24:40
% DurationCPUTime: 0.43s
% Computational Cost: add. (570->89), mult. (1209->175), div. (0->0), fcn. (1347->8), ass. (0->65)
t52 = sin(pkin(9));
t53 = cos(pkin(9));
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t34 = t52 * t58 + t53 * t55;
t56 = sin(qJ(2));
t28 = t34 * t56;
t66 = t58 * t56;
t69 = t55 * t56;
t29 = -t52 * t69 + t53 * t66;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t11 = -t54 * t28 + t57 * t29;
t76 = -0.2e1 * t11;
t33 = -t52 * t55 + t53 * t58;
t46 = -t58 * pkin(3) - pkin(2);
t27 = -t33 * pkin(4) + t46;
t75 = 0.2e1 * t27;
t74 = -0.2e1 * t56;
t59 = cos(qJ(2));
t73 = 0.2e1 * t59;
t72 = pkin(2) * t58;
t71 = pkin(3) * t52;
t70 = pkin(6) * t55;
t68 = t55 * t58;
t67 = t55 * t59;
t65 = t58 * t59;
t64 = -qJ(4) - pkin(7);
t40 = -t59 * pkin(2) - t56 * pkin(7) - pkin(1);
t35 = t58 * t40;
t63 = qJ(4) * t56;
t18 = -t58 * t63 + t35 + (-pkin(3) - t70) * t59;
t61 = pkin(6) * t65;
t21 = t61 + (t40 - t63) * t55;
t9 = t52 * t18 + t53 * t21;
t41 = t64 * t55;
t42 = t64 * t58;
t23 = t52 * t41 - t53 * t42;
t47 = t56 * pkin(6);
t39 = pkin(3) * t69 + t47;
t62 = t56 * t73;
t8 = t53 * t18 - t52 * t21;
t4 = -t59 * pkin(4) - t29 * pkin(8) + t8;
t7 = -t28 * pkin(8) + t9;
t1 = t57 * t4 - t54 * t7;
t22 = t53 * t41 + t52 * t42;
t2 = t54 * t4 + t57 * t7;
t51 = t59 ^ 2;
t50 = t58 ^ 2;
t49 = t56 ^ 2;
t48 = t55 ^ 2;
t45 = t53 * pkin(3) + pkin(4);
t31 = t54 * t45 + t57 * t71;
t30 = t57 * t45 - t54 * t71;
t26 = t55 * t40 + t61;
t25 = -pkin(6) * t67 + t35;
t19 = t28 * pkin(4) + t39;
t17 = t54 * t33 + t57 * t34;
t16 = -t57 * t33 + t54 * t34;
t13 = t33 * pkin(8) + t23;
t12 = -t34 * pkin(8) + t22;
t10 = t57 * t28 + t54 * t29;
t6 = t54 * t12 + t57 * t13;
t5 = t57 * t12 - t54 * t13;
t3 = [1, 0, 0, t49, t62, 0, 0, 0, pkin(1) * t73, pkin(1) * t74, t50 * t49, -0.2e1 * t49 * t68, t65 * t74, t55 * t62, t51, -0.2e1 * t25 * t59 + 0.2e1 * t49 * t70, 0.2e1 * t49 * pkin(6) * t58 + 0.2e1 * t26 * t59, -0.2e1 * t9 * t28 - 0.2e1 * t8 * t29, t39 ^ 2 + t8 ^ 2 + t9 ^ 2, t11 ^ 2, t10 * t76, t59 * t76, t10 * t73, t51, -0.2e1 * t1 * t59 + 0.2e1 * t19 * t10, 0.2e1 * t19 * t11 + 0.2e1 * t2 * t59; 0, 0, 0, 0, 0, t56, t59, 0, -t47, -t59 * pkin(6), t55 * t66, (-t48 + t50) * t56, -t67, -t65, 0, -pkin(6) * t66 + (-pkin(2) * t56 + pkin(7) * t59) * t55, pkin(7) * t65 + (t70 - t72) * t56, -t22 * t29 - t23 * t28 + t9 * t33 - t8 * t34, t8 * t22 + t9 * t23 + t39 * t46, t11 * t17, -t17 * t10 - t11 * t16, -t17 * t59, t16 * t59, 0, t27 * t10 + t19 * t16 - t5 * t59, t27 * t11 + t19 * t17 + t6 * t59; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t48, 0.2e1 * t68, 0, 0, 0, 0.2e1 * t72, -0.2e1 * pkin(2) * t55, -0.2e1 * t22 * t34 + 0.2e1 * t23 * t33, t22 ^ 2 + t23 ^ 2 + t46 ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0, 0, 0, t16 * t75, t17 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t69, -t59, t25, -t26, (-t28 * t52 - t29 * t53) * pkin(3), (t52 * t9 + t53 * t8) * pkin(3), 0, 0, t11, -t10, -t59, -t30 * t59 + t1, t31 * t59 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t58, 0, -t55 * pkin(7), -t58 * pkin(7), (t33 * t52 - t34 * t53) * pkin(3), (t22 * t53 + t23 * t52) * pkin(3), 0, 0, t17, -t16, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t52 ^ 2 + t53 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, -t59, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
