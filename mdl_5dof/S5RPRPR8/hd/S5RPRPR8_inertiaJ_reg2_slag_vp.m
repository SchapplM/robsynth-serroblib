% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:11
% EndTime: 2019-12-31 18:22:14
% DurationCPUTime: 0.70s
% Computational Cost: add. (474->90), mult. (953->165), div. (0->0), fcn. (989->8), ass. (0->68)
t40 = sin(pkin(9));
t42 = cos(pkin(9));
t44 = sin(qJ(5));
t71 = cos(qJ(5));
t78 = -t44 * t40 + t71 * t42;
t35 = -t42 * pkin(4) - pkin(3);
t77 = 0.2e1 * t35;
t45 = sin(qJ(3));
t76 = 0.2e1 * t45;
t46 = cos(qJ(3));
t75 = -0.2e1 * t46;
t41 = sin(pkin(8));
t74 = t41 * pkin(1);
t43 = cos(pkin(8));
t73 = t43 * pkin(1);
t72 = t46 * pkin(3);
t21 = t71 * t40 + t44 * t42;
t12 = t21 * t45;
t70 = t21 * t12;
t69 = t21 * t46;
t36 = t40 ^ 2;
t68 = t36 * t45;
t33 = pkin(6) + t74;
t38 = t45 ^ 2;
t67 = t38 * t33;
t66 = t40 * t42;
t65 = t40 * t45;
t64 = t40 * t46;
t63 = t42 * t45;
t62 = t42 * t46;
t26 = t45 * t33;
t60 = t45 * t46;
t59 = t46 * t78;
t58 = t46 * t33;
t57 = pkin(7) + qJ(4);
t34 = -pkin(2) - t73;
t18 = -t45 * qJ(4) + t34 - t72;
t7 = t40 * t18 + t42 * t58;
t37 = t42 ^ 2;
t56 = t36 + t37;
t39 = t46 ^ 2;
t55 = t38 + t39;
t54 = 0.2e1 * t60;
t53 = t40 * t63;
t51 = t56 * qJ(4);
t16 = t42 * t18;
t6 = -t40 * t58 + t16;
t50 = -t6 * t40 + t7 * t42;
t49 = -pkin(3) * t45 + qJ(4) * t46;
t31 = t37 * t45;
t30 = t37 * t38;
t29 = t36 * t38;
t28 = t33 ^ 2;
t25 = t38 * t28;
t24 = t57 * t42;
t23 = t57 * t40;
t17 = pkin(4) * t65 + t26;
t14 = t78 * t45;
t11 = t14 ^ 2;
t10 = t12 ^ 2;
t9 = -t44 * t23 + t71 * t24;
t8 = -t71 * t23 - t44 * t24;
t5 = t14 * t78;
t4 = -pkin(7) * t65 + t7;
t3 = -pkin(7) * t63 + t16 + (-t33 * t40 - pkin(4)) * t46;
t2 = t44 * t3 + t71 * t4;
t1 = t71 * t3 - t44 * t4;
t13 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t73, -0.2e1 * t74, 0, (t41 ^ 2 + t43 ^ 2) * pkin(1) ^ 2, t38, t54, 0, t39, 0, 0, t34 * t75, t34 * t76, 0.2e1 * t55 * t33, t39 * t28 + t34 ^ 2 + t25, t30, -0.2e1 * t38 * t66, -0.2e1 * t42 * t60, t29, t40 * t54, t39, 0.2e1 * t40 * t67 - 0.2e1 * t6 * t46, 0.2e1 * t42 * t67 + 0.2e1 * t7 * t46, (-t40 * t7 - t42 * t6) * t76, t6 ^ 2 + t7 ^ 2 + t25, t11, -0.2e1 * t14 * t12, t14 * t75, t10, -t12 * t75, t39, -0.2e1 * t1 * t46 + 0.2e1 * t17 * t12, 0.2e1 * t17 * t14 + 0.2e1 * t2 * t46, -0.2e1 * t1 * t14 - 0.2e1 * t2 * t12, t1 ^ 2 + t17 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t50 - t58) * t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t12 + t2 * t14 - t17 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 + t29 + t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t10 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t46, 0, -t26, -t58, 0, 0, t53, t31 - t68, -t64, -t53, -t62, 0, -t42 * t26 + t49 * t40, t40 * t26 + t49 * t42, t50, -pkin(3) * t26 + t50 * qJ(4), t14 * t21, t5 - t70, -t69, -t12 * t78, -t59, 0, t35 * t12 - t17 * t78 - t8 * t46, t35 * t14 + t17 * t21 + t9 * t46, -t1 * t21 - t9 * t12 - t8 * t14 + t2 * t78, t1 * t8 + t17 * t35 + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t64, t31 + t68, t45 * t51 + t72, 0, 0, 0, 0, 0, 0, t59, -t69, t5 + t70, -t12 * t8 + t14 * t9 - t46 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, 0.2e1 * t66, 0, t37, 0, 0, 0.2e1 * pkin(3) * t42, -0.2e1 * pkin(3) * t40, 0.2e1 * t51, t56 * qJ(4) ^ 2 + pkin(3) ^ 2, t21 ^ 2, 0.2e1 * t21 * t78, 0, t78 ^ 2, 0, 0, -t78 * t77, t21 * t77, -0.2e1 * t8 * t21 + 0.2e1 * t78 * t9, t35 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t63, 0, t26, 0, 0, 0, 0, 0, 0, t12, t14, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t40, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t78, t21, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t12, -t46, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, t78, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t13;
