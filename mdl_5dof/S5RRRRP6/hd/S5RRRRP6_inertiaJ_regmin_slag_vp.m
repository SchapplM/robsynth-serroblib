% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:24
% EndTime: 2021-01-16 00:10:26
% DurationCPUTime: 0.50s
% Computational Cost: add. (539->93), mult. (1028->164), div. (0->0), fcn. (1139->6), ass. (0->66)
t60 = cos(qJ(2));
t36 = -t60 * pkin(2) - pkin(1);
t69 = 0.2e1 * t36;
t39 = sin(qJ(4));
t68 = 0.2e1 * t39;
t42 = cos(qJ(4));
t67 = -0.2e1 * t42;
t40 = sin(qJ(3));
t41 = sin(qJ(2));
t59 = cos(qJ(3));
t23 = t40 * t41 - t59 * t60;
t66 = t23 * pkin(4);
t65 = t39 * pkin(4);
t64 = t40 * pkin(2);
t63 = t42 * pkin(4);
t29 = (-pkin(7) - pkin(6)) * t41;
t48 = t60 * pkin(6);
t30 = t60 * pkin(7) + t48;
t12 = -t59 * t29 + t40 * t30;
t24 = t40 * t60 + t59 * t41;
t57 = t39 * t24;
t7 = pkin(4) * t57 + t12;
t62 = t7 * t42;
t47 = t59 * pkin(2);
t34 = -t47 - pkin(3);
t61 = pkin(3) - t34;
t58 = t12 * t42;
t56 = t39 * t42;
t13 = t40 * t29 + t59 * t30;
t55 = t42 * t13;
t18 = t42 * t24;
t54 = -qJ(5) - pkin(8);
t26 = t34 - t63;
t35 = -pkin(3) - t63;
t53 = t26 + t35;
t52 = qJ(5) * t24;
t33 = pkin(8) + t64;
t51 = qJ(5) + t33;
t50 = 0.2e1 * t60;
t49 = -0.2e1 * t24 * t23;
t9 = t23 * pkin(3) - t24 * pkin(8) + t36;
t4 = -t39 * t13 + t42 * t9;
t43 = -t42 * t52 + t4;
t1 = t43 + t66;
t3 = t55 + (t9 - t52) * t39;
t46 = -t1 * t39 + t3 * t42;
t45 = -pkin(3) * t24 - pkin(8) * t23;
t44 = -t23 * t33 + t24 * t34;
t38 = t42 ^ 2;
t37 = t39 ^ 2;
t31 = 0.2e1 * t56;
t28 = t54 * t42;
t27 = t54 * t39;
t22 = t28 * t42;
t21 = t24 ^ 2;
t20 = t51 * t42;
t19 = t51 * t39;
t17 = t42 * t23;
t16 = t39 * t23;
t15 = t20 * t42;
t14 = t39 * t18;
t11 = t12 * t39;
t10 = (-t37 + t38) * t24;
t6 = t7 * t39;
t5 = t39 * t9 + t55;
t2 = [1, 0, 0, t41 ^ 2, t41 * t50, 0, 0, 0, pkin(1) * t50, -0.2e1 * pkin(1) * t41, t21, t49, 0, 0, 0, t23 * t69, t24 * t69, t38 * t21, -0.2e1 * t21 * t56, 0.2e1 * t23 * t18, t39 * t49, t23 ^ 2, 0.2e1 * t12 * t57 + 0.2e1 * t4 * t23, 0.2e1 * t12 * t18 - 0.2e1 * t5 * t23, 0.2e1 * t1 * t23 + 0.2e1 * t7 * t57, 0.2e1 * t7 * t18 - 0.2e1 * t3 * t23, 0.2e1 * (-t1 * t42 - t3 * t39) * t24, t1 ^ 2 + t3 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, t41, t60, 0, -t41 * pkin(6), -t48, 0, 0, t24, -t23, 0, -t12, -t13, t14, t10, t16, t17, 0, t44 * t39 - t58, t44 * t42 + t11, -t19 * t23 + t26 * t57 - t62, t26 * t18 - t20 * t23 + t6, (t19 * t42 - t20 * t39) * t24 + t46, -t1 * t19 + t3 * t20 + t7 * t26; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t47, -0.2e1 * t64, t37, t31, 0, 0, 0, t34 * t67, t34 * t68, t26 * t67, t26 * t68, 0.2e1 * t19 * t39 + 0.2e1 * t15, t19 ^ 2 + t20 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t12, -t13, t14, t10, t16, t17, 0, t45 * t39 - t58, t45 * t42 + t11, t27 * t23 + t35 * t57 - t62, t35 * t18 + t28 * t23 + t6, (-t27 * t42 + t28 * t39) * t24 + t46, t1 * t27 - t3 * t28 + t7 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t47, -t64, t37, t31, 0, 0, 0, t61 * t42, -t61 * t39, -t53 * t42, t53 * t39, t15 - t22 + (t19 - t27) * t39, -t19 * t27 - t20 * t28 + t26 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t37, t31, 0, 0, 0, 0.2e1 * pkin(3) * t42, -0.2e1 * pkin(3) * t39, t35 * t67, t35 * t68, -0.2e1 * t27 * t39 - 0.2e1 * t22, t27 ^ 2 + t28 ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t57, t23, t4, -t5, t43 + 0.2e1 * t66, -t3, -pkin(4) * t18, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42, 0, -t39 * t33, -t42 * t33, -t19, -t20, -t65, -t19 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42, 0, -t39 * pkin(8), -t42 * pkin(8), t27, t28, -t65, t27 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t18, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t39, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t39, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
