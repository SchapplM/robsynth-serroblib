% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:06
% EndTime: 2019-12-31 20:11:08
% DurationCPUTime: 0.66s
% Computational Cost: add. (316->81), mult. (582->132), div. (0->0), fcn. (534->4), ass. (0->62)
t42 = sin(qJ(2));
t37 = t42 ^ 2;
t44 = cos(qJ(2));
t39 = t44 ^ 2;
t68 = t37 + t39;
t41 = sin(qJ(4));
t26 = t41 * pkin(4) + qJ(3);
t67 = 0.2e1 * t26;
t66 = -0.2e1 * t42;
t65 = 0.2e1 * t44;
t64 = 0.2e1 * qJ(3);
t45 = -pkin(2) - pkin(7);
t63 = t42 * pkin(4);
t43 = cos(qJ(4));
t62 = t43 * pkin(4);
t32 = t42 * pkin(6);
t19 = t42 * pkin(3) + t32;
t61 = t41 * t19;
t60 = t41 * t42;
t59 = t41 * t44;
t58 = t42 * t44;
t57 = t43 * t41;
t56 = t43 * t44;
t55 = t68 * pkin(6) ^ 2;
t34 = t44 * pkin(6);
t20 = t44 * pkin(3) + t34;
t36 = t41 ^ 2;
t38 = t43 ^ 2;
t22 = t36 + t38;
t54 = qJ(5) * t44;
t53 = t44 * qJ(3);
t52 = -0.2e1 * t58;
t51 = -t42 * qJ(3) - pkin(1);
t8 = t45 * t44 + t51;
t4 = t43 * t19 - t41 * t8;
t50 = t41 * t54 + t4;
t5 = t43 * t8 + t61;
t1 = t4 * t43 + t5 * t41;
t49 = -t42 * pkin(2) + t53;
t48 = t42 * t45 + t53;
t46 = qJ(3) ^ 2;
t30 = t43 * t45;
t29 = t43 * t42;
t28 = t38 * t39;
t27 = t36 * t39;
t25 = -0.2e1 * t57;
t24 = 0.2e1 * t58;
t21 = t41 * t56;
t18 = t43 * t52;
t17 = t41 * t52;
t16 = 0.2e1 * t39 * t57;
t15 = -t44 * pkin(2) + t51;
t14 = -t43 * qJ(5) + t30;
t13 = (-qJ(5) + t45) * t41;
t12 = 0.2e1 * t68 * pkin(6);
t11 = t22 * t45;
t9 = (t36 - t38) * t44;
t7 = pkin(4) * t56 + t20;
t6 = t13 * t41 + t14 * t43;
t3 = t61 + (t8 - t54) * t43;
t2 = t50 + t63;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t24, 0, t39, 0, 0, pkin(1) * t65, pkin(1) * t66, t12, pkin(1) ^ 2 + t55, 0, 0, 0, t37, t24, t39, t12, t15 * t65, t15 * t66, t15 ^ 2 + t55, t27, t16, t17, t28, t18, t37, 0.2e1 * t20 * t56 + 0.2e1 * t4 * t42, -0.2e1 * t20 * t59 - 0.2e1 * t5 * t42, (t4 * t41 - t43 * t5) * t65, t20 ^ 2 + t4 ^ 2 + t5 ^ 2, t27, t16, t17, t28, t18, t37, 0.2e1 * t2 * t42 + 0.2e1 * t7 * t56, -0.2e1 * t3 * t42 - 0.2e1 * t7 * t59, (t2 * t41 - t3 * t43) * t65, t2 ^ 2 + t3 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t44, 0, -t32, -t34, 0, 0, 0, -t42, -t44, 0, 0, 0, t49, t32, t34, t49 * pkin(6), -t21, t9, t29, t21, -t60, 0, t20 * t41 + t48 * t43, t20 * t43 - t48 * t41, -t1, t20 * qJ(3) + t1 * t45, -t21, t9, t29, t21, -t60, 0, t14 * t42 + t26 * t56 + t7 * t41, -t13 * t42 - t26 * t59 + t7 * t43, (-t13 * t44 - t2) * t43 + (t14 * t44 - t3) * t41, t3 * t13 + t2 * t14 + t7 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t64, pkin(2) ^ 2 + t46, t38, t25, 0, t36, 0, 0, t41 * t64, t43 * t64, -0.2e1 * t11, t22 * t45 ^ 2 + t46, t38, t25, 0, t36, 0, 0, t41 * t67, t43 * t67, -0.2e1 * t6, t13 ^ 2 + t14 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, t32, 0, 0, 0, 0, 0, 0, t29, -t60, 0, t1, 0, 0, 0, 0, 0, 0, t29, -t60, 0, t2 * t43 + t3 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t22, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, -t56, t42, t4, -t5, 0, 0, 0, 0, -t59, 0, -t56, t42, t50 + 0.2e1 * t63, -t3, pkin(4) * t59, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, -t41, 0, t30, -t41 * t45, 0, 0, 0, 0, t43, 0, -t41, 0, t14, -t13, -t62, t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, 0, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t59, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t43, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t10;
