% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:22
% EndTime: 2019-12-31 19:36:24
% DurationCPUTime: 0.62s
% Computational Cost: add. (462->64), mult. (866->124), div. (0->0), fcn. (965->6), ass. (0->58)
t35 = sin(pkin(8));
t36 = cos(pkin(8));
t38 = sin(qJ(2));
t40 = cos(qJ(2));
t17 = t35 * t38 - t36 * t40;
t19 = t35 * t40 + t36 * t38;
t30 = -t40 * pkin(2) - pkin(1);
t44 = -t19 * qJ(4) + t30;
t7 = t17 * pkin(3) + t44;
t68 = -0.2e1 * t7;
t15 = t17 ^ 2;
t16 = t19 ^ 2;
t63 = t35 * pkin(2);
t26 = qJ(4) + t63;
t67 = t26 ^ 2;
t66 = 0.2e1 * t26;
t65 = 0.2e1 * t30;
t64 = 0.2e1 * t40;
t62 = t36 * pkin(2);
t61 = t19 * t17;
t60 = t26 * t17;
t37 = sin(qJ(5));
t31 = t37 ^ 2;
t59 = t31 * t17;
t58 = t37 * t17;
t57 = t37 * t19;
t39 = cos(qJ(5));
t56 = t39 * t17;
t14 = t39 * t19;
t55 = t39 * t37;
t54 = -qJ(3) - pkin(6);
t33 = t39 ^ 2;
t53 = t31 + t33;
t32 = t38 ^ 2;
t34 = t40 ^ 2;
t52 = t32 + t34;
t51 = -0.2e1 * t61;
t50 = 0.2e1 * t61;
t49 = t17 * t55;
t22 = t54 * t40;
t47 = t54 * t38;
t10 = -t36 * t22 + t35 * t47;
t8 = -t35 * t22 - t36 * t47;
t48 = t10 ^ 2 + t8 ^ 2;
t29 = -pkin(3) - t62;
t4 = (pkin(3) + pkin(7)) * t17 + t44;
t5 = t19 * pkin(4) + t8;
t2 = -t37 * t4 + t39 * t5;
t3 = t37 * t5 + t39 * t4;
t1 = t2 * t39 + t3 * t37;
t46 = -t2 * t37 + t3 * t39;
t25 = -pkin(7) + t29;
t45 = -t19 * t25 + t60;
t43 = -0.2e1 * t10 * t17 + 0.2e1 * t8 * t19;
t13 = t33 * t17;
t12 = t53 * t25;
t6 = -t17 * pkin(4) + t10;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t32, t38 * t64, 0, t34, 0, 0, pkin(1) * t64, -0.2e1 * pkin(1) * t38, 0.2e1 * t52 * pkin(6), t52 * pkin(6) ^ 2 + pkin(1) ^ 2, t16, t51, 0, t15, 0, 0, t17 * t65, t19 * t65, t43, t30 ^ 2 + t48, 0, 0, 0, t16, t51, t15, t43, t17 * t68, t19 * t68, t7 ^ 2 + t48, t31 * t15, 0.2e1 * t15 * t55, t37 * t50, t33 * t15, t39 * t50, t16, 0.2e1 * t2 * t19 - 0.2e1 * t6 * t56, -0.2e1 * t3 * t19 + 0.2e1 * t6 * t58, 0.2e1 * t46 * t17, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t40, 0, -t38 * pkin(6), -t40 * pkin(6), 0, 0, 0, 0, t19, 0, -t17, 0, -t8, -t10, (-t17 * t35 - t19 * t36) * pkin(2), (t10 * t35 - t36 * t8) * pkin(2), 0, -t19, t17, 0, 0, 0, t29 * t19 - t60, t8, t10, t10 * t26 + t8 * t29, t49, t13 - t59, t14, -t49, -t57, 0, t6 * t37 - t39 * t45, t37 * t45 + t6 * t39, -t1, t1 * t25 + t6 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t62, -0.2e1 * t63, 0, (t35 ^ 2 + t36 ^ 2) * pkin(2) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t29, t66, t29 ^ 2 + t67, t33, -0.2e1 * t55, 0, t31, 0, 0, t37 * t66, t39 * t66, -0.2e1 * t12, t53 * t25 ^ 2 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t19, 0, t30, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, t7, 0, 0, 0, 0, 0, 0, -t57, -t14, t13 + t59, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, t8, 0, 0, 0, 0, 0, 0, t14, -t57, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t56, t19, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t37, 0, t39 * t25, -t37 * t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t9;
