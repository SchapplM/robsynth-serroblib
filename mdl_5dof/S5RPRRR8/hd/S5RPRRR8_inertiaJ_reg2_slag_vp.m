% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:06
% EndTime: 2019-12-31 19:06:08
% DurationCPUTime: 0.68s
% Computational Cost: add. (441->91), mult. (718->139), div. (0->0), fcn. (750->6), ass. (0->59)
t32 = sin(qJ(5));
t33 = sin(qJ(4));
t35 = cos(qJ(5));
t36 = cos(qJ(4));
t10 = t32 * t33 - t35 * t36;
t12 = t32 * t36 + t35 * t33;
t68 = (t10 * t32 + t12 * t35) * pkin(4);
t28 = t33 ^ 2;
t30 = t36 ^ 2;
t45 = t28 + t30;
t67 = t45 * pkin(7);
t34 = sin(qJ(3));
t37 = cos(qJ(3));
t38 = -pkin(1) - pkin(2);
t20 = t37 * qJ(2) + t34 * t38;
t17 = -pkin(7) + t20;
t43 = t45 * t17;
t18 = t34 * qJ(2) - t37 * t38;
t16 = pkin(3) + t18;
t57 = t36 * pkin(4);
t9 = t16 + t57;
t66 = -0.2e1 * t9;
t65 = t10 ^ 2;
t64 = t12 ^ 2;
t26 = -pkin(3) - t57;
t63 = 0.2e1 * t26;
t62 = -0.2e1 * t33;
t61 = 0.2e1 * t36;
t60 = -pkin(8) - pkin(7);
t59 = t32 * pkin(4);
t58 = t35 * pkin(4);
t56 = pkin(3) + t16;
t55 = pkin(8) - t17;
t54 = -t26 + t9;
t53 = t12 * t10;
t52 = t33 * t34;
t51 = t33 * t36;
t50 = t36 * t34;
t49 = t37 * t10;
t48 = t37 * t12;
t47 = t37 * t33;
t46 = t37 * t36;
t44 = -0.2e1 * t53;
t14 = t45 * t34;
t7 = t12 * t34;
t8 = -t32 * t52 + t35 * t50;
t42 = t8 * t10 - t7 * t12;
t31 = t37 ^ 2;
t29 = t34 ^ 2;
t24 = 0.2e1 * t51;
t22 = t60 * t36;
t21 = t60 * t33;
t6 = t55 * t36;
t5 = t55 * t33;
t4 = t32 * t21 - t35 * t22;
t3 = t35 * t21 + t32 * t22;
t2 = t32 * t5 - t35 * t6;
t1 = t32 * t6 + t35 * t5;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 * pkin(1), 0, 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t18, 0.2e1 * t20, 0, t18 ^ 2 + t20 ^ 2, t28, t24, 0, t30, 0, 0, t16 * t61, t16 * t62, -0.2e1 * t43, t45 * t17 ^ 2 + t16 ^ 2, t64, t44, 0, t65, 0, 0, t10 * t66, t12 * t66, 0.2e1 * t1 * t12 + 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t37, t34, 0, -t18 * t37 + t20 * t34, 0, 0, 0, 0, 0, 0, -t46, t47, -t14, t17 * t14 - t16 * t37, 0, 0, 0, 0, 0, 0, t49, t48, t42, -t1 * t7 + t2 * t8 - t9 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 + t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t29 + t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t18, -t20, 0, 0, -t28, -0.2e1 * t51, 0, -t30, 0, 0, -t56 * t36, t56 * t33, t43 - t67, -t16 * pkin(3) + pkin(7) * t43, -t64, 0.2e1 * t53, 0, -t65, 0, 0, t54 * t10, t54 * t12, (-t1 + t3) * t12 + (-t2 + t4) * t10, t1 * t3 + t2 * t4 + t9 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t34, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t47, t14, t37 * pkin(3) + pkin(7) * t14, 0, 0, 0, 0, 0, 0, -t49, -t48, -t42, -t37 * t26 - t7 * t3 + t8 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t28, t24, 0, t30, 0, 0, pkin(3) * t61, pkin(3) * t62, 0.2e1 * t67, t45 * pkin(7) ^ 2 + pkin(3) ^ 2, t64, t44, 0, t65, 0, 0, t10 * t63, t12 * t63, -0.2e1 * t4 * t10 - 0.2e1 * t3 * t12, t26 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, -t36, 0, -t33 * t17, -t36 * t17, 0, 0, 0, 0, -t12, 0, t10, 0, t1, -t2, t68, (t1 * t35 + t2 * t32) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t50, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, (t32 * t8 - t35 * t7) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t36, 0, -t33 * pkin(7), -t36 * pkin(7), 0, 0, 0, 0, t12, 0, -t10, 0, t3, -t4, -t68, (t3 * t35 + t32 * t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t59, 0, (t32 ^ 2 + t35 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, t10, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t11;
