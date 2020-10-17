% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:46
% EndTime: 2019-12-31 20:55:48
% DurationCPUTime: 0.70s
% Computational Cost: add. (599->69), mult. (1162->135), div. (0->0), fcn. (1308->6), ass. (0->47)
t39 = sin(qJ(3));
t40 = sin(qJ(2));
t41 = cos(qJ(3));
t42 = cos(qJ(2));
t21 = t39 * t40 - t41 * t42;
t23 = t39 * t42 + t41 * t40;
t37 = sin(pkin(8));
t38 = cos(pkin(8));
t10 = t38 * t21 + t37 * t23;
t62 = t10 ^ 2;
t32 = -t42 * pkin(2) - pkin(1);
t15 = t21 * pkin(3) + t32;
t61 = 0.2e1 * t15;
t60 = 0.2e1 * t32;
t59 = 0.2e1 * t42;
t58 = -pkin(7) - pkin(6);
t57 = t38 * pkin(3);
t56 = t39 * pkin(2);
t12 = -t37 * t21 + t38 * t23;
t55 = t12 * t10;
t34 = t41 * pkin(2);
t31 = t34 + pkin(3);
t19 = t37 * t31 + t38 * t56;
t35 = t40 ^ 2;
t36 = t42 ^ 2;
t54 = t35 + t36;
t50 = t58 * t40;
t51 = t58 * t42;
t13 = t39 * t51 + t41 * t50;
t46 = -t23 * qJ(4) + t13;
t14 = t39 * t50 - t41 * t51;
t8 = -t21 * qJ(4) + t14;
t4 = t37 * t8 - t38 * t46;
t6 = t37 * t46 + t38 * t8;
t53 = t4 ^ 2 + t6 ^ 2;
t33 = t37 * pkin(3);
t52 = -t33 - t19;
t49 = -t38 * t31 + t37 * t56;
t48 = -0.2e1 * t6 * t10 + 0.2e1 * t4 * t12;
t47 = -t49 + t57;
t28 = pkin(4) + t57;
t27 = t33 + qJ(5);
t17 = -pkin(4) + t49;
t16 = qJ(5) + t19;
t9 = t12 ^ 2;
t2 = t10 * pkin(4) - t12 * qJ(5) + t15;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t35, t40 * t59, 0, t36, 0, 0, pkin(1) * t59, -0.2e1 * pkin(1) * t40, 0.2e1 * t54 * pkin(6), t54 * pkin(6) ^ 2 + pkin(1) ^ 2, t23 ^ 2, -0.2e1 * t23 * t21, 0, t21 ^ 2, 0, 0, t21 * t60, t23 * t60, -0.2e1 * t13 * t23 - 0.2e1 * t14 * t21, t13 ^ 2 + t14 ^ 2 + t32 ^ 2, t9, -0.2e1 * t55, 0, t62, 0, 0, t10 * t61, t12 * t61, t48, t15 ^ 2 + t53, t9, 0, 0.2e1 * t55, 0, 0, t62, 0.2e1 * t2 * t10, t48, -0.2e1 * t2 * t12, t2 ^ 2 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t42, 0, -t40 * pkin(6), -t42 * pkin(6), 0, 0, 0, 0, t23, 0, -t21, 0, t13, -t14, (-t21 * t39 - t23 * t41) * pkin(2), (t13 * t41 + t14 * t39) * pkin(2), 0, 0, t12, 0, -t10, 0, -t4, -t6, -t19 * t10 + t12 * t49, t6 * t19 + t4 * t49, 0, t12, 0, 0, t10, 0, -t4, -t16 * t10 + t17 * t12, t6, t6 * t16 + t4 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t56, 0, (t39 ^ 2 + t41 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t49, -0.2e1 * t19, 0, t19 ^ 2 + t49 ^ 2, 0, 0, 0, 1, 0, 0, -0.2e1 * t17, 0, 0.2e1 * t16, t16 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t21, 0, t13, -t14, 0, 0, 0, 0, t12, 0, -t10, 0, -t4, -t6, (-t10 * t37 - t12 * t38) * pkin(3), (t37 * t6 - t38 * t4) * pkin(3), 0, t12, 0, 0, t10, 0, -t4, -t27 * t10 - t28 * t12, t6, t6 * t27 - t4 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t56, 0, 0, 0, 0, 0, 0, 0, 1, t47, t52, 0, (t19 * t37 - t38 * t49) * pkin(3), 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) + t47, 0, 0.2e1 * qJ(5) - t52, t16 * t27 - t17 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t33, 0, (t37 ^ 2 + t38 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t28, 0, 0.2e1 * t27, t27 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, 0, t15, 0, 0, 0, 0, 0, 0, t10, 0, -t12, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
