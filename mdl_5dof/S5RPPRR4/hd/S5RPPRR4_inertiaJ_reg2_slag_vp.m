% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:22
% EndTime: 2020-01-03 11:31:26
% DurationCPUTime: 0.73s
% Computational Cost: add. (674->74), mult. (1446->158), div. (0->0), fcn. (1587->8), ass. (0->53)
t38 = sin(pkin(9));
t40 = cos(pkin(9));
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t25 = -t43 * t38 + t45 * t40;
t39 = sin(pkin(8));
t64 = -0.2e1 * t39;
t41 = cos(pkin(8));
t63 = -0.2e1 * t41;
t62 = 0.2e1 * t41;
t61 = t41 * pkin(4);
t42 = sin(qJ(5));
t60 = t42 * pkin(4);
t44 = cos(qJ(5));
t59 = t44 * pkin(4);
t26 = t38 * t45 + t40 * t43;
t20 = t26 * t39;
t28 = -pkin(2) * t41 - qJ(3) * t39 - pkin(1);
t24 = t40 * t28;
t56 = t40 * t39;
t12 = -pkin(6) * t56 + t24 + (-qJ(2) * t38 - pkin(3)) * t41;
t51 = qJ(2) * t41;
t18 = t38 * t28 + t40 * t51;
t57 = t38 * t39;
t16 = -pkin(6) * t57 + t18;
t7 = t12 * t43 + t16 * t45;
t5 = -pkin(7) * t20 + t7;
t58 = t44 * t5;
t55 = t40 * t41;
t33 = t39 * qJ(2);
t27 = pkin(3) * t57 + t33;
t34 = t38 ^ 2;
t36 = t40 ^ 2;
t52 = t34 + t36;
t35 = t39 ^ 2;
t50 = t35 * qJ(2);
t49 = t39 * t62;
t22 = t25 * t39;
t6 = t45 * t12 - t16 * t43;
t4 = -pkin(7) * t22 + t6 - t61;
t1 = t44 * t4 - t42 * t5;
t17 = -t38 * t51 + t24;
t48 = t17 * t40 + t18 * t38;
t46 = qJ(2) ^ 2;
t37 = t41 ^ 2;
t32 = t35 * t46;
t15 = t25 * t42 + t26 * t44;
t14 = t25 * t44 - t26 * t42;
t13 = pkin(4) * t20 + t27;
t10 = -t20 * t42 + t22 * t44;
t8 = t44 * t20 + t22 * t42;
t2 = t4 * t42 + t58;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t35, t49, 0, t37, 0, 0, pkin(1) * t62, pkin(1) * t64, 0.2e1 * (t35 + t37) * qJ(2), pkin(1) ^ 2 + t37 * t46 + t32, t36 * t35, -0.2e1 * t40 * t35 * t38, t55 * t64, t34 * t35, t38 * t49, t37, -0.2e1 * t17 * t41 + 0.2e1 * t38 * t50, 0.2e1 * t18 * t41 + 0.2e1 * t40 * t50, t48 * t64, t17 ^ 2 + t18 ^ 2 + t32, t22 ^ 2, -0.2e1 * t22 * t20, t22 * t63, t20 ^ 2, -t20 * t63, t37, 0.2e1 * t20 * t27 - 0.2e1 * t41 * t6, 0.2e1 * t22 * t27 + 0.2e1 * t41 * t7, -0.2e1 * t20 * t7 - 0.2e1 * t22 * t6, t27 ^ 2 + t6 ^ 2 + t7 ^ 2, t10 ^ 2, -0.2e1 * t10 * t8, t10 * t63, t8 ^ 2, t8 * t62, t37, -0.2e1 * t1 * t41 + 0.2e1 * t13 * t8, 0.2e1 * t10 * t13 + 0.2e1 * t2 * t41, -0.2e1 * t1 * t10 - 0.2e1 * t2 * t8, t1 ^ 2 + t13 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t39, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t55, t38 * t41, -t52 * t39, t48, 0, 0, 0, 0, 0, 0, -t25 * t41, t26 * t41, -t20 * t26 - t22 * t25, t25 * t6 + t26 * t7, 0, 0, 0, 0, 0, 0, -t14 * t41, t15 * t41, -t10 * t14 - t15 * t8, t1 * t14 + t15 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 ^ 2 + t26 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t56, 0, t33, 0, 0, 0, 0, 0, 0, t20, t22, 0, t27, 0, 0, 0, 0, 0, 0, t8, t10, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t20, -t41, t6, -t7, 0, 0, 0, 0, t10, 0, -t8, -t41, -t41 * t59 + t1, -t58 + (-t4 + t61) * t42, (-t10 * t44 - t42 * t8) * pkin(4), (t1 * t44 + t2 * t42) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t26, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, (t14 * t44 + t15 * t42) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t60, 0, (t42 ^ 2 + t44 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t8, -t41, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
