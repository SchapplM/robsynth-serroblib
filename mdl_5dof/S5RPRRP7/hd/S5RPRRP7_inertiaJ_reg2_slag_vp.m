% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:37
% EndTime: 2019-12-31 18:45:39
% DurationCPUTime: 0.66s
% Computational Cost: add. (300->76), mult. (604->136), div. (0->0), fcn. (548->6), ass. (0->64)
t40 = cos(qJ(4));
t34 = t40 ^ 2;
t39 = sin(qJ(3));
t27 = t34 * t39;
t38 = sin(qJ(4));
t32 = t38 ^ 2;
t60 = t32 * t39;
t11 = t27 + t60;
t71 = t32 + t34;
t45 = -pkin(4) * t40 - qJ(5) * t38;
t14 = -pkin(3) + t45;
t70 = -0.2e1 * t14;
t69 = 0.2e1 * t39;
t68 = pkin(3) * t38;
t67 = pkin(3) * t40;
t36 = sin(pkin(8));
t66 = t36 * pkin(1);
t37 = cos(pkin(8));
t65 = t37 * pkin(1);
t64 = t38 * pkin(7);
t63 = t40 * pkin(7);
t41 = cos(qJ(3));
t62 = t41 * pkin(3);
t22 = pkin(6) + t66;
t54 = t41 * t22;
t23 = -pkin(2) - t65;
t9 = -pkin(7) * t39 + t23 - t62;
t4 = t38 * t9 + t40 * t54;
t61 = t22 * t38;
t33 = t39 ^ 2;
t59 = t33 * t22;
t25 = t38 * t39;
t58 = t38 * t40;
t57 = t38 * t41;
t56 = t39 * t22;
t55 = t39 * t41;
t28 = t40 * t39;
t29 = t40 * t41;
t53 = t11 * pkin(7);
t52 = t71 * pkin(7) ^ 2;
t35 = t41 ^ 2;
t51 = t33 + t35;
t50 = t41 * qJ(5);
t49 = t38 * t55;
t48 = t33 * t58;
t1 = -t50 + t4;
t7 = t40 * t9;
t2 = -t7 + (pkin(4) + t61) * t41;
t47 = t1 * t40 + t2 * t38;
t3 = -t38 * t54 + t7;
t46 = -t3 * t38 + t4 * t40;
t44 = -pkin(4) * t38 + qJ(5) * t40;
t26 = t34 * t33;
t24 = t32 * t33;
t21 = t22 ^ 2;
t19 = pkin(7) * t57;
t17 = t38 * t28;
t16 = -0.2e1 * t39 * t29;
t15 = t33 * t21;
t13 = 0.2e1 * t71 * pkin(7);
t10 = -t27 + t60;
t8 = t26 + t24 + t35;
t5 = (t22 - t44) * t39;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t66, 0, (t36 ^ 2 + t37 ^ 2) * pkin(1) ^ 2, t33, 0.2e1 * t55, 0, t35, 0, 0, -0.2e1 * t23 * t41, t23 * t69, 0.2e1 * t51 * t22, t21 * t35 + t23 ^ 2 + t15, t26, -0.2e1 * t48, t16, t24, 0.2e1 * t49, t35, -0.2e1 * t3 * t41 + 0.2e1 * t38 * t59, 0.2e1 * t4 * t41 + 0.2e1 * t40 * t59, (-t3 * t40 - t38 * t4) * t69, t3 ^ 2 + t4 ^ 2 + t15, t26, t16, 0.2e1 * t48, t35, -0.2e1 * t49, t24, 0.2e1 * t2 * t41 + 0.2e1 * t25 * t5, (-t1 * t38 + t2 * t40) * t69, -0.2e1 * t1 * t41 - 0.2e1 * t28 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t46 - t54) * t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t47 - t5 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t41, 0, -t56, -t54, 0, 0, t17, -t10, -t57, -t17, -t29, 0, t19 + (-t22 * t40 - t68) * t39, pkin(7) * t29 + (t61 - t67) * t39, t46, -pkin(3) * t56 + pkin(7) * t46, t17, -t57, t10, 0, t29, -t17, t14 * t25 - t40 * t5 + t19, t47, -t5 * t38 + (-pkin(7) * t41 - t14 * t39) * t40, pkin(7) * t47 + t5 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t57, t11, t53 + t62, 0, 0, 0, 0, 0, 0, t29, t11, t57, -t14 * t41 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t32, 0.2e1 * t58, 0, t34, 0, 0, 0.2e1 * t67, -0.2e1 * t68, t13, pkin(3) ^ 2 + t52, t32, 0, -0.2e1 * t58, 0, 0, t34, t40 * t70, t13, t38 * t70, t14 ^ 2 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t25, -t41, t3, -t4, 0, 0, 0, t28, 0, -t41, t25, 0, t7 + (-0.2e1 * pkin(4) - t61) * t41, t45 * t39, -0.2e1 * t50 + t4, -pkin(4) * t2 + qJ(5) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t28, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, t28, t44 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t40, 0, -t64, -t63, 0, 0, 0, t38, 0, 0, -t40, 0, -t64, t44, t63, t44 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t28, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
