% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR15_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:24
% EndTime: 2019-12-31 18:37:27
% DurationCPUTime: 0.73s
% Computational Cost: add. (437->90), mult. (866->173), div. (0->0), fcn. (905->6), ass. (0->57)
t41 = cos(qJ(3));
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t38 = sin(qJ(5));
t40 = cos(qJ(5));
t69 = -t38 * t36 + t40 * t37;
t13 = t69 * t41;
t29 = -t37 * pkin(4) - pkin(3);
t68 = 0.2e1 * t29;
t67 = 0.2e1 * t37;
t39 = sin(qJ(3));
t66 = 0.2e1 * t39;
t65 = 2 * qJ(2);
t64 = t41 * pkin(3);
t63 = t69 * t39;
t19 = t40 * t36 + t38 * t37;
t62 = t19 * t39;
t34 = t41 ^ 2;
t42 = -pkin(1) - pkin(6);
t61 = t34 * t42;
t60 = t36 * t37;
t59 = t36 * t41;
t58 = t36 * t42;
t27 = t37 * t41;
t56 = t39 * t42;
t54 = t41 * t19;
t53 = t41 * t39;
t52 = t41 * t42;
t51 = pkin(7) + qJ(4);
t21 = t39 * pkin(3) - t41 * qJ(4) + qJ(2);
t8 = t36 * t21 + t37 * t56;
t31 = t36 ^ 2;
t32 = t37 ^ 2;
t50 = t31 + t32;
t33 = t39 ^ 2;
t26 = t33 + t34;
t49 = -0.2e1 * t53;
t48 = t36 * t27;
t47 = t50 * qJ(4);
t15 = t37 * t21;
t7 = -t36 * t56 + t15;
t46 = -t7 * t36 + t8 * t37;
t45 = -qJ(4) * t39 - t64;
t44 = qJ(2) ^ 2;
t35 = t42 ^ 2;
t30 = t34 * t35;
t23 = t51 * t37;
t22 = t51 * t36;
t20 = t26 * t42;
t16 = (pkin(4) * t36 - t42) * t41;
t6 = -t38 * t22 + t40 * t23;
t5 = -t40 * t22 - t38 * t23;
t4 = -pkin(7) * t59 + t8;
t3 = -pkin(7) * t27 + t15 + (pkin(4) - t58) * t39;
t2 = t38 * t3 + t40 * t4;
t1 = t40 * t3 - t38 * t4;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t65, pkin(1) ^ 2 + t44, t34, t49, 0, t33, 0, 0, t39 * t65, t41 * t65, -0.2e1 * t20, t33 * t35 + t30 + t44, t32 * t34, -0.2e1 * t34 * t60, t53 * t67, t31 * t34, t36 * t49, t33, -0.2e1 * t34 * t58 + 0.2e1 * t7 * t39, -0.2e1 * t37 * t61 - 0.2e1 * t8 * t39, 0.2e1 * (-t36 * t8 - t37 * t7) * t41, t7 ^ 2 + t8 ^ 2 + t30, t13 ^ 2, -0.2e1 * t13 * t54, t13 * t66, t54 ^ 2, -t54 * t66, t33, 0.2e1 * t1 * t39 + 0.2e1 * t16 * t54, 0.2e1 * t16 * t13 - 0.2e1 * t2 * t39, -0.2e1 * t1 * t13 - 0.2e1 * t2 * t54, t1 ^ 2 + t16 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t26, t20, 0, 0, 0, 0, 0, 0, -t26 * t36, -t26 * t37, 0, t46 * t39 + t61, 0, 0, 0, 0, 0, 0, -t39 * t62 - t41 * t54, -t41 * t13 - t39 * t63, t13 * t62 - t54 * t63, -t1 * t62 - t16 * t41 + t2 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t33 + t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 ^ 2 + t63 ^ 2 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t39, 0, t52, -t56, 0, 0, t48, (-t31 + t32) * t41, t36 * t39, -t48, t37 * t39, 0, t45 * t36 + t37 * t52, -t36 * t52 + t45 * t37, t46, pkin(3) * t52 + t46 * qJ(4), t13 * t19, t13 * t69 - t19 * t54, t62, -t54 * t69, t63, 0, -t16 * t69 + t29 * t54 + t5 * t39, t29 * t13 + t16 * t19 - t6 * t39, -t1 * t19 - t5 * t13 + t2 * t69 - t54 * t6, t1 * t5 + t16 * t29 + t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t59, t50 * t39, t39 * t47 + t64, 0, 0, 0, 0, 0, 0, t13, -t54, t19 * t62 + t63 * t69, -t41 * t29 - t5 * t62 + t6 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t31, 0.2e1 * t60, 0, t32, 0, 0, pkin(3) * t67, -0.2e1 * pkin(3) * t36, 0.2e1 * t47, t50 * qJ(4) ^ 2 + pkin(3) ^ 2, t19 ^ 2, 0.2e1 * t19 * t69, 0, t69 ^ 2, 0, 0, -t69 * t68, t19 * t68, -0.2e1 * t5 * t19 + 0.2e1 * t6 * t69, t29 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t27, 0, -t52, 0, 0, 0, 0, 0, 0, t54, t13, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t69, t19, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t54, t39, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t63, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t69, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t9;
