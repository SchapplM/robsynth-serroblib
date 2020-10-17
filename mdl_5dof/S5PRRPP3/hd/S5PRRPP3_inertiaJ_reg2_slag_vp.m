% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:32
% EndTime: 2019-12-05 16:13:36
% DurationCPUTime: 0.69s
% Computational Cost: add. (225->81), mult. (567->152), div. (0->0), fcn. (565->6), ass. (0->66)
t40 = sin(pkin(8));
t34 = t40 ^ 2;
t41 = cos(pkin(8));
t35 = t41 ^ 2;
t74 = t34 + t35;
t73 = -0.2e1 * t41;
t42 = sin(qJ(3));
t72 = 0.2e1 * t42;
t71 = pkin(3) * t40;
t36 = t42 ^ 2;
t70 = t36 * pkin(6);
t33 = t42 * pkin(6);
t45 = cos(qJ(2));
t43 = sin(qJ(2));
t44 = cos(qJ(3));
t63 = t44 * t43;
t12 = -t40 * t45 + t41 * t63;
t8 = t12 * t41;
t69 = t36 * t43;
t68 = t40 * t41;
t25 = t40 * t42;
t67 = t40 * t44;
t17 = -pkin(3) * t44 - qJ(4) * t42 - pkin(2);
t66 = t41 * t17;
t27 = t41 * t42;
t65 = t41 * t44;
t29 = t42 * t43;
t64 = t42 * t44;
t5 = pkin(6) * t65 + t40 * t17;
t62 = t74 * qJ(4) ^ 2;
t38 = t44 ^ 2;
t61 = t36 + t38;
t60 = qJ(4) * t44;
t59 = t40 * qJ(4);
t58 = t36 * t68;
t57 = t40 * t29;
t56 = t41 * t29;
t55 = t40 * t64;
t10 = t40 * t63 + t41 * t45;
t54 = t10 * t40 + t8;
t37 = t43 ^ 2;
t28 = t36 * t37;
t53 = t10 ^ 2 + t12 ^ 2 + t28;
t52 = t10 * t44 + t40 * t69;
t51 = qJ(4) * t8 + t10 * t59;
t2 = -qJ(5) * t44 + t5;
t3 = -t66 + (pkin(6) * t40 + pkin(4)) * t44;
t50 = t2 * t41 + t3 * t40;
t4 = -pkin(6) * t67 + t66;
t49 = -t4 * t40 + t41 * t5;
t48 = (t10 * t41 - t12 * t40) * t42;
t47 = pkin(6) ^ 2;
t39 = t45 ^ 2;
t32 = t36 * t47;
t26 = t35 * t36;
t24 = t34 * t36;
t23 = pkin(6) * t69;
t21 = t44 * t59;
t19 = t40 * t27;
t18 = t64 * t73;
t16 = -pkin(4) * t41 - qJ(5) * t40 - pkin(3);
t15 = 0.2e1 * t74 * qJ(4);
t14 = (t34 - t35) * t42;
t7 = t33 + (pkin(4) * t40 - qJ(5) * t41) * t42;
t1 = t12 * t44 + t41 * t69;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 + t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t38 + t28 + t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t43, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t44, -t45 * t42, t61 * t43, pkin(6) * t38 * t43 + pkin(2) * t45 + t23, 0, 0, 0, 0, 0, 0, t52, t1, t48, -t10 * t4 + t12 * t5 + t23, 0, 0, 0, 0, 0, 0, t52, t48, -t1, t10 * t3 + t12 * t2 + t29 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, 0.2e1 * t64, 0, t38, 0, 0, 0.2e1 * pkin(2) * t44, -0.2e1 * pkin(2) * t42, 0.2e1 * t61 * pkin(6), pkin(2) ^ 2 + t38 * t47 + t32, t26, -0.2e1 * t58, t18, t24, 0.2e1 * t55, t38, -0.2e1 * t4 * t44 + 0.2e1 * t40 * t70, 0.2e1 * t41 * t70 + 0.2e1 * t44 * t5, (-t4 * t41 - t40 * t5) * t72, t4 ^ 2 + t5 ^ 2 + t32, t26, t18, 0.2e1 * t58, t38, -0.2e1 * t55, t24, 0.2e1 * t25 * t7 + 0.2e1 * t3 * t44, (-t2 * t40 + t3 * t41) * t72, -0.2e1 * t2 * t44 - 0.2e1 * t27 * t7, t2 ^ 2 + t3 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t63, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t57, t54, -pkin(3) * t29 + t51, 0, 0, 0, 0, 0, 0, -t56, t54, -t57, t16 * t29 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t44, 0, -t33, -t44 * pkin(6), 0, 0, t19, -t14, -t67, -t19, -t65, 0, t21 + (-pkin(6) * t41 - t71) * t42, pkin(6) * t25 + (-pkin(3) * t42 + t60) * t41, t49, -pkin(3) * t33 + qJ(4) * t49, t19, -t67, t14, 0, t65, -t19, t16 * t25 - t41 * t7 + t21, t50, -t7 * t40 + (-t16 * t42 - t60) * t41, qJ(4) * t50 + t7 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t34, 0.2e1 * t68, 0, t35, 0, 0, 0.2e1 * pkin(3) * t41, -0.2e1 * t71, t15, pkin(3) ^ 2 + t62, t34, 0, -0.2e1 * t68, 0, 0, t35, t16 * t73, t15, -0.2e1 * t16 * t40, t16 ^ 2 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t27, 0, t33, 0, 0, 0, 0, 0, 0, t25, 0, -t27, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t40, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t41, 0, -t40, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t27, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
