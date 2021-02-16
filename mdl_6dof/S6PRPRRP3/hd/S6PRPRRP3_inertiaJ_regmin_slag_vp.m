% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:31
% EndTime: 2021-01-16 01:38:34
% DurationCPUTime: 0.60s
% Computational Cost: add. (588->93), mult. (1264->171), div. (0->0), fcn. (1550->10), ass. (0->63)
t37 = sin(pkin(11));
t39 = cos(pkin(11));
t42 = sin(qJ(4));
t65 = cos(qJ(4));
t24 = t65 * t37 + t42 * t39;
t69 = -0.2e1 * t24;
t31 = -t39 * pkin(3) - pkin(2);
t68 = 0.2e1 * t31;
t23 = t42 * t37 - t65 * t39;
t67 = t23 * pkin(5);
t44 = cos(qJ(5));
t66 = t44 * pkin(5);
t40 = cos(pkin(6));
t38 = sin(pkin(6));
t63 = t38 * sin(qJ(2));
t17 = -t37 * t63 + t39 * t40;
t18 = t40 * t37 + t39 * t63;
t11 = -t65 * t17 + t42 * t18;
t64 = t11 * t44;
t45 = cos(qJ(2));
t62 = t38 * t45;
t41 = sin(qJ(5));
t19 = t41 * t23;
t61 = t41 * t24;
t60 = t41 * t44;
t58 = qJ(3) + pkin(8);
t25 = t58 * t37;
t26 = t58 * t39;
t16 = -t42 * t25 + t65 * t26;
t59 = t44 * t16;
t21 = t44 * t24;
t57 = -qJ(6) - pkin(9);
t56 = t37 ^ 2 + t39 ^ 2;
t35 = t41 ^ 2;
t36 = t44 ^ 2;
t55 = t35 + t36;
t54 = qJ(6) * t24;
t53 = t23 * t69;
t14 = t23 * pkin(4) - t24 * pkin(9) + t31;
t5 = t44 * t14 - t41 * t16;
t52 = -pkin(4) * t24 - pkin(9) * t23;
t47 = -t44 * t54 + t5;
t3 = t47 + t67;
t4 = t59 + (t14 - t54) * t41;
t51 = t3 * t44 + t4 * t41;
t12 = t42 * t17 + t65 * t18;
t7 = -t41 * t12 - t44 * t62;
t8 = t44 * t12 - t41 * t62;
t50 = t8 * t41 + t7 * t44;
t49 = -t17 * t37 + t18 * t39;
t27 = t57 * t41;
t28 = t57 * t44;
t48 = t44 * t27 - t41 * t28;
t15 = t65 * t25 + t42 * t26;
t32 = -pkin(4) - t66;
t22 = t24 ^ 2;
t20 = t44 * t23;
t10 = t11 * t41;
t9 = pkin(5) * t61 + t15;
t6 = t41 * t14 + t59;
t2 = t11 * t21 - t8 * t23;
t1 = t11 * t61 + t7 * t23;
t13 = [1, 0, 0, 0, 0, 0, t38 ^ 2 * t45 ^ 2 + t17 ^ 2 + t18 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, t62, -t63, t39 * t62, t49, pkin(2) * t62 + t49 * qJ(3), 0, 0, 0, 0, 0, -t23 * t62, -t24 * t62, 0, 0, 0, 0, 0, t1, t2, t1, t2, -t50 * t24, t11 * t9 + t7 * t3 + t8 * t4; 0, 1, 0, 0, 0.2e1 * pkin(2) * t39, 0.2e1 * t56 * qJ(3), t56 * qJ(3) ^ 2 + pkin(2) ^ 2, t22, t53, 0, 0, 0, t23 * t68, t24 * t68, t36 * t22, -0.2e1 * t22 * t60, 0.2e1 * t23 * t21, t41 * t53, t23 ^ 2, 0.2e1 * t15 * t61 + 0.2e1 * t5 * t23, 0.2e1 * t15 * t21 - 0.2e1 * t6 * t23, 0.2e1 * t3 * t23 + 0.2e1 * t9 * t61, 0.2e1 * t9 * t21 - 0.2e1 * t4 * t23, t51 * t69, t3 ^ 2 + t4 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, -t39, 0, -pkin(2), 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, t20, -t19, t20, -t19, -t55 * t24, t51; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, -t64, t10, -t64, t10, -t7 * t41 + t8 * t44, t11 * t32 + t7 * t27 - t8 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t15, -t16, t41 * t21, (-t35 + t36) * t24, t19, t20, 0, -t15 * t44 + t52 * t41, t15 * t41 + t52 * t44, t27 * t23 + t32 * t61 - t9 * t44, t32 * t21 + t28 * t23 + t9 * t41, -t48 * t24 - t3 * t41 + t4 * t44, t3 * t27 - t4 * t28 + t9 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t35, 0.2e1 * t60, 0, 0, 0, 0.2e1 * pkin(4) * t44, -0.2e1 * pkin(4) * t41, -0.2e1 * t32 * t44, 0.2e1 * t32 * t41, -0.2e1 * t27 * t41 - 0.2e1 * t28 * t44, t27 ^ 2 + t28 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, t7, -t8, 0, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t61, t23, t5, -t6, t47 + 0.2e1 * t67, -t4, -pkin(5) * t21, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, t44, -t41, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t44, 0, -t41 * pkin(9), -t44 * pkin(9), t27, t28, -t41 * pkin(5), t27 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t21, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t41, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t13;
