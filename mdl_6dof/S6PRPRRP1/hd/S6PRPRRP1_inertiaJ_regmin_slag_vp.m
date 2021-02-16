% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRP1
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
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:18
% EndTime: 2021-01-16 01:24:20
% DurationCPUTime: 0.49s
% Computational Cost: add. (376->97), mult. (851->167), div. (0->0), fcn. (985->10), ass. (0->59)
t37 = sin(qJ(4));
t61 = 0.2e1 * t37;
t39 = cos(qJ(5));
t60 = pkin(4) * t39;
t36 = sin(qJ(5));
t59 = t36 * pkin(5);
t32 = sin(pkin(11));
t33 = sin(pkin(6));
t34 = cos(pkin(11));
t38 = sin(qJ(2));
t41 = cos(qJ(2));
t14 = (t32 * t41 + t34 * t38) * t33;
t35 = cos(pkin(6));
t40 = cos(qJ(4));
t10 = t37 * t14 - t35 * t40;
t58 = t10 * t39;
t21 = t32 * pkin(2) + pkin(8);
t57 = t21 * t36;
t28 = t36 ^ 2;
t56 = t28 * t37;
t55 = t33 * t38;
t54 = t33 * t41;
t53 = t36 * t37;
t52 = t36 * t39;
t51 = t36 * t40;
t25 = t39 * t37;
t26 = t39 * t40;
t50 = t40 * t21;
t49 = -qJ(6) - pkin(9);
t48 = qJ(6) * t37;
t47 = t40 * t61;
t46 = t39 * t50;
t22 = -t34 * pkin(2) - pkin(3);
t17 = -t40 * pkin(4) - t37 * pkin(9) + t22;
t15 = t39 * t17;
t45 = -t39 * t48 + t15;
t11 = t40 * t14 + t35 * t37;
t12 = t32 * t55 - t34 * t54;
t3 = -t36 * t11 + t39 * t12;
t4 = t39 * t11 + t36 * t12;
t44 = -t3 * t36 + t4 * t39;
t18 = t49 * t36;
t19 = t49 * t39;
t43 = -t18 * t36 - t19 * t39;
t31 = t40 ^ 2;
t30 = t39 ^ 2;
t29 = t37 ^ 2;
t27 = -t39 * pkin(5) - pkin(4);
t24 = t30 * t37;
t23 = t30 * t29;
t16 = (t21 + t59) * t37;
t9 = t36 * t17 + t46;
t8 = -t36 * t50 + t15;
t7 = t10 * t36;
t6 = t46 + (t17 - t48) * t36;
t5 = (-pkin(5) - t57) * t40 + t45;
t2 = t10 * t53 - t3 * t40;
t1 = t10 * t25 + t4 * t40;
t13 = [1, 0, 0, 0, t12 ^ 2 + t14 ^ 2 + t35 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, t54, -t55, (-t12 * t34 + t14 * t32) * pkin(2), 0, 0, 0, 0, 0, -t12 * t40, t12 * t37, 0, 0, 0, 0, 0, t2, t1, t2, t1, (-t3 * t39 - t36 * t4) * t37, t10 * t16 + t3 * t5 + t4 * t6; 0, 1, 0, 0, (t32 ^ 2 + t34 ^ 2) * pkin(2) ^ 2, t29, t47, 0, 0, 0, -0.2e1 * t22 * t40, t22 * t61, t23, -0.2e1 * t29 * t52, -0.2e1 * t37 * t26, t36 * t47, t31, 0.2e1 * t29 * t57 - 0.2e1 * t8 * t40, 0.2e1 * t29 * t21 * t39 + 0.2e1 * t9 * t40, 0.2e1 * t16 * t53 - 0.2e1 * t5 * t40, 0.2e1 * t16 * t25 + 0.2e1 * t6 * t40, (-t36 * t6 - t39 * t5) * t61, t16 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t40 + t44 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t40 + (-t5 * t36 + t6 * t39) * t37; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t29 + t23 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t58, t7, -t58, t7, t44, t10 * t27 + t3 * t18 - t4 * t19; 0, 0, 0, 0, 0, 0, 0, t37, t40, 0, -t37 * t21, -t50, t36 * t25, t24 - t56, -t51, -t26, 0, -t21 * t25 + (-pkin(4) * t37 + pkin(9) * t40) * t36, pkin(9) * t26 + (t57 - t60) * t37, -t16 * t39 - t18 * t40 + t27 * t53, t16 * t36 - t19 * t40 + t27 * t25, (-t18 * t37 + t6) * t39 + (t19 * t37 - t5) * t36, t16 * t27 + t5 * t18 - t6 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t37, 0, 0, 0, 0, 0, t26, -t51, t26, -t51, t24 + t56, -t40 * t27 + t43 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t28, 0.2e1 * t52, 0, 0, 0, 0.2e1 * t60, -0.2e1 * pkin(4) * t36, -0.2e1 * t27 * t39, 0.2e1 * t27 * t36, 0.2e1 * t43, t18 ^ 2 + t19 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, t3, -t4, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t53, -t40, t8, -t9, (-0.2e1 * pkin(5) - t57) * t40 + t45, -t6, -pkin(5) * t25, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t25, -t53, -t25, 0, -pkin(5) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t39, 0, -t36 * pkin(9), -t39 * pkin(9), t18, t19, -t59, t18 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t25, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t36, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t13;
