% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:35
% EndTime: 2019-12-05 15:51:38
% DurationCPUTime: 0.60s
% Computational Cost: add. (277->77), mult. (691->139), div. (0->0), fcn. (794->10), ass. (0->59)
t27 = sin(pkin(10));
t28 = sin(pkin(5));
t29 = cos(pkin(10));
t33 = sin(qJ(2));
t36 = cos(qJ(2));
t10 = (t27 * t36 + t29 * t33) * t28;
t30 = cos(pkin(5));
t32 = sin(qJ(4));
t35 = cos(qJ(4));
t5 = t10 * t32 - t30 * t35;
t65 = t5 ^ 2;
t54 = t28 * t36;
t55 = t28 * t33;
t8 = t27 * t55 - t29 * t54;
t64 = t8 ^ 2;
t63 = 0.2e1 * t32;
t34 = cos(qJ(5));
t62 = pkin(4) * t34;
t61 = t27 * pkin(2);
t60 = t29 * pkin(2);
t59 = t35 * pkin(4);
t58 = t5 * t35;
t31 = sin(qJ(5));
t23 = t31 ^ 2;
t57 = t23 * t32;
t16 = pkin(7) + t61;
t24 = t32 ^ 2;
t56 = t24 * t16;
t53 = t31 * t32;
t52 = t31 * t34;
t51 = t31 * t35;
t50 = t32 * t16;
t49 = t34 * t32;
t48 = t34 * t35;
t47 = t35 * t16;
t25 = t34 ^ 2;
t46 = t23 + t25;
t26 = t35 ^ 2;
t45 = t24 + t26;
t44 = t35 * t63;
t43 = t31 * t49;
t17 = -pkin(3) - t60;
t42 = t46 * pkin(8);
t7 = t10 * t35 + t30 * t32;
t1 = -t7 * t31 + t8 * t34;
t2 = t8 * t31 + t7 * t34;
t41 = -t1 * t31 + t2 * t34;
t11 = -t32 * pkin(8) + t17 - t59;
t3 = t34 * t11 - t31 * t47;
t4 = t31 * t11 + t34 * t47;
t40 = -t3 * t31 + t4 * t34;
t39 = t5 * t32 + t7 * t35;
t22 = t30 ^ 2;
t20 = t25 * t32;
t19 = t25 * t24;
t18 = t23 * t24;
t14 = t16 ^ 2;
t12 = t24 * t14;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 + (t33 ^ 2 + t36 ^ 2) * t28 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t22 + t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t64 + t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t55, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t10, 0, (t10 * t27 - t29 * t8) * pkin(2), 0, 0, 0, 0, 0, 0, -t8 * t35, t8 * t32, t39, t16 * t39 + t8 * t17, 0, 0, 0, 0, 0, 0, -t1 * t35 + t5 * t53, t2 * t35 + t5 * t49, (-t1 * t34 - t2 * t31) * t32, t1 * t3 + t2 * t4 + t5 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t60, -0.2e1 * t61, 0, (t27 ^ 2 + t29 ^ 2) * pkin(2) ^ 2, t24, t44, 0, t26, 0, 0, -0.2e1 * t17 * t35, t17 * t63, 0.2e1 * t45 * t16, t26 * t14 + t17 ^ 2 + t12, t19, -0.2e1 * t24 * t52, -0.2e1 * t32 * t48, t18, t31 * t44, t26, -0.2e1 * t3 * t35 + 0.2e1 * t31 * t56, 0.2e1 * t34 * t56 + 0.2e1 * t4 * t35, (-t3 * t34 - t31 * t4) * t63, t3 ^ 2 + t4 ^ 2 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t32 - t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t41 - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t40 - t47) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 + t18 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t34, t5 * t31, t41, -t5 * pkin(4) + pkin(8) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t35, 0, -t50, -t47, 0, 0, t43, t20 - t57, -t51, -t43, -t48, 0, -t16 * t49 + (-pkin(4) * t32 + pkin(8) * t35) * t31, pkin(8) * t48 + (t16 * t31 - t62) * t32, t40, -pkin(4) * t50 + pkin(8) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t51, t20 + t57, t32 * t42 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t23, 0.2e1 * t52, 0, t25, 0, 0, 0.2e1 * t62, -0.2e1 * pkin(4) * t31, 0.2e1 * t42, t46 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, -t53, -t35, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t34, 0, -t31 * pkin(8), -t34 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
