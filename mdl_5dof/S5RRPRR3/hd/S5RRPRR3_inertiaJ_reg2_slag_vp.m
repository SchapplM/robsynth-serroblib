% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:20
% EndTime: 2022-01-20 10:34:23
% DurationCPUTime: 0.58s
% Computational Cost: add. (376->54), mult. (682->92), div. (0->0), fcn. (596->8), ass. (0->43)
t35 = sin(qJ(5));
t31 = t35 ^ 2;
t38 = cos(qJ(5));
t32 = t38 ^ 2;
t47 = t31 + t32;
t40 = cos(qJ(2));
t29 = t40 * pkin(1);
t25 = t29 + pkin(2);
t33 = sin(pkin(9));
t34 = cos(pkin(9));
t37 = sin(qJ(2));
t52 = t37 * pkin(1);
t17 = t34 * t25 - t33 * t52;
t14 = pkin(3) + t17;
t46 = t34 * t52;
t19 = t33 * t25 + t46;
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t8 = t36 * t14 + t39 * t19;
t6 = pkin(8) + t8;
t55 = t47 * t6;
t26 = t34 * pkin(2);
t24 = t26 + pkin(3);
t53 = t33 * pkin(2);
t20 = t36 * t24 + t39 * t53;
t16 = pkin(8) + t20;
t56 = t47 * t16;
t54 = pkin(4) * t35;
t11 = t39 * t14;
t7 = -t36 * t19 + t11;
t5 = -pkin(4) - t7;
t51 = t5 * t38;
t22 = t39 * t24;
t18 = -t36 * t53 + t22;
t15 = -pkin(4) - t18;
t49 = t15 * t38;
t48 = pkin(8) * t47;
t45 = -t19 - t53;
t30 = pkin(4) * t38;
t23 = 0.2e1 * t35 * t38;
t12 = t15 * t35;
t3 = t5 * t35;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t52, 0, (t37 ^ 2 + t40 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t17, -0.2e1 * t19, 0, t17 ^ 2 + t19 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t7, -0.2e1 * t8, 0, t7 ^ 2 + t8 ^ 2, t31, t23, 0, t32, 0, 0, -0.2e1 * t51, 0.2e1 * t3, 0.2e1 * t55, t47 * t6 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t52, 0, 0, 0, 0, 0, 0, 0, 1, t17 + t26, -t46 + (-pkin(2) - t25) * t33, 0, (t17 * t34 + t19 * t33) * pkin(2), 0, 0, 0, 0, 0, 1, t45 * t36 + t11 + t22, t45 * t39 + (-t14 - t24) * t36, 0, t7 * t18 + t8 * t20, t31, t23, 0, t32, 0, 0, (-t15 - t5) * t38, t12 + t3, t56 + t55, t5 * t15 + t56 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t53, 0, (t33 ^ 2 + t34 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t18, -0.2e1 * t20, 0, t18 ^ 2 + t20 ^ 2, t31, t23, 0, t32, 0, 0, -0.2e1 * t49, 0.2e1 * t12, 0.2e1 * t56, t47 * t16 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t7, -t8, 0, 0, t31, t23, 0, t32, 0, 0, t30 - t51, t3 - t54, t48 + t55, -t5 * pkin(4) + pkin(8) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t18, -t20, 0, 0, t31, t23, 0, t32, 0, 0, t30 - t49, t12 - t54, t48 + t56, -t15 * pkin(4) + pkin(8) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t31, t23, 0, t32, 0, 0, 0.2e1 * t30, -0.2e1 * t54, 0.2e1 * t48, t47 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t35 * t6, -t38 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t35 * t16, -t38 * t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t35 * pkin(8), -t38 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
