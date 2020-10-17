% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:44:52
% EndTime: 2019-05-05 14:44:54
% DurationCPUTime: 0.37s
% Computational Cost: add. (441->64), mult. (805->116), div. (0->0), fcn. (911->8), ass. (0->55)
t33 = sin(pkin(10));
t35 = cos(pkin(10));
t38 = sin(qJ(4));
t55 = cos(qJ(4));
t20 = t55 * t33 + t38 * t35;
t60 = -0.2e1 * t20;
t36 = cos(pkin(9));
t27 = -t36 * pkin(1) - pkin(2);
t21 = -t35 * pkin(3) + t27;
t59 = 0.2e1 * t21;
t39 = cos(qJ(5));
t58 = t39 * pkin(5);
t34 = sin(pkin(9));
t25 = t34 * pkin(1) + qJ(3);
t56 = pkin(7) + t25;
t15 = t56 * t33;
t16 = t56 * t35;
t9 = -t38 * t15 + t55 * t16;
t57 = t39 * t9;
t37 = sin(qJ(5));
t31 = t37 ^ 2;
t54 = t31 * t20;
t19 = t38 * t33 - t55 * t35;
t11 = t37 * t19;
t53 = t37 * t20;
t52 = t37 * t39;
t51 = t39 * t20;
t50 = -qJ(6) - pkin(8);
t49 = t33 ^ 2 + t35 ^ 2;
t48 = qJ(6) * t20;
t47 = t19 * t60;
t46 = pkin(5) * t53;
t7 = t19 * pkin(4) - t20 * pkin(8) + t21;
t3 = -t37 * t9 + t39 * t7;
t45 = -pkin(4) * t20 - pkin(8) * t19;
t1 = t19 * pkin(5) - t39 * t48 + t3;
t2 = t57 + (t7 - t48) * t37;
t44 = t1 * t39 + t2 * t37;
t43 = -t1 * t37 + t2 * t39;
t22 = t50 * t37;
t23 = t50 * t39;
t42 = t39 * t22 - t37 * t23;
t41 = -t22 * t37 - t23 * t39;
t8 = t55 * t15 + t38 * t16;
t32 = t39 ^ 2;
t28 = -pkin(4) - t58;
t18 = t20 ^ 2;
t17 = t19 ^ 2;
t14 = t39 * t19;
t13 = t32 * t20;
t12 = t32 * t18;
t10 = -t13 - t54;
t5 = t8 + t46;
t4 = t37 * t7 + t57;
t6 = [1, 0, 0 (t34 ^ 2 + t36 ^ 2) * pkin(1) ^ 2, -0.2e1 * t27 * t35, 0.2e1 * t27 * t33, 0.2e1 * t49 * t25, t49 * t25 ^ 2 + t27 ^ 2, t18, t47, 0, 0, 0, t19 * t59, t20 * t59, t12, -0.2e1 * t18 * t52, 0.2e1 * t19 * t51, t37 * t47, t17, 0.2e1 * t3 * t19 + 0.2e1 * t8 * t53, -0.2e1 * t4 * t19 + 0.2e1 * t8 * t51, t44 * t60, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t19 + t43 * t20; 0, 0, 0, 1, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t18 + t12 + t17; 0, 0, 0, 0, -t35, t33, 0, t27, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, t14, -t11, t10, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, -t8, -t9, t37 * t51, t13 - t54, t11, t14, 0, t45 * t37 - t8 * t39, t8 * t37 + t45 * t39, -t42 * t20 + t43, t1 * t22 - t2 * t23 + t5 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, 0, 0, 0, 0, -t14, t11, -t10, t19 * t28 + t41 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t31, 0.2e1 * t52, 0, 0, 0, 0.2e1 * pkin(4) * t39, -0.2e1 * pkin(4) * t37, 0.2e1 * t41, t22 ^ 2 + t23 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t53, t19, t3, -t4, -pkin(5) * t51, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t51, 0, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, -t37 * pkin(8), -t39 * pkin(8), -t37 * pkin(5), t22 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
