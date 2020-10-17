% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:28
% EndTime: 2019-12-31 20:07:30
% DurationCPUTime: 0.52s
% Computational Cost: add. (480->85), mult. (991->163), div. (0->0), fcn. (1041->6), ass. (0->47)
t35 = sin(pkin(8));
t36 = cos(pkin(8));
t37 = sin(qJ(4));
t53 = cos(qJ(4));
t61 = -t37 * t35 + t53 * t36;
t38 = sin(qJ(2));
t17 = t61 * t38;
t60 = -0.2e1 * t17;
t22 = t53 * t35 + t37 * t36;
t59 = -0.2e1 * t22;
t30 = -t36 * pkin(3) - pkin(2);
t58 = 0.2e1 * t30;
t39 = cos(qJ(2));
t57 = 0.2e1 * t39;
t56 = pkin(6) * t35;
t31 = t38 * pkin(6);
t55 = t39 * pkin(4);
t54 = t39 * pkin(6);
t24 = -t39 * pkin(2) - t38 * qJ(3) - pkin(1);
t15 = t35 * t24 + t36 * t54;
t50 = t35 * t38;
t11 = -pkin(7) * t50 + t15;
t19 = t36 * t24;
t49 = t36 * t38;
t9 = -pkin(7) * t49 + t19 + (-pkin(3) - t56) * t39;
t4 = t53 * t11 + t37 * t9;
t47 = pkin(7) + qJ(3);
t25 = t47 * t36;
t43 = t47 * t35;
t12 = t37 * t25 + t53 * t43;
t52 = t12 * t39;
t13 = t53 * t25 - t37 * t43;
t51 = t13 * t39;
t23 = pkin(3) * t50 + t31;
t46 = t35 ^ 2 + t36 ^ 2;
t45 = t39 * qJ(5);
t3 = -t37 * t11 + t53 * t9;
t42 = -pkin(2) * t38 + qJ(3) * t39;
t14 = -t35 * t54 + t19;
t41 = -t14 * t35 + t15 * t36;
t34 = t38 ^ 2;
t16 = t22 * t38;
t8 = -pkin(4) * t61 - t22 * qJ(5) + t30;
t5 = t16 * pkin(4) - t17 * qJ(5) + t23;
t2 = -t3 + t55;
t1 = -t45 + t4;
t6 = [1, 0, 0, t34, t38 * t57, 0, 0, 0, pkin(1) * t57, -0.2e1 * pkin(1) * t38, -0.2e1 * t14 * t39 + 0.2e1 * t34 * t56, 0.2e1 * t34 * pkin(6) * t36 + 0.2e1 * t15 * t39, 0.2e1 * (-t14 * t36 - t15 * t35) * t38, t34 * pkin(6) ^ 2 + t14 ^ 2 + t15 ^ 2, t17 ^ 2, t16 * t60, t39 * t60, t16 * t57, t39 ^ 2, 0.2e1 * t23 * t16 - 0.2e1 * t3 * t39, 0.2e1 * t23 * t17 + 0.2e1 * t4 * t39, 0.2e1 * t5 * t16 + 0.2e1 * t2 * t39, -0.2e1 * t1 * t16 + 0.2e1 * t2 * t17, -0.2e1 * t1 * t39 - 0.2e1 * t5 * t17, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t38, t39, 0, -t31, -t54, -pkin(6) * t49 + t42 * t35, pkin(6) * t50 + t42 * t36, t41, -pkin(2) * t31 + t41 * qJ(3), t17 * t22, -t22 * t16 + t17 * t61, -t22 * t39, -t61 * t39, 0, t30 * t16 - t23 * t61 + t52, t30 * t17 + t23 * t22 + t51, t8 * t16 - t5 * t61 + t52, t1 * t61 + t12 * t17 - t13 * t16 + t2 * t22, -t8 * t17 - t5 * t22 - t51, t1 * t13 + t2 * t12 + t5 * t8; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t36, -0.2e1 * pkin(2) * t35, 0.2e1 * t46 * qJ(3), t46 * qJ(3) ^ 2 + pkin(2) ^ 2, t22 ^ 2, -t61 * t59, 0, 0, 0, -t61 * t58, t22 * t58, -0.2e1 * t8 * t61, 0.2e1 * t12 * t22 + 0.2e1 * t13 * t61, t8 * t59, t12 ^ 2 + t13 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t49, 0, t31, 0, 0, 0, 0, 0, t16, t17, t16, 0, -t17, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t35, 0, -pkin(2), 0, 0, 0, 0, 0, -t61, t22, -t61, 0, -t22, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t39, t3, -t4, t3 - 0.2e1 * t55, -pkin(4) * t17 - t16 * qJ(5), -0.2e1 * t45 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t61, 0, -t12, -t13, -t12, -pkin(4) * t22 + qJ(5) * t61, t13, -t12 * pkin(4) + t13 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
