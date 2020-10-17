% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:41
% EndTime: 2019-12-31 19:44:42
% DurationCPUTime: 0.43s
% Computational Cost: add. (264->75), mult. (566->144), div. (0->0), fcn. (572->6), ass. (0->53)
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t62 = t38 ^ 2 + t39 ^ 2;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t18 = t42 * t38 - t40 * t39;
t41 = sin(qJ(2));
t13 = t18 * t41;
t61 = 0.2e1 * t13;
t47 = t38 * qJ(4) + pkin(2);
t56 = pkin(3) + pkin(4);
t15 = t56 * t39 + t47;
t60 = 0.2e1 * t15;
t20 = -t39 * pkin(3) - t47;
t59 = -0.2e1 * t20;
t58 = 0.2e1 * t41;
t43 = cos(qJ(2));
t57 = 0.2e1 * t43;
t55 = pkin(2) * t38;
t54 = pkin(6) * t39;
t33 = t41 * pkin(6);
t53 = t43 * pkin(6);
t29 = t39 * t41;
t28 = t41 * t38;
t21 = -t43 * pkin(2) - t41 * qJ(3) - pkin(1);
t11 = t38 * t21 + t39 * t53;
t50 = t62 * qJ(3) ^ 2;
t49 = qJ(3) * t43;
t48 = qJ(4) * t39;
t31 = t38 * qJ(3);
t26 = t38 * t53;
t10 = t39 * t21 - t26;
t8 = -t43 * qJ(4) + t11;
t34 = t43 * pkin(3);
t9 = -t10 + t34;
t46 = t9 * t38 + t8 * t39;
t45 = -t10 * t38 + t11 * t39;
t17 = t40 * t38 + t42 * t39;
t37 = t41 ^ 2;
t25 = t43 * t31;
t23 = (-pkin(7) + qJ(3)) * t39;
t22 = -t38 * pkin(7) + t31;
t19 = 0.2e1 * t62 * qJ(3);
t14 = t17 * t41;
t12 = t33 + (pkin(3) * t38 - t48) * t41;
t7 = -t33 + (-t56 * t38 + t48) * t41;
t6 = t40 * t22 + t42 * t23;
t5 = t42 * t22 - t40 * t23;
t4 = pkin(7) * t28 + t8;
t3 = t43 * pkin(4) + t26 + t34 + (-pkin(7) * t41 - t21) * t39;
t2 = t40 * t3 + t42 * t4;
t1 = t42 * t3 - t40 * t4;
t16 = [1, 0, 0, t37, t41 * t57, 0, 0, 0, pkin(1) * t57, -0.2e1 * pkin(1) * t41, 0.2e1 * t37 * pkin(6) * t38 - 0.2e1 * t10 * t43, 0.2e1 * t11 * t43 + 0.2e1 * t37 * t54, (-t10 * t39 - t11 * t38) * t58, t37 * pkin(6) ^ 2 + t10 ^ 2 + t11 ^ 2, 0.2e1 * t12 * t28 + 0.2e1 * t9 * t43, (-t38 * t8 + t39 * t9) * t58, -0.2e1 * t12 * t29 - 0.2e1 * t8 * t43, t12 ^ 2 + t8 ^ 2 + t9 ^ 2, t14 ^ 2, t14 * t61, t14 * t57, t43 * t61, t43 ^ 2, 0.2e1 * t1 * t43 - 0.2e1 * t7 * t13, 0.2e1 * t7 * t14 - 0.2e1 * t2 * t43; 0, 0, 0, 0, 0, t41, t43, 0, -t33, -t53, t25 + (-t54 - t55) * t41, pkin(6) * t28 + (-pkin(2) * t41 + t49) * t39, t45, -pkin(2) * t33 + t45 * qJ(3), -t12 * t39 + t20 * t28 + t25, t46, -t12 * t38 + (-t20 * t41 - t49) * t39, t46 * qJ(3) + t12 * t20, t14 * t18, t18 * t13 - t14 * t17, t18 * t43, -t17 * t43, 0, -t15 * t13 + t7 * t17 + t5 * t43, t15 * t14 + t7 * t18 - t6 * t43; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t39, -0.2e1 * t55, t19, pkin(2) ^ 2 + t50, t39 * t59, t19, t38 * t59, t20 ^ 2 + t50, t18 ^ 2, -0.2e1 * t18 * t17, 0, 0, 0, t17 * t60, t18 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t29, 0, t33, t28, 0, -t29, t12, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, -pkin(2), -t39, 0, -t38, t20, 0, 0, 0, 0, 0, -t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t29, 0, t9, 0, 0, 0, 0, 0, t42 * t43, -t40 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, t43, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t16;
