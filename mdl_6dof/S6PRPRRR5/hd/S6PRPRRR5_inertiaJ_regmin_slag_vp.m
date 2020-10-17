% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:18:29
% EndTime: 2019-05-05 01:18:31
% DurationCPUTime: 0.50s
% Computational Cost: add. (261->79), mult. (551->120), div. (0->0), fcn. (708->10), ass. (0->59)
t40 = sin(qJ(5));
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t57 = cos(qJ(5));
t24 = t40 * t44 + t41 * t57;
t21 = t24 ^ 2;
t23 = t40 * t41 - t44 * t57;
t22 = t23 ^ 2;
t64 = -t21 - t22;
t63 = -0.2e1 * t23;
t62 = 0.2e1 * t24;
t61 = 2 * qJ(3);
t60 = t40 * pkin(4);
t43 = cos(qJ(6));
t38 = cos(pkin(6));
t37 = sin(pkin(6));
t45 = cos(qJ(2));
t53 = t37 * t45;
t15 = -t38 * t41 - t44 * t53;
t16 = t38 * t44 - t41 * t53;
t6 = -t15 * t57 + t16 * t40;
t59 = t6 * t43;
t50 = t57 * pkin(4);
t31 = -t50 - pkin(5);
t58 = pkin(5) - t31;
t46 = -pkin(2) - pkin(8);
t26 = (-pkin(9) + t46) * t41;
t32 = t44 * t46;
t49 = -pkin(9) * t44 + t32;
t11 = t26 * t40 - t49 * t57;
t56 = t11 * t43;
t39 = sin(qJ(6));
t19 = t23 * t39;
t55 = t23 * t43;
t42 = sin(qJ(2));
t54 = t37 * t42;
t52 = t39 * t43;
t18 = t43 * t24;
t28 = pkin(4) * t41 + qJ(3);
t51 = t23 * t62;
t48 = pkin(5) * t23 - pkin(10) * t24;
t30 = pkin(10) + t60;
t47 = -t23 * t31 - t24 * t30;
t36 = t43 ^ 2;
t35 = t39 ^ 2;
t27 = 0.2e1 * t52;
t17 = t39 * t24;
t14 = t23 * t52;
t12 = t26 * t57 + t40 * t49;
t10 = t11 * t39;
t9 = (t35 - t36) * t23;
t8 = pkin(5) * t24 + pkin(10) * t23 + t28;
t7 = t15 * t40 + t16 * t57;
t5 = t6 * t39;
t4 = t39 * t54 + t43 * t7;
t3 = -t39 * t7 + t43 * t54;
t2 = t12 * t43 + t39 * t8;
t1 = -t12 * t39 + t43 * t8;
t13 = [1, 0, 0, 0, 0, 0, t38 ^ 2 + (t42 ^ 2 + t45 ^ 2) * t37 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t53, -t54, -t53, t54 (pkin(2) * t45 + qJ(3) * t42) * t37, 0, 0, 0, 0, 0, t41 * t54, t44 * t54, 0, 0, 0, 0, 0, t24 * t54, -t23 * t54, 0, 0, 0, 0, 0, -t19 * t6 + t24 * t3, -t24 * t4 - t55 * t6; 0, 1, 0, 0, -0.2e1 * pkin(2), t61, pkin(2) ^ 2 + (qJ(3) ^ 2) t44 ^ 2, -0.2e1 * t44 * t41, 0, 0, 0, t41 * t61, t44 * t61, t22, t51, 0, 0, 0, t28 * t62, t28 * t63, t36 * t22, -0.2e1 * t22 * t52, t18 * t63, t39 * t51, t21, 0.2e1 * t1 * t24 - 0.2e1 * t11 * t19, -0.2e1 * t11 * t55 - 0.2e1 * t2 * t24; 0, 0, 0, 0, 0, 0, -t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t39, t64 * t43; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t59, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, t32, -t41 * t46, 0, 0, -t23, -t24, 0, -t11, -t12, -t14, t9, t17, t18, 0, t39 * t47 - t56, t43 * t47 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, 0, 0, 0, 0, -t23, -t24, 0, 0, 0, 0, 0, -t55, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t60, t35, t27, 0, 0, 0, -0.2e1 * t31 * t43, 0.2e1 * t31 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t59, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, -t11, -t12, -t14, t9, t17, t18, 0, t39 * t48 - t56, t43 * t48 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, 0, 0, 0, 0, -t55, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t60, t35, t27, 0, 0, 0, t58 * t43, -t58 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t35, t27, 0, 0, 0, 0.2e1 * pkin(5) * t43, -0.2e1 * pkin(5) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t19, t24, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t43, 0, -t39 * t30, -t43 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t43, 0, -t39 * pkin(10), -t43 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
