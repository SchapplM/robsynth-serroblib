% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:36:19
% EndTime: 2019-05-04 22:36:20
% DurationCPUTime: 0.35s
% Computational Cost: add. (177->61), mult. (420->113), div. (0->0), fcn. (501->10), ass. (0->53)
t26 = sin(pkin(11));
t28 = cos(pkin(11));
t27 = sin(pkin(6));
t35 = cos(qJ(2));
t51 = t27 * t35;
t32 = sin(qJ(2));
t52 = t27 * t32;
t8 = t26 * t52 - t28 * t51;
t57 = t8 ^ 2;
t56 = 2 * qJ(5);
t36 = -pkin(4) - pkin(9);
t34 = cos(qJ(4));
t55 = t34 * pkin(4);
t31 = sin(qJ(4));
t54 = t8 * t31;
t53 = t8 * t34;
t30 = sin(qJ(6));
t50 = t30 * t31;
t49 = t30 * t34;
t48 = t31 * t34;
t33 = cos(qJ(6));
t47 = t33 * t30;
t46 = t33 * t34;
t23 = t31 ^ 2;
t25 = t34 ^ 2;
t45 = t23 + t25;
t44 = qJ(5) * t34;
t43 = t31 * qJ(5);
t42 = -0.2e1 * t48;
t20 = -t28 * pkin(2) - pkin(3);
t10 = (t26 * t35 + t28 * t32) * t27;
t29 = cos(pkin(6));
t5 = t10 * t31 - t29 * t34;
t6 = t10 * t34 + t29 * t31;
t41 = t5 * t31 + t6 * t34;
t40 = -pkin(4) * t31 + t44;
t39 = t31 * t36 + t44;
t38 = t20 - t43;
t24 = t33 ^ 2;
t22 = t30 ^ 2;
t21 = t33 * t31;
t19 = t26 * pkin(2) + pkin(8);
t16 = t34 * t19;
t15 = t31 * t19;
t14 = t34 * pkin(5) + t16;
t13 = t31 * pkin(5) + t15;
t12 = t38 - t55;
t11 = t36 * t34 + t38;
t4 = t33 * t11 + t30 * t13;
t3 = -t30 * t11 + t33 * t13;
t2 = t5 * t30 + t8 * t33;
t1 = -t8 * t30 + t5 * t33;
t7 = [1, 0, 0, 0, t10 ^ 2 + t29 ^ 2 + t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, t51, -t52 (t10 * t26 - t28 * t8) * pkin(2), 0, 0, 0, 0, 0, -t53, t54, t41, t53, -t54, t8 * t12 + t41 * t19, 0, 0, 0, 0, 0, t1 * t31 + t6 * t46, -t2 * t31 - t6 * t49; 0, 1, 0, 0 (t26 ^ 2 + t28 ^ 2) * pkin(2) ^ 2, t23, 0.2e1 * t48, 0, 0, 0, -0.2e1 * t20 * t34, 0.2e1 * t20 * t31, 0.2e1 * t45 * t19, 0.2e1 * t12 * t34, -0.2e1 * t12 * t31, t45 * t19 ^ 2 + t12 ^ 2, t22 * t25, 0.2e1 * t25 * t47, t30 * t42, t33 * t42, t23, 0.2e1 * t14 * t46 + 0.2e1 * t3 * t31, -0.2e1 * t14 * t49 - 0.2e1 * t4 * t31; 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t31 - t5 * t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t5, t6, -t5 * pkin(4) + t6 * qJ(5), 0, 0, 0, 0, 0, t6 * t30, t6 * t33; 0, 0, 0, 0, 0, 0, 0, t31, t34, 0, -t15, -t16, t40, t15, t16, t40 * t19, -t30 * t46 (t22 - t24) * t34, t21, -t50, 0, t14 * t30 + t39 * t33, t14 * t33 - t39 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31, 0, -t34, t31, t43 + t55, 0, 0, 0, 0, 0, t50, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t56, pkin(4) ^ 2 + (qJ(5) ^ 2) t24, -0.2e1 * t47, 0, 0, 0, t30 * t56, t33 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, t15, 0, 0, 0, 0, 0, t21, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t46, t31, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t30, 0, t33 * t36, -t30 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t7;
