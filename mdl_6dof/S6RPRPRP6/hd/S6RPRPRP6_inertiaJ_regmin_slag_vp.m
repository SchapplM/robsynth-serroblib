% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:54:19
% EndTime: 2019-05-05 17:54:20
% DurationCPUTime: 0.44s
% Computational Cost: add. (490->73), mult. (885->124), div. (0->0), fcn. (1018->6), ass. (0->59)
t34 = sin(pkin(9));
t35 = cos(pkin(9));
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t19 = t37 * t34 - t39 * t35;
t20 = t39 * t34 + t37 * t35;
t28 = -t35 * pkin(2) - pkin(1);
t42 = -t20 * qJ(4) + t28;
t10 = t19 * pkin(3) + t42;
t64 = -0.2e1 * t10;
t63 = 0.2e1 * t28;
t62 = 0.2e1 * qJ(4);
t61 = pkin(3) + pkin(8);
t36 = sin(qJ(5));
t60 = t36 * pkin(5);
t51 = pkin(7) + qJ(2);
t23 = t51 * t34;
t24 = t51 * t35;
t12 = t39 * t23 + t37 * t24;
t8 = t20 * pkin(4) + t12;
t59 = t36 * t8;
t38 = cos(qJ(5));
t58 = t38 * pkin(5);
t57 = t20 * t19;
t32 = t36 ^ 2;
t56 = t32 * t19;
t55 = t36 * t19;
t54 = t36 * t20;
t53 = t38 * t19;
t16 = t38 * t20;
t52 = t38 * t36;
t50 = t34 ^ 2 + t35 ^ 2;
t49 = qJ(4) * t19;
t48 = 0.2e1 * t57;
t6 = t61 * t19 + t42;
t47 = qJ(6) * t19 + t6;
t7 = t38 * t8;
t1 = t20 * pkin(5) - t47 * t36 + t7;
t2 = t47 * t38 + t59;
t46 = t1 * t38 + t2 * t36;
t45 = -t1 * t36 + t2 * t38;
t21 = (-qJ(6) - t61) * t36;
t29 = t38 * t61;
t22 = -t38 * qJ(6) - t29;
t44 = t38 * t21 - t36 * t22;
t13 = -t37 * t23 + t39 * t24;
t43 = t20 * t61 + t49;
t33 = t38 ^ 2;
t27 = qJ(4) + t60;
t25 = -t32 - t33;
t18 = t20 ^ 2;
t17 = t19 ^ 2;
t15 = t33 * t19;
t11 = t21 * t36 + t22 * t38;
t9 = -t19 * pkin(4) + t13;
t5 = (-pkin(4) - t58) * t19 + t13;
t4 = t38 * t6 + t59;
t3 = -t36 * t6 + t7;
t14 = [1, 0, 0, 0.2e1 * pkin(1) * t35, -0.2e1 * pkin(1) * t34, 0.2e1 * t50 * qJ(2), t50 * qJ(2) ^ 2 + pkin(1) ^ 2, t18, -0.2e1 * t57, 0, 0, 0, t19 * t63, t20 * t63, 0.2e1 * t12 * t20 - 0.2e1 * t13 * t19, t19 * t64, t20 * t64, t10 ^ 2 + t12 ^ 2 + t13 ^ 2, t32 * t17, 0.2e1 * t17 * t52, t36 * t48, t38 * t48, t18, 0.2e1 * t3 * t20 - 0.2e1 * t9 * t53, -0.2e1 * t4 * t20 + 0.2e1 * t9 * t55, 0.2e1 * t45 * t19, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t35, t34, 0, -pkin(1), 0, 0, 0, 0, 0, t19, t20, 0, -t19, -t20, t10, 0, 0, 0, 0, 0, -t54, -t16, t15 + t56, t45; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, -t12, -t13, -pkin(3) * t20 - t49, t12, t13, -t12 * pkin(3) + t13 * qJ(4), t19 * t52, t15 - t56, t16, -t54, 0, t9 * t36 - t43 * t38, t43 * t36 + t9 * t38, t44 * t19 - t46, t1 * t22 + t2 * t21 + t5 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t62, pkin(3) ^ 2 + qJ(4) ^ 2, t33, -0.2e1 * t52, 0, 0, 0, t36 * t62, t38 * t62, -0.2e1 * t11, t21 ^ 2 + t22 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, t12, 0, 0, 0, 0, 0, t16, -t54, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, t25, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t53, t20, t3, -t4, -pkin(5) * t55, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t38, 0, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0, -t29, t36 * t61, -t58, t22 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t14;
