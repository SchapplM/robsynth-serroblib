% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:04:30
% EndTime: 2019-05-05 04:04:31
% DurationCPUTime: 0.42s
% Computational Cost: add. (274->78), mult. (564->143), div. (0->0), fcn. (618->8), ass. (0->60)
t34 = cos(pkin(6));
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t33 = sin(pkin(6));
t60 = t33 * sin(qJ(2));
t13 = t34 * t36 + t39 * t60;
t65 = t13 ^ 2;
t64 = -0.2e1 * t36;
t63 = 0.2e1 * t39;
t62 = 2 * qJ(4);
t41 = -pkin(3) - pkin(9);
t38 = cos(qJ(5));
t61 = t38 * pkin(5);
t40 = cos(qJ(2));
t59 = t33 * t40;
t26 = t36 * pkin(8);
t20 = t36 * pkin(4) + t26;
t35 = sin(qJ(5));
t58 = t35 * t20;
t57 = t35 * t36;
t56 = t35 * t39;
t55 = t36 * t39;
t54 = t38 * t35;
t53 = t38 * t39;
t27 = t39 * pkin(8);
t21 = t39 * pkin(4) + t27;
t30 = t36 ^ 2;
t32 = t39 ^ 2;
t52 = t30 + t32;
t51 = qJ(4) * t39;
t50 = -0.2e1 * t55;
t49 = t36 * t59;
t48 = t39 * t59;
t47 = -t36 * qJ(4) - pkin(2);
t15 = t41 * t39 + t47;
t46 = qJ(6) * t39 - t15;
t45 = -pkin(3) * t36 + t51;
t12 = -t34 * t39 + t36 * t60;
t44 = t12 * t36 + t13 * t39;
t43 = t36 * t41 + t51;
t31 = t38 ^ 2;
t29 = t35 ^ 2;
t25 = t38 * t41;
t24 = t38 * t36;
t23 = t35 * pkin(5) + qJ(4);
t22 = t29 + t31;
t19 = -t39 * pkin(3) + t47;
t18 = -t38 * qJ(6) + t25;
t17 = (-qJ(6) + t41) * t35;
t16 = t38 * t20;
t11 = pkin(5) * t53 + t21;
t8 = -t12 * t35 + t38 * t59;
t7 = t12 * t38 + t35 * t59;
t6 = t17 * t35 + t18 * t38;
t5 = t38 * t15 + t58;
t4 = -t35 * t15 + t16;
t3 = -t46 * t38 + t58;
t2 = t36 * pkin(5) + t46 * t35 + t16;
t1 = -t8 * t35 + t7 * t38;
t9 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 ^ 2 * t40 ^ 2 + t12 ^ 2 + t65, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t65; 0, 0, t59, -t60, 0, 0, 0, 0, 0, t48, -t49, t44, -t48, t49, t44 * pkin(8) - t19 * t59, 0, 0, 0, 0, 0, t13 * t53 + t7 * t36, -t13 * t56 + t8 * t36 (t35 * t7 + t38 * t8) * t39, t13 * t11 + t7 * t2 - t8 * t3; 0, 1, 0, 0, t30, 0.2e1 * t55, 0, 0, 0, pkin(2) * t63, pkin(2) * t64, 0.2e1 * t52 * pkin(8), t19 * t63, t19 * t64, t52 * pkin(8) ^ 2 + t19 ^ 2, t29 * t32, 0.2e1 * t32 * t54, t35 * t50, t38 * t50, t30, 0.2e1 * t21 * t53 + 0.2e1 * t4 * t36, -0.2e1 * t21 * t56 - 0.2e1 * t5 * t36 (t2 * t35 - t3 * t38) * t63, t11 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, t12, t13, -t12 * pkin(3) + t13 * qJ(4), 0, 0, 0, 0, 0, t13 * t35, t13 * t38, -t1, t13 * t23 - t8 * t17 + t7 * t18; 0, 0, 0, 0, 0, 0, t36, t39, 0, -t26, -t27, t45, t26, t27, t45 * pkin(8), -t35 * t53 (t29 - t31) * t39, t24, -t57, 0, t21 * t35 + t43 * t38, t21 * t38 - t43 * t35 (-t17 * t39 - t2) * t38 + (t18 * t39 - t3) * t35, t11 * t23 + t3 * t17 + t2 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t62, pkin(3) ^ 2 + (qJ(4) ^ 2) t31, -0.2e1 * t54, 0, 0, 0, t35 * t62, t38 * t62, -0.2e1 * t6, t17 ^ 2 + t18 ^ 2 + t23 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, t26, 0, 0, 0, 0, 0, t24, -t57, 0, t2 * t38 + t3 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, -t22, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t53, t36, t4, -t5, pkin(5) * t56, t2 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, t25, -t35 * t41, -t61, t18 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t9;
