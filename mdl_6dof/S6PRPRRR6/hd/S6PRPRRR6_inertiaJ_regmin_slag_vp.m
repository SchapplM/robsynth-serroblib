% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRR6
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:34:32
% EndTime: 2019-05-05 01:34:33
% DurationCPUTime: 0.60s
% Computational Cost: add. (263->92), mult. (622->160), div. (0->0), fcn. (758->10), ass. (0->66)
t38 = sin(qJ(6));
t39 = sin(qJ(5));
t42 = cos(qJ(6));
t43 = cos(qJ(5));
t23 = t38 * t43 + t42 * t39;
t44 = cos(qJ(4));
t15 = t44 * t23;
t71 = -0.2e1 * t15;
t30 = -t43 * pkin(5) - pkin(4);
t70 = 0.2e1 * t30;
t69 = 0.2e1 * t43;
t68 = 2 * qJ(3);
t67 = pkin(9) + pkin(10);
t66 = t38 * pkin(5);
t65 = t42 * pkin(5);
t40 = sin(qJ(4));
t24 = t40 * pkin(4) - t44 * pkin(9) + qJ(3);
t46 = -pkin(2) - pkin(8);
t54 = t43 * t46;
t48 = t40 * t54;
t7 = t48 + (-pkin(10) * t44 + t24) * t39;
t64 = t42 * t7;
t63 = t23 * t40;
t36 = sin(pkin(6));
t41 = sin(qJ(2));
t62 = t36 * t41;
t45 = cos(qJ(2));
t61 = t36 * t45;
t60 = t39 * t40;
t59 = t39 * t43;
t58 = t39 * t44;
t57 = t39 * t46;
t56 = t40 * t46;
t55 = t43 * t40;
t29 = t43 * t44;
t22 = t38 * t39 - t42 * t43;
t53 = t44 * t22;
t52 = t44 * t40;
t51 = t44 * t46;
t33 = t40 ^ 2;
t35 = t44 ^ 2;
t50 = -t33 - t35;
t49 = -0.2e1 * t52;
t20 = t43 * t24;
t6 = -pkin(10) * t29 + t20 + (pkin(5) - t57) * t40;
t1 = -t38 * t7 + t42 * t6;
t47 = -pkin(4) * t44 - pkin(9) * t40;
t37 = cos(pkin(6));
t34 = t43 ^ 2;
t32 = t39 ^ 2;
t26 = t67 * t43;
t25 = t67 * t39;
t21 = (pkin(5) * t39 - t46) * t44;
t19 = t37 * t44 - t40 * t61;
t18 = t37 * t40 + t44 * t61;
t16 = -t38 * t60 + t42 * t55;
t13 = t39 * t24 + t48;
t12 = -t39 * t56 + t20;
t11 = -t38 * t25 + t42 * t26;
t10 = -t42 * t25 - t38 * t26;
t9 = t19 * t43 + t39 * t62;
t8 = -t19 * t39 + t43 * t62;
t4 = t38 * t8 + t42 * t9;
t3 = -t38 * t9 + t42 * t8;
t2 = t38 * t6 + t64;
t5 = [1, 0, 0, 0, 0, 0, t37 ^ 2 + (t41 ^ 2 + t45 ^ 2) * t36 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t61, -t62, -t61, t62 (pkin(2) * t45 + qJ(3) * t41) * t36, 0, 0, 0, 0, 0, t40 * t62, t44 * t62, 0, 0, 0, 0, 0, t18 * t58 + t8 * t40, t18 * t29 - t9 * t40, 0, 0, 0, 0, 0, t18 * t15 + t3 * t40, -t18 * t53 - t4 * t40; 0, 1, 0, 0, -0.2e1 * pkin(2), t68, pkin(2) ^ 2 + (qJ(3) ^ 2) t35, t49, 0, 0, 0, t40 * t68, t44 * t68, t34 * t35, -0.2e1 * t35 * t59, t52 * t69, t39 * t49, t33, 0.2e1 * t12 * t40 - 0.2e1 * t35 * t57, -0.2e1 * t13 * t40 - 0.2e1 * t35 * t54, t53 ^ 2, -t53 * t71, -0.2e1 * t53 * t40, t40 * t71, t33, 0.2e1 * t1 * t40 + 0.2e1 * t21 * t15, -0.2e1 * t2 * t40 - 0.2e1 * t21 * t53; 0, 0, 0, 0, 0, 0, -t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t39, t50 * t43, 0, 0, 0, 0, 0, -t44 * t15 - t40 * t63, -t16 * t40 + t44 * t53; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, 0, 0, 0, 0, -t18 * t43, t18 * t39, 0, 0, 0, 0, 0, t18 * t22, t18 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t40, 0, t51, -t56, t39 * t29 (-t32 + t34) * t44, t60, t55, 0, t47 * t39 + t43 * t51, -t39 * t51 + t47 * t43, -t53 * t23, -t23 * t15 + t22 * t53, t63, -t22 * t40, 0, t10 * t40 + t30 * t15 + t21 * t22, -t11 * t40 + t21 * t23 - t30 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t40, 0, 0, 0, 0, 0, t29, -t58, 0, 0, 0, 0, 0, -t53, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, 0.2e1 * t59, 0, 0, 0, pkin(4) * t69, -0.2e1 * pkin(4) * t39, t23 ^ 2, -0.2e1 * t23 * t22, 0, 0, 0, t22 * t70, t23 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t58, t40, t12, -t13, 0, 0, -t53, -t15, t40, t40 * t65 + t1, -t64 + (-t40 * pkin(5) - t6) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t55, 0, 0, 0, 0, 0, -t63, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t43, 0, -t39 * pkin(9), -t43 * pkin(9), 0, 0, t23, -t22, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t15, t40, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
