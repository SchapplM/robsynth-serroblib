% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:10:55
% EndTime: 2019-05-04 22:10:56
% DurationCPUTime: 0.34s
% Computational Cost: add. (346->68), mult. (778->139), div. (0->0), fcn. (971->12), ass. (0->50)
t34 = sin(pkin(11));
t37 = cos(pkin(11));
t35 = sin(pkin(6));
t43 = cos(qJ(2));
t54 = t35 * t43;
t41 = sin(qJ(2));
t55 = t35 * t41;
t14 = t34 * t55 - t37 * t54;
t59 = t14 ^ 2;
t33 = sin(pkin(12));
t36 = cos(pkin(12));
t40 = sin(qJ(4));
t56 = cos(qJ(4));
t22 = t33 * t40 - t36 * t56;
t58 = t22 ^ 2;
t57 = 0.2e1 * t40;
t39 = sin(qJ(6));
t18 = t39 * t22;
t24 = t33 * t56 + t36 * t40;
t53 = t39 * t24;
t42 = cos(qJ(6));
t52 = t39 * t42;
t51 = t42 * t24;
t30 = -t37 * pkin(2) - pkin(3);
t50 = t34 * pkin(2) + pkin(8);
t28 = t33 * pkin(4) + pkin(9);
t29 = -t36 * pkin(4) - pkin(5);
t49 = -t22 * t28 + t24 * t29;
t48 = (-qJ(5) - t50) * t40;
t47 = t56 * t50;
t16 = (t34 * t43 + t37 * t41) * t35;
t38 = cos(pkin(6));
t46 = -t16 * t40 + t38 * t56;
t25 = -t56 * pkin(4) + t30;
t32 = t42 ^ 2;
t31 = t39 ^ 2;
t21 = t24 ^ 2;
t20 = t56 * qJ(5) + t47;
t19 = t42 * t22;
t13 = t16 * t56 + t38 * t40;
t11 = t36 * t20 + t33 * t48;
t9 = t33 * t20 - t36 * t48;
t8 = t22 * pkin(5) - t24 * pkin(9) + t25;
t7 = t36 * t13 + t33 * t46;
t5 = t33 * t13 - t36 * t46;
t4 = t42 * t11 + t39 * t8;
t3 = -t39 * t11 + t42 * t8;
t2 = t14 * t39 + t42 * t7;
t1 = t14 * t42 - t39 * t7;
t6 = [1, 0, 0, 0, t16 ^ 2 + t38 ^ 2 + t59, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t7 ^ 2 + t59, 0, 0, 0, 0, 0, 0, 0; 0, 0, t54, -t55 (-t14 * t37 + t16 * t34) * pkin(2), 0, 0, 0, 0, 0, -t14 * t56, t14 * t40, -t7 * t22 + t5 * t24, t7 * t11 + t14 * t25 + t5 * t9, 0, 0, 0, 0, 0, t1 * t22 + t5 * t53, -t2 * t22 + t5 * t51; 0, 1, 0, 0 (t34 ^ 2 + t37 ^ 2) * pkin(2) ^ 2, t40 ^ 2, t56 * t57, 0, 0, 0, -0.2e1 * t30 * t56, t30 * t57, -0.2e1 * t11 * t22 + 0.2e1 * t9 * t24, t11 ^ 2 + t25 ^ 2 + t9 ^ 2, t32 * t21, -0.2e1 * t21 * t52, 0.2e1 * t22 * t51, -0.2e1 * t22 * t53, t58, 0.2e1 * t3 * t22 + 0.2e1 * t9 * t53, -0.2e1 * t4 * t22 + 0.2e1 * t9 * t51; 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t22 + t7 * t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t24 + t9 * t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t21 + t58, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t13, 0 (t33 * t7 - t36 * t5) * pkin(4), 0, 0, 0, 0, 0, -t5 * t42, t5 * t39; 0, 0, 0, 0, 0, 0, 0, t40, t56, 0, -t40 * t50, -t47 (-t22 * t33 - t24 * t36) * pkin(4) (t11 * t33 - t36 * t9) * pkin(4), t39 * t51 (-t31 + t32) * t24, t18, t19, 0, t49 * t39 - t9 * t42, t9 * t39 + t49 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t40, 0 (-t22 * t36 + t24 * t33) * pkin(4), 0, 0, 0, 0, 0, -t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t33 ^ 2 + t36 ^ 2) * pkin(4) ^ 2, t31, 0.2e1 * t52, 0, 0, 0, -0.2e1 * t29 * t42, 0.2e1 * t29 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, t19, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t53, t22, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42, 0, -t39 * t28, -t42 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
