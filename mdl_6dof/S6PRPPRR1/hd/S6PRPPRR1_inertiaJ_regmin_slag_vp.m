% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:42:06
% EndTime: 2019-05-04 21:42:08
% DurationCPUTime: 0.34s
% Computational Cost: add. (263->60), mult. (619->119), div. (0->0), fcn. (785->12), ass. (0->50)
t33 = sin(pkin(11));
t36 = cos(pkin(11));
t34 = sin(pkin(6));
t42 = cos(qJ(2));
t51 = t34 * t42;
t40 = sin(qJ(2));
t52 = t34 * t40;
t12 = t33 * t52 - t36 * t51;
t56 = t12 ^ 2;
t27 = -t36 * pkin(2) - pkin(3);
t35 = cos(pkin(12));
t22 = -t35 * pkin(4) + t27;
t55 = 0.2e1 * t22;
t25 = t33 * pkin(2) + qJ(4);
t54 = pkin(8) + t25;
t53 = cos(qJ(5));
t32 = sin(pkin(12));
t39 = sin(qJ(5));
t20 = t39 * t32 - t53 * t35;
t38 = sin(qJ(6));
t15 = t38 * t20;
t21 = t53 * t32 + t39 * t35;
t50 = t38 * t21;
t41 = cos(qJ(6));
t49 = t38 * t41;
t48 = t41 * t21;
t47 = t32 ^ 2 + t35 ^ 2;
t46 = -0.2e1 * t21 * t20;
t45 = -pkin(5) * t21 - pkin(9) * t20;
t14 = (t33 * t42 + t36 * t40) * t34;
t37 = cos(pkin(6));
t10 = -t14 * t32 + t37 * t35;
t11 = t14 * t35 + t37 * t32;
t44 = -t10 * t32 + t11 * t35;
t31 = t41 ^ 2;
t30 = t38 ^ 2;
t19 = t21 ^ 2;
t18 = t54 * t35;
t17 = t54 * t32;
t16 = t41 * t20;
t9 = -t39 * t17 + t53 * t18;
t8 = t53 * t17 + t39 * t18;
t7 = t20 * pkin(5) - t21 * pkin(9) + t22;
t6 = t39 * t10 + t53 * t11;
t5 = -t53 * t10 + t39 * t11;
t4 = t38 * t7 + t41 * t9;
t3 = -t38 * t9 + t41 * t7;
t2 = t12 * t38 + t41 * t6;
t1 = t12 * t41 - t38 * t6;
t13 = [1, 0, 0, 0, t14 ^ 2 + t37 ^ 2 + t56, 0, 0, 0, t10 ^ 2 + t11 ^ 2 + t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t51, -t52 (-t12 * t36 + t14 * t33) * pkin(2), -t12 * t35, t12 * t32, t44, t12 * t27 + t44 * t25, 0, 0, 0, 0, 0, t12 * t20, t12 * t21, 0, 0, 0, 0, 0, t1 * t20 + t5 * t50, -t2 * t20 + t5 * t48; 0, 1, 0, 0 (t33 ^ 2 + t36 ^ 2) * pkin(2) ^ 2, -0.2e1 * t27 * t35, 0.2e1 * t27 * t32, 0.2e1 * t47 * t25, t47 * t25 ^ 2 + t27 ^ 2, t19, t46, 0, 0, 0, t20 * t55, t21 * t55, t31 * t19, -0.2e1 * t19 * t49, 0.2e1 * t20 * t48, t38 * t46, t20 ^ 2, 0.2e1 * t3 * t20 + 0.2e1 * t8 * t50, -0.2e1 * t4 * t20 + 0.2e1 * t8 * t48; 0, 0, 0, 0, t37, 0, 0, 0, t10 * t35 + t11 * t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t35, t32, 0, t27, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, t16, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t41, t5 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t8, -t9, t38 * t48 (-t30 + t31) * t21, t15, t16, 0, t45 * t38 - t8 * t41, t8 * t38 + t45 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, 0, 0, 0, 0, -t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t30, 0.2e1 * t49, 0, 0, 0, 0.2e1 * pkin(5) * t41, -0.2e1 * pkin(5) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t50, t20, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * pkin(9), -t41 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
