% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:41
% EndTime: 2019-12-05 17:16:43
% DurationCPUTime: 0.33s
% Computational Cost: add. (203->58), mult. (485->106), div. (0->0), fcn. (613->10), ass. (0->54)
t38 = cos(qJ(3));
t28 = -t38 * pkin(3) - pkin(2);
t58 = 0.2e1 * t28;
t57 = 0.2e1 * t38;
t56 = pkin(7) + pkin(8);
t34 = sin(qJ(4));
t55 = t34 * pkin(3);
t37 = cos(qJ(5));
t32 = cos(pkin(5));
t35 = sin(qJ(3));
t31 = sin(pkin(5));
t50 = t31 * sin(qJ(2));
t15 = t32 * t38 - t35 * t50;
t16 = t32 * t35 + t38 * t50;
t52 = cos(qJ(4));
t6 = -t52 * t15 + t34 * t16;
t54 = t6 * t37;
t43 = t52 * pkin(3);
t27 = -t43 - pkin(4);
t53 = pkin(4) - t27;
t23 = t56 * t38;
t42 = t52 * t35;
t11 = t34 * t23 + t56 * t42;
t51 = t11 * t37;
t49 = t31 * cos(qJ(2));
t21 = t34 * t38 + t42;
t33 = sin(qJ(5));
t48 = t33 * t21;
t47 = t33 * t37;
t46 = t34 * t35;
t45 = t37 * t21;
t20 = -t52 * t38 + t46;
t44 = -0.2e1 * t21 * t20;
t41 = -pkin(4) * t21 - pkin(9) * t20;
t26 = pkin(9) + t55;
t40 = -t20 * t26 + t21 * t27;
t30 = t37 ^ 2;
t29 = t33 ^ 2;
t24 = 0.2e1 * t47;
t19 = t21 ^ 2;
t18 = t37 * t20;
t17 = t33 * t20;
t14 = t33 * t45;
t12 = t52 * t23 - t56 * t46;
t10 = t11 * t33;
t9 = (-t29 + t30) * t21;
t8 = t20 * pkin(4) - t21 * pkin(9) + t28;
t7 = t34 * t15 + t52 * t16;
t5 = t6 * t33;
t4 = -t33 * t49 + t37 * t7;
t3 = -t33 * t7 - t37 * t49;
t2 = t37 * t12 + t33 * t8;
t1 = -t33 * t12 + t37 * t8;
t13 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t49, -t50, 0, 0, 0, 0, 0, t38 * t49, -t35 * t49, 0, 0, 0, 0, 0, -t20 * t49, -t21 * t49, 0, 0, 0, 0, 0, t3 * t20 + t6 * t48, -t4 * t20 + t6 * t45; 0, 1, 0, 0, t35 ^ 2, t35 * t57, 0, 0, 0, pkin(2) * t57, -0.2e1 * pkin(2) * t35, t19, t44, 0, 0, 0, t20 * t58, t21 * t58, t30 * t19, -0.2e1 * t19 * t47, 0.2e1 * t20 * t45, t33 * t44, t20 ^ 2, 0.2e1 * t1 * t20 + 0.2e1 * t11 * t48, 0.2e1 * t11 * t45 - 0.2e1 * t2 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t54, t5; 0, 0, 0, 0, 0, 0, t35, t38, 0, -t35 * pkin(7), -t38 * pkin(7), 0, 0, t21, -t20, 0, -t11, -t12, t14, t9, t17, t18, 0, t40 * t33 - t51, t40 * t37 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t55, t29, t24, 0, 0, 0, -0.2e1 * t27 * t37, 0.2e1 * t27 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t54, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t11, -t12, t14, t9, t17, t18, 0, t41 * t33 - t51, t41 * t37 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t43, -t55, t29, t24, 0, 0, 0, t53 * t37, -t53 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t29, t24, 0, 0, 0, 0.2e1 * pkin(4) * t37, -0.2e1 * pkin(4) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t48, t20, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t37, 0, -t33 * t26, -t37 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t37, 0, -t33 * pkin(9), -t37 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t13;
