% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:01
% DurationCPUTime: 0.61s
% Computational Cost: add. (456->86), mult. (1015->151), div. (0->0), fcn. (1176->8), ass. (0->47)
t36 = sin(qJ(4));
t37 = sin(qJ(3));
t40 = cos(qJ(4));
t41 = cos(qJ(3));
t18 = t36 * t37 - t40 * t41;
t28 = -t41 * pkin(3) - pkin(2);
t12 = t18 * pkin(4) + t28;
t56 = 0.2e1 * t12;
t55 = 0.2e1 * t28;
t54 = 0.2e1 * t41;
t53 = -pkin(7) - pkin(6);
t35 = sin(qJ(5));
t52 = t35 * pkin(4);
t51 = t36 * pkin(3);
t38 = sin(qJ(2));
t50 = t37 * t38;
t49 = t41 * t38;
t31 = t37 ^ 2;
t33 = t41 ^ 2;
t48 = t31 + t33;
t39 = cos(qJ(5));
t47 = t39 * t51;
t46 = t48 * t38;
t22 = t53 * t37;
t23 = t53 * t41;
t10 = t40 * t22 + t36 * t23;
t30 = t40 * pkin(3);
t27 = t30 + pkin(4);
t15 = t39 * t27 - t35 * t51;
t11 = t36 * t22 - t40 * t23;
t20 = t36 * t41 + t40 * t37;
t42 = cos(qJ(2));
t34 = t42 ^ 2;
t32 = t38 ^ 2;
t29 = t39 * pkin(4);
t16 = t35 * t27 + t47;
t14 = -t36 * t50 + t40 * t49;
t13 = t20 * t38;
t9 = -t35 * t18 + t39 * t20;
t7 = t39 * t18 + t35 * t20;
t6 = -t18 * pkin(8) + t11;
t5 = -t20 * pkin(8) + t10;
t4 = -t35 * t13 + t39 * t14;
t3 = -t39 * t13 - t35 * t14;
t2 = t35 * t5 + t39 * t6;
t1 = -t35 * t6 + t39 * t5;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 + t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t32 + t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t14 ^ 2 + t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t38, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t41, -t42 * t37, t46, t42 * pkin(2) + pkin(6) * t46, 0, 0, 0, 0, 0, 0, -t42 * t18, -t42 * t20, t13 * t20 - t14 * t18, -t13 * t10 + t14 * t11 - t42 * t28, 0, 0, 0, 0, 0, 0, -t42 * t7, -t42 * t9, -t3 * t9 - t4 * t7, t3 * t1 - t42 * t12 + t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t31, t37 * t54, 0, t33, 0, 0, pkin(2) * t54, -0.2e1 * pkin(2) * t37, 0.2e1 * t48 * pkin(6), t48 * pkin(6) ^ 2 + pkin(2) ^ 2, t20 ^ 2, -0.2e1 * t20 * t18, 0, t18 ^ 2, 0, 0, t18 * t55, t20 * t55, -0.2e1 * t10 * t20 - 0.2e1 * t11 * t18, t10 ^ 2 + t11 ^ 2 + t28 ^ 2, t9 ^ 2, -0.2e1 * t9 * t7, 0, t7 ^ 2, 0, 0, t7 * t56, t9 * t56, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t49, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, (-t13 * t40 + t14 * t36) * pkin(3), 0, 0, 0, 0, 0, 0, t3, -t4, 0, t3 * t15 + t4 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t41, 0, -t37 * pkin(6), -t41 * pkin(6), 0, 0, 0, 0, t20, 0, -t18, 0, t10, -t11, (-t18 * t36 - t20 * t40) * pkin(3), (t10 * t40 + t11 * t36) * pkin(3), 0, 0, t9, 0, -t7, 0, t1, -t2, -t15 * t9 - t16 * t7, t1 * t15 + t2 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t51, 0, (t36 ^ 2 + t40 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t15, -0.2e1 * t16, 0, t15 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, (t3 * t39 + t35 * t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, 0, t10, -t11, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, (-t35 * t7 - t39 * t9) * pkin(4), (t1 * t39 + t2 * t35) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t30, -t51, 0, 0, 0, 0, 0, 0, 0, 1, t15 + t29, -t47 + (-pkin(4) - t27) * t35, 0, (t15 * t39 + t16 * t35) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t52, 0, (t35 ^ 2 + t39 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t15, -t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t52, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
