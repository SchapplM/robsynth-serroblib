% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:02
% EndTime: 2020-01-03 11:43:03
% DurationCPUTime: 0.34s
% Computational Cost: add. (204->74), mult. (271->110), div. (0->0), fcn. (285->10), ass. (0->60)
t36 = cos(pkin(8));
t40 = cos(qJ(3));
t41 = cos(qJ(1));
t46 = t41 * t40;
t38 = sin(qJ(3));
t39 = sin(qJ(1));
t53 = t39 * t38;
t12 = -t36 * t53 - t46;
t47 = t41 * t38;
t52 = t39 * t40;
t14 = t36 * t47 - t52;
t35 = sin(pkin(8));
t64 = g(1) * t35;
t65 = -g(2) * t12 - g(3) * t14 + t38 * t64;
t61 = t38 * pkin(3);
t60 = t35 * t39;
t59 = t36 * t39;
t34 = qJ(3) + pkin(9);
t26 = sin(t34);
t18 = pkin(4) * t26 + t61;
t58 = t39 * t18;
t28 = qJ(5) + t34;
t23 = sin(t28);
t57 = t39 * t23;
t24 = cos(t28);
t56 = t39 * t24;
t55 = t39 * t26;
t27 = cos(t34);
t54 = t39 * t27;
t51 = t41 * t23;
t50 = t41 * t24;
t49 = t41 * t26;
t48 = t41 * t27;
t37 = -qJ(4) - pkin(6);
t31 = t40 * pkin(3);
t19 = pkin(4) * t27 + t31;
t45 = t41 * pkin(1) + t39 * qJ(2);
t30 = t39 * pkin(1);
t43 = -t41 * qJ(2) + t30;
t42 = pkin(2) * t36 + pkin(6) * t35;
t21 = g(2) * t41 + g(3) * t39;
t20 = g(2) * t39 - g(3) * t41;
t33 = -pkin(7) + t37;
t25 = t31 + pkin(2);
t17 = pkin(2) + t19;
t16 = t21 * t35;
t15 = t36 * t46 + t53;
t13 = t36 * t52 - t47;
t11 = g(1) * t36 - t20 * t35;
t10 = t36 * t48 + t55;
t9 = t36 * t49 - t54;
t8 = t36 * t54 - t49;
t7 = -t36 * t55 - t48;
t6 = t36 * t50 + t57;
t5 = t36 * t51 - t56;
t4 = t36 * t56 - t51;
t3 = -t36 * t57 - t50;
t2 = g(2) * t4 - g(3) * t6 + t24 * t64;
t1 = -g(2) * t3 - g(3) * t5 + t23 * t64;
t22 = [0, 0, 0, 0, 0, 0, -t21, t20, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t36, t16, -t20, -g(2) * t45 - g(3) * t43, 0, 0, 0, 0, 0, 0, -g(2) * t15 - g(3) * t13, g(2) * t14 - g(3) * t12, -t16, -g(2) * (t42 * t41 + t45) - g(3) * (t42 * t39 + t43), 0, 0, 0, 0, 0, 0, -g(2) * t10 - g(3) * t8, g(2) * t9 - g(3) * t7, -t16, -g(2) * (pkin(3) * t53 + t45) - g(3) * (t25 * t59 - t37 * t60 + t30) + (-g(2) * (t25 * t36 - t35 * t37) - g(3) * (-qJ(2) - t61)) * t41, 0, 0, 0, 0, 0, 0, -g(2) * t6 - g(3) * t4, g(2) * t5 - g(3) * t3, -t16, -g(2) * (t45 + t58) - g(3) * (t17 * t59 - t33 * t60 + t30) + (-g(2) * (t17 * t36 - t33 * t35) - g(3) * (-qJ(2) - t18)) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, g(2) * t13 - g(3) * t15 + t40 * t64, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t7 - g(3) * t9 + t26 * t64, g(2) * t8 - g(3) * t10 + t27 * t64, 0, t65 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, t18 * t64 - g(2) * (-t41 * t19 - t36 * t58) - g(3) * (t41 * t36 * t18 - t39 * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t22;
