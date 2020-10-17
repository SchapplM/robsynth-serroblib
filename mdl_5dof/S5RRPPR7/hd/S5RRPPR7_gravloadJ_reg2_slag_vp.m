% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:22
% EndTime: 2019-12-31 19:36:23
% DurationCPUTime: 0.28s
% Computational Cost: add. (180->65), mult. (242->81), div. (0->0), fcn. (233->8), ass. (0->43)
t24 = qJ(2) + pkin(8);
t20 = sin(t24);
t21 = cos(t24);
t54 = t21 * pkin(3) + t20 * qJ(4);
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t13 = g(1) * t31 + g(2) * t28;
t53 = t13 * t20;
t2 = g(3) * t20 + t13 * t21;
t27 = sin(qJ(2));
t50 = pkin(2) * t27;
t46 = g(3) * t21;
t45 = t21 * pkin(7);
t25 = -qJ(3) - pkin(6);
t44 = pkin(4) - t25;
t26 = sin(qJ(5));
t43 = t28 * t26;
t29 = cos(qJ(5));
t42 = t28 * t29;
t41 = t31 * t25;
t40 = t31 * t26;
t39 = t31 * t29;
t38 = qJ(4) * t21;
t30 = cos(qJ(2));
t22 = t30 * pkin(2);
t37 = t22 + t54;
t19 = t22 + pkin(1);
t15 = t31 * t19;
t36 = g(2) * (t54 * t31 + t15);
t35 = -pkin(3) * t20 - t50;
t12 = g(1) * t28 - g(2) * t31;
t34 = -t19 - t54;
t32 = -g(3) * t30 + t13 * t27;
t11 = t31 * t38;
t9 = t28 * t38;
t8 = -t20 * t43 + t39;
t7 = t20 * t42 + t40;
t6 = t20 * t40 + t42;
t5 = t20 * t39 - t43;
t4 = t12 * t21;
t3 = t12 * t20;
t1 = -t46 + t53;
t10 = [0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t30, -t12 * t27, -t13, -g(1) * (-t28 * pkin(1) + t31 * pkin(6)) - g(2) * (t31 * pkin(1) + t28 * pkin(6)), 0, 0, 0, 0, 0, 0, t4, -t3, -t13, -g(1) * (-t28 * t19 - t41) - g(2) * (-t28 * t25 + t15), 0, 0, 0, 0, 0, 0, -t13, -t4, t3, g(1) * t41 - t36 + (-g(1) * t34 + g(2) * t25) * t28, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5, t4, -t36 + (-g(1) * t44 - g(2) * t45) * t31 + (-g(1) * (t34 - t45) - g(2) * t44) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t27 + t13 * t30, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t32 * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (t31 * t35 + t11) - g(2) * (t28 * t35 + t9) - g(3) * t37, 0, 0, 0, 0, 0, 0, -t2 * t26, -t2 * t29, t1, -g(1) * (-t31 * t50 + t11) - g(2) * (-t28 * t50 + t9) - g(3) * (t37 + t45) + (pkin(3) + pkin(7)) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t29 * t46, g(1) * t6 - g(2) * t8 - t26 * t46, 0, 0;];
taug_reg = t10;
