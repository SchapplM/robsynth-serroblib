% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:29
% EndTime: 2020-01-03 12:09:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (251->52), mult. (193->60), div. (0->0), fcn. (173->10), ass. (0->40)
t36 = qJ(3) + pkin(9);
t28 = cos(t36);
t41 = cos(qJ(3));
t33 = t41 * pkin(3);
t49 = pkin(4) * t28 + t33;
t13 = pkin(2) + t49;
t37 = qJ(1) + qJ(2);
t30 = sin(t37);
t31 = cos(t37);
t38 = -qJ(4) - pkin(7);
t35 = -pkin(8) + t38;
t51 = t30 * t13 + t31 * t35;
t26 = t33 + pkin(2);
t50 = t30 * t26 + t31 * t38;
t48 = t31 * pkin(2) + t30 * pkin(7);
t47 = t30 * pkin(2) - t31 * pkin(7);
t46 = t31 * t13 - t30 * t35;
t45 = t31 * t26 - t30 * t38;
t12 = g(2) * t31 + g(3) * t30;
t11 = g(2) * t30 - g(3) * t31;
t40 = sin(qJ(1));
t42 = cos(qJ(1));
t44 = -g(2) * t42 - g(3) * t40;
t39 = sin(qJ(3));
t43 = -g(1) * t41 + t11 * t39;
t34 = t42 * pkin(1);
t32 = t40 * pkin(1);
t29 = qJ(5) + t36;
t27 = sin(t36);
t21 = cos(t29);
t20 = sin(t29);
t8 = t12 * t41;
t7 = t12 * t39;
t6 = t12 * t28;
t5 = t12 * t27;
t4 = t12 * t21;
t3 = t12 * t20;
t2 = -g(1) * t21 + t11 * t20;
t1 = g(1) * t20 + t11 * t21;
t9 = [0, 0, 0, 0, 0, 0, t44, g(2) * t40 - g(3) * t42, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, t44 * pkin(1), 0, 0, 0, 0, 0, 0, -t8, t7, -t11, -g(2) * (t34 + t48) - g(3) * (t32 + t47), 0, 0, 0, 0, 0, 0, -t6, t5, -t11, -g(2) * (t34 + t45) - g(3) * (t32 + t50), 0, 0, 0, 0, 0, 0, -t4, t3, -t11, -g(2) * (t34 + t46) - g(3) * (t32 + t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, -t11, -g(2) * t48 - g(3) * t47, 0, 0, 0, 0, 0, 0, -t6, t5, -t11, -g(2) * t45 - g(3) * t50, 0, 0, 0, 0, 0, 0, -t4, t3, -t11, -g(2) * t46 - g(3) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, g(1) * t39 + t11 * t41, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t28 + t11 * t27, g(1) * t27 + t11 * t28, 0, t43 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t49 - t11 * (-t39 * pkin(3) - pkin(4) * t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t9;
