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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:43:15
% EndTime: 2022-01-20 11:43:16
% DurationCPUTime: 0.19s
% Computational Cost: add. (251->52), mult. (193->60), div. (0->0), fcn. (173->10), ass. (0->40)
t36 = sin(qJ(1));
t48 = t36 * pkin(1);
t34 = -qJ(4) - pkin(7);
t32 = qJ(3) + pkin(9);
t25 = cos(t32);
t37 = cos(qJ(3));
t29 = t37 * pkin(3);
t47 = pkin(4) * t25 + t29;
t33 = qJ(1) + qJ(2);
t27 = sin(t33);
t28 = cos(t33);
t46 = t28 * pkin(2) + t27 * pkin(7);
t45 = -t27 * pkin(2) + t28 * pkin(7);
t12 = pkin(2) + t47;
t31 = pkin(8) - t34;
t44 = t28 * t12 + t27 * t31;
t43 = -t27 * t12 + t31 * t28;
t23 = t29 + pkin(2);
t42 = t28 * t23 - t27 * t34;
t11 = g(1) * t28 + g(2) * t27;
t10 = g(1) * t27 - g(2) * t28;
t38 = cos(qJ(1));
t41 = g(1) * t36 - g(2) * t38;
t40 = -t27 * t23 - t28 * t34;
t35 = sin(qJ(3));
t39 = -g(3) * t37 + t11 * t35;
t30 = t38 * pkin(1);
t26 = qJ(5) + t32;
t24 = sin(t32);
t18 = cos(t26);
t17 = sin(t26);
t8 = t10 * t37;
t7 = t10 * t35;
t6 = t10 * t25;
t5 = t10 * t24;
t4 = t10 * t18;
t3 = t10 * t17;
t2 = g(3) * t17 + t11 * t18;
t1 = -g(3) * t18 + t11 * t17;
t9 = [0, 0, 0, 0, 0, 0, t41, g(1) * t38 + g(2) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, t41 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * (t45 - t48) - g(2) * (t30 + t46), 0, 0, 0, 0, 0, 0, t6, -t5, -t11, -g(1) * (t40 - t48) - g(2) * (t30 + t42), 0, 0, 0, 0, 0, 0, t4, -t3, -t11, -g(1) * (t43 - t48) - g(2) * (t30 + t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * t45 - g(2) * t46, 0, 0, 0, 0, 0, 0, t6, -t5, -t11, -g(1) * t40 - g(2) * t42, 0, 0, 0, 0, 0, 0, t4, -t3, -t11, -g(1) * t43 - g(2) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, g(3) * t35 + t11 * t37, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t25 + t11 * t24, g(3) * t24 + t11 * t25, 0, t39 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t47 - t11 * (-t35 * pkin(3) - pkin(4) * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
