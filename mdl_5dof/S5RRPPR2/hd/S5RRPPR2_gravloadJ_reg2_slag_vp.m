% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:51
% EndTime: 2022-01-20 10:05:52
% DurationCPUTime: 0.19s
% Computational Cost: add. (252->49), mult. (170->62), div. (0->0), fcn. (168->10), ass. (0->43)
t27 = qJ(1) + qJ(2);
t24 = sin(t27);
t48 = pkin(2) * t24;
t29 = cos(pkin(9));
t47 = pkin(4) * t29;
t23 = pkin(8) + t27;
t20 = sin(t23);
t46 = g(1) * t20;
t28 = sin(pkin(9));
t45 = g(3) * t28;
t31 = sin(qJ(1));
t44 = t31 * pkin(1);
t21 = cos(t23);
t43 = t21 * t28;
t30 = sin(qJ(5));
t42 = t29 * t30;
t32 = cos(qJ(5));
t41 = t29 * t32;
t25 = cos(t27);
t22 = pkin(2) * t25;
t40 = t21 * pkin(3) + t20 * qJ(4) + t22;
t17 = t21 * qJ(4);
t39 = t17 - t48;
t38 = pkin(7) * t43 + t21 * t47 + t40;
t37 = -t44 - t48;
t9 = -g(2) * t21 + t46;
t11 = g(1) * t24 - g(2) * t25;
t33 = cos(qJ(1));
t36 = g(1) * t31 - g(2) * t33;
t35 = -t20 * pkin(3) + t39;
t34 = (-pkin(7) * t28 - pkin(3) - t47) * t46;
t26 = t33 * pkin(1);
t12 = g(1) * t25 + g(2) * t24;
t10 = g(1) * t21 + g(2) * t20;
t8 = t9 * t29;
t7 = -g(2) * t43 + t28 * t46;
t6 = t20 * t30 + t21 * t41;
t5 = t20 * t32 - t21 * t42;
t4 = -t20 * t41 + t21 * t30;
t3 = t20 * t42 + t21 * t32;
t2 = -g(1) * t4 - g(2) * t6;
t1 = -g(1) * t3 - g(2) * t5;
t13 = [0, 0, 0, 0, 0, 0, t36, g(1) * t33 + g(2) * t31, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t36 * pkin(1), 0, 0, 0, 0, 0, 0, t9, t10, 0, -g(1) * t37 - g(2) * (t22 + t26), 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (t35 - t44) - g(2) * (t26 + t40), 0, 0, 0, 0, 0, 0, t2, t1, t7, -g(1) * (t17 + t37) - g(2) * (t26 + t38) - t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t11 * pkin(2), 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * t35 - g(2) * t40, 0, 0, 0, 0, 0, 0, t2, t1, t7, -g(1) * t39 - g(2) * t38 - t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t30 * t45, g(1) * t6 - g(2) * t4 + t32 * t45, 0, 0;];
taug_reg = t13;
