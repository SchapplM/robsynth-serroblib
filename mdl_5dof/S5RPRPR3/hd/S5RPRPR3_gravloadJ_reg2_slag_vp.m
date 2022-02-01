% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR3
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:55
% EndTime: 2022-01-23 09:20:56
% DurationCPUTime: 0.20s
% Computational Cost: add. (236->48), mult. (158->61), div. (0->0), fcn. (158->10), ass. (0->39)
t25 = qJ(1) + pkin(8);
t23 = qJ(3) + t25;
t20 = cos(t23);
t26 = sin(pkin(9));
t41 = t20 * t26;
t27 = cos(pkin(9));
t45 = pkin(4) * t27;
t46 = pkin(7) * t41 + t20 * t45;
t19 = sin(t23);
t44 = g(1) * t19;
t43 = g(3) * t26;
t42 = t19 * pkin(3);
t28 = sin(qJ(5));
t40 = t27 * t28;
t30 = cos(qJ(5));
t39 = t27 * t30;
t38 = t20 * pkin(3) + t19 * qJ(4);
t22 = cos(t25);
t31 = cos(qJ(1));
t37 = t31 * pkin(1) + pkin(2) * t22;
t36 = t37 + t38;
t21 = sin(t25);
t29 = sin(qJ(1));
t35 = -t29 * pkin(1) - pkin(2) * t21;
t9 = -g(2) * t20 + t44;
t34 = g(1) * t29 - g(2) * t31;
t15 = t20 * qJ(4);
t33 = t15 + t35;
t32 = (-pkin(7) * t26 - pkin(3) - t45) * t44;
t10 = g(1) * t20 + g(2) * t19;
t8 = t9 * t27;
t7 = -g(2) * t41 + t26 * t44;
t6 = t19 * t28 + t20 * t39;
t5 = t19 * t30 - t20 * t40;
t4 = -t19 * t39 + t20 * t28;
t3 = t19 * t40 + t20 * t30;
t2 = -g(1) * t4 - g(2) * t6;
t1 = -g(1) * t3 - g(2) * t5;
t11 = [0, 0, 0, 0, 0, 0, t34, g(1) * t31 + g(2) * t29, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t21 - g(2) * t22, g(1) * t22 + g(2) * t21, 0, t34 * pkin(1), 0, 0, 0, 0, 0, 0, t9, t10, 0, -g(1) * t35 - g(2) * t37, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (t33 - t42) - g(2) * t36, 0, 0, 0, 0, 0, 0, t2, t1, t7, -g(1) * t33 - g(2) * (t36 + t46) - t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (t15 - t42) - g(2) * t38, 0, 0, 0, 0, 0, 0, t2, t1, t7, -g(1) * t15 - g(2) * (t38 + t46) - t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t28 * t43, g(1) * t6 - g(2) * t4 + t30 * t43, 0, 0;];
taug_reg = t11;
