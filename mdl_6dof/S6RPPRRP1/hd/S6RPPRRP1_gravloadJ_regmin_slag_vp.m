% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:44:52
% EndTime: 2019-05-05 14:44:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (192->48), mult. (169->63), div. (0->0), fcn. (165->10), ass. (0->35)
t17 = qJ(1) + pkin(9);
t12 = sin(t17);
t14 = cos(t17);
t31 = g(1) * t14 + g(2) * t12;
t22 = sin(qJ(5));
t16 = pkin(10) + qJ(4);
t13 = cos(t16);
t24 = cos(qJ(5));
t34 = t14 * t24;
t37 = t12 * t22;
t4 = t13 * t37 + t34;
t11 = sin(t16);
t40 = g(3) * t11;
t35 = t14 * t22;
t36 = t12 * t24;
t6 = -t13 * t35 + t36;
t45 = -g(1) * t6 + g(2) * t4 + t22 * t40;
t1 = -g(3) * t13 + t31 * t11;
t23 = sin(qJ(1));
t38 = t23 * pkin(1);
t32 = pkin(5) * t22 + pkin(7) + qJ(3);
t30 = g(1) * t12 - g(2) * t14;
t25 = cos(qJ(1));
t29 = g(1) * t23 - g(2) * t25;
t10 = t24 * pkin(5) + pkin(4);
t20 = -qJ(6) - pkin(8);
t28 = t13 * t10 - t11 * t20;
t19 = cos(pkin(10));
t26 = t19 * pkin(3) + pkin(2) + t28;
t15 = t25 * pkin(1);
t7 = t13 * t34 + t37;
t5 = -t13 * t36 + t35;
t3 = t30 * t11;
t2 = t31 * t13 + t40;
t8 = [0, t29, g(1) * t25 + g(2) * t23, t29 * pkin(1), t30 * t19, -t30 * sin(pkin(10)) -t31, -g(1) * (-t12 * pkin(2) + t14 * qJ(3) - t38) - g(2) * (t14 * pkin(2) + t12 * qJ(3) + t15) 0, 0, 0, 0, 0, t30 * t13, -t3, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, -g(1) * t4 - g(2) * t6, t3, g(1) * t38 - g(2) * t15 + (-g(1) * t32 - g(2) * t26) * t14 + (g(1) * t26 - g(2) * t32) * t12; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, t1 * t24, -t1 * t22, -t2, -g(3) * t28 + t31 * (t10 * t11 + t13 * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(1) * t7 - g(2) * t5 + t24 * t40, 0, t45 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t8;
