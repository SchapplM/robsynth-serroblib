% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRP2
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
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:48:06
% EndTime: 2019-05-05 14:48:06
% DurationCPUTime: 0.19s
% Computational Cost: add. (284->59), mult. (257->77), div. (0->0), fcn. (271->10), ass. (0->38)
t21 = qJ(1) + pkin(9);
t16 = sin(t21);
t18 = cos(t21);
t35 = g(1) * t18 + g(2) * t16;
t20 = pkin(10) + qJ(4);
t15 = sin(t20);
t42 = g(3) * t15;
t26 = sin(qJ(1));
t41 = t26 * pkin(1);
t25 = sin(qJ(5));
t40 = t16 * t25;
t27 = cos(qJ(5));
t39 = t16 * t27;
t38 = t18 * t25;
t37 = t18 * t27;
t17 = cos(t20);
t7 = t17 * t40 + t37;
t9 = t17 * t38 - t39;
t36 = g(1) * t7 - g(2) * t9;
t34 = g(1) * t16 - g(2) * t18;
t28 = cos(qJ(1));
t33 = g(1) * t26 - g(2) * t28;
t23 = cos(pkin(10));
t32 = t23 * pkin(3) + t17 * pkin(4) + t15 * pkin(8) + pkin(2);
t31 = pkin(5) * t27 + qJ(6) * t25 + pkin(4);
t1 = g(1) * t9 + g(2) * t7 + t25 * t42;
t10 = t17 * t37 + t40;
t8 = t17 * t39 - t38;
t30 = g(1) * t10 + g(2) * t8 + t27 * t42;
t29 = -g(3) * t17 + t35 * t15;
t24 = -pkin(7) - qJ(3);
t19 = t28 * pkin(1);
t6 = t34 * t15;
t5 = t35 * t17 + t42;
t4 = t29 * t27;
t3 = t29 * t25;
t2 = g(1) * t8 - g(2) * t10;
t11 = [0, t33, g(1) * t28 + g(2) * t26, t33 * pkin(1), t34 * t23, -t34 * sin(pkin(10)) -t35, -g(1) * (-t16 * pkin(2) + t18 * qJ(3) - t41) - g(2) * (t18 * pkin(2) + t16 * qJ(3) + t19) 0, 0, 0, 0, 0, t34 * t17, -t6, 0, 0, 0, 0, 0, t2, -t36, t2, t6, t36, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) - t41) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t19) + (g(1) * t24 - g(2) * t32) * t18 + (g(1) * t32 + g(2) * t24) * t16; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t5, t3 (-t35 * pkin(8) - g(3) * t31) * t17 + (-g(3) * pkin(8) + t35 * t31) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t30, t1, 0, -t30, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (-pkin(5) * t25 + qJ(6) * t27) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t11;
