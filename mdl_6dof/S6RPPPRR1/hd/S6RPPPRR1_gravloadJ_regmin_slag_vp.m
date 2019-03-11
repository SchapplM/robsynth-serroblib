% Calculate minimal parameter regressor of gravitation load for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(qJ(5));
t19 = cos(qJ(5));
t14 = qJ(1) + pkin(9);
t11 = sin(t14);
t12 = cos(t14);
t23 = g(1) * t12 + g(2) * t11;
t30 = -g(3) * t16 + t23 * t19;
t28 = g(3) * t19;
t15 = sin(qJ(6));
t27 = t15 * t16;
t18 = cos(qJ(6));
t26 = t16 * t18;
t20 = cos(qJ(1));
t25 = t20 * pkin(1) + t12 * pkin(2) + t11 * qJ(3);
t17 = sin(qJ(1));
t24 = -t17 * pkin(1) + t12 * qJ(3);
t5 = g(1) * t11 - g(2) * t12;
t22 = g(1) * t17 - g(2) * t20;
t4 = -t11 * t15 + t12 * t26;
t3 = -t11 * t18 - t12 * t27;
t2 = -t11 * t26 - t12 * t15;
t1 = t11 * t27 - t12 * t18;
t6 = [0, t22, g(1) * t20 + g(2) * t17, t22 * pkin(1), -t5, -t23, -g(1) * (-t11 * pkin(2) + t24) - g(2) * t25, -t23, t5, -g(1) * ((-pkin(2) - qJ(4)) * t11 + t24) - g(2) * (t12 * qJ(4) + t25) 0, 0, 0, 0, 0, t5 * t16, t5 * t19, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t5, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t23 * t16 + t28, 0, 0, 0, 0, 0, -t30 * t18, t30 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t15 * t28, g(1) * t4 - g(2) * t2 + t18 * t28;];
taug_reg  = t6;
