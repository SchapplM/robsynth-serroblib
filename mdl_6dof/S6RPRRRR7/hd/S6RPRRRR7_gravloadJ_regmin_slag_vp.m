% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x34]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t11 = g(1) * t22 - g(2) * t25;
t19 = qJ(3) + qJ(4);
t17 = qJ(5) + t19;
t13 = sin(t17);
t14 = cos(t17);
t32 = -g(3) * t13 + t11 * t14;
t30 = g(3) * t14;
t20 = sin(qJ(6));
t29 = t22 * t20;
t23 = cos(qJ(6));
t28 = t22 * t23;
t27 = t25 * t20;
t26 = t25 * t23;
t12 = g(1) * t25 + g(2) * t22;
t24 = cos(qJ(3));
t21 = sin(qJ(3));
t16 = cos(t19);
t15 = sin(t19);
t10 = t13 * t26 - t29;
t9 = t13 * t27 + t28;
t8 = t13 * t28 + t27;
t7 = -t13 * t29 + t26;
t6 = g(3) * t15 - t11 * t16;
t5 = g(3) * t16 + t11 * t15;
t3 = t11 * t13 + t30;
t2 = t32 * t23;
t1 = t32 * t20;
t4 = [0, t11, t12, -t11, -t12, -g(1) * (-t22 * pkin(1) + t25 * qJ(2)) - g(2) * (t25 * pkin(1) + t22 * qJ(2)) 0, 0, 0, 0, 0, -t12 * t21, -t12 * t24, 0, 0, 0, 0, 0, -t12 * t15, -t12 * t16, 0, 0, 0, 0, 0, -t12 * t13, -t12 * t14, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7; 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t21 - t11 * t24, g(3) * t24 + t11 * t21, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, -t32, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, -t32, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t20 * t30, g(1) * t8 - g(2) * t10 + t23 * t30;];
taug_reg  = t4;
