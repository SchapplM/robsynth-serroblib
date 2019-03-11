% Calculate inertial parameters regressor of gravitation load for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t15 = qJ(1) + pkin(7);
t12 = cos(t15);
t17 = cos(qJ(1));
t20 = t17 * pkin(1) + pkin(2) * t12;
t13 = qJ(3) + t15;
t8 = sin(t13);
t9 = cos(t13);
t3 = g(1) * t8 - g(2) * t9;
t11 = sin(t15);
t16 = sin(qJ(1));
t19 = -t16 * pkin(1) - pkin(2) * t11;
t18 = g(1) * t16 - g(2) * t17;
t10 = qJ(4) + t13;
t6 = cos(t10);
t5 = sin(t10);
t4 = g(1) * t9 + g(2) * t8;
t2 = g(1) * t6 + g(2) * t5;
t1 = g(1) * t5 - g(2) * t6;
t7 = [0, 0, 0, 0, 0, 0, t18, g(1) * t17 + g(2) * t16, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t12, g(1) * t12 + g(2) * t11, 0, t18 * pkin(1), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t19 - g(2) * t20, 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-pkin(3) * t8 + t19) - g(2) * (pkin(3) * t9 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t7;
