% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t21 = qJ(2) + qJ(3);
t20 = qJ(4) + t21;
t17 = qJ(5) + t20;
t13 = sin(t17);
t14 = cos(t17);
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t29 = g(1) * t27 + g(2) * t24;
t3 = -g(3) * t14 + t29 * t13;
t35 = g(3) * t13;
t22 = sin(qJ(6));
t33 = t24 * t22;
t25 = cos(qJ(6));
t32 = t24 * t25;
t31 = t27 * t22;
t30 = t27 * t25;
t28 = g(1) * t24 - g(2) * t27;
t26 = cos(qJ(2));
t23 = sin(qJ(2));
t19 = cos(t21);
t18 = sin(t21);
t16 = cos(t20);
t15 = sin(t20);
t12 = t14 * t30 + t33;
t11 = -t14 * t31 + t32;
t10 = -t14 * t32 + t31;
t9 = t14 * t33 + t30;
t8 = g(3) * t18 + t29 * t19;
t7 = -g(3) * t19 + t29 * t18;
t6 = g(3) * t15 + t29 * t16;
t5 = -g(3) * t16 + t29 * t15;
t4 = t29 * t14 + t35;
t2 = t3 * t25;
t1 = t3 * t22;
t34 = [0, t28, t29, 0, 0, 0, 0, 0, t28 * t26, -t28 * t23, 0, 0, 0, 0, 0, t28 * t19, -t28 * t18, 0, 0, 0, 0, 0, t28 * t16, -t28 * t15, 0, 0, 0, 0, 0, t28 * t14, -t28 * t13, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 + t29 * t23, g(3) * t23 + t29 * t26, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t22 * t35, g(1) * t12 - g(2) * t10 + t25 * t35;];
taug_reg  = t34;
