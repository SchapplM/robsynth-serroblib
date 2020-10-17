% Calculate minimal parameter regressor of gravitation load for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:12
% EndTime: 2019-12-05 17:10:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (139->28), mult. (154->42), div. (0->0), fcn. (166->10), ass. (0->31)
t14 = qJ(2) + qJ(3);
t10 = sin(t14);
t12 = cos(t14);
t15 = sin(pkin(9));
t16 = cos(pkin(9));
t22 = g(1) * t16 + g(2) * t15;
t7 = -g(3) * t12 + t22 * t10;
t32 = g(3) * t10;
t13 = qJ(4) + qJ(5);
t9 = sin(t13);
t30 = t15 * t9;
t29 = t16 * t9;
t11 = cos(t13);
t28 = t15 * t11;
t17 = sin(qJ(4));
t27 = t15 * t17;
t19 = cos(qJ(4));
t26 = t15 * t19;
t25 = t16 * t11;
t24 = t16 * t17;
t23 = t16 * t19;
t20 = cos(qJ(2));
t18 = sin(qJ(2));
t8 = t22 * t12 + t32;
t6 = t7 * t19;
t5 = t7 * t17;
t4 = t7 * t11;
t3 = t7 * t9;
t2 = -g(1) * (-t12 * t25 - t30) - g(2) * (-t12 * t28 + t29) + t11 * t32;
t1 = -g(1) * (-t12 * t29 + t28) - g(2) * (-t12 * t30 - t25) + t9 * t32;
t21 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(3) * t20 + t22 * t18, g(3) * t18 + t22 * t20, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t24 + t26) - g(2) * (-t12 * t27 - t23) + t17 * t32, -g(1) * (-t12 * t23 - t27) - g(2) * (-t12 * t26 + t24) + t19 * t32, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t21;
