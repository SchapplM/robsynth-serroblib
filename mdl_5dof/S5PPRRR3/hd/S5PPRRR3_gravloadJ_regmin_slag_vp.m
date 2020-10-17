% Calculate minimal parameter regressor of gravitation load for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:59
% EndTime: 2019-12-05 15:17:00
% DurationCPUTime: 0.12s
% Computational Cost: add. (99->34), mult. (182->62), div. (0->0), fcn. (224->10), ass. (0->26)
t10 = sin(pkin(9));
t11 = sin(pkin(8));
t27 = t10 * t11;
t13 = cos(pkin(8));
t26 = t10 * t13;
t14 = sin(qJ(4));
t25 = t10 * t14;
t16 = cos(qJ(4));
t24 = t10 * t16;
t17 = cos(qJ(3));
t23 = t10 * t17;
t15 = sin(qJ(3));
t22 = t11 * t15;
t21 = t11 * t17;
t20 = t13 * t15;
t19 = t13 * t17;
t12 = cos(pkin(9));
t18 = g(3) * t10 * t15 - g(1) * (-t12 * t20 + t21) - g(2) * (-t12 * t22 - t19);
t9 = qJ(4) + qJ(5);
t8 = cos(t9);
t7 = sin(t9);
t6 = t12 * t19 + t22;
t4 = t12 * t21 - t20;
t2 = -g(1) * (-t7 * t26 - t6 * t8) - g(2) * (-t7 * t27 - t4 * t8) - g(3) * (t12 * t7 - t8 * t23);
t1 = -g(1) * (t8 * t26 - t6 * t7) - g(2) * (t8 * t27 - t4 * t7) - g(3) * (-t12 * t8 - t7 * t23);
t3 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(1) * t11 + g(2) * t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t18, g(1) * t6 + g(2) * t4 + g(3) * t23, 0, 0, 0, 0, 0, t18 * t16, -t18 * t14, 0, 0, 0, 0, 0, t18 * t8, -t18 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t24 - t6 * t14) - g(2) * (t11 * t24 - t4 * t14) - g(3) * (-t12 * t16 - t14 * t23), -g(1) * (-t13 * t25 - t6 * t16) - g(2) * (-t11 * t25 - t4 * t16) - g(3) * (t12 * t14 - t16 * t23), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t3;
