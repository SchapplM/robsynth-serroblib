% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% taug_reg [5x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_gravloadJ_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = sin(qJ(1));
t17 = cos(qJ(1));
t37 = -g(1) * t17 - g(2) * t13;
t12 = sin(qJ(3));
t16 = cos(qJ(3));
t18 = -g(3) * t16 - t12 * t37;
t34 = g(3) * t12;
t10 = sin(qJ(5));
t32 = t12 * t10;
t14 = cos(qJ(5));
t31 = t12 * t14;
t15 = cos(qJ(4));
t30 = t12 * t15;
t29 = t12 * t17;
t28 = t13 * t16;
t27 = t16 * t10;
t26 = t16 * t14;
t11 = sin(qJ(4));
t25 = t17 * t11;
t24 = t17 * t15;
t7 = g(1) * t13 - g(2) * t17;
t4 = t15 * t28 - t25;
t23 = -t4 * t10 + t13 * t31;
t22 = t13 * t32 + t4 * t14;
t21 = t14 * t30 - t27;
t20 = t10 * t30 + t26;
t3 = t11 * t28 + t24;
t5 = t13 * t15 - t16 * t25;
t19 = -g(1) * t5 + g(2) * t3 + t11 * t34;
t6 = t13 * t11 + t16 * t24;
t2 = t10 * t29 + t6 * t14;
t1 = -t6 * t10 + t14 * t29;
t8 = [0, t7, -t37, t7, t37, t37 * qJ(2), 0, 0, 0, 0, 0, t7 * t16, -t7 * t12, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t22 - g(2) * t2, g(1) * t23 - g(2) * t1; 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t16 * t37 + t34, 0, 0, 0, 0, 0, t18 * t15, -t18 * t11, 0, 0, 0, 0, 0, -g(3) * (t15 * t26 + t32) - t37 * t21, -g(3) * (-t15 * t27 + t31) + t37 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(1) * t6 + g(2) * t4 + g(3) * t30, 0, 0, 0, 0, 0, t19 * t14, -t19 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t23 + g(3) * t20, g(1) * t2 + g(2) * t22 + g(3) * t21;];
taug_reg = t8;
