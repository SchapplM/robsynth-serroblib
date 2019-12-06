% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t18 = qJ(1) + qJ(2);
t16 = sin(t18);
t34 = pkin(2) * t16;
t17 = cos(t18);
t33 = pkin(2) * t17;
t19 = sin(pkin(9));
t32 = g(1) * t19;
t22 = sin(qJ(1));
t31 = t22 * pkin(1);
t24 = cos(qJ(1));
t30 = t24 * pkin(1);
t20 = cos(pkin(9));
t21 = sin(qJ(5));
t29 = t20 * t21;
t23 = cos(qJ(5));
t28 = t20 * t23;
t15 = pkin(8) + t18;
t13 = sin(t15);
t14 = cos(t15);
t27 = g(2) * t14 + g(3) * t13;
t11 = g(2) * t17 + g(3) * t16;
t26 = -t13 * pkin(3) + t14 * qJ(4) - t34;
t25 = -t14 * pkin(3) - t13 * qJ(4) - t33;
t10 = -g(2) * t16 + g(3) * t17;
t9 = g(2) * t13 - g(3) * t14;
t8 = t27 * t20;
t7 = t27 * t19;
t6 = -t13 * t21 - t14 * t28;
t5 = -t13 * t23 + t14 * t29;
t4 = t13 * t28 - t14 * t21;
t3 = t13 * t29 + t14 * t23;
t2 = -g(2) * t6 + g(3) * t4;
t1 = -g(2) * t5 - g(3) * t3;
t12 = [0, g(2) * t24 + g(3) * t22, -g(2) * t22 + g(3) * t24, 0, t11, t10, -g(2) * (-t30 - t33) - g(3) * (-t31 - t34), t8, -t7, t9, -g(2) * (t25 - t30) - g(3) * (t26 - t31), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, t11, t10, t11 * pkin(2), t8, -t7, t9, -g(2) * t25 - g(3) * t26, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t3 + g(3) * t5 + t21 * t32, -g(2) * t4 - g(3) * t6 + t23 * t32;];
taug_reg = t12;
