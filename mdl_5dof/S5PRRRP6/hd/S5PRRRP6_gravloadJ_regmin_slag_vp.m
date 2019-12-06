% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t32 = g(1) * t21 + g(2) * t20;
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t30 = -g(3) * t25 + t32 * t23;
t40 = g(3) * t23;
t19 = qJ(3) + qJ(4);
t17 = sin(t19);
t38 = t17 * t23;
t18 = cos(t19);
t37 = t18 * t23;
t36 = t20 * t25;
t35 = t21 * t25;
t22 = sin(qJ(3));
t34 = t22 * t25;
t24 = cos(qJ(3));
t33 = t24 * t25;
t31 = t24 * pkin(3) + pkin(4) * t18 + qJ(5) * t17 + pkin(2);
t10 = t17 * t35 - t20 * t18;
t8 = t17 * t36 + t21 * t18;
t1 = g(1) * t10 + g(2) * t8 + g(3) * t38;
t11 = t20 * t17 + t18 * t35;
t9 = -t21 * t17 + t18 * t36;
t3 = g(1) * t11 + g(2) * t9 + g(3) * t37;
t28 = -g(1) * (-t10 * pkin(4) + t11 * qJ(5)) - g(2) * (-t8 * pkin(4) + t9 * qJ(5)) - g(3) * (-pkin(4) * t38 + qJ(5) * t37);
t27 = -g(1) * (t20 * t24 - t21 * t34) - g(2) * (-t20 * t34 - t21 * t24) + t22 * t40;
t26 = -pkin(7) - pkin(6);
t12 = t32 * t25 + t40;
t5 = t30 * t18;
t4 = t30 * t17;
t2 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t30, t12, 0, 0, 0, 0, 0, t30 * t24, -t30 * t22, 0, 0, 0, 0, 0, t5, -t4, t5, -t12, t4, -g(3) * (-t23 * t26 + t31 * t25) + t32 * (t31 * t23 + t25 * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -g(1) * (-t20 * t22 - t21 * t33) - g(2) * (-t20 * t33 + t21 * t22) + t24 * t40, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t27 * pkin(3) + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t2;
