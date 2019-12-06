% Calculate minimal parameter regressor of gravitation load for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t21 = cos(pkin(8));
t24 = sin(qJ(3));
t22 = cos(pkin(7));
t26 = cos(qJ(3));
t32 = t22 * t26;
t33 = t22 * t24;
t20 = sin(pkin(7));
t34 = t20 * t26;
t35 = t20 * t24;
t19 = sin(pkin(8));
t39 = g(3) * t19;
t28 = t24 * t39 - g(1) * (-t21 * t33 + t34) - g(2) * (-t21 * t35 - t32);
t23 = sin(qJ(4));
t38 = t19 * t23;
t25 = cos(qJ(4));
t37 = t19 * t25;
t36 = t19 * t26;
t11 = t21 * t32 + t35;
t9 = t21 * t34 - t33;
t31 = -g(1) * t11 - g(2) * t9;
t12 = t21 * t25 + t23 * t36;
t4 = -t20 * t37 + t9 * t23;
t6 = t11 * t23 - t22 * t37;
t1 = g(1) * t6 + g(2) * t4 + g(3) * t12;
t13 = -t21 * t23 + t25 * t36;
t5 = t20 * t38 + t9 * t25;
t7 = t11 * t25 + t22 * t38;
t29 = g(1) * t7 + g(2) * t5 + g(3) * t13;
t27 = g(3) * t36 - t31;
t14 = -g(1) * t20 + g(2) * t22;
t3 = t28 * t25;
t2 = t28 * t23;
t8 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, t28, t27, 0, 0, 0, 0, 0, t3, -t2, t3, -t27, t2, (-t26 * t39 + t31) * pkin(6) + t28 * (pkin(4) * t25 + qJ(5) * t23 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t29, t1, 0, -t29, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t12 * pkin(4) + t13 * qJ(5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t8;
