% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = qJ(3) + qJ(4);
t18 = sin(t23);
t25 = cos(pkin(8));
t19 = cos(t23);
t29 = cos(qJ(1));
t34 = t29 * t19;
t27 = sin(qJ(1));
t39 = t27 * t18;
t3 = t25 * t39 + t34;
t24 = sin(pkin(8));
t41 = g(1) * t24;
t35 = t29 * t18;
t38 = t27 * t19;
t5 = t25 * t35 - t38;
t1 = -g(2) * t3 + g(3) * t5 + t18 * t41;
t26 = sin(qJ(3));
t13 = t26 * pkin(3) + pkin(4) * t18;
t40 = t13 * t25;
t37 = t27 * t26;
t28 = cos(qJ(3));
t36 = t27 * t28;
t33 = t29 * t26;
t32 = t29 * t28;
t14 = t28 * pkin(3) + pkin(4) * t19;
t16 = g(2) * t29 + g(3) * t27;
t15 = g(2) * t27 - g(3) * t29;
t30 = (pkin(2) + t14) * t25 - (-qJ(5) - pkin(7) - pkin(6)) * t24 + pkin(1);
t20 = t29 * qJ(2);
t11 = t16 * t24;
t10 = -t25 * t32 - t37;
t9 = t25 * t33 - t36;
t8 = t25 * t36 - t33;
t7 = t25 * t37 + t32;
t6 = -t25 * t34 - t39;
t4 = t25 * t38 - t35;
t2 = -g(2) * t4 - g(3) * t6 + t19 * t41;
t12 = [0, t16, -t15, t16 * t25, -t11, t15, -g(2) * (-t29 * pkin(1) - t27 * qJ(2)) - g(3) * (-t27 * pkin(1) + t20), 0, 0, 0, 0, 0, -g(2) * t10 + g(3) * t8, -g(2) * t9 - g(3) * t7, 0, 0, 0, 0, 0, -g(2) * t6 + g(3) * t4, -g(2) * t5 - g(3) * t3, t11, -g(3) * t20 + (g(2) * t30 - g(3) * t13) * t29 + (-g(2) * (-qJ(2) - t13) + g(3) * t30) * t27; 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t7 + g(3) * t9 + t26 * t41, -g(2) * t8 - g(3) * t10 + t28 * t41, 0, 0, 0, 0, 0, t1, t2, 0, t13 * t41 - g(2) * (t29 * t14 + t27 * t40) - g(3) * (t27 * t14 - t29 * t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t25 + t15 * t24;];
taug_reg = t12;
