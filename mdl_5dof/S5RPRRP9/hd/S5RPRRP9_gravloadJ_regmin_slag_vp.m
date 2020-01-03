% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP9
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
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = pkin(8) + qJ(3);
t16 = cos(t20);
t17 = qJ(4) + t20;
t13 = sin(t17);
t14 = cos(t17);
t28 = t14 * pkin(4) + t13 * qJ(5);
t30 = pkin(3) * t16 + t28;
t29 = pkin(4) * t13;
t27 = qJ(5) * t14;
t15 = sin(t20);
t26 = -pkin(3) * t15 - t29;
t23 = sin(qJ(1));
t24 = cos(qJ(1));
t9 = g(1) * t24 + g(2) * t23;
t8 = g(1) * t23 - g(2) * t24;
t22 = cos(pkin(8));
t25 = t22 * pkin(2) + pkin(1) + t30;
t19 = -pkin(7) - pkin(6) - qJ(2);
t7 = t24 * t27;
t6 = t23 * t27;
t4 = t8 * t14;
t3 = t8 * t13;
t2 = g(3) * t13 + t9 * t14;
t1 = -g(3) * t14 + t9 * t13;
t5 = [0, t8, t9, t8 * t22, -t8 * sin(pkin(8)), -t9, -g(1) * (-t23 * pkin(1) + t24 * qJ(2)) - g(2) * (t24 * pkin(1) + t23 * qJ(2)), 0, 0, 0, 0, 0, t8 * t16, -t8 * t15, 0, 0, 0, 0, 0, t4, -t3, t4, -t9, t3, (g(1) * t19 - g(2) * t25) * t24 + (g(1) * t25 + g(2) * t19) * t23; 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t16 + t9 * t15, g(3) * t15 + t9 * t16, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (t26 * t24 + t7) - g(2) * (t26 * t23 + t6) - g(3) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t24 * t29 + t7) - g(2) * (-t23 * t29 + t6) - g(3) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t5;
