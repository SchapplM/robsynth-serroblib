% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP7
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
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = qJ(1) + pkin(8);
t13 = sin(t15);
t14 = cos(t15);
t37 = -g(1) * t14 - g(2) * t13;
t17 = sin(qJ(3));
t34 = g(3) * t17;
t33 = t17 * pkin(7);
t16 = sin(qJ(4));
t20 = cos(qJ(3));
t32 = t16 * t20;
t19 = cos(qJ(4));
t31 = t19 * t20;
t6 = t13 * t32 + t14 * t19;
t8 = -t13 * t19 + t14 * t32;
t30 = g(1) * t6 - g(2) * t8;
t28 = g(1) * t13 - g(2) * t14;
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t27 = g(1) * t18 - g(2) * t21;
t26 = t20 * pkin(3) + pkin(2) + t33;
t25 = pkin(4) * t19 + qJ(5) * t16 + pkin(3);
t1 = g(1) * t8 + g(2) * t6 + t16 * t34;
t7 = t13 * t31 - t14 * t16;
t9 = t13 * t16 + t14 * t31;
t24 = g(1) * t9 + g(2) * t7 + t19 * t34;
t23 = -g(3) * t20 - t37 * t17;
t10 = t28 * t17;
t5 = -t20 * t37 + t34;
t4 = t23 * t19;
t3 = t23 * t16;
t2 = g(1) * t7 - g(2) * t9;
t11 = [0, t27, g(1) * t21 + g(2) * t18, t27 * pkin(1), 0, 0, 0, 0, 0, t28 * t20, -t10, 0, 0, 0, 0, 0, t2, -t30, t2, t10, t30, -g(1) * (-t18 * pkin(1) - t7 * pkin(4) - t6 * qJ(5)) - g(2) * (t21 * pkin(1) + t9 * pkin(4) + t8 * qJ(5)) + (-g(1) * pkin(6) - g(2) * t26) * t14 + (-g(2) * pkin(6) + g(1) * t26) * t13; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t5, t3, -g(3) * (t25 * t20 + t33) + t37 * (pkin(7) * t20 - t25 * t17); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t24, t1, 0, -t24, -g(1) * (-t8 * pkin(4) + t9 * qJ(5)) - g(2) * (-t6 * pkin(4) + t7 * qJ(5)) - (-pkin(4) * t16 + qJ(5) * t19) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
