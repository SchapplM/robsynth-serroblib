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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
t25 = qJ(3) + qJ(4);
t19 = cos(t25);
t30 = cos(qJ(3));
t14 = t30 * pkin(3) + pkin(4) * t19;
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t48 = -(-qJ(5) - pkin(7) - pkin(6)) * t26 + (pkin(2) + t14) * t27;
t18 = sin(t25);
t31 = cos(qJ(1));
t36 = t31 * t19;
t29 = sin(qJ(1));
t41 = t29 * t18;
t3 = -t27 * t41 - t36;
t45 = g(1) * t26;
t37 = t31 * t18;
t40 = t29 * t19;
t5 = t27 * t37 - t40;
t1 = -g(2) * t3 - g(3) * t5 + t18 * t45;
t28 = sin(qJ(3));
t13 = t28 * pkin(3) + pkin(4) * t18;
t42 = t29 * t13;
t39 = t29 * t28;
t38 = t29 * t30;
t35 = t31 * t28;
t34 = t31 * t30;
t33 = t31 * pkin(1) + t29 * qJ(2);
t16 = g(2) * t31 + g(3) * t29;
t15 = g(2) * t29 - g(3) * t31;
t21 = t29 * pkin(1);
t11 = t16 * t26;
t10 = t27 * t34 + t39;
t9 = t27 * t35 - t38;
t8 = t27 * t38 - t35;
t7 = -t27 * t39 - t34;
t6 = t27 * t36 + t41;
t4 = t27 * t40 - t37;
t2 = g(2) * t4 - g(3) * t6 + t19 * t45;
t12 = [0, -t16, t15, -t16 * t27, t11, -t15, -g(2) * t33 - g(3) * (-t31 * qJ(2) + t21), 0, 0, 0, 0, 0, -g(2) * t10 - g(3) * t8, g(2) * t9 - g(3) * t7, 0, 0, 0, 0, 0, -g(2) * t6 - g(3) * t4, g(2) * t5 - g(3) * t3, -t11, -g(2) * (t33 + t42) - g(3) * (t48 * t29 + t21) + (-g(2) * t48 - g(3) * (-qJ(2) - t13)) * t31; 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t7 - g(3) * t9 + t28 * t45, g(2) * t8 - g(3) * t10 + t30 * t45, 0, 0, 0, 0, 0, t1, t2, 0, t13 * t45 - g(2) * (-t31 * t14 - t27 * t42) - g(3) * (t31 * t27 * t13 - t29 * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t27 - t15 * t26;];
taug_reg = t12;
