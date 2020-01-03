% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t38 = g(1) * t31 + g(2) * t28;
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t42 = t31 * t29;
t30 = cos(qJ(2));
t46 = t28 * t30;
t14 = t26 * t46 + t42;
t43 = t31 * t26;
t16 = t28 * t29 - t30 * t43;
t27 = sin(qJ(2));
t50 = g(3) * t27;
t55 = -g(1) * t16 + g(2) * t14 + t26 * t50;
t34 = -g(3) * t30 + t38 * t27;
t25 = qJ(3) + qJ(4);
t23 = sin(t25);
t48 = t23 * t27;
t24 = cos(t25);
t47 = t24 * t27;
t45 = t31 * t23;
t44 = t31 * t24;
t40 = pkin(3) * t26 + pkin(6);
t11 = -t28 * t24 + t30 * t45;
t9 = t23 * t46 + t44;
t39 = g(1) * t9 - g(2) * t11;
t37 = g(1) * t28 - g(2) * t31;
t22 = t29 * pkin(3) + pkin(2);
t32 = -pkin(8) - pkin(7);
t36 = t30 * t22 - t27 * t32 + pkin(1);
t35 = pkin(4) * t24 + qJ(5) * t23 + t22;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t48;
t10 = t24 * t46 - t45;
t12 = t28 * t23 + t30 * t44;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t47;
t33 = -g(1) * (-t11 * pkin(4) + t12 * qJ(5)) - g(2) * (-t9 * pkin(4) + t10 * qJ(5)) - g(3) * (-pkin(4) * t48 + qJ(5) * t47);
t18 = t37 * t27;
t17 = t28 * t26 + t30 * t42;
t15 = -t29 * t46 + t43;
t13 = t38 * t30 + t50;
t6 = t34 * t24;
t5 = t34 * t23;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, t37, t38, 0, 0, 0, 0, 0, t37 * t30, -t18, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, 0, 0, 0, 0, 0, t4, -t39, t4, t18, t39, -g(1) * (-t10 * pkin(4) - t9 * qJ(5)) - g(2) * (t12 * pkin(4) + t11 * qJ(5)) + (-g(1) * t40 - g(2) * t36) * t31 + (g(1) * t36 - g(2) * t40) * t28; 0, 0, 0, 0, 0, 0, 0, 0, t34, t13, 0, 0, 0, 0, 0, t34 * t29, -t34 * t26, 0, 0, 0, 0, 0, t6, -t5, t6, -t13, t5, (-g(3) * t35 + t38 * t32) * t30 + (g(3) * t32 + t38 * t35) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, g(1) * t17 - g(2) * t15 + t29 * t50, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t55 * pkin(3) + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t2;
