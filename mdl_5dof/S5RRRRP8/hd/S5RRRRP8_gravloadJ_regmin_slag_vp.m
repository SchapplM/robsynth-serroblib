% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP8
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
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t34 = g(1) * t28 + g(2) * t25;
t22 = qJ(3) + qJ(4);
t18 = sin(t22);
t19 = cos(t22);
t38 = t28 * t19;
t27 = cos(qJ(2));
t40 = t25 * t27;
t3 = t18 * t40 + t38;
t24 = sin(qJ(2));
t43 = g(3) * t24;
t39 = t28 * t18;
t5 = t25 * t19 - t27 * t39;
t1 = -g(1) * t5 + g(2) * t3 + t18 * t43;
t7 = -g(3) * t27 + t34 * t24;
t23 = sin(qJ(3));
t15 = t23 * pkin(3) + pkin(4) * t18;
t41 = pkin(6) + t15;
t37 = t28 * t23;
t26 = cos(qJ(3));
t36 = t28 * t26;
t16 = t26 * pkin(3) + pkin(4) * t19;
t33 = g(1) * t25 - g(2) * t28;
t14 = pkin(2) + t16;
t21 = -qJ(5) - pkin(8) - pkin(7);
t32 = t27 * t14 - t24 * t21;
t30 = pkin(1) + t32;
t13 = t33 * t24;
t12 = t25 * t23 + t27 * t36;
t11 = t25 * t26 - t27 * t37;
t10 = -t26 * t40 + t37;
t9 = t23 * t40 + t36;
t8 = t34 * t27 + t43;
t6 = t25 * t18 + t27 * t38;
t4 = -t19 * t40 + t39;
t2 = g(1) * t6 - g(2) * t4 + t19 * t43;
t17 = [0, t33, t34, 0, 0, 0, 0, 0, t33 * t27, -t13, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t13, (-g(1) * t41 - g(2) * t30) * t28 + (g(1) * t30 - g(2) * t41) * t25; 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t7 * t26, -t7 * t23, 0, 0, 0, 0, 0, t7 * t19, -t7 * t18, -t8, -g(3) * t32 + t34 * (t14 * t24 + t21 * t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t23 * t43, g(1) * t12 - g(2) * t10 + t26 * t43, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t28 * t27 * t15 + t25 * t16) - g(2) * (-t15 * t40 - t28 * t16) + t15 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t17;
