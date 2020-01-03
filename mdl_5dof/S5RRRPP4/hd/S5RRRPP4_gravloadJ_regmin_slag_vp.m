% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = qJ(2) + qJ(3);
t17 = pkin(8) + t23;
t14 = sin(t17);
t15 = cos(t17);
t28 = t15 * pkin(4) + t14 * qJ(5);
t18 = sin(t23);
t35 = pkin(3) * t18;
t34 = pkin(4) * t14;
t19 = cos(t23);
t16 = pkin(3) * t19;
t26 = cos(qJ(2));
t20 = t26 * pkin(2);
t33 = t16 + t20;
t32 = qJ(5) * t15;
t31 = t16 + t28;
t24 = sin(qJ(2));
t7 = -t24 * pkin(2) - t35;
t30 = t7 - t34;
t29 = -t34 - t35;
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t11 = g(1) * t27 + g(2) * t25;
t10 = g(1) * t25 - g(2) * t27;
t3 = -g(3) * t19 + t11 * t18;
t22 = -qJ(4) - pkin(7) - pkin(6);
t9 = t27 * t32;
t8 = t25 * t32;
t6 = pkin(1) + t33;
t5 = t27 * t6;
t4 = g(3) * t18 + t11 * t19;
t2 = -g(3) * t14 - t11 * t15;
t1 = -g(3) * t15 + t11 * t14;
t12 = [0, t10, t11, 0, 0, 0, 0, 0, t10 * t26, -t10 * t24, 0, 0, 0, 0, 0, t10 * t19, -t10 * t18, -t11, -g(1) * (-t27 * t22 - t25 * t6) - g(2) * (-t25 * t22 + t5), t10 * t15, -t11, t10 * t14, -g(2) * t5 + (g(1) * t22 - g(2) * t28) * t27 + (-g(1) * (-t28 - t6) + g(2) * t22) * t25; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 + t11 * t24, g(3) * t24 + t11 * t26, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t33 - t11 * t7, t1, 0, t2, -g(1) * (t30 * t27 + t9) - g(2) * (t30 * t25 + t8) - g(3) * (t20 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t3 * pkin(3), t1, 0, t2, -g(1) * (t29 * t27 + t9) - g(2) * (t29 * t25 + t8) - g(3) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
