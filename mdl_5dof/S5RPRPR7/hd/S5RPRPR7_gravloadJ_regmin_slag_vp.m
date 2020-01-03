% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t11 = qJ(1) + pkin(8);
t7 = sin(t11);
t9 = cos(t11);
t22 = g(1) * t9 + g(2) * t7;
t10 = qJ(3) + pkin(9);
t6 = sin(t10);
t8 = cos(t10);
t29 = -g(3) * t8 + t22 * t6;
t28 = g(3) * t6;
t13 = sin(qJ(5));
t26 = t7 * t13;
t16 = cos(qJ(5));
t25 = t7 * t16;
t24 = t9 * t13;
t23 = t9 * t16;
t21 = g(1) * t7 - g(2) * t9;
t15 = sin(qJ(1));
t18 = cos(qJ(1));
t20 = g(1) * t15 - g(2) * t18;
t14 = sin(qJ(3));
t17 = cos(qJ(3));
t19 = -g(3) * t17 + t22 * t14;
t12 = -qJ(4) - pkin(6);
t5 = t17 * pkin(3) + pkin(2);
t4 = t8 * t23 + t26;
t3 = -t8 * t24 + t25;
t2 = -t8 * t25 + t24;
t1 = t8 * t26 + t23;
t27 = [0, t20, g(1) * t18 + g(2) * t15, t20 * pkin(1), 0, 0, 0, 0, 0, t21 * t17, -t21 * t14, -t22, -g(1) * (-t15 * pkin(1) - t9 * t12 - t7 * t5) - g(2) * (t18 * pkin(1) - t7 * t12 + t9 * t5), 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(3) * t14 + t22 * t17, 0, t19 * pkin(3), 0, 0, 0, 0, 0, t29 * t16, -t29 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t13 * t28, g(1) * t4 - g(2) * t2 + t16 * t28;];
taug_reg = t27;
