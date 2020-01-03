% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t12 = g(1) * t27 + g(2) * t24;
t21 = qJ(2) + qJ(3);
t16 = pkin(9) + t21;
t13 = sin(t16);
t14 = cos(t16);
t35 = -g(3) * t14 + t12 * t13;
t34 = g(3) * t13;
t22 = sin(qJ(5));
t32 = t24 * t22;
t25 = cos(qJ(5));
t31 = t24 * t25;
t30 = t27 * t22;
t29 = t27 * t25;
t18 = cos(t21);
t26 = cos(qJ(2));
t28 = t26 * pkin(2) + pkin(3) * t18;
t11 = g(1) * t24 - g(2) * t27;
t17 = sin(t21);
t3 = -g(3) * t18 + t12 * t17;
t23 = sin(qJ(2));
t20 = -qJ(4) - pkin(7) - pkin(6);
t9 = pkin(1) + t28;
t8 = t14 * t29 + t32;
t7 = -t14 * t30 + t31;
t6 = -t14 * t31 + t30;
t5 = t14 * t32 + t29;
t4 = g(3) * t17 + t12 * t18;
t2 = t35 * t25;
t1 = t35 * t22;
t10 = [0, t11, t12, 0, 0, 0, 0, 0, t11 * t26, -t11 * t23, 0, 0, 0, 0, 0, t11 * t18, -t11 * t17, -t12, -g(1) * (-t27 * t20 - t24 * t9) - g(2) * (-t24 * t20 + t27 * t9), 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 + t12 * t23, g(3) * t23 + t12 * t26, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t28 - t12 * (-t23 * pkin(2) - pkin(3) * t17), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t3 * pkin(3), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t22 * t34, g(1) * t8 - g(2) * t6 + t25 * t34;];
taug_reg = t10;
