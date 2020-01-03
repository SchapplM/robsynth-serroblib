% Calculate minimal parameter regressor of gravitation load for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t19 = sin(qJ(2));
t20 = sin(qJ(1));
t23 = cos(qJ(2));
t24 = cos(qJ(1));
t30 = cos(pkin(4));
t28 = t24 * t30;
t10 = t20 * t19 - t23 * t28;
t17 = sin(qJ(4));
t21 = cos(qJ(4));
t11 = t19 * t28 + t20 * t23;
t18 = sin(qJ(3));
t22 = cos(qJ(3));
t16 = sin(pkin(4));
t34 = t16 * t24;
t4 = t11 * t22 - t18 * t34;
t42 = -t10 * t21 + t4 * t17;
t41 = t10 * t17 + t4 * t21;
t40 = g(3) * t16;
t37 = t16 * t19;
t36 = t16 * t22;
t35 = t16 * t23;
t33 = t17 * t22;
t32 = t21 * t22;
t31 = t22 * t23;
t29 = t20 * t30;
t27 = t11 * t18 + t22 * t34;
t13 = -t19 * t29 + t24 * t23;
t6 = -t13 * t18 + t20 * t36;
t26 = g(1) * t6 - g(2) * t27 + g(3) * (-t18 * t37 + t30 * t22);
t12 = t24 * t19 + t23 * t29;
t25 = -g(1) * t12 - g(2) * t10 + g(3) * t35;
t9 = t30 * t18 + t19 * t36;
t7 = t20 * t16 * t18 + t13 * t22;
t2 = t12 * t17 + t7 * t21;
t1 = t12 * t21 - t7 * t17;
t3 = [0, g(1) * t20 - g(2) * t24, g(1) * t24 + g(2) * t20, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t13, -g(1) * t10 + g(2) * t12, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t27 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t41 - g(2) * t2, -g(1) * t42 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t25, g(1) * t13 + g(2) * t11 + g(3) * t37, 0, 0, 0, 0, 0, -t25 * t22, t25 * t18, 0, 0, 0, 0, 0, -g(1) * (-t12 * t32 + t13 * t17) - g(2) * (-t10 * t32 + t11 * t17) - (t17 * t19 + t21 * t31) * t40, -g(1) * (t12 * t33 + t13 * t21) - g(2) * (t10 * t33 + t11 * t21) - (-t17 * t31 + t19 * t21) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, g(1) * t7 + g(2) * t4 + g(3) * t9, 0, 0, 0, 0, 0, -t26 * t21, t26 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t42 - g(3) * (-t9 * t17 - t21 * t35), g(1) * t2 + g(2) * t41 - g(3) * (t17 * t35 - t9 * t21);];
taug_reg = t3;
