% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR12
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
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = sin(qJ(1));
t24 = cos(qJ(1));
t9 = g(1) * t24 + g(2) * t23;
t17 = pkin(8) + qJ(3);
t12 = sin(t17);
t14 = cos(t17);
t1 = -g(3) * t14 + t9 * t12;
t38 = g(3) * t12;
t16 = pkin(9) + qJ(5);
t11 = sin(t16);
t36 = t23 * t11;
t13 = cos(t16);
t35 = t23 * t13;
t18 = sin(pkin(9));
t34 = t23 * t18;
t20 = cos(pkin(9));
t33 = t23 * t20;
t32 = t24 * t11;
t31 = t24 * t13;
t30 = t24 * t18;
t29 = t24 * t20;
t8 = g(1) * t23 - g(2) * t24;
t28 = t14 * pkin(3) + t12 * qJ(4);
t21 = cos(pkin(8));
t26 = t21 * pkin(2) + pkin(1) + t28;
t22 = -pkin(6) - qJ(2);
t7 = t8 * t12;
t6 = t14 * t31 + t36;
t5 = -t14 * t32 + t35;
t4 = -t14 * t35 + t32;
t3 = t14 * t36 + t31;
t2 = t9 * t14 + t38;
t10 = [0, t8, t9, t8 * t21, -t8 * sin(pkin(8)), -t9, -g(1) * (-t23 * pkin(1) + t24 * qJ(2)) - g(2) * (t24 * pkin(1) + t23 * qJ(2)), 0, 0, 0, 0, 0, t8 * t14, -t7, -g(1) * (-t14 * t33 + t30) - g(2) * (t14 * t29 + t34), -g(1) * (t14 * t34 + t29) - g(2) * (-t14 * t30 + t33), t7, (g(1) * t22 - g(2) * t26) * t24 + (g(1) * t26 + g(2) * t22) * t23, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1 * t20, -t1 * t18, -t2, -g(3) * t28 + t9 * (pkin(3) * t12 - qJ(4) * t14), 0, 0, 0, 0, 0, t1 * t13, -t1 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t11 * t38, g(1) * t6 - g(2) * t4 + t13 * t38;];
taug_reg = t10;
