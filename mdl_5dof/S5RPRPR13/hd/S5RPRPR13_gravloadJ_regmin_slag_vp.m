% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t10 = g(1) * t23 + g(2) * t21;
t16 = pkin(8) + qJ(3);
t13 = sin(t16);
t14 = cos(t16);
t2 = g(3) * t13 + t10 * t14;
t32 = g(3) * t14;
t20 = sin(qJ(5));
t31 = t21 * t20;
t22 = cos(qJ(5));
t30 = t21 * t22;
t29 = t23 * t20;
t28 = t23 * t22;
t9 = g(1) * t21 - g(2) * t23;
t27 = t14 * pkin(3) + t13 * qJ(4);
t18 = cos(pkin(8));
t25 = t18 * pkin(2) + pkin(1) + t27;
t19 = -pkin(6) - qJ(2);
t8 = -t13 * t31 + t28;
t7 = t13 * t30 + t29;
t6 = t13 * t29 + t30;
t5 = t13 * t28 - t31;
t4 = t9 * t14;
t3 = t9 * t13;
t1 = t10 * t13 - t32;
t11 = [0, t9, t10, t9 * t18, -t9 * sin(pkin(8)), -t10, -g(1) * (-t21 * pkin(1) + t23 * qJ(2)) - g(2) * (t23 * pkin(1) + t21 * qJ(2)), 0, 0, 0, 0, 0, t4, -t3, -t10, -t4, t3, (g(1) * t19 - g(2) * t25) * t23 + (g(1) * t25 + g(2) * t19) * t21, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5; 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1, -t2, -g(3) * t27 + t10 * (pkin(3) * t13 - qJ(4) * t14), 0, 0, 0, 0, 0, -t2 * t20, -t2 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t22 * t32, g(1) * t6 - g(2) * t8 - t20 * t32;];
taug_reg = t11;
