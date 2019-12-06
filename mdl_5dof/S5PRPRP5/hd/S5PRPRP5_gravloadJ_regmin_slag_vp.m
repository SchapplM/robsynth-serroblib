% Calculate minimal parameter regressor of gravitation load for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t18 = sin(pkin(7));
t20 = cos(pkin(7));
t28 = g(1) * t20 + g(2) * t18;
t22 = sin(qJ(2));
t23 = cos(qJ(2));
t8 = -g(3) * t23 + t28 * t22;
t32 = g(3) * t22;
t30 = t18 * t23;
t29 = t20 * t23;
t16 = pkin(8) + qJ(4);
t13 = sin(t16);
t14 = cos(t16);
t19 = cos(pkin(8));
t26 = t19 * pkin(3) + pkin(4) * t14 + qJ(5) * t13 + pkin(2);
t4 = t13 * t30 + t20 * t14;
t6 = t13 * t29 - t18 * t14;
t1 = g(1) * t6 + g(2) * t4 + t13 * t32;
t5 = -t20 * t13 + t14 * t30;
t7 = t18 * t13 + t14 * t29;
t25 = g(1) * t7 + g(2) * t5 + t14 * t32;
t21 = -pkin(6) - qJ(3);
t9 = t28 * t23 + t32;
t3 = t8 * t14;
t2 = t8 * t13;
t10 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t8, t9, t8 * t19, -t8 * sin(pkin(8)), -t9, -g(3) * (t23 * pkin(2) + t22 * qJ(3)) + t28 * (pkin(2) * t22 - qJ(3) * t23), 0, 0, 0, 0, 0, t3, -t2, t3, -t9, t2, -g(3) * (-t22 * t21 + t26 * t23) + t28 * (t21 * t23 + t26 * t22); 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t25, t1, 0, -t25, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - (-pkin(4) * t13 + qJ(5) * t14) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
