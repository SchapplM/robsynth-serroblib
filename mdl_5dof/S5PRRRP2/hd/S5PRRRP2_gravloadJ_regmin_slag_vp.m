% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t19 = sin(qJ(4));
t20 = cos(qJ(4));
t23 = t20 * pkin(4) + t19 * qJ(5);
t31 = -pkin(3) - t23;
t18 = pkin(8) + qJ(2);
t17 = qJ(3) + t18;
t14 = cos(t17);
t13 = sin(t17);
t30 = g(1) * t13;
t5 = -g(2) * t14 + t30;
t6 = g(1) * t14 + g(2) * t13;
t24 = t13 * pkin(7) - t31 * t14;
t21 = t31 * t30;
t16 = cos(t18);
t15 = sin(t18);
t11 = t14 * pkin(7);
t4 = t5 * t20;
t3 = t5 * t19;
t2 = g(3) * t19 + t6 * t20;
t1 = -g(3) * t20 + t6 * t19;
t7 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, g(1) * t15 - g(2) * t16, g(1) * t16 + g(2) * t15, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (-pkin(2) * t15 + t11) - g(2) * (pkin(2) * t16 + t24) - t21; 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t11 - g(2) * t24 - t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t23 + t6 * (pkin(4) * t19 - qJ(5) * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t7;
