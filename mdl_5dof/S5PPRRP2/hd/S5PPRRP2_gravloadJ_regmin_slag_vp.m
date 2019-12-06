% Calculate minimal parameter regressor of gravitation load for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = sin(pkin(7));
t16 = cos(pkin(7));
t22 = g(1) * t16 + g(2) * t15;
t14 = pkin(8) + qJ(3);
t12 = sin(t14);
t27 = g(3) * t12;
t17 = sin(qJ(4));
t26 = t15 * t17;
t18 = cos(qJ(4));
t25 = t15 * t18;
t24 = t16 * t17;
t23 = t16 * t18;
t21 = pkin(4) * t18 + qJ(5) * t17 + pkin(3);
t13 = cos(t14);
t5 = t13 * t26 + t23;
t7 = t13 * t24 - t25;
t1 = g(1) * t7 + g(2) * t5 + t17 * t27;
t6 = t13 * t25 - t24;
t8 = t13 * t23 + t26;
t20 = g(1) * t8 + g(2) * t6 + t18 * t27;
t19 = -g(3) * t13 + t22 * t12;
t9 = -g(1) * t15 + g(2) * t16;
t4 = t22 * t13 + t27;
t3 = t19 * t18;
t2 = t19 * t17;
t10 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, t19, t4, 0, 0, 0, 0, 0, t3, -t2, t3, -t4, t2, (-t22 * pkin(6) - g(3) * t21) * t13 + (-g(3) * pkin(6) + t22 * t21) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t20, t1, 0, -t20, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - (-pkin(4) * t17 + qJ(5) * t18) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
