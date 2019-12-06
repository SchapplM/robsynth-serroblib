% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP2
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
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = qJ(3) + qJ(4);
t16 = cos(t20);
t24 = cos(qJ(3));
t28 = t24 * pkin(3) + pkin(4) * t16;
t21 = qJ(1) + qJ(2);
t15 = sin(t21);
t17 = cos(t21);
t19 = -qJ(5) - pkin(8) - pkin(7);
t9 = pkin(2) + t28;
t27 = t15 * t19 - t17 * t9;
t8 = g(2) * t17 + g(3) * t15;
t7 = g(2) * t15 - g(3) * t17;
t26 = -t15 * t9 - t17 * t19;
t14 = sin(t20);
t2 = -g(1) * t16 - t7 * t14;
t25 = cos(qJ(1));
t23 = sin(qJ(1));
t22 = sin(qJ(3));
t6 = t8 * t24;
t5 = t8 * t22;
t4 = t8 * t16;
t3 = t8 * t14;
t1 = g(1) * t14 - t7 * t16;
t10 = [0, g(2) * t25 + g(3) * t23, -g(2) * t23 + g(3) * t25, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * (-t25 * pkin(1) + t27) - g(3) * (-t23 * pkin(1) + t26); 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * t27 - g(3) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t24 - t7 * t22, g(1) * t22 - t7 * t24, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t28 + t7 * (-t22 * pkin(3) - pkin(4) * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg = t10;
