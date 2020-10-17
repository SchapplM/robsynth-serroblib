% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP1
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
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:00
% EndTime: 2019-12-05 18:46:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (154->26), mult. (137->39), div. (0->0), fcn. (126->8), ass. (0->22)
t19 = qJ(2) + qJ(3);
t17 = qJ(4) + t19;
t13 = cos(t17);
t15 = cos(t19);
t26 = pkin(3) * t15 + pkin(4) * t13;
t22 = cos(qJ(2));
t25 = t22 * pkin(2) + t26;
t12 = sin(t17);
t14 = sin(t19);
t24 = -pkin(3) * t14 - pkin(4) * t12;
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t9 = g(1) * t23 + g(2) * t21;
t8 = g(1) * t21 - g(2) * t23;
t1 = -g(3) * t13 + t9 * t12;
t20 = sin(qJ(2));
t16 = -qJ(5) - pkin(8) - pkin(7) - pkin(6);
t5 = pkin(1) + t25;
t4 = g(3) * t14 + t9 * t15;
t3 = -g(3) * t15 + t9 * t14;
t2 = g(3) * t12 + t9 * t13;
t6 = [0, t8, t9, 0, 0, 0, 0, 0, t8 * t22, -t8 * t20, 0, 0, 0, 0, 0, t8 * t15, -t8 * t14, 0, 0, 0, 0, 0, t8 * t13, -t8 * t12, -t9, -g(1) * (-t23 * t16 - t21 * t5) - g(2) * (-t21 * t16 + t23 * t5); 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t22 + t9 * t20, g(3) * t20 + t9 * t22, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t25 - t9 * (-t20 * pkin(2) + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t26 - t9 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg = t6;
