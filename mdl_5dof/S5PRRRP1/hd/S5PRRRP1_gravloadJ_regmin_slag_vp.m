% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP1
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
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:12
% EndTime: 2019-12-05 16:40:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (115->23), mult. (70->26), div. (0->0), fcn. (63->6), ass. (0->18)
t12 = pkin(8) + qJ(2);
t13 = -qJ(5) - pkin(7);
t11 = qJ(3) + t12;
t6 = sin(t11);
t7 = cos(t11);
t15 = cos(qJ(4));
t8 = t15 * pkin(4) + pkin(3);
t18 = -t6 * t13 + t7 * t8;
t4 = g(1) * t7 + g(2) * t6;
t3 = g(1) * t6 - g(2) * t7;
t17 = -t7 * t13 - t6 * t8;
t14 = sin(qJ(4));
t16 = -g(3) * t15 + t4 * t14;
t10 = cos(t12);
t9 = sin(t12);
t2 = t3 * t15;
t1 = t3 * t14;
t5 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, g(1) * t9 - g(2) * t10, g(1) * t10 + g(2) * t9, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-pkin(2) * t9 + t17) - g(2) * (pkin(2) * t10 + t18); 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t17 - g(2) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, g(3) * t14 + t4 * t15, 0, t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg = t5;
