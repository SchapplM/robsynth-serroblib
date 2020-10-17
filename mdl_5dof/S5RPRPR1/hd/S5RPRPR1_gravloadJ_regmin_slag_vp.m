% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:43
% EndTime: 2019-12-05 17:47:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (70->26), mult. (88->30), div. (0->0), fcn. (81->6), ass. (0->17)
t13 = sin(qJ(3));
t19 = pkin(3) * t13;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t18 = t16 * pkin(1) + t14 * qJ(2);
t4 = g(1) * t16 + g(2) * t14;
t3 = g(1) * t14 - g(2) * t16;
t15 = cos(qJ(3));
t17 = g(3) * t13 - t3 * t15;
t12 = -qJ(4) - pkin(6);
t9 = t16 * qJ(2);
t7 = qJ(3) + pkin(8) + qJ(5);
t6 = cos(t7);
t5 = sin(t7);
t2 = g(3) * t5 - t3 * t6;
t1 = g(3) * t6 + t3 * t5;
t8 = [0, t3, t4, -t3, -t4, -g(1) * (-t14 * pkin(1) + t9) - g(2) * t18, 0, 0, 0, 0, 0, -t4 * t13, -t4 * t15, t3, -g(1) * (t16 * t19 + t9 + (-pkin(1) + t12) * t14) - g(2) * (-t16 * t12 + t14 * t19 + t18), 0, 0, 0, 0, 0, -t4 * t5, -t4 * t6; 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, g(3) * t15 + t3 * t13, 0, t17 * pkin(3), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t8;
