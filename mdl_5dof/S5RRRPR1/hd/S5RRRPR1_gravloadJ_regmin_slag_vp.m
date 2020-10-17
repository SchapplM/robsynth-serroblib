% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:51
% EndTime: 2019-12-05 18:38:51
% DurationCPUTime: 0.09s
% Computational Cost: add. (142->23), mult. (121->35), div. (0->0), fcn. (113->8), ass. (0->20)
t17 = qJ(2) + qJ(3);
t14 = cos(t17);
t20 = cos(qJ(2));
t22 = t20 * pkin(2) + pkin(3) * t14;
t19 = sin(qJ(1));
t21 = cos(qJ(1));
t8 = g(1) * t21 + g(2) * t19;
t7 = g(1) * t19 - g(2) * t21;
t13 = sin(t17);
t3 = -g(3) * t14 + t8 * t13;
t18 = sin(qJ(2));
t16 = -qJ(4) - pkin(7) - pkin(6);
t12 = pkin(9) + qJ(5) + t17;
t10 = cos(t12);
t9 = sin(t12);
t5 = pkin(1) + t22;
t4 = g(3) * t13 + t8 * t14;
t2 = g(3) * t9 + t8 * t10;
t1 = -g(3) * t10 + t8 * t9;
t6 = [0, t7, t8, 0, 0, 0, 0, 0, t7 * t20, -t7 * t18, 0, 0, 0, 0, 0, t7 * t14, -t7 * t13, -t8, -g(1) * (-t21 * t16 - t19 * t5) - g(2) * (-t19 * t16 + t21 * t5), 0, 0, 0, 0, 0, t7 * t10, -t7 * t9; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t20 + t8 * t18, g(3) * t18 + t8 * t20, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t22 - t8 * (-t18 * pkin(2) - pkin(3) * t13), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t3 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t6;
