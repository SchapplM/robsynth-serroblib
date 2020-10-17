% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:51
% EndTime: 2019-12-05 17:04:52
% DurationCPUTime: 0.14s
% Computational Cost: add. (140->31), mult. (96->33), div. (0->0), fcn. (84->8), ass. (0->24)
t16 = qJ(2) + qJ(3);
t12 = sin(t16);
t24 = pkin(3) * t12;
t20 = cos(qJ(2));
t13 = cos(t16);
t9 = pkin(3) * t13;
t23 = t20 * pkin(2) + t9;
t18 = sin(qJ(2));
t22 = -t18 * pkin(2) - t24;
t14 = qJ(4) + t16;
t10 = sin(t14);
t11 = cos(t14);
t4 = g(1) * t11 + g(2) * t10;
t3 = g(1) * t10 - g(2) * t11;
t5 = g(1) * t12 - g(2) * t13;
t21 = g(1) * t18 - g(2) * t20;
t19 = cos(qJ(5));
t17 = sin(qJ(5));
t8 = t11 * pkin(6);
t7 = t10 * pkin(6);
t6 = g(1) * t13 + g(2) * t12;
t2 = t3 * t19;
t1 = t3 * t17;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, g(1) * t20 + g(2) * t18, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t21 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t22 - g(2) * t23, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t22 + t8) - g(2) * (t7 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t8 - t24) - g(2) * (t7 + t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -t4 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t19 + t4 * t17, g(3) * t17 + t4 * t19, 0, 0;];
taug_reg = t15;
