% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:56
% EndTime: 2019-12-05 17:39:57
% DurationCPUTime: 0.14s
% Computational Cost: add. (111->45), mult. (130->43), div. (0->0), fcn. (119->8), ass. (0->21)
t17 = sin(pkin(8));
t24 = t17 * pkin(3);
t19 = -pkin(6) - qJ(3);
t20 = sin(qJ(1));
t21 = cos(qJ(1));
t23 = t21 * pkin(1) + t20 * qJ(2);
t16 = pkin(8) + qJ(4);
t5 = g(1) * t21 + g(2) * t20;
t4 = g(1) * t20 - g(2) * t21;
t8 = sin(t16);
t9 = cos(t16);
t22 = g(3) * t8 - t4 * t9;
t15 = -pkin(7) + t19;
t12 = t21 * qJ(2);
t10 = qJ(5) + t16;
t7 = cos(t10);
t6 = sin(t10);
t3 = pkin(4) * t8 + t24;
t2 = g(3) * t6 - t4 * t7;
t1 = g(3) * t7 + t4 * t6;
t11 = [0, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5, -g(1) * (-t20 * pkin(1) + t12) - g(2) * t23, 0, 0, 0, 0, 0, 0, -t5 * t17, -t5 * cos(pkin(8)), t4, -g(1) * (t12 + (-pkin(1) - qJ(3)) * t20) - g(2) * (t21 * qJ(3) + t23), 0, 0, 0, 0, 0, 0, -t5 * t8, -t5 * t9, t4, -g(1) * (t21 * t24 + t12 + (-pkin(1) + t19) * t20) - g(2) * (-t21 * t19 + t20 * t24 + t23), 0, 0, 0, 0, 0, 0, -t5 * t6, -t5 * t7, t4, -g(1) * (t21 * t3 + t12 + (-pkin(1) + t15) * t20) - g(2) * (-t21 * t15 + t20 * t3 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t9 + t4 * t8, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t22 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t11;
