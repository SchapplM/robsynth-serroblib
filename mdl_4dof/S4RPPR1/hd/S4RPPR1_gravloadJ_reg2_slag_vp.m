% Calculate inertial parameters regressor of gravitation load for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:08:32
% EndTime: 2019-05-04 19:08:32
% DurationCPUTime: 0.09s
% Computational Cost: add. (74->26), mult. (62->29), div. (0->0), fcn. (64->6), ass. (0->17)
t24 = cos(qJ(4));
t23 = sin(qJ(4));
t15 = qJ(1) + pkin(6);
t12 = sin(t15);
t13 = cos(t15);
t17 = cos(qJ(1));
t22 = t17 * pkin(1) + t13 * pkin(2) + t12 * qJ(3);
t16 = sin(qJ(1));
t21 = -t16 * pkin(1) + t13 * qJ(3);
t1 = -t12 * t23 - t13 * t24;
t2 = -t12 * t24 + t13 * t23;
t20 = g(1) * t2 - g(2) * t1;
t19 = g(1) * t1 + g(2) * t2;
t18 = g(1) * t16 - g(2) * t17;
t4 = g(1) * t13 + g(2) * t12;
t3 = g(1) * t12 - g(2) * t13;
t5 = [0, 0, 0, 0, 0, 0, t18, g(1) * t17 + g(2) * t16, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t18 * pkin(1), 0, 0, 0, 0, 0, 0, t3, 0, -t4, -g(1) * (-t12 * pkin(2) + t21) - g(2) * t22, 0, 0, 0, 0, 0, 0, -t20, t19, 0, -g(1) * ((-pkin(2) - pkin(3)) * t12 + t21) - g(2) * (t13 * pkin(3) + t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, 0;];
taug_reg  = t5;
