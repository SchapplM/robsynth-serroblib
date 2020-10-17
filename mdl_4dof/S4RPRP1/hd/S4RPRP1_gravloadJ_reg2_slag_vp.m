% Calculate inertial parameters regressor of gravitation load for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:14:05
% EndTime: 2019-05-04 19:14:06
% DurationCPUTime: 0.08s
% Computational Cost: add. (90->23), mult. (52->27), div. (0->0), fcn. (44->6), ass. (0->16)
t13 = qJ(1) + pkin(6);
t11 = qJ(3) + t13;
t7 = sin(t11);
t8 = cos(t11);
t20 = t8 * pkin(3) + t7 * qJ(4);
t10 = cos(t13);
t15 = cos(qJ(1));
t19 = t15 * pkin(1) + pkin(2) * t10;
t18 = -t7 * pkin(3) + t8 * qJ(4);
t14 = sin(qJ(1));
t9 = sin(t13);
t17 = -t14 * pkin(1) - pkin(2) * t9;
t16 = g(1) * t14 - g(2) * t15;
t2 = g(1) * t8 + g(2) * t7;
t1 = g(1) * t7 - g(2) * t8;
t3 = [0, 0, 0, 0, 0, 0, t16, g(1) * t15 + g(2) * t14, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t9 - g(2) * t10, g(1) * t10 + g(2) * t9, 0, t16 * pkin(1), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * t17 - g(2) * t19, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t17 + t18) - g(2) * (t19 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * t18 - g(2) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t3;
