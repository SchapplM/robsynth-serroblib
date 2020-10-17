% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:50:36
% EndTime: 2019-05-05 15:50:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (93->34), mult. (146->46), div. (0->0), fcn. (150->8), ass. (0->27)
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t10 = g(1) * t23 + g(2) * t20;
t17 = qJ(4) + qJ(5);
t11 = sin(t17);
t12 = cos(t17);
t31 = -g(3) * t11 + t10 * t12;
t29 = g(3) * t12;
t18 = sin(qJ(6));
t28 = t20 * t18;
t21 = cos(qJ(6));
t27 = t20 * t21;
t26 = t23 * t18;
t25 = t23 * t21;
t24 = t23 * pkin(1) + t20 * qJ(2);
t9 = g(1) * t20 - g(2) * t23;
t22 = cos(qJ(4));
t19 = sin(qJ(4));
t14 = t23 * qJ(2);
t8 = t11 * t25 - t28;
t7 = -t11 * t26 - t27;
t6 = -t11 * t27 - t26;
t5 = t11 * t28 - t25;
t3 = t10 * t11 + t29;
t2 = t31 * t21;
t1 = t31 * t18;
t4 = [0, t9, t10, -t9, -t10, -g(1) * (-t20 * pkin(1) + t14) - g(2) * t24, -t10, t9, -g(1) * (t14 + (-pkin(1) - qJ(3)) * t20) - g(2) * (t23 * qJ(3) + t24) 0, 0, 0, 0, 0, t9 * t19, t9 * t22, 0, 0, 0, 0, 0, t9 * t11, t9 * t12, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, -t9, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t19 - t10 * t22, g(3) * t22 + t10 * t19, 0, 0, 0, 0, 0, -t31, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t18 * t29, g(1) * t8 - g(2) * t6 + t21 * t29;];
taug_reg  = t4;
