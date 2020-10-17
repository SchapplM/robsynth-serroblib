% Calculate minimal parameter regressor of gravitation load for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:36:14
% EndTime: 2019-05-05 13:36:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (130->35), mult. (116->49), div. (0->0), fcn. (118->10), ass. (0->26)
t16 = pkin(10) + qJ(5);
t11 = sin(t16);
t13 = cos(t16);
t17 = qJ(1) + pkin(9);
t12 = sin(t17);
t14 = cos(t17);
t5 = g(1) * t12 - g(2) * t14;
t34 = -g(3) * t11 + t5 * t13;
t32 = g(3) * t13;
t20 = sin(qJ(6));
t31 = t12 * t20;
t22 = cos(qJ(6));
t30 = t12 * t22;
t29 = t14 * t20;
t28 = t14 * t22;
t23 = cos(qJ(1));
t27 = t23 * pkin(1) + t14 * pkin(2) + t12 * qJ(3);
t21 = sin(qJ(1));
t26 = -t21 * pkin(1) + t14 * qJ(3);
t6 = -g(1) * t14 - g(2) * t12;
t25 = g(1) * t21 - g(2) * t23;
t4 = t11 * t28 - t31;
t3 = t11 * t29 + t30;
t2 = t11 * t30 + t29;
t1 = -t11 * t31 + t28;
t7 = [0, t25, g(1) * t23 + g(2) * t21, t25 * pkin(1), -t5, t6, -g(1) * (-t12 * pkin(2) + t26) - g(2) * t27, t6 * sin(pkin(10)) t6 * cos(pkin(10)) t5, -g(1) * ((-pkin(2) - qJ(4)) * t12 + t26) - g(2) * (t14 * qJ(4) + t27) 0, 0, 0, 0, 0, t6 * t11, t6 * t13, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t5 * t11 + t32, 0, 0, 0, 0, 0, -t34 * t22, t34 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t20 * t32, g(1) * t2 - g(2) * t4 + t22 * t32;];
taug_reg  = t7;
