% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRR6
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:58:07
% EndTime: 2019-05-05 15:58:07
% DurationCPUTime: 0.15s
% Computational Cost: add. (101->41), mult. (174->60), div. (0->0), fcn. (190->8), ass. (0->34)
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t12 = g(1) * t25 + g(2) * t22;
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t27 = -g(3) * t21 + t12 * t24;
t37 = g(3) * t24;
t19 = qJ(5) + qJ(6);
t13 = sin(t19);
t36 = t22 * t13;
t14 = cos(t19);
t35 = t22 * t14;
t20 = sin(qJ(5));
t34 = t22 * t20;
t23 = cos(qJ(5));
t33 = t22 * t23;
t32 = t25 * t13;
t31 = t25 * t14;
t30 = t25 * t20;
t29 = t25 * t23;
t28 = t25 * pkin(1) + t22 * qJ(2);
t11 = g(1) * t22 - g(2) * t25;
t16 = t25 * qJ(2);
t10 = t21 * t29 - t34;
t9 = -t21 * t30 - t33;
t8 = -t21 * t33 - t30;
t7 = t21 * t34 - t29;
t6 = t21 * t31 - t36;
t5 = -t21 * t32 - t35;
t4 = -t21 * t35 - t32;
t3 = t21 * t36 - t31;
t2 = g(1) * t6 - g(2) * t4 + t14 * t37;
t1 = -g(1) * t5 + g(2) * t3 + t13 * t37;
t15 = [0, t11, t12, -t11, -t12, -g(1) * (-t22 * pkin(1) + t16) - g(2) * t28, -t12, t11, -g(1) * (t16 + (-pkin(1) - qJ(3)) * t22) - g(2) * (t25 * qJ(3) + t28) 0, 0, 0, 0, 0, t11 * t21, t11 * t24, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, -t11, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t12 * t21 + t37, 0, 0, 0, 0, 0, -t27 * t23, t27 * t20, 0, 0, 0, 0, 0, -t27 * t14, t27 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t20 * t37, g(1) * t10 - g(2) * t8 + t23 * t37, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t15;
