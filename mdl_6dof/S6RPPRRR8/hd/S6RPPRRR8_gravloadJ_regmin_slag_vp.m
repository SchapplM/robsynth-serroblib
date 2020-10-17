% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:14:52
% EndTime: 2019-05-05 16:14:52
% DurationCPUTime: 0.17s
% Computational Cost: add. (150->45), mult. (180->62), div. (0->0), fcn. (196->10), ass. (0->35)
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t11 = g(1) * t26 - g(2) * t28;
t21 = pkin(10) + qJ(4);
t13 = sin(t21);
t14 = cos(t21);
t30 = -g(3) * t13 + t11 * t14;
t40 = g(3) * t14;
t22 = qJ(5) + qJ(6);
t15 = sin(t22);
t39 = t26 * t15;
t16 = cos(t22);
t38 = t26 * t16;
t25 = sin(qJ(5));
t37 = t26 * t25;
t27 = cos(qJ(5));
t36 = t26 * t27;
t35 = t28 * t15;
t34 = t28 * t16;
t33 = t28 * t25;
t32 = t28 * t27;
t31 = t28 * pkin(1) + t26 * qJ(2);
t12 = g(1) * t28 + g(2) * t26;
t18 = t28 * qJ(2);
t10 = t13 * t32 - t37;
t9 = t13 * t33 + t36;
t8 = t13 * t36 + t33;
t7 = -t13 * t37 + t32;
t6 = t13 * t34 - t39;
t5 = t13 * t35 + t38;
t4 = t13 * t38 + t35;
t3 = -t13 * t39 + t34;
t2 = g(1) * t4 - g(2) * t6 + t16 * t40;
t1 = -g(1) * t3 - g(2) * t5 + t15 * t40;
t17 = [0, t11, t12, -t11, -t12, -g(1) * (-t26 * pkin(1) + t18) - g(2) * t31, -t12 * sin(pkin(10)) -t12 * cos(pkin(10)) t11, -g(1) * (t18 + (-pkin(1) - qJ(3)) * t26) - g(2) * (t28 * qJ(3) + t31) 0, 0, 0, 0, 0, -t12 * t13, -t12 * t14, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, -t11, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t11 * t13 + t40, 0, 0, 0, 0, 0, -t30 * t27, t30 * t25, 0, 0, 0, 0, 0, -t30 * t16, t30 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t25 * t40, g(1) * t8 - g(2) * t10 + t27 * t40, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t17;
