% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:31:48
% EndTime: 2019-05-05 15:31:48
% DurationCPUTime: 0.16s
% Computational Cost: add. (163->40), mult. (164->67), div. (0->0), fcn. (180->10), ass. (0->32)
t19 = sin(qJ(4));
t22 = cos(qJ(4));
t16 = qJ(1) + pkin(10);
t12 = sin(t16);
t13 = cos(t16);
t27 = g(1) * t12 - g(2) * t13;
t25 = -g(3) * t19 + t27 * t22;
t33 = g(3) * t22;
t17 = qJ(5) + qJ(6);
t14 = sin(t17);
t32 = t14 * t19;
t15 = cos(t17);
t31 = t15 * t19;
t18 = sin(qJ(5));
t30 = t18 * t19;
t21 = cos(qJ(5));
t29 = t19 * t21;
t28 = -g(1) * t13 - g(2) * t12;
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t26 = g(1) * t20 - g(2) * t23;
t10 = -t12 * t18 + t13 * t29;
t9 = t12 * t21 + t13 * t30;
t8 = t12 * t29 + t13 * t18;
t7 = -t12 * t30 + t13 * t21;
t6 = -t12 * t14 + t13 * t31;
t5 = t12 * t15 + t13 * t32;
t4 = t12 * t31 + t13 * t14;
t3 = -t12 * t32 + t13 * t15;
t2 = g(1) * t4 - g(2) * t6 + t15 * t33;
t1 = -g(1) * t3 - g(2) * t5 + t14 * t33;
t11 = [0, t26, g(1) * t23 + g(2) * t20, t26 * pkin(1), -t27, t28, -g(1) * (-t20 * pkin(1) - t12 * pkin(2) + t13 * qJ(3)) - g(2) * (t23 * pkin(1) + t13 * pkin(2) + t12 * qJ(3)) 0, 0, 0, 0, 0, t28 * t19, t28 * t22, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t27 * t19 + t33, 0, 0, 0, 0, 0, -t25 * t21, t25 * t18, 0, 0, 0, 0, 0, -t25 * t15, t25 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t18 * t33, g(1) * t8 - g(2) * t10 + t21 * t33, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t11;
