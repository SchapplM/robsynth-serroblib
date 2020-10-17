% Calculate minimal parameter regressor of gravitation load for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:02:48
% EndTime: 2019-05-05 14:02:49
% DurationCPUTime: 0.19s
% Computational Cost: add. (182->49), mult. (156->60), div. (0->0), fcn. (155->10), ass. (0->33)
t19 = qJ(1) + pkin(9);
t14 = sin(t19);
t16 = cos(t19);
t33 = g(1) * t16 + g(2) * t14;
t18 = pkin(10) + qJ(4);
t13 = sin(t18);
t15 = cos(t18);
t2 = g(3) * t13 + t33 * t15;
t39 = g(3) * t15;
t24 = sin(qJ(1));
t38 = t24 * pkin(1);
t23 = sin(qJ(6));
t37 = t14 * t23;
t25 = cos(qJ(6));
t36 = t14 * t25;
t35 = t16 * t23;
t34 = t16 * t25;
t32 = g(1) * t14 - g(2) * t16;
t26 = cos(qJ(1));
t31 = g(1) * t24 - g(2) * t26;
t30 = t15 * pkin(4) + t13 * qJ(5);
t21 = cos(pkin(10));
t28 = t21 * pkin(3) + pkin(2) + t30;
t22 = -pkin(7) - qJ(3);
t17 = t26 * pkin(1);
t8 = -t13 * t37 + t34;
t7 = t13 * t36 + t35;
t6 = t13 * t35 + t36;
t5 = t13 * t34 - t37;
t4 = t32 * t15;
t3 = t32 * t13;
t1 = t33 * t13 - t39;
t9 = [0, t31, g(1) * t26 + g(2) * t24, t31 * pkin(1), t32 * t21, -t32 * sin(pkin(10)) -t33, -g(1) * (-t14 * pkin(2) + t16 * qJ(3) - t38) - g(2) * (t16 * pkin(2) + t14 * qJ(3) + t17) 0, 0, 0, 0, 0, t4, -t3, -t33, -t4, t3, g(1) * t38 - g(2) * t17 + (g(1) * t22 - g(2) * t28) * t16 + (g(1) * t28 + g(2) * t22) * t14, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1, -t2, -g(3) * t30 + t33 * (pkin(4) * t13 - qJ(5) * t15) 0, 0, 0, 0, 0, -t2 * t23, -t2 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t25 * t39, g(1) * t6 - g(2) * t8 - t23 * t39;];
taug_reg  = t9;
