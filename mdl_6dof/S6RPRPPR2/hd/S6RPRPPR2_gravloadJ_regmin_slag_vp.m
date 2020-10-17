% Calculate minimal parameter regressor of gravitation load for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:34:18
% EndTime: 2019-05-05 16:34:18
% DurationCPUTime: 0.17s
% Computational Cost: add. (179->46), mult. (163->64), div. (0->0), fcn. (159->10), ass. (0->34)
t16 = qJ(3) + pkin(10);
t10 = sin(t16);
t12 = cos(t16);
t28 = t12 * pkin(4) + t10 * qJ(5);
t17 = qJ(1) + pkin(9);
t11 = sin(t17);
t13 = cos(t17);
t32 = g(1) * t13 + g(2) * t11;
t40 = g(3) * t12;
t24 = cos(qJ(1));
t23 = cos(qJ(3));
t14 = t23 * pkin(3);
t9 = t14 + pkin(2);
t38 = t24 * pkin(1) + t13 * t9;
t19 = sin(qJ(6));
t37 = t11 * t19;
t22 = cos(qJ(6));
t36 = t11 * t22;
t35 = t13 * t19;
t34 = t13 * t22;
t31 = g(1) * t11 - g(2) * t13;
t21 = sin(qJ(1));
t30 = g(1) * t21 - g(2) * t24;
t18 = -qJ(4) - pkin(7);
t29 = -t21 * pkin(1) - t13 * t18;
t26 = -g(3) * t10 - t32 * t12;
t20 = sin(qJ(3));
t25 = -g(3) * t23 + t32 * t20;
t5 = -t10 * t37 + t34;
t4 = t10 * t36 + t35;
t3 = t10 * t35 + t36;
t2 = t10 * t34 - t37;
t1 = -t32 * t10 + t40;
t6 = [0, t30, g(1) * t24 + g(2) * t21, t30 * pkin(1), 0, 0, 0, 0, 0, t31 * t23, -t31 * t20, -t32, -g(1) * (-t11 * t9 + t29) - g(2) * (-t11 * t18 + t38) -t32, -t31 * t12, t31 * t10, -g(1) * t29 - g(2) * (t28 * t13 + t38) + (-g(1) * (-t28 - t9) + g(2) * t18) * t11, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t3, g(1) * t4 - g(2) * t2; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t20 + t32 * t23, 0, t25 * pkin(3), 0, t1, t26, -g(3) * (t14 + t28) + t32 * (pkin(3) * t20 + pkin(4) * t10 - qJ(5) * t12) 0, 0, 0, 0, 0, t26 * t19, t26 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, 0, 0, -t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4 + t22 * t40, g(1) * t3 - g(2) * t5 - t19 * t40;];
taug_reg  = t6;
