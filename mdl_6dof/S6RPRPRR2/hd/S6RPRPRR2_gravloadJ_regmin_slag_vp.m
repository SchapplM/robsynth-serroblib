% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:24:30
% EndTime: 2019-05-05 18:24:30
% DurationCPUTime: 0.16s
% Computational Cost: add. (206->44), mult. (172->68), div. (0->0), fcn. (185->12), ass. (0->42)
t18 = qJ(3) + pkin(11);
t12 = sin(t18);
t14 = cos(t18);
t19 = qJ(1) + pkin(10);
t13 = sin(t19);
t15 = cos(t19);
t33 = g(1) * t15 + g(2) * t13;
t30 = -g(3) * t14 + t33 * t12;
t43 = g(3) * t12;
t20 = qJ(5) + qJ(6);
t16 = sin(t20);
t41 = t13 * t16;
t17 = cos(t20);
t40 = t13 * t17;
t22 = sin(qJ(5));
t39 = t13 * t22;
t25 = cos(qJ(5));
t38 = t13 * t25;
t37 = t15 * t16;
t36 = t15 * t17;
t35 = t15 * t22;
t34 = t15 * t25;
t32 = g(1) * t13 - g(2) * t15;
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t31 = g(1) * t24 - g(2) * t27;
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t28 = -g(3) * t26 + t33 * t23;
t21 = -qJ(4) - pkin(7);
t11 = t26 * pkin(3) + pkin(2);
t10 = t14 * t34 + t39;
t9 = -t14 * t35 + t38;
t8 = -t14 * t38 + t35;
t7 = t14 * t39 + t34;
t6 = t14 * t36 + t41;
t5 = -t14 * t37 + t40;
t4 = -t14 * t40 + t37;
t3 = t14 * t41 + t36;
t2 = g(1) * t6 - g(2) * t4 + t17 * t43;
t1 = -g(1) * t5 + g(2) * t3 + t16 * t43;
t29 = [0, t31, g(1) * t27 + g(2) * t24, t31 * pkin(1), 0, 0, 0, 0, 0, t32 * t26, -t32 * t23, -t33, -g(1) * (-t24 * pkin(1) - t13 * t11 - t15 * t21) - g(2) * (t27 * pkin(1) + t15 * t11 - t13 * t21) 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, g(3) * t23 + t33 * t26, 0, t28 * pkin(3), 0, 0, 0, 0, 0, t30 * t25, -t30 * t22, 0, 0, 0, 0, 0, t30 * t17, -t30 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t22 * t43, g(1) * t10 - g(2) * t8 + t25 * t43, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t29;
