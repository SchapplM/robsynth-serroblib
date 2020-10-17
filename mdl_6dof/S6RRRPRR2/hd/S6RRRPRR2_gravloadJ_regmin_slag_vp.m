% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:58:56
% EndTime: 2019-05-07 09:58:56
% DurationCPUTime: 0.20s
% Computational Cost: add. (254->46), mult. (239->72), div. (0->0), fcn. (251->12), ass. (0->47)
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t20 = g(1) * t38 + g(2) * t35;
t32 = qJ(2) + qJ(3);
t24 = pkin(11) + t32;
t21 = sin(t24);
t22 = cos(t24);
t40 = -g(3) * t22 + t20 * t21;
t51 = g(3) * t21;
t31 = qJ(5) + qJ(6);
t25 = sin(t31);
t49 = t35 * t25;
t27 = cos(t31);
t48 = t35 * t27;
t33 = sin(qJ(5));
t47 = t35 * t33;
t36 = cos(qJ(5));
t46 = t35 * t36;
t45 = t38 * t25;
t44 = t38 * t27;
t43 = t38 * t33;
t42 = t38 * t36;
t28 = cos(t32);
t37 = cos(qJ(2));
t41 = t37 * pkin(2) + pkin(3) * t28;
t19 = g(1) * t35 - g(2) * t38;
t26 = sin(t32);
t7 = -g(3) * t28 + t20 * t26;
t34 = sin(qJ(2));
t30 = -qJ(4) - pkin(8) - pkin(7);
t17 = pkin(1) + t41;
t16 = t22 * t42 + t47;
t15 = -t22 * t43 + t46;
t14 = -t22 * t46 + t43;
t13 = t22 * t47 + t42;
t12 = t22 * t44 + t49;
t11 = -t22 * t45 + t48;
t10 = -t22 * t48 + t45;
t9 = t22 * t49 + t44;
t8 = g(3) * t26 + t20 * t28;
t6 = t40 * t36;
t5 = t40 * t33;
t4 = t40 * t27;
t3 = t40 * t25;
t2 = g(1) * t12 - g(2) * t10 + t27 * t51;
t1 = -g(1) * t11 + g(2) * t9 + t25 * t51;
t18 = [0, t19, t20, 0, 0, 0, 0, 0, t19 * t37, -t19 * t34, 0, 0, 0, 0, 0, t19 * t28, -t19 * t26, -t20, -g(1) * (-t35 * t17 - t38 * t30) - g(2) * (t38 * t17 - t35 * t30) 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t37 + t20 * t34, g(3) * t34 + t20 * t37, 0, 0, 0, 0, 0, t7, t8, 0, -g(3) * t41 - t20 * (-t34 * pkin(2) - pkin(3) * t26) 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t7 * pkin(3), 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t33 * t51, g(1) * t16 - g(2) * t14 + t36 * t51, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t18;
