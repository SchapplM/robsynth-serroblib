% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:10:28
% EndTime: 2019-05-05 23:10:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (162->48), mult. (227->72), div. (0->0), fcn. (237->8), ass. (0->41)
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t52 = -g(1) * t27 + g(2) * t30;
t26 = sin(qJ(3));
t25 = sin(qJ(4));
t38 = t30 * t25;
t28 = cos(qJ(4));
t41 = t27 * t28;
t11 = t26 * t38 + t41;
t29 = cos(qJ(3));
t45 = g(3) * t29;
t37 = t30 * t28;
t42 = t27 * t25;
t9 = -t26 * t42 + t37;
t51 = -g(1) * t9 - g(2) * t11 + t25 * t45;
t8 = -g(3) * t26 - t29 * t52;
t19 = qJ(4) + pkin(10) + qJ(6);
t16 = sin(t19);
t44 = t27 * t16;
t17 = cos(t19);
t43 = t27 * t17;
t40 = t30 * t16;
t39 = t30 * t17;
t35 = pkin(4) * t25 + pkin(7);
t34 = g(2) * (t30 * pkin(1) + t27 * qJ(2));
t15 = g(1) * t30 + g(2) * t27;
t18 = t28 * pkin(4) + pkin(3);
t24 = -qJ(5) - pkin(8);
t32 = t26 * t18 + t29 * t24;
t21 = t30 * qJ(2);
t13 = t15 * t29;
t12 = t26 * t37 - t42;
t10 = t26 * t41 + t38;
t7 = -t26 * t52 + t45;
t6 = t26 * t39 - t44;
t5 = t26 * t40 + t43;
t4 = t26 * t43 + t40;
t3 = -t26 * t44 + t39;
t2 = g(1) * t4 - g(2) * t6 + t17 * t45;
t1 = -g(1) * t3 - g(2) * t5 + t16 * t45;
t14 = [0, -t52, t15, t52, -t15, -g(1) * (-t27 * pkin(1) + t21) - t34, 0, 0, 0, 0, 0, -t15 * t26, -t13, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t13, -g(1) * t21 - t34 + (-g(1) * t32 - g(2) * t35) * t30 + (-g(1) * (-pkin(1) - t35) - g(2) * t32) * t27, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, -t8 * t28, t8 * t25, -t7, g(3) * t32 + t52 * (t18 * t29 - t24 * t26) 0, 0, 0, 0, 0, -t8 * t17, t8 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, g(1) * t10 - g(2) * t12 + t28 * t45, 0, t51 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t14;
