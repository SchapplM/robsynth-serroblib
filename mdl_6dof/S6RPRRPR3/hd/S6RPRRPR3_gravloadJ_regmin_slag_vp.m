% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:17:58
% EndTime: 2019-05-05 22:17:59
% DurationCPUTime: 0.28s
% Computational Cost: add. (295->64), mult. (399->93), div. (0->0), fcn. (461->10), ass. (0->42)
t23 = sin(qJ(3));
t27 = cos(qJ(3));
t20 = qJ(1) + pkin(10);
t18 = sin(t20);
t19 = cos(t20);
t39 = g(1) * t19 + g(2) * t18;
t59 = g(3) * t27 - t23 * t39;
t26 = cos(qJ(4));
t22 = sin(qJ(4));
t43 = t22 * t27;
t11 = -t18 * t26 + t19 * t43;
t42 = t26 * t27;
t12 = t18 * t22 + t19 * t42;
t21 = sin(qJ(6));
t25 = cos(qJ(6));
t2 = t11 * t25 - t12 * t21;
t35 = t21 * t26 - t22 * t25;
t10 = t18 * t42 - t19 * t22;
t9 = t18 * t43 + t19 * t26;
t41 = t10 * t21 - t25 * t9;
t47 = g(3) * t23;
t58 = -g(1) * t2 + g(2) * t41 + t35 * t47;
t3 = t11 * t21 + t12 * t25;
t34 = t21 * t22 + t25 * t26;
t36 = t10 * t25 + t21 * t9;
t54 = g(1) * t3 + g(2) * t36 + t34 * t47;
t45 = t23 * pkin(8);
t40 = g(1) * t9 - g(2) * t11;
t38 = g(1) * t18 - g(2) * t19;
t24 = sin(qJ(1));
t28 = cos(qJ(1));
t37 = g(1) * t24 - g(2) * t28;
t33 = pkin(3) * t27 + pkin(2) + t45;
t32 = pkin(4) * t26 + qJ(5) * t22 + pkin(3);
t1 = g(1) * t11 + g(2) * t9 + t22 * t47;
t31 = g(1) * t12 + g(2) * t10 + t26 * t47;
t13 = t38 * t23;
t7 = t27 * t39 + t47;
t6 = t59 * t26;
t5 = t59 * t22;
t4 = g(1) * t10 - g(2) * t12;
t8 = [0, t37, g(1) * t28 + g(2) * t24, t37 * pkin(1), 0, 0, 0, 0, 0, t38 * t27, -t13, 0, 0, 0, 0, 0, t4, -t40, t4, t13, t40, -g(1) * (-t24 * pkin(1) - t10 * pkin(4) - t9 * qJ(5)) - g(2) * (t28 * pkin(1) + t12 * pkin(4) + t11 * qJ(5)) + (-g(1) * pkin(7) - g(2) * t33) * t19 + (-g(2) * pkin(7) + g(1) * t33) * t18, 0, 0, 0, 0, 0, g(1) * t36 - g(2) * t3, -g(1) * t41 - g(2) * t2; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t7, 0, 0, 0, 0, 0, -t6, t5, -t6, -t7, -t5, -g(3) * (t27 * t32 + t45) - t39 * (pkin(8) * t27 - t23 * t32) 0, 0, 0, 0, 0, -t59 * t34, t59 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t31, t1, 0, -t31, -g(1) * (-pkin(4) * t11 + qJ(5) * t12) - g(2) * (-pkin(4) * t9 + qJ(5) * t10) - (-pkin(4) * t22 + qJ(5) * t26) * t47, 0, 0, 0, 0, 0, -t58, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t54;];
taug_reg  = t8;
