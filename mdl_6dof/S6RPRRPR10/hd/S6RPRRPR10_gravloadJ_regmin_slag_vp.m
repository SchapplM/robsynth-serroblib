% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% taug_reg [6x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:55:59
% EndTime: 2019-05-05 23:56:01
% DurationCPUTime: 0.33s
% Computational Cost: add. (155->62), mult. (409->89), div. (0->0), fcn. (471->8), ass. (0->42)
t26 = sin(qJ(3));
t30 = cos(qJ(3));
t27 = sin(qJ(1));
t31 = cos(qJ(1));
t61 = -g(1) * t27 + g(2) * t31;
t64 = g(3) * t26 + t61 * t30;
t29 = cos(qJ(4));
t41 = t31 * t29;
t25 = sin(qJ(4));
t44 = t27 * t25;
t11 = t26 * t44 - t41;
t42 = t31 * t25;
t43 = t27 * t29;
t12 = t26 * t43 + t42;
t24 = sin(qJ(6));
t28 = cos(qJ(6));
t2 = t11 * t28 - t12 * t24;
t36 = t24 * t29 - t25 * t28;
t13 = t26 * t42 + t43;
t14 = t26 * t41 - t44;
t39 = -t13 * t28 + t14 * t24;
t48 = g(3) * t30;
t65 = -g(1) * t2 - g(2) * t39 + t36 * t48;
t46 = t30 * pkin(8);
t62 = t26 * pkin(3) - t46;
t3 = t11 * t24 + t12 * t28;
t35 = t24 * t25 + t28 * t29;
t37 = t13 * t24 + t14 * t28;
t57 = -g(1) * t3 + g(2) * t37 - t35 * t48;
t40 = t31 * pkin(1) + t27 * qJ(2);
t38 = g(1) * t13 + g(2) * t11;
t17 = g(1) * t31 + g(2) * t27;
t34 = pkin(4) * t29 + qJ(5) * t25 + pkin(3);
t1 = g(1) * t11 - g(2) * t13 + t25 * t48;
t33 = g(1) * t12 - g(2) * t14 + t29 * t48;
t21 = t31 * qJ(2);
t15 = t17 * t30;
t7 = -t26 * t61 + t48;
t6 = t64 * t29;
t5 = t64 * t25;
t4 = -g(1) * t14 - g(2) * t12;
t8 = [0, -t61, t17, t61, -t17, -g(1) * (-t27 * pkin(1) + t21) - g(2) * t40, 0, 0, 0, 0, 0, -t17 * t26, -t15, 0, 0, 0, 0, 0, t4, t38, t4, t15, -t38, -g(1) * (t14 * pkin(4) + t13 * qJ(5) + t62 * t31 + t21) - g(2) * (t12 * pkin(4) + t31 * pkin(7) + t11 * qJ(5) + t40) + (-g(1) * (-pkin(1) - pkin(7)) - g(2) * t62) * t27, 0, 0, 0, 0, 0, -g(1) * t37 - g(2) * t3, g(1) * t39 - g(2) * t2; 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t7, 0, 0, 0, 0, 0, t6, -t5, t6, -t7, t5, -g(3) * (-t34 * t26 + t46) + t61 * (pkin(8) * t26 + t34 * t30) 0, 0, 0, 0, 0, t64 * t35, -t64 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t33, t1, 0, -t33, -g(1) * (-t11 * pkin(4) + t12 * qJ(5)) - g(2) * (t13 * pkin(4) - t14 * qJ(5)) - (-pkin(4) * t25 + qJ(5) * t29) * t48, 0, 0, 0, 0, 0, -t65, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t57;];
taug_reg  = t8;
