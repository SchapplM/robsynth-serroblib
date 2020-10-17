% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x34]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:28:52
% EndTime: 2019-05-06 04:28:53
% DurationCPUTime: 0.21s
% Computational Cost: add. (214->48), mult. (250->76), div. (0->0), fcn. (284->10), ass. (0->45)
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t17 = g(1) * t28 - g(2) * t31;
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t33 = -g(3) * t27 + t17 * t30;
t46 = g(3) * t30;
t25 = qJ(4) + qJ(5);
t23 = qJ(6) + t25;
t19 = sin(t23);
t45 = t28 * t19;
t20 = cos(t23);
t44 = t28 * t20;
t21 = sin(t25);
t43 = t28 * t21;
t22 = cos(t25);
t42 = t28 * t22;
t26 = sin(qJ(4));
t41 = t28 * t26;
t29 = cos(qJ(4));
t40 = t28 * t29;
t39 = t31 * t19;
t38 = t31 * t20;
t37 = t31 * t21;
t36 = t31 * t22;
t35 = t31 * t26;
t34 = t31 * t29;
t18 = g(1) * t31 + g(2) * t28;
t16 = t27 * t34 - t41;
t15 = t27 * t35 + t40;
t14 = t27 * t40 + t35;
t13 = -t27 * t41 + t34;
t12 = t27 * t36 - t43;
t11 = t27 * t37 + t42;
t10 = t27 * t42 + t37;
t9 = -t27 * t43 + t36;
t8 = t27 * t38 - t45;
t7 = t27 * t39 + t44;
t6 = t27 * t44 + t39;
t5 = -t27 * t45 + t38;
t4 = g(1) * t10 - g(2) * t12 + t22 * t46;
t3 = -g(1) * t9 - g(2) * t11 + t21 * t46;
t2 = g(1) * t6 - g(2) * t8 + t20 * t46;
t1 = -g(1) * t5 - g(2) * t7 + t19 * t46;
t24 = [0, t17, t18, -t17, -t18, -g(1) * (-t28 * pkin(1) + t31 * qJ(2)) - g(2) * (t31 * pkin(1) + t28 * qJ(2)) 0, 0, 0, 0, 0, -t18 * t27, -t18 * t30, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14, g(1) * t15 - g(2) * t13, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5; 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t17 * t27 + t46, 0, 0, 0, 0, 0, -t33 * t29, t33 * t26, 0, 0, 0, 0, 0, -t33 * t22, t33 * t21, 0, 0, 0, 0, 0, -t33 * t20, t33 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15 + t26 * t46, g(1) * t14 - g(2) * t16 + t29 * t46, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t24;
