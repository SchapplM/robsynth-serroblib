% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t27 = qJ(2) + qJ(3);
t25 = qJ(4) + t27;
t19 = sin(t25);
t20 = cos(t25);
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t36 = g(1) * t33 + g(2) * t30;
t7 = -g(3) * t20 + t19 * t36;
t46 = g(3) * t19;
t26 = qJ(5) + qJ(6);
t21 = sin(t26);
t44 = t30 * t21;
t23 = cos(t26);
t43 = t30 * t23;
t28 = sin(qJ(5));
t42 = t30 * t28;
t31 = cos(qJ(5));
t41 = t30 * t31;
t40 = t33 * t21;
t39 = t33 * t23;
t38 = t33 * t28;
t37 = t33 * t31;
t35 = g(1) * t30 - g(2) * t33;
t32 = cos(qJ(2));
t29 = sin(qJ(2));
t24 = cos(t27);
t22 = sin(t27);
t18 = t20 * t37 + t42;
t17 = -t20 * t38 + t41;
t16 = -t20 * t41 + t38;
t15 = t20 * t42 + t37;
t14 = t20 * t39 + t44;
t13 = -t20 * t40 + t43;
t12 = -t20 * t43 + t40;
t11 = t20 * t44 + t39;
t10 = g(3) * t22 + t24 * t36;
t9 = -g(3) * t24 + t22 * t36;
t8 = t20 * t36 + t46;
t6 = t7 * t31;
t5 = t7 * t28;
t4 = t7 * t23;
t3 = t7 * t21;
t2 = g(1) * t14 - g(2) * t12 + t23 * t46;
t1 = -g(1) * t13 + g(2) * t11 + t21 * t46;
t34 = [0, t35, t36, 0, 0, 0, 0, 0, t35 * t32, -t35 * t29, 0, 0, 0, 0, 0, t35 * t24, -t35 * t22, 0, 0, 0, 0, 0, t35 * t20, -t35 * t19, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t32 + t29 * t36, g(3) * t29 + t32 * t36, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t17 + g(2) * t15 + t28 * t46, g(1) * t18 - g(2) * t16 + t31 * t46, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t34;
