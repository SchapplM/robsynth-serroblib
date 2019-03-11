% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = qJ(1) + pkin(10);
t15 = sin(t20);
t16 = cos(t20);
t35 = g(1) * t16 + g(2) * t15;
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t8 = g(3) * t23 + t35 * t26;
t40 = g(3) * t26;
t21 = qJ(5) + qJ(6);
t17 = sin(t21);
t39 = t17 * t23;
t18 = cos(t21);
t38 = t18 * t23;
t22 = sin(qJ(5));
t37 = t22 * t23;
t25 = cos(qJ(5));
t36 = t23 * t25;
t34 = g(1) * t15 - g(2) * t16;
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t33 = g(1) * t24 - g(2) * t27;
t32 = t26 * pkin(3) + t23 * qJ(4);
t30 = pkin(2) + t32;
t29 = t33 * pkin(1);
t14 = t34 * t26;
t13 = t34 * t23;
t12 = -t15 * t37 + t16 * t25;
t11 = t15 * t36 + t16 * t22;
t10 = t15 * t25 + t16 * t37;
t9 = -t15 * t22 + t16 * t36;
t7 = t35 * t23 - t40;
t6 = -t15 * t39 + t16 * t18;
t5 = t15 * t38 + t16 * t17;
t4 = t15 * t18 + t16 * t39;
t3 = -t15 * t17 + t16 * t38;
t2 = g(1) * t4 - g(2) * t6 - t17 * t40;
t1 = -g(1) * t3 - g(2) * t5 + t18 * t40;
t19 = [0, t33, g(1) * t27 + g(2) * t24, t29, 0, 0, 0, 0, 0, t14, -t13, -t35, -t14, t13, t29 + (-g(1) * pkin(7) - g(2) * t30) * t16 + (-g(2) * pkin(7) + g(1) * t30) * t15, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, -t7, -t8, -g(3) * t32 + t35 * (pkin(3) * t23 - qJ(4) * t26) 0, 0, 0, 0, 0, -t8 * t22, -t8 * t25, 0, 0, 0, 0, 0, -t8 * t17, -t8 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11 + t25 * t40, g(1) * t10 - g(2) * t12 - t22 * t40, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t19;
