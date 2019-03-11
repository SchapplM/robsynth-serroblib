% Calculate minimal parameter regressor of gravitation load for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t19 = qJ(4) + pkin(10);
t12 = sin(t19);
t13 = cos(t19);
t33 = sin(pkin(9));
t34 = cos(pkin(9));
t38 = sin(qJ(1));
t39 = cos(qJ(1));
t4 = -t38 * t33 - t39 * t34;
t5 = t39 * t33 - t38 * t34;
t30 = g(1) * t4 + g(2) * t5;
t42 = -g(3) * t13 + t30 * t12;
t41 = g(3) * t12;
t21 = sin(qJ(6));
t37 = t13 * t21;
t23 = cos(qJ(6));
t36 = t13 * t23;
t35 = t39 * pkin(1) + t38 * qJ(2);
t32 = t39 * pkin(2) + t35;
t31 = g(1) * t5 - g(2) * t4;
t29 = -t38 * pkin(1) + t39 * qJ(2);
t28 = t4 * t21 + t5 * t36;
t27 = -t4 * t23 + t5 * t37;
t26 = -t38 * pkin(2) + t29;
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t25 = g(3) * t24 - t30 * t22;
t20 = -qJ(5) - pkin(7);
t11 = t24 * pkin(4) + pkin(3);
t7 = g(1) * t39 + g(2) * t38;
t6 = g(1) * t38 - g(2) * t39;
t2 = t5 * t21 - t4 * t36;
t1 = t5 * t23 + t4 * t37;
t3 = [0, t6, t7, t6, -t7, -g(1) * t29 - g(2) * t35, -t31, t30, -g(1) * t26 - g(2) * t32, 0, 0, 0, 0, 0, -t31 * t24, t31 * t22, -t30, -g(1) * (t5 * t11 - t4 * t20 + t26) - g(2) * (-t4 * t11 - t5 * t20 + t32) 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t2, g(1) * t27 - g(2) * t1; 0, 0, 0, 0, 0, -t6, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -g(3) * t22 - t30 * t24, 0, t25 * pkin(4), 0, 0, 0, 0, 0, -t42 * t23, t42 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t27 - t21 * t41, g(1) * t2 - g(2) * t28 - t23 * t41;];
taug_reg  = t3;
