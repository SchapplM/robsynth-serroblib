% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR7
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t11 = g(1) * t26 - g(2) * t29;
t23 = qJ(3) + qJ(4);
t15 = pkin(10) + t23;
t13 = sin(t15);
t14 = cos(t15);
t37 = -g(3) * t13 + t11 * t14;
t35 = g(3) * t14;
t24 = sin(qJ(6));
t34 = t26 * t24;
t27 = cos(qJ(6));
t33 = t26 * t27;
t32 = t29 * t24;
t31 = t29 * t27;
t30 = t29 * pkin(1) + t26 * qJ(2);
t12 = g(1) * t29 + g(2) * t26;
t16 = sin(t23);
t17 = cos(t23);
t4 = g(3) * t16 - t11 * t17;
t28 = cos(qJ(3));
t25 = sin(qJ(3));
t22 = -qJ(5) - pkin(8) - pkin(7);
t19 = t29 * qJ(2);
t9 = t25 * pkin(3) + pkin(4) * t16;
t8 = t13 * t31 - t34;
t7 = t13 * t32 + t33;
t6 = t13 * t33 + t32;
t5 = -t13 * t34 + t31;
t3 = g(3) * t17 + t11 * t16;
t2 = t37 * t27;
t1 = t37 * t24;
t10 = [0, t11, t12, -t11, -t12, -g(1) * (-t26 * pkin(1) + t19) - g(2) * t30, 0, 0, 0, 0, 0, -t12 * t25, -t12 * t28, 0, 0, 0, 0, 0, -t12 * t16, -t12 * t17, t11, -g(1) * (t29 * t9 + t19 + (-pkin(1) + t22) * t26) - g(2) * (-t29 * t22 + t26 * t9 + t30) 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5; 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t25 - t11 * t28, g(3) * t28 + t11 * t25, 0, 0, 0, 0, 0, t4, t3, 0, g(3) * t9 - t11 * (t28 * pkin(3) + pkin(4) * t17) 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, t4 * pkin(4), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t24 * t35, g(1) * t6 - g(2) * t8 + t27 * t35;];
taug_reg  = t10;
