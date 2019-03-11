% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t44 = -g(1) * t23 + g(2) * t25;
t22 = sin(qJ(5));
t17 = pkin(9) + qJ(4);
t12 = cos(t17);
t37 = g(3) * t12;
t11 = sin(t17);
t24 = cos(qJ(5));
t33 = t25 * t24;
t36 = t23 * t22;
t4 = -t11 * t36 + t33;
t34 = t25 * t22;
t35 = t23 * t24;
t6 = t11 * t34 + t35;
t43 = -g(1) * t4 - g(2) * t6 + t22 * t37;
t2 = -g(3) * t11 - t12 * t44;
t32 = t25 * pkin(1) + t23 * qJ(2);
t30 = g(2) * t32;
t29 = pkin(5) * t22 + pkin(7) + qJ(3);
t9 = g(1) * t25 + g(2) * t23;
t10 = t24 * pkin(5) + pkin(4);
t20 = -qJ(6) - pkin(8);
t27 = -t11 * t10 - t12 * t20;
t18 = sin(pkin(9));
t26 = pkin(3) * t18 - t27;
t14 = t25 * qJ(2);
t7 = t11 * t33 - t36;
t5 = t11 * t35 + t34;
t3 = t9 * t12;
t1 = -t11 * t44 + t37;
t8 = [0, -t44, t9, t44, -t9, -g(1) * (-t23 * pkin(1) + t14) - t30, -t9 * t18, -t9 * cos(pkin(9)) -t44, -g(1) * (t14 + (-pkin(1) - qJ(3)) * t23) - g(2) * (t25 * qJ(3) + t32) 0, 0, 0, 0, 0, -t9 * t11, -t3, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t5, g(1) * t6 - g(2) * t4, t3, -g(1) * t14 - t30 + (-g(1) * t26 - g(2) * t29) * t25 + (-g(1) * (-pkin(1) - t29) - g(2) * t26) * t23; 0, 0, 0, 0, 0, t44, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, -t2 * t24, t2 * t22, -t1, -g(3) * t27 + t44 * (t10 * t12 - t11 * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, g(1) * t5 - g(2) * t7 + t24 * t37, 0, t43 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg  = t8;
