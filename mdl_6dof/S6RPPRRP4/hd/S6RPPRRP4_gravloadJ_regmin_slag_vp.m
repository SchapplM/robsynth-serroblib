% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRP4
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
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t49 = sin(qJ(1));
t50 = cos(qJ(1));
t16 = -t49 * t44 - t50 * t45;
t17 = t50 * t44 - t49 * t45;
t39 = g(1) * t16 + g(2) * t17;
t29 = sin(qJ(4));
t52 = g(3) * t29;
t51 = t29 * pkin(8);
t28 = sin(qJ(5));
t31 = cos(qJ(4));
t48 = t28 * t31;
t30 = cos(qJ(5));
t47 = t30 * t31;
t46 = t50 * pkin(1) + t49 * qJ(2);
t43 = t50 * pkin(2) + t46;
t10 = -t16 * t48 - t17 * t30;
t6 = -t16 * t30 + t17 * t48;
t42 = g(1) * t6 + g(2) * t10;
t41 = -t49 * pkin(1) + t50 * qJ(2);
t40 = g(1) * t17 - g(2) * t16;
t38 = t31 * pkin(4) + pkin(3) + t51;
t37 = pkin(5) * t30 + qJ(6) * t28 + pkin(4);
t7 = t16 * t28 + t17 * t47;
t36 = -g(1) * t10 + g(2) * t6 + t28 * t52;
t11 = -t16 * t47 + t17 * t28;
t35 = -g(1) * t11 + g(2) * t7 + t30 * t52;
t34 = -t49 * pkin(2) + t41;
t33 = g(3) * t31 - t29 * t39;
t19 = g(1) * t50 + g(2) * t49;
t18 = g(1) * t49 - g(2) * t50;
t12 = t40 * t29;
t5 = -t31 * t39 - t52;
t4 = t33 * t30;
t3 = t33 * t28;
t2 = -g(1) * t7 - g(2) * t11;
t1 = [0, t18, t19, t18, -t19, -g(1) * t41 - g(2) * t46, -t40, t39, -g(1) * t34 - g(2) * t43, 0, 0, 0, 0, 0, -t40 * t31, t12, 0, 0, 0, 0, 0, t2, t42, t2, -t12, -t42, -g(1) * (t7 * pkin(5) + t6 * qJ(6) + t34) - g(2) * (t11 * pkin(5) + t10 * qJ(6) + t43) + (-g(2) * pkin(7) - g(1) * t38) * t17 + (-g(1) * pkin(7) + g(2) * t38) * t16; 0, 0, 0, 0, 0, -t18, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t5, t3, -g(3) * (-t31 * t37 - t51) + t39 * (pkin(8) * t31 - t29 * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t35, -t36, 0, t35, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (pkin(5) * t6 - qJ(6) * t7) - (pkin(5) * t28 - qJ(6) * t30) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36;];
taug_reg  = t1;
