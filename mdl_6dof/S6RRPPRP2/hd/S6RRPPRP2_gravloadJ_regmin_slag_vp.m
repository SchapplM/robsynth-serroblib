% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t24 = qJ(2) + pkin(9);
t20 = sin(t24);
t21 = cos(t24);
t58 = t21 * pkin(3) + t20 * qJ(4);
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t56 = -g(1) * t32 - g(2) * t29;
t30 = cos(qJ(5));
t43 = t32 * t30;
t27 = sin(qJ(5));
t47 = t29 * t27;
t4 = t20 * t43 - t47;
t49 = g(3) * t21;
t44 = t32 * t27;
t46 = t29 * t30;
t6 = t20 * t46 + t44;
t55 = -g(1) * t4 - g(2) * t6 + t30 * t49;
t28 = sin(qJ(2));
t52 = pkin(2) * t28;
t48 = t21 * t27;
t26 = -qJ(3) - pkin(7);
t45 = t32 * t26;
t42 = t30 * pkin(5) + pkin(4) - t26;
t41 = qJ(4) * t21;
t31 = cos(qJ(2));
t22 = t31 * pkin(2);
t39 = t22 + t58;
t19 = t22 + pkin(1);
t14 = t32 * t19;
t38 = g(2) * (t58 * t32 + t14);
t37 = -pkin(3) * t20 - t52;
t11 = g(1) * t29 - g(2) * t32;
t25 = -qJ(6) - pkin(8);
t36 = t20 * t27 * pkin(5) - t21 * t25;
t35 = -t19 - t58;
t2 = -g(3) * t20 + t21 * t56;
t34 = -g(3) * t31 - t28 * t56;
t10 = t32 * t41;
t8 = t29 * t41;
t7 = -t20 * t47 + t43;
t5 = t20 * t44 + t46;
t3 = t11 * t21;
t1 = -t20 * t56 - t49;
t9 = [0, t11, -t56, 0, 0, 0, 0, 0, t11 * t31, -t11 * t28, t56, -g(1) * (-t29 * t19 - t45) - g(2) * (-t29 * t26 + t14) t56, -t3, t11 * t20, g(1) * t45 - t38 + (-g(1) * t35 + g(2) * t26) * t29, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t5, g(1) * t6 - g(2) * t4, t3, -t38 + (-g(1) * t42 - g(2) * t36) * t32 + (-g(1) * (t35 - t36) - g(2) * t42) * t29; 0, 0, 0, 0, 0, 0, 0, 0, t34, g(3) * t28 - t31 * t56, 0, t34 * pkin(2), 0, -t1, t2, -g(1) * (t37 * t32 + t10) - g(2) * (t37 * t29 + t8) - g(3) * t39, 0, 0, 0, 0, 0, t2 * t27, t2 * t30, t1, -g(1) * t10 - g(2) * t8 - g(3) * (t36 + t39) + t56 * (pkin(5) * t48 - t52 + (-pkin(3) + t25) * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, g(1) * t5 - g(2) * t7 - g(3) * t48, 0, t55 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg  = t9;
