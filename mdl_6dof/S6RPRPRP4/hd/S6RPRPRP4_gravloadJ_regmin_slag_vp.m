% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t29 = sin(qJ(3));
t23 = t29 * qJ(4);
t32 = cos(qJ(3));
t45 = t32 * pkin(3) + t23;
t27 = qJ(1) + pkin(9);
t21 = sin(t27);
t22 = cos(t27);
t56 = -g(1) * t22 - g(2) * t21;
t6 = g(3) * t29 - t56 * t32;
t55 = pkin(3) * t29;
t54 = g(1) * t21;
t50 = g(3) * t32;
t49 = t32 * pkin(8);
t48 = t22 * t32;
t28 = sin(qJ(5));
t47 = t28 * t29;
t31 = cos(qJ(5));
t46 = t29 * t31;
t44 = qJ(4) * t32;
t30 = sin(qJ(1));
t43 = -t30 * pkin(1) + t22 * pkin(7);
t7 = t21 * t28 - t22 * t46;
t9 = t21 * t46 + t22 * t28;
t42 = g(1) * t9 + g(2) * t7;
t33 = cos(qJ(1));
t41 = t33 * pkin(1) + pkin(3) * t48 + t21 * pkin(7) + (pkin(2) + t23) * t22;
t39 = -g(2) * t22 + t54;
t38 = g(1) * t30 - g(2) * t33;
t37 = pkin(5) * t28 - qJ(6) * t31;
t36 = -pkin(2) - t45;
t1 = g(1) * t7 - g(2) * t9 + t31 * t50;
t10 = -t21 * t47 + t22 * t31;
t8 = t21 * t31 + t22 * t47;
t35 = -g(1) * t8 + g(2) * t10 + t28 * t50;
t16 = t22 * t44;
t14 = t21 * t44;
t12 = t39 * t32;
t11 = t39 * t29;
t5 = -t29 * t56 - t50;
t4 = t6 * t31;
t3 = t6 * t28;
t2 = -g(1) * t10 - g(2) * t8;
t13 = [0, t38, g(1) * t33 + g(2) * t30, t38 * pkin(1), 0, 0, 0, 0, 0, t12, -t11, t56, -t12, t11, -g(1) * t43 - g(2) * t41 - t36 * t54, 0, 0, 0, 0, 0, t2, t42, t2, t12, -t42, -g(1) * (t22 * pkin(4) + t10 * pkin(5) + t9 * qJ(6) + t43) - g(2) * (t8 * pkin(5) + pkin(8) * t48 + t7 * qJ(6) + t41) + (-g(1) * (t36 - t49) - g(2) * pkin(4)) * t21; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, -t5, -t6, -g(1) * (-t22 * t55 + t16) - g(2) * (-t21 * t55 + t14) - g(3) * t45, 0, 0, 0, 0, 0, -t3, -t4, -t3, t5, t4, -g(1) * t16 - g(2) * t14 - g(3) * (t37 * t29 + t45 + t49) + t56 * (t37 * t32 + (-pkin(3) - pkin(8)) * t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t35, t1, 0, t35, -g(1) * (-t7 * pkin(5) + t8 * qJ(6)) - g(2) * (t9 * pkin(5) - t10 * qJ(6)) - (-pkin(5) * t31 - qJ(6) * t28) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
