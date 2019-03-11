% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t36 = cos(qJ(2));
t30 = qJ(2) + qJ(3);
t27 = sin(t30);
t28 = cos(t30);
t49 = t28 * pkin(3) + t27 * qJ(4);
t55 = t36 * pkin(2) + t49;
t31 = sin(qJ(6));
t37 = cos(qJ(1));
t32 = sin(qJ(5));
t51 = cos(qJ(5));
t54 = t27 * t51 - t28 * t32;
t11 = t54 * t37;
t41 = t27 * t32 + t28 * t51;
t34 = sin(qJ(1));
t9 = t54 * t34;
t40 = -g(1) * t11 - g(2) * t9 + g(3) * t41;
t1 = t40 * t31;
t35 = cos(qJ(6));
t2 = t40 * t35;
t53 = pkin(3) * t27;
t52 = g(3) * t54;
t48 = qJ(4) * t28;
t33 = sin(qJ(2));
t46 = -pkin(2) * t33 - t53;
t19 = g(1) * t37 + g(2) * t34;
t45 = g(1) * t34 - g(2) * t37;
t10 = t41 * t34;
t44 = t10 * t35 + t37 * t31;
t43 = t10 * t31 - t37 * t35;
t42 = pkin(1) + t55;
t12 = t41 * t37;
t39 = g(1) * t12 + g(2) * t10 + t52;
t38 = -pkin(8) - pkin(7);
t21 = t37 * t48;
t20 = t34 * t48;
t14 = t45 * t28;
t13 = t45 * t27;
t8 = g(3) * t27 + t19 * t28;
t7 = -g(3) * t28 + t19 * t27;
t6 = t12 * t35 - t34 * t31;
t5 = -t12 * t31 - t34 * t35;
t3 = [0, t45, t19, 0, 0, 0, 0, 0, t45 * t36, -t45 * t33, 0, 0, 0, 0, 0, t14, -t13, t14, -t19, t13 (g(1) * t38 - g(2) * t42) * t37 + (g(1) * t42 + g(2) * t38) * t34, 0, 0, 0, 0, 0, g(1) * t10 - g(2) * t12, g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t6, -g(1) * t43 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t36 + t19 * t33, g(3) * t33 + t19 * t36, 0, 0, 0, 0, 0, t7, t8, t7, 0, -t8, -g(1) * (t46 * t37 + t21) - g(2) * (t46 * t34 + t20) - g(3) * t55, 0, 0, 0, 0, 0, -t40, -t39, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, t7, 0, -t8, -g(1) * (-t37 * t53 + t21) - g(2) * (-t34 * t53 + t20) - g(3) * t49, 0, 0, 0, 0, 0, -t40, -t39, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t43 + t31 * t52, g(1) * t6 + g(2) * t44 + t35 * t52;];
taug_reg  = t3;
