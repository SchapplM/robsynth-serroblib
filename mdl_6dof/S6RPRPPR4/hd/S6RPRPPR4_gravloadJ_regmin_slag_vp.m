% Calculate minimal parameter regressor of gravitation load for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = cos(pkin(9));
t29 = pkin(9) + qJ(3);
t26 = sin(t29);
t27 = cos(t29);
t50 = pkin(3) * t27 + qJ(4) * t26;
t62 = pkin(2) * t33 + pkin(1) + t50;
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t18 = g(1) * t38 + g(2) * t36;
t61 = t18 * t26;
t6 = -g(3) * t27 + t61;
t60 = pkin(3) * t26;
t57 = g(3) * t26;
t30 = sin(pkin(10));
t55 = t36 * t30;
t32 = cos(pkin(10));
t54 = t36 * t32;
t53 = t38 * t30;
t52 = t38 * t32;
t34 = -pkin(7) - qJ(2);
t51 = t38 * t34;
t49 = qJ(4) * t27;
t48 = t62 * t38;
t11 = t27 * t53 - t54;
t9 = t27 * t55 + t52;
t47 = g(1) * t9 - g(2) * t11;
t17 = g(1) * t36 - g(2) * t38;
t10 = t27 * t54 - t53;
t35 = sin(qJ(6));
t37 = cos(qJ(6));
t46 = t10 * t37 + t35 * t9;
t45 = t10 * t35 - t37 * t9;
t44 = pkin(4) * t32 + qJ(5) * t30;
t43 = t30 * t37 - t32 * t35;
t42 = t30 * t35 + t32 * t37;
t40 = g(3) * t43;
t39 = (g(1) * t62 + g(2) * t34) * t36;
t15 = t38 * t49;
t13 = t36 * t49;
t12 = t27 * t52 + t55;
t8 = t17 * t26;
t7 = t18 * t27 + t57;
t5 = t6 * t32;
t4 = t6 * t30;
t3 = g(1) * t10 - g(2) * t12;
t2 = t11 * t35 + t12 * t37;
t1 = t11 * t37 - t12 * t35;
t14 = [0, t17, t18, t17 * t33, -t17 * sin(pkin(9)) -t18, -g(1) * (-pkin(1) * t36 + qJ(2) * t38) - g(2) * (pkin(1) * t38 + qJ(2) * t36) 0, 0, 0, 0, 0, t17 * t27, -t8, t3, -t47, t8, g(1) * t51 - g(2) * t48 + t39, t3, t8, t47, -g(1) * (-t10 * pkin(4) - t9 * qJ(5) - t51) - g(2) * (pkin(4) * t12 + qJ(5) * t11 + t48) + t39, 0, 0, 0, 0, 0, g(1) * t46 - g(2) * t2, -g(1) * t45 - g(2) * t1; 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, t5, -t4, -t7, -g(1) * (-t38 * t60 + t15) - g(2) * (-t36 * t60 + t13) - g(3) * t50, t5, -t7, t4, -g(1) * t15 - g(2) * t13 - g(3) * (t27 * t44 + t50) + (pkin(3) + t44) * t61, 0, 0, 0, 0, 0, t6 * t42, -t27 * t40 + t43 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t9 - t30 * t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t45 - t26 * t40, g(1) * t2 + g(2) * t46 + t42 * t57;];
taug_reg  = t14;
