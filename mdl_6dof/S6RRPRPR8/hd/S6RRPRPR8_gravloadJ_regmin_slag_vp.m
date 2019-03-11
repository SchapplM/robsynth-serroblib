% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t44 = g(1) * t35 + g(2) * t32;
t72 = -g(3) * t34 + t44 * t31;
t26 = pkin(10) + qJ(4);
t21 = cos(t26);
t20 = sin(t26);
t51 = t35 * t20;
t12 = -t32 * t21 + t34 * t51;
t50 = t35 * t21;
t13 = t32 * t20 + t34 * t50;
t30 = sin(qJ(6));
t33 = cos(qJ(6));
t2 = t12 * t33 - t13 * t30;
t39 = t20 * t33 - t21 * t30;
t53 = t32 * t34;
t10 = t20 * t53 + t50;
t11 = t21 * t53 - t51;
t46 = -t10 * t33 + t11 * t30;
t57 = g(3) * t31;
t71 = -g(1) * t2 + g(2) * t46 - t39 * t57;
t28 = cos(pkin(10));
t19 = t28 * pkin(3) + pkin(2);
t29 = -pkin(8) - qJ(3);
t68 = t34 * t19 - t31 * t29;
t3 = t12 * t30 + t13 * t33;
t38 = t20 * t30 + t21 * t33;
t40 = t10 * t30 + t11 * t33;
t66 = g(1) * t3 + g(2) * t40 + t38 * t57;
t60 = g(1) * t32;
t27 = sin(pkin(10));
t49 = t35 * t27;
t48 = t35 * t28;
t47 = t35 * pkin(1) + t32 * pkin(7);
t45 = g(1) * t10 - g(2) * t12;
t43 = -g(2) * t35 + t60;
t42 = t34 * pkin(2) + t31 * qJ(3);
t37 = pkin(4) * t21 + qJ(5) * t20 + t19;
t1 = g(1) * t12 + g(2) * t10 + t20 * t57;
t36 = g(1) * t13 + g(2) * t11 + t21 * t57;
t23 = t35 * pkin(7);
t16 = t43 * t31;
t15 = t44 * t34 + t57;
t6 = t72 * t21;
t5 = t72 * t20;
t4 = g(1) * t11 - g(2) * t13;
t7 = [0, t43, t44, 0, 0, 0, 0, 0, t43 * t34, -t16, -g(1) * (-t28 * t53 + t49) - g(2) * (t32 * t27 + t34 * t48) -g(1) * (t27 * t53 + t48) - g(2) * (t32 * t28 - t34 * t49) t16, -g(1) * t23 - g(2) * (t42 * t35 + t47) - (-pkin(1) - t42) * t60, 0, 0, 0, 0, 0, t4, -t45, t4, t16, t45, -g(1) * (pkin(3) * t49 - t11 * pkin(4) - t10 * qJ(5) + t23) - g(2) * (t13 * pkin(4) + t12 * qJ(5) + t68 * t35 + t47) + (-g(1) * (-pkin(1) - t68) - g(2) * pkin(3) * t27) * t32, 0, 0, 0, 0, 0, g(1) * t40 - g(2) * t3, -g(1) * t46 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, t72, t15, t72 * t28, -t72 * t27, -t15, -g(3) * t42 + t44 * (pkin(2) * t31 - qJ(3) * t34) 0, 0, 0, 0, 0, t6, -t5, t6, -t15, t5 (-g(3) * t37 + t44 * t29) * t34 + (g(3) * t29 + t44 * t37) * t31, 0, 0, 0, 0, 0, t72 * t38, t72 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t36, t1, 0, -t36, -g(1) * (-t12 * pkin(4) + t13 * qJ(5)) - g(2) * (-t10 * pkin(4) + t11 * qJ(5)) - (-pkin(4) * t20 + qJ(5) * t21) * t57, 0, 0, 0, 0, 0, -t71, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t66;];
taug_reg  = t7;
