% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t36 = cos(qJ(2));
t31 = qJ(2) + qJ(3);
t26 = cos(t31);
t27 = qJ(4) + t31;
t21 = sin(t27);
t22 = cos(t27);
t46 = t22 * pkin(4) + t21 * qJ(5);
t44 = pkin(3) * t26 + t46;
t59 = t36 * pkin(2) + t44;
t25 = sin(t31);
t57 = pkin(4) * t21;
t42 = -pkin(3) * t25 - t57;
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t41 = g(1) * t37 + g(2) * t35;
t5 = -g(3) * t22 + t41 * t21;
t56 = g(3) * t21;
t29 = pkin(11) + qJ(6);
t23 = sin(t29);
t54 = t35 * t23;
t24 = cos(t29);
t53 = t35 * t24;
t32 = sin(pkin(11));
t52 = t35 * t32;
t33 = cos(pkin(11));
t51 = t35 * t33;
t50 = t37 * t23;
t49 = t37 * t24;
t48 = t37 * t32;
t47 = t37 * t33;
t45 = qJ(5) * t22;
t34 = sin(qJ(2));
t43 = -t34 * pkin(2) + t42;
t40 = g(1) * t35 - g(2) * t37;
t39 = pkin(1) + t59;
t30 = -pkin(9) - pkin(8) - pkin(7);
t17 = t37 * t45;
t16 = t35 * t45;
t13 = t40 * t21;
t12 = t22 * t49 + t54;
t11 = -t22 * t50 + t53;
t10 = -t22 * t53 + t50;
t9 = t22 * t54 + t49;
t8 = g(3) * t25 + t41 * t26;
t7 = -g(3) * t26 + t41 * t25;
t6 = t41 * t22 + t56;
t4 = t5 * t33;
t3 = t5 * t32;
t2 = t5 * t24;
t1 = t5 * t23;
t14 = [0, t40, t41, 0, 0, 0, 0, 0, t40 * t36, -t40 * t34, 0, 0, 0, 0, 0, t40 * t26, -t40 * t25, 0, 0, 0, 0, 0, t40 * t22, -t13, -g(1) * (-t22 * t51 + t48) - g(2) * (t22 * t47 + t52) -g(1) * (t22 * t52 + t47) - g(2) * (-t22 * t48 + t51) t13 (g(1) * t30 - g(2) * t39) * t37 + (g(1) * t39 + g(2) * t30) * t35, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t36 + t41 * t34, g(3) * t34 + t41 * t36, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (t43 * t37 + t17) - g(2) * (t43 * t35 + t16) - g(3) * t59, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (t42 * t37 + t17) - g(2) * (t42 * t35 + t16) - g(3) * t44, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (-t37 * t57 + t17) - g(2) * (-t35 * t57 + t16) - g(3) * t46, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t23 * t56, g(1) * t12 - g(2) * t10 + t24 * t56;];
taug_reg  = t14;
