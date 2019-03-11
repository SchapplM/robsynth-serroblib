% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t50 = qJ(4) + qJ(5);
t29 = sin(t50);
t32 = sin(qJ(2));
t36 = cos(qJ(2));
t49 = cos(t50);
t38 = t32 * t29 + t36 * t49;
t33 = sin(qJ(1));
t45 = t32 * t49;
t57 = -t36 * t29 + t45;
t7 = t57 * t33;
t37 = cos(qJ(1));
t51 = t36 * t37;
t9 = t29 * t51 - t37 * t45;
t3 = g(1) * t9 - g(2) * t7 + g(3) * t38;
t30 = sin(qJ(6));
t1 = t3 * t30;
t34 = cos(qJ(6));
t2 = t3 * t34;
t25 = g(1) * t37 + g(2) * t33;
t54 = g(3) * t57;
t35 = cos(qJ(4));
t53 = t32 * t35;
t48 = g(1) * t33 - g(2) * t37;
t8 = t38 * t33;
t47 = t8 * t30 - t37 * t34;
t46 = t37 * t30 + t8 * t34;
t44 = t36 * pkin(2) + t32 * qJ(3);
t31 = sin(qJ(4));
t42 = t36 * t31 - t53;
t21 = t32 * t31 + t36 * t35;
t41 = pkin(1) + t44;
t10 = t38 * t37;
t4 = g(1) * t10 + g(2) * t8 + t54;
t13 = t42 * t33;
t15 = t31 * t51 - t37 * t53;
t40 = g(1) * t15 + g(2) * t13 + g(3) * t21;
t14 = t21 * t33;
t16 = t21 * t37;
t39 = g(1) * t16 + g(2) * t14 - g(3) * t42;
t20 = t48 * t36;
t19 = t48 * t32;
t12 = g(3) * t32 + t25 * t36;
t11 = -g(3) * t36 + t25 * t32;
t6 = t10 * t34 - t33 * t30;
t5 = -t10 * t30 - t33 * t34;
t17 = [0, t48, t25, 0, 0, 0, 0, 0, t20, -t19, t20, -t25, t19 (-g(1) * pkin(7) - g(2) * t41) * t37 + (-g(2) * pkin(7) + g(1) * t41) * t33, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t16, -g(1) * t13 + g(2) * t15, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t10, g(1) * t7 + g(2) * t9, 0, 0, 0, 0, 0, g(1) * t46 - g(2) * t6, -g(1) * t47 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, t11, 0, -t12, -g(3) * t44 + t25 * (pkin(2) * t32 - qJ(3) * t36) 0, 0, 0, 0, 0, -t40, -t39, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t47 + t30 * t54, g(1) * t6 + g(2) * t46 + t34 * t54;];
taug_reg  = t17;
