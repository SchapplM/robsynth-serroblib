% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = sin(qJ(3));
t65 = pkin(3) + pkin(8);
t67 = t65 * t34;
t32 = qJ(1) + pkin(9);
t25 = sin(t32);
t26 = cos(t32);
t66 = -g(1) * t26 - g(2) * t25;
t37 = cos(qJ(3));
t6 = g(3) * t34 - t66 * t37;
t64 = pkin(3) * t34;
t63 = g(1) * t25;
t59 = g(3) * t37;
t29 = t37 * pkin(3);
t58 = t26 * t37;
t33 = sin(qJ(5));
t57 = t33 * t34;
t36 = cos(qJ(5));
t56 = t34 * t36;
t27 = t34 * qJ(4);
t55 = t29 + t27;
t54 = qJ(4) * t37;
t38 = cos(qJ(1));
t53 = t38 * pkin(1) + t26 * pkin(2) + t25 * pkin(7);
t52 = t37 * pkin(8) + t55;
t35 = sin(qJ(1));
t51 = -t35 * pkin(1) + t26 * pkin(7);
t50 = -pkin(2) - t27;
t49 = t26 * pkin(4) + t51;
t7 = t25 * t33 - t26 * t56;
t9 = t25 * t56 + t26 * t33;
t48 = g(1) * t9 + g(2) * t7;
t47 = pkin(3) * t58 + t26 * t27 + t53;
t15 = t25 * t54;
t17 = t26 * t54;
t46 = -g(1) * t17 - g(2) * t15;
t45 = -g(2) * t26 + t63;
t44 = g(1) * t35 - g(2) * t38;
t43 = pkin(5) * t33 - qJ(6) * t36;
t42 = t25 * pkin(4) + pkin(8) * t58 + t47;
t1 = g(1) * t7 - g(2) * t9 + t36 * t59;
t10 = -t25 * t57 + t26 * t36;
t8 = t25 * t36 + t26 * t57;
t41 = -g(1) * t8 + g(2) * t10 + t33 * t59;
t40 = (-t65 * t37 + t50) * t63;
t12 = t45 * t37;
t11 = t45 * t34;
t5 = -t34 * t66 - t59;
t4 = t6 * t36;
t3 = t6 * t33;
t2 = -g(1) * t10 - g(2) * t8;
t13 = [0, 0, 0, 0, 0, 0, t44, g(1) * t38 + g(2) * t35, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t66, 0, t44 * pkin(1), 0, 0, 0, 0, 0, 0, t12, -t11, t66, -g(1) * (-t25 * pkin(2) + t51) - g(2) * t53, 0, 0, 0, 0, 0, 0, t66, -t12, t11, -g(1) * t51 - g(2) * t47 - (t50 - t29) * t63, 0, 0, 0, 0, 0, 0, t2, t48, t12, -g(1) * t49 - g(2) * t42 - t40, 0, 0, 0, 0, 0, 0, t2, t12, -t48, -g(1) * (t10 * pkin(5) + t9 * qJ(6) + t49) - g(2) * (t8 * pkin(5) + t7 * qJ(6) + t42) - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, -g(1) * (-t26 * t64 + t17) - g(2) * (-t25 * t64 + t15) - g(3) * t55, 0, 0, 0, 0, 0, 0, -t3, -t4, t5, -g(3) * t52 - t66 * t67 + t46, 0, 0, 0, 0, 0, 0, -t3, t5, t4, -g(3) * (t43 * t34 + t52) + t46 + t66 * (t43 * t37 - t67); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t41, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, t41, -g(1) * (-t7 * pkin(5) + t8 * qJ(6)) - g(2) * (t9 * pkin(5) - t10 * qJ(6)) - (-pkin(5) * t36 - qJ(6) * t33) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
