% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t33 = cos(pkin(10));
t23 = t33 * pkin(4) + pkin(3);
t34 = -pkin(8) - qJ(4);
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t50 = t37 * t23 - t35 * t34;
t31 = qJ(1) + pkin(9);
t25 = sin(t31);
t27 = cos(t31);
t13 = g(1) * t27 + g(2) * t25;
t9 = -g(3) * t37 + t13 * t35;
t65 = g(1) * t25;
t62 = g(3) * t35;
t32 = sin(pkin(10));
t60 = t25 * t32;
t59 = t25 * t37;
t58 = t27 * t32;
t57 = t27 * t37;
t56 = t32 * t37;
t55 = t33 * t37;
t54 = t34 * t37;
t38 = cos(qJ(1));
t52 = t38 * pkin(1) + t27 * pkin(2) + t25 * pkin(7);
t36 = sin(qJ(1));
t51 = -t36 * pkin(1) + t27 * pkin(7);
t30 = pkin(10) + qJ(5);
t24 = sin(t30);
t26 = cos(t30);
t5 = t24 * t59 + t27 * t26;
t7 = t24 * t57 - t25 * t26;
t49 = g(1) * t5 - g(2) * t7;
t48 = -g(2) * t27 + t65;
t47 = g(1) * t36 - g(2) * t38;
t46 = t37 * pkin(3) + t35 * qJ(4);
t44 = pkin(5) * t26 + qJ(6) * t24;
t1 = g(1) * t7 + g(2) * t5 + t24 * t62;
t6 = -t27 * t24 + t26 * t59;
t8 = t25 * t24 + t26 * t57;
t42 = g(1) * t8 + g(2) * t6 + t26 * t62;
t41 = pkin(4) * t60 + t27 * t50 + t52;
t39 = pkin(4) * t58 + t51 + (-pkin(2) - t50) * t25;
t11 = t48 * t35;
t10 = t13 * t37 + t62;
t4 = t9 * t26;
t3 = t9 * t24;
t2 = g(1) * t6 - g(2) * t8;
t12 = [0, 0, 0, 0, 0, 0, t47, g(1) * t38 + g(2) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t48, t13, 0, t47 * pkin(1), 0, 0, 0, 0, 0, 0, t48 * t37, -t11, -t13, -g(1) * (-t25 * pkin(2) + t51) - g(2) * t52, 0, 0, 0, 0, 0, 0, -g(1) * (-t25 * t55 + t58) - g(2) * (t27 * t55 + t60) -g(1) * (t25 * t56 + t27 * t33) - g(2) * (t25 * t33 - t27 * t56) t11, -g(1) * t51 - g(2) * (t46 * t27 + t52) - (-pkin(2) - t46) * t65, 0, 0, 0, 0, 0, 0, t2, -t49, t11, -g(1) * t39 - g(2) * t41, 0, 0, 0, 0, 0, 0, t2, t11, t49, -g(1) * (-t6 * pkin(5) - t5 * qJ(6) + t39) - g(2) * (t8 * pkin(5) + t7 * qJ(6) + t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t33, -t9 * t32, -t10, -g(3) * t46 + t13 * (pkin(3) * t35 - qJ(4) * t37) 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(3) * t50 + t13 * (t23 * t35 + t54) 0, 0, 0, 0, 0, 0, t4, -t10, t3, -g(3) * (t44 * t37 + t50) + t13 * (t54 - (-t23 - t44) * t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t42, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t42, -g(1) * (-t7 * pkin(5) + t8 * qJ(6)) - g(2) * (-t5 * pkin(5) + t6 * qJ(6)) - (-pkin(5) * t24 + qJ(6) * t26) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
