% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = cos(qJ(1));
t67 = g(2) * t38;
t36 = sin(qJ(1));
t68 = g(1) * t36;
t70 = t67 - t68;
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t10 = -g(3) * t35 - t37 * t70;
t69 = -pkin(1) - pkin(7);
t65 = g(3) * t37;
t64 = t35 * t36;
t63 = t35 * t38;
t31 = pkin(9) + qJ(5);
t24 = sin(t31);
t62 = t36 * t24;
t25 = cos(t31);
t61 = t36 * t25;
t32 = sin(pkin(9));
t60 = t36 * t32;
t33 = cos(pkin(9));
t59 = t36 * t33;
t34 = -pkin(8) - qJ(4);
t58 = t37 * t34;
t57 = t37 * t38;
t56 = t38 * t24;
t55 = t38 * t25;
t54 = t38 * t32;
t53 = t38 * t33;
t52 = t38 * pkin(1) + t36 * qJ(2);
t51 = t37 * qJ(4);
t23 = t33 * pkin(4) + pkin(3);
t27 = t38 * qJ(2);
t50 = t23 * t63 + t34 * t57 + t27;
t49 = t38 * pkin(7) + t52;
t48 = g(2) * t49;
t5 = t35 * t62 - t55;
t7 = t35 * t56 + t61;
t47 = g(1) * t7 + g(2) * t5;
t16 = g(1) * t38 + g(2) * t36;
t45 = t35 * pkin(3) - t51;
t44 = pkin(5) * t25 + qJ(6) * t24;
t43 = pkin(4) * t54 + t23 * t64 + t36 * t58 + t49;
t42 = t23 + t44;
t41 = (-pkin(4) * t32 + t69) * t68;
t1 = g(1) * t5 - g(2) * t7 + t24 * t65;
t6 = t35 * t61 + t56;
t8 = t35 * t55 - t62;
t40 = g(1) * t6 - g(2) * t8 + t25 * t65;
t18 = t34 * t63;
t13 = t36 * t37 * t23;
t11 = t16 * t37;
t9 = g(1) * t64 - g(2) * t63 + t65;
t4 = t10 * t25;
t3 = t10 * t24;
t2 = -g(1) * t8 - g(2) * t6;
t12 = [0, 0, 0, 0, 0, 0, -t70, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t16, -g(1) * (-t36 * pkin(1) + t27) - g(2) * t52, 0, 0, 0, 0, 0, 0, -t16 * t35, -t11, -t70, -g(1) * (t69 * t36 + t27) - t48, 0, 0, 0, 0, 0, 0, -g(1) * (t35 * t53 - t60) - g(2) * (t35 * t59 + t54) -g(1) * (-t35 * t54 - t59) - g(2) * (-t35 * t60 + t53) t11, -g(1) * (pkin(3) * t63 - t38 * t51 + t27) - t48 + (-g(1) * t69 - g(2) * t45) * t36, 0, 0, 0, 0, 0, 0, t2, t47, t11, -g(1) * t50 - g(2) * t43 - t41, 0, 0, 0, 0, 0, 0, t2, t11, -t47, -g(1) * (t8 * pkin(5) + t7 * qJ(6) + t50) - g(2) * (t6 * pkin(5) + t5 * qJ(6) + t43) - t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t33, t10 * t32, -t9, g(3) * t45 + t70 * (pkin(3) * t37 + qJ(4) * t35) 0, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(1) * (-t34 * t64 + t13) - g(2) * (-t23 * t57 + t18) - g(3) * (-t35 * t23 - t58) 0, 0, 0, 0, 0, 0, -t4, -t9, -t3, -g(1) * t13 - g(2) * t18 + (g(3) * t42 + t34 * t68) * t35 + (g(3) * t34 + t42 * t67 - t44 * t68) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t40, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t40, -g(1) * (-t5 * pkin(5) + t6 * qJ(6)) - g(2) * (t7 * pkin(5) - t8 * qJ(6)) - (-pkin(5) * t24 + qJ(6) * t25) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
