% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t19 = g(1) * t44 + g(2) * t42;
t41 = sin(qJ(2));
t87 = t19 * t41;
t40 = -pkin(8) - qJ(4);
t43 = cos(qJ(2));
t38 = sin(pkin(9));
t84 = pkin(4) * t38;
t61 = t43 * t84;
t78 = t41 * t42;
t86 = t40 * t78 + t42 * t61;
t77 = t41 * t44;
t85 = t40 * t77 + t44 * t61;
t10 = g(3) * t41 + t19 * t43;
t83 = g(1) * t42;
t79 = g(3) * t43;
t33 = t43 * pkin(2);
t37 = pkin(9) + qJ(5);
t29 = sin(t37);
t76 = t42 * t29;
t30 = cos(t37);
t75 = t42 * t30;
t74 = t42 * t38;
t39 = cos(pkin(9));
t73 = t42 * t39;
t72 = t43 * t40;
t71 = t43 * t44;
t70 = t44 * t29;
t69 = t44 * t30;
t68 = t44 * t38;
t67 = t44 * t39;
t31 = t41 * qJ(3);
t65 = t33 + t31;
t64 = t44 * pkin(1) + t42 * pkin(7);
t63 = qJ(3) * t43;
t62 = t43 * qJ(4);
t60 = t41 * t68;
t28 = t39 * pkin(4) + pkin(3);
t34 = t44 * pkin(7);
t59 = t44 * t28 + t42 * t72 + t34;
t58 = -pkin(1) - t33;
t57 = pkin(2) * t71 + t44 * t31 + t64;
t24 = t42 * t63;
t56 = -pkin(2) * t78 + t24;
t26 = t44 * t63;
t55 = -pkin(2) * t77 + t26;
t5 = -t41 * t69 + t76;
t7 = t41 * t75 + t70;
t54 = g(1) * t7 + g(2) * t5;
t53 = -g(2) * t44 + t83;
t52 = pkin(5) * t29 - qJ(6) * t30;
t51 = t41 * t84 + t65 - t72;
t50 = t58 - t31;
t6 = t41 * t70 + t75;
t8 = -t41 * t76 + t69;
t49 = -g(1) * t6 + g(2) * t8 + t29 * t79;
t1 = g(1) * t5 - g(2) * t7 + t30 * t79;
t48 = pkin(4) * t60 + t42 * t28 - t40 * t71 + t57;
t45 = ((-qJ(3) - t84) * t41 + t58) * t83;
t12 = t53 * t43;
t11 = t53 * t41;
t9 = -t79 + t87;
t4 = t10 * t30;
t3 = t10 * t29;
t2 = -g(1) * t8 - g(2) * t6;
t13 = [0, 0, 0, 0, 0, 0, t53, t19, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, -t19, -g(1) * (-t42 * pkin(1) + t34) - g(2) * t64, 0, 0, 0, 0, 0, 0, -t19, -t12, t11, -g(1) * t34 - g(2) * t57 - t50 * t83, 0, 0, 0, 0, 0, 0, -g(1) * (-t41 * t74 + t67) - g(2) * (t60 + t73) -g(1) * (-t41 * t73 - t68) - g(2) * (t41 * t67 - t74) t12, -g(1) * (t44 * pkin(3) + t34) - g(2) * (t44 * t62 + t57) + (-g(1) * (t50 - t62) - g(2) * pkin(3)) * t42, 0, 0, 0, 0, 0, 0, t2, t54, t12, -g(1) * t59 - g(2) * t48 - t45, 0, 0, 0, 0, 0, 0, t2, t12, -t54, -g(1) * (t8 * pkin(5) + t7 * qJ(6) + t59) - g(2) * (t6 * pkin(5) + t5 * qJ(6) + t48) - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -g(1) * t55 - g(2) * t56 - g(3) * t65, 0, 0, 0, 0, 0, 0, -t10 * t38, -t10 * t39, t9, -g(1) * t26 - g(2) * t24 - g(3) * (t62 + t65) + (pkin(2) + qJ(4)) * t87, 0, 0, 0, 0, 0, 0, -t3, -t4, t9, -g(1) * (t55 + t85) - g(2) * (t56 + t86) - g(3) * t51, 0, 0, 0, 0, 0, 0, -t3, t9, t4, -g(1) * (t26 + t85) - g(2) * (t24 + t86) - g(3) * (t41 * t52 + t51) + t19 * (pkin(2) * t41 - t43 * t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t49, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, t49, -g(1) * (-t5 * pkin(5) + t6 * qJ(6)) - g(2) * (t7 * pkin(5) - t8 * qJ(6)) - (-pkin(5) * t30 - qJ(6) * t29) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
