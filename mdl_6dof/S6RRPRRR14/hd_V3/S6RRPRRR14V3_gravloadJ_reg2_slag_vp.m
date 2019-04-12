% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR14V3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_gravloadJ_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = sin(qJ(1));
t43 = cos(qJ(1));
t29 = g(1) * t43 + g(2) * t38;
t36 = sin(qJ(4));
t41 = cos(qJ(4));
t56 = t43 * t41;
t42 = cos(qJ(2));
t61 = t38 * t42;
t21 = t36 * t61 + t56;
t34 = sin(qJ(6));
t39 = cos(qJ(6));
t57 = t43 * t36;
t22 = t41 * t61 - t57;
t40 = cos(qJ(5));
t35 = sin(qJ(5));
t37 = sin(qJ(2));
t65 = t37 * t35;
t5 = t22 * t40 + t38 * t65;
t75 = -t21 * t39 + t5 * t34;
t74 = t21 * t34 + t5 * t39;
t71 = g(3) * t37;
t17 = t29 * t42 + t71;
t68 = t34 * t40;
t67 = t36 * t37;
t66 = t36 * t42;
t64 = t37 * t40;
t63 = t37 * t41;
t62 = t37 * t43;
t60 = t39 * t40;
t59 = t42 * t35;
t58 = t42 * t40;
t55 = t34 * t67;
t54 = t39 * t67;
t50 = -t22 * t35 + t38 * t64;
t26 = t38 * t36 + t42 * t56;
t8 = t26 * t35 - t40 * t62;
t53 = g(1) * t50 + g(2) * t8;
t25 = -t38 * t41 + t42 * t57;
t52 = g(1) * t21 - g(2) * t25;
t51 = g(1) * t38 - g(2) * t43;
t20 = t40 * t63 - t59;
t49 = t35 * t63 + t58;
t27 = t51 * t37;
t48 = g(1) * t8 - g(2) * t50 + g(3) * t49;
t9 = t26 * t40 + t35 * t62;
t47 = g(1) * t9 + g(2) * t5 + g(3) * t20;
t46 = -g(3) * (t41 * t59 - t64) + t29 * t49;
t45 = g(1) * t25 + g(2) * t21 + g(3) * t67;
t44 = g(1) * t26 + g(2) * t22 + g(3) * t63;
t16 = -g(3) * t42 + t29 * t37;
t28 = t51 * t42;
t24 = t41 * t58 + t65;
t18 = qJ(3) * t27;
t15 = t20 * t43;
t13 = t20 * t38;
t11 = t17 * qJ(3);
t10 = t16 * t36;
t3 = t45 * t35;
t2 = t25 * t34 + t9 * t39;
t1 = t25 * t39 - t9 * t34;
t4 = [0, 0, 0, 0, 0, 0, t51, t29, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, -t29, 0, 0, 0, 0, 0, 0, 0, t28, -t29, t27, t18, 0, 0, 0, 0, 0, 0, g(1) * t22 - g(2) * t26, -t52, t27, t18, 0, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t53, t52, t18, 0, 0, 0, 0, 0, 0, g(1) * t74 - g(2) * t2, -g(1) * t75 - g(2) * t1, -t53, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t17, -t11, 0, 0, 0, 0, 0, 0, t16 * t41, -t10, -t17, -t11, 0, 0, 0, 0, 0, 0, g(1) * t15 + g(2) * t13 - g(3) * t24, -t46, t10, -t11, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t39 - t43 * t55) - g(2) * (-t13 * t39 - t38 * t55) - g(3) * (t24 * t39 + t34 * t66) -g(1) * (t15 * t34 - t43 * t54) - g(2) * (t13 * t34 - t38 * t54) - g(3) * (-t24 * t34 + t39 * t66) t46, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t44, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t40, -t3, -t44, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t25 * t60 + t26 * t34) - g(2) * (-t21 * t60 + t22 * t34) - (t34 * t41 - t36 * t60) * t71, -g(1) * (t25 * t68 + t26 * t39) - g(2) * (t21 * t68 + t22 * t39) - (t36 * t68 + t39 * t41) * t71, t3, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t47, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t39, -t48 * t34, -t47, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t75 - g(3) * (-t20 * t34 + t54) g(1) * t2 + g(2) * t74 - g(3) * (-t20 * t39 - t55) 0, 0;];
taug_reg  = t4;
