% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t40 = sin(qJ(2));
t41 = sin(qJ(1));
t44 = cos(qJ(2));
t45 = cos(qJ(1));
t59 = cos(pkin(6));
t54 = t45 * t59;
t20 = t41 * t40 - t44 * t54;
t35 = qJ(4) + pkin(11);
t32 = sin(t35);
t33 = cos(t35);
t36 = sin(pkin(6));
t63 = t36 * t45;
t10 = -t20 * t32 + t33 * t63;
t21 = t40 * t54 + t41 * t44;
t38 = sin(qJ(6));
t42 = cos(qJ(6));
t76 = t10 * t38 + t21 * t42;
t75 = t10 * t42 - t21 * t38;
t55 = t41 * t59;
t22 = t45 * t40 + t44 * t55;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t65 = t36 * t41;
t12 = t22 * t43 - t39 * t65;
t50 = t20 * t43 + t39 * t63;
t64 = t36 * t44;
t74 = -g(3) * (-t59 * t39 - t43 * t64) - g(2) * t50 - g(1) * t12;
t73 = pkin(4) * t39;
t71 = g(3) * t36;
t68 = t32 * t38;
t67 = t32 * t42;
t66 = t36 * t40;
t62 = t38 * t40;
t61 = t40 * t42;
t60 = pkin(2) * t64 + qJ(3) * t66;
t23 = -t40 * t55 + t45 * t44;
t58 = t45 * pkin(1) + t23 * pkin(2) + pkin(8) * t65;
t57 = qJ(3) + t73;
t56 = t20 * t39 - t43 * t63;
t53 = -t41 * pkin(1) - t21 * pkin(2) + pkin(8) * t63;
t52 = g(1) * t20 - g(2) * t22;
t6 = g(1) * t21 - g(2) * t23;
t51 = g(1) * t45 + g(2) * t41;
t49 = g(1) * (t22 * t33 - t32 * t65) + g(2) * (t20 * t33 + t32 * t63) + g(3) * (-t59 * t32 - t33 * t64);
t3 = -g(1) * t22 - g(2) * t20 + g(3) * t64;
t47 = g(1) * t23 + g(2) * t21 + g(3) * t66;
t37 = -qJ(5) - pkin(9);
t31 = t43 * pkin(4) + pkin(3);
t18 = t22 * pkin(2);
t16 = t20 * pkin(2);
t15 = -t32 * t64 + t59 * t33;
t13 = t22 * t39 + t43 * t65;
t8 = t22 * t32 + t33 * t65;
t2 = t23 * t38 + t8 * t42;
t1 = t23 * t42 - t8 * t38;
t4 = [0, g(1) * t41 - g(2) * t45, t51, 0, 0, 0, 0, 0, t6, -t52, -t51 * t36, -t6, t52, -g(1) * (-t20 * qJ(3) + t53) - g(2) * (t22 * qJ(3) + t58) 0, 0, 0, 0, 0, g(1) * t56 - g(2) * t13, g(1) * t50 - g(2) * t12, t6, -g(1) * (-t57 * t20 + t21 * t37 + t31 * t63 + t53) - g(2) * (t57 * t22 - t23 * t37 + t31 * t65 + t58) 0, 0, 0, 0, 0, -g(1) * t75 - g(2) * t2, g(1) * t76 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t47, 0, t3, -t47, -g(1) * (t23 * qJ(3) - t18) - g(2) * (t21 * qJ(3) - t16) - g(3) * t60, 0, 0, 0, 0, 0, -t47 * t39, -t47 * t43, -t3, -g(1) * (t22 * t37 + t57 * t23 - t18) - g(2) * (t20 * t37 + t57 * t21 - t16) - g(3) * ((-t37 * t44 + t40 * t73) * t36 + t60) 0, 0, 0, 0, 0, -g(1) * (-t22 * t38 + t23 * t67) - g(2) * (-t20 * t38 + t21 * t67) - (t32 * t61 + t38 * t44) * t71, -g(1) * (-t22 * t42 - t23 * t68) - g(2) * (-t20 * t42 - t21 * t68) - (-t32 * t62 + t42 * t44) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, g(1) * t13 + g(2) * t56 - g(3) * (t39 * t64 - t59 * t43) 0, t74 * pkin(4), 0, 0, 0, 0, 0, -t49 * t42, t49 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t76 - g(3) * (-t15 * t38 + t36 * t61) g(1) * t2 - g(2) * t75 - g(3) * (-t15 * t42 - t36 * t62);];
taug_reg  = t4;
