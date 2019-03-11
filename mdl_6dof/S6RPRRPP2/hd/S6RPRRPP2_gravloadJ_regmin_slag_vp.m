% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t39 = sin(qJ(3));
t37 = qJ(1) + pkin(9);
t32 = sin(t37);
t33 = cos(t37);
t50 = g(1) * t33 + g(2) * t32;
t68 = t50 * t39;
t67 = -pkin(4) - pkin(5);
t66 = g(1) * t32;
t34 = t39 * pkin(8);
t42 = cos(qJ(3));
t35 = t42 * pkin(3);
t63 = t33 * t39;
t62 = t33 * t42;
t38 = sin(qJ(4));
t61 = t38 * t39;
t60 = t38 * t42;
t41 = cos(qJ(4));
t59 = t39 * t41;
t58 = t41 * t42;
t57 = qJ(5) * t38;
t56 = qJ(6) * t42;
t55 = -pkin(2) - t35;
t54 = -pkin(3) - t57;
t14 = t32 * t60 + t33 * t41;
t15 = t32 * t58 - t33 * t38;
t53 = -t14 * pkin(4) + t15 * qJ(5);
t16 = -t32 * t41 + t33 * t60;
t17 = t32 * t38 + t33 * t58;
t52 = -t16 * pkin(4) + t17 * qJ(5);
t51 = pkin(4) * t58 + t42 * t57 + t34 + t35;
t4 = g(1) * t14 - g(2) * t16;
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t49 = g(1) * t40 - g(2) * t43;
t47 = -t40 * pkin(1) - t15 * pkin(4) + t33 * pkin(7) - t14 * qJ(5);
t2 = g(1) * t16 + g(2) * t14 + g(3) * t61;
t45 = g(1) * t17 + g(2) * t15 + g(3) * t59;
t44 = t43 * pkin(1) + t33 * pkin(2) + pkin(3) * t62 + t17 * pkin(4) + t32 * pkin(7) + pkin(8) * t63 + t16 * qJ(5);
t8 = -g(3) * t42 + t68;
t26 = qJ(5) * t59;
t21 = pkin(8) * t62;
t19 = t32 * t42 * pkin(8);
t18 = -g(2) * t63 + t39 * t66;
t9 = g(3) * t39 + t50 * t42;
t7 = t8 * t41;
t6 = t8 * t38;
t5 = g(1) * t15 - g(2) * t17;
t1 = [0, t49, g(1) * t43 + g(2) * t40, t49 * pkin(1), 0, 0, 0, 0, 0 (-g(2) * t33 + t66) * t42, -t18, 0, 0, 0, 0, 0, t5, -t4, t5, t18, t4, -g(1) * t47 - g(2) * t44 - (t55 - t34) * t66, t5, t4, -t18, -g(1) * (-t15 * pkin(5) + t47) - g(2) * (t17 * pkin(5) - qJ(6) * t63 + t44) - ((-pkin(8) + qJ(6)) * t39 + t55) * t66; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t6, t7, -t9, t6, -g(1) * t21 - g(2) * t19 - g(3) * t51 + (pkin(4) * t41 - t54) * t68, t7, t6, t9, -g(1) * (-t33 * t56 + t21) - g(2) * (-t32 * t56 + t19) - g(3) * (pkin(5) * t58 + t51) + (g(3) * qJ(6) + t50 * (-t67 * t41 - t54)) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t45, t2, 0, -t45, -g(1) * t52 - g(2) * t53 - g(3) * (-pkin(4) * t61 + t26) t2, -t45, 0, -g(1) * (-t16 * pkin(5) + t52) - g(2) * (-t14 * pkin(5) + t53) - g(3) * (t67 * t61 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t1;
