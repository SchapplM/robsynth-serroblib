% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t21 = g(1) * t44 + g(2) * t41;
t37 = qJ(2) + pkin(9);
t33 = sin(t37);
t77 = t21 * t33;
t28 = t33 * pkin(8);
t34 = cos(t37);
t76 = -t34 * pkin(3) - t28;
t75 = -pkin(4) - pkin(5);
t40 = sin(qJ(2));
t74 = pkin(2) * t40;
t38 = -qJ(3) - pkin(7);
t72 = g(2) * t38;
t39 = sin(qJ(4));
t70 = t33 * t39;
t42 = cos(qJ(4));
t69 = t33 * t42;
t68 = t34 * t42;
t67 = t34 * t44;
t66 = t41 * t39;
t65 = t41 * t42;
t64 = t44 * t38;
t63 = t44 * t39;
t62 = t44 * t42;
t61 = qJ(5) * t39;
t60 = qJ(6) * t34;
t59 = t33 * qJ(6);
t58 = -pkin(3) - t61;
t14 = t34 * t66 + t62;
t15 = t34 * t65 - t63;
t57 = -pkin(4) * t14 + qJ(5) * t15;
t16 = t34 * t63 - t65;
t17 = t34 * t62 + t66;
t56 = -pkin(4) * t16 + qJ(5) * t17;
t55 = (pkin(8) * t34 - t74) * t41;
t54 = pkin(8) * t67 - t44 * t74;
t43 = cos(qJ(2));
t35 = t43 * pkin(2);
t53 = pkin(4) * t68 + t34 * t61 + t35 - t76;
t4 = g(1) * t14 - g(2) * t16;
t20 = g(1) * t41 - g(2) * t44;
t32 = t35 + pkin(1);
t52 = -t32 + t76;
t50 = -pkin(4) * t15 - t14 * qJ(5) - t64;
t27 = t44 * t32;
t49 = pkin(3) * t67 + pkin(4) * t17 + qJ(5) * t16 + t28 * t44 + t27;
t2 = g(1) * t16 + g(2) * t14 + g(3) * t70;
t47 = g(1) * t17 + g(2) * t15 + g(3) * t69;
t46 = -g(3) * t34 + t77;
t45 = -g(3) * t43 + t21 * t40;
t18 = qJ(5) * t69;
t13 = t20 * t33;
t8 = g(3) * t33 + t21 * t34;
t7 = t46 * t42;
t6 = t46 * t39;
t5 = g(1) * t15 - g(2) * t17;
t1 = [0, t20, t21, 0, 0, 0, 0, 0, t20 * t43, -t20 * t40, -t21, -g(1) * (-t41 * t32 - t64) - g(2) * (-t38 * t41 + t27) 0, 0, 0, 0, 0, t5, -t4, t5, t13, t4, -g(1) * t50 - g(2) * t49 + (-g(1) * t52 + t72) * t41, t5, t4, -t13, -g(1) * (-t15 * pkin(5) + t50) - g(2) * (pkin(5) * t17 - t44 * t59 + t49) + (-g(1) * (t52 + t59) + t72) * t41; 0, 0, 0, 0, 0, 0, 0, 0, t45, g(3) * t40 + t21 * t43, 0, t45 * pkin(2), 0, 0, 0, 0, 0, t7, -t6, t7, -t8, t6, -g(1) * t54 - g(2) * t55 - g(3) * t53 + (pkin(4) * t42 - t58) * t77, t7, t6, t8, -g(1) * (-t44 * t60 + t54) - g(2) * (-t41 * t60 + t55) - g(3) * (pkin(5) * t68 + t53) + (g(3) * qJ(6) + t21 * (-t42 * t75 - t58)) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t47, t2, 0, -t47, -g(1) * t56 - g(2) * t57 - g(3) * (-pkin(4) * t70 + t18) t2, -t47, 0, -g(1) * (-pkin(5) * t16 + t56) - g(2) * (-pkin(5) * t14 + t57) - g(3) * (t70 * t75 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46;];
taug_reg  = t1;
