% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = qJ(1) + pkin(9);
t33 = sin(t38);
t34 = cos(t38);
t19 = g(1) * t34 + g(2) * t33;
t40 = sin(qJ(3));
t75 = t19 * t40;
t74 = -pkin(4) - pkin(5);
t73 = g(1) * t33;
t35 = t40 * pkin(8);
t43 = cos(qJ(3));
t36 = t43 * pkin(3);
t70 = t33 * t40;
t69 = t34 * t40;
t68 = t34 * t43;
t39 = sin(qJ(4));
t67 = t39 * t40;
t66 = t39 * t43;
t42 = cos(qJ(4));
t65 = t40 * t42;
t64 = t42 * t43;
t63 = t36 + t35;
t62 = qJ(5) * t39;
t61 = qJ(6) * t43;
t44 = cos(qJ(1));
t60 = t44 * pkin(1) + t34 * pkin(2) + t33 * pkin(7);
t59 = -pkin(2) - t36;
t41 = sin(qJ(1));
t58 = -t41 * pkin(1) + t34 * pkin(7);
t57 = -pkin(3) - t62;
t14 = t33 * t66 + t34 * t42;
t15 = t33 * t64 - t34 * t39;
t56 = -t14 * pkin(4) + t15 * qJ(5);
t16 = -t33 * t42 + t34 * t66;
t17 = t33 * t39 + t34 * t64;
t55 = -t16 * pkin(4) + t17 * qJ(5);
t54 = pkin(4) * t64 + t43 * t62 + t63;
t53 = pkin(3) * t68 + pkin(8) * t69 + t60;
t4 = g(1) * t14 - g(2) * t16;
t52 = -g(2) * t34 + t73;
t51 = g(1) * t41 - g(2) * t44;
t49 = -t15 * pkin(4) - t14 * qJ(5) + t58;
t47 = (t59 - t35) * t73;
t2 = g(1) * t16 + g(2) * t14 + g(3) * t67;
t46 = g(1) * t17 + g(2) * t15 + g(3) * t65;
t45 = t17 * pkin(4) + t16 * qJ(5) + t53;
t8 = -g(3) * t43 + t75;
t27 = qJ(5) * t65;
t22 = pkin(8) * t68;
t20 = t33 * t43 * pkin(8);
t18 = g(1) * t70 - g(2) * t69;
t9 = g(3) * t40 + t19 * t43;
t7 = t8 * t42;
t6 = t8 * t39;
t5 = g(1) * t15 - g(2) * t17;
t1 = [0, 0, 0, 0, 0, 0, t51, g(1) * t44 + g(2) * t41, 0, 0, 0, 0, 0, 0, 0, 0, t52, t19, 0, t51 * pkin(1), 0, 0, 0, 0, 0, 0, t52 * t43, -t18, -t19, -g(1) * (-t33 * pkin(2) + t58) - g(2) * t60, 0, 0, 0, 0, 0, 0, t5, -t4, t18, -g(1) * t58 - g(2) * t53 - t47, 0, 0, 0, 0, 0, 0, t5, t18, t4, -g(1) * t49 - g(2) * t45 - t47, 0, 0, 0, 0, 0, 0, t5, t4, -t18, -g(1) * (-t15 * pkin(5) + t49) - g(2) * (t17 * pkin(5) - qJ(6) * t69 + t45) - ((-pkin(8) + qJ(6)) * t40 + t59) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (-pkin(3) * t69 + t22) - g(2) * (-pkin(3) * t70 + t20) - g(3) * t63, 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * t22 - g(2) * t20 - g(3) * t54 + (pkin(4) * t42 - t57) * t75, 0, 0, 0, 0, 0, 0, t7, t6, t9, -g(1) * (-t34 * t61 + t22) - g(2) * (-t33 * t61 + t20) - g(3) * (pkin(5) * t64 + t54) + (g(3) * qJ(6) + t19 * (-t74 * t42 - t57)) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t46, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t46, -g(1) * t55 - g(2) * t56 - g(3) * (-pkin(4) * t67 + t27) 0, 0, 0, 0, 0, 0, t2, -t46, 0, -g(1) * (-t16 * pkin(5) + t55) - g(2) * (-t14 * pkin(5) + t56) - g(3) * (t74 * t67 + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t1;
