% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPP5
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
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t43 = sin(qJ(1));
t45 = cos(qJ(1));
t23 = g(1) * t45 + g(2) * t43;
t38 = pkin(9) + qJ(3);
t35 = sin(t38);
t74 = t23 * t35;
t36 = cos(t38);
t73 = -t36 * pkin(3) - t35 * pkin(8);
t72 = -pkin(4) - pkin(5);
t71 = g(1) * t43;
t41 = -pkin(7) - qJ(2);
t69 = g(2) * t41;
t42 = sin(qJ(4));
t67 = t35 * t42;
t44 = cos(qJ(4));
t66 = t35 * t44;
t65 = t35 * t45;
t64 = t36 * t44;
t63 = t36 * t45;
t62 = t43 * t42;
t61 = t43 * t44;
t60 = t45 * t42;
t59 = t45 * t44;
t58 = qJ(5) * t42;
t57 = qJ(6) * t36;
t56 = t35 * qJ(6);
t55 = -pkin(3) - t58;
t15 = t36 * t62 + t59;
t16 = t36 * t61 - t60;
t54 = -t15 * pkin(4) + t16 * qJ(5);
t17 = t36 * t60 - t61;
t18 = t36 * t59 + t62;
t53 = -t17 * pkin(4) + t18 * qJ(5);
t52 = pkin(4) * t64 + t36 * t58 - t73;
t4 = g(1) * t15 - g(2) * t17;
t22 = -g(2) * t45 + t71;
t40 = cos(pkin(9));
t32 = t40 * pkin(2) + pkin(1);
t51 = -t32 + t73;
t49 = -t16 * pkin(4) - t15 * qJ(5) - t45 * t41;
t48 = pkin(3) * t63 + t18 * pkin(4) + pkin(8) * t65 + t17 * qJ(5) + t45 * t32;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t67;
t46 = g(1) * t18 + g(2) * t16 + g(3) * t66;
t8 = -g(3) * t36 + t74;
t27 = pkin(8) * t63;
t24 = t43 * t36 * pkin(8);
t19 = qJ(5) * t66;
t14 = -g(2) * t65 + t35 * t71;
t9 = g(3) * t35 + t23 * t36;
t7 = t8 * t44;
t6 = t8 * t42;
t5 = g(1) * t16 - g(2) * t18;
t1 = [0, t22, t23, t22 * t40, -t22 * sin(pkin(9)) -t23, -g(1) * (-t43 * pkin(1) + t45 * qJ(2)) - g(2) * (t45 * pkin(1) + t43 * qJ(2)) 0, 0, 0, 0, 0, t22 * t36, -t14, 0, 0, 0, 0, 0, t5, -t4, t5, t14, t4, -g(1) * t49 - g(2) * t48 + (-g(1) * t51 + t69) * t43, t5, t4, -t14, -g(1) * (-t16 * pkin(5) + t49) - g(2) * (t18 * pkin(5) - t45 * t56 + t48) + (-g(1) * (t51 + t56) + t69) * t43; 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t6, t7, -t9, t6, -g(1) * t27 - g(2) * t24 - g(3) * t52 + (pkin(4) * t44 - t55) * t74, t7, t6, t9, -g(1) * (-t45 * t57 + t27) - g(2) * (-t43 * t57 + t24) - g(3) * (pkin(5) * t64 + t52) + (g(3) * qJ(6) + t23 * (-t72 * t44 - t55)) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t46, t2, 0, -t46, -g(1) * t53 - g(2) * t54 - g(3) * (-pkin(4) * t67 + t19) t2, -t46, 0, -g(1) * (-t17 * pkin(5) + t53) - g(2) * (-t15 * pkin(5) + t54) - g(3) * (t72 * t67 + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t1;
