% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t37 = sin(qJ(3));
t36 = sin(qJ(4));
t41 = cos(qJ(1));
t56 = t41 * t36;
t38 = sin(qJ(1));
t39 = cos(qJ(4));
t61 = t38 * t39;
t9 = t37 * t56 + t61;
t55 = t41 * t39;
t62 = t38 * t36;
t7 = -t37 * t62 + t55;
t69 = g(2) * t41;
t70 = g(1) * t38;
t16 = -t69 + t70;
t40 = cos(qJ(3));
t6 = -g(3) * t37 + t16 * t40;
t71 = pkin(4) * t36;
t67 = g(3) * t40;
t66 = t37 * t38;
t65 = t37 * t41;
t34 = qJ(4) + pkin(9);
t27 = sin(t34);
t64 = t38 * t27;
t28 = cos(t34);
t63 = t38 * t28;
t35 = -qJ(5) - pkin(8);
t60 = t40 * t35;
t59 = t40 * t41;
t58 = t41 * t27;
t57 = t41 * t28;
t54 = t9 * pkin(4);
t53 = t41 * pkin(1) + t38 * qJ(2);
t52 = t36 * t67;
t26 = t39 * pkin(4) + pkin(3);
t30 = t41 * qJ(2);
t49 = t26 * t65 + t35 * t59 + t30;
t17 = g(1) * t41 + g(2) * t38;
t48 = pkin(5) * t28 + qJ(6) * t27;
t47 = t7 * pkin(4);
t46 = pkin(4) * t56 + t41 * pkin(7) + t26 * t66 + t38 * t60 + t53;
t45 = t26 + t48;
t44 = (-pkin(1) - pkin(7) - t71) * t70;
t1 = t37 * t64 - t57;
t3 = t37 * t58 + t63;
t43 = g(1) * t1 - g(2) * t3 + t27 * t67;
t19 = t35 * t65;
t13 = t38 * t40 * t26;
t11 = t17 * t40;
t10 = t37 * t55 - t62;
t8 = t37 * t61 + t56;
t5 = g(1) * t66 - g(2) * t65 + t67;
t4 = t37 * t57 - t64;
t2 = t37 * t63 + t58;
t12 = [0, t16, t17, -t16, -t17, -g(1) * (-t38 * pkin(1) + t30) - g(2) * t53, 0, 0, 0, 0, 0, -t17 * t37, -t11, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t11, -g(1) * t49 - g(2) * t46 - t44, -g(1) * t4 - g(2) * t2, t11, -g(1) * t3 - g(2) * t1, -g(1) * (t4 * pkin(5) + t3 * qJ(6) + t49) - g(2) * (t2 * pkin(5) + t1 * qJ(6) + t46) - t44; 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t6 * t39, t6 * t36, -t5, -g(1) * (-t35 * t66 + t13) - g(2) * (-t26 * t59 + t19) - g(3) * (-t37 * t26 - t60) -t6 * t28, -t5, -t6 * t27, -g(1) * t13 - g(2) * t19 + (g(3) * t45 + t35 * t70) * t37 + (g(3) * t35 + t45 * t69 - t48 * t70) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t52, g(1) * t8 - g(2) * t10 + t39 * t67, 0, pkin(4) * t52 - g(1) * t47 - g(2) * t54, t43, 0, -g(1) * t2 + g(2) * t4 - t28 * t67, -g(1) * (-t1 * pkin(5) + t2 * qJ(6) + t47) - g(2) * (t3 * pkin(5) - t4 * qJ(6) + t54) - (-pkin(5) * t27 + qJ(6) * t28 - t71) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43;];
taug_reg  = t12;
