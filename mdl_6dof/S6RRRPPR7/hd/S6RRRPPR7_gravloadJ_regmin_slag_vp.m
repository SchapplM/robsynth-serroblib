% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t46 = sin(qJ(1));
t47 = cos(qJ(3));
t48 = cos(qJ(2));
t44 = sin(qJ(3));
t49 = cos(qJ(1));
t77 = t49 * t44;
t21 = -t46 * t47 + t48 * t77;
t76 = t49 * t47;
t22 = t46 * t44 + t48 * t76;
t41 = pkin(10) + qJ(6);
t34 = sin(t41);
t35 = cos(t41);
t1 = t21 * t35 - t22 * t34;
t58 = t34 * t47 - t35 * t44;
t45 = sin(qJ(2));
t86 = g(3) * t45;
t80 = t46 * t48;
t19 = t44 * t80 + t76;
t79 = t47 * t48;
t20 = t46 * t79 - t77;
t99 = t19 * t35 - t20 * t34;
t100 = g(1) * t1 + g(2) * t99 - t58 * t86;
t63 = g(1) * t49 + g(2) * t46;
t97 = t63 * t45;
t95 = g(3) * t48 - t97;
t2 = t21 * t34 + t22 * t35;
t57 = t34 * t44 + t35 * t47;
t61 = t19 * t34 + t20 * t35;
t96 = g(1) * t2 + g(2) * t61 + t57 * t86;
t90 = -pkin(3) - pkin(4);
t89 = g(1) * t46;
t36 = t45 * pkin(8);
t38 = t48 * pkin(2);
t83 = t44 * t45;
t82 = t45 * t47;
t81 = t45 * t49;
t78 = t48 * t49;
t75 = qJ(4) * t44;
t74 = qJ(5) * t49;
t73 = -pkin(1) - t38;
t72 = -pkin(2) - t75;
t42 = sin(pkin(10));
t43 = cos(pkin(10));
t70 = -t19 * t43 + t20 * t42;
t69 = t21 * t42 + t22 * t43;
t68 = -t19 * pkin(3) + t20 * qJ(4);
t67 = -t21 * pkin(3) + t22 * qJ(4);
t66 = pkin(3) * t79 + t48 * t75 + t36 + t38;
t65 = -t20 * pkin(3) + t49 * pkin(7) - t19 * qJ(4);
t64 = g(1) * t19 - g(2) * t21;
t62 = -g(2) * t49 + t89;
t60 = t19 * t42 + t20 * t43;
t59 = t21 * t43 - t22 * t42;
t56 = t42 * t47 - t43 * t44;
t55 = t42 * t44 + t43 * t47;
t53 = g(3) * t55;
t51 = t49 * pkin(1) + pkin(2) * t78 + t22 * pkin(3) + t46 * pkin(7) + pkin(8) * t81 + t21 * qJ(4);
t4 = g(1) * t21 + g(2) * t19 + g(3) * t83;
t50 = g(1) * t22 + g(2) * t20 + g(3) * t82;
t29 = pkin(8) * t78;
t26 = pkin(8) * t80;
t24 = qJ(4) * t82;
t23 = -g(2) * t81 + t45 * t89;
t14 = t63 * t48 + t86;
t7 = t95 * t47;
t6 = t95 * t44;
t5 = g(1) * t20 - g(2) * t22;
t3 = [0, t62, t63, 0, 0, 0, 0, 0, t62 * t48, -t23, 0, 0, 0, 0, 0, t5, -t64, t5, t23, t64, -g(1) * t65 - g(2) * t51 - (t73 - t36) * t89, g(1) * t60 - g(2) * t69, -g(1) * t70 - g(2) * t59, -t23, -g(1) * (-t20 * pkin(4) + t65) - g(2) * (t22 * pkin(4) - t45 * t74 + t51) - ((-pkin(8) + qJ(5)) * t45 + t73) * t89, 0, 0, 0, 0, 0, g(1) * t61 - g(2) * t2, g(1) * t99 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t95, t14, 0, 0, 0, 0, 0, -t7, t6, -t7, -t14, -t6, -g(1) * t29 - g(2) * t26 - g(3) * t66 + (pkin(3) * t47 - t72) * t97, -t48 * t53 + t55 * t97, t95 * t56, t14, -g(1) * (-t48 * t74 + t29) - g(2) * (-qJ(5) * t80 + t26) - g(3) * (pkin(4) * t79 + t66) + (g(3) * qJ(5) + t63 * (-t90 * t47 - t72)) * t45, 0, 0, 0, 0, 0, -t95 * t57, t95 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t50, t4, 0, -t50, -g(1) * t67 - g(2) * t68 - g(3) * (-pkin(3) * t83 + t24) g(1) * t59 - g(2) * t70 - t56 * t86, -g(1) * t69 - g(2) * t60 - t45 * t53, 0, -g(1) * (-t21 * pkin(4) + t67) - g(2) * (-t19 * pkin(4) + t68) - g(3) * (t90 * t83 + t24) 0, 0, 0, 0, 0, t100, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t96;];
taug_reg  = t3;
