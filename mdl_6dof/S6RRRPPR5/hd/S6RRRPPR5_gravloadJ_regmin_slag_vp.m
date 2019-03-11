% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t52 = sin(qJ(2));
t53 = sin(qJ(1));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t72 = cos(pkin(6));
t65 = t56 * t72;
t22 = t53 * t52 - t55 * t65;
t45 = pkin(12) + qJ(6);
t40 = sin(t45);
t42 = cos(t45);
t23 = t52 * t65 + t53 * t55;
t46 = qJ(3) + pkin(11);
t41 = sin(t46);
t43 = cos(t46);
t48 = sin(pkin(6));
t76 = t48 * t56;
t5 = t23 * t43 - t41 * t76;
t91 = -t22 * t42 + t5 * t40;
t90 = t22 * t40 + t5 * t42;
t89 = g(3) * t48;
t66 = t53 * t72;
t25 = -t52 * t66 + t56 * t55;
t51 = sin(qJ(3));
t86 = t25 * t51;
t85 = t40 * t43;
t84 = t42 * t43;
t47 = sin(pkin(12));
t83 = t43 * t47;
t49 = cos(pkin(12));
t82 = t43 * t49;
t81 = t43 * t55;
t80 = t48 * t52;
t79 = t48 * t53;
t54 = cos(qJ(3));
t78 = t48 * t54;
t77 = t48 * t55;
t50 = -qJ(4) - pkin(9);
t75 = t50 * t52;
t39 = t54 * pkin(3) + pkin(2);
t74 = -t22 * t39 - t23 * t50;
t24 = t56 * t52 + t55 * t66;
t73 = -t24 * t39 - t25 * t50;
t71 = t51 * t80;
t70 = t51 * t79;
t69 = t53 * t78;
t33 = t51 * t76;
t68 = t54 * t76;
t67 = t23 * t54 - t33;
t64 = t72 * t54;
t63 = t56 * pkin(1) + pkin(3) * t70 + pkin(8) * t79 - t24 * t50 + t25 * t39;
t62 = g(1) * t22 - g(2) * t24;
t61 = pkin(4) * t43 + qJ(5) * t41;
t4 = t23 * t41 + t43 * t76;
t60 = t23 * t51 + t68;
t59 = -t53 * pkin(1) + pkin(3) * t33 + pkin(8) * t76 + t22 * t50 - t23 * t39;
t16 = t41 * t80 - t72 * t43;
t8 = t25 * t41 - t43 * t79;
t58 = g(1) * t8 + g(2) * t4 + g(3) * t16;
t3 = -g(1) * t24 - g(2) * t22 + g(3) * t77;
t57 = g(1) * t25 + g(2) * t23 + g(3) * t80;
t38 = pkin(3) * t64;
t30 = pkin(3) * t69;
t26 = t39 * t77;
t17 = t72 * t41 + t43 * t80;
t11 = t25 * t54 + t70;
t10 = t69 - t86;
t9 = t25 * t43 + t41 * t79;
t2 = t24 * t40 + t9 * t42;
t1 = t24 * t42 - t9 * t40;
t6 = [0, g(1) * t53 - g(2) * t56, g(1) * t56 + g(2) * t53, 0, 0, 0, 0, 0, g(1) * t23 - g(2) * t25, -t62, 0, 0, 0, 0, 0, g(1) * t67 - g(2) * t11, -g(1) * t60 - g(2) * t10, t62, -g(1) * t59 - g(2) * t63, -g(1) * (-t22 * t47 - t49 * t5) - g(2) * (t24 * t47 + t9 * t49) -g(1) * (-t22 * t49 + t47 * t5) - g(2) * (t24 * t49 - t9 * t47) g(1) * t4 - g(2) * t8, -g(1) * (-pkin(4) * t5 - qJ(5) * t4 + t59) - g(2) * (t9 * pkin(4) + t8 * qJ(5) + t63) 0, 0, 0, 0, 0, g(1) * t90 - g(2) * t2, -g(1) * t91 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t57, 0, 0, 0, 0, 0, -t3 * t54, t3 * t51, -t57, -g(1) * t73 - g(2) * t74 - g(3) * (-t48 * t75 + t26) -g(1) * (-t24 * t82 + t25 * t47) - g(2) * (-t22 * t82 + t23 * t47) - (t47 * t52 + t49 * t81) * t89, -g(1) * (t24 * t83 + t25 * t49) - g(2) * (t22 * t83 + t23 * t49) - (-t47 * t81 + t49 * t52) * t89, -t3 * t41, -g(1) * (-t24 * t61 + t73) - g(2) * (-t22 * t61 + t74) - g(3) * t26 - (t55 * t61 - t75) * t89, 0, 0, 0, 0, 0, -g(1) * (-t24 * t84 + t25 * t40) - g(2) * (-t22 * t84 + t23 * t40) - (t40 * t52 + t42 * t81) * t89, -g(1) * (t24 * t85 + t25 * t42) - g(2) * (t22 * t85 + t23 * t42) - (-t40 * t81 + t42 * t52) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t60 - g(3) * (t64 - t71) g(1) * t11 + g(2) * t67 - g(3) * (-t72 * t51 - t52 * t78) 0, -g(1) * t30 - g(3) * t38 + (g(2) * t68 + t57 * t51) * pkin(3), t58 * t49, -t58 * t47, -g(1) * t9 - g(2) * t5 - g(3) * t17, -g(1) * (-pkin(3) * t86 - t8 * pkin(4) + t9 * qJ(5) + t30) - g(2) * (-pkin(3) * t60 - t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-pkin(3) * t71 - t16 * pkin(4) + t17 * qJ(5) + t38) 0, 0, 0, 0, 0, t58 * t42, -t58 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t91 - g(3) * (-t17 * t40 - t42 * t77) g(1) * t2 + g(2) * t90 - g(3) * (-t17 * t42 + t40 * t77);];
taug_reg  = t6;
