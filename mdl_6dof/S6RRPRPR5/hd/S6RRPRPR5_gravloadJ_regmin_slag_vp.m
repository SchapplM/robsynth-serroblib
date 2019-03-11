% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t44 = sin(pkin(11));
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t70 = cos(pkin(11));
t29 = -t52 * t44 - t49 * t70;
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t47 = cos(pkin(6));
t58 = -t49 * t44 + t52 * t70;
t54 = t58 * t47;
t11 = t50 * t29 + t53 * t54;
t42 = pkin(12) + qJ(6);
t40 = sin(t42);
t41 = cos(t42);
t71 = t29 * t47;
t10 = -t50 * t58 + t53 * t71;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t45 = sin(pkin(6));
t81 = t45 * t53;
t5 = -t10 * t51 - t48 * t81;
t96 = t11 * t41 + t5 * t40;
t95 = -t11 * t40 + t5 * t41;
t73 = t53 * t49;
t76 = t50 * t52;
t26 = -t47 * t76 - t73;
t82 = t45 * t52;
t94 = -g(1) * t26 - g(3) * t82;
t14 = t53 * t29 - t50 * t54;
t21 = t58 * t45;
t93 = -g(1) * t14 - g(2) * t11 - g(3) * t21;
t15 = t50 * t71 + t53 * t58;
t86 = t40 * t51;
t85 = t41 * t51;
t43 = sin(pkin(12));
t84 = t43 * t51;
t83 = t45 * t50;
t46 = cos(pkin(12));
t80 = t46 * t51;
t77 = t50 * t49;
t72 = t53 * t52;
t67 = t47 * t72;
t23 = t47 * t49 * pkin(2) + (-pkin(8) - qJ(3)) * t45;
t39 = t52 * pkin(2) + pkin(1);
t66 = -t50 * t23 + t53 * t39;
t4 = -t10 * t48 + t51 * t81;
t8 = t15 * t48 - t51 * t83;
t63 = -g(1) * t4 + g(2) * t8;
t62 = g(1) * t53 + g(2) * t50;
t61 = g(1) * t50 - g(2) * t53;
t60 = -t53 * t23 - t50 * t39;
t22 = t29 * t45;
t16 = -t22 * t48 - t47 * t51;
t57 = g(1) * t8 + g(2) * t4 + g(3) * t16;
t17 = -t22 * t51 + t47 * t48;
t9 = t15 * t51 + t48 * t83;
t56 = g(1) * t9 + g(2) * t5 + g(3) * t17;
t30 = pkin(2) * t67;
t27 = -t47 * t77 + t72;
t25 = -t47 * t73 - t76;
t24 = -t67 + t77;
t20 = -g(3) * t47 - t61 * t45;
t3 = -t14 * t40 + t9 * t41;
t2 = -t14 * t41 - t9 * t40;
t1 = t93 * t48;
t6 = [0, t61, t62, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t27, -g(1) * t24 - g(2) * t26, -t62 * t45, -g(1) * t60 - g(2) * t66, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t63, -g(1) * (t11 * t43 - t46 * t5) - g(2) * (-t14 * t43 + t9 * t46) -g(1) * (t11 * t46 + t43 * t5) - g(2) * (-t14 * t46 - t9 * t43) -t63, -g(1) * (pkin(3) * t10 - pkin(4) * t5 + t11 * pkin(9) - qJ(5) * t4 + t60) - g(2) * (t15 * pkin(3) + t9 * pkin(4) - t14 * pkin(9) + t8 * qJ(5) + t66) 0, 0, 0, 0, 0, g(1) * t95 - g(2) * t3, -g(1) * t96 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t24 + t94, g(3) * t45 * t49 + g(1) * t27 - g(2) * t25, 0, -g(2) * t30 + (g(2) * t77 + t94) * pkin(2), 0, 0, 0, 0, 0, t93 * t51, -t1, -g(1) * (t14 * t80 + t15 * t43) - g(2) * (-t10 * t43 + t11 * t80) - g(3) * (t21 * t80 - t22 * t43) -g(1) * (-t14 * t84 + t15 * t46) - g(2) * (-t10 * t46 - t11 * t84) - g(3) * (-t21 * t84 - t22 * t46) t1, -g(1) * (t26 * pkin(2) + pkin(9) * t15) - g(2) * (-pkin(2) * t77 - t10 * pkin(9) + t30) - g(3) * (pkin(2) * t82 - t22 * pkin(9)) + t93 * (pkin(4) * t51 + qJ(5) * t48 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (t14 * t85 + t15 * t40) - g(2) * (-t10 * t40 + t11 * t85) - g(3) * (t21 * t85 - t22 * t40) -g(1) * (-t14 * t86 + t15 * t41) - g(2) * (-t10 * t41 - t11 * t86) - g(3) * (-t21 * t86 - t22 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t56, t57 * t46, -t57 * t43, -t56, -g(1) * (-t8 * pkin(4) + t9 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t16 * pkin(4) + t17 * qJ(5)) 0, 0, 0, 0, 0, t57 * t41, -t57 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t96 - g(3) * (-t17 * t40 - t21 * t41) g(1) * t3 + g(2) * t95 - g(3) * (-t17 * t41 + t21 * t40);];
taug_reg  = t6;
