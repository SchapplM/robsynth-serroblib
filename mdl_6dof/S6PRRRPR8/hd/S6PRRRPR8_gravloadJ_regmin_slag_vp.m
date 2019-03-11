% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t43 = sin(pkin(7));
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t45 = cos(pkin(12));
t79 = cos(pkin(6));
t73 = t45 * t79;
t77 = sin(pkin(12));
t58 = t77 * t49 - t52 * t73;
t78 = cos(pkin(7));
t44 = sin(pkin(6));
t84 = t44 * t45;
t93 = t43 * t84 + t58 * t78;
t68 = t79 * t77;
t59 = t45 * t49 + t52 * t68;
t74 = t44 * t77;
t92 = -t43 * t74 + t59 * t78;
t35 = t49 * t73 + t77 * t52;
t48 = sin(qJ(3));
t87 = cos(qJ(3));
t12 = t35 * t48 + t93 * t87;
t36 = t45 * t52 - t49 * t68;
t14 = t36 * t48 + t92 * t87;
t69 = t78 * t87;
t71 = t79 * t43;
t82 = t44 * t52;
t83 = t44 * t49;
t26 = t48 * t83 - t69 * t82 - t87 * t71;
t62 = g(1) * t14 + g(2) * t12 + g(3) * t26;
t91 = pkin(9) * t43;
t47 = sin(qJ(4));
t86 = t43 * t47;
t51 = cos(qJ(4));
t85 = t43 * t51;
t46 = sin(qJ(6));
t81 = t46 * t47;
t50 = cos(qJ(6));
t80 = t47 * t50;
t75 = t43 * t83;
t72 = t48 * t78;
t27 = t48 * t71 + (t87 * t49 + t52 * t72) * t44;
t57 = -t43 * t82 + t79 * t78;
t16 = t27 * t47 - t57 * t51;
t13 = t35 * t87 - t93 * t48;
t54 = t58 * t43 - t78 * t84;
t4 = t13 * t47 - t54 * t51;
t15 = t36 * t87 - t92 * t48;
t53 = t59 * t43 + t78 * t74;
t6 = t15 * t47 - t53 * t51;
t66 = g(1) * t6 + g(2) * t4 + g(3) * t16;
t17 = t27 * t51 + t57 * t47;
t5 = t13 * t51 + t54 * t47;
t7 = t15 * t51 + t53 * t47;
t65 = g(1) * t7 + g(2) * t5 + g(3) * t17;
t21 = -t36 * t72 - t59 * t87;
t10 = t21 * t47 - t36 * t85;
t33 = (-t49 * t72 + t87 * t52) * t44;
t22 = t33 * t47 - t51 * t75;
t19 = -t35 * t72 - t58 * t87;
t8 = t19 * t47 - t35 * t85;
t64 = g(1) * t10 + g(2) * t8 + g(3) * t22;
t11 = t21 * t51 + t36 * t86;
t23 = t33 * t51 + t47 * t75;
t9 = t19 * t51 + t35 * t86;
t63 = g(1) * t11 + g(2) * t9 + g(3) * t23;
t61 = g(1) * t15 + g(2) * t13 + g(3) * t27;
t18 = t35 * t69 - t58 * t48;
t20 = t36 * t69 - t59 * t48;
t32 = (t48 * t52 + t49 * t69) * t44;
t60 = g(1) * t20 + g(2) * t18 + g(3) * t32;
t3 = t62 * t51;
t2 = t62 * t47;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t59 + g(2) * t58 - g(3) * t82, g(1) * t36 + g(2) * t35 + g(3) * t83, 0, 0, 0, 0, 0, -g(1) * t21 - g(2) * t19 - g(3) * t33, t60, 0, 0, 0, 0, 0, -t63, t64, -t60, t63, -t64, -g(1) * (-pkin(2) * t59 + t21 * pkin(3) + t11 * pkin(4) + t20 * pkin(10) + t10 * qJ(5) + t36 * t91) - g(2) * (-pkin(2) * t58 + t19 * pkin(3) + t9 * pkin(4) + t18 * pkin(10) + t8 * qJ(5) + t35 * t91) - g(3) * (t33 * pkin(3) + t23 * pkin(4) + t32 * pkin(10) + t22 * qJ(5) + (pkin(2) * t52 + t49 * t91) * t44) 0, 0, 0, 0, 0, -g(1) * (t10 * t46 + t20 * t50) - g(2) * (t18 * t50 + t8 * t46) - g(3) * (t22 * t46 + t32 * t50) -g(1) * (t10 * t50 - t20 * t46) - g(2) * (-t18 * t46 + t8 * t50) - g(3) * (t22 * t50 - t32 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t61, 0, 0, 0, 0, 0, t3, -t2, -t61, -t3, t2, -pkin(10) * t61 + t62 * (pkin(4) * t51 + qJ(5) * t47 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t14 * t81 + t15 * t50) - g(2) * (-t12 * t81 + t13 * t50) - g(3) * (-t26 * t81 + t27 * t50) -g(1) * (-t14 * t80 - t15 * t46) - g(2) * (-t12 * t80 - t13 * t46) - g(3) * (-t26 * t80 - t27 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, -t66, -t65, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t16 * pkin(4) + t17 * qJ(5)) 0, 0, 0, 0, 0, -t65 * t46, -t65 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t46 + t6 * t50) - g(2) * (-t12 * t46 + t4 * t50) - g(3) * (t16 * t50 - t26 * t46) -g(1) * (-t14 * t50 - t6 * t46) - g(2) * (-t12 * t50 - t4 * t46) - g(3) * (-t16 * t46 - t26 * t50);];
taug_reg  = t1;
