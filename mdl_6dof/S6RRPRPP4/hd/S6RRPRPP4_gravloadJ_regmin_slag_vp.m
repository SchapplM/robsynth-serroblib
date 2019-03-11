% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = -qJ(5) - pkin(8);
t46 = cos(qJ(2));
t42 = sin(qJ(4));
t44 = sin(qJ(1));
t77 = t44 * t42;
t65 = pkin(4) * t77;
t43 = sin(qJ(2));
t81 = t43 * t44;
t90 = t41 * t81 + t46 * t65;
t47 = cos(qJ(1));
t71 = t47 * t42;
t29 = pkin(4) * t71;
t80 = t43 * t47;
t89 = t46 * t29 + t41 * t80;
t19 = g(1) * t47 + g(2) * t44;
t6 = g(3) * t43 + t19 * t46;
t88 = pkin(4) * t42;
t87 = g(1) * t44;
t83 = g(3) * t46;
t45 = cos(qJ(4));
t82 = t45 * pkin(4);
t36 = t46 * pkin(2);
t40 = qJ(4) + pkin(9);
t32 = sin(t40);
t79 = t44 * t32;
t33 = cos(t40);
t78 = t44 * t33;
t76 = t44 * t45;
t75 = t46 * t41;
t74 = t46 * t47;
t73 = t47 * t32;
t72 = t47 * t33;
t70 = t47 * t45;
t64 = t43 * t76;
t69 = pkin(4) * t64 + t29;
t34 = t43 * qJ(3);
t68 = t36 + t34;
t67 = qJ(3) * t46;
t66 = t45 * t83;
t63 = t43 * t71;
t62 = t43 * t70;
t31 = pkin(3) + t82;
t37 = t47 * pkin(7);
t61 = t47 * t31 + t44 * t75 + t37;
t60 = -pkin(1) - t36;
t59 = pkin(2) * t74 + t44 * pkin(7) + (pkin(1) + t34) * t47;
t25 = t44 * t67;
t58 = -pkin(2) * t81 + t25;
t27 = t47 * t67;
t57 = -pkin(2) * t80 + t27;
t56 = pkin(4) * t62 - t65;
t55 = -g(2) * t47 + t87;
t54 = pkin(5) * t32 - qJ(6) * t33;
t53 = t43 * t88 + t68 - t75;
t1 = -t43 * t72 + t79;
t3 = t43 * t78 + t73;
t52 = g(1) * t1 - g(2) * t3 + t33 * t83;
t51 = pkin(4) * t63 + t44 * t31 - t41 * t74 + t59;
t48 = ((-qJ(3) - t88) * t43 + t60) * t87;
t12 = t55 * t46;
t11 = t55 * t43;
t10 = -t43 * t77 + t70;
t9 = t64 + t71;
t8 = t63 + t76;
t7 = t62 - t77;
t5 = t19 * t43 - t83;
t4 = -t43 * t79 + t72;
t2 = t43 * t73 + t78;
t13 = [0, t55, t19, 0, 0, 0, 0, 0, t12, -t11, -t19, -t12, t11, -g(1) * t37 - g(2) * t59 - (t60 - t34) * t87, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t12, -g(1) * t61 - g(2) * t51 - t48, -g(1) * t4 - g(2) * t2, t12, -g(1) * t3 - g(2) * t1, -g(1) * (t4 * pkin(5) + t3 * qJ(6) + t61) - g(2) * (t2 * pkin(5) + t1 * qJ(6) + t51) - t48; 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, -t5, -t6, -g(1) * t57 - g(2) * t58 - g(3) * t68, 0, 0, 0, 0, 0, -t6 * t42, -t6 * t45, t5, -g(1) * (t57 + t89) - g(2) * (t58 + t90) - g(3) * t53, -t6 * t32, t5, t6 * t33, -g(1) * (t27 + t89) - g(2) * (t25 + t90) - g(3) * (t54 * t43 + t53) + t19 * (pkin(2) * t43 - t54 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t66, g(1) * t8 - g(2) * t10 - t42 * t83, 0, pkin(4) * t66 - g(1) * t56 - g(2) * t69, t52, 0, -g(1) * t2 + g(2) * t4 + t32 * t83, -g(1) * (-t1 * pkin(5) + t2 * qJ(6) + t56) - g(2) * (t3 * pkin(5) - t4 * qJ(6) + t69) - (-pkin(5) * t33 - qJ(6) * t32 - t82) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52;];
taug_reg  = t13;
