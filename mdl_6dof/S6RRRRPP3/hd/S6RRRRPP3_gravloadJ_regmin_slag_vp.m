% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = qJ(2) + qJ(3);
t38 = sin(t41);
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t60 = g(1) * t47 + g(2) * t44;
t86 = t60 * t38;
t39 = cos(t41);
t85 = -t39 * pkin(3) - t38 * pkin(9);
t8 = -g(3) * t39 + t86;
t43 = sin(qJ(2));
t83 = pkin(2) * t43;
t48 = -pkin(8) - pkin(7);
t80 = g(2) * t48;
t33 = t38 * pkin(5);
t42 = sin(qJ(4));
t78 = t38 * t42;
t45 = cos(qJ(4));
t77 = t38 * t45;
t76 = t38 * t47;
t75 = t39 * t44;
t74 = t39 * t45;
t73 = t39 * t47;
t72 = t44 * t42;
t71 = t44 * t45;
t70 = t47 * t42;
t69 = t47 * t45;
t68 = -pkin(4) - qJ(6);
t67 = qJ(5) * t42;
t66 = -pkin(3) - t67;
t15 = t39 * t72 + t69;
t16 = t39 * t71 - t70;
t65 = -t15 * pkin(4) + t16 * qJ(5);
t17 = t39 * t70 - t71;
t18 = t39 * t69 + t72;
t64 = -t17 * pkin(4) + t18 * qJ(5);
t63 = pkin(4) * t74 + t39 * t67 - t85;
t23 = pkin(9) * t75;
t62 = -t44 * t83 + t23;
t27 = pkin(9) * t73;
t61 = -t47 * t83 + t27;
t4 = g(1) * t15 - g(2) * t17;
t5 = g(1) * t16 - g(2) * t18;
t59 = g(1) * t44 - g(2) * t47;
t58 = qJ(6) * t74 + t33 + t63;
t46 = cos(qJ(2));
t40 = t46 * pkin(2);
t37 = t40 + pkin(1);
t57 = -t37 + t85;
t55 = -t16 * pkin(4) - t15 * qJ(5) - t47 * t48;
t53 = pkin(3) * t73 + t18 * pkin(4) + pkin(9) * t76 + t17 * qJ(5) + t47 * t37;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t78;
t51 = g(1) * t18 + g(2) * t16 + g(3) * t77;
t50 = (pkin(4) * t45 - t66) * t86;
t49 = (-t68 * t45 - t66) * t86;
t28 = pkin(5) * t73;
t24 = pkin(5) * t75;
t19 = qJ(5) * t77;
t14 = t59 * t38;
t9 = g(3) * t38 + t60 * t39;
t7 = -g(3) * t74 + t45 * t86;
t6 = t8 * t42;
t1 = [0, t59, t60, 0, 0, 0, 0, 0, t59 * t46, -t59 * t43, 0, 0, 0, 0, 0, t59 * t39, -t14, 0, 0, 0, 0, 0, t5, -t4, t14, -t5, t4, -g(1) * t55 - g(2) * t53 + (-g(1) * t57 + t80) * t44, t14, t4, t5, -g(1) * (-t16 * qJ(6) + t55) - g(2) * (pkin(5) * t76 + t18 * qJ(6) + t53) + (-g(1) * (t57 - t33) + t80) * t44; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t46 + t60 * t43, g(3) * t43 + t60 * t46, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t6, -t9, -t7, t6, -g(1) * t61 - g(2) * t62 - g(3) * (t40 + t63) + t50, -t9, t6, t7, -g(1) * (t28 + t61) - g(2) * (t24 + t62) - g(3) * (t40 + t58) + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t6, -t9, -t7, t6, -g(1) * t27 - g(2) * t23 - g(3) * t63 + t50, -t9, t6, t7, -g(1) * (t27 + t28) - g(2) * (t23 + t24) - g(3) * t58 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t51, 0, -t2, -t51, -g(1) * t64 - g(2) * t65 - g(3) * (-pkin(4) * t78 + t19) 0, -t51, t2, -g(1) * (-t17 * qJ(6) + t64) - g(2) * (-t15 * qJ(6) + t65) - g(3) * (t68 * t78 + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51;];
taug_reg  = t1;
