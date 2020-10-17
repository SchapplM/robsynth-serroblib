% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP2
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:04:16
% EndTime: 2019-05-07 18:04:17
% DurationCPUTime: 0.32s
% Computational Cost: add. (435->88), mult. (591->111), div. (0->0), fcn. (617->8), ass. (0->64)
t38 = qJ(2) + qJ(3);
t35 = sin(t38);
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t56 = g(1) * t44 + g(2) * t41;
t51 = t56 * t35;
t36 = cos(t38);
t85 = -t36 * pkin(3) - t35 * pkin(9);
t8 = -g(3) * t36 + t51;
t84 = -pkin(4) - pkin(5);
t40 = sin(qJ(2));
t83 = pkin(2) * t40;
t82 = g(1) * t41;
t45 = -pkin(8) - pkin(7);
t79 = g(2) * t45;
t39 = sin(qJ(4));
t77 = t35 * t39;
t42 = cos(qJ(4));
t76 = t35 * t42;
t75 = t35 * t44;
t74 = t36 * t42;
t73 = t36 * t44;
t72 = t41 * t39;
t71 = t41 * t42;
t70 = t44 * t39;
t69 = t44 * t42;
t68 = qJ(5) * t39;
t67 = qJ(6) * t36;
t66 = t35 * qJ(6);
t65 = t41 * t83;
t64 = t44 * t83;
t63 = -pkin(3) - t68;
t15 = t36 * t72 + t69;
t16 = t36 * t71 - t70;
t62 = -t15 * pkin(4) + t16 * qJ(5);
t17 = t36 * t70 - t71;
t18 = t36 * t69 + t72;
t61 = -t17 * pkin(4) + t18 * qJ(5);
t60 = pkin(4) * t74 + t36 * t68 - t85;
t22 = t41 * t36 * pkin(9);
t59 = -t41 * t67 + t22;
t26 = pkin(9) * t73;
t58 = -t44 * t67 + t26;
t43 = cos(qJ(2));
t37 = t43 * pkin(2);
t57 = t37 + t60;
t4 = g(1) * t15 - g(2) * t17;
t55 = -g(2) * t44 + t82;
t34 = t37 + pkin(1);
t54 = -t34 + t85;
t52 = -t16 * pkin(4) - t15 * qJ(5) - t44 * t45;
t50 = pkin(3) * t73 + t18 * pkin(4) + pkin(9) * t75 + t17 * qJ(5) + t44 * t34;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t77;
t48 = g(1) * t18 + g(2) * t16 + g(3) * t76;
t47 = (pkin(4) * t42 - t63) * t51;
t46 = (g(3) * qJ(6) + t56 * (-t84 * t42 - t63)) * t35;
t23 = pkin(5) * t74;
t19 = qJ(5) * t76;
t14 = -g(2) * t75 + t35 * t82;
t9 = g(3) * t35 + t56 * t36;
t7 = t8 * t42;
t6 = t8 * t39;
t5 = g(1) * t16 - g(2) * t18;
t1 = [0, t55, t56, 0, 0, 0, 0, 0, t55 * t43, -t55 * t40, 0, 0, 0, 0, 0, t55 * t36, -t14, 0, 0, 0, 0, 0, t5, -t4, t5, t14, t4, -g(1) * t52 - g(2) * t50 + (-g(1) * t54 + t79) * t41, t5, t4, -t14, -g(1) * (-t16 * pkin(5) + t52) - g(2) * (t18 * pkin(5) - t44 * t66 + t50) + (-g(1) * (t54 + t66) + t79) * t41; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t43 + t56 * t40, g(3) * t40 + t56 * t43, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t6, t7, -t9, t6, -g(1) * (t26 - t64) - g(2) * (t22 - t65) - g(3) * t57 + t47, t7, t6, t9, -g(1) * (t58 - t64) - g(2) * (t59 - t65) - g(3) * (t23 + t57) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t6, t7, -t9, t6, -g(1) * t26 - g(2) * t22 - g(3) * t60 + t47, t7, t6, t9, -g(1) * t58 - g(2) * t59 - g(3) * (t23 + t60) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t48, t2, 0, -t48, -g(1) * t61 - g(2) * t62 - g(3) * (-pkin(4) * t77 + t19) t2, -t48, 0, -g(1) * (-t17 * pkin(5) + t61) - g(2) * (-t15 * pkin(5) + t62) - g(3) * (t84 * t77 + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t1;
