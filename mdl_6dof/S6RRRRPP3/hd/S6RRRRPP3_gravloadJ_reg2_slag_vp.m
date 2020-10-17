% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:10:48
% EndTime: 2019-05-07 18:10:50
% DurationCPUTime: 0.47s
% Computational Cost: add. (504->116), mult. (683->133), div. (0->0), fcn. (696->8), ass. (0->69)
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t19 = g(1) * t48 + g(2) * t45;
t42 = qJ(2) + qJ(3);
t39 = sin(t42);
t92 = t19 * t39;
t40 = cos(t42);
t72 = t40 * pkin(3) + t39 * pkin(9);
t8 = -g(3) * t40 + t92;
t44 = sin(qJ(2));
t90 = pkin(2) * t44;
t89 = pkin(3) * t39;
t49 = -pkin(8) - pkin(7);
t86 = g(2) * t49;
t34 = t39 * pkin(5);
t43 = sin(qJ(4));
t84 = t39 * t43;
t46 = cos(qJ(4));
t83 = t39 * t46;
t82 = t39 * t48;
t81 = t40 * t45;
t80 = t40 * t46;
t79 = t40 * t48;
t78 = t45 * t43;
t77 = t45 * t46;
t76 = t48 * t43;
t75 = t48 * t46;
t74 = t48 * t49;
t73 = -pkin(4) - qJ(6);
t71 = qJ(5) * t43;
t47 = cos(qJ(2));
t41 = t47 * pkin(2);
t38 = t41 + pkin(1);
t23 = t48 * t38;
t70 = pkin(3) * t79 + pkin(9) * t82 + t23;
t69 = -pkin(3) - t71;
t15 = t40 * t78 + t75;
t16 = t40 * t77 - t76;
t68 = -t15 * pkin(4) + t16 * qJ(5);
t17 = t40 * t76 - t77;
t18 = t40 * t75 + t78;
t67 = -t17 * pkin(4) + t18 * qJ(5);
t66 = pkin(4) * t80 + t40 * t71 + t72;
t24 = pkin(9) * t81;
t65 = -t45 * t90 + t24;
t28 = pkin(9) * t79;
t64 = -t48 * t90 + t28;
t63 = -t89 - t90;
t4 = g(1) * t15 - g(2) * t17;
t5 = g(1) * t16 - g(2) * t18;
t62 = g(1) * t45 - g(2) * t48;
t61 = qJ(6) * t80 + t34 + t66;
t60 = -t38 - t72;
t58 = -t16 * pkin(4) - t15 * qJ(5) - t74;
t56 = t18 * pkin(4) + t17 * qJ(5) + t70;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t84;
t54 = g(1) * t18 + g(2) * t16 + g(3) * t83;
t53 = -g(3) * t47 + t19 * t44;
t52 = (-g(1) * t60 + t86) * t45;
t51 = (pkin(4) * t46 - t69) * t92;
t50 = (-t73 * t46 - t69) * t92;
t29 = pkin(5) * t79;
t25 = pkin(5) * t81;
t20 = qJ(5) * t83;
t14 = t62 * t39;
t9 = g(3) * t39 + t19 * t40;
t7 = -g(3) * t80 + t46 * t92;
t6 = t8 * t43;
t1 = [0, 0, 0, 0, 0, 0, t62, t19, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t47, -t62 * t44, -t19, -g(1) * (-t45 * pkin(1) + t48 * pkin(7)) - g(2) * (t48 * pkin(1) + t45 * pkin(7)) 0, 0, 0, 0, 0, 0, t62 * t40, -t14, -t19, -g(1) * (-t45 * t38 - t74) - g(2) * (-t45 * t49 + t23) 0, 0, 0, 0, 0, 0, t5, -t4, t14, g(1) * t74 - g(2) * t70 + t52, 0, 0, 0, 0, 0, 0, t14, -t5, t4, -g(1) * t58 - g(2) * t56 + t52, 0, 0, 0, 0, 0, 0, t14, t4, t5, -g(1) * (-t16 * qJ(6) + t58) - g(2) * (pkin(5) * t82 + t18 * qJ(6) + t56) + (-g(1) * (t60 - t34) + t86) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(3) * t44 + t19 * t47, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t53 * pkin(2), 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (t48 * t63 + t28) - g(2) * (t45 * t63 + t24) - g(3) * (t41 + t72) 0, 0, 0, 0, 0, 0, -t9, -t7, t6, -g(1) * t64 - g(2) * t65 - g(3) * (t41 + t66) + t51, 0, 0, 0, 0, 0, 0, -t9, t6, t7, -g(1) * (t29 + t64) - g(2) * (t25 + t65) - g(3) * (t41 + t61) + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (-pkin(3) * t82 + t28) - g(2) * (-t45 * t89 + t24) - g(3) * t72, 0, 0, 0, 0, 0, 0, -t9, -t7, t6, -g(1) * t28 - g(2) * t24 - g(3) * t66 + t51, 0, 0, 0, 0, 0, 0, -t9, t6, t7, -g(1) * (t28 + t29) - g(2) * (t24 + t25) - g(3) * t61 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t54, -g(1) * t67 - g(2) * t68 - g(3) * (-pkin(4) * t84 + t20) 0, 0, 0, 0, 0, 0, 0, -t54, t2, -g(1) * (-t17 * qJ(6) + t67) - g(2) * (-t15 * qJ(6) + t68) - g(3) * (t73 * t84 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54;];
taug_reg  = t1;
