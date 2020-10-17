% Calculate inertial parameters regressor of gravitation load for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:17:33
% EndTime: 2019-05-04 20:17:34
% DurationCPUTime: 0.47s
% Computational Cost: add. (676->103), mult. (1892->152), div. (0->0), fcn. (2447->14), ass. (0->71)
t76 = sin(pkin(11));
t79 = cos(pkin(12));
t62 = t76 * t79;
t75 = sin(pkin(12));
t80 = cos(pkin(11));
t65 = t80 * t75;
t82 = cos(pkin(6));
t36 = t82 * t65 + t62;
t44 = sin(qJ(3));
t89 = cos(qJ(3));
t60 = t76 * t75;
t67 = t80 * t79;
t52 = -t82 * t67 + t60;
t77 = sin(pkin(7));
t78 = sin(pkin(6));
t64 = t78 * t77;
t81 = cos(pkin(7));
t93 = t52 * t81 + t80 * t64;
t19 = t36 * t44 + t93 * t89;
t43 = sin(qJ(4));
t83 = qJ(5) * t43;
t46 = cos(qJ(4));
t88 = t19 * t46;
t96 = -pkin(4) * t88 - t19 * t83;
t37 = -t82 * t60 + t67;
t53 = t82 * t62 + t65;
t61 = t76 * t78;
t92 = t53 * t81 - t77 * t61;
t21 = t37 * t44 + t92 * t89;
t87 = t21 * t46;
t95 = -pkin(4) * t87 - t21 * t83;
t63 = t78 * t75;
t91 = t79 * t81 * t78 + t82 * t77;
t30 = t44 * t63 - t91 * t89;
t86 = t30 * t46;
t94 = -pkin(4) * t86 - t30 * t83;
t90 = pkin(5) + pkin(9);
t42 = sin(qJ(6));
t85 = t42 * t43;
t45 = cos(qJ(6));
t84 = t43 * t45;
t16 = t19 * pkin(3);
t20 = t36 * t89 - t93 * t44;
t74 = t20 * pkin(9) - t16;
t17 = t21 * pkin(3);
t22 = t37 * t89 - t92 * t44;
t73 = t22 * pkin(9) - t17;
t29 = t30 * pkin(3);
t31 = t91 * t44 + t89 * t63;
t72 = t31 * pkin(9) - t29;
t66 = t80 * t78;
t47 = t52 * t77 - t81 * t66;
t8 = t20 * t43 - t47 * t46;
t9 = t20 * t46 + t47 * t43;
t71 = -t8 * pkin(4) + t9 * qJ(5);
t48 = t53 * t77 + t81 * t61;
t10 = t22 * t43 - t48 * t46;
t11 = t22 * t46 + t48 * t43;
t70 = -t10 * pkin(4) + t11 * qJ(5);
t51 = -t79 * t64 + t82 * t81;
t23 = t31 * t43 - t51 * t46;
t24 = t31 * t46 + t51 * t43;
t69 = -t23 * pkin(4) + t24 * qJ(5);
t2 = g(1) * t10 + g(2) * t8 + g(3) * t23;
t59 = g(1) * t11 + g(2) * t9 + g(3) * t24;
t58 = g(1) * t21 + g(2) * t19 + g(3) * t30;
t57 = g(1) * t22 + g(2) * t20 + g(3) * t31;
t34 = -g(1) * t61 + g(2) * t66 - g(3) * t82;
t4 = t58 * t46;
t3 = t58 * t43;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t57, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t57, -g(1) * t73 - g(2) * t74 - g(3) * t72, 0, 0, 0, 0, 0, 0, -t57, -t4, t3, -g(1) * (t73 + t95) - g(2) * (t74 + t96) - g(3) * (t72 + t94) 0, 0, 0, 0, 0, 0, -g(1) * (-t21 * t85 + t22 * t45) - g(2) * (-t19 * t85 + t20 * t45) - g(3) * (-t30 * t85 + t31 * t45) -g(1) * (-t21 * t84 - t22 * t42) - g(2) * (-t19 * t84 - t20 * t42) - g(3) * (-t30 * t84 - t31 * t42) t4, -g(1) * (-pkin(10) * t87 + t90 * t22 - t17 + t95) - g(2) * (-pkin(10) * t88 + t90 * t20 - t16 + t96) - g(3) * (-pkin(10) * t86 + t90 * t31 - t29 + t94); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t59, -g(1) * t70 - g(2) * t71 - g(3) * t69, 0, 0, 0, 0, 0, 0, -t59 * t42, -t59 * t45, t2, -g(1) * (-t10 * pkin(10) + t70) - g(2) * (-t8 * pkin(10) + t71) - g(3) * (-t23 * pkin(10) + t69); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t45 - t21 * t42) - g(2) * (-t19 * t42 + t8 * t45) - g(3) * (t23 * t45 - t30 * t42) -g(1) * (-t10 * t42 - t21 * t45) - g(2) * (-t19 * t45 - t8 * t42) - g(3) * (-t23 * t42 - t30 * t45) 0, 0;];
taug_reg  = t1;
