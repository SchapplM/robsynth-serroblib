% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t49 = sin(qJ(2));
t50 = sin(qJ(1));
t53 = cos(qJ(2));
t54 = cos(qJ(1));
t79 = cos(pkin(6));
t75 = t54 * t79;
t30 = t50 * t49 - t53 * t75;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t46 = sin(pkin(6));
t87 = t46 * t54;
t10 = t30 * t52 + t48 * t87;
t31 = t49 * t75 + t50 * t53;
t47 = sin(qJ(6));
t51 = cos(qJ(6));
t98 = t10 * t47 - t31 * t51;
t97 = t10 * t51 + t31 * t47;
t76 = t50 * t79;
t33 = -t49 * t76 + t54 * t53;
t96 = -g(1) * t33 - g(2) * t31;
t93 = g(3) * t46;
t90 = t46 * t49;
t89 = t46 * t50;
t88 = t46 * t53;
t86 = t47 * t52;
t85 = t51 * t52;
t84 = pkin(2) * t88 + qJ(3) * t90;
t83 = t54 * pkin(1) + pkin(8) * t89;
t82 = qJ(4) * t46;
t81 = t31 * qJ(3);
t80 = t33 * qJ(3);
t78 = pkin(3) * t88 + t84;
t77 = -t50 * pkin(1) + pkin(8) * t87;
t18 = t30 * pkin(2);
t74 = -t18 + t81;
t32 = t54 * t49 + t53 * t76;
t24 = t32 * pkin(2);
t73 = -t24 + t80;
t72 = pkin(4) * t90 + pkin(9) * t88 + t78;
t71 = pkin(5) * t52 + pkin(10) * t48;
t13 = t32 * t48 + t52 * t89;
t63 = -t30 * t48 + t52 * t87;
t70 = g(1) * t63 + g(2) * t13;
t7 = g(1) * t30 - g(2) * t32;
t8 = g(1) * t31 - g(2) * t33;
t69 = g(1) * t54 + g(2) * t50;
t68 = g(1) * t50 - g(2) * t54;
t17 = t30 * pkin(3);
t67 = t31 * pkin(4) - t30 * pkin(9) - t17 - t18;
t23 = t32 * pkin(3);
t66 = t33 * pkin(4) - t32 * pkin(9) - t23 - t24;
t65 = t33 * pkin(2) + t32 * qJ(3) + t83;
t28 = t48 * t88 - t79 * t52;
t62 = g(1) * t13 - g(2) * t63 - g(3) * t28;
t61 = -t31 * pkin(2) - t30 * qJ(3) + t77;
t14 = t32 * t52 - t48 * t89;
t29 = t79 * t48 + t52 * t88;
t60 = g(1) * t14 + g(2) * t10 - g(3) * t29;
t4 = -g(1) * t32 - g(2) * t30 + g(3) * t88;
t59 = g(3) * t90 - t96;
t58 = t33 * pkin(3) - t50 * t82 + t65;
t57 = -t31 * pkin(3) - t54 * t82 + t61;
t56 = t32 * pkin(4) + t33 * pkin(9) + t58;
t55 = -t30 * pkin(4) - t31 * pkin(9) + t57;
t34 = t69 * t46;
t15 = g(3) * t79 + t68 * t46;
t3 = t59 * t48;
t2 = t14 * t51 + t33 * t47;
t1 = -t14 * t47 + t33 * t51;
t5 = [0, 0, 0, 0, 0, 0, t68, t69, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t34, -g(1) * t77 - g(2) * t83, 0, 0, 0, 0, 0, 0, t8, -t34, t7, -g(1) * t61 - g(2) * t65, 0, 0, 0, 0, 0, 0, t7, -t8, t34, -g(1) * t57 - g(2) * t58, 0, 0, 0, 0, 0, 0, g(1) * t10 - g(2) * t14, t70, t8, -g(1) * t55 - g(2) * t56, 0, 0, 0, 0, 0, 0, g(1) * t97 - g(2) * t2, -g(1) * t98 - g(2) * t1, -t70, -g(1) * (-pkin(5) * t10 + pkin(10) * t63 + t55) - g(2) * (t14 * pkin(5) + t13 * pkin(10) + t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t59, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, -t59, -g(1) * t73 - g(2) * t74 - g(3) * t84, 0, 0, 0, 0, 0, 0, -t59, t4, 0, -g(1) * (-t23 + t73) - g(2) * (-t17 + t74) - g(3) * t78, 0, 0, 0, 0, 0, 0, -t59 * t52, t3, -t4, -g(1) * (t66 + t80) - g(2) * (t67 + t81) - g(3) * t72, 0, 0, 0, 0, 0, 0, -g(1) * (-t32 * t47 + t33 * t85) - g(2) * (-t30 * t47 + t31 * t85) - (t47 * t53 + t49 * t85) * t93, -g(1) * (-t32 * t51 - t33 * t86) - g(2) * (-t30 * t51 - t31 * t86) - (-t49 * t86 + t51 * t53) * t93, -t3, -g(1) * t66 - g(2) * t67 - g(3) * (t71 * t90 + t72) + t96 * (qJ(3) + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t60, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t51, -t62 * t47, -t60, -g(1) * (-t13 * pkin(5) + t14 * pkin(10)) - g(2) * (pkin(5) * t63 + t10 * pkin(10)) - g(3) * (t28 * pkin(5) - t29 * pkin(10)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t98 - g(3) * (t29 * t47 + t51 * t90) g(1) * t2 + g(2) * t97 - g(3) * (t29 * t51 - t47 * t90) 0, 0;];
taug_reg  = t5;
