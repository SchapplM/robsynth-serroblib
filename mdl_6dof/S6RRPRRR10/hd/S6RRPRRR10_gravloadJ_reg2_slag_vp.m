% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t57 = sin(qJ(2));
t58 = sin(qJ(1));
t60 = cos(qJ(2));
t101 = cos(qJ(1));
t79 = cos(pkin(6));
t69 = t79 * t101;
t30 = t57 * t69 + t58 * t60;
t50 = pkin(12) + qJ(4);
t45 = sin(t50);
t46 = cos(t50);
t53 = sin(pkin(6));
t77 = t53 * t101;
t14 = t30 * t46 - t45 * t77;
t29 = t58 * t57 - t60 * t69;
t51 = qJ(5) + qJ(6);
t47 = sin(t51);
t48 = cos(t51);
t109 = t14 * t47 - t29 * t48;
t56 = sin(qJ(5));
t59 = cos(qJ(5));
t108 = t14 * t56 - t29 * t59;
t98 = t29 * t56;
t107 = t14 * t59 + t98;
t106 = t14 * t48 + t29 * t47;
t88 = t53 * t57;
t24 = t79 * t45 + t46 * t88;
t74 = t58 * t79;
t32 = t101 * t60 - t57 * t74;
t87 = t53 * t58;
t18 = t32 * t46 + t45 * t87;
t31 = t101 * t57 + t60 * t74;
t8 = -t18 * t56 + t31 * t59;
t83 = t59 * t60;
t105 = g(2) * t108 - g(3) * (-t24 * t56 - t53 * t83) - g(1) * t8;
t54 = cos(pkin(12));
t43 = t54 * pkin(3) + pkin(2);
t86 = t53 * t60;
t33 = t43 * t86;
t103 = g(3) * t33;
t102 = g(3) * t53;
t96 = t30 * t56;
t95 = t31 * t56;
t94 = t32 * t56;
t93 = t46 * t47;
t92 = t46 * t48;
t91 = t46 * t56;
t90 = t46 * t59;
t89 = t46 * t60;
t55 = -pkin(9) - qJ(3);
t85 = t55 * t57;
t84 = t56 * t60;
t82 = -t29 * t43 - t30 * t55;
t81 = -t31 * t43 - t32 * t55;
t80 = t101 * pkin(1) + pkin(8) * t87;
t52 = sin(pkin(12));
t78 = t52 * t87;
t76 = -t58 * pkin(1) + pkin(8) * t77;
t75 = -t30 * t45 - t46 * t77;
t73 = t52 * t77;
t72 = pkin(3) * t78 - t31 * t55 + t32 * t43 + t80;
t71 = pkin(4) * t46 + pkin(10) * t45;
t17 = t32 * t45 - t46 * t87;
t70 = g(1) * t75 + g(2) * t17;
t12 = g(1) * t29 - g(2) * t31;
t44 = t59 * pkin(5) + pkin(4);
t61 = -pkin(11) - pkin(10);
t68 = t44 * t46 - t45 * t61;
t67 = g(1) * t101 + g(2) * t58;
t66 = pkin(3) * t73 + t29 * t55 - t30 * t43 + t76;
t23 = -t45 * t88 + t79 * t46;
t65 = g(1) * t17 - g(2) * t75 - g(3) * t23;
t64 = g(1) * t18 + g(2) * t14 + g(3) * t24;
t10 = -g(1) * t31 - g(2) * t29 + g(3) * t86;
t63 = g(1) * t32 + g(2) * t30 + g(3) * t88;
t9 = t18 * t59 + t95;
t7 = t10 * t45;
t6 = t18 * t48 + t31 * t47;
t5 = -t18 * t47 + t31 * t48;
t2 = g(1) * t6 + g(2) * t106 - g(3) * (-t24 * t48 + t47 * t86);
t1 = -g(1) * t5 + g(2) * t109 - g(3) * (-t24 * t47 - t48 * t86);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t58 - g(2) * t101, t67, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t30 - g(2) * t32, -t12, -t67 * t53, -g(1) * t76 - g(2) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t54 + t73) - g(2) * (t32 * t54 + t78) -g(1) * (t30 * t52 + t54 * t77) - g(2) * (-t32 * t52 + t54 * t87) t12, -g(1) * (-t30 * pkin(2) - t29 * qJ(3) + t76) - g(2) * (t32 * pkin(2) + t31 * qJ(3) + t80) 0, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t70, t12, -g(1) * t66 - g(2) * t72, 0, 0, 0, 0, 0, 0, g(1) * t107 - g(2) * t9, -g(1) * t108 - g(2) * t8, -t70, -g(1) * (-pkin(4) * t14 + pkin(10) * t75 + t66) - g(2) * (t18 * pkin(4) + t17 * pkin(10) + t72) 0, 0, 0, 0, 0, 0, g(1) * t106 - g(2) * t6, -g(1) * t109 - g(2) * t5, -t70, -g(1) * (-pkin(5) * t98 - t14 * t44 - t61 * t75 + t66) - g(2) * (pkin(5) * t95 - t17 * t61 + t18 * t44 + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t63, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t54, t10 * t52, -t63, -g(1) * (-t31 * pkin(2) + t32 * qJ(3)) - g(2) * (-t29 * pkin(2) + t30 * qJ(3)) - (pkin(2) * t60 + qJ(3) * t57) * t102, 0, 0, 0, 0, 0, 0, -t10 * t46, t7, -t63, -g(1) * t81 - g(2) * t82 - g(3) * (-t53 * t85 + t33) 0, 0, 0, 0, 0, 0, -g(1) * (-t31 * t90 + t94) - g(2) * (-t29 * t90 + t96) - (t46 * t83 + t56 * t57) * t102, -g(1) * (t31 * t91 + t32 * t59) - g(2) * (t29 * t91 + t30 * t59) - (-t46 * t84 + t57 * t59) * t102, -t7, -g(1) * (-t71 * t31 + t81) - g(2) * (-t71 * t29 + t82) - t103 - (t71 * t60 - t85) * t102, 0, 0, 0, 0, 0, 0, -g(1) * (-t31 * t92 + t32 * t47) - g(2) * (-t29 * t92 + t30 * t47) - (t47 * t57 + t48 * t89) * t102, -g(1) * (t31 * t93 + t32 * t48) - g(2) * (t29 * t93 + t30 * t48) - (-t47 * t89 + t48 * t57) * t102, -t7, -g(1) * (pkin(5) * t94 - t68 * t31 + t81) - g(2) * (pkin(5) * t96 - t68 * t29 + t82) - t103 - (t68 * t60 + (pkin(5) * t56 - t55) * t57) * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t64, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t59, -t65 * t56, -t64, -g(1) * (-t17 * pkin(4) + t18 * pkin(10)) - g(2) * (pkin(4) * t75 + t14 * pkin(10)) - g(3) * (t23 * pkin(4) + t24 * pkin(10)) 0, 0, 0, 0, 0, 0, t65 * t48, -t65 * t47, -t64, -g(1) * (-t17 * t44 - t18 * t61) - g(2) * (-t14 * t61 + t44 * t75) - g(3) * (t23 * t44 - t24 * t61); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, g(1) * t9 + g(2) * t107 - g(3) * (-t24 * t59 + t53 * t84) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t105 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
