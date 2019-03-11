% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t36 = g(1) * t63 + g(2) * t60;
t59 = sin(qJ(2));
t105 = t36 * t59;
t56 = sin(pkin(9));
t57 = cos(pkin(9));
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t75 = t56 * t58 + t57 * t61;
t104 = t75 * t59;
t50 = t59 * qJ(3);
t62 = cos(qJ(2));
t89 = t62 * pkin(2) + t50;
t98 = t56 * t61;
t103 = t59 * (t57 * t58 - t98);
t102 = g(1) * t60;
t99 = g(3) * t59;
t96 = t57 * t62;
t95 = t59 * t60;
t94 = t59 * t63;
t93 = t60 * t62;
t92 = t62 * t63;
t91 = t63 * t56;
t90 = t63 * t57;
t88 = t63 * pkin(1) + t60 * pkin(7);
t87 = qJ(3) * t62;
t86 = qJ(4) * t56;
t85 = -pkin(2) - t86;
t84 = pkin(3) * t96 + t62 * t86 + t89;
t83 = pkin(2) * t92 + t63 * t50 + t88;
t40 = t60 * t87;
t82 = -pkin(8) * t93 + t40;
t44 = t63 * t87;
t81 = -pkin(8) * t92 + t44;
t31 = t56 * t93 + t90;
t32 = t57 * t93 - t91;
t7 = t31 * t61 - t32 * t58;
t33 = -t60 * t57 + t62 * t91;
t34 = t60 * t56 + t62 * t90;
t9 = -t33 * t61 + t34 * t58;
t80 = g(1) * t7 + g(2) * t9;
t53 = t63 * pkin(7);
t79 = -t32 * pkin(3) - t31 * qJ(4) + t53;
t78 = pkin(4) * t96 + t84;
t77 = g(1) * t31 - g(2) * t33;
t76 = -g(2) * t63 + t102;
t6 = t31 * t58 + t32 * t61;
t73 = -t32 * pkin(4) + pkin(8) * t95 + t79;
t1 = g(1) * t9 - g(2) * t7 + g(3) * t103;
t10 = t33 * t58 + t34 * t61;
t70 = g(1) * t10 + g(2) * t6 + g(3) * t104;
t14 = t60 * t103;
t16 = t63 * t103;
t25 = t58 * t96 - t62 * t98;
t69 = g(1) * t16 + g(2) * t14 - g(3) * t25;
t67 = t34 * pkin(3) + t33 * qJ(4) + t83;
t66 = (-pkin(1) - t89) * t102;
t21 = -g(3) * t62 + t105;
t65 = t34 * pkin(4) - pkin(8) * t94 + t67;
t64 = (g(3) * pkin(8) + t36 * (-(-pkin(3) - pkin(4)) * t57 - t85)) * t59;
t35 = g(1) * t95 - g(2) * t94;
t26 = t75 * t62;
t22 = t36 * t62 + t99;
t17 = t63 * t104;
t15 = t60 * t104;
t13 = t21 * t57;
t12 = t21 * t56;
t11 = g(1) * t32 - g(2) * t34;
t4 = -g(1) * t33 - g(2) * t31 - t56 * t99;
t3 = g(1) * t17 + g(2) * t15 - g(3) * t26;
t2 = g(1) * t6 - g(2) * t10;
t5 = [0, 0, 0, 0, 0, 0, t76, t36, 0, 0, 0, 0, 0, 0, 0, 0, t76 * t62, -t35, -t36, -g(1) * (-t60 * pkin(1) + t53) - g(2) * t88, 0, 0, 0, 0, 0, 0, t11, -t77, t35, -g(1) * t53 - g(2) * t83 - t66, 0, 0, 0, 0, 0, 0, t11, t35, t77, -g(1) * t79 - g(2) * t67 - t66, 0, 0, 0, 0, 0, 0, t2, t80, -t35, -g(1) * t73 - g(2) * t65 - t66, 0, 0, 0, 0, 0, 0, t2, -t35, -t80, -g(1) * (-pkin(5) * t6 + t7 * qJ(6) + t73) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t65) - t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t22, -g(1) * (-pkin(2) * t94 + t44) - g(2) * (-pkin(2) * t95 + t40) - g(3) * t89, 0, 0, 0, 0, 0, 0, t13, -t22, t12, -g(1) * t44 - g(2) * t40 - g(3) * t84 + (pkin(3) * t57 - t85) * t105, 0, 0, 0, 0, 0, 0, t3, -t69, t22, -g(1) * t81 - g(2) * t82 - g(3) * t78 + t64, 0, 0, 0, 0, 0, 0, t3, t22, t69, -g(1) * (-t17 * pkin(5) - t16 * qJ(6) + t81) - g(2) * (-t15 * pkin(5) - t14 * qJ(6) + t82) - g(3) * (t26 * pkin(5) + t25 * qJ(6) + t78) + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t70, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t70, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (pkin(5) * t7 + t6 * qJ(6)) - g(3) * (-pkin(5) * t103 + qJ(6) * t104); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
