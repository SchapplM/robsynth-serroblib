% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:33
% EndTime: 2019-12-31 21:46:35
% DurationCPUTime: 0.56s
% Computational Cost: add. (364->117), mult. (943->173), div. (0->0), fcn. (1144->10), ass. (0->71)
t51 = sin(qJ(2));
t52 = sin(qJ(1));
t55 = cos(qJ(2));
t76 = cos(pkin(5));
t93 = cos(qJ(1));
t62 = t76 * t93;
t32 = t51 * t62 + t52 * t55;
t50 = sin(qJ(3));
t54 = cos(qJ(3));
t48 = sin(pkin(5));
t74 = t48 * t93;
t14 = t32 * t50 + t54 * t74;
t31 = t52 * t51 - t55 * t62;
t49 = sin(qJ(5));
t53 = cos(qJ(5));
t101 = t14 * t49 + t31 * t53;
t100 = t14 * t53 - t31 * t49;
t77 = qJ(4) * t50;
t90 = t31 * t54;
t99 = -pkin(3) * t90 - t31 * t77;
t70 = t52 * t76;
t33 = t51 * t93 + t55 * t70;
t89 = t33 * t54;
t98 = -pkin(3) * t89 - t33 * t77;
t97 = pkin(4) + pkin(8);
t96 = g(3) * t48;
t95 = t31 * pkin(8);
t94 = t33 * pkin(8);
t88 = t48 * t51;
t87 = t48 * t52;
t86 = t48 * t54;
t85 = t48 * t55;
t84 = t49 * t50;
t83 = t49 * t55;
t82 = t50 * t53;
t81 = t53 * t55;
t80 = t54 * t55;
t79 = pkin(2) * t85 + pkin(8) * t88;
t78 = t93 * pkin(1) + pkin(7) * t87;
t34 = -t51 * t70 + t55 * t93;
t75 = t34 * pkin(2) + t78;
t73 = -t52 * pkin(1) + pkin(7) * t74;
t25 = t31 * pkin(2);
t72 = t32 * pkin(8) - t25;
t27 = t33 * pkin(2);
t71 = t34 * pkin(8) - t27;
t15 = t32 * t54 - t50 * t74;
t69 = -t14 * pkin(3) + t15 * qJ(4);
t18 = t34 * t50 - t52 * t86;
t19 = t34 * t54 + t50 * t87;
t68 = -t18 * pkin(3) + t19 * qJ(4);
t29 = t50 * t88 - t54 * t76;
t30 = t50 * t76 + t51 * t86;
t67 = -t29 * pkin(3) + t30 * qJ(4);
t66 = t48 * pkin(3) * t80 + t77 * t85 + t79;
t65 = -t32 * pkin(2) + t73;
t64 = -g(1) * t14 + g(2) * t18;
t63 = -g(1) * t15 + g(2) * t19;
t9 = g(1) * t31 - g(2) * t33;
t61 = g(1) * t93 + g(2) * t52;
t60 = t19 * pkin(3) + t18 * qJ(4) + t75;
t2 = g(1) * t18 + g(2) * t14 + g(3) * t29;
t59 = g(1) * t19 + g(2) * t15 + g(3) * t30;
t58 = -pkin(3) * t15 - qJ(4) * t14 + t65;
t57 = -g(1) * t33 - g(2) * t31 + g(3) * t85;
t56 = g(1) * t34 + g(2) * t32 + g(3) * t88;
t7 = t57 * t54;
t6 = t57 * t50;
t5 = t18 * t49 + t33 * t53;
t4 = t18 * t53 - t33 * t49;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t52 - g(2) * t93, t61, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t32 - g(2) * t34, -t9, -t61 * t48, -g(1) * t73 - g(2) * t78, 0, 0, 0, 0, 0, 0, -t63, t64, t9, -g(1) * (t65 - t95) - g(2) * (t75 + t94), 0, 0, 0, 0, 0, 0, t9, t63, -t64, -g(1) * (t58 - t95) - g(2) * (t60 + t94), 0, 0, 0, 0, 0, 0, g(1) * t101 - g(2) * t5, g(1) * t100 - g(2) * t4, -t63, -g(1) * (-pkin(9) * t15 - t31 * t97 + t58) - g(2) * (t19 * pkin(9) + t33 * t97 + t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t56, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t56, -g(1) * t71 - g(2) * t72 - g(3) * t79, 0, 0, 0, 0, 0, 0, -t56, t7, -t6, -g(1) * (t71 + t98) - g(2) * (t72 + t99) - g(3) * t66, 0, 0, 0, 0, 0, 0, -g(1) * (-t33 * t84 + t34 * t53) - g(2) * (-t31 * t84 + t32 * t53) - (t50 * t83 + t51 * t53) * t96, -g(1) * (-t33 * t82 - t34 * t49) - g(2) * (-t31 * t82 - t32 * t49) - (-t49 * t51 + t50 * t81) * t96, -t7, -g(1) * (-pkin(9) * t89 + t34 * t97 - t27 + t98) - g(2) * (-pkin(9) * t90 + t32 * t97 - t25 + t99) - g(3) * ((pkin(4) * t51 + pkin(9) * t80) * t48 + t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t59, -g(1) * t68 - g(2) * t69 - g(3) * t67, 0, 0, 0, 0, 0, 0, -t59 * t49, -t59 * t53, t2, -g(1) * (-t18 * pkin(9) + t68) - g(2) * (-t14 * pkin(9) + t69) - g(3) * (-t29 * pkin(9) + t67); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t100 - g(3) * (t29 * t53 + t48 * t83), g(1) * t5 + g(2) * t101 - g(3) * (-t29 * t49 + t48 * t81), 0, 0;];
taug_reg = t1;
