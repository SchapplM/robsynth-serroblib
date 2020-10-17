% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:16:27
% EndTime: 2019-05-07 06:16:29
% DurationCPUTime: 0.51s
% Computational Cost: add. (360->109), mult. (939->162), div. (0->0), fcn. (1152->10), ass. (0->59)
t52 = sin(qJ(3));
t94 = qJ(4) * t52 + pkin(2);
t53 = sin(qJ(2));
t54 = sin(qJ(1));
t57 = cos(qJ(2));
t74 = cos(pkin(6));
t90 = cos(qJ(1));
t62 = t74 * t90;
t35 = t53 * t62 + t54 * t57;
t56 = cos(qJ(3));
t50 = sin(pkin(6));
t71 = t50 * t90;
t17 = t35 * t52 + t56 * t71;
t34 = t54 * t53 - t57 * t62;
t51 = sin(qJ(6));
t55 = cos(qJ(6));
t93 = t17 * t51 + t34 * t55;
t92 = t17 * t55 - t34 * t51;
t91 = g(3) * t50;
t87 = t34 * t56;
t70 = t54 * t74;
t36 = t90 * t53 + t57 * t70;
t86 = t36 * t56;
t85 = t50 * t53;
t84 = t50 * t54;
t83 = t50 * t56;
t82 = t50 * t57;
t81 = t51 * t52;
t80 = t51 * t57;
t79 = t52 * t55;
t78 = t55 * t57;
t77 = t56 * t57;
t76 = pkin(9) - qJ(5);
t73 = -pkin(3) * t87 - t94 * t34;
t72 = -pkin(3) * t86 - t94 * t36;
t18 = t35 * t56 - t52 * t71;
t69 = -t17 * pkin(3) + t18 * qJ(4);
t37 = -t53 * t70 + t90 * t57;
t21 = t37 * t52 - t54 * t83;
t22 = t37 * t56 + t52 * t84;
t68 = -t21 * pkin(3) + t22 * qJ(4);
t32 = t52 * t85 - t74 * t56;
t33 = t74 * t52 + t53 * t83;
t67 = -t32 * pkin(3) + t33 * qJ(4);
t66 = t50 * pkin(3) * t77 + pkin(9) * t85 + t94 * t82;
t65 = -g(1) * t17 + g(2) * t21;
t64 = -g(1) * t18 + g(2) * t22;
t63 = g(1) * t34 - g(2) * t36;
t61 = t90 * pkin(1) + t37 * pkin(2) + t22 * pkin(3) + pkin(8) * t84 + t21 * qJ(4);
t2 = g(1) * t21 + g(2) * t17 + g(3) * t32;
t60 = g(1) * t22 + g(2) * t18 + g(3) * t33;
t59 = -t54 * pkin(1) - t35 * pkin(2) - pkin(3) * t18 + pkin(8) * t71 - qJ(4) * t17;
t58 = -g(1) * t36 - g(2) * t34 + g(3) * t82;
t11 = g(1) * t37 + g(2) * t35 + g(3) * t85;
t9 = t58 * t56;
t8 = t58 * t52;
t7 = t21 * t55 - t36 * t51;
t6 = -t21 * t51 - t36 * t55;
t1 = [0, g(1) * t54 - g(2) * t90, g(1) * t90 + g(2) * t54, 0, 0, 0, 0, 0, g(1) * t35 - g(2) * t37, -t63, 0, 0, 0, 0, 0, -t64, t65, -t64, t63, -t65, -g(1) * (-t34 * pkin(9) + t59) - g(2) * (t36 * pkin(9) + t61) -t65, t64, -t63, -g(1) * (-pkin(4) * t18 - t76 * t34 + t59) - g(2) * (t22 * pkin(4) + t76 * t36 + t61) 0, 0, 0, 0, 0, g(1) * t92 - g(2) * t7, -g(1) * t93 - g(2) * t6; 0, 0, 0, 0, 0, 0, 0, 0, -t58, t11, 0, 0, 0, 0, 0, -t9, t8, -t9, -t11, -t8, -g(1) * (t37 * pkin(9) + t72) - g(2) * (t35 * pkin(9) + t73) - g(3) * t66, -t8, t9, t11, -g(1) * (-pkin(4) * t86 + t76 * t37 + t72) - g(2) * (-pkin(4) * t87 + t76 * t35 + t73) - g(3) * ((pkin(4) * t77 - qJ(5) * t53) * t50 + t66) 0, 0, 0, 0, 0, -g(1) * (-t36 * t79 - t37 * t51) - g(2) * (-t34 * t79 - t35 * t51) - (-t51 * t53 + t52 * t78) * t91, -g(1) * (t36 * t81 - t37 * t55) - g(2) * (t34 * t81 - t35 * t55) - (-t52 * t80 - t53 * t55) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t60, t2, 0, -t60, -g(1) * t68 - g(2) * t69 - g(3) * t67, -t60, -t2, 0, -g(1) * (-t21 * pkin(4) + t68) - g(2) * (-t17 * pkin(4) + t69) - g(3) * (-t32 * pkin(4) + t67) 0, 0, 0, 0, 0, -t60 * t55, t60 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t93 - g(3) * (-t32 * t51 + t50 * t78) g(1) * t7 + g(2) * t92 - g(3) * (-t32 * t55 - t50 * t80);];
taug_reg  = t1;
