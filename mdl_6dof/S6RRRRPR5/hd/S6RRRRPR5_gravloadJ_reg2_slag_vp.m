% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:30:58
% EndTime: 2019-05-07 20:31:01
% DurationCPUTime: 0.74s
% Computational Cost: add. (579->131), mult. (809->160), div. (0->0), fcn. (864->10), ass. (0->82)
t44 = qJ(2) + qJ(3);
t42 = cos(t44);
t50 = cos(qJ(4));
t52 = cos(qJ(1));
t85 = t52 * t50;
t46 = sin(qJ(4));
t48 = sin(qJ(1));
t88 = t48 * t46;
t20 = t42 * t88 + t85;
t86 = t52 * t46;
t87 = t48 * t50;
t21 = t42 * t87 - t86;
t45 = sin(qJ(6));
t49 = cos(qJ(6));
t110 = t20 * t49 - t21 * t45;
t41 = sin(t44);
t112 = g(3) * t41;
t22 = t42 * t86 - t87;
t23 = t42 * t85 + t88;
t5 = t22 * t49 - t23 * t45;
t66 = t45 * t50 - t46 * t49;
t114 = g(1) * t5 + g(2) * t110 - t66 * t112;
t24 = g(1) * t52 + g(2) * t48;
t111 = t24 * t41;
t113 = g(3) * t42 - t111;
t83 = t42 * pkin(3) + t41 * pkin(9);
t6 = t22 * t45 + t23 * t49;
t65 = t45 * t46 + t49 * t50;
t67 = t20 * t45 + t21 * t49;
t107 = g(1) * t6 + g(2) * t67 + t65 * t112;
t104 = -pkin(4) - pkin(5);
t47 = sin(qJ(2));
t103 = pkin(2) * t47;
t53 = -pkin(8) - pkin(7);
t100 = g(2) * t53;
t95 = t41 * t46;
t94 = t41 * t48;
t93 = t41 * t50;
t92 = t41 * t52;
t91 = t42 * t48;
t90 = t42 * t50;
t89 = t42 * t52;
t84 = t52 * t53;
t82 = qJ(5) * t46;
t81 = t48 * t103;
t80 = t52 * t103;
t51 = cos(qJ(2));
t43 = t51 * pkin(2);
t40 = t43 + pkin(1);
t27 = t52 * t40;
t79 = pkin(3) * t89 + pkin(9) * t92 + t27;
t78 = -pkin(3) - t82;
t76 = -t20 * pkin(4) + t21 * qJ(5);
t75 = -t22 * pkin(4) + t23 * qJ(5);
t74 = pkin(4) * t90 + t42 * t82 + t83;
t28 = pkin(9) * t91;
t73 = -pkin(10) * t91 + t28;
t32 = pkin(9) * t89;
t72 = -pkin(10) * t89 + t32;
t71 = t43 + t74;
t70 = -pkin(3) * t41 - t103;
t69 = g(1) * t20 - g(2) * t22;
t68 = g(1) * t48 - g(2) * t52;
t64 = -t40 - t83;
t62 = -t21 * pkin(4) - t20 * qJ(5) - t84;
t60 = t23 * pkin(4) + t22 * qJ(5) + t79;
t4 = g(1) * t22 + g(2) * t20 + g(3) * t95;
t58 = g(1) * t23 + g(2) * t21 + g(3) * t93;
t57 = -g(3) * t51 + t24 * t47;
t56 = (-g(1) * t64 + t100) * t48;
t55 = (pkin(4) * t50 - t78) * t111;
t54 = (g(3) * pkin(10) + t24 * (-t104 * t50 - t78)) * t41;
t29 = pkin(5) * t90;
t25 = qJ(5) * t93;
t19 = g(1) * t94 - g(2) * t92;
t11 = t24 * t42 + t112;
t9 = t113 * t50;
t8 = t113 * t46;
t7 = g(1) * t21 - g(2) * t23;
t2 = t113 * t65;
t1 = t113 * t66;
t3 = [0, 0, 0, 0, 0, 0, t68, t24, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t51, -t68 * t47, -t24, -g(1) * (-t48 * pkin(1) + t52 * pkin(7)) - g(2) * (t52 * pkin(1) + t48 * pkin(7)) 0, 0, 0, 0, 0, 0, t68 * t42, -t19, -t24, -g(1) * (-t48 * t40 - t84) - g(2) * (-t48 * t53 + t27) 0, 0, 0, 0, 0, 0, t7, -t69, t19, g(1) * t84 - g(2) * t79 + t56, 0, 0, 0, 0, 0, 0, t7, t19, t69, -g(1) * t62 - g(2) * t60 + t56, 0, 0, 0, 0, 0, 0, g(1) * t67 - g(2) * t6, g(1) * t110 - g(2) * t5, -t19, -g(1) * (-t21 * pkin(5) + t62) - g(2) * (t23 * pkin(5) - pkin(10) * t92 + t60) + (-g(1) * (t41 * pkin(10) + t64) + t100) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, g(3) * t47 + t24 * t51, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t11, 0, t57 * pkin(2), 0, 0, 0, 0, 0, 0, -t9, t8, -t11, -g(1) * (t70 * t52 + t32) - g(2) * (t70 * t48 + t28) - g(3) * (t43 + t83) 0, 0, 0, 0, 0, 0, -t9, -t11, -t8, -g(1) * (t32 - t80) - g(2) * (t28 - t81) - g(3) * t71 + t55, 0, 0, 0, 0, 0, 0, -t2, t1, t11, -g(1) * (t72 - t80) - g(2) * (t73 - t81) - g(3) * (t29 + t71) + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, -t11, -g(1) * (-pkin(3) * t92 + t32) - g(2) * (-pkin(3) * t94 + t28) - g(3) * t83, 0, 0, 0, 0, 0, 0, -t9, -t11, -t8, -g(1) * t32 - g(2) * t28 - g(3) * t74 + t55, 0, 0, 0, 0, 0, 0, -t2, t1, t11, -g(1) * t72 - g(2) * t73 - g(3) * (t29 + t74) + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t58, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, -t58, -g(1) * t75 - g(2) * t76 - g(3) * (-pkin(4) * t95 + t25) 0, 0, 0, 0, 0, 0, t114, -t107, 0, -g(1) * (-t22 * pkin(5) + t75) - g(2) * (-t20 * pkin(5) + t76) - g(3) * (t104 * t95 + t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t107, 0, 0;];
taug_reg  = t3;
