% Calculate inertial parameters regressor of potential energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:00
% EndTime: 2019-03-08 19:20:00
% DurationCPUTime: 0.25s
% Computational Cost: add. (312->85), mult. (712->134), div. (0->0), fcn. (881->12), ass. (0->53)
t82 = sin(pkin(11));
t85 = cos(pkin(11));
t90 = sin(qJ(2));
t93 = cos(qJ(2));
t72 = -t90 * t82 + t93 * t85;
t83 = sin(pkin(10));
t84 = sin(pkin(6));
t119 = t83 * t84;
t86 = cos(pkin(10));
t118 = t84 * t86;
t89 = sin(qJ(5));
t117 = t84 * t89;
t116 = t84 * t90;
t92 = cos(qJ(5));
t115 = t84 * t92;
t87 = cos(pkin(6));
t114 = t87 * t90;
t113 = t87 * t93;
t70 = pkin(2) * t114 + (-pkin(7) - qJ(3)) * t84;
t78 = t93 * pkin(2) + pkin(1);
t110 = t86 * t70 + t83 * t78;
t109 = t87 * pkin(7) + qJ(1);
t108 = -t83 * t70 + t86 * t78;
t107 = pkin(2) * t116 + t87 * qJ(3) + t109;
t106 = g(1) * t83 - g(2) * t86;
t105 = t93 * t82 + t90 * t85;
t102 = t72 * t87;
t56 = t86 * t102 - t105 * t83;
t103 = t105 * t87;
t57 = t86 * t103 + t83 * t72;
t104 = t57 * pkin(3) - t56 * qJ(4) + t110;
t58 = -t83 * t102 - t105 * t86;
t49 = t83 * t117 + t58 * t92;
t51 = t86 * t117 - t56 * t92;
t68 = t72 * t84;
t60 = t68 * t92 + t87 * t89;
t101 = g(1) * t49 - g(2) * t51 + g(3) * t60;
t100 = g(1) * t58 + g(2) * t56 + g(3) * t68;
t59 = -t83 * t103 + t86 * t72;
t69 = t105 * t84;
t99 = g(1) * t59 + g(2) * t57 + g(3) * t69;
t98 = t59 * pkin(3) - t58 * qJ(4) + t108;
t97 = t69 * pkin(3) - t68 * qJ(4) + t107;
t96 = pkin(4) * t119 + t59 * pkin(8) + t98;
t95 = t87 * pkin(4) + t69 * pkin(8) + t97;
t94 = -pkin(4) * t118 + t57 * pkin(8) + t104;
t91 = cos(qJ(6));
t88 = sin(qJ(6));
t67 = -g(3) * t87 - t106 * t84;
t61 = -t68 * t89 + t87 * t92;
t52 = -t86 * t115 - t56 * t89;
t50 = t83 * t115 - t58 * t89;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t86 - g(2) * t83, t106, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t83 * t114 + t86 * t93) - g(2) * (t86 * t114 + t83 * t93) - g(3) * t116, -g(1) * (-t83 * t113 - t86 * t90) - g(2) * (t86 * t113 - t83 * t90) - g(3) * t84 * t93, t67, -g(1) * (t86 * pkin(1) + pkin(7) * t119) - g(2) * (t83 * pkin(1) - pkin(7) * t118) - g(3) * t109, 0, 0, 0, 0, 0, 0, -t99, -t100, t67, -g(1) * t108 - g(2) * t110 - g(3) * t107, 0, 0, 0, 0, 0, 0, t67, t99, t100, -g(1) * t98 - g(2) * t104 - g(3) * t97, 0, 0, 0, 0, 0, 0, -g(1) * t50 - g(2) * t52 - g(3) * t61, t101, -t99, -g(1) * t96 - g(2) * t94 - g(3) * t95, 0, 0, 0, 0, 0, 0, -g(1) * (t50 * t91 + t59 * t88) - g(2) * (t52 * t91 + t57 * t88) - g(3) * (t61 * t91 + t69 * t88) -g(1) * (-t50 * t88 + t59 * t91) - g(2) * (-t52 * t88 + t57 * t91) - g(3) * (-t61 * t88 + t69 * t91) -t101, -g(1) * (t50 * pkin(5) + t49 * pkin(9) + t96) - g(2) * (t52 * pkin(5) - t51 * pkin(9) + t94) - g(3) * (t61 * pkin(5) + t60 * pkin(9) + t95);];
U_reg  = t1;
