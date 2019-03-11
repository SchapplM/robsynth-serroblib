% Calculate inertial parameters regressor of potential energy for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:56
% EndTime: 2019-03-08 21:10:56
% DurationCPUTime: 0.21s
% Computational Cost: add. (254->84), mult. (578->108), div. (0->0), fcn. (697->10), ass. (0->51)
t90 = sin(pkin(6));
t96 = cos(qJ(2));
t120 = t90 * t96;
t123 = cos(qJ(3));
t111 = t90 * t123;
t113 = cos(pkin(6));
t93 = sin(qJ(3));
t94 = sin(qJ(2));
t78 = t94 * t111 + t113 * t93;
t125 = t78 * pkin(4) + qJ(5) * t120;
t124 = pkin(7) * t90;
t122 = t90 * t93;
t121 = t90 * t94;
t119 = pkin(5) + qJ(4);
t118 = pkin(8) - qJ(5);
t89 = sin(pkin(10));
t91 = cos(pkin(10));
t117 = t91 * pkin(1) + t89 * t124;
t110 = t94 * t113;
t74 = t91 * t110 + t89 * t96;
t65 = t91 * t111 + t74 * t93;
t116 = qJ(4) * t65;
t76 = -t89 * t110 + t91 * t96;
t67 = -t89 * t111 + t76 * t93;
t115 = qJ(4) * t67;
t114 = t113 * pkin(7) + qJ(1);
t112 = t76 * pkin(2) + t117;
t109 = t96 * t113;
t108 = t89 * pkin(1) - t91 * t124;
t107 = g(1) * t89 - g(2) * t91;
t106 = t74 * pkin(2) + t108;
t75 = t89 * t109 + t91 * t94;
t105 = pkin(8) * t75 + t112;
t104 = pkin(2) * t121 - pkin(8) * t120 + t114;
t77 = -t113 * t123 + t93 * t121;
t103 = g(1) * t67 + g(2) * t65 + g(3) * t77;
t66 = -t91 * t122 + t74 * t123;
t68 = t89 * t122 + t76 * t123;
t102 = g(1) * t68 + g(2) * t66 + g(3) * t78;
t101 = t78 * pkin(3) + t104;
t73 = -t91 * t109 + t89 * t94;
t100 = pkin(8) * t73 + t106;
t59 = -g(1) * t75 - g(2) * t73 + g(3) * t120;
t64 = t68 * pkin(3);
t99 = t68 * pkin(4) + t118 * t75 + t112 + t64;
t98 = t77 * qJ(4) + t101;
t62 = t66 * pkin(3);
t97 = t66 * pkin(4) + t118 * t73 + t106 + t62;
t95 = cos(qJ(6));
t92 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t89, t107, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t74 - g(3) * t121, -t59, -g(3) * t113 - t107 * t90, -g(1) * t117 - g(2) * t108 - g(3) * t114, 0, 0, 0, 0, 0, 0, -t102, t103, t59, -g(1) * t105 - g(2) * t100 - g(3) * t104, 0, 0, 0, 0, 0, 0, -t102, t59, -t103, -g(1) * (t105 + t64 + t115) - g(2) * (t100 + t62 + t116) - g(3) * t98, 0, 0, 0, 0, 0, 0, -t103, t102, -t59, -g(1) * (t99 + t115) - g(2) * (t97 + t116) - g(3) * (t98 + t125) 0, 0, 0, 0, 0, 0, -g(1) * (t67 * t95 - t75 * t92) - g(2) * (t65 * t95 - t73 * t92) - g(3) * (t92 * t120 + t77 * t95) -g(1) * (-t67 * t92 - t75 * t95) - g(2) * (-t65 * t92 - t73 * t95) - g(3) * (t95 * t120 - t77 * t92) -t102, -g(1) * (pkin(9) * t68 + t119 * t67 + t99) - g(2) * (pkin(9) * t66 + t119 * t65 + t97) - g(3) * (t78 * pkin(9) + t119 * t77 + t101 + t125);];
U_reg  = t1;
