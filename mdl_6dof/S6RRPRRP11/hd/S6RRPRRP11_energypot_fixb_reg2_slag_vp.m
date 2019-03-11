% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:35
% EndTime: 2019-03-09 12:48:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (152->81), mult. (225->93), div. (0->0), fcn. (221->8), ass. (0->46)
t89 = sin(qJ(2));
t101 = qJ(3) * t89;
t90 = sin(qJ(1));
t92 = cos(qJ(2));
t109 = t90 * t92;
t119 = pkin(2) * t109 + t90 * t101;
t118 = g(3) * pkin(6);
t94 = -pkin(9) - pkin(8);
t117 = g(3) * t92;
t88 = sin(qJ(4));
t116 = t88 * pkin(4);
t115 = t89 * pkin(2) + pkin(6);
t91 = cos(qJ(4));
t77 = t91 * pkin(4) + pkin(3);
t87 = qJ(4) + qJ(5);
t78 = sin(t87);
t71 = pkin(5) * t78 + t116;
t114 = t71 * t89;
t113 = t90 * t78;
t79 = cos(t87);
t112 = t90 * t79;
t111 = t90 * t88;
t110 = t90 * t91;
t93 = cos(qJ(1));
t108 = t92 * t93;
t107 = t92 * t94;
t106 = t93 * t78;
t105 = t93 * t79;
t104 = t93 * t88;
t103 = t93 * t91;
t102 = t93 * pkin(1) + t90 * pkin(7);
t100 = t89 * t111;
t82 = t90 * pkin(1);
t99 = t82 + t119;
t98 = -t93 * pkin(7) + t82;
t97 = pkin(2) * t108 + t93 * t101 + t102;
t96 = -t92 * qJ(3) + t115;
t95 = g(1) * t93 + g(2) * t90;
t86 = -qJ(6) + t94;
t72 = g(1) * t90 - g(2) * t93;
t70 = pkin(5) * t79 + t77;
t69 = g(3) * t89 + t95 * t92;
t68 = t95 * t89 - t117;
t67 = -g(1) * (t89 * t106 + t112) - g(2) * (t89 * t113 - t105) + t78 * t117;
t66 = -g(1) * (t89 * t105 - t113) - g(2) * (t89 * t112 + t106) + t79 * t117;
t1 = [0, 0, 0, 0, 0, 0, -t95, t72, -g(3), -t118, 0, 0, 0, 0, 0, 0, -t69, t68, -t72, -g(1) * t102 - g(2) * t98 - t118, 0, 0, 0, 0, 0, 0, -t72, t69, -t68, -g(1) * t97 - g(2) * (t98 + t119) - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t104 + t110) - g(2) * (t100 - t103) + t88 * t117, -g(1) * (t89 * t103 - t111) - g(2) * (t89 * t110 + t104) + t91 * t117, -t69, -g(1) * (t90 * pkin(3) + pkin(8) * t108 + t97) - g(2) * (pkin(8) * t109 + (-pkin(3) - pkin(7)) * t93 + t99) - g(3) * (t89 * pkin(8) + t96) 0, 0, 0, 0, 0, 0, t67, t66, -t69, -g(1) * (t90 * t77 + t97) - g(2) * (pkin(4) * t100 - t90 * t107 + t99) - g(3) * (-t89 * t94 + (-qJ(3) - t116) * t92 + t115) + (-g(1) * (t89 * t116 - t107) - g(2) * (-pkin(7) - t77)) * t93, 0, 0, 0, 0, 0, 0, t67, t66, -t69, -g(1) * (t90 * t70 + t97) - g(2) * (-t86 * t109 + t90 * t114 + t99) - g(3) * (-t89 * t86 + (-qJ(3) - t71) * t92 + t115) + (-g(1) * (-t86 * t92 + t114) - g(2) * (-pkin(7) - t70)) * t93;];
U_reg  = t1;
