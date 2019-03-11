% Calculate inertial parameters regressor of potential energy for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:28
% EndTime: 2019-03-09 09:52:28
% DurationCPUTime: 0.16s
% Computational Cost: add. (199->60), mult. (236->72), div. (0->0), fcn. (242->8), ass. (0->40)
t91 = qJ(2) + pkin(9);
t87 = sin(t91);
t88 = cos(t91);
t121 = pkin(3) * t88 + pkin(8) * t87;
t96 = cos(qJ(4));
t114 = t87 * t96;
t93 = sin(qJ(4));
t115 = t87 * t93;
t120 = pkin(4) * t114 + qJ(5) * t115;
t119 = g(3) * pkin(6);
t94 = sin(qJ(2));
t116 = t94 * pkin(2) + pkin(6);
t95 = sin(qJ(1));
t113 = t95 * t93;
t112 = t95 * t96;
t98 = cos(qJ(1));
t111 = t98 * t93;
t110 = t98 * t96;
t97 = cos(qJ(2));
t86 = t97 * pkin(2) + pkin(1);
t92 = -pkin(7) - qJ(3);
t109 = t95 * t86 + t98 * t92;
t108 = qJ(6) * t87;
t107 = t87 * pkin(3) + t116;
t106 = t98 * t86 - t95 * t92;
t105 = t121 * t95 + t109;
t104 = g(1) * t98 + g(2) * t95;
t103 = -t88 * pkin(8) + t107;
t102 = t121 * t98 + t106;
t69 = t88 * t113 + t110;
t70 = t88 * t112 - t111;
t101 = t70 * pkin(4) + t69 * qJ(5) + t105;
t71 = t88 * t111 - t112;
t100 = g(1) * t71 + g(2) * t69 + g(3) * t115;
t72 = t88 * t110 + t113;
t99 = t72 * pkin(4) + t71 * qJ(5) + t102;
t74 = g(1) * t95 - g(2) * t98;
t66 = -g(3) * t88 + t104 * t87;
t65 = -g(1) * t72 - g(2) * t70 - g(3) * t114;
t1 = [0, 0, 0, 0, 0, 0, -t104, t74, -g(3), -t119, 0, 0, 0, 0, 0, 0, -g(3) * t94 - t104 * t97, -g(3) * t97 + t104 * t94, -t74, -g(1) * (t98 * pkin(1) + t95 * pkin(7)) - g(2) * (t95 * pkin(1) - t98 * pkin(7)) - t119, 0, 0, 0, 0, 0, 0, -g(3) * t87 - t104 * t88, t66, -t74, -g(1) * t106 - g(2) * t109 - g(3) * t116, 0, 0, 0, 0, 0, 0, t65, t100, -t66, -g(1) * t102 - g(2) * t105 - g(3) * t103, 0, 0, 0, 0, 0, 0, t65, -t66, -t100, -g(1) * t99 - g(2) * t101 - g(3) * (t103 + t120) 0, 0, 0, 0, 0, 0, t65, -t100, t66, -g(1) * (t72 * pkin(5) - t98 * t108 + t99) - g(2) * (t70 * pkin(5) - t95 * t108 + t101) - g(3) * (pkin(5) * t114 + (-pkin(8) + qJ(6)) * t88 + t107 + t120);];
U_reg  = t1;
