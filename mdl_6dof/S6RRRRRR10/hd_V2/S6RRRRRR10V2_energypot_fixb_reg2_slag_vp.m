% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR10V2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:51:47
% EndTime: 2019-04-11 14:51:47
% DurationCPUTime: 0.19s
% Computational Cost: add. (199->59), mult. (265->88), div. (0->0), fcn. (295->12), ass. (0->45)
t96 = cos(qJ(2));
t82 = pkin(2) * t96 + pkin(1);
t87 = qJ(2) + qJ(3);
t84 = cos(t87);
t115 = pkin(3) * t84 + t82;
t114 = g(3) * pkin(4);
t91 = sin(qJ(2));
t112 = pkin(2) * t91 + pkin(4);
t83 = sin(t87);
t90 = sin(qJ(4));
t111 = t83 * t90;
t92 = sin(qJ(1));
t110 = t83 * t92;
t95 = cos(qJ(4));
t109 = t83 * t95;
t97 = cos(qJ(1));
t108 = t83 * t97;
t107 = t90 * t97;
t106 = t92 * t90;
t105 = t92 * t95;
t104 = t95 * t97;
t103 = pkin(5) * t110 + t115 * t92;
t102 = pkin(5) * t108 + t115 * t97;
t101 = g(1) * t97 + g(2) * t92;
t100 = pkin(3) * t83 - pkin(5) * t84 + t112;
t69 = t105 * t84 - t107;
t89 = sin(qJ(5));
t94 = cos(qJ(5));
t61 = -t110 * t94 + t69 * t89;
t71 = t104 * t84 + t106;
t63 = -t108 * t94 + t71 * t89;
t66 = t109 * t89 + t84 * t94;
t99 = g(1) * t63 + g(2) * t61 + g(3) * t66;
t68 = t106 * t84 + t104;
t70 = t107 * t84 - t105;
t98 = g(1) * t70 + g(2) * t68 + g(3) * t111;
t93 = cos(qJ(6));
t88 = sin(qJ(6));
t74 = g(1) * t92 - g(2) * t97;
t67 = t109 * t94 - t84 * t89;
t65 = -g(3) * t84 + t101 * t83;
t64 = t108 * t89 + t71 * t94;
t62 = t110 * t89 + t69 * t94;
t60 = -g(1) * t102 - g(2) * t103 - g(3) * t100;
t1 = [0, 0, 0, 0, 0, 0, -t101, t74, -g(3), -t114, 0, 0, 0, 0, 0, 0, -g(3) * t91 - t101 * t96, -g(3) * t96 + t101 * t91, -t74, -pkin(1) * t101 - t114, 0, 0, 0, 0, 0, 0, -g(3) * t83 - t101 * t84, t65, -t74, -g(3) * t112 - t101 * t82, 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 - g(3) * t109, t98, -t65, t60, 0, 0, 0, 0, 0, 0, -g(1) * t64 - g(2) * t62 - g(3) * t67, t99, -t98, t60, 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t93 + t70 * t88) - g(2) * (t62 * t93 + t68 * t88) - g(3) * (t111 * t88 + t67 * t93) -g(1) * (-t64 * t88 + t70 * t93) - g(2) * (-t62 * t88 + t68 * t93) - g(3) * (t111 * t93 - t67 * t88) -t99, -g(1) * (pkin(6) * t63 + t102) - g(2) * (pkin(6) * t61 + t103) - g(3) * (pkin(6) * t66 + t100);];
U_reg  = t1;
