% Calculate potential energy for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR13_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR13_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:13
% EndTime: 2019-12-31 21:43:13
% DurationCPUTime: 0.26s
% Computational Cost: add. (196->85), mult. (413->108), div. (0->0), fcn. (481->10), ass. (0->45)
t116 = rSges(5,1) + pkin(8);
t115 = rSges(4,3) + pkin(8);
t85 = cos(pkin(5));
t114 = t85 * pkin(7) + pkin(6);
t113 = pkin(9) + rSges(6,3);
t84 = sin(pkin(5));
t88 = sin(qJ(2));
t112 = t84 * t88;
t89 = sin(qJ(1));
t111 = t84 * t89;
t91 = cos(qJ(3));
t110 = t84 * t91;
t92 = cos(qJ(2));
t109 = t84 * t92;
t93 = cos(qJ(1));
t108 = t84 * t93;
t107 = t89 * t88;
t106 = t89 * t92;
t105 = t93 * t88;
t104 = t93 * t92;
t103 = t93 * pkin(1) + pkin(7) * t111;
t102 = rSges(5,3) + qJ(4);
t101 = pkin(2) * t112 + t114;
t74 = -t85 * t107 + t104;
t100 = t74 * pkin(2) + t103;
t87 = sin(qJ(3));
t70 = t88 * t110 + t85 * t87;
t99 = t70 * pkin(3) + t101;
t65 = t87 * t111 + t74 * t91;
t98 = t65 * pkin(3) + t100;
t72 = t85 * t105 + t106;
t82 = t89 * pkin(1);
t97 = t72 * pkin(2) - pkin(7) * t108 + t82;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t96 = t86 * rSges(6,1) + t90 * rSges(6,2) + qJ(4);
t63 = -t87 * t108 + t72 * t91;
t95 = t63 * pkin(3) + t97;
t94 = t90 * rSges(6,1) - t86 * rSges(6,2) + pkin(4) + pkin(8);
t73 = t85 * t106 + t105;
t71 = -t85 * t104 + t107;
t69 = t87 * t112 - t85 * t91;
t64 = -t89 * t110 + t74 * t87;
t62 = t91 * t108 + t72 * t87;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t93 * rSges(2,1) - t89 * rSges(2,2)) + g(2) * (t89 * rSges(2,1) + t93 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t74 * rSges(3,1) - t73 * rSges(3,2) + t103) + g(2) * (t72 * rSges(3,1) - t71 * rSges(3,2) + t82) + g(3) * (t85 * rSges(3,3) + t114) + (g(1) * rSges(3,3) * t89 + g(3) * (rSges(3,1) * t88 + rSges(3,2) * t92) + g(2) * (-rSges(3,3) - pkin(7)) * t93) * t84) - m(4) * (g(1) * (t65 * rSges(4,1) - t64 * rSges(4,2) + t115 * t73 + t100) + g(2) * (t63 * rSges(4,1) - t62 * rSges(4,2) + t115 * t71 + t97) + g(3) * (t70 * rSges(4,1) - t69 * rSges(4,2) - t115 * t109 + t101)) - m(5) * (g(1) * (-t65 * rSges(5,2) + t102 * t64 + t116 * t73 + t98) + g(2) * (-t63 * rSges(5,2) + t102 * t62 + t116 * t71 + t95) + g(3) * (-t70 * rSges(5,2) + t102 * t69 - t116 * t109 + t99)) - m(6) * (g(1) * (t113 * t65 + t96 * t64 + t94 * t73 + t98) + g(2) * (t113 * t63 + t96 * t62 + t94 * t71 + t95) + g(3) * (-t94 * t109 + t113 * t70 + t96 * t69 + t99));
U = t1;
