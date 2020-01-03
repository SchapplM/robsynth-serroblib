% Calculate potential energy for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:31:54
% EndTime: 2019-12-31 22:31:54
% DurationCPUTime: 0.36s
% Computational Cost: add. (227->103), mult. (371->138), div. (0->0), fcn. (422->12), ass. (0->44)
t92 = sin(qJ(3));
t116 = pkin(3) * t92;
t90 = cos(pkin(5));
t115 = t90 * pkin(7) + pkin(6);
t114 = pkin(8) + rSges(4,3);
t113 = pkin(10) + rSges(6,3);
t89 = sin(pkin(5));
t93 = sin(qJ(2));
t112 = t89 * t93;
t94 = sin(qJ(1));
t111 = t89 * t94;
t97 = cos(qJ(2));
t110 = t89 * t97;
t98 = cos(qJ(1));
t109 = t89 * t98;
t108 = t93 * t94;
t107 = t93 * t98;
t106 = t94 * t97;
t105 = t97 * t98;
t104 = t98 * pkin(1) + pkin(7) * t111;
t103 = t92 * t111;
t96 = cos(qJ(3));
t82 = pkin(3) * t96 + pkin(2);
t99 = -pkin(9) - pkin(8);
t102 = t99 * t110 + t82 * t112 + t90 * t116 + t115;
t72 = t90 * t106 + t107;
t73 = -t90 * t108 + t105;
t101 = pkin(3) * t103 - t72 * t99 + t73 * t82 + t104;
t70 = -t90 * t105 + t108;
t71 = t90 * t107 + t106;
t86 = t94 * pkin(1);
t100 = t71 * t82 + t86 + (-pkin(7) - t116) * t109 - t70 * t99;
t95 = cos(qJ(5));
t91 = sin(qJ(5));
t88 = qJ(3) + qJ(4);
t84 = cos(t88);
t83 = sin(t88);
t67 = t84 * t112 + t83 * t90;
t66 = t83 * t112 - t90 * t84;
t63 = t83 * t111 + t73 * t84;
t62 = -t84 * t111 + t73 * t83;
t61 = -t83 * t109 + t71 * t84;
t60 = t84 * t109 + t71 * t83;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t98 - rSges(2,2) * t94) + g(2) * (rSges(2,1) * t94 + rSges(2,2) * t98) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t73 - rSges(3,2) * t72 + t104) + g(2) * (rSges(3,1) * t71 - rSges(3,2) * t70 + t86) + g(3) * (rSges(3,3) * t90 + t115) + (g(1) * rSges(3,3) * t94 + g(3) * (rSges(3,1) * t93 + rSges(3,2) * t97) + g(2) * (-rSges(3,3) - pkin(7)) * t98) * t89) - m(4) * (g(1) * (t73 * pkin(2) + (t73 * t96 + t103) * rSges(4,1) + (t96 * t111 - t73 * t92) * rSges(4,2) + t114 * t72 + t104) + g(2) * (t71 * pkin(2) + t86 - pkin(7) * t109 + (-t92 * t109 + t71 * t96) * rSges(4,1) + (-t96 * t109 - t71 * t92) * rSges(4,2) + t114 * t70) + g(3) * ((t92 * rSges(4,1) + t96 * rSges(4,2)) * t90 + (-t114 * t97 + (t96 * rSges(4,1) - t92 * rSges(4,2) + pkin(2)) * t93) * t89 + t115)) - m(5) * (g(1) * (rSges(5,1) * t63 - rSges(5,2) * t62 + rSges(5,3) * t72 + t101) + g(2) * (rSges(5,1) * t61 - rSges(5,2) * t60 + rSges(5,3) * t70 + t100) + g(3) * (rSges(5,1) * t67 - rSges(5,2) * t66 - rSges(5,3) * t110 + t102)) - m(6) * (g(1) * (t63 * pkin(4) + (t63 * t95 + t72 * t91) * rSges(6,1) + (-t63 * t91 + t72 * t95) * rSges(6,2) + t113 * t62 + t101) + g(2) * (t61 * pkin(4) + (t61 * t95 + t70 * t91) * rSges(6,1) + (-t61 * t91 + t70 * t95) * rSges(6,2) + t113 * t60 + t100) + g(3) * (t67 * pkin(4) + (-t91 * t110 + t67 * t95) * rSges(6,1) + (-t95 * t110 - t67 * t91) * rSges(6,2) + t113 * t66 + t102));
U = t1;
