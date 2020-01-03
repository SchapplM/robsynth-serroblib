% Calculate potential energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:04
% EndTime: 2019-12-31 22:08:04
% DurationCPUTime: 0.29s
% Computational Cost: add. (207->96), mult. (435->127), div. (0->0), fcn. (511->10), ass. (0->47)
t114 = rSges(4,3) + pkin(8);
t113 = rSges(5,3) + pkin(9);
t87 = cos(pkin(5));
t112 = t87 * pkin(7) + pkin(6);
t86 = sin(pkin(5));
t91 = sin(qJ(2));
t111 = t86 * t91;
t92 = sin(qJ(1));
t110 = t86 * t92;
t94 = cos(qJ(3));
t109 = t86 * t94;
t95 = cos(qJ(2));
t108 = t86 * t95;
t96 = cos(qJ(1));
t107 = t86 * t96;
t106 = t92 * t91;
t105 = t92 * t95;
t104 = t96 * t91;
t103 = t96 * t95;
t102 = rSges(6,3) + qJ(5) + pkin(9);
t101 = t96 * pkin(1) + pkin(7) * t110;
t100 = pkin(2) * t111 + t112;
t76 = -t87 * t106 + t103;
t99 = t76 * pkin(2) + t101;
t89 = sin(qJ(4));
t98 = pkin(4) * t89 + pkin(8);
t74 = t87 * t104 + t105;
t84 = t92 * pkin(1);
t97 = t74 * pkin(2) - pkin(7) * t107 + t84;
t93 = cos(qJ(4));
t90 = sin(qJ(3));
t82 = t93 * pkin(4) + pkin(3);
t75 = t87 * t105 + t104;
t73 = -t87 * t103 + t106;
t72 = t91 * t109 + t87 * t90;
t71 = t90 * t111 - t87 * t94;
t68 = t90 * t110 + t76 * t94;
t67 = -t92 * t109 + t76 * t90;
t66 = -t90 * t107 + t74 * t94;
t65 = t94 * t107 + t74 * t90;
t64 = -t89 * t108 + t72 * t93;
t63 = -t93 * t108 - t72 * t89;
t62 = t68 * t93 + t75 * t89;
t61 = -t68 * t89 + t75 * t93;
t60 = t66 * t93 + t73 * t89;
t59 = -t66 * t89 + t73 * t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t96 * rSges(2,1) - t92 * rSges(2,2)) + g(2) * (t92 * rSges(2,1) + t96 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t76 * rSges(3,1) - t75 * rSges(3,2) + t101) + g(2) * (t74 * rSges(3,1) - t73 * rSges(3,2) + t84) + g(3) * (t87 * rSges(3,3) + t112) + (g(1) * rSges(3,3) * t92 + g(3) * (rSges(3,1) * t91 + rSges(3,2) * t95) + g(2) * (-rSges(3,3) - pkin(7)) * t96) * t86) - m(4) * (g(1) * (t68 * rSges(4,1) - t67 * rSges(4,2) + t114 * t75 + t99) + g(2) * (t66 * rSges(4,1) - t65 * rSges(4,2) + t114 * t73 + t97) + g(3) * (t72 * rSges(4,1) - t71 * rSges(4,2) - t114 * t108 + t100)) - m(5) * (g(1) * (t62 * rSges(5,1) + t61 * rSges(5,2) + t68 * pkin(3) + t75 * pkin(8) + t113 * t67 + t99) + g(2) * (t60 * rSges(5,1) + t59 * rSges(5,2) + t66 * pkin(3) + t73 * pkin(8) + t113 * t65 + t97) + g(3) * (t64 * rSges(5,1) + t63 * rSges(5,2) + t72 * pkin(3) - pkin(8) * t108 + t113 * t71 + t100)) - m(6) * (g(1) * (t62 * rSges(6,1) + t61 * rSges(6,2) + t102 * t67 + t68 * t82 + t98 * t75 + t99) + g(2) * (t60 * rSges(6,1) + t59 * rSges(6,2) + t102 * t65 + t66 * t82 + t98 * t73 + t97) + g(3) * (t64 * rSges(6,1) + t63 * rSges(6,2) + t102 * t71 - t98 * t108 + t72 * t82 + t100));
U = t1;
