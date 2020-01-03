% Calculate potential energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:17
% EndTime: 2019-12-31 19:23:17
% DurationCPUTime: 0.33s
% Computational Cost: add. (168->86), mult. (367->107), div. (0->0), fcn. (406->8), ass. (0->43)
t104 = rSges(6,2) + qJ(4);
t102 = rSges(6,3) + qJ(5);
t90 = cos(pkin(5));
t92 = sin(qJ(1));
t110 = t92 * t90;
t88 = sin(pkin(5));
t94 = cos(qJ(1));
t115 = t88 * t94;
t91 = sin(qJ(2));
t120 = t91 * t110 + t115;
t119 = rSges(6,1) + pkin(4);
t118 = t91 * pkin(2) + pkin(6);
t117 = t88 * t92;
t93 = cos(qJ(2));
t116 = t88 * t93;
t114 = t90 * t93;
t87 = sin(pkin(8));
t113 = t91 * t87;
t89 = cos(pkin(8));
t112 = t91 * t89;
t111 = t91 * t92;
t109 = t92 * t93;
t108 = t93 * t94;
t107 = t94 * t90;
t106 = t94 * pkin(1) + t92 * pkin(7);
t105 = qJ(3) * t90;
t103 = rSges(5,3) + qJ(4);
t70 = t114 * t87 + t112;
t101 = t70 * pkin(3) + t118;
t99 = qJ(3) * t88 * t91;
t98 = pkin(2) * t108 + t92 * t105 + t94 * t99 + t106;
t67 = -t107 * t113 + t108 * t89 + t117 * t87;
t97 = t67 * pkin(3) + t98;
t85 = t92 * pkin(1);
t96 = t92 * t99 + pkin(2) * t109 + t85 + (-pkin(7) - t105) * t94;
t65 = t109 * t89 - t120 * t87;
t95 = t65 * pkin(3) + t96;
t72 = t115 * t91 + t110;
t71 = t111 * t88 - t107;
t69 = -t114 * t89 + t113;
t66 = -t89 * t117 + (t112 * t90 + t87 * t93) * t94;
t64 = t87 * t109 + t120 * t89;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t94 * rSges(2,1) - t92 * rSges(2,2)) + g(2) * (t92 * rSges(2,1) + t94 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t92 * rSges(3,3) + t106) + g(2) * (rSges(3,1) * t109 - rSges(3,2) * t111 + t85) + g(3) * (t91 * rSges(3,1) + t93 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t93 - rSges(3,2) * t91) + g(2) * (-rSges(3,3) - pkin(7))) * t94) - m(4) * (g(1) * (t67 * rSges(4,1) - t66 * rSges(4,2) + t72 * rSges(4,3) + t98) + g(2) * (t65 * rSges(4,1) - t64 * rSges(4,2) + t71 * rSges(4,3) + t96) + g(3) * (t70 * rSges(4,1) - t69 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t116 + t118)) - m(5) * (g(1) * (t72 * rSges(5,1) - t67 * rSges(5,2) + t103 * t66 + t97) + g(2) * (t71 * rSges(5,1) - t65 * rSges(5,2) + t103 * t64 + t95) + g(3) * (-t70 * rSges(5,2) + (-rSges(5,1) - qJ(3)) * t116 + t103 * t69 + t101)) - m(6) * (g(1) * (t102 * t67 + t104 * t66 + t119 * t72 + t97) + g(2) * (t102 * t65 + t104 * t64 + t119 * t71 + t95) + (t101 + (-qJ(3) - t119) * t116 + t102 * t70 + t104 * t69) * g(3));
U = t1;
