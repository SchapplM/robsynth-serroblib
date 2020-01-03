% Calculate potential energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR16_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR16_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:37
% EndTime: 2019-12-31 20:44:38
% DurationCPUTime: 0.39s
% Computational Cost: add. (175->97), mult. (358->130), div. (0->0), fcn. (405->10), ass. (0->43)
t85 = sin(qJ(1));
t110 = g(1) * t85;
t89 = cos(qJ(1));
t109 = g(2) * t89;
t81 = cos(pkin(5));
t108 = t81 * pkin(7) + pkin(6);
t107 = pkin(9) + rSges(6,3);
t80 = sin(pkin(5));
t84 = sin(qJ(2));
t106 = t80 * t84;
t105 = t80 * t85;
t88 = cos(qJ(2));
t104 = t80 * t88;
t103 = t80 * t89;
t102 = t84 * t85;
t101 = t84 * t89;
t100 = t85 * t88;
t99 = t88 * t89;
t98 = t89 * pkin(1) + pkin(7) * t105;
t97 = t88 * qJ(3);
t96 = pkin(2) * t106 + t108;
t95 = (-pkin(3) - pkin(7)) * t89;
t94 = t81 * pkin(3) + pkin(8) * t106 + t96;
t66 = -t81 * t99 + t102;
t67 = t81 * t101 + t100;
t78 = t85 * pkin(1);
t93 = t67 * pkin(2) + t66 * qJ(3) + t78;
t68 = t81 * t100 + t101;
t69 = -t81 * t102 + t99;
t92 = t69 * pkin(2) + t68 * qJ(3) + t98;
t91 = pkin(3) * t105 + t92;
t90 = t67 * pkin(8) + t93;
t87 = cos(qJ(4));
t86 = cos(qJ(5));
t83 = sin(qJ(4));
t82 = sin(qJ(5));
t65 = -t83 * t104 + t81 * t87;
t64 = t87 * t104 + t81 * t83;
t60 = -t87 * t103 + t66 * t83;
t59 = t83 * t103 + t66 * t87;
t58 = t87 * t105 + t68 * t83;
t57 = t83 * t105 - t68 * t87;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t89 - t85 * rSges(2,2)) + g(2) * (t85 * rSges(2,1) + rSges(2,2) * t89) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t69 - rSges(3,2) * t68 + t98) + g(2) * (t67 * rSges(3,1) - t66 * rSges(3,2) + t78) + g(3) * (rSges(3,3) * t81 + t108) + (rSges(3,3) * t110 + g(3) * (rSges(3,1) * t84 + rSges(3,2) * t88) + (-rSges(3,3) - pkin(7)) * t109) * t80) - m(4) * (g(1) * (-rSges(4,2) * t69 + rSges(4,3) * t68 + t92) + g(2) * (-t67 * rSges(4,2) + t66 * rSges(4,3) + t93) + g(3) * (rSges(4,1) * t81 + t96) + (rSges(4,1) * t110 + g(3) * (-rSges(4,2) * t84 - rSges(4,3) * t88 - t97) + (-rSges(4,1) - pkin(7)) * t109) * t80) - m(5) * (g(1) * (rSges(5,1) * t58 - rSges(5,2) * t57 + (rSges(5,3) + pkin(8)) * t69 + t91) + g(2) * (t60 * rSges(5,1) + t59 * rSges(5,2) + t67 * rSges(5,3) + t90) + g(3) * (rSges(5,1) * t65 - rSges(5,2) * t64 + t94) + (g(3) * (rSges(5,3) * t84 - t97) + g(2) * t95) * t80) - m(6) * (g(1) * (t58 * pkin(4) + t69 * pkin(8) + (t58 * t86 + t69 * t82) * rSges(6,1) + (-t58 * t82 + t69 * t86) * rSges(6,2) + t107 * t57 + t91) + g(2) * (t60 * pkin(4) + (t60 * t86 + t67 * t82) * rSges(6,1) + (-t60 * t82 + t67 * t86) * rSges(6,2) + t80 * t95 - t107 * t59 + t90) + g(3) * (t65 * pkin(4) - t80 * t97 + (t82 * t106 + t65 * t86) * rSges(6,1) + (t86 * t106 - t65 * t82) * rSges(6,2) + t107 * t64 + t94));
U = t1;
