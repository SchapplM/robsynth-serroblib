% Calculate potential energy for
% S6RRPRRP7
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
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:39
% EndTime: 2019-03-09 12:16:40
% DurationCPUTime: 0.33s
% Computational Cost: add. (186->98), mult. (351->120), div. (0->0), fcn. (377->8), ass. (0->39)
t111 = rSges(7,2) + pkin(9);
t112 = rSges(7,1) + pkin(5);
t110 = rSges(6,3) + pkin(9);
t88 = sin(qJ(2));
t109 = t88 * pkin(2) + pkin(6);
t89 = sin(qJ(1));
t108 = t88 * t89;
t91 = cos(qJ(4));
t107 = t88 * t91;
t92 = cos(qJ(2));
t106 = t89 * t92;
t93 = cos(qJ(1));
t105 = t92 * t93;
t104 = t93 * pkin(1) + t89 * pkin(7);
t103 = qJ(3) * t88;
t102 = rSges(7,3) + qJ(6);
t83 = t89 * pkin(1);
t101 = pkin(2) * t106 + t89 * t103 + t83;
t100 = pkin(2) * t105 + t93 * t103 + t104;
t99 = pkin(3) * t106 + t93 * pkin(8) + t101;
t98 = pkin(3) * t105 + t100;
t87 = sin(qJ(4));
t69 = t87 * t88 + t91 * t92;
t97 = t88 * pkin(3) - qJ(3) * t92 + t109;
t70 = -t87 * t92 + t107;
t96 = t70 * pkin(4) + t97;
t65 = t69 * t89;
t95 = t65 * pkin(4) - pkin(7) * t93 + t99;
t67 = t69 * t93;
t94 = t67 * pkin(4) - pkin(8) * t89 + t98;
t90 = cos(qJ(5));
t86 = sin(qJ(5));
t66 = t87 * t105 - t93 * t107;
t64 = t87 * t106 - t89 * t107;
t61 = t67 * t90 - t86 * t89;
t60 = t67 * t86 + t89 * t90;
t59 = t65 * t90 + t86 * t93;
t58 = t65 * t86 - t93 * t90;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t93 - t89 * rSges(2,2)) + g(2) * (t89 * rSges(2,1) + rSges(2,2) * t93) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t89 * rSges(3,3) + t104) + g(2) * (rSges(3,1) * t106 - rSges(3,2) * t108 + t83) + g(3) * (rSges(3,1) * t88 + rSges(3,2) * t92 + pkin(6)) + (g(1) * (rSges(3,1) * t92 - rSges(3,2) * t88) + g(2) * (-rSges(3,3) - pkin(7))) * t93) - m(4) * (g(1) * (t89 * rSges(4,2) + t100) + g(2) * (rSges(4,1) * t106 + rSges(4,3) * t108 + t101) + g(3) * (rSges(4,1) * t88 + (-rSges(4,3) - qJ(3)) * t92 + t109) + (g(1) * (rSges(4,1) * t92 + rSges(4,3) * t88) + g(2) * (-rSges(4,2) - pkin(7))) * t93) - m(5) * (g(1) * (rSges(5,1) * t67 - rSges(5,2) * t66 + (-rSges(5,3) - pkin(8)) * t89 + t98) + g(2) * (t65 * rSges(5,1) - t64 * rSges(5,2) + (rSges(5,3) - pkin(7)) * t93 + t99) + g(3) * (rSges(5,1) * t70 - rSges(5,2) * t69 + t97)) - m(6) * (g(1) * (rSges(6,1) * t61 - rSges(6,2) * t60 + t110 * t66 + t94) + g(2) * (t59 * rSges(6,1) - t58 * rSges(6,2) + t110 * t64 + t95) + g(3) * ((rSges(6,1) * t90 - rSges(6,2) * t86) * t70 + t110 * t69 + t96)) - m(7) * (g(1) * (t102 * t60 + t111 * t66 + t112 * t61 + t94) + g(2) * (t102 * t58 + t111 * t64 + t112 * t59 + t95) + (t96 + (t102 * t86 + t112 * t90) * t70 + t111 * t69) * g(3));
U  = t1;
