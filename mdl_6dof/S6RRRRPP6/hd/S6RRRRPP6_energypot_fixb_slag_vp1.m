% Calculate potential energy for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:04
% EndTime: 2019-03-09 21:11:04
% DurationCPUTime: 0.48s
% Computational Cost: add. (224->99), mult. (280->122), div. (0->0), fcn. (278->8), ass. (0->35)
t108 = -rSges(7,1) - pkin(5);
t107 = rSges(4,3) + pkin(8);
t106 = rSges(7,3) + qJ(6);
t80 = sin(qJ(1));
t83 = cos(qJ(1));
t105 = g(1) * t83 + g(2) * t80;
t79 = sin(qJ(2));
t101 = rSges(3,2) * t79;
t78 = sin(qJ(3));
t100 = t78 * t80;
t99 = t78 * t83;
t82 = cos(qJ(2));
t98 = t80 * t82;
t97 = t82 * t83;
t77 = qJ(3) + qJ(4);
t72 = sin(t77);
t96 = t83 * t72;
t93 = pkin(1) * t83 + pkin(7) * t80;
t81 = cos(qJ(3));
t70 = pkin(3) * t81 + pkin(2);
t84 = -pkin(9) - pkin(8);
t91 = t70 * t79 + t82 * t84 + pkin(6);
t75 = t80 * pkin(1);
t90 = -t83 * pkin(7) + t75;
t89 = pkin(3) * t100 + t70 * t97 + t93;
t73 = cos(t77);
t88 = t91 + (pkin(4) * t73 + qJ(5) * t72) * t79;
t60 = -t73 * t80 + t82 * t96;
t61 = t72 * t80 + t73 * t97;
t87 = pkin(4) * t61 + t60 * qJ(5) + t89;
t86 = -pkin(3) * t99 + t70 * t98 + t90;
t58 = t72 * t98 + t73 * t83;
t59 = t73 * t98 - t96;
t85 = pkin(4) * t59 + t58 * qJ(5) + t86;
t1 = -m(1) * (rSges(1,1) * g(1) + g(2) * rSges(1,2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t83 - rSges(2,2) * t80) + g(2) * (rSges(2,1) * t80 + rSges(2,2) * t83) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t80 + t93) + g(2) * (rSges(3,1) * t98 - t80 * t101 + t75) + g(3) * (rSges(3,1) * t79 + rSges(3,2) * t82 + pkin(6)) + (g(1) * (rSges(3,1) * t82 - t101) + g(2) * (-rSges(3,3) - pkin(7))) * t83) - m(4) * (g(1) * (pkin(2) * t97 + (t81 * t97 + t100) * rSges(4,1) + (-t78 * t97 + t80 * t81) * rSges(4,2) + t93) + g(2) * (pkin(2) * t98 + (t81 * t98 - t99) * rSges(4,1) + (-t78 * t98 - t81 * t83) * rSges(4,2) + t90) + g(3) * (-t107 * t82 + pkin(6)) + (g(3) * (rSges(4,1) * t81 - rSges(4,2) * t78 + pkin(2)) + t105 * t107) * t79) - m(5) * (g(1) * (t61 * rSges(5,1) - t60 * rSges(5,2) + t89) + g(2) * (t59 * rSges(5,1) - t58 * rSges(5,2) + t86) + g(3) * (-rSges(5,3) * t82 + t91) + (g(3) * (rSges(5,1) * t73 - rSges(5,2) * t72) + t105 * (rSges(5,3) - t84)) * t79) - m(6) * (g(1) * (-t61 * rSges(6,2) + t60 * rSges(6,3) + t87) + g(2) * (-t59 * rSges(6,2) + t58 * rSges(6,3) + t85) + g(3) * (-rSges(6,1) * t82 + t88) + (g(3) * (-rSges(6,2) * t73 + rSges(6,3) * t72) + t105 * (rSges(6,1) - t84)) * t79) - m(7) * (g(1) * (t60 * rSges(7,2) + t106 * t61 + t87) + g(2) * (t58 * rSges(7,2) + t106 * t59 + t85) + g(3) * (t108 * t82 + t88) + (g(3) * (rSges(7,2) * t72 + t106 * t73) + t105 * (-t84 - t108)) * t79);
U  = t1;
