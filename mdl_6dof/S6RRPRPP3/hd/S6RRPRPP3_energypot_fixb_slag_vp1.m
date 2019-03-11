% Calculate potential energy for
% S6RRPRPP3
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:09
% EndTime: 2019-03-09 09:54:09
% DurationCPUTime: 0.48s
% Computational Cost: add. (224->99), mult. (280->122), div. (0->0), fcn. (278->8), ass. (0->34)
t107 = -rSges(7,1) - pkin(5);
t106 = rSges(4,3) + qJ(3);
t105 = rSges(7,3) + qJ(6);
t82 = sin(qJ(1));
t84 = cos(qJ(1));
t104 = g(1) * t84 + g(2) * t82;
t81 = sin(qJ(2));
t101 = rSges(3,2) * t81;
t78 = sin(pkin(9));
t100 = t78 * t84;
t99 = t82 * t78;
t83 = cos(qJ(2));
t98 = t82 * t83;
t97 = t84 * t83;
t94 = t84 * pkin(1) + t82 * pkin(7);
t79 = cos(pkin(9));
t70 = pkin(3) * t79 + pkin(2);
t80 = -pkin(8) - qJ(3);
t91 = t81 * t70 + t83 * t80 + pkin(6);
t75 = t82 * pkin(1);
t90 = -t84 * pkin(7) + t75;
t89 = pkin(3) * t99 + t70 * t97 + t94;
t77 = pkin(9) + qJ(4);
t72 = sin(t77);
t73 = cos(t77);
t88 = t91 + (pkin(4) * t73 + qJ(5) * t72) * t81;
t60 = t72 * t97 - t82 * t73;
t61 = t82 * t72 + t73 * t97;
t87 = t61 * pkin(4) + t60 * qJ(5) + t89;
t86 = -pkin(3) * t100 + t70 * t98 + t90;
t58 = t72 * t98 + t73 * t84;
t59 = -t84 * t72 + t73 * t98;
t85 = t59 * pkin(4) + t58 * qJ(5) + t86;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t84 - t82 * rSges(2,2)) + g(2) * (t82 * rSges(2,1) + rSges(2,2) * t84) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t82 * rSges(3,3) + t94) + g(2) * (rSges(3,1) * t98 - t82 * t101 + t75) + g(3) * (rSges(3,1) * t81 + rSges(3,2) * t83 + pkin(6)) + (g(1) * (rSges(3,1) * t83 - t101) + g(2) * (-rSges(3,3) - pkin(7))) * t84) - m(4) * (g(1) * (pkin(2) * t97 + (t79 * t97 + t99) * rSges(4,1) + (-t78 * t97 + t82 * t79) * rSges(4,2) + t94) + g(2) * (pkin(2) * t98 + (t79 * t98 - t100) * rSges(4,1) + (-t78 * t98 - t79 * t84) * rSges(4,2) + t90) + g(3) * (-t106 * t83 + pkin(6)) + (g(3) * (rSges(4,1) * t79 - rSges(4,2) * t78 + pkin(2)) + t104 * t106) * t81) - m(5) * (g(1) * (t61 * rSges(5,1) - t60 * rSges(5,2) + t89) + g(2) * (t59 * rSges(5,1) - t58 * rSges(5,2) + t86) + g(3) * (-rSges(5,3) * t83 + t91) + (g(3) * (rSges(5,1) * t73 - rSges(5,2) * t72) + t104 * (rSges(5,3) - t80)) * t81) - m(6) * (g(1) * (-t61 * rSges(6,2) + t60 * rSges(6,3) + t87) + g(2) * (-t59 * rSges(6,2) + t58 * rSges(6,3) + t85) + g(3) * (-rSges(6,1) * t83 + t88) + (g(3) * (-rSges(6,2) * t73 + rSges(6,3) * t72) + t104 * (rSges(6,1) - t80)) * t81) - m(7) * (g(1) * (t60 * rSges(7,2) + t105 * t61 + t87) + g(2) * (t58 * rSges(7,2) + t105 * t59 + t85) + g(3) * (t107 * t83 + t88) + (g(3) * (rSges(7,2) * t72 + t105 * t73) + t104 * (-t80 - t107)) * t81);
U  = t1;
