% Calculate potential energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:20
% EndTime: 2019-03-09 08:41:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (183->101), mult. (252->122), div. (0->0), fcn. (240->8), ass. (0->39)
t110 = rSges(7,1) + pkin(5);
t109 = rSges(5,3) + qJ(4);
t108 = rSges(7,3) + qJ(6);
t81 = sin(qJ(1));
t83 = cos(qJ(1));
t107 = g(1) * t83 + g(2) * t81;
t80 = sin(qJ(2));
t104 = t80 * pkin(2) + pkin(6);
t103 = t80 * t81;
t102 = t80 * t83;
t76 = pkin(9) + qJ(5);
t71 = cos(t76);
t101 = t81 * t71;
t77 = sin(pkin(9));
t100 = t81 * t77;
t78 = cos(pkin(9));
t99 = t81 * t78;
t82 = cos(qJ(2));
t98 = t81 * t82;
t95 = t83 * pkin(1) + t81 * pkin(7);
t94 = qJ(3) * t80;
t92 = t77 * t102;
t91 = t80 * t100;
t74 = t81 * pkin(1);
t90 = pkin(2) * t98 + t81 * t94 + t74;
t89 = -pkin(4) * t77 - qJ(3);
t88 = t95 + (pkin(2) * t82 + t94) * t83;
t79 = -pkin(8) - qJ(4);
t87 = -t80 * t79 + t104;
t86 = -t83 * pkin(7) + t90;
t69 = pkin(4) * t78 + pkin(3);
t85 = pkin(4) * t92 + t81 * t69 + t88;
t84 = pkin(4) * t91 - t83 * t69 + t86;
t70 = sin(t76);
t60 = t70 * t103 - t71 * t83;
t59 = t80 * t101 + t70 * t83;
t58 = t70 * t102 + t101;
t57 = -t71 * t102 + t70 * t81;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t83 - t81 * rSges(2,2)) + g(2) * (t81 * rSges(2,1) + rSges(2,2) * t83) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t81 * rSges(3,3) + t95) + g(2) * (rSges(3,1) * t98 - rSges(3,2) * t103 + t74) + g(3) * (rSges(3,1) * t80 + rSges(3,2) * t82 + pkin(6)) + (g(1) * (rSges(3,1) * t82 - rSges(3,2) * t80) + g(2) * (-rSges(3,3) - pkin(7))) * t83) - m(4) * (g(1) * (t81 * rSges(4,1) + t88) + g(2) * (-rSges(4,2) * t98 + rSges(4,3) * t103 + t90) + g(3) * (-rSges(4,2) * t80 + (-rSges(4,3) - qJ(3)) * t82 + t104) + (g(1) * (-rSges(4,2) * t82 + rSges(4,3) * t80) + g(2) * (-rSges(4,1) - pkin(7))) * t83) - m(5) * (g(1) * (t81 * pkin(3) + (t92 + t99) * rSges(5,1) + (t78 * t102 - t100) * rSges(5,2) + t88) + g(2) * (-t83 * pkin(3) + (-t78 * t83 + t91) * rSges(5,1) + (t77 * t83 + t80 * t99) * rSges(5,2) + t86) + g(3) * (t109 * t80 + t104) + (g(3) * (-rSges(5,1) * t77 - rSges(5,2) * t78 - qJ(3)) + t107 * t109) * t82) - m(6) * (g(1) * (t58 * rSges(6,1) - t57 * rSges(6,2) + t85) + g(2) * (t60 * rSges(6,1) + t59 * rSges(6,2) + t84) + g(3) * (t80 * rSges(6,3) + t87) + (g(3) * (-rSges(6,1) * t70 - rSges(6,2) * t71 + t89) + t107 * (rSges(6,3) - t79)) * t82) - m(7) * (g(1) * (t108 * t57 + t110 * t58 + t85) + g(2) * (-t108 * t59 + t110 * t60 + t84) + g(3) * (rSges(7,2) * t80 + t87) + (g(3) * (t108 * t71 - t110 * t70 + t89) + t107 * (rSges(7,2) - t79)) * t82);
U  = t1;
