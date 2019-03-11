% Calculate potential energy for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:12:57
% EndTime: 2019-03-09 09:12:57
% DurationCPUTime: 0.52s
% Computational Cost: add. (208->105), mult. (290->127), div. (0->0), fcn. (293->10), ass. (0->44)
t114 = pkin(9) + rSges(7,3);
t90 = sin(qJ(2));
t116 = t90 * pkin(2) + pkin(6);
t88 = -pkin(8) - qJ(4);
t115 = -pkin(7) - t88;
t85 = pkin(10) + qJ(5);
t80 = cos(t85);
t113 = t80 * t90;
t86 = sin(pkin(10));
t112 = t86 * t90;
t93 = cos(qJ(2));
t111 = t86 * t93;
t91 = sin(qJ(1));
t110 = t90 * t91;
t109 = t91 * t93;
t94 = cos(qJ(1));
t108 = t93 * t94;
t107 = t94 * pkin(1) + t91 * pkin(7);
t106 = qJ(3) * t90;
t105 = -qJ(4) - rSges(5,3);
t104 = pkin(4) * t112;
t83 = t91 * pkin(1);
t103 = pkin(2) * t109 + t91 * t106 + t83;
t102 = pkin(2) * t108 + t94 * t106 + t107;
t101 = -t93 * qJ(3) + t116;
t87 = cos(pkin(10));
t77 = pkin(4) * t87 + pkin(3);
t100 = t91 * t104 + t77 * t109 + t103;
t79 = sin(t85);
t64 = t79 * t90 + t80 * t93;
t99 = t87 * t90 - t111;
t98 = t87 * t93 + t112;
t89 = sin(qJ(6));
t92 = cos(qJ(6));
t97 = rSges(7,1) * t92 - rSges(7,2) * t89 + pkin(5);
t96 = t94 * t104 + t77 * t108 + t91 * t88 + t102;
t95 = rSges(5,1) * t98 + rSges(5,2) * t99 + t93 * pkin(3);
t72 = t90 * t77;
t65 = -t79 * t93 + t113;
t63 = t64 * t94;
t62 = t108 * t79 - t113 * t94;
t61 = t64 * t91;
t60 = t109 * t79 - t110 * t80;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t94 - t91 * rSges(2,2)) + g(2) * (t91 * rSges(2,1) + rSges(2,2) * t94) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t91 * rSges(3,3) + t107) + g(2) * (rSges(3,1) * t109 - rSges(3,2) * t110 + t83) + g(3) * (rSges(3,1) * t90 + rSges(3,2) * t93 + pkin(6)) + (g(1) * (rSges(3,1) * t93 - rSges(3,2) * t90) + g(2) * (-rSges(3,3) - pkin(7))) * t94) - m(4) * (g(1) * (t91 * rSges(4,2) + t102) + g(2) * (rSges(4,1) * t109 + rSges(4,3) * t110 + t103) + g(3) * (rSges(4,1) * t90 + (-rSges(4,3) - qJ(3)) * t93 + t116) + (g(1) * (rSges(4,1) * t93 + rSges(4,3) * t90) + g(2) * (-rSges(4,2) - pkin(7))) * t94) - m(5) * (g(1) * t102 + g(2) * t103 + g(3) * (rSges(5,1) * t99 - rSges(5,2) * t98 + t90 * pkin(3) + t101) + (g(1) * t105 + g(2) * t95) * t91 + (g(1) * t95 + g(2) * (-pkin(7) - t105)) * t94) - m(6) * (g(1) * (rSges(6,1) * t63 - rSges(6,2) * t62 - rSges(6,3) * t91 + t96) + g(3) * (t65 * rSges(6,1) - t64 * rSges(6,2) + t72 + (-pkin(4) * t86 - qJ(3)) * t93 + t116) + (t61 * rSges(6,1) - t60 * rSges(6,2) + t100 + (rSges(6,3) + t115) * t94) * g(2)) - m(7) * (g(1) * (t63 * pkin(5) + (t63 * t92 - t89 * t91) * rSges(7,1) + (-t63 * t89 - t91 * t92) * rSges(7,2) + t114 * t62 + t96) + g(2) * (t97 * t61 + (t89 * rSges(7,1) + t92 * rSges(7,2) + t115) * t94 + t114 * t60 + t100) + (-pkin(4) * t111 + t114 * t64 + t97 * t65 + t101 + t72) * g(3));
U  = t1;
