% Calculate potential energy for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:11
% EndTime: 2019-03-09 10:45:12
% DurationCPUTime: 0.49s
% Computational Cost: add. (208->105), mult. (290->127), div. (0->0), fcn. (293->10), ass. (0->44)
t113 = pkin(9) + rSges(7,3);
t89 = sin(qJ(2));
t116 = t89 * pkin(2) + pkin(6);
t86 = -qJ(5) - pkin(8);
t115 = -pkin(7) - t86;
t114 = -pkin(8) - rSges(5,3);
t85 = qJ(4) + pkin(10);
t80 = cos(t85);
t112 = t80 * t89;
t88 = sin(qJ(4));
t111 = t88 * t89;
t93 = cos(qJ(2));
t110 = t88 * t93;
t90 = sin(qJ(1));
t109 = t89 * t90;
t108 = t90 * t93;
t94 = cos(qJ(1));
t107 = t93 * t94;
t106 = t94 * pkin(1) + t90 * pkin(7);
t105 = qJ(3) * t89;
t104 = pkin(4) * t111;
t83 = t90 * pkin(1);
t103 = pkin(2) * t108 + t90 * t105 + t83;
t102 = pkin(2) * t107 + t94 * t105 + t106;
t101 = -t93 * qJ(3) + t116;
t92 = cos(qJ(4));
t78 = pkin(4) * t92 + pkin(3);
t100 = t90 * t104 + t78 * t108 + t103;
t79 = sin(t85);
t64 = t79 * t89 + t80 * t93;
t99 = t89 * t92 - t110;
t98 = t92 * t93 + t111;
t87 = sin(qJ(6));
t91 = cos(qJ(6));
t97 = rSges(7,1) * t91 - rSges(7,2) * t87 + pkin(5);
t96 = t94 * t104 + t78 * t107 + t90 * t86 + t102;
t95 = rSges(5,1) * t98 + rSges(5,2) * t99 + t93 * pkin(3);
t72 = t89 * t78;
t65 = -t79 * t93 + t112;
t63 = t64 * t94;
t62 = t107 * t79 - t112 * t94;
t61 = t64 * t90;
t60 = t108 * t79 - t109 * t80;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t94 - t90 * rSges(2,2)) + g(2) * (t90 * rSges(2,1) + rSges(2,2) * t94) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t90 * rSges(3,3) + t106) + g(2) * (rSges(3,1) * t108 - rSges(3,2) * t109 + t83) + g(3) * (rSges(3,1) * t89 + rSges(3,2) * t93 + pkin(6)) + (g(1) * (rSges(3,1) * t93 - rSges(3,2) * t89) + g(2) * (-rSges(3,3) - pkin(7))) * t94) - m(4) * (g(1) * (t90 * rSges(4,2) + t102) + g(2) * (rSges(4,1) * t108 + rSges(4,3) * t109 + t103) + g(3) * (rSges(4,1) * t89 + (-rSges(4,3) - qJ(3)) * t93 + t116) + (g(1) * (rSges(4,1) * t93 + rSges(4,3) * t89) + g(2) * (-rSges(4,2) - pkin(7))) * t94) - m(5) * (g(1) * t102 + g(2) * t103 + g(3) * (rSges(5,1) * t99 - rSges(5,2) * t98 + t89 * pkin(3) + t101) + (g(1) * t114 + g(2) * t95) * t90 + (g(1) * t95 + g(2) * (-pkin(7) - t114)) * t94) - m(6) * (g(1) * (rSges(6,1) * t63 - rSges(6,2) * t62 - rSges(6,3) * t90 + t96) + g(3) * (rSges(6,1) * t65 - rSges(6,2) * t64 + t72 + (-pkin(4) * t88 - qJ(3)) * t93 + t116) + (t61 * rSges(6,1) - t60 * rSges(6,2) + t100 + (rSges(6,3) + t115) * t94) * g(2)) - m(7) * (g(1) * (t63 * pkin(5) + (t63 * t91 - t87 * t90) * rSges(7,1) + (-t63 * t87 - t90 * t91) * rSges(7,2) + t113 * t62 + t96) + g(2) * (t97 * t61 + (t87 * rSges(7,1) + t91 * rSges(7,2) + t115) * t94 + t113 * t60 + t100) + (-pkin(4) * t110 + t113 * t64 + t97 * t65 + t101 + t72) * g(3));
U  = t1;
