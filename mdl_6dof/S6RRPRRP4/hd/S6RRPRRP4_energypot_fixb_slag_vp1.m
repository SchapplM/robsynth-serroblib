% Calculate potential energy for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:28
% EndTime: 2019-03-09 11:52:29
% DurationCPUTime: 0.44s
% Computational Cost: add. (239->95), mult. (227->115), div. (0->0), fcn. (215->10), ass. (0->45)
t118 = rSges(7,1) + pkin(5);
t117 = rSges(5,3) + pkin(8);
t116 = rSges(7,3) + qJ(6);
t86 = sin(qJ(1));
t89 = cos(qJ(1));
t115 = g(1) * t89 + g(2) * t86;
t112 = rSges(3,3) + pkin(7);
t85 = sin(qJ(2));
t110 = t85 * pkin(2) + pkin(6);
t81 = qJ(2) + pkin(10);
t76 = sin(t81);
t109 = rSges(4,2) * t76;
t77 = cos(t81);
t108 = t86 * t77;
t82 = qJ(4) + qJ(5);
t78 = sin(t82);
t107 = t86 * t78;
t79 = cos(t82);
t106 = t86 * t79;
t84 = sin(qJ(4));
t105 = t86 * t84;
t87 = cos(qJ(4));
t104 = t86 * t87;
t103 = t89 * t77;
t102 = t89 * t78;
t101 = t89 * t79;
t100 = t89 * t84;
t99 = t89 * t87;
t88 = cos(qJ(2));
t75 = t88 * pkin(2) + pkin(1);
t83 = -qJ(3) - pkin(7);
t96 = t86 * t75 + t89 * t83;
t69 = t89 * t75;
t95 = -t86 * t83 + t69;
t74 = t87 * pkin(4) + pkin(3);
t90 = -pkin(9) - pkin(8);
t94 = t76 * t74 + t77 * t90 + t110;
t93 = pkin(4) * t105 + t74 * t103 + t95;
t92 = rSges(3,1) * t88 - rSges(3,2) * t85 + pkin(1);
t91 = -pkin(4) * t100 + t74 * t108 + t96;
t64 = t77 * t101 + t107;
t63 = t77 * t102 - t106;
t62 = t77 * t106 - t102;
t61 = t77 * t107 + t101;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t89 * rSges(2,1) - t86 * rSges(2,2)) + g(2) * (t86 * rSges(2,1) + t89 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t85 * rSges(3,1) + t88 * rSges(3,2) + pkin(6)) + (g(1) * t92 - g(2) * t112) * t89 + (g(1) * t112 + g(2) * t92) * t86) - m(4) * (g(1) * (rSges(4,1) * t103 - t89 * t109 + t69) + g(2) * (-t89 * rSges(4,3) + t96) + g(3) * (t76 * rSges(4,1) + t77 * rSges(4,2) + t110) + (g(1) * (rSges(4,3) - t83) + g(2) * (rSges(4,1) * t77 - t109)) * t86) - m(5) * (g(1) * (pkin(3) * t103 + (t77 * t99 + t105) * rSges(5,1) + (-t77 * t100 + t104) * rSges(5,2) + t95) + g(2) * (pkin(3) * t108 + (t77 * t104 - t100) * rSges(5,1) + (-t77 * t105 - t99) * rSges(5,2) + t96) + g(3) * (-t117 * t77 + t110) + (g(3) * (rSges(5,1) * t87 - rSges(5,2) * t84 + pkin(3)) + t115 * t117) * t76) - m(6) * (g(1) * (t64 * rSges(6,1) - t63 * rSges(6,2) + t93) + g(2) * (t62 * rSges(6,1) - t61 * rSges(6,2) + t91) + g(3) * (-t77 * rSges(6,3) + t94) + (g(3) * (rSges(6,1) * t79 - rSges(6,2) * t78) + t115 * (rSges(6,3) - t90)) * t76) - m(7) * (g(1) * (t116 * t63 + t118 * t64 + t93) + g(2) * (t116 * t61 + t118 * t62 + t91) + g(3) * (-t77 * rSges(7,2) + t94) + (g(3) * (t116 * t78 + t118 * t79) + t115 * (rSges(7,2) - t90)) * t76);
U  = t1;
