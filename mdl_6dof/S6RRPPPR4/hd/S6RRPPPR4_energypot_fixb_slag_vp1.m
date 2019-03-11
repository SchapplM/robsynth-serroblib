% Calculate potential energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:10
% EndTime: 2019-03-09 08:17:11
% DurationCPUTime: 0.46s
% Computational Cost: add. (165->103), mult. (293->130), div. (0->0), fcn. (297->8), ass. (0->32)
t100 = -pkin(8) - rSges(7,3);
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t84 = g(1) * t79 + g(2) * t76;
t75 = sin(qJ(2));
t97 = t75 * pkin(2) + pkin(6);
t95 = t75 * t76;
t94 = t75 * t79;
t73 = cos(pkin(9));
t93 = t76 * t73;
t78 = cos(qJ(2));
t92 = t76 * t78;
t91 = t79 * pkin(1) + t76 * pkin(7);
t90 = qJ(3) * t75;
t89 = qJ(4) * t78;
t88 = t75 * qJ(4) + t97;
t70 = t76 * pkin(1);
t87 = pkin(2) * t92 + t76 * t90 + t70;
t86 = t78 * t73 * qJ(5) + t88;
t85 = t91 + (pkin(2) * t78 + t90) * t79;
t83 = t76 * pkin(3) + t79 * t89 + t85;
t82 = t76 * t89 + t87 + (-pkin(3) - pkin(7)) * t79;
t72 = sin(pkin(9));
t54 = t72 * t76 - t73 * t94;
t55 = t72 * t94 + t93;
t81 = t55 * pkin(4) + t54 * qJ(5) + t83;
t56 = t72 * t79 + t75 * t93;
t57 = t72 * t95 - t73 * t79;
t80 = t57 * pkin(4) - t56 * qJ(5) + t82;
t77 = cos(qJ(6));
t74 = sin(qJ(6));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - t76 * rSges(2,2)) + g(2) * (t76 * rSges(2,1) + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t76 * rSges(3,3) + t91) + g(2) * (rSges(3,1) * t92 - rSges(3,2) * t95 + t70) + g(3) * (rSges(3,1) * t75 + rSges(3,2) * t78 + pkin(6)) + (g(1) * (rSges(3,1) * t78 - rSges(3,2) * t75) + g(2) * (-rSges(3,3) - pkin(7))) * t79) - m(4) * (g(1) * (t76 * rSges(4,1) + t85) + g(2) * (-rSges(4,2) * t92 + rSges(4,3) * t95 + t87) + g(3) * (-rSges(4,2) * t75 + (-rSges(4,3) - qJ(3)) * t78 + t97) + (g(1) * (-rSges(4,2) * t78 + rSges(4,3) * t75) + g(2) * (-rSges(4,1) - pkin(7))) * t79) - m(5) * (g(1) * (t55 * rSges(5,1) - t54 * rSges(5,2) + t83) + g(2) * (t57 * rSges(5,1) + t56 * rSges(5,2) + t82) + g(3) * (rSges(5,3) * t75 + t88) + (g(3) * (-rSges(5,1) * t72 - rSges(5,2) * t73 - qJ(3)) + t84 * rSges(5,3)) * t78) - m(6) * (g(1) * (t55 * rSges(6,1) + t54 * rSges(6,3) + t81) + g(2) * (t57 * rSges(6,1) - t56 * rSges(6,3) + t80) + g(3) * (rSges(6,2) * t75 + t86) + (g(3) * (rSges(6,3) * t73 - qJ(3) + (-rSges(6,1) - pkin(4)) * t72) + t84 * rSges(6,2)) * t78) - m(7) * (g(1) * (t55 * pkin(5) + (t54 * t74 + t55 * t77) * rSges(7,1) + (t54 * t77 - t55 * t74) * rSges(7,2) + t81) + g(2) * (t57 * pkin(5) + (-t56 * t74 + t57 * t77) * rSges(7,1) + (-t56 * t77 - t57 * t74) * rSges(7,2) + t80) + g(3) * (t100 * t75 + t86) + (g(3) * (-qJ(3) + (t74 * rSges(7,1) + t77 * rSges(7,2)) * t73 + (-t77 * rSges(7,1) + t74 * rSges(7,2) - pkin(4) - pkin(5)) * t72) + t84 * t100) * t78);
U  = t1;
