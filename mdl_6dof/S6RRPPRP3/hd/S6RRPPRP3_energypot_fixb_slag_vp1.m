% Calculate potential energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:29
% EndTime: 2019-03-09 08:33:29
% DurationCPUTime: 0.42s
% Computational Cost: add. (142->98), mult. (226->118), div. (0->0), fcn. (206->6), ass. (0->35)
t102 = rSges(6,3) + pkin(8);
t101 = rSges(7,3) + qJ(6) + pkin(8);
t73 = sin(qJ(1));
t76 = cos(qJ(1));
t100 = g(1) * t76 + g(2) * t73;
t72 = sin(qJ(2));
t96 = t72 * pkin(2) + pkin(6);
t95 = t72 * t73;
t94 = t72 * t76;
t71 = sin(qJ(5));
t93 = t73 * t71;
t74 = cos(qJ(5));
t92 = t73 * t74;
t75 = cos(qJ(2));
t91 = t73 * t75;
t90 = t75 * t76;
t89 = t76 * t71;
t88 = t76 * t74;
t86 = t76 * pkin(1) + t73 * pkin(7);
t85 = qJ(3) * t72;
t84 = t72 * pkin(3) + t96;
t68 = t73 * pkin(1);
t83 = pkin(2) * t91 + t73 * t85 + t68;
t82 = pkin(2) * t90 + t76 * t85 + t86;
t81 = pkin(3) * t91 + t76 * qJ(4) + t83;
t80 = pkin(3) * t90 + t82;
t79 = rSges(5,1) * t72 - rSges(5,2) * t75;
t78 = -t76 * pkin(7) + t81;
t77 = -t73 * qJ(4) + t80;
t63 = t74 * pkin(5) + pkin(4);
t56 = t72 * t88 - t93;
t55 = -t72 * t89 - t92;
t54 = t72 * t92 + t89;
t53 = -t72 * t93 + t88;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t76 * rSges(2,1) - t73 * rSges(2,2)) + g(2) * (t73 * rSges(2,1) + t76 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t73 * rSges(3,3) + t86) + g(2) * (rSges(3,1) * t91 - rSges(3,2) * t95 + t68) + g(3) * (t72 * rSges(3,1) + t75 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t75 - rSges(3,2) * t72) + g(2) * (-rSges(3,3) - pkin(7))) * t76) - m(4) * (g(1) * (t73 * rSges(4,2) + t82) + g(2) * (rSges(4,1) * t91 + rSges(4,3) * t95 + t83) + g(3) * (t72 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t75 + t96) + (g(1) * (rSges(4,1) * t75 + rSges(4,3) * t72) + g(2) * (-rSges(4,2) - pkin(7))) * t76) - m(5) * (g(1) * t80 + g(2) * t81 + g(3) * (-t72 * rSges(5,2) + (-rSges(5,1) - qJ(3)) * t75 + t84) + (g(1) * t79 + g(2) * (rSges(5,3) - pkin(7))) * t76 + (g(1) * (-rSges(5,3) - qJ(4)) + g(2) * t79) * t73) - m(6) * (g(1) * (t56 * rSges(6,1) + t55 * rSges(6,2) + pkin(4) * t94 + t77) + g(2) * (t54 * rSges(6,1) + t53 * rSges(6,2) + pkin(4) * t95 + t78) + g(3) * (t102 * t72 + t84) + (g(3) * (-rSges(6,1) * t74 + rSges(6,2) * t71 - pkin(4) - qJ(3)) + t100 * t102) * t75) - m(7) * (g(1) * (t56 * rSges(7,1) + t55 * rSges(7,2) - pkin(5) * t93 + t63 * t94 + t77) + g(2) * (t54 * rSges(7,1) + t53 * rSges(7,2) + pkin(5) * t89 + t63 * t95 + t78) + g(3) * (t101 * t72 + t84) + (g(3) * (-rSges(7,1) * t74 + rSges(7,2) * t71 - qJ(3) - t63) + t100 * t101) * t75);
U  = t1;
