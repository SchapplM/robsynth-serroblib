% Calculate potential energy for
% S6RRRRPP3
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:38
% EndTime: 2019-03-09 20:53:39
% DurationCPUTime: 0.41s
% Computational Cost: add. (223->92), mult. (248->110), div. (0->0), fcn. (242->8), ass. (0->38)
t104 = rSges(7,2) + qJ(5);
t103 = rSges(7,3) + qJ(6);
t102 = rSges(7,1) + pkin(5);
t101 = rSges(3,3) + pkin(7);
t78 = sin(qJ(2));
t100 = t78 * pkin(2) + pkin(6);
t76 = qJ(2) + qJ(3);
t73 = sin(t76);
t79 = sin(qJ(1));
t99 = t73 * t79;
t82 = cos(qJ(1));
t98 = t73 * t82;
t74 = cos(t76);
t97 = t74 * t82;
t77 = sin(qJ(4));
t96 = t77 * t79;
t80 = cos(qJ(4));
t95 = t79 * t80;
t94 = t80 * t82;
t93 = t82 * t77;
t81 = cos(qJ(2));
t71 = pkin(2) * t81 + pkin(1);
t83 = -pkin(8) - pkin(7);
t92 = t79 * t71 + t82 * t83;
t91 = rSges(6,3) + qJ(5);
t90 = t73 * pkin(3) + t100;
t89 = t79 * t74 * pkin(3) + pkin(9) * t99 + t92;
t88 = t90 + (pkin(4) * t80 + qJ(5) * t77) * t73;
t57 = t74 * t95 - t93;
t87 = t57 * pkin(4) + t89;
t62 = t82 * t71;
t86 = pkin(3) * t97 + pkin(9) * t98 - t79 * t83 + t62;
t85 = rSges(3,1) * t81 - rSges(3,2) * t78 + pkin(1);
t59 = t74 * t94 + t96;
t84 = t59 * pkin(4) + t86;
t58 = t74 * t93 - t95;
t56 = t74 * t96 + t94;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t82 - rSges(2,2) * t79) + g(2) * (rSges(2,1) * t79 + rSges(2,2) * t82) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t78 + rSges(3,2) * t81 + pkin(6)) + (g(1) * t85 - g(2) * t101) * t82 + (g(1) * t101 + g(2) * t85) * t79) - m(4) * (g(1) * (rSges(4,1) * t97 - rSges(4,2) * t98 + t62) + g(2) * (-rSges(4,3) * t82 + t92) + g(3) * (rSges(4,1) * t73 + rSges(4,2) * t74 + t100) + (g(1) * (rSges(4,3) - t83) + g(2) * (rSges(4,1) * t74 - rSges(4,2) * t73)) * t79) - m(5) * (g(1) * (t59 * rSges(5,1) - t58 * rSges(5,2) + rSges(5,3) * t98 + t86) + g(2) * (rSges(5,1) * t57 - rSges(5,2) * t56 + rSges(5,3) * t99 + t89) + g(3) * ((-rSges(5,3) - pkin(9)) * t74 + (rSges(5,1) * t80 - rSges(5,2) * t77) * t73 + t90)) - m(6) * (g(1) * (rSges(6,1) * t98 - t59 * rSges(6,2) + t91 * t58 + t84) + g(2) * (rSges(6,1) * t99 - rSges(6,2) * t57 + t91 * t56 + t87) + g(3) * ((-rSges(6,1) - pkin(9)) * t74 + (-rSges(6,2) * t80 + rSges(6,3) * t77) * t73 + t88)) - m(7) * (g(1) * (t103 * t59 + t104 * t58 + t84) + g(2) * (t103 * t57 + t104 * t56 + t87) + (g(1) * t82 + g(2) * t79) * t73 * t102 + (t88 + (-pkin(9) - t102) * t74 + (rSges(7,2) * t77 + t103 * t80) * t73) * g(3));
U  = t1;
