% Calculate potential energy for
% S6RRRRPP2
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:48:52
% EndTime: 2019-03-09 20:48:52
% DurationCPUTime: 0.41s
% Computational Cost: add. (223->92), mult. (248->110), div. (0->0), fcn. (242->8), ass. (0->38)
t102 = rSges(7,1) + pkin(5);
t101 = rSges(7,2) + qJ(5);
t100 = rSges(3,3) + pkin(7);
t76 = sin(qJ(2));
t99 = t76 * pkin(2) + pkin(6);
t74 = qJ(2) + qJ(3);
t71 = sin(t74);
t77 = sin(qJ(1));
t98 = t71 * t77;
t80 = cos(qJ(1));
t97 = t71 * t80;
t72 = cos(t74);
t96 = t72 * t80;
t75 = sin(qJ(4));
t95 = t75 * t77;
t94 = t75 * t80;
t78 = cos(qJ(4));
t93 = t77 * t78;
t92 = t78 * t80;
t79 = cos(qJ(2));
t69 = pkin(2) * t79 + pkin(1);
t81 = -pkin(8) - pkin(7);
t91 = t77 * t69 + t80 * t81;
t90 = rSges(6,3) + qJ(5);
t89 = -rSges(7,3) - qJ(6);
t88 = t71 * pkin(3) + t99;
t87 = t77 * t72 * pkin(3) + pkin(9) * t98 + t91;
t86 = t88 + (pkin(4) * t78 + qJ(5) * t75) * t71;
t56 = t72 * t93 - t94;
t85 = t56 * pkin(4) + t87;
t61 = t80 * t69;
t84 = pkin(3) * t96 + pkin(9) * t97 - t77 * t81 + t61;
t83 = rSges(3,1) * t79 - rSges(3,2) * t76 + pkin(1);
t58 = t72 * t92 + t95;
t82 = t58 * pkin(4) + t84;
t57 = t72 * t94 - t93;
t55 = t72 * t95 + t92;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t80 * rSges(2,1) - t77 * rSges(2,2)) + g(2) * (rSges(2,1) * t77 + rSges(2,2) * t80) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t76 + rSges(3,2) * t79 + pkin(6)) + (g(1) * t83 - g(2) * t100) * t80 + (g(1) * t100 + g(2) * t83) * t77) - m(4) * (g(1) * (rSges(4,1) * t96 - rSges(4,2) * t97 + t61) + g(2) * (-rSges(4,3) * t80 + t91) + g(3) * (rSges(4,1) * t71 + rSges(4,2) * t72 + t99) + (g(1) * (rSges(4,3) - t81) + g(2) * (rSges(4,1) * t72 - rSges(4,2) * t71)) * t77) - m(5) * (g(1) * (t58 * rSges(5,1) - t57 * rSges(5,2) + rSges(5,3) * t97 + t84) + g(2) * (rSges(5,1) * t56 - rSges(5,2) * t55 + rSges(5,3) * t98 + t87) + g(3) * ((-rSges(5,3) - pkin(9)) * t72 + (rSges(5,1) * t78 - rSges(5,2) * t75) * t71 + t88)) - m(6) * (g(1) * (t58 * rSges(6,1) + rSges(6,2) * t97 + t90 * t57 + t82) + g(2) * (rSges(6,1) * t56 + rSges(6,2) * t98 + t90 * t55 + t85) + g(3) * ((-rSges(6,2) - pkin(9)) * t72 + (rSges(6,1) * t78 + rSges(6,3) * t75) * t71 + t86)) - m(7) * (g(1) * (t101 * t57 + t102 * t58 + t82) + g(2) * (t101 * t55 + t102 * t56 + t85) + (g(1) * t80 + g(2) * t77) * t71 * t89 + (t86 + (-pkin(9) - t89) * t72 + (rSges(7,2) * t75 + t102 * t78) * t71) * g(3));
U  = t1;
