% Calculate potential energy for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:32
% EndTime: 2019-03-09 04:42:33
% DurationCPUTime: 0.40s
% Computational Cost: add. (223->92), mult. (248->110), div. (0->0), fcn. (242->8), ass. (0->38)
t102 = rSges(7,1) + pkin(5);
t101 = rSges(7,2) + qJ(5);
t75 = sin(pkin(9));
t100 = t75 * pkin(2) + pkin(6);
t74 = pkin(9) + qJ(3);
t71 = sin(t74);
t79 = sin(qJ(1));
t99 = t71 * t79;
t81 = cos(qJ(1));
t98 = t71 * t81;
t72 = cos(t74);
t97 = t72 * t81;
t78 = sin(qJ(4));
t96 = t78 * t81;
t95 = t79 * t78;
t80 = cos(qJ(4));
t94 = t79 * t80;
t93 = t80 * t81;
t76 = cos(pkin(9));
t68 = pkin(2) * t76 + pkin(1);
t77 = -pkin(7) - qJ(2);
t92 = t79 * t68 + t81 * t77;
t91 = rSges(3,3) + qJ(2);
t90 = rSges(6,3) + qJ(5);
t89 = -rSges(7,3) - qJ(6);
t88 = t71 * pkin(3) + t100;
t87 = t79 * t72 * pkin(3) + pkin(8) * t99 + t92;
t86 = t88 + (pkin(4) * t80 + qJ(5) * t78) * t71;
t56 = t72 * t94 - t96;
t85 = t56 * pkin(4) + t87;
t61 = t81 * t68;
t84 = pkin(3) * t97 + pkin(8) * t98 - t79 * t77 + t61;
t83 = rSges(3,1) * t76 - rSges(3,2) * t75 + pkin(1);
t58 = t72 * t93 + t95;
t82 = t58 * pkin(4) + t84;
t57 = t72 * t96 - t94;
t55 = t72 * t95 + t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t81 - t79 * rSges(2,2)) + g(2) * (t79 * rSges(2,1) + rSges(2,2) * t81) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t75 + rSges(3,2) * t76 + pkin(6)) + (g(1) * t83 - g(2) * t91) * t81 + (g(1) * t91 + g(2) * t83) * t79) - m(4) * (g(1) * (rSges(4,1) * t97 - rSges(4,2) * t98 + t61) + g(2) * (-rSges(4,3) * t81 + t92) + g(3) * (rSges(4,1) * t71 + rSges(4,2) * t72 + t100) + (g(1) * (rSges(4,3) - t77) + g(2) * (rSges(4,1) * t72 - rSges(4,2) * t71)) * t79) - m(5) * (g(1) * (t58 * rSges(5,1) - t57 * rSges(5,2) + rSges(5,3) * t98 + t84) + g(2) * (rSges(5,1) * t56 - rSges(5,2) * t55 + rSges(5,3) * t99 + t87) + g(3) * ((-rSges(5,3) - pkin(8)) * t72 + (rSges(5,1) * t80 - rSges(5,2) * t78) * t71 + t88)) - m(6) * (g(1) * (t58 * rSges(6,1) + rSges(6,2) * t98 + t90 * t57 + t82) + g(2) * (rSges(6,1) * t56 + rSges(6,2) * t99 + t90 * t55 + t85) + g(3) * ((-rSges(6,2) - pkin(8)) * t72 + (rSges(6,1) * t80 + rSges(6,3) * t78) * t71 + t86)) - m(7) * (g(1) * (t101 * t57 + t102 * t58 + t82) + g(2) * (t101 * t55 + t102 * t56 + t85) + (g(1) * t81 + g(2) * t79) * t71 * t89 + (t86 + (-pkin(8) - t89) * t72 + (rSges(7,2) * t78 + t102 * t80) * t71) * g(3));
U  = t1;
