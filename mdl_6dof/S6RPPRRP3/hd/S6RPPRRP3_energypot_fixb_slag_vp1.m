% Calculate potential energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:38
% EndTime: 2019-03-09 02:02:39
% DurationCPUTime: 0.33s
% Computational Cost: add. (191->84), mult. (174->100), div. (0->0), fcn. (158->8), ass. (0->34)
t91 = rSges(7,1) + pkin(5);
t90 = rSges(7,3) + qJ(6);
t68 = sin(qJ(4));
t89 = -pkin(4) * t68 - qJ(3);
t66 = qJ(1) + pkin(9);
t60 = sin(t66);
t87 = g(1) * t60;
t61 = cos(t66);
t86 = g(2) * t61;
t71 = cos(qJ(4));
t85 = rSges(5,2) * t71;
t84 = t60 * t68;
t67 = sin(qJ(5));
t83 = t67 * t68;
t70 = cos(qJ(5));
t82 = t68 * t70;
t81 = pkin(6) + qJ(2);
t69 = sin(qJ(1));
t63 = t69 * pkin(1);
t80 = t60 * pkin(2) + t63;
t79 = pkin(3) + t81;
t72 = cos(qJ(1));
t65 = t72 * pkin(1);
t78 = t61 * pkin(2) + t60 * qJ(3) + t65;
t77 = t60 * pkin(7) + t80;
t76 = t61 * t71 * pkin(8) + t77;
t75 = t61 * pkin(7) + t78;
t74 = t71 * pkin(4) + t68 * pkin(8) + t79;
t73 = pkin(4) * t84 + t75;
t51 = t60 * t67 - t61 * t82;
t50 = t60 * t70 + t61 * t83;
t49 = t60 * t82 + t61 * t67;
t48 = t60 * t83 - t61 * t70;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t72 * rSges(2,1) - t69 * rSges(2,2)) + g(2) * (t69 * rSges(2,1) + t72 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t61 * rSges(3,1) - t60 * rSges(3,2) + t65) + g(2) * (t60 * rSges(3,1) + t61 * rSges(3,2) + t63) + g(3) * (rSges(3,3) + t81)) - m(4) * (g(1) * (-t61 * rSges(4,2) + t60 * rSges(4,3) + t78) + g(2) * (-t60 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t61 + t80) + g(3) * (rSges(4,1) + t81)) - m(5) * (g(1) * (rSges(5,1) * t84 + t60 * t85 + t75) + g(2) * (t60 * rSges(5,3) + t77) + g(3) * (t71 * rSges(5,1) - t68 * rSges(5,2) + t79) + (g(1) * rSges(5,3) + g(2) * (-rSges(5,1) * t68 - qJ(3) - t85)) * t61) - m(6) * (g(1) * (t49 * rSges(6,1) - t48 * rSges(6,2) + t73) + g(2) * (t51 * rSges(6,1) + t50 * rSges(6,2) + t76) + g(3) * (t68 * rSges(6,3) + t74) + (g(3) * (rSges(6,1) * t70 - rSges(6,2) * t67) + (-rSges(6,3) - pkin(8)) * t87) * t71 + (rSges(6,3) * t71 + t89) * t86) - m(7) * (g(1) * (t90 * t48 + t91 * t49 + t73) + g(2) * (-t90 * t50 + t91 * t51 + t89 * t61 + t76) + g(3) * (t68 * rSges(7,2) + t74) + (rSges(7,2) * t86 + g(3) * (t90 * t67 + t91 * t70) + (-rSges(7,2) - pkin(8)) * t87) * t71);
U  = t1;
