% Calculate potential energy for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:20
% EndTime: 2019-03-09 05:04:21
% DurationCPUTime: 0.47s
% Computational Cost: add. (250->99), mult. (260->125), div. (0->0), fcn. (264->10), ass. (0->34)
t95 = -rSges(7,3) - pkin(9);
t71 = qJ(1) + pkin(10);
t66 = sin(t71);
t74 = sin(qJ(3));
t94 = t66 * t74;
t78 = cos(qJ(3));
t93 = t66 * t78;
t67 = cos(t71);
t92 = t67 * t74;
t73 = sin(qJ(4));
t91 = t73 * t78;
t77 = cos(qJ(4));
t90 = t77 * t78;
t89 = pkin(6) + qJ(2);
t75 = sin(qJ(1));
t69 = t75 * pkin(1);
t88 = t66 * pkin(2) + t69;
t87 = rSges(6,3) + qJ(5);
t86 = t74 * pkin(3) + t89;
t79 = cos(qJ(1));
t70 = t79 * pkin(1);
t85 = t67 * pkin(2) + t66 * pkin(7) + t70;
t84 = t86 + (pkin(4) * t77 + qJ(5) * t73) * t74;
t83 = t67 * t78 * pkin(3) + pkin(8) * t92 + t85;
t55 = t66 * t73 + t67 * t90;
t82 = t55 * pkin(4) + t83;
t81 = pkin(3) * t93 - t67 * pkin(7) + pkin(8) * t94 + t88;
t53 = t66 * t90 - t67 * t73;
t80 = t53 * pkin(4) + t81;
t76 = cos(qJ(6));
t72 = sin(qJ(6));
t54 = -t66 * t77 + t67 * t91;
t52 = t66 * t91 + t67 * t77;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - t75 * rSges(2,2)) + g(2) * (t75 * rSges(2,1) + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t67 - rSges(3,2) * t66 + t70) + g(2) * (rSges(3,1) * t66 + rSges(3,2) * t67 + t69) + g(3) * (rSges(3,3) + t89)) - m(4) * (g(1) * (rSges(4,3) * t66 + t85) + g(2) * (rSges(4,1) * t93 - rSges(4,2) * t94 + t88) + g(3) * (rSges(4,1) * t74 + rSges(4,2) * t78 + t89) + (g(1) * (rSges(4,1) * t78 - rSges(4,2) * t74) + g(2) * (-rSges(4,3) - pkin(7))) * t67) - m(5) * (g(1) * (rSges(5,1) * t55 - rSges(5,2) * t54 + rSges(5,3) * t92 + t83) + g(2) * (rSges(5,1) * t53 - rSges(5,2) * t52 + rSges(5,3) * t94 + t81) + g(3) * ((-rSges(5,3) - pkin(8)) * t78 + (rSges(5,1) * t77 - rSges(5,2) * t73) * t74 + t86)) - m(6) * (g(1) * (rSges(6,1) * t55 + rSges(6,2) * t92 + t87 * t54 + t82) + g(2) * (rSges(6,1) * t53 + rSges(6,2) * t94 + t87 * t52 + t80) + g(3) * ((-rSges(6,2) - pkin(8)) * t78 + (rSges(6,1) * t77 + rSges(6,3) * t73) * t74 + t84)) - m(7) * (g(1) * (t55 * pkin(5) + t54 * qJ(5) + (t54 * t72 + t55 * t76) * rSges(7,1) + (t54 * t76 - t55 * t72) * rSges(7,2) + t82) + g(2) * (t53 * pkin(5) + t52 * qJ(5) + (t52 * t72 + t53 * t76) * rSges(7,1) + (t52 * t76 - t53 * t72) * rSges(7,2) + t80) + (g(1) * t67 + g(2) * t66) * t74 * t95 + (t84 + (-pkin(8) - t95) * t78 + (t77 * pkin(5) + (t72 * t73 + t76 * t77) * rSges(7,1) + (-t72 * t77 + t73 * t76) * rSges(7,2)) * t74) * g(3));
U  = t1;
