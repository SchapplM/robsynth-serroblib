% Calculate potential energy for
% S6RPRRPP2
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:23
% EndTime: 2019-03-09 04:31:24
% DurationCPUTime: 0.39s
% Computational Cost: add. (232->90), mult. (234->108), div. (0->0), fcn. (228->8), ass. (0->34)
t95 = rSges(7,1) + pkin(5);
t94 = rSges(7,2) + qJ(5);
t71 = qJ(1) + pkin(9);
t66 = sin(t71);
t73 = sin(qJ(3));
t93 = t66 * t73;
t76 = cos(qJ(3));
t92 = t66 * t76;
t67 = cos(t71);
t91 = t67 * t73;
t72 = sin(qJ(4));
t90 = t72 * t76;
t75 = cos(qJ(4));
t89 = t75 * t76;
t88 = pkin(6) + qJ(2);
t74 = sin(qJ(1));
t69 = t74 * pkin(1);
t87 = t66 * pkin(2) + t69;
t86 = rSges(6,3) + qJ(5);
t85 = -rSges(7,3) - qJ(6);
t84 = t73 * pkin(3) + t88;
t77 = cos(qJ(1));
t70 = t77 * pkin(1);
t83 = t67 * pkin(2) + t66 * pkin(7) + t70;
t82 = t84 + (pkin(4) * t75 + qJ(5) * t72) * t73;
t81 = t67 * t76 * pkin(3) + pkin(8) * t91 + t83;
t55 = t66 * t72 + t67 * t89;
t80 = t55 * pkin(4) + t81;
t79 = pkin(3) * t92 - pkin(7) * t67 + pkin(8) * t93 + t87;
t53 = t66 * t89 - t67 * t72;
t78 = t53 * pkin(4) + t79;
t54 = -t66 * t75 + t67 * t90;
t52 = t66 * t90 + t67 * t75;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t77 - t74 * rSges(2,2)) + g(2) * (t74 * rSges(2,1) + rSges(2,2) * t77) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t67 - rSges(3,2) * t66 + t70) + g(2) * (rSges(3,1) * t66 + rSges(3,2) * t67 + t69) + g(3) * (rSges(3,3) + t88)) - m(4) * (g(1) * (rSges(4,3) * t66 + t83) + g(2) * (rSges(4,1) * t92 - rSges(4,2) * t93 + t87) + g(3) * (rSges(4,1) * t73 + rSges(4,2) * t76 + t88) + (g(1) * (rSges(4,1) * t76 - rSges(4,2) * t73) + g(2) * (-rSges(4,3) - pkin(7))) * t67) - m(5) * (g(1) * (rSges(5,1) * t55 - rSges(5,2) * t54 + rSges(5,3) * t91 + t81) + g(2) * (rSges(5,1) * t53 - rSges(5,2) * t52 + rSges(5,3) * t93 + t79) + g(3) * ((-rSges(5,3) - pkin(8)) * t76 + (rSges(5,1) * t75 - rSges(5,2) * t72) * t73 + t84)) - m(6) * (g(1) * (rSges(6,1) * t55 + rSges(6,2) * t91 + t86 * t54 + t80) + g(2) * (rSges(6,1) * t53 + rSges(6,2) * t93 + t86 * t52 + t78) + g(3) * ((-rSges(6,2) - pkin(8)) * t76 + (rSges(6,1) * t75 + rSges(6,3) * t72) * t73 + t82)) - m(7) * (g(1) * (t94 * t54 + t95 * t55 + t80) + g(2) * (t94 * t52 + t95 * t53 + t78) + (g(1) * t67 + g(2) * t66) * t73 * t85 + (t82 + (-pkin(8) - t85) * t76 + (rSges(7,2) * t72 + t95 * t75) * t73) * g(3));
U  = t1;
