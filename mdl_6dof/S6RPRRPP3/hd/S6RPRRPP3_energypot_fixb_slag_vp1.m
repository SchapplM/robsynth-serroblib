% Calculate potential energy for
% S6RPRRPP3
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:34:57
% EndTime: 2019-03-09 04:34:58
% DurationCPUTime: 0.39s
% Computational Cost: add. (232->90), mult. (234->108), div. (0->0), fcn. (228->8), ass. (0->34)
t97 = rSges(7,2) + qJ(5);
t96 = rSges(7,3) + qJ(6);
t95 = rSges(7,1) + pkin(5);
t73 = qJ(1) + pkin(9);
t68 = sin(t73);
t75 = sin(qJ(3));
t94 = t68 * t75;
t78 = cos(qJ(3));
t93 = t68 * t78;
t69 = cos(t73);
t92 = t69 * t75;
t74 = sin(qJ(4));
t91 = t74 * t78;
t77 = cos(qJ(4));
t90 = t77 * t78;
t89 = pkin(6) + qJ(2);
t76 = sin(qJ(1));
t71 = t76 * pkin(1);
t88 = t68 * pkin(2) + t71;
t87 = rSges(6,3) + qJ(5);
t86 = t75 * pkin(3) + t89;
t79 = cos(qJ(1));
t72 = t79 * pkin(1);
t85 = t69 * pkin(2) + t68 * pkin(7) + t72;
t84 = t86 + (pkin(4) * t77 + qJ(5) * t74) * t75;
t83 = t69 * t78 * pkin(3) + pkin(8) * t92 + t85;
t56 = t68 * t74 + t69 * t90;
t82 = t56 * pkin(4) + t83;
t81 = pkin(3) * t93 - pkin(7) * t69 + pkin(8) * t94 + t88;
t54 = t68 * t90 - t69 * t74;
t80 = t54 * pkin(4) + t81;
t55 = -t68 * t77 + t69 * t91;
t53 = t68 * t91 + t69 * t77;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - t76 * rSges(2,2)) + g(2) * (t76 * rSges(2,1) + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t69 - rSges(3,2) * t68 + t72) + g(2) * (rSges(3,1) * t68 + rSges(3,2) * t69 + t71) + g(3) * (rSges(3,3) + t89)) - m(4) * (g(1) * (rSges(4,3) * t68 + t85) + g(2) * (rSges(4,1) * t93 - rSges(4,2) * t94 + t88) + g(3) * (rSges(4,1) * t75 + rSges(4,2) * t78 + t89) + (g(1) * (rSges(4,1) * t78 - rSges(4,2) * t75) + g(2) * (-rSges(4,3) - pkin(7))) * t69) - m(5) * (g(1) * (rSges(5,1) * t56 - rSges(5,2) * t55 + rSges(5,3) * t92 + t83) + g(2) * (rSges(5,1) * t54 - rSges(5,2) * t53 + rSges(5,3) * t94 + t81) + g(3) * ((-rSges(5,3) - pkin(8)) * t78 + (rSges(5,1) * t77 - rSges(5,2) * t74) * t75 + t86)) - m(6) * (g(1) * (rSges(6,1) * t92 - rSges(6,2) * t56 + t87 * t55 + t82) + g(2) * (rSges(6,1) * t94 - rSges(6,2) * t54 + t87 * t53 + t80) + g(3) * ((-rSges(6,1) - pkin(8)) * t78 + (-rSges(6,2) * t77 + rSges(6,3) * t74) * t75 + t84)) - m(7) * (g(1) * (t55 * t97 + t96 * t56 + t82) + g(2) * (t53 * t97 + t96 * t54 + t80) + (g(1) * t69 + g(2) * t68) * t75 * t95 + (t84 + (-pkin(8) - t95) * t78 + (rSges(7,2) * t74 + t77 * t96) * t75) * g(3));
U  = t1;
