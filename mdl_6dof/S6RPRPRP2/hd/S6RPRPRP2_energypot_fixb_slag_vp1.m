% Calculate potential energy for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:26
% EndTime: 2019-03-09 03:04:26
% DurationCPUTime: 0.32s
% Computational Cost: add. (238->87), mult. (185->103), div. (0->0), fcn. (169->10), ass. (0->39)
t99 = rSges(7,1) + pkin(5);
t98 = rSges(7,3) + qJ(6);
t97 = rSges(4,3) + pkin(7);
t73 = qJ(3) + pkin(10);
t66 = sin(t73);
t74 = qJ(1) + pkin(9);
t67 = sin(t74);
t96 = t66 * t67;
t69 = cos(t74);
t95 = t66 * t69;
t76 = sin(qJ(5));
t94 = t67 * t76;
t79 = cos(qJ(5));
t93 = t67 * t79;
t68 = cos(t73);
t92 = t68 * t69;
t91 = t69 * t76;
t90 = t69 * t79;
t89 = pkin(6) + qJ(2);
t80 = cos(qJ(3));
t65 = t80 * pkin(3) + pkin(2);
t81 = cos(qJ(1));
t72 = t81 * pkin(1);
t88 = t69 * t65 + t72;
t77 = sin(qJ(3));
t87 = t77 * pkin(3) + t89;
t78 = sin(qJ(1));
t71 = t78 * pkin(1);
t75 = -qJ(4) - pkin(7);
t86 = t67 * t65 + t69 * t75 + t71;
t85 = t66 * pkin(4) + t87;
t84 = t67 * t68 * pkin(4) + pkin(8) * t96 + t86;
t83 = rSges(4,1) * t80 - rSges(4,2) * t77 + pkin(2);
t82 = pkin(4) * t92 + pkin(8) * t95 - t67 * t75 + t88;
t55 = t68 * t90 + t94;
t54 = t68 * t91 - t93;
t53 = t68 * t93 - t91;
t52 = t68 * t94 + t90;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t81 * rSges(2,1) - t78 * rSges(2,2)) + g(2) * (t78 * rSges(2,1) + t81 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t69 * rSges(3,1) - t67 * rSges(3,2) + t72) + g(2) * (t67 * rSges(3,1) + t69 * rSges(3,2) + t71) + g(3) * (rSges(3,3) + t89)) - m(4) * (g(1) * t72 + g(2) * t71 + g(3) * (t77 * rSges(4,1) + t80 * rSges(4,2) + t89) + (g(1) * t83 - g(2) * t97) * t69 + (g(1) * t97 + g(2) * t83) * t67) - m(5) * (g(1) * (rSges(5,1) * t92 - rSges(5,2) * t95 + t88) + g(2) * (-t69 * rSges(5,3) + t86) + g(3) * (t66 * rSges(5,1) + t68 * rSges(5,2) + t87) + (g(1) * (rSges(5,3) - t75) + g(2) * (rSges(5,1) * t68 - rSges(5,2) * t66)) * t67) - m(6) * (g(1) * (t55 * rSges(6,1) - t54 * rSges(6,2) + rSges(6,3) * t95 + t82) + g(2) * (t53 * rSges(6,1) - t52 * rSges(6,2) + rSges(6,3) * t96 + t84) + g(3) * ((-rSges(6,3) - pkin(8)) * t68 + (rSges(6,1) * t79 - rSges(6,2) * t76) * t66 + t85)) - m(7) * (g(1) * (t98 * t54 + t99 * t55 + t82) + g(2) * (t98 * t52 + t99 * t53 + t84) + g(3) * (t85 + (-rSges(7,2) - pkin(8)) * t68) + (g(3) * (t98 * t76 + t99 * t79) + (g(1) * t69 + g(2) * t67) * rSges(7,2)) * t66);
U  = t1;
