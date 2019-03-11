% Calculate potential energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:02
% EndTime: 2019-03-10 01:18:03
% DurationCPUTime: 0.54s
% Computational Cost: add. (235->110), mult. (240->134), div. (0->0), fcn. (228->10), ass. (0->37)
t100 = rSges(4,3) + pkin(8);
t82 = -pkin(9) - pkin(8);
t99 = rSges(5,3) - t82;
t74 = -pkin(10) + t82;
t98 = rSges(6,3) - t74;
t97 = rSges(7,3) + qJ(6) - t74;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t96 = g(1) * t81 + g(2) * t78;
t79 = cos(qJ(3));
t64 = t79 * pkin(3) + pkin(2);
t77 = sin(qJ(2));
t92 = rSges(3,2) * t77;
t76 = sin(qJ(3));
t91 = t76 * t78;
t90 = t76 * t81;
t80 = cos(qJ(2));
t89 = t78 * t80;
t88 = t80 * t81;
t75 = qJ(3) + qJ(4);
t65 = sin(t75);
t59 = t76 * pkin(3) + pkin(4) * t65;
t84 = t81 * pkin(1) + t78 * pkin(7);
t66 = cos(t75);
t58 = pkin(4) * t66 + t64;
t71 = t78 * pkin(1);
t83 = -t81 * pkin(7) + t71;
t68 = qJ(5) + t75;
t63 = cos(t68);
t62 = sin(t68);
t57 = pkin(5) * t62 + t59;
t56 = pkin(5) * t63 + t58;
t55 = t62 * t78 + t63 * t88;
t54 = -t62 * t88 + t63 * t78;
t53 = -t62 * t81 + t63 * t89;
t52 = -t62 * t89 - t63 * t81;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t81 - rSges(2,2) * t78) + g(2) * (rSges(2,1) * t78 + rSges(2,2) * t81) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t78 + t84) + g(2) * (rSges(3,1) * t89 - t78 * t92 + t71) + g(3) * (rSges(3,1) * t77 + rSges(3,2) * t80 + pkin(6)) + (g(1) * (rSges(3,1) * t80 - t92) + g(2) * (-rSges(3,3) - pkin(7))) * t81) - m(4) * (g(1) * (pkin(2) * t88 + (t79 * t88 + t91) * rSges(4,1) + (-t76 * t88 + t78 * t79) * rSges(4,2) + t84) + g(2) * (pkin(2) * t89 + (t79 * t89 - t90) * rSges(4,1) + (-t76 * t89 - t79 * t81) * rSges(4,2) + t83) + g(3) * (-t100 * t80 + pkin(6)) + (g(3) * (rSges(4,1) * t79 - rSges(4,2) * t76 + pkin(2)) + t96 * t100) * t77) - m(5) * (g(1) * (t64 * t88 + pkin(3) * t91 + (t65 * t78 + t66 * t88) * rSges(5,1) + (-t65 * t88 + t66 * t78) * rSges(5,2) + t84) + g(2) * (t64 * t89 - pkin(3) * t90 + (-t65 * t81 + t66 * t89) * rSges(5,1) + (-t65 * t89 - t66 * t81) * rSges(5,2) + t83) + g(3) * (-t99 * t80 + pkin(6)) + (g(3) * (rSges(5,1) * t66 - rSges(5,2) * t65 + t64) + t96 * t99) * t77) - m(6) * (g(1) * (rSges(6,1) * t55 + rSges(6,2) * t54 + t58 * t88 + t59 * t78 + t84) + g(2) * (rSges(6,1) * t53 + rSges(6,2) * t52 + t58 * t89 - t59 * t81 + t83) + g(3) * (-t98 * t80 + pkin(6)) + (g(3) * (rSges(6,1) * t63 - rSges(6,2) * t62 + t58) + t96 * t98) * t77) - m(7) * (g(1) * (rSges(7,1) * t55 + rSges(7,2) * t54 + t56 * t88 + t57 * t78 + t84) + g(2) * (rSges(7,1) * t53 + rSges(7,2) * t52 + t56 * t89 - t57 * t81 + t83) + g(3) * (-t97 * t80 + pkin(6)) + (g(3) * (rSges(7,1) * t63 - rSges(7,2) * t62 + t56) + t96 * t97) * t77);
U  = t1;
