% Calculate potential energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10V2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:31
% EndTime: 2019-04-11 14:41:32
% DurationCPUTime: 0.25s
% Computational Cost: add. (223->90), mult. (277->123), div. (0->0), fcn. (295->12), ass. (0->44)
t78 = qJ(2) + qJ(3);
t76 = cos(t78);
t105 = pkin(3) * t76;
t82 = sin(qJ(2));
t104 = t82 * pkin(2) + pkin(4);
t103 = pkin(6) + rSges(7,3);
t75 = sin(t78);
t81 = sin(qJ(4));
t102 = t75 * t81;
t83 = sin(qJ(1));
t101 = t75 * t83;
t86 = cos(qJ(4));
t100 = t75 * t86;
t88 = cos(qJ(1));
t99 = t75 * t88;
t98 = t83 * t81;
t97 = t83 * t86;
t96 = t88 * t81;
t95 = t88 * t86;
t94 = t75 * pkin(3) + t104;
t87 = cos(qJ(2));
t74 = t87 * pkin(2) + pkin(1);
t67 = t83 * t74;
t93 = pkin(5) * t101 + t83 * t105 + t67;
t68 = t88 * t74;
t92 = pkin(5) * t99 + t88 * t105 + t68;
t91 = rSges(4,1) * t76 - rSges(4,2) * t75;
t90 = -t76 * pkin(5) + t94;
t89 = rSges(3,1) * t87 - rSges(3,2) * t82 + pkin(1);
t85 = cos(qJ(5));
t84 = cos(qJ(6));
t80 = sin(qJ(5));
t79 = sin(qJ(6));
t64 = t76 * t95 + t98;
t63 = t76 * t96 - t97;
t62 = t76 * t97 - t96;
t61 = t76 * t98 + t95;
t60 = t100 * t85 - t76 * t80;
t59 = t100 * t80 + t76 * t85;
t58 = t64 * t85 + t80 * t99;
t57 = t64 * t80 - t85 * t99;
t56 = t101 * t80 + t62 * t85;
t55 = -t101 * t85 + t62 * t80;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t88 * rSges(2,1) - t83 * rSges(2,2)) + g(2) * (t83 * rSges(2,1) + t88 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(3) * (t82 * rSges(3,1) + t87 * rSges(3,2) + pkin(4)) + (-g(2) * rSges(3,3) + g(1) * t89) * t88 + (g(1) * rSges(3,3) + g(2) * t89) * t83) - m(4) * (g(1) * (t83 * rSges(4,3) + t88 * t91 + t68) + g(2) * (-t88 * rSges(4,3) + t83 * t91 + t67) + g(3) * (t75 * rSges(4,1) + t76 * rSges(4,2) + t104)) - m(5) * (g(1) * (t64 * rSges(5,1) - t63 * rSges(5,2) + rSges(5,3) * t99 + t92) + g(2) * (t62 * rSges(5,1) - t61 * rSges(5,2) + rSges(5,3) * t101 + t93) + g(3) * ((-rSges(5,3) - pkin(5)) * t76 + (rSges(5,1) * t86 - rSges(5,2) * t81) * t75 + t94)) - m(6) * (g(1) * (t58 * rSges(6,1) - t57 * rSges(6,2) + t63 * rSges(6,3) + t92) + g(2) * (t56 * rSges(6,1) - t55 * rSges(6,2) + t61 * rSges(6,3) + t93) + g(3) * (t60 * rSges(6,1) - t59 * rSges(6,2) + rSges(6,3) * t102 + t90)) - m(7) * (g(1) * ((t58 * t84 + t63 * t79) * rSges(7,1) + (-t58 * t79 + t63 * t84) * rSges(7,2) + t103 * t57 + t92) + g(2) * ((t56 * t84 + t61 * t79) * rSges(7,1) + (-t56 * t79 + t61 * t84) * rSges(7,2) + t103 * t55 + t93) + g(3) * ((t102 * t79 + t60 * t84) * rSges(7,1) + (t102 * t84 - t60 * t79) * rSges(7,2) + t103 * t59 + t90));
U  = t1;
