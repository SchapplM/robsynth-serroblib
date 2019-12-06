% Calculate potential energy for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:12
% EndTime: 2019-12-05 16:25:13
% DurationCPUTime: 0.37s
% Computational Cost: add. (227->103), mult. (371->142), div. (0->0), fcn. (422->12), ass. (0->44)
t95 = sin(qJ(3));
t116 = pkin(3) * t95;
t115 = pkin(7) + rSges(4,3);
t114 = pkin(8) + rSges(6,3);
t89 = sin(pkin(9));
t90 = sin(pkin(5));
t113 = t89 * t90;
t91 = cos(pkin(9));
t112 = t90 * t91;
t111 = t90 * t95;
t96 = sin(qJ(2));
t110 = t90 * t96;
t98 = cos(qJ(3));
t109 = t90 * t98;
t99 = cos(qJ(2));
t108 = t90 * t99;
t92 = cos(pkin(5));
t107 = t92 * t96;
t106 = t92 * t99;
t105 = t91 * pkin(1) + pkin(6) * t113;
t104 = t92 * pkin(6) + qJ(1);
t103 = t89 * t111;
t72 = t106 * t89 + t91 * t96;
t73 = -t107 * t89 + t91 * t99;
t82 = pkin(3) * t98 + pkin(2);
t93 = -qJ(4) - pkin(7);
t102 = pkin(3) * t103 - t72 * t93 + t73 * t82 + t105;
t101 = t93 * t108 + t82 * t110 + t92 * t116 + t104;
t70 = -t106 * t91 + t89 * t96;
t71 = t107 * t91 + t89 * t99;
t85 = t89 * pkin(1);
t100 = t71 * t82 + t85 + (-pkin(6) - t116) * t112 - t70 * t93;
t97 = cos(qJ(5));
t94 = sin(qJ(5));
t88 = qJ(3) + pkin(10);
t84 = cos(t88);
t83 = sin(t88);
t67 = t110 * t84 + t83 * t92;
t66 = t110 * t83 - t92 * t84;
t63 = t113 * t83 + t73 * t84;
t62 = -t113 * t84 + t73 * t83;
t61 = -t112 * t83 + t71 * t84;
t60 = t112 * t84 + t71 * t83;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t91 - rSges(2,2) * t89) + g(2) * (rSges(2,1) * t89 + rSges(2,2) * t91) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t73 - rSges(3,2) * t72 + t105) + g(2) * (rSges(3,1) * t71 - rSges(3,2) * t70 + t85) + g(3) * (t92 * rSges(3,3) + t104) + (g(1) * rSges(3,3) * t89 + g(3) * (rSges(3,1) * t96 + rSges(3,2) * t99) + g(2) * (-rSges(3,3) - pkin(6)) * t91) * t90) - m(4) * (g(1) * (t73 * pkin(2) + (t73 * t98 + t103) * rSges(4,1) + (t109 * t89 - t73 * t95) * rSges(4,2) + t115 * t72 + t105) + g(2) * (t71 * pkin(2) + t85 - pkin(6) * t112 + (-t111 * t91 + t71 * t98) * rSges(4,1) + (-t109 * t91 - t71 * t95) * rSges(4,2) + t115 * t70) + g(3) * ((rSges(4,1) * t95 + rSges(4,2) * t98) * t92 + (-t115 * t99 + (rSges(4,1) * t98 - rSges(4,2) * t95 + pkin(2)) * t96) * t90 + t104)) - m(5) * (g(1) * (rSges(5,1) * t63 - rSges(5,2) * t62 + rSges(5,3) * t72 + t102) + g(2) * (rSges(5,1) * t61 - rSges(5,2) * t60 + rSges(5,3) * t70 + t100) + g(3) * (rSges(5,1) * t67 - rSges(5,2) * t66 - rSges(5,3) * t108 + t101)) - m(6) * (g(1) * (t63 * pkin(4) + (t63 * t97 + t72 * t94) * rSges(6,1) + (-t63 * t94 + t72 * t97) * rSges(6,2) + t114 * t62 + t102) + g(2) * (t61 * pkin(4) + (t61 * t97 + t70 * t94) * rSges(6,1) + (-t61 * t94 + t70 * t97) * rSges(6,2) + t114 * t60 + t100) + g(3) * (t67 * pkin(4) + (-t108 * t94 + t67 * t97) * rSges(6,1) + (-t108 * t97 - t67 * t94) * rSges(6,2) + t114 * t66 + t101));
U = t1;
