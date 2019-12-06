% Calculate potential energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:01
% EndTime: 2019-12-05 16:54:02
% DurationCPUTime: 0.29s
% Computational Cost: add. (207->96), mult. (435->129), div. (0->0), fcn. (511->10), ass. (0->45)
t85 = sin(pkin(5));
t110 = pkin(6) * t85;
t109 = rSges(4,3) + pkin(7);
t108 = rSges(5,3) + pkin(8);
t90 = sin(qJ(3));
t107 = t85 * t90;
t91 = sin(qJ(2));
t106 = t85 * t91;
t93 = cos(qJ(3));
t105 = t85 * t93;
t94 = cos(qJ(2));
t104 = t85 * t94;
t87 = cos(pkin(5));
t103 = t87 * t91;
t102 = t87 * t94;
t101 = rSges(6,3) + qJ(5) + pkin(8);
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t100 = t86 * pkin(1) + t84 * t110;
t99 = t87 * pkin(6) + qJ(1);
t72 = -t84 * t103 + t86 * t94;
t98 = t72 * pkin(2) + t100;
t97 = pkin(2) * t106 + t99;
t89 = sin(qJ(4));
t96 = pkin(4) * t89 + pkin(7);
t70 = t86 * t103 + t84 * t94;
t81 = t84 * pkin(1);
t95 = t70 * pkin(2) - t86 * t110 + t81;
t92 = cos(qJ(4));
t80 = t92 * pkin(4) + pkin(3);
t74 = t91 * t105 + t87 * t90;
t73 = t90 * t106 - t87 * t93;
t71 = t84 * t102 + t86 * t91;
t69 = -t86 * t102 + t84 * t91;
t66 = -t89 * t104 + t74 * t92;
t65 = -t92 * t104 - t74 * t89;
t64 = t84 * t107 + t72 * t93;
t63 = -t84 * t105 + t72 * t90;
t62 = -t86 * t107 + t70 * t93;
t61 = t86 * t105 + t70 * t90;
t60 = t64 * t92 + t71 * t89;
t59 = -t64 * t89 + t71 * t92;
t58 = t62 * t92 + t69 * t89;
t57 = -t62 * t89 + t69 * t92;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t86 * rSges(2,1) - t84 * rSges(2,2)) + g(2) * (t84 * rSges(2,1) + t86 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t72 * rSges(3,1) - t71 * rSges(3,2) + t100) + g(2) * (t70 * rSges(3,1) - t69 * rSges(3,2) + t81) + g(3) * (t87 * rSges(3,3) + t99) + (g(1) * rSges(3,3) * t84 + g(3) * (rSges(3,1) * t91 + rSges(3,2) * t94) + g(2) * (-rSges(3,3) - pkin(6)) * t86) * t85) - m(4) * (g(1) * (t64 * rSges(4,1) - t63 * rSges(4,2) + t109 * t71 + t98) + g(2) * (t62 * rSges(4,1) - t61 * rSges(4,2) + t109 * t69 + t95) + g(3) * (t74 * rSges(4,1) - t73 * rSges(4,2) - t109 * t104 + t97)) - m(5) * (g(1) * (t60 * rSges(5,1) + t59 * rSges(5,2) + t64 * pkin(3) + t71 * pkin(7) + t108 * t63 + t98) + g(2) * (t58 * rSges(5,1) + t57 * rSges(5,2) + t62 * pkin(3) + t69 * pkin(7) + t108 * t61 + t95) + g(3) * (t66 * rSges(5,1) + t65 * rSges(5,2) + t74 * pkin(3) - pkin(7) * t104 + t108 * t73 + t97)) - m(6) * (g(1) * (t60 * rSges(6,1) + t59 * rSges(6,2) + t101 * t63 + t64 * t80 + t96 * t71 + t98) + g(2) * (t58 * rSges(6,1) + t57 * rSges(6,2) + t101 * t61 + t62 * t80 + t96 * t69 + t95) + g(3) * (t66 * rSges(6,1) + t65 * rSges(6,2) + t101 * t73 - t96 * t104 + t74 * t80 + t97));
U = t1;
