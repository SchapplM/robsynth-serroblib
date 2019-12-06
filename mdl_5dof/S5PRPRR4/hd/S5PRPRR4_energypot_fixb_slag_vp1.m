% Calculate potential energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:31
% EndTime: 2019-12-05 15:49:32
% DurationCPUTime: 0.46s
% Computational Cost: add. (253->100), mult. (549->143), div. (0->0), fcn. (668->12), ass. (0->48)
t121 = rSges(5,3) + pkin(7);
t120 = pkin(8) + rSges(6,3);
t94 = sin(pkin(5));
t96 = cos(pkin(9));
t119 = t96 * t94;
t100 = sin(qJ(2));
t97 = cos(pkin(5));
t80 = t97 * t100 * pkin(2) + (-pkin(6) - qJ(3)) * t94;
t103 = cos(qJ(2));
t89 = t103 * pkin(2) + pkin(1);
t93 = sin(pkin(9));
t118 = t96 * t80 + t93 * t89;
t117 = t100 * t94;
t102 = cos(qJ(4));
t116 = t102 * t94;
t95 = cos(pkin(10));
t115 = t103 * t95;
t114 = t93 * t100;
t113 = t93 * t103;
t112 = t96 * t100;
t111 = t96 * t103;
t110 = t97 * pkin(6) + qJ(1);
t92 = sin(pkin(10));
t105 = t100 * t95 + t103 * t92;
t79 = t105 * t97;
t82 = -t100 * t92 + t115;
t69 = t96 * t79 + t93 * t82;
t109 = t69 * pkin(3) + t118;
t108 = pkin(2) * t117 + t97 * qJ(3) + t110;
t71 = -t93 * t79 + t96 * t82;
t84 = t96 * t89;
t107 = t71 * pkin(3) - t93 * t80 + t84;
t78 = t105 * t94;
t106 = t78 * pkin(3) + t108;
t104 = t82 * t97;
t101 = cos(qJ(5));
t99 = sin(qJ(4));
t98 = sin(qJ(5));
t77 = -t115 * t94 + t117 * t92;
t73 = t78 * t102 + t97 * t99;
t72 = -t97 * t102 + t78 * t99;
t70 = -t104 * t93 - t105 * t96;
t68 = t104 * t96 - t105 * t93;
t65 = t93 * t94 * t99 + t71 * t102;
t64 = -t116 * t93 + t71 * t99;
t63 = t69 * t102 - t119 * t99;
t62 = t116 * t96 + t69 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t96 * rSges(2,1) - t93 * rSges(2,2)) + g(2) * (t93 * rSges(2,1) + t96 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t96 * pkin(1) + (-t114 * t97 + t111) * rSges(3,1) + (-t113 * t97 - t112) * rSges(3,2)) + g(2) * (t93 * pkin(1) + (t112 * t97 + t113) * rSges(3,1) + (t111 * t97 - t114) * rSges(3,2)) + g(3) * (t97 * rSges(3,3) + t110) + (g(3) * (rSges(3,1) * t100 + rSges(3,2) * t103) + (g(1) * t93 - g(2) * t96) * (rSges(3,3) + pkin(6))) * t94) - m(4) * (g(1) * (t71 * rSges(4,1) + t70 * rSges(4,2) + t84 + (rSges(4,3) * t94 - t80) * t93) + g(2) * (t69 * rSges(4,1) + t68 * rSges(4,2) - rSges(4,3) * t119 + t118) + g(3) * (t78 * rSges(4,1) - t77 * rSges(4,2) + t97 * rSges(4,3) + t108)) - m(5) * (g(1) * (t65 * rSges(5,1) - t64 * rSges(5,2) - t121 * t70 + t107) + g(2) * (t63 * rSges(5,1) - t62 * rSges(5,2) - t121 * t68 + t109) + g(3) * (t73 * rSges(5,1) - t72 * rSges(5,2) + t121 * t77 + t106)) - m(6) * (g(1) * (t65 * pkin(4) - t70 * pkin(7) + (t65 * t101 - t70 * t98) * rSges(6,1) + (-t70 * t101 - t65 * t98) * rSges(6,2) + t120 * t64 + t107) + g(2) * (t63 * pkin(4) - t68 * pkin(7) + (t63 * t101 - t68 * t98) * rSges(6,1) + (-t68 * t101 - t63 * t98) * rSges(6,2) + t120 * t62 + t109) + g(3) * (t73 * pkin(4) + t77 * pkin(7) + (t73 * t101 + t77 * t98) * rSges(6,1) + (t77 * t101 - t73 * t98) * rSges(6,2) + t120 * t72 + t106));
U = t1;
