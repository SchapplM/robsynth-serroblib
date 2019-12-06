% Calculate potential energy for
% S5PRPRR6
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:21
% EndTime: 2019-12-05 15:56:22
% DurationCPUTime: 0.36s
% Computational Cost: add. (227->103), mult. (371->140), div. (0->0), fcn. (422->12), ass. (0->42)
t87 = sin(pkin(10));
t112 = pkin(3) * t87;
t111 = pkin(8) + rSges(6,3);
t88 = sin(pkin(9));
t89 = sin(pkin(5));
t110 = t88 * t89;
t91 = cos(pkin(9));
t109 = t89 * t91;
t95 = sin(qJ(2));
t108 = t89 * t95;
t97 = cos(qJ(2));
t107 = t89 * t97;
t92 = cos(pkin(5));
t106 = t92 * t95;
t105 = t92 * t97;
t104 = t91 * pkin(1) + pkin(6) * t110;
t103 = t92 * pkin(6) + qJ(1);
t102 = qJ(3) + rSges(4,3);
t101 = t87 * t110;
t70 = t105 * t88 + t91 * t95;
t71 = -t106 * t88 + t91 * t97;
t90 = cos(pkin(10));
t80 = pkin(3) * t90 + pkin(2);
t93 = -pkin(7) - qJ(3);
t100 = pkin(3) * t101 - t70 * t93 + t71 * t80 + t104;
t99 = t93 * t107 + t80 * t108 + t92 * t112 + t103;
t68 = -t105 * t91 + t88 * t95;
t69 = t106 * t91 + t88 * t97;
t83 = t88 * pkin(1);
t98 = t69 * t80 + t83 + (-pkin(6) - t112) * t109 - t68 * t93;
t96 = cos(qJ(5));
t94 = sin(qJ(5));
t86 = pkin(10) + qJ(4);
t82 = cos(t86);
t81 = sin(t86);
t65 = t108 * t82 + t81 * t92;
t64 = t108 * t81 - t92 * t82;
t61 = t110 * t81 + t71 * t82;
t60 = -t110 * t82 + t71 * t81;
t59 = -t109 * t81 + t69 * t82;
t58 = t109 * t82 + t69 * t81;
t1 = -m(1) * (rSges(1,1) * g(1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t91 - rSges(2,2) * t88) + g(2) * (rSges(2,1) * t88 + rSges(2,2) * t91) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t71 - rSges(3,2) * t70 + t104) + g(2) * (rSges(3,1) * t69 - rSges(3,2) * t68 + t83) + g(3) * (t92 * rSges(3,3) + t103) + (g(1) * rSges(3,3) * t88 + g(3) * (rSges(3,1) * t95 + rSges(3,2) * t97) + g(2) * (-rSges(3,3) - pkin(6)) * t91) * t89) - m(4) * (g(1) * (t71 * pkin(2) + (t71 * t90 + t101) * rSges(4,1) + (t110 * t90 - t71 * t87) * rSges(4,2) + t102 * t70 + t104) + g(2) * (t69 * pkin(2) + t83 - pkin(6) * t109 + (-t109 * t87 + t69 * t90) * rSges(4,1) + (-t109 * t90 - t69 * t87) * rSges(4,2) + t102 * t68) + g(3) * ((rSges(4,1) * t87 + rSges(4,2) * t90) * t92 + (-t102 * t97 + (rSges(4,1) * t90 - rSges(4,2) * t87 + pkin(2)) * t95) * t89 + t103)) - m(5) * (g(1) * (rSges(5,1) * t61 - rSges(5,2) * t60 + rSges(5,3) * t70 + t100) + g(2) * (rSges(5,1) * t59 - rSges(5,2) * t58 + rSges(5,3) * t68 + t98) + g(3) * (rSges(5,1) * t65 - rSges(5,2) * t64 - rSges(5,3) * t107 + t99)) - m(6) * (g(1) * (t61 * pkin(4) + (t61 * t96 + t70 * t94) * rSges(6,1) + (-t61 * t94 + t70 * t96) * rSges(6,2) + t111 * t60 + t100) + g(2) * (t59 * pkin(4) + (t59 * t96 + t68 * t94) * rSges(6,1) + (-t59 * t94 + t68 * t96) * rSges(6,2) + t111 * t58 + t98) + g(3) * (t65 * pkin(4) + (-t107 * t94 + t65 * t96) * rSges(6,1) + (-t107 * t96 - t65 * t94) * rSges(6,2) + t111 * t64 + t99));
U = t1;
