% Calculate potential energy for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:29
% EndTime: 2019-12-05 15:10:29
% DurationCPUTime: 0.23s
% Computational Cost: add. (145->83), mult. (279->108), div. (0->0), fcn. (301->8), ass. (0->42)
t112 = rSges(6,1) + pkin(4);
t111 = rSges(6,2) + pkin(6);
t110 = rSges(5,3) + pkin(6);
t83 = sin(pkin(8));
t84 = sin(pkin(7));
t109 = t83 * t84;
t86 = cos(pkin(7));
t108 = t83 * t86;
t87 = sin(qJ(4));
t107 = t83 * t87;
t88 = sin(qJ(3));
t106 = t83 * t88;
t89 = cos(qJ(4));
t105 = t83 * t89;
t90 = cos(qJ(3));
t104 = t83 * t90;
t85 = cos(pkin(8));
t103 = t84 * t85;
t102 = t84 * t88;
t101 = t84 * t90;
t100 = t86 * t88;
t99 = t86 * t90;
t98 = t86 * pkin(1) + t84 * qJ(2);
t97 = rSges(6,3) + qJ(5);
t96 = t83 * pkin(2) + qJ(1);
t95 = t86 * t85 * pkin(2) + pkin(5) * t108 + t98;
t67 = t85 * t99 + t102;
t94 = t67 * pkin(3) + t95;
t81 = t84 * pkin(1);
t93 = pkin(2) * t103 + pkin(5) * t109 - t86 * qJ(2) + t81;
t92 = pkin(3) * t104 - t85 * pkin(5) + pkin(6) * t106 + t96;
t65 = t85 * t101 - t100;
t91 = t65 * pkin(3) + t93;
t69 = t89 * t104 - t85 * t87;
t68 = t87 * t104 + t85 * t89;
t66 = t85 * t100 - t101;
t64 = t85 * t102 + t99;
t61 = t86 * t107 + t67 * t89;
t60 = -t86 * t105 + t67 * t87;
t59 = t84 * t107 + t65 * t89;
t58 = -t84 * t105 + t65 * t87;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t86 * rSges(2,1) - t84 * rSges(2,2)) + g(2) * (t84 * rSges(2,1) + t86 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t84 * rSges(3,3) + t98) + g(2) * (rSges(3,1) * t103 - rSges(3,2) * t109 + t81) + g(3) * (t83 * rSges(3,1) + t85 * rSges(3,2) + qJ(1)) + (g(1) * (rSges(3,1) * t85 - rSges(3,2) * t83) + g(2) * (-rSges(3,3) - qJ(2))) * t86) - m(4) * (g(1) * (t67 * rSges(4,1) - t66 * rSges(4,2) + rSges(4,3) * t108 + t95) + g(2) * (t65 * rSges(4,1) - t64 * rSges(4,2) + rSges(4,3) * t109 + t93) + g(3) * ((-rSges(4,3) - pkin(5)) * t85 + (rSges(4,1) * t90 - rSges(4,2) * t88) * t83 + t96)) - m(5) * (g(1) * (t61 * rSges(5,1) - t60 * rSges(5,2) + t110 * t66 + t94) + g(2) * (t59 * rSges(5,1) - t58 * rSges(5,2) + t110 * t64 + t91) + g(3) * (t69 * rSges(5,1) - t68 * rSges(5,2) + rSges(5,3) * t106 + t92)) - m(6) * (g(1) * (t111 * t66 + t112 * t61 + t97 * t60 + t94) + g(2) * (t111 * t64 + t112 * t59 + t97 * t58 + t91) + g(3) * (rSges(6,2) * t106 + t112 * t69 + t97 * t68 + t92));
U = t1;
