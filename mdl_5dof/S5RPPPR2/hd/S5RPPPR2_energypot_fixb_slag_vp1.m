% Calculate potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:53
% EndTime: 2019-12-05 17:30:54
% DurationCPUTime: 0.42s
% Computational Cost: add. (158->94), mult. (313->124), div. (0->0), fcn. (348->10), ass. (0->43)
t106 = pkin(6) + rSges(6,3);
t85 = sin(qJ(1));
t108 = g(2) * t85;
t80 = sin(pkin(7));
t107 = t80 * pkin(2) + pkin(5);
t79 = sin(pkin(8));
t105 = t79 * t80;
t82 = cos(pkin(8));
t104 = t80 * t82;
t87 = cos(qJ(1));
t103 = t80 * t87;
t83 = cos(pkin(7));
t102 = t83 * t87;
t101 = t85 * t79;
t100 = t85 * t80;
t99 = t85 * t82;
t98 = t87 * t79;
t97 = t87 * t82;
t96 = t87 * pkin(1) + t85 * qJ(2);
t95 = t80 * qJ(3);
t94 = -rSges(4,3) - qJ(3);
t93 = -t83 * pkin(2) - pkin(1);
t92 = pkin(2) * t102 + t87 * t95 + t96;
t64 = t83 * t101 + t97;
t65 = -t83 * t99 + t98;
t76 = t87 * qJ(2);
t91 = t65 * pkin(3) - t64 * qJ(4) + t76;
t67 = t83 * t97 + t101;
t90 = t67 * pkin(3) + t92;
t89 = pkin(3) * t104 - t83 * qJ(3) + qJ(4) * t105 + t107;
t88 = (t93 - t95) * t108;
t86 = cos(qJ(5));
t84 = sin(qJ(5));
t81 = cos(pkin(9));
t78 = sin(pkin(9));
t66 = t83 * t98 - t99;
t63 = t81 * t104 - t83 * t78;
t62 = t78 * t104 + t83 * t81;
t59 = t78 * t103 + t67 * t81;
t58 = -t81 * t103 + t67 * t78;
t57 = -t78 * t100 + t65 * t81;
t56 = t81 * t100 + t65 * t78;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t85 * rSges(2,1) - t87 * rSges(2,2)) + g(3) * (t87 * rSges(2,1) - t85 * rSges(2,2))) - m(3) * (g(1) * (t80 * rSges(3,1) + t83 * rSges(3,2) + pkin(5)) + g(2) * (t87 * rSges(3,3) + t76) + g(3) * (rSges(3,1) * t102 - rSges(3,2) * t103 + t96) + (g(2) * (-rSges(3,1) * t83 + rSges(3,2) * t80 - pkin(1)) + g(3) * rSges(3,3)) * t85) - m(4) * (g(1) * (t94 * t83 + (rSges(4,1) * t82 - rSges(4,2) * t79) * t80 + t107) + g(2) * (t65 * rSges(4,1) + t64 * rSges(4,2) + t76) + g(3) * (t67 * rSges(4,1) - t66 * rSges(4,2) + rSges(4,3) * t103 + t92) + (t94 * t80 + t93) * t108) - m(5) * (g(1) * (t63 * rSges(5,1) - t62 * rSges(5,2) + rSges(5,3) * t105 + t89) + g(2) * (t57 * rSges(5,1) - t56 * rSges(5,2) - t64 * rSges(5,3) + t91) + g(3) * (t59 * rSges(5,1) - t58 * rSges(5,2) + (rSges(5,3) + qJ(4)) * t66 + t90) + t88) - m(6) * (g(1) * (t63 * pkin(4) + (t84 * t105 + t63 * t86) * rSges(6,1) + (t86 * t105 - t63 * t84) * rSges(6,2) + t106 * t62 + t89) + g(2) * (t57 * pkin(4) + (t57 * t86 - t64 * t84) * rSges(6,1) + (-t57 * t84 - t64 * t86) * rSges(6,2) + t91 + t106 * t56) + g(3) * (t59 * pkin(4) + t66 * qJ(4) + (t59 * t86 + t66 * t84) * rSges(6,1) + (-t59 * t84 + t66 * t86) * rSges(6,2) + t106 * t58 + t90) + t88);
U = t1;
