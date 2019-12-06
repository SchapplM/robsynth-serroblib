% Calculate potential energy for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:45
% EndTime: 2019-12-05 14:57:46
% DurationCPUTime: 0.49s
% Computational Cost: add. (162->83), mult. (216->105), div. (0->0), fcn. (216->10), ass. (0->34)
t98 = pkin(6) + rSges(6,3);
t75 = sin(qJ(5));
t76 = cos(qJ(5));
t97 = -t75 * rSges(6,1) - t76 * rSges(6,2);
t70 = sin(pkin(7));
t73 = cos(pkin(7));
t96 = g(1) * t73 + g(2) * t70;
t95 = rSges(4,3) + qJ(3);
t69 = sin(pkin(8));
t92 = rSges(3,2) * t69;
t68 = sin(pkin(9));
t91 = t68 * t70;
t90 = t68 * t73;
t72 = cos(pkin(8));
t89 = t70 * t72;
t88 = t72 * t73;
t84 = t73 * pkin(1) + t70 * qJ(2);
t71 = cos(pkin(9));
t61 = pkin(3) * t71 + pkin(2);
t74 = -pkin(5) - qJ(3);
t82 = t69 * t61 + t72 * t74 + qJ(1);
t65 = t70 * pkin(1);
t81 = -t73 * qJ(2) + t65;
t80 = pkin(3) * t91 + t61 * t88 + t84;
t79 = t76 * rSges(6,1) - t75 * rSges(6,2) + pkin(4);
t77 = -pkin(3) * t90 + t61 * t89 + t81;
t67 = pkin(9) + qJ(4);
t63 = cos(t67);
t62 = sin(t67);
t54 = t62 * t70 + t63 * t88;
t53 = t62 * t88 - t70 * t63;
t52 = -t62 * t73 + t63 * t89;
t51 = t62 * t89 + t63 * t73;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t73 - rSges(2,2) * t70) + g(2) * (rSges(2,1) * t70 + rSges(2,2) * t73) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t70 + t84) + g(2) * (rSges(3,1) * t89 - t70 * t92 + t65) + g(3) * (rSges(3,1) * t69 + rSges(3,2) * t72 + qJ(1)) + (g(1) * (rSges(3,1) * t72 - t92) + g(2) * (-rSges(3,3) - qJ(2))) * t73) - m(4) * (g(1) * (pkin(2) * t88 + (t71 * t88 + t91) * rSges(4,1) + (-t68 * t88 + t70 * t71) * rSges(4,2) + t84) + g(2) * (pkin(2) * t89 + (t71 * t89 - t90) * rSges(4,1) + (-t68 * t89 - t71 * t73) * rSges(4,2) + t81) + g(3) * (-t95 * t72 + qJ(1)) + (g(3) * (rSges(4,1) * t71 - rSges(4,2) * t68 + pkin(2)) + t96 * t95) * t69) - m(5) * (g(1) * (rSges(5,1) * t54 - rSges(5,2) * t53 + t80) + g(2) * (rSges(5,1) * t52 - rSges(5,2) * t51 + t77) + g(3) * (-rSges(5,3) * t72 + t82) + (g(3) * (rSges(5,1) * t63 - rSges(5,2) * t62) + t96 * (rSges(5,3) - t74)) * t69) - m(6) * (g(3) * (t97 * t72 + t82) + (g(3) * (t98 * t62 + t79 * t63) + t96 * (-t74 - t97)) * t69 + (t98 * t51 + t79 * t52 + t77) * g(2) + (t98 * t53 + t79 * t54 + t80) * g(1));
U = t1;
