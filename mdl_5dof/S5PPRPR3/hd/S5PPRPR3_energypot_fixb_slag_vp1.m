% Calculate potential energy for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:41
% EndTime: 2019-12-05 15:04:42
% DurationCPUTime: 0.50s
% Computational Cost: add. (162->83), mult. (216->105), div. (0->0), fcn. (216->10), ass. (0->36)
t102 = rSges(4,3) + pkin(5);
t101 = pkin(6) + rSges(6,3);
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t100 = -t75 * rSges(6,1) - t77 * rSges(6,2);
t71 = sin(pkin(7));
t73 = cos(pkin(7));
t99 = g(1) * t73 + g(2) * t71;
t70 = sin(pkin(8));
t95 = rSges(3,2) * t70;
t72 = cos(pkin(8));
t94 = t71 * t72;
t76 = sin(qJ(3));
t93 = t71 * t76;
t78 = cos(qJ(3));
t92 = t71 * t78;
t91 = t72 * t73;
t90 = t73 * t76;
t89 = t73 * t78;
t85 = t73 * pkin(1) + t71 * qJ(2);
t63 = pkin(3) * t78 + pkin(2);
t74 = -qJ(4) - pkin(5);
t84 = t70 * t63 + t72 * t74 + qJ(1);
t67 = t71 * pkin(1);
t83 = -t73 * qJ(2) + t67;
t82 = pkin(3) * t93 + t63 * t91 + t85;
t81 = rSges(6,1) * t77 - rSges(6,2) * t75 + pkin(4);
t79 = -pkin(3) * t90 + t63 * t94 + t83;
t69 = qJ(3) + pkin(9);
t65 = cos(t69);
t64 = sin(t69);
t56 = t64 * t71 + t65 * t91;
t55 = t64 * t91 - t71 * t65;
t54 = -t64 * t73 + t65 * t94;
t53 = t64 * t94 + t65 * t73;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t73 - rSges(2,2) * t71) + g(2) * (rSges(2,1) * t71 + rSges(2,2) * t73) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t71 + t85) + g(2) * (rSges(3,1) * t94 - t71 * t95 + t67) + g(3) * (rSges(3,1) * t70 + rSges(3,2) * t72 + qJ(1)) + (g(1) * (rSges(3,1) * t72 - t95) + g(2) * (-rSges(3,3) - qJ(2))) * t73) - m(4) * (g(1) * (pkin(2) * t91 + (t72 * t89 + t93) * rSges(4,1) + (-t72 * t90 + t92) * rSges(4,2) + t85) + g(2) * (pkin(2) * t94 + (t72 * t92 - t90) * rSges(4,1) + (-t72 * t93 - t89) * rSges(4,2) + t83) + g(3) * (-t102 * t72 + qJ(1)) + (g(3) * (rSges(4,1) * t78 - rSges(4,2) * t76 + pkin(2)) + t99 * t102) * t70) - m(5) * (g(1) * (rSges(5,1) * t56 - rSges(5,2) * t55 + t82) + g(2) * (rSges(5,1) * t54 - rSges(5,2) * t53 + t79) + g(3) * (-rSges(5,3) * t72 + t84) + (g(3) * (rSges(5,1) * t65 - rSges(5,2) * t64) + t99 * (rSges(5,3) - t74)) * t70) - m(6) * (g(3) * (t100 * t72 + t84) + (g(3) * (t101 * t64 + t81 * t65) + t99 * (-t74 - t100)) * t70 + (t101 * t53 + t81 * t54 + t79) * g(2) + (t101 * t55 + t81 * t56 + t82) * g(1));
U = t1;
