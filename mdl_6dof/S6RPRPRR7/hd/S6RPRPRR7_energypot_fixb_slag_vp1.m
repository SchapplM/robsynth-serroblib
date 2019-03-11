% Calculate potential energy for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:52
% EndTime: 2019-03-09 03:54:53
% DurationCPUTime: 0.30s
% Computational Cost: add. (172->80), mult. (156->86), div. (0->0), fcn. (132->10), ass. (0->36)
t84 = rSges(7,3) + pkin(9);
t57 = qJ(3) + pkin(10);
t51 = qJ(5) + t57;
t47 = sin(t51);
t59 = sin(qJ(6));
t62 = cos(qJ(6));
t69 = rSges(7,1) * t62 - rSges(7,2) * t59 + pkin(5);
t83 = t47 * t69;
t82 = pkin(2) + pkin(6);
t60 = sin(qJ(3));
t81 = pkin(3) * t60;
t61 = sin(qJ(1));
t53 = t61 * pkin(1);
t80 = g(2) * t53;
t79 = rSges(4,3) + pkin(7);
t58 = -qJ(4) - pkin(7);
t77 = rSges(5,3) - t58;
t56 = -pkin(8) + t58;
t76 = rSges(6,3) - t56;
t64 = cos(qJ(1));
t75 = t64 * pkin(1) + t61 * qJ(2);
t49 = sin(t57);
t45 = pkin(4) * t49 + t81;
t74 = -qJ(2) - t45;
t63 = cos(qJ(3));
t73 = t63 * pkin(3) + t82;
t50 = cos(t57);
t72 = pkin(4) * t50 + t73;
t71 = rSges(4,1) * t60 + rSges(4,2) * t63;
t48 = cos(t51);
t70 = rSges(6,1) * t47 + rSges(6,2) * t48;
t68 = rSges(7,1) * t59 + rSges(7,2) * t62 - t56;
t67 = g(1) * t75 + t80;
t66 = g(1) * (t61 * t45 + t75) + t80;
t65 = rSges(5,1) * t49 + rSges(5,2) * t50 + t81;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t64 - t61 * rSges(2,2)) + g(2) * (t61 * rSges(2,1) + rSges(2,2) * t64) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t64 + t61 * rSges(3,3) + t75) + g(2) * (-t61 * rSges(3,2) + t53 + (-rSges(3,3) - qJ(2)) * t64) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(3) * (rSges(4,1) * t63 - rSges(4,2) * t60 + t82) + (g(1) * t71 + g(2) * t79) * t61 + (g(1) * t79 + g(2) * (-qJ(2) - t71)) * t64 + t67) - m(5) * (g(3) * (rSges(5,1) * t50 - rSges(5,2) * t49 + t73) + (g(1) * t65 + g(2) * t77) * t61 + (g(1) * t77 + g(2) * (-qJ(2) - t65)) * t64 + t67) - m(6) * (g(3) * (rSges(6,1) * t48 - rSges(6,2) * t47 + t72) + (g(1) * t70 + g(2) * t76) * t61 + (g(1) * t76 + g(2) * (-t70 + t74)) * t64 + t66) - m(7) * (g(3) * (t84 * t47 + t72) + (g(1) * t83 + g(2) * t68) * t61 + (g(1) * t68 + (t74 - t83) * g(2)) * t64 + (g(3) * t69 + (-g(1) * t61 + g(2) * t64) * t84) * t48 + t66);
U  = t1;
