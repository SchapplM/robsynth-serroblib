% Calculate potential energy for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:07
% EndTime: 2019-03-09 03:21:07
% DurationCPUTime: 0.40s
% Computational Cost: add. (163->92), mult. (180->107), div. (0->0), fcn. (160->8), ass. (0->37)
t85 = rSges(6,3) + pkin(8);
t84 = rSges(7,3) + qJ(6) + pkin(8);
t83 = pkin(2) + pkin(6);
t60 = sin(qJ(3));
t82 = pkin(3) * t60;
t61 = sin(qJ(1));
t81 = g(1) * t61;
t53 = t61 * pkin(1);
t80 = g(2) * t53;
t79 = rSges(4,3) + pkin(7);
t56 = qJ(3) + pkin(9);
t50 = sin(t56);
t77 = t50 * t61;
t59 = sin(qJ(5));
t64 = cos(qJ(1));
t76 = t59 * t64;
t75 = t61 * t59;
t62 = cos(qJ(5));
t74 = t61 * t62;
t73 = t62 * t64;
t58 = -qJ(4) - pkin(7);
t72 = rSges(5,3) - t58;
t71 = t64 * pkin(1) + t61 * qJ(2);
t63 = cos(qJ(3));
t70 = t63 * pkin(3) + t83;
t69 = t61 * t82 + t71;
t68 = -qJ(2) - t82;
t67 = -t61 * t58 + t53;
t66 = rSges(4,1) * t60 + rSges(4,2) * t63;
t51 = cos(t56);
t65 = rSges(5,1) * t50 + rSges(5,2) * t51;
t49 = pkin(5) * t62 + pkin(4);
t47 = -t50 * t73 + t75;
t46 = t50 * t76 + t74;
t45 = t50 * t74 + t76;
t44 = -t50 * t75 + t73;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t64 - t61 * rSges(2,2)) + g(2) * (t61 * rSges(2,1) + rSges(2,2) * t64) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t64 + t61 * rSges(3,3) + t71) + g(2) * (-t61 * rSges(3,2) + t53 + (-rSges(3,3) - qJ(2)) * t64) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * t71 + t80 + g(3) * (rSges(4,1) * t63 - rSges(4,2) * t60 + t83) + (g(1) * t66 + g(2) * t79) * t61 + (g(1) * t79 + g(2) * (-qJ(2) - t66)) * t64) - m(5) * (g(1) * t69 + t80 + g(3) * (rSges(5,1) * t51 - rSges(5,2) * t50 + t70) + (g(1) * t65 + g(2) * t72) * t61 + (g(1) * t72 + g(2) * (-t65 + t68)) * t64) - m(6) * (g(1) * (t45 * rSges(6,1) + t44 * rSges(6,2) + pkin(4) * t77 + t69) + g(2) * (t47 * rSges(6,1) + t46 * rSges(6,2) + t67) + g(3) * (t85 * t50 + t70) + (-g(1) * t58 + g(2) * (-pkin(4) * t50 + t68)) * t64 + (g(3) * (rSges(6,1) * t62 - rSges(6,2) * t59 + pkin(4)) + (g(2) * t64 - t81) * t85) * t51) - m(7) * (g(1) * (t45 * rSges(7,1) + t44 * rSges(7,2) + t49 * t77 + t69) + g(2) * (t47 * rSges(7,1) + t46 * rSges(7,2) + pkin(5) * t75 + t67) + g(3) * (t84 * t50 + t70) + (g(3) * (rSges(7,1) * t62 - rSges(7,2) * t59 + t49) - t84 * t81) * t51 + (g(1) * (pkin(5) * t59 - t58) + g(2) * (-t49 * t50 + t84 * t51 + t68)) * t64);
U  = t1;
