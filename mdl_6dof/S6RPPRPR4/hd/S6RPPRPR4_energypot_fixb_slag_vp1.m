% Calculate potential energy for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:11
% EndTime: 2019-03-09 01:46:11
% DurationCPUTime: 0.36s
% Computational Cost: add. (175->81), mult. (229->99), div. (0->0), fcn. (247->10), ass. (0->31)
t87 = rSges(7,3) + pkin(8);
t62 = qJ(4) + pkin(10);
t56 = cos(t62);
t86 = t56 * pkin(5);
t85 = rSges(5,3) + pkin(7);
t83 = cos(qJ(1));
t82 = sin(qJ(1));
t64 = sin(qJ(6));
t81 = t56 * t64;
t66 = cos(qJ(6));
t80 = t56 * t66;
t79 = pkin(6) - qJ(3);
t78 = t83 * pkin(1) + t82 * qJ(2);
t77 = cos(pkin(9));
t76 = sin(pkin(9));
t75 = t83 * pkin(2) + t78;
t65 = sin(qJ(4));
t74 = -pkin(4) * t65 + t79;
t49 = -t82 * t76 - t83 * t77;
t50 = t83 * t76 - t82 * t77;
t67 = cos(qJ(4));
t54 = pkin(4) * t67 + pkin(3);
t63 = -qJ(5) - pkin(7);
t73 = -t49 * t54 - t50 * t63 + t75;
t72 = t82 * pkin(1) - t83 * qJ(2);
t55 = sin(t62);
t71 = -rSges(6,1) * t56 + rSges(6,2) * t55;
t70 = -rSges(5,1) * t67 + rSges(5,2) * t65 - pkin(3);
t69 = t82 * pkin(2) + t72;
t68 = t49 * t63 - t50 * t54 + t69;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t83 * rSges(2,1) - t82 * rSges(2,2)) + g(2) * (t82 * rSges(2,1) + t83 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t83 * rSges(3,1) + t82 * rSges(3,3) + t78) + g(2) * (t82 * rSges(3,1) - t83 * rSges(3,3) + t72) + g(3) * (pkin(6) + rSges(3,2))) - m(4) * (g(1) * (-rSges(4,1) * t49 - rSges(4,2) * t50 + t75) + g(2) * (-t50 * rSges(4,1) + t49 * rSges(4,2) + t69) + g(3) * (-rSges(4,3) + t79)) - m(5) * (g(1) * t75 + g(2) * t69 + g(3) * (-t65 * rSges(5,1) - rSges(5,2) * t67 + t79) + (g(1) * t85 + g(2) * t70) * t50 + (g(1) * t70 - g(2) * t85) * t49) - m(6) * (g(1) * (rSges(6,3) * t50 + t71 * t49 + t73) + g(2) * (-t49 * rSges(6,3) + t71 * t50 + t68) + g(3) * (-rSges(6,1) * t55 - rSges(6,2) * t56 + t74)) - m(7) * (g(1) * (-t49 * t86 + (-t49 * t80 + t50 * t64) * rSges(7,1) + (t49 * t81 + t50 * t66) * rSges(7,2) + t73) + g(2) * (-t50 * t86 + (-t49 * t64 - t50 * t80) * rSges(7,1) + (-t49 * t66 + t50 * t81) * rSges(7,2) + t68) + g(3) * (t87 * t56 + t74) + (g(3) * (-rSges(7,1) * t66 + rSges(7,2) * t64 - pkin(5)) - (g(1) * t49 + g(2) * t50) * t87) * t55);
U  = t1;
