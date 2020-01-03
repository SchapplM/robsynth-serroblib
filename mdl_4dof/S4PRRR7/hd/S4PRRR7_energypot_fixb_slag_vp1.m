% Calculate potential energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR7_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:00
% EndTime: 2019-12-31 16:36:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (130->75), mult. (272->108), div. (0->0), fcn. (311->10), ass. (0->36)
t68 = sin(pkin(4));
t90 = pkin(5) * t68;
t89 = rSges(4,3) + pkin(6);
t88 = pkin(7) + rSges(5,3);
t72 = sin(qJ(3));
t87 = t68 * t72;
t73 = sin(qJ(2));
t86 = t68 * t73;
t75 = cos(qJ(3));
t85 = t68 * t75;
t76 = cos(qJ(2));
t84 = t68 * t76;
t70 = cos(pkin(4));
t83 = t70 * t73;
t82 = t70 * t76;
t67 = sin(pkin(8));
t69 = cos(pkin(8));
t81 = t69 * pkin(1) + t67 * t90;
t80 = t70 * pkin(5) + qJ(1);
t56 = -t67 * t83 + t69 * t76;
t79 = t56 * pkin(2) + t81;
t78 = pkin(2) * t86 + t80;
t54 = t67 * t76 + t69 * t83;
t64 = t67 * pkin(1);
t77 = t54 * pkin(2) - t69 * t90 + t64;
t74 = cos(qJ(4));
t71 = sin(qJ(4));
t58 = t70 * t72 + t73 * t85;
t57 = -t70 * t75 + t72 * t86;
t55 = t67 * t82 + t69 * t73;
t53 = t67 * t73 - t69 * t82;
t50 = t56 * t75 + t67 * t87;
t49 = t56 * t72 - t67 * t85;
t48 = t54 * t75 - t69 * t87;
t47 = t54 * t72 + t69 * t85;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t69 * rSges(2,1) - t67 * rSges(2,2)) + g(2) * (t67 * rSges(2,1) + t69 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t56 * rSges(3,1) - t55 * rSges(3,2) + t81) + g(2) * (t54 * rSges(3,1) - t53 * rSges(3,2) + t64) + g(3) * (t70 * rSges(3,3) + t80) + (g(1) * rSges(3,3) * t67 + g(3) * (rSges(3,1) * t73 + rSges(3,2) * t76) + g(2) * (-rSges(3,3) - pkin(5)) * t69) * t68) - m(4) * (g(1) * (t50 * rSges(4,1) - t49 * rSges(4,2) + t89 * t55 + t79) + g(2) * (t48 * rSges(4,1) - t47 * rSges(4,2) + t89 * t53 + t77) + g(3) * (t58 * rSges(4,1) - t57 * rSges(4,2) - t89 * t84 + t78)) - m(5) * (g(1) * (t50 * pkin(3) + t55 * pkin(6) + (t50 * t74 + t55 * t71) * rSges(5,1) + (-t50 * t71 + t55 * t74) * rSges(5,2) + t88 * t49 + t79) + g(2) * (t48 * pkin(3) + t53 * pkin(6) + (t48 * t74 + t53 * t71) * rSges(5,1) + (-t48 * t71 + t53 * t74) * rSges(5,2) + t88 * t47 + t77) + g(3) * (t58 * pkin(3) - pkin(6) * t84 + (t58 * t74 - t71 * t84) * rSges(5,1) + (-t58 * t71 - t74 * t84) * rSges(5,2) + t88 * t57 + t78));
U = t1;
