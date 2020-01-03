% Calculate potential energy for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:32
% EndTime: 2019-12-31 17:47:33
% DurationCPUTime: 0.36s
% Computational Cost: add. (116->85), mult. (203->110), div. (0->0), fcn. (199->8), ass. (0->31)
t68 = sin(qJ(1));
t70 = cos(qJ(1));
t61 = t68 * pkin(1);
t64 = sin(pkin(7));
t77 = qJ(3) * t64;
t66 = cos(pkin(7));
t84 = t66 * t68;
t74 = pkin(2) * t84 + t68 * t77 + t61;
t76 = qJ(4) * t66;
t88 = t68 * t76 + t74 + (-pkin(3) - qJ(2)) * t70;
t87 = t64 * pkin(2) + pkin(5);
t86 = pkin(6) + rSges(6,3);
t85 = t64 * t68;
t83 = t66 * t70;
t63 = sin(pkin(8));
t82 = t68 * t63;
t65 = cos(pkin(8));
t81 = t68 * t65;
t80 = t70 * t63;
t79 = t70 * t65;
t78 = t70 * pkin(1) + t68 * qJ(2);
t75 = t64 * qJ(4) + t87;
t72 = pkin(2) * t83 + t70 * t77 + t78;
t71 = t68 * pkin(3) + t70 * t76 + t72;
t69 = cos(qJ(5));
t67 = sin(qJ(5));
t49 = t64 * t82 - t79;
t48 = t64 * t81 + t80;
t47 = t64 * t80 + t81;
t46 = -t64 * t79 + t82;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t70 * rSges(2,1) - t68 * rSges(2,2)) + g(2) * (t68 * rSges(2,1) + t70 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t68 * rSges(3,3) + t78) + g(2) * (rSges(3,1) * t84 - rSges(3,2) * t85 + t61) + g(3) * (t64 * rSges(3,1) + t66 * rSges(3,2) + pkin(5)) + (g(1) * (rSges(3,1) * t66 - rSges(3,2) * t64) + g(2) * (-rSges(3,3) - qJ(2))) * t70) - m(4) * (g(1) * (t68 * rSges(4,1) + t72) + g(2) * (-rSges(4,2) * t84 + rSges(4,3) * t85 + t74) + g(3) * (-t64 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t66 + t87) + (g(1) * (-rSges(4,2) * t66 + rSges(4,3) * t64) + g(2) * (-rSges(4,1) - qJ(2))) * t70) - m(5) * (g(1) * (t47 * rSges(5,1) - t46 * rSges(5,2) + t71) + g(2) * (t49 * rSges(5,1) + t48 * rSges(5,2) + t88) + g(3) * (t64 * rSges(5,3) + t75) + (g(3) * (-rSges(5,1) * t63 - rSges(5,2) * t65 - qJ(3)) + (g(1) * t70 + g(2) * t68) * rSges(5,3)) * t66) - m(6) * (g(1) * (t47 * pkin(4) + (t47 * t69 + t67 * t83) * rSges(6,1) + (-t47 * t67 + t69 * t83) * rSges(6,2) + t86 * t46 + t71) + g(2) * (t49 * pkin(4) + (t49 * t69 + t67 * t84) * rSges(6,1) + (-t49 * t67 + t69 * t84) * rSges(6,2) - t86 * t48 + t88) + g(3) * ((t67 * rSges(6,1) + t69 * rSges(6,2)) * t64 + (-qJ(3) + t86 * t65 + (-t69 * rSges(6,1) + t67 * rSges(6,2) - pkin(4)) * t63) * t66 + t75));
U = t1;
