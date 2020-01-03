% Calculate potential energy for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP11_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:17
% EndTime: 2019-12-31 20:12:18
% DurationCPUTime: 0.34s
% Computational Cost: add. (108->76), mult. (182->95), div. (0->0), fcn. (170->6), ass. (0->29)
t83 = rSges(6,1) + pkin(4);
t82 = rSges(6,3) + qJ(5);
t62 = sin(qJ(2));
t81 = t62 * pkin(2) + pkin(5);
t61 = sin(qJ(4));
t66 = cos(qJ(1));
t80 = t61 * t66;
t63 = sin(qJ(1));
t79 = t62 * t63;
t64 = cos(qJ(4));
t78 = t63 * t64;
t65 = cos(qJ(2));
t77 = t63 * t65;
t76 = t64 * t66;
t75 = t65 * t66;
t74 = t66 * pkin(1) + t63 * pkin(6);
t73 = qJ(3) * t62;
t72 = t62 * pkin(7) + t81;
t59 = t63 * pkin(1);
t71 = pkin(2) * t77 + t63 * t73 + t59;
t70 = pkin(2) * t75 + t66 * t73 + t74;
t69 = g(1) * t66 + g(2) * t63;
t68 = t63 * pkin(3) + pkin(7) * t75 + t70;
t67 = pkin(7) * t77 + t71 + (-pkin(3) - pkin(6)) * t66;
t47 = t61 * t79 - t76;
t46 = t62 * t78 + t80;
t45 = t62 * t80 + t78;
t44 = t61 * t63 - t62 * t76;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t66 - t63 * rSges(2,2)) + g(2) * (t63 * rSges(2,1) + rSges(2,2) * t66) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t63 * rSges(3,3) + t74) + g(2) * (rSges(3,1) * t77 - rSges(3,2) * t79 + t59) + g(3) * (rSges(3,1) * t62 + rSges(3,2) * t65 + pkin(5)) + (g(1) * (rSges(3,1) * t65 - rSges(3,2) * t62) + g(2) * (-rSges(3,3) - pkin(6))) * t66) - m(4) * (g(1) * (t63 * rSges(4,1) + t70) + g(2) * (-rSges(4,2) * t77 + rSges(4,3) * t79 + t71) + g(3) * (-t62 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t65 + t81) + (g(1) * (-rSges(4,2) * t65 + rSges(4,3) * t62) + g(2) * (-rSges(4,1) - pkin(6))) * t66) - m(5) * (g(1) * (t45 * rSges(5,1) - t44 * rSges(5,2) + t68) + g(2) * (t47 * rSges(5,1) + t46 * rSges(5,2) + t67) + g(3) * (t62 * rSges(5,3) + t72) + (g(3) * (-rSges(5,1) * t61 - rSges(5,2) * t64 - qJ(3)) + t69 * rSges(5,3)) * t65) - m(6) * (g(1) * (t82 * t44 + t83 * t45 + t68) + g(2) * (-t82 * t46 + t83 * t47 + t67) + g(3) * (t62 * rSges(6,2) + t72) + (g(3) * (-t83 * t61 + t82 * t64 - qJ(3)) + t69 * rSges(6,2)) * t65);
U = t1;
