% Calculate potential energy for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:28
% EndTime: 2019-03-09 01:41:28
% DurationCPUTime: 0.36s
% Computational Cost: add. (214->88), mult. (159->103), div. (0->0), fcn. (135->10), ass. (0->34)
t89 = rSges(7,3) + pkin(8);
t63 = pkin(10) + qJ(4);
t56 = sin(t63);
t64 = qJ(1) + pkin(9);
t59 = cos(t64);
t87 = t56 * t59;
t57 = sin(t64);
t68 = sin(qJ(6));
t86 = t57 * t68;
t70 = cos(qJ(6));
t85 = t57 * t70;
t58 = cos(t63);
t84 = t58 * t59;
t83 = t59 * t68;
t82 = t59 * t70;
t81 = pkin(6) + qJ(2);
t66 = cos(pkin(10));
t55 = pkin(3) * t66 + pkin(2);
t71 = cos(qJ(1));
t62 = t71 * pkin(1);
t80 = t59 * t55 + t62;
t79 = qJ(5) * t56;
t78 = rSges(4,3) + qJ(3);
t65 = sin(pkin(10));
t77 = t65 * pkin(3) + t81;
t69 = sin(qJ(1));
t61 = t69 * pkin(1);
t67 = -pkin(7) - qJ(3);
t76 = t57 * t55 + t59 * t67 + t61;
t75 = t56 * pkin(4) + t77;
t74 = pkin(4) * t84 + t59 * t79 + t80;
t73 = t76 + (pkin(4) * t58 + t79) * t57;
t72 = rSges(4,1) * t66 - rSges(4,2) * t65 + pkin(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t71 - t69 * rSges(2,2)) + g(2) * (t69 * rSges(2,1) + rSges(2,2) * t71) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t59 - rSges(3,2) * t57 + t62) + g(2) * (rSges(3,1) * t57 + rSges(3,2) * t59 + t61) + g(3) * (rSges(3,3) + t81)) - m(4) * (g(1) * t62 + g(2) * t61 + g(3) * (t65 * rSges(4,1) + t66 * rSges(4,2) + t81) + (g(1) * t72 - g(2) * t78) * t59 + (g(1) * t78 + g(2) * t72) * t57) - m(5) * (g(1) * (rSges(5,1) * t84 - rSges(5,2) * t87 + t80) + g(2) * (-rSges(5,3) * t59 + t76) + g(3) * (rSges(5,1) * t56 + rSges(5,2) * t58 + t77) + (g(1) * (rSges(5,3) - t67) + g(2) * (rSges(5,1) * t58 - rSges(5,2) * t56)) * t57) - m(6) * (g(1) * (-rSges(6,2) * t84 + rSges(6,3) * t87 + t74) + g(2) * (-rSges(6,1) * t59 + t73) + g(3) * (-rSges(6,2) * t56 + (-rSges(6,3) - qJ(5)) * t58 + t75) + (g(1) * (rSges(6,1) - t67) + g(2) * (-rSges(6,2) * t58 + rSges(6,3) * t56)) * t57) - m(7) * (g(1) * ((t56 * t83 + t85) * rSges(7,1) + (t56 * t82 - t86) * rSges(7,2) + t74 + (pkin(5) - t67) * t57) + g(2) * (-t59 * pkin(5) + (t56 * t86 - t82) * rSges(7,1) + (t56 * t85 + t83) * rSges(7,2) + t73) + g(3) * (t89 * t56 + t75) + (g(3) * (-rSges(7,1) * t68 - rSges(7,2) * t70 - qJ(5)) + (g(1) * t59 + g(2) * t57) * t89) * t58);
U  = t1;
