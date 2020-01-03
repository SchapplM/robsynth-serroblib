% Calculate potential energy for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:51
% EndTime: 2019-12-31 19:59:51
% DurationCPUTime: 0.27s
% Computational Cost: add. (149->73), mult. (167->91), div. (0->0), fcn. (155->8), ass. (0->33)
t84 = rSges(6,1) + pkin(4);
t83 = rSges(6,3) + qJ(5);
t82 = rSges(3,3) + pkin(6);
t64 = sin(qJ(2));
t81 = t64 * pkin(2) + pkin(5);
t61 = qJ(2) + pkin(8);
t58 = sin(t61);
t65 = sin(qJ(1));
t80 = t58 * t65;
t68 = cos(qJ(1));
t79 = t58 * t68;
t59 = cos(t61);
t78 = t59 * t68;
t63 = sin(qJ(4));
t77 = t65 * t63;
t66 = cos(qJ(4));
t76 = t65 * t66;
t75 = t68 * t63;
t74 = t68 * t66;
t67 = cos(qJ(2));
t57 = t67 * pkin(2) + pkin(1);
t62 = -qJ(3) - pkin(6);
t73 = t65 * t57 + t68 * t62;
t72 = t58 * pkin(3) + t81;
t71 = t65 * t59 * pkin(3) + pkin(7) * t80 + t73;
t53 = t68 * t57;
t70 = pkin(3) * t78 + pkin(7) * t79 - t65 * t62 + t53;
t69 = rSges(3,1) * t67 - rSges(3,2) * t64 + pkin(1);
t47 = t59 * t74 + t77;
t46 = t59 * t75 - t76;
t45 = t59 * t76 - t75;
t44 = t59 * t77 + t74;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t68 * rSges(2,1) - t65 * rSges(2,2)) + g(2) * (t65 * rSges(2,1) + t68 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (t64 * rSges(3,1) + t67 * rSges(3,2) + pkin(5)) + (g(1) * t69 - g(2) * t82) * t68 + (g(1) * t82 + g(2) * t69) * t65) - m(4) * (g(1) * (rSges(4,1) * t78 - rSges(4,2) * t79 + t53) + g(2) * (-t68 * rSges(4,3) + t73) + g(3) * (t58 * rSges(4,1) + t59 * rSges(4,2) + t81) + (g(1) * (rSges(4,3) - t62) + g(2) * (rSges(4,1) * t59 - rSges(4,2) * t58)) * t65) - m(5) * (g(1) * (t47 * rSges(5,1) - t46 * rSges(5,2) + rSges(5,3) * t79 + t70) + g(2) * (t45 * rSges(5,1) - t44 * rSges(5,2) + rSges(5,3) * t80 + t71) + g(3) * ((-rSges(5,3) - pkin(7)) * t59 + (rSges(5,1) * t66 - rSges(5,2) * t63) * t58 + t72)) - m(6) * (g(1) * (t83 * t46 + t84 * t47 + t70) + g(2) * (t83 * t44 + t84 * t45 + t71) + g(3) * (t72 + (-rSges(6,2) - pkin(7)) * t59) + (g(3) * (t83 * t63 + t84 * t66) + (g(1) * t68 + g(2) * t65) * rSges(6,2)) * t58);
U = t1;
