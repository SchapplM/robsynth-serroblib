% Calculate potential energy for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR11_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:20
% EndTime: 2019-12-31 21:32:21
% DurationCPUTime: 0.42s
% Computational Cost: add. (131->87), mult. (242->113), div. (0->0), fcn. (250->8), ass. (0->28)
t80 = -rSges(6,3) - pkin(8);
t61 = sin(qJ(2));
t79 = t61 * pkin(2) + pkin(5);
t62 = sin(qJ(1));
t78 = t62 * t61;
t65 = cos(qJ(2));
t77 = t62 * t65;
t60 = sin(qJ(3));
t66 = cos(qJ(1));
t76 = t66 * t60;
t75 = t66 * t61;
t64 = cos(qJ(3));
t74 = t66 * t64;
t73 = t66 * pkin(1) + t62 * pkin(6);
t72 = rSges(5,3) + qJ(4);
t71 = t79 + (pkin(3) * t64 + qJ(4) * t60) * t61;
t70 = t66 * t65 * pkin(2) + pkin(7) * t75 + t73;
t47 = t62 * t60 + t65 * t74;
t69 = t47 * pkin(3) + t70;
t57 = t62 * pkin(1);
t68 = pkin(2) * t77 - t66 * pkin(6) + pkin(7) * t78 + t57;
t45 = t64 * t77 - t76;
t67 = t45 * pkin(3) + t68;
t63 = cos(qJ(5));
t59 = sin(qJ(5));
t46 = -t62 * t64 + t65 * t76;
t44 = t60 * t77 + t74;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t66 * rSges(2,1) - t62 * rSges(2,2)) + g(2) * (t62 * rSges(2,1) + t66 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t62 * rSges(3,3) + t73) + g(2) * (rSges(3,1) * t77 - rSges(3,2) * t78 + t57) + g(3) * (t61 * rSges(3,1) + t65 * rSges(3,2) + pkin(5)) + (g(1) * (rSges(3,1) * t65 - rSges(3,2) * t61) + g(2) * (-rSges(3,3) - pkin(6))) * t66) - m(4) * (g(1) * (t47 * rSges(4,1) - t46 * rSges(4,2) + rSges(4,3) * t75 + t70) + g(2) * (t45 * rSges(4,1) - t44 * rSges(4,2) + rSges(4,3) * t78 + t68) + g(3) * ((-rSges(4,3) - pkin(7)) * t65 + (rSges(4,1) * t64 - rSges(4,2) * t60) * t61 + t79)) - m(5) * (g(1) * (t47 * rSges(5,1) + rSges(5,2) * t75 + t46 * t72 + t69) + g(2) * (t45 * rSges(5,1) + rSges(5,2) * t78 + t44 * t72 + t67) + g(3) * ((-rSges(5,2) - pkin(7)) * t65 + (rSges(5,1) * t64 + rSges(5,3) * t60) * t61 + t71)) - m(6) * (g(1) * (t47 * pkin(4) + t46 * qJ(4) + (t46 * t59 + t47 * t63) * rSges(6,1) + (t46 * t63 - t47 * t59) * rSges(6,2) + t69) + g(2) * (t45 * pkin(4) + t44 * qJ(4) + (t44 * t59 + t45 * t63) * rSges(6,1) + (t44 * t63 - t45 * t59) * rSges(6,2) + t67) + (g(1) * t66 + g(2) * t62) * t61 * t80 + (t71 + (-pkin(7) - t80) * t65 + (t64 * pkin(4) + (t59 * t60 + t63 * t64) * rSges(6,1) + (-t59 * t64 + t60 * t63) * rSges(6,2)) * t61) * g(3));
U = t1;
