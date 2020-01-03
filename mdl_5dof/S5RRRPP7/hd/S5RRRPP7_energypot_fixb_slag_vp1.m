% Calculate potential energy for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:53
% EndTime: 2019-12-31 21:03:53
% DurationCPUTime: 0.36s
% Computational Cost: add. (121->78), mult. (216->96), div. (0->0), fcn. (214->6), ass. (0->27)
t79 = rSges(6,1) + pkin(4);
t78 = rSges(6,2) + qJ(4);
t60 = sin(qJ(2));
t77 = t60 * pkin(2) + pkin(5);
t61 = sin(qJ(1));
t76 = t60 * t61;
t64 = cos(qJ(1));
t75 = t60 * t64;
t63 = cos(qJ(2));
t74 = t61 * t63;
t73 = t63 * t64;
t72 = t64 * pkin(1) + t61 * pkin(6);
t71 = rSges(5,3) + qJ(4);
t70 = -rSges(6,3) - qJ(5);
t59 = sin(qJ(3));
t62 = cos(qJ(3));
t69 = t77 + (pkin(3) * t62 + qJ(4) * t59) * t60;
t68 = pkin(2) * t73 + pkin(7) * t75 + t72;
t47 = t61 * t59 + t62 * t73;
t67 = t47 * pkin(3) + t68;
t57 = t61 * pkin(1);
t66 = pkin(2) * t74 - pkin(6) * t64 + pkin(7) * t76 + t57;
t45 = -t59 * t64 + t62 * t74;
t65 = t45 * pkin(3) + t66;
t46 = t59 * t73 - t61 * t62;
t44 = t59 * t74 + t62 * t64;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t64 - t61 * rSges(2,2)) + g(2) * (t61 * rSges(2,1) + rSges(2,2) * t64) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t61 * rSges(3,3) + t72) + g(2) * (rSges(3,1) * t74 - rSges(3,2) * t76 + t57) + g(3) * (rSges(3,1) * t60 + rSges(3,2) * t63 + pkin(5)) + (g(1) * (rSges(3,1) * t63 - rSges(3,2) * t60) + g(2) * (-rSges(3,3) - pkin(6))) * t64) - m(4) * (g(1) * (t47 * rSges(4,1) - t46 * rSges(4,2) + rSges(4,3) * t75 + t68) + g(2) * (t45 * rSges(4,1) - t44 * rSges(4,2) + rSges(4,3) * t76 + t66) + g(3) * ((-rSges(4,3) - pkin(7)) * t63 + (rSges(4,1) * t62 - rSges(4,2) * t59) * t60 + t77)) - m(5) * (g(1) * (t47 * rSges(5,1) + rSges(5,2) * t75 + t71 * t46 + t67) + g(2) * (t45 * rSges(5,1) + rSges(5,2) * t76 + t71 * t44 + t65) + g(3) * ((-rSges(5,2) - pkin(7)) * t63 + (rSges(5,1) * t62 + rSges(5,3) * t59) * t60 + t69)) - m(6) * (g(1) * (t78 * t46 + t47 * t79 + t67) + g(2) * (t78 * t44 + t45 * t79 + t65) + (g(1) * t64 + g(2) * t61) * t60 * t70 + (t69 + (-pkin(7) - t70) * t63 + (rSges(6,2) * t59 + t62 * t79) * t60) * g(3));
U = t1;
