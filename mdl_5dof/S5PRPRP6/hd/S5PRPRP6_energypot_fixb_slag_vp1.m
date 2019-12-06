% Calculate potential energy for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:16
% EndTime: 2019-12-05 15:40:16
% DurationCPUTime: 0.34s
% Computational Cost: add. (108->76), mult. (182->97), div. (0->0), fcn. (170->6), ass. (0->28)
t80 = rSges(6,1) + pkin(4);
t79 = rSges(6,3) + qJ(5);
t59 = sin(pkin(7));
t62 = sin(qJ(2));
t78 = t59 * t62;
t64 = cos(qJ(2));
t77 = t59 * t64;
t60 = cos(pkin(7));
t76 = t60 * t64;
t61 = sin(qJ(4));
t75 = t61 * t62;
t63 = cos(qJ(4));
t74 = t62 * t63;
t73 = t60 * pkin(1) + t59 * pkin(5);
t72 = qJ(3) * t62;
t71 = t62 * pkin(2) + qJ(1);
t55 = t59 * pkin(1);
t70 = pkin(2) * t77 + t59 * t72 + t55;
t69 = t62 * pkin(6) + t71;
t68 = pkin(2) * t76 + t60 * t72 + t73;
t67 = g(1) * t60 + g(2) * t59;
t66 = t59 * pkin(3) + pkin(6) * t76 + t68;
t65 = pkin(6) * t77 + t70 + (-pkin(3) - pkin(5)) * t60;
t45 = t59 * t75 - t60 * t63;
t44 = t59 * t74 + t60 * t61;
t43 = t59 * t63 + t60 * t75;
t42 = t59 * t61 - t60 * t74;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t60 - rSges(2,2) * t59) + g(2) * (rSges(2,1) * t59 + rSges(2,2) * t60) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t59 * rSges(3,3) + t73) + g(2) * (rSges(3,1) * t77 - rSges(3,2) * t78 + t55) + g(3) * (rSges(3,1) * t62 + rSges(3,2) * t64 + qJ(1)) + (g(1) * (rSges(3,1) * t64 - rSges(3,2) * t62) + g(2) * (-rSges(3,3) - pkin(5))) * t60) - m(4) * (g(1) * (t59 * rSges(4,1) + t68) + g(2) * (-rSges(4,2) * t77 + rSges(4,3) * t78 + t70) + g(3) * (-t62 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t64 + t71) + (g(1) * (-rSges(4,2) * t64 + rSges(4,3) * t62) + g(2) * (-rSges(4,1) - pkin(5))) * t60) - m(5) * (g(1) * (rSges(5,1) * t43 - rSges(5,2) * t42 + t66) + g(2) * (t45 * rSges(5,1) + t44 * rSges(5,2) + t65) + g(3) * (t62 * rSges(5,3) + t69) + (g(3) * (-rSges(5,1) * t61 - rSges(5,2) * t63 - qJ(3)) + t67 * rSges(5,3)) * t64) - m(6) * (g(1) * (t42 * t79 + t43 * t80 + t66) + g(2) * (-t44 * t79 + t45 * t80 + t65) + g(3) * (t62 * rSges(6,2) + t69) + (g(3) * (-t61 * t80 + t63 * t79 - qJ(3)) + t67 * rSges(6,2)) * t64);
U = t1;
