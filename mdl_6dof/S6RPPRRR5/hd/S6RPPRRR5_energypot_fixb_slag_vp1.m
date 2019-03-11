% Calculate potential energy for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:11
% EndTime: 2019-03-09 02:28:11
% DurationCPUTime: 0.34s
% Computational Cost: add. (130->83), mult. (143->95), div. (0->0), fcn. (119->8), ass. (0->31)
t80 = rSges(7,3) + pkin(9);
t60 = -pkin(8) - pkin(7);
t79 = -t60 - qJ(2);
t78 = pkin(2) + pkin(6);
t55 = sin(qJ(4));
t77 = pkin(4) * t55;
t76 = -rSges(5,3) - pkin(7);
t53 = qJ(4) + qJ(5);
t46 = cos(t53);
t74 = rSges(6,2) * t46;
t45 = sin(t53);
t56 = sin(qJ(1));
t73 = t45 * t56;
t54 = sin(qJ(6));
t72 = t54 * t56;
t59 = cos(qJ(1));
t71 = t54 * t59;
t57 = cos(qJ(6));
t70 = t56 * t57;
t69 = t57 * t59;
t50 = t56 * pkin(1);
t68 = t56 * qJ(3) + t50;
t67 = t59 * pkin(1) + t56 * qJ(2);
t66 = pkin(3) + t78;
t65 = t56 * t77 + t68;
t64 = t59 * qJ(3) + t67;
t58 = cos(qJ(4));
t63 = t58 * pkin(4) + t66;
t62 = t56 * t60 + t59 * t77 + t64;
t61 = rSges(5,1) * t55 + rSges(5,2) * t58;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t59 - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 + rSges(2,2) * t59) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t59 + rSges(3,3) * t56 + t67) + g(2) * (-rSges(3,2) * t56 + t50 + (-rSges(3,3) - qJ(2)) * t59) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,2) * t56 + rSges(4,3) * t59 + t64) + g(2) * (rSges(4,3) * t56 + (-rSges(4,2) - qJ(2)) * t59 + t68) + g(3) * (rSges(4,1) + t78)) - m(5) * (g(1) * t64 + g(2) * t68 + g(3) * (rSges(5,1) * t58 - rSges(5,2) * t55 + t66) + (g(1) * t76 + g(2) * t61) * t56 + (g(1) * t61 + g(2) * (-qJ(2) - t76)) * t59) - m(6) * (g(1) * (-rSges(6,3) * t56 + t62) + g(2) * (rSges(6,1) * t73 + t56 * t74 + t65) + g(3) * (rSges(6,1) * t46 - rSges(6,2) * t45 + t63) + (g(1) * (rSges(6,1) * t45 + t74) + g(2) * (rSges(6,3) + t79)) * t59) - m(7) * (g(1) * (t59 * t45 * pkin(5) + (t45 * t69 - t72) * rSges(7,1) + (-t45 * t71 - t70) * rSges(7,2) + t62) + g(2) * (pkin(5) * t73 + (t45 * t70 + t71) * rSges(7,1) + (-t45 * t72 + t69) * rSges(7,2) + t65 + t79 * t59) + g(3) * (t80 * t45 + t63) + (g(3) * (rSges(7,1) * t57 - rSges(7,2) * t54 + pkin(5)) - (g(1) * t59 + g(2) * t56) * t80) * t46);
U  = t1;
