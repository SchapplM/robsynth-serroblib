% Calculate potential energy for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:06
% EndTime: 2019-03-09 04:18:06
% DurationCPUTime: 0.42s
% Computational Cost: add. (138->92), mult. (192->112), div. (0->0), fcn. (172->8), ass. (0->30)
t78 = rSges(6,3) + pkin(8);
t81 = pkin(2) + pkin(6);
t58 = sin(qJ(1));
t80 = g(1) * t58;
t61 = cos(qJ(1));
t79 = g(2) * t61;
t57 = sin(qJ(3));
t77 = t57 * t58;
t60 = cos(qJ(3));
t76 = t58 * t60;
t75 = t60 * t61;
t74 = rSges(7,3) + pkin(9) + pkin(8);
t73 = t61 * pkin(1) + t58 * qJ(2);
t51 = t58 * pkin(1);
t72 = t58 * pkin(7) + t51;
t71 = t60 * qJ(4);
t70 = t61 * t71 + t72;
t69 = t61 * pkin(7) + t73;
t68 = t60 * pkin(3) + t57 * qJ(4) + t81;
t67 = pkin(3) * t77 + t69;
t66 = -rSges(5,2) * t57 - rSges(5,3) * t60;
t55 = qJ(5) + qJ(6);
t46 = sin(t55);
t47 = cos(t55);
t59 = cos(qJ(5));
t65 = rSges(7,1) * t47 - rSges(7,2) * t46 + pkin(5) * t59 + pkin(4);
t56 = sin(qJ(5));
t64 = rSges(7,1) * t46 + rSges(7,2) * t47 + pkin(5) * t56;
t63 = g(1) * t67 + g(2) * t70;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t61 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t61) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t61 + rSges(3,3) * t58 + t73) + g(2) * (-rSges(3,2) * t58 + t51 + (-rSges(3,3) - qJ(2)) * t61) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t77 + rSges(4,2) * t76 + t69) + g(2) * (rSges(4,3) * t58 + t72) + g(3) * (rSges(4,1) * t60 - rSges(4,2) * t57 + t81) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t57 - rSges(4,2) * t60 - qJ(2))) * t61) - m(5) * (g(3) * (-rSges(5,2) * t60 + rSges(5,3) * t57 + t68) + (g(1) * (t66 - t71) + g(2) * rSges(5,1)) * t58 + (g(1) * rSges(5,1) + g(2) * (-t57 * pkin(3) - qJ(2) - t66)) * t61 + t63) - m(6) * (g(1) * (t61 * pkin(4) - t58 * t71 + (-t56 * t76 + t59 * t61) * rSges(6,1) + (-t56 * t61 - t59 * t76) * rSges(6,2) + t67) + g(2) * (t58 * pkin(4) - t61 * qJ(2) + (t56 * t75 + t58 * t59) * rSges(6,1) + (-t56 * t58 + t59 * t75) * rSges(6,2) + t70) + g(3) * (t78 * t60 + t68) + (g(3) * (rSges(6,1) * t56 + rSges(6,2) * t59) + t78 * t80 + (-pkin(3) - t78) * t79) * t57) - m(7) * (g(3) * t68 + (g(1) * t65 - g(2) * qJ(2)) * t61 + g(2) * t65 * t58 + (g(3) * t74 + t64 * t79 + (-qJ(4) - t64) * t80) * t60 + (g(3) * t64 + t74 * t80 + (-pkin(3) - t74) * t79) * t57 + t63);
U  = t1;
