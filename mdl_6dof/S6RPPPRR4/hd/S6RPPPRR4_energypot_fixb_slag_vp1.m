% Calculate potential energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:14
% EndTime: 2019-03-09 01:35:14
% DurationCPUTime: 0.41s
% Computational Cost: add. (149->72), mult. (223->81), div. (0->0), fcn. (243->8), ass. (0->26)
t55 = sin(qJ(6));
t57 = cos(qJ(6));
t83 = rSges(7,1) * t57 - rSges(7,2) * t55 + pkin(5);
t80 = rSges(7,3) + pkin(8);
t79 = -t55 * rSges(7,1) - t57 * rSges(7,2) - pkin(7);
t56 = sin(qJ(5));
t78 = t83 * t56 + qJ(4);
t76 = -rSges(6,3) - pkin(7);
t74 = cos(qJ(1));
t73 = sin(qJ(1));
t70 = pkin(6) - qJ(3);
t69 = t74 * pkin(1) + t73 * qJ(2);
t68 = rSges(5,3) + qJ(4);
t67 = cos(pkin(9));
t66 = sin(pkin(9));
t65 = -pkin(4) + t70;
t64 = t74 * pkin(2) + t69;
t45 = -t66 * t73 - t67 * t74;
t63 = -t45 * pkin(3) + t64;
t62 = t73 * pkin(1) - qJ(2) * t74;
t61 = t73 * pkin(2) + t62;
t58 = cos(qJ(5));
t60 = rSges(6,1) * t56 + rSges(6,2) * t58 + qJ(4);
t46 = t66 * t74 - t67 * t73;
t59 = -t46 * pkin(3) + t61;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t74 - rSges(2,2) * t73) + g(2) * (rSges(2,1) * t73 + rSges(2,2) * t74) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t74 + rSges(3,3) * t73 + t69) + g(2) * (rSges(3,1) * t73 - rSges(3,3) * t74 + t62) + g(3) * (pkin(6) + rSges(3,2))) - m(4) * (g(1) * (-t45 * rSges(4,1) - t46 * rSges(4,2) + t64) + g(2) * (-t46 * rSges(4,1) + t45 * rSges(4,2) + t61) + g(3) * (-rSges(4,3) + t70)) - m(5) * (g(1) * (t45 * rSges(5,2) + t46 * t68 + t63) + g(2) * (t46 * rSges(5,2) - t45 * t68 + t59) + g(3) * (-rSges(5,1) + t70)) - m(6) * (g(1) * t63 + g(2) * t59 + g(3) * (-t58 * rSges(6,1) + t56 * rSges(6,2) + t65) + (g(1) * t60 + g(2) * t76) * t46 + (g(1) * t76 - g(2) * t60) * t45) - m(7) * (g(1) * (t79 * t45 + t78 * t46 + t63) + g(2) * (-t78 * t45 + t79 * t46 + t59) + g(3) * (-t80 * t56 + t65) + (-g(3) * t83 + (-g(1) * t46 + g(2) * t45) * t80) * t58);
U  = t1;
