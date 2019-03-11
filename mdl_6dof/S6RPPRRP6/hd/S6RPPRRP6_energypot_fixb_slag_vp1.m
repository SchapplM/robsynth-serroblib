% Calculate potential energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:07
% EndTime: 2019-03-09 02:10:08
% DurationCPUTime: 0.32s
% Computational Cost: add. (120->81), mult. (176->93), div. (0->0), fcn. (160->6), ass. (0->30)
t88 = rSges(7,1) + pkin(5);
t63 = sin(qJ(1));
t66 = cos(qJ(1));
t87 = g(1) * t66 + g(2) * t63;
t86 = rSges(7,3) + qJ(6);
t85 = pkin(2) + pkin(6);
t62 = sin(qJ(4));
t84 = pkin(4) * t62;
t61 = sin(qJ(5));
t79 = t63 * t61;
t64 = cos(qJ(5));
t78 = t63 * t64;
t77 = t66 * t61;
t76 = t66 * t64;
t57 = t63 * pkin(1);
t75 = t63 * qJ(3) + t57;
t74 = t66 * pkin(1) + t63 * qJ(2);
t73 = pkin(3) + t85;
t72 = t66 * pkin(7) + t75;
t71 = t66 * qJ(3) + t74;
t65 = cos(qJ(4));
t70 = t65 * pkin(4) + t62 * pkin(8) + t73;
t69 = rSges(5,1) * t62 + rSges(5,2) * t65;
t68 = -t63 * pkin(7) + t66 * t84 + t71;
t67 = -t66 * qJ(2) + t63 * t84 + t72;
t49 = t62 * t76 - t79;
t48 = t62 * t77 + t78;
t47 = t62 * t78 + t77;
t46 = t62 * t79 - t76;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t66 * rSges(2,1) - t63 * rSges(2,2)) + g(2) * (t63 * rSges(2,1) + t66 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-t66 * rSges(3,2) + t63 * rSges(3,3) + t74) + g(2) * (-t63 * rSges(3,2) + t57 + (-rSges(3,3) - qJ(2)) * t66) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (t63 * rSges(4,2) + t66 * rSges(4,3) + t71) + g(2) * (t63 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t66 + t75) + g(3) * (rSges(4,1) + t85)) - m(5) * (g(1) * t71 + g(2) * t72 + g(3) * (t65 * rSges(5,1) - t62 * rSges(5,2) + t73) + (g(1) * t69 + g(2) * (rSges(5,3) - qJ(2))) * t66 + (g(1) * (-rSges(5,3) - pkin(7)) + g(2) * t69) * t63) - m(6) * (g(1) * (t49 * rSges(6,1) - t48 * rSges(6,2) + t68) + g(2) * (t47 * rSges(6,1) - t46 * rSges(6,2) + t67) + g(3) * (t62 * rSges(6,3) + t70) + (g(3) * (rSges(6,1) * t64 - rSges(6,2) * t61) + t87 * (-rSges(6,3) - pkin(8))) * t65) - m(7) * (g(1) * (t86 * t48 + t88 * t49 + t68) + g(2) * (t86 * t46 + t88 * t47 + t67) + g(3) * (t62 * rSges(7,2) + t70) + (g(3) * (t86 * t61 + t88 * t64) + t87 * (-rSges(7,2) - pkin(8))) * t65);
U  = t1;
