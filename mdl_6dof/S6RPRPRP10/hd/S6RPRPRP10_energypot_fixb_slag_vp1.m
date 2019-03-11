% Calculate potential energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:52
% EndTime: 2019-03-09 03:30:52
% DurationCPUTime: 0.36s
% Computational Cost: add. (130->89), mult. (202->108), div. (0->0), fcn. (186->6), ass. (0->32)
t91 = rSges(7,1) + pkin(5);
t90 = rSges(7,3) + qJ(6);
t89 = pkin(2) + pkin(6);
t88 = -pkin(3) - pkin(8);
t68 = sin(qJ(1));
t87 = g(1) * t68;
t71 = cos(qJ(1));
t86 = g(2) * t71;
t67 = sin(qJ(3));
t85 = t67 * t68;
t70 = cos(qJ(3));
t84 = t68 * t70;
t83 = t70 * t71;
t82 = t71 * pkin(1) + t68 * qJ(2);
t60 = t68 * pkin(1);
t81 = t68 * pkin(7) + t60;
t80 = qJ(4) * t70;
t79 = t71 * t80 + t81;
t78 = t71 * pkin(7) + t82;
t77 = t70 * pkin(3) + t67 * qJ(4) + t89;
t76 = pkin(3) * t85 + t78;
t75 = t70 * pkin(8) + t77;
t74 = -rSges(5,2) * t67 - rSges(5,3) * t70;
t73 = t68 * pkin(4) - qJ(2) * t71 + t79;
t72 = t71 * pkin(4) + pkin(8) * t85 - t68 * t80 + t76;
t69 = cos(qJ(5));
t66 = sin(qJ(5));
t51 = -t66 * t84 + t69 * t71;
t50 = t66 * t71 + t69 * t84;
t49 = t66 * t83 + t68 * t69;
t48 = t66 * t68 - t69 * t83;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t71 - t68 * rSges(2,2)) + g(2) * (t68 * rSges(2,1) + rSges(2,2) * t71) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t71 + t68 * rSges(3,3) + t82) + g(2) * (-t68 * rSges(3,2) + t60 + (-rSges(3,3) - qJ(2)) * t71) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t85 + rSges(4,2) * t84 + t78) + g(2) * (t68 * rSges(4,3) + t81) + g(3) * (rSges(4,1) * t70 - rSges(4,2) * t67 + t89) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t67 - rSges(4,2) * t70 - qJ(2))) * t71) - m(5) * (g(1) * t76 + g(2) * t79 + g(3) * (-rSges(5,2) * t70 + rSges(5,3) * t67 + t77) + (g(1) * (t74 - t80) + g(2) * rSges(5,1)) * t68 + (g(1) * rSges(5,1) + g(2) * (-pkin(3) * t67 - qJ(2) - t74)) * t71) - m(6) * (g(1) * (rSges(6,1) * t51 - rSges(6,2) * t50 + t72) + g(2) * (t49 * rSges(6,1) - t48 * rSges(6,2) + t73) + g(3) * (rSges(6,3) * t70 + t75) + (rSges(6,3) * t87 + g(3) * (rSges(6,1) * t66 + rSges(6,2) * t69) + (-rSges(6,3) + t88) * t86) * t67) - m(7) * (g(1) * (t90 * t50 + t91 * t51 + t72) + g(2) * (t90 * t48 + t91 * t49 + t73) + g(3) * (rSges(7,2) * t70 + t75) + (rSges(7,2) * t87 + g(3) * (t91 * t66 - t90 * t69) + (-rSges(7,2) + t88) * t86) * t67);
U  = t1;
