% Calculate potential energy for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:47:51
% EndTime: 2019-03-09 03:47:52
% DurationCPUTime: 0.44s
% Computational Cost: add. (231->94), mult. (247->114), div. (0->0), fcn. (244->10), ass. (0->39)
t104 = pkin(9) + rSges(7,3);
t80 = sin(pkin(10));
t106 = t80 * pkin(2) + pkin(6);
t82 = -pkin(7) - qJ(2);
t105 = -pkin(8) - t82;
t79 = pkin(10) + qJ(3);
t75 = sin(t79);
t87 = cos(qJ(5));
t103 = t75 * t87;
t88 = cos(qJ(1));
t102 = t75 * t88;
t76 = cos(t79);
t85 = sin(qJ(1));
t101 = t76 * t85;
t100 = t76 * t88;
t81 = cos(pkin(10));
t73 = pkin(2) * t81 + pkin(1);
t99 = t85 * t73 + t88 * t82;
t98 = qJ(4) * t75;
t97 = rSges(3,3) + qJ(2);
t96 = t75 * pkin(3) + t106;
t66 = t88 * t73;
t95 = pkin(3) * t100 + t88 * t98 + t66;
t94 = pkin(3) * t101 + t85 * t98 + t99;
t93 = pkin(4) * t100 + t95;
t84 = sin(qJ(5));
t59 = t75 * t84 + t76 * t87;
t92 = pkin(4) * t101 + t88 * pkin(8) + t94;
t91 = rSges(3,1) * t81 - rSges(3,2) * t80 + pkin(1);
t83 = sin(qJ(6));
t86 = cos(qJ(6));
t90 = rSges(7,1) * t86 - rSges(7,2) * t83 + pkin(5);
t89 = t75 * pkin(4) - qJ(4) * t76 + t96;
t60 = -t76 * t84 + t103;
t58 = t59 * t88;
t57 = t84 * t100 - t87 * t102;
t56 = t59 * t85;
t55 = t84 * t101 - t85 * t103;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t88 - t85 * rSges(2,2)) + g(2) * (t85 * rSges(2,1) + rSges(2,2) * t88) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t80 + t81 * rSges(3,2) + pkin(6)) + (g(1) * t91 - g(2) * t97) * t88 + (g(1) * t97 + g(2) * t91) * t85) - m(4) * (g(1) * (rSges(4,1) * t100 - rSges(4,2) * t102 + t66) + g(2) * (-rSges(4,3) * t88 + t99) + g(3) * (rSges(4,1) * t75 + rSges(4,2) * t76 + t106) + (g(1) * (rSges(4,3) - t82) + g(2) * (rSges(4,1) * t76 - rSges(4,2) * t75)) * t85) - m(5) * (g(1) * (rSges(5,1) * t100 + rSges(5,3) * t102 + t95) + g(2) * (-rSges(5,2) * t88 + t94) + g(3) * (rSges(5,1) * t75 + (-rSges(5,3) - qJ(4)) * t76 + t96) + (g(1) * (rSges(5,2) - t82) + g(2) * (rSges(5,1) * t76 + rSges(5,3) * t75)) * t85) - m(6) * (g(2) * (t56 * rSges(6,1) - t55 * rSges(6,2) + rSges(6,3) * t88 + t92) + g(3) * (rSges(6,1) * t60 - rSges(6,2) * t59 + t89) + (rSges(6,1) * t58 - rSges(6,2) * t57 + t93 + (-rSges(6,3) + t105) * t85) * g(1)) - m(7) * (g(1) * (t90 * t58 + (-rSges(7,1) * t83 - rSges(7,2) * t86 + t105) * t85 + t104 * t57 + t93) + g(2) * (t56 * pkin(5) + (t56 * t86 + t83 * t88) * rSges(7,1) + (-t56 * t83 + t86 * t88) * rSges(7,2) + t104 * t55 + t92) + (t104 * t59 + t90 * t60 + t89) * g(3));
U  = t1;
