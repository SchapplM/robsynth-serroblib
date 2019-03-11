% Calculate potential energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:15
% EndTime: 2019-03-09 08:21:16
% DurationCPUTime: 0.51s
% Computational Cost: add. (181->103), mult. (335->124), div. (0->0), fcn. (351->8), ass. (0->40)
t81 = sin(qJ(6));
t84 = cos(qJ(6));
t112 = t81 * rSges(7,1) + t84 * rSges(7,2) + qJ(5);
t111 = t84 * rSges(7,1) - t81 * rSges(7,2) + pkin(5) + qJ(4);
t85 = cos(qJ(2));
t110 = g(3) * t85;
t109 = rSges(7,3) + pkin(8);
t82 = sin(qJ(2));
t108 = t82 * pkin(2) + pkin(6);
t79 = sin(pkin(9));
t107 = t79 * t82;
t80 = cos(pkin(9));
t106 = t80 * t82;
t83 = sin(qJ(1));
t105 = t82 * t83;
t86 = cos(qJ(1));
t104 = t82 * t86;
t103 = t83 * t85;
t102 = t85 * t86;
t101 = t86 * t79;
t100 = -pkin(4) - qJ(3);
t99 = t86 * pkin(1) + t83 * pkin(7);
t98 = qJ(3) * t82;
t97 = rSges(6,1) + qJ(4);
t96 = rSges(5,3) + qJ(4);
t95 = rSges(6,3) + qJ(5);
t94 = pkin(3) * t106 + qJ(4) * t107 + t108;
t93 = pkin(2) * t102 + t86 * t98 + t99;
t92 = qJ(5) * t106 + t94;
t63 = t80 * t102 + t83 * t79;
t91 = t63 * pkin(3) + t93;
t77 = t83 * pkin(1);
t90 = pkin(2) * t103 - t86 * pkin(7) + t83 * t98 + t77;
t89 = pkin(4) * t104 + t91;
t61 = t80 * t103 - t101;
t88 = t61 * pkin(3) + t90;
t87 = pkin(4) * t105 + t88;
t62 = t85 * t101 - t83 * t80;
t60 = t79 * t103 + t80 * t86;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t86 - t83 * rSges(2,2)) + g(2) * (t83 * rSges(2,1) + rSges(2,2) * t86) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t83 * rSges(3,3) + t99) + g(2) * (rSges(3,1) * t103 - rSges(3,2) * t105 + t77) + g(3) * (rSges(3,1) * t82 + rSges(3,2) * t85 + pkin(6)) + (g(1) * (rSges(3,1) * t85 - rSges(3,2) * t82) + g(2) * (-rSges(3,3) - pkin(7))) * t86) - m(4) * (g(1) * (t63 * rSges(4,1) - t62 * rSges(4,2) + rSges(4,3) * t104 + t93) + g(2) * (t61 * rSges(4,1) - t60 * rSges(4,2) + rSges(4,3) * t105 + t90) + g(3) * ((-rSges(4,3) - qJ(3)) * t85 + (rSges(4,1) * t80 - rSges(4,2) * t79) * t82 + t108)) - m(5) * (g(1) * (rSges(5,1) * t104 - t63 * rSges(5,2) + t96 * t62 + t91) + g(2) * (rSges(5,1) * t105 - t61 * rSges(5,2) + t96 * t60 + t88) + g(3) * ((-rSges(5,1) - qJ(3)) * t85 + (-rSges(5,2) * t80 + rSges(5,3) * t79) * t82 + t94)) - m(6) * (g(1) * (-rSges(6,2) * t104 + t97 * t62 + t95 * t63 + t89) + g(2) * (-rSges(6,2) * t105 + t97 * t60 + t95 * t61 + t87) + g(3) * (rSges(6,1) * t107 + rSges(6,3) * t106 + t92) + (rSges(6,2) + t100) * t110) - m(7) * (g(1) * (t111 * t62 + t112 * t63 + t89) + g(2) * (t111 * t60 + t112 * t61 + t87) + g(3) * t92 + (t100 - t109) * t110 + (g(3) * (t79 * pkin(5) + (t79 * t84 + t80 * t81) * rSges(7,1) + (-t79 * t81 + t80 * t84) * rSges(7,2)) + (g(1) * t86 + g(2) * t83) * t109) * t82);
U  = t1;
