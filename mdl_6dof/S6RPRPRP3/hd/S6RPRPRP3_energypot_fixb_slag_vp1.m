% Calculate potential energy for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:43
% EndTime: 2019-03-09 03:07:44
% DurationCPUTime: 0.46s
% Computational Cost: add. (250->93), mult. (213->113), div. (0->0), fcn. (201->10), ass. (0->39)
t107 = rSges(7,1) + pkin(5);
t106 = rSges(5,3) + qJ(4);
t105 = rSges(7,3) + qJ(6);
t77 = qJ(1) + pkin(9);
t71 = sin(t77);
t73 = cos(t77);
t104 = g(1) * t73 + g(2) * t71;
t81 = sin(qJ(3));
t101 = rSges(4,2) * t81;
t78 = sin(pkin(10));
t100 = t71 * t78;
t83 = cos(qJ(3));
t99 = t71 * t83;
t98 = t73 * t78;
t97 = t73 * t83;
t96 = t78 * t83;
t79 = cos(pkin(10));
t95 = t79 * t83;
t92 = pkin(6) + qJ(2);
t82 = sin(qJ(1));
t74 = t82 * pkin(1);
t91 = t71 * pkin(2) + t74;
t84 = cos(qJ(1));
t75 = t84 * pkin(1);
t89 = t73 * pkin(2) + t71 * pkin(7) + t75;
t68 = pkin(4) * t79 + pkin(3);
t80 = -pkin(8) - qJ(4);
t88 = t81 * t68 + t83 * t80 + t92;
t87 = -t73 * pkin(7) + t91;
t86 = pkin(4) * t100 + t68 * t97 + t89;
t85 = -pkin(4) * t98 + t68 * t99 + t87;
t76 = pkin(10) + qJ(5);
t72 = cos(t76);
t70 = sin(t76);
t59 = t70 * t71 + t72 * t97;
t58 = t70 * t97 - t71 * t72;
t57 = -t70 * t73 + t72 * t99;
t56 = t70 * t99 + t72 * t73;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t84 - t82 * rSges(2,2)) + g(2) * (t82 * rSges(2,1) + rSges(2,2) * t84) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t73 - rSges(3,2) * t71 + t75) + g(2) * (rSges(3,1) * t71 + rSges(3,2) * t73 + t74) + g(3) * (rSges(3,3) + t92)) - m(4) * (g(1) * (t71 * rSges(4,3) + t89) + g(2) * (rSges(4,1) * t99 - t71 * t101 + t91) + g(3) * (rSges(4,1) * t81 + rSges(4,2) * t83 + t92) + (g(1) * (rSges(4,1) * t83 - t101) + g(2) * (-rSges(4,3) - pkin(7))) * t73) - m(5) * (g(1) * (pkin(3) * t97 + (t73 * t95 + t100) * rSges(5,1) + (t71 * t79 - t73 * t96) * rSges(5,2) + t89) + g(2) * (pkin(3) * t99 + (t71 * t95 - t98) * rSges(5,1) + (-t71 * t96 - t73 * t79) * rSges(5,2) + t87) + g(3) * (-t106 * t83 + t92) + (g(3) * (rSges(5,1) * t79 - rSges(5,2) * t78 + pkin(3)) + t104 * t106) * t81) - m(6) * (g(1) * (t59 * rSges(6,1) - t58 * rSges(6,2) + t86) + g(2) * (t57 * rSges(6,1) - t56 * rSges(6,2) + t85) + g(3) * (-t83 * rSges(6,3) + t88) + (g(3) * (rSges(6,1) * t72 - rSges(6,2) * t70) + t104 * (rSges(6,3) - t80)) * t81) - m(7) * (g(1) * (t105 * t58 + t107 * t59 + t86) + g(2) * (t105 * t56 + t107 * t57 + t85) + g(3) * (-t83 * rSges(7,2) + t88) + (g(3) * (t105 * t70 + t107 * t72) + t104 * (rSges(7,2) - t80)) * t81);
U  = t1;
