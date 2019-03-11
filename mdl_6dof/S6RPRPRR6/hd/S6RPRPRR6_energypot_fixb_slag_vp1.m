% Calculate potential energy for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:07
% EndTime: 2019-03-09 03:51:08
% DurationCPUTime: 0.52s
% Computational Cost: add. (236->104), mult. (212->128), div. (0->0), fcn. (196->12), ass. (0->43)
t76 = -pkin(8) - qJ(4);
t104 = rSges(6,3) - t76;
t103 = rSges(7,3) + pkin(9) - t76;
t102 = rSges(5,3) + qJ(4);
t78 = sin(qJ(1));
t79 = cos(qJ(1));
t101 = g(1) * t79 + g(2) * t78;
t73 = sin(pkin(10));
t98 = t73 * pkin(2) + pkin(6);
t74 = cos(pkin(11));
t59 = t74 * pkin(4) + pkin(3);
t71 = pkin(10) + qJ(3);
t63 = sin(t71);
t97 = rSges(4,2) * t63;
t65 = cos(t71);
t96 = t65 * t78;
t95 = t65 * t79;
t72 = sin(pkin(11));
t94 = t72 * t79;
t93 = t74 * t79;
t70 = pkin(11) + qJ(5);
t66 = qJ(6) + t70;
t57 = sin(t66);
t92 = t78 * t57;
t58 = cos(t66);
t91 = t78 * t58;
t62 = sin(t70);
t90 = t78 * t62;
t64 = cos(t70);
t89 = t78 * t64;
t88 = t78 * t72;
t87 = t78 * t74;
t75 = cos(pkin(10));
t60 = pkin(2) * t75 + pkin(1);
t77 = -pkin(7) - qJ(2);
t84 = t78 * t60 + t79 * t77;
t83 = rSges(3,3) + qJ(2);
t56 = t79 * t60;
t81 = -t78 * t77 + t56;
t80 = rSges(3,1) * t75 - rSges(3,2) * t73 + pkin(1);
t54 = t72 * pkin(4) + pkin(5) * t62;
t53 = pkin(5) * t64 + t59;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - t78 * rSges(2,2)) + g(2) * (t78 * rSges(2,1) + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t73 * rSges(3,1) + t75 * rSges(3,2) + pkin(6)) + (g(1) * t80 - g(2) * t83) * t79 + (g(1) * t83 + g(2) * t80) * t78) - m(4) * (g(1) * (rSges(4,1) * t95 - t79 * t97 + t56) + g(2) * (-rSges(4,3) * t79 + t84) + g(3) * (rSges(4,1) * t63 + rSges(4,2) * t65 + t98) + (g(1) * (rSges(4,3) - t77) + g(2) * (rSges(4,1) * t65 - t97)) * t78) - m(5) * (g(1) * (pkin(3) * t95 + (t65 * t93 + t88) * rSges(5,1) + (-t65 * t94 + t87) * rSges(5,2) + t81) + g(2) * (pkin(3) * t96 + (t65 * t87 - t94) * rSges(5,1) + (-t65 * t88 - t93) * rSges(5,2) + t84) + g(3) * (-t102 * t65 + t98) + (g(3) * (rSges(5,1) * t74 - rSges(5,2) * t72 + pkin(3)) + t101 * t102) * t63) - m(6) * (g(1) * (t59 * t95 + pkin(4) * t88 + (t64 * t95 + t90) * rSges(6,1) + (-t62 * t95 + t89) * rSges(6,2) + t81) + g(2) * (t59 * t96 - pkin(4) * t94 + (-t62 * t79 + t65 * t89) * rSges(6,1) + (-t64 * t79 - t65 * t90) * rSges(6,2) + t84) + g(3) * (-t104 * t65 + t98) + (g(3) * (rSges(6,1) * t64 - rSges(6,2) * t62 + t59) + t101 * t104) * t63) - m(7) * (g(1) * (t53 * t95 + t78 * t54 + (t58 * t95 + t92) * rSges(7,1) + (-t57 * t95 + t91) * rSges(7,2) + t81) + g(2) * (t53 * t96 - t79 * t54 + (-t57 * t79 + t65 * t91) * rSges(7,1) + (-t58 * t79 - t65 * t92) * rSges(7,2) + t84) + g(3) * (-t103 * t65 + t98) + (g(3) * (rSges(7,1) * t58 - rSges(7,2) * t57 + t53) + t101 * t103) * t63);
U  = t1;
