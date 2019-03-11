% Calculate potential energy for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:38:56
% EndTime: 2019-03-09 01:38:56
% DurationCPUTime: 0.41s
% Computational Cost: add. (234->93), mult. (172->111), div. (0->0), fcn. (152->12), ass. (0->40)
t97 = rSges(7,3) + pkin(8) + qJ(5);
t66 = qJ(1) + pkin(9);
t57 = sin(t66);
t60 = cos(t66);
t96 = g(1) * t60 + g(2) * t57;
t95 = rSges(6,3) + qJ(5);
t65 = pkin(10) + qJ(4);
t56 = sin(t65);
t92 = rSges(5,2) * t56;
t59 = cos(t65);
t91 = t57 * t59;
t67 = sin(pkin(11));
t90 = t57 * t67;
t69 = cos(pkin(11));
t89 = t57 * t69;
t88 = t59 * t60;
t64 = pkin(11) + qJ(6);
t55 = sin(t64);
t87 = t60 * t55;
t58 = cos(t64);
t86 = t60 * t58;
t85 = t60 * t67;
t84 = t60 * t69;
t82 = pkin(6) + qJ(2);
t70 = cos(pkin(10));
t54 = pkin(3) * t70 + pkin(2);
t74 = cos(qJ(1));
t63 = t74 * pkin(1);
t81 = t54 * t60 + t63;
t80 = rSges(4,3) + qJ(3);
t68 = sin(pkin(10));
t78 = pkin(3) * t68 + t82;
t73 = sin(qJ(1));
t62 = t73 * pkin(1);
t72 = -pkin(7) - qJ(3);
t77 = t54 * t57 + t60 * t72 + t62;
t76 = -t57 * t72 + t81;
t75 = rSges(4,1) * t70 - rSges(4,2) * t68 + pkin(2);
t53 = pkin(5) * t69 + pkin(4);
t1 = -m(1) * (rSges(1,1) * g(1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t74 - rSges(2,2) * t73) + g(2) * (rSges(2,1) * t73 + rSges(2,2) * t74) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t60 - rSges(3,2) * t57 + t63) + g(2) * (rSges(3,1) * t57 + rSges(3,2) * t60 + t62) + g(3) * (rSges(3,3) + t82)) - m(4) * (g(1) * t63 + g(2) * t62 + g(3) * (rSges(4,1) * t68 + rSges(4,2) * t70 + t82) + (g(1) * t75 - g(2) * t80) * t60 + (g(1) * t80 + g(2) * t75) * t57) - m(5) * (g(1) * (rSges(5,1) * t88 - t60 * t92 + t81) + g(2) * (-t60 * rSges(5,3) + t77) + g(3) * (rSges(5,1) * t56 + rSges(5,2) * t59 + t78) + (g(1) * (rSges(5,3) - t72) + g(2) * (rSges(5,1) * t59 - t92)) * t57) - m(6) * (g(1) * (pkin(4) * t88 + (t59 * t84 + t90) * rSges(6,1) + (-t59 * t85 + t89) * rSges(6,2) + t76) + g(2) * (pkin(4) * t91 + (t59 * t89 - t85) * rSges(6,1) + (-t59 * t90 - t84) * rSges(6,2) + t77) + g(3) * (-t95 * t59 + t78) + (g(3) * (rSges(6,1) * t69 - rSges(6,2) * t67 + pkin(4)) + t96 * t95) * t56) - m(7) * (g(1) * (t53 * t88 + pkin(5) * t90 + (t55 * t57 + t59 * t86) * rSges(7,1) + (t57 * t58 - t59 * t87) * rSges(7,2) + t76) + g(2) * (t53 * t91 - pkin(5) * t85 + (t58 * t91 - t87) * rSges(7,1) + (-t55 * t91 - t86) * rSges(7,2) + t77) + g(3) * (-t97 * t59 + t78) + (g(3) * (rSges(7,1) * t58 - rSges(7,2) * t55 + t53) + t96 * t97) * t56);
U  = t1;
