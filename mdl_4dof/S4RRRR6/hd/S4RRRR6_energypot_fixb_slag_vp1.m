% Calculate potential energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR6_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:18
% EndTime: 2019-12-31 17:29:19
% DurationCPUTime: 0.20s
% Computational Cost: add. (130->75), mult. (272->106), div. (0->0), fcn. (311->10), ass. (0->38)
t94 = rSges(4,3) + pkin(7);
t70 = cos(pkin(4));
t93 = t70 * pkin(6) + pkin(5);
t92 = pkin(8) + rSges(5,3);
t69 = sin(pkin(4));
t73 = sin(qJ(2));
t91 = t69 * t73;
t74 = sin(qJ(1));
t90 = t69 * t74;
t76 = cos(qJ(3));
t89 = t69 * t76;
t77 = cos(qJ(2));
t88 = t69 * t77;
t78 = cos(qJ(1));
t87 = t69 * t78;
t86 = t74 * t73;
t85 = t74 * t77;
t84 = t78 * t73;
t83 = t78 * t77;
t82 = t78 * pkin(1) + pkin(6) * t90;
t81 = pkin(2) * t91 + t93;
t60 = -t70 * t86 + t83;
t80 = t60 * pkin(2) + t82;
t58 = t70 * t84 + t85;
t67 = t74 * pkin(1);
t79 = t58 * pkin(2) - pkin(6) * t87 + t67;
t75 = cos(qJ(4));
t72 = sin(qJ(3));
t71 = sin(qJ(4));
t59 = t70 * t85 + t84;
t57 = -t70 * t83 + t86;
t56 = t70 * t72 + t73 * t89;
t55 = -t70 * t76 + t72 * t91;
t52 = t60 * t76 + t72 * t90;
t51 = t60 * t72 - t74 * t89;
t50 = t58 * t76 - t72 * t87;
t49 = t58 * t72 + t76 * t87;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t78 * rSges(2,1) - t74 * rSges(2,2)) + g(2) * (t74 * rSges(2,1) + t78 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t60 * rSges(3,1) - t59 * rSges(3,2) + t82) + g(2) * (t58 * rSges(3,1) - t57 * rSges(3,2) + t67) + g(3) * (t70 * rSges(3,3) + t93) + (g(1) * rSges(3,3) * t74 + g(3) * (rSges(3,1) * t73 + rSges(3,2) * t77) + g(2) * (-rSges(3,3) - pkin(6)) * t78) * t69) - m(4) * (g(1) * (t52 * rSges(4,1) - t51 * rSges(4,2) + t94 * t59 + t80) + g(2) * (t50 * rSges(4,1) - t49 * rSges(4,2) + t94 * t57 + t79) + g(3) * (t56 * rSges(4,1) - t55 * rSges(4,2) - t94 * t88 + t81)) - m(5) * (g(1) * (t52 * pkin(3) + t59 * pkin(7) + (t52 * t75 + t59 * t71) * rSges(5,1) + (-t52 * t71 + t59 * t75) * rSges(5,2) + t92 * t51 + t80) + g(2) * (t50 * pkin(3) + t57 * pkin(7) + (t50 * t75 + t57 * t71) * rSges(5,1) + (-t50 * t71 + t57 * t75) * rSges(5,2) + t92 * t49 + t79) + g(3) * (t56 * pkin(3) - pkin(7) * t88 + (t56 * t75 - t71 * t88) * rSges(5,1) + (-t56 * t71 - t75 * t88) * rSges(5,2) + t92 * t55 + t81));
U = t1;
