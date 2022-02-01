% Calculate potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:58:59
% EndTime: 2022-01-23 08:59:00
% DurationCPUTime: 0.34s
% Computational Cost: add. (158->103), mult. (263->147), div. (0->0), fcn. (280->10), ass. (0->47)
t74 = sin(pkin(7));
t98 = t74 * qJ(3) + pkin(1);
t73 = sin(pkin(8));
t97 = t73 * qJ(4) + pkin(2);
t78 = sin(qJ(5));
t96 = t73 * t78;
t80 = cos(qJ(5));
t95 = t73 * t80;
t76 = cos(pkin(8));
t94 = t74 * t76;
t79 = sin(qJ(1));
t93 = t74 * t79;
t81 = cos(qJ(1));
t92 = t74 * t81;
t75 = cos(pkin(9));
t77 = cos(pkin(7));
t91 = t77 * t75;
t90 = t79 * t73;
t89 = t79 * t76;
t88 = t81 * t73;
t87 = t81 * t76;
t86 = -t77 * qJ(3) + pkin(5);
t85 = qJ(4) * t76 - qJ(2);
t84 = rSges(3,1) * t77 - rSges(3,2) * t74 + pkin(1);
t72 = sin(pkin(9));
t83 = t72 * pkin(4) - t75 * pkin(6) + qJ(3);
t68 = t75 * pkin(4) + t72 * pkin(6) + pkin(3);
t82 = -t68 * t73 + t85;
t71 = t79 * qJ(2);
t67 = t76 * pkin(3) + t97;
t66 = pkin(2) * t77 + t98;
t65 = -t73 * pkin(3) + t85;
t64 = t77 * t87 + t90;
t63 = t77 * t88 - t89;
t62 = t77 * t89 - t88;
t61 = t77 * t90 + t87;
t60 = t75 * t96 + t80 * t76;
t59 = t74 * t72 + t76 * t91;
t58 = -t77 * t72 + t75 * t94;
t57 = t72 * t94 + t91;
t56 = t68 * t76 + t97;
t55 = t67 * t77 + t98;
t54 = t64 * t72 - t75 * t92;
t53 = t62 * t72 - t75 * t93;
t52 = -t59 * t78 + t77 * t95;
t51 = t56 * t77 + t83 * t74 + pkin(1);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t81 * rSges(2,1) - t79 * rSges(2,2)) + g(2) * (t79 * rSges(2,1) + t81 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * t71 + g(3) * (t74 * rSges(3,1) + t77 * rSges(3,2) + pkin(5)) + (g(1) * rSges(3,3) + g(2) * t84) * t79 + (g(1) * t84 + g(2) * (-rSges(3,3) - qJ(2))) * t81) - m(4) * (g(1) * (t64 * rSges(4,1) - t63 * rSges(4,2) + t66 * t81 + t71) + g(2) * (t62 * rSges(4,1) - t61 * rSges(4,2) - t81 * qJ(2) + t66 * t79) + g(3) * (-t77 * rSges(4,3) + t86) + (g(3) * (rSges(4,1) * t76 - rSges(4,2) * t73 + pkin(2)) + (g(1) * t81 + g(2) * t79) * rSges(4,3)) * t74) - m(5) * (g(1) * (t55 * t81 - t65 * t79 + (t64 * t75 + t72 * t92) * rSges(5,1) - t54 * rSges(5,2) + t63 * rSges(5,3)) + g(2) * (t55 * t79 + t65 * t81 + (t62 * t75 + t72 * t93) * rSges(5,1) - t53 * rSges(5,2) + t61 * rSges(5,3)) + g(3) * (t58 * rSges(5,1) - t57 * rSges(5,2) + (rSges(5,3) * t73 + t67) * t74 + t86)) - m(6) * (g(1) * (t51 * t81 - t82 * t79 + ((t59 * t80 + t77 * t96) * t81 + t79 * (t75 * t95 - t78 * t76)) * rSges(6,1) + (t52 * t81 - t79 * t60) * rSges(6,2) + t54 * rSges(6,3)) + g(2) * (t51 * t79 + t82 * t81 + ((t59 * t79 - t75 * t88) * t80 + t61 * t78) * rSges(6,1) + (t52 * t79 + t81 * t60) * rSges(6,2) + t53 * rSges(6,3)) + g(3) * (t56 * t74 - t83 * t77 + pkin(5) + (t58 * t80 + t74 * t96) * rSges(6,1) + (-t58 * t78 + t74 * t95) * rSges(6,2) + t57 * rSges(6,3)));
U = t1;
