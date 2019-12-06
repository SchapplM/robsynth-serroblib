% Calculate potential energy for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:28
% EndTime: 2019-12-05 15:08:28
% DurationCPUTime: 0.25s
% Computational Cost: add. (149->73), mult. (167->91), div. (0->0), fcn. (155->8), ass. (0->33)
t84 = rSges(6,1) + pkin(4);
t83 = rSges(6,3) + qJ(5);
t61 = pkin(8) + qJ(3);
t58 = sin(t61);
t63 = sin(pkin(7));
t82 = t58 * t63;
t65 = cos(pkin(7));
t81 = t58 * t65;
t59 = cos(t61);
t80 = t59 * t65;
t67 = sin(qJ(4));
t79 = t63 * t67;
t68 = cos(qJ(4));
t78 = t63 * t68;
t77 = t65 * t67;
t76 = t65 * t68;
t64 = cos(pkin(8));
t57 = t64 * pkin(2) + pkin(1);
t66 = -pkin(5) - qJ(2);
t75 = t63 * t57 + t65 * t66;
t74 = rSges(3,3) + qJ(2);
t62 = sin(pkin(8));
t73 = t62 * pkin(2) + qJ(1);
t72 = t58 * pkin(3) + t73;
t71 = t63 * t59 * pkin(3) + pkin(6) * t82 + t75;
t49 = t65 * t57;
t70 = pkin(3) * t80 + pkin(6) * t81 - t63 * t66 + t49;
t69 = rSges(3,1) * t64 - rSges(3,2) * t62 + pkin(1);
t47 = t59 * t76 + t79;
t46 = t59 * t77 - t78;
t45 = t59 * t78 - t77;
t44 = t59 * t79 + t76;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t65 * rSges(2,1) - t63 * rSges(2,2)) + g(2) * (t63 * rSges(2,1) + t65 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (t62 * rSges(3,1) + t64 * rSges(3,2) + qJ(1)) + (g(1) * t69 - g(2) * t74) * t65 + (g(1) * t74 + g(2) * t69) * t63) - m(4) * (g(1) * (rSges(4,1) * t80 - rSges(4,2) * t81 + t49) + g(2) * (-t65 * rSges(4,3) + t75) + g(3) * (t58 * rSges(4,1) + t59 * rSges(4,2) + t73) + (g(1) * (rSges(4,3) - t66) + g(2) * (rSges(4,1) * t59 - rSges(4,2) * t58)) * t63) - m(5) * (g(1) * (t47 * rSges(5,1) - t46 * rSges(5,2) + rSges(5,3) * t81 + t70) + g(2) * (t45 * rSges(5,1) - t44 * rSges(5,2) + rSges(5,3) * t82 + t71) + g(3) * ((-rSges(5,3) - pkin(6)) * t59 + (rSges(5,1) * t68 - rSges(5,2) * t67) * t58 + t72)) - m(6) * (g(1) * (t83 * t46 + t84 * t47 + t70) + g(2) * (t83 * t44 + t84 * t45 + t71) + g(3) * (t72 + (-rSges(6,2) - pkin(6)) * t59) + (g(3) * (t83 * t67 + t84 * t68) + (g(1) * t65 + g(2) * t63) * rSges(6,2)) * t58);
U = t1;
