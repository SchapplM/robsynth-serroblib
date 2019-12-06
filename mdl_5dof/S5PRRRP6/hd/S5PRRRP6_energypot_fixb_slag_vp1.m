% Calculate potential energy for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:41
% EndTime: 2019-12-05 16:50:42
% DurationCPUTime: 0.41s
% Computational Cost: add. (149->81), mult. (195->103), div. (0->0), fcn. (187->8), ass. (0->33)
t92 = rSges(6,1) + pkin(4);
t91 = rSges(4,3) + pkin(6);
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t90 = g(1) * t66 + g(2) * t65;
t89 = rSges(6,3) + qJ(5);
t68 = sin(qJ(2));
t85 = rSges(3,2) * t68;
t67 = sin(qJ(3));
t84 = t65 * t67;
t70 = cos(qJ(2));
t83 = t65 * t70;
t82 = t66 * t67;
t81 = t66 * t70;
t80 = t67 * t70;
t69 = cos(qJ(3));
t79 = t69 * t70;
t76 = t66 * pkin(1) + t65 * pkin(5);
t57 = t69 * pkin(3) + pkin(2);
t71 = -pkin(7) - pkin(6);
t75 = t68 * t57 + t70 * t71 + qJ(1);
t62 = t65 * pkin(1);
t74 = -t66 * pkin(5) + t62;
t73 = pkin(3) * t84 + t57 * t81 + t76;
t72 = -pkin(3) * t82 + t57 * t83 + t74;
t64 = qJ(3) + qJ(4);
t60 = cos(t64);
t59 = sin(t64);
t51 = t65 * t59 + t60 * t81;
t50 = t59 * t81 - t65 * t60;
t49 = -t66 * t59 + t60 * t83;
t48 = t59 * t83 + t66 * t60;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t66 * rSges(2,1) - t65 * rSges(2,2)) + g(2) * (t65 * rSges(2,1) + t66 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t65 * rSges(3,3) + t76) + g(2) * (rSges(3,1) * t83 - t65 * t85 + t62) + g(3) * (t68 * rSges(3,1) + t70 * rSges(3,2) + qJ(1)) + (g(1) * (rSges(3,1) * t70 - t85) + g(2) * (-rSges(3,3) - pkin(5))) * t66) - m(4) * (g(1) * (pkin(2) * t81 + (t66 * t79 + t84) * rSges(4,1) + (t65 * t69 - t66 * t80) * rSges(4,2) + t76) + g(2) * (pkin(2) * t83 + (t65 * t79 - t82) * rSges(4,1) + (-t65 * t80 - t66 * t69) * rSges(4,2) + t74) + g(3) * (-t91 * t70 + qJ(1)) + (g(3) * (rSges(4,1) * t69 - rSges(4,2) * t67 + pkin(2)) + t90 * t91) * t68) - m(5) * (g(1) * (t51 * rSges(5,1) - t50 * rSges(5,2) + t73) + g(2) * (t49 * rSges(5,1) - t48 * rSges(5,2) + t72) + g(3) * (-t70 * rSges(5,3) + t75) + (g(3) * (rSges(5,1) * t60 - rSges(5,2) * t59) + t90 * (rSges(5,3) - t71)) * t68) - m(6) * (g(1) * (t89 * t50 + t92 * t51 + t73) + g(2) * (t89 * t48 + t92 * t49 + t72) + g(3) * (-t70 * rSges(6,2) + t75) + (g(3) * (t89 * t59 + t92 * t60) + t90 * (rSges(6,2) - t71)) * t68);
U = t1;
