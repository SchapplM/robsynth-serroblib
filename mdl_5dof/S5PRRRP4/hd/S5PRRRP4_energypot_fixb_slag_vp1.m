% Calculate potential energy for
% S5PRRRP4
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:32
% EndTime: 2019-12-05 16:45:32
% DurationCPUTime: 0.26s
% Computational Cost: add. (149->73), mult. (167->91), div. (0->0), fcn. (155->8), ass. (0->33)
t84 = rSges(6,1) + pkin(4);
t83 = rSges(6,3) + qJ(5);
t82 = rSges(3,3) + pkin(5);
t61 = qJ(2) + qJ(3);
t58 = sin(t61);
t62 = sin(pkin(8));
t81 = t58 * t62;
t63 = cos(pkin(8));
t80 = t58 * t63;
t59 = cos(t61);
t79 = t59 * t63;
t64 = sin(qJ(4));
t78 = t62 * t64;
t66 = cos(qJ(4));
t77 = t62 * t66;
t76 = t63 * t64;
t75 = t63 * t66;
t67 = cos(qJ(2));
t57 = t67 * pkin(2) + pkin(1);
t68 = -pkin(6) - pkin(5);
t74 = t62 * t57 + t63 * t68;
t65 = sin(qJ(2));
t73 = t65 * pkin(2) + qJ(1);
t72 = t58 * pkin(3) + t73;
t71 = t62 * t59 * pkin(3) + pkin(7) * t81 + t74;
t49 = t63 * t57;
t70 = pkin(3) * t79 + pkin(7) * t80 - t62 * t68 + t49;
t69 = rSges(3,1) * t67 - rSges(3,2) * t65 + pkin(1);
t47 = t59 * t75 + t78;
t46 = t59 * t76 - t77;
t45 = t59 * t77 - t76;
t44 = t59 * t78 + t75;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t63 * rSges(2,1) - t62 * rSges(2,2)) + g(2) * (t62 * rSges(2,1) + t63 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (t65 * rSges(3,1) + t67 * rSges(3,2) + qJ(1)) + (g(1) * t69 - g(2) * t82) * t63 + (g(1) * t82 + g(2) * t69) * t62) - m(4) * (g(1) * (rSges(4,1) * t79 - rSges(4,2) * t80 + t49) + g(2) * (-t63 * rSges(4,3) + t74) + g(3) * (t58 * rSges(4,1) + t59 * rSges(4,2) + t73) + (g(1) * (rSges(4,3) - t68) + g(2) * (rSges(4,1) * t59 - rSges(4,2) * t58)) * t62) - m(5) * (g(1) * (t47 * rSges(5,1) - t46 * rSges(5,2) + rSges(5,3) * t80 + t70) + g(2) * (t45 * rSges(5,1) - t44 * rSges(5,2) + rSges(5,3) * t81 + t71) + g(3) * ((-rSges(5,3) - pkin(7)) * t59 + (rSges(5,1) * t66 - rSges(5,2) * t64) * t58 + t72)) - m(6) * (g(1) * (t83 * t46 + t47 * t84 + t70) + g(2) * (t83 * t44 + t84 * t45 + t71) + g(3) * (t72 + (-rSges(6,2) - pkin(7)) * t59) + (g(3) * (t83 * t64 + t66 * t84) + (g(1) * t63 + g(2) * t62) * rSges(6,2)) * t58);
U = t1;
